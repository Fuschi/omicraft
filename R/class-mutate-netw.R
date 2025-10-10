#' Mutate Network/Community Slots in an `omic` Object
#'
#' Applies user-defined expressions to compute new columns for the `taxa` slot,
#' using the network (`netw`) and community (`comm`) slots of an `omic` object.
#' Expressions are evaluated in a data mask where the symbols
#' `netw` and `comm` are bound to the current object's slots—or, if the object
#' is grouped via `group_omic()`, to the corresponding slots of each subgroup-
#' specific subset.
#'
#' @param object An `omic` object.
#'   This function targets the `netw` (network graph) and `comm` (community/cluster)
#'   slots for computations, and writes results as new columns into `taxa(object)`.
#' @param ... One or more expressions to be evaluated. Within each expression, the
#'   symbols `netw` and `comm` refer to the network and community slots respectively.
#'
#' @details
#' - **Grouping semantics:** This function respects the grouping state set with
#'   `group_omic()`. The taxa-level grouping variables are discovered through
#'   `taxa_vars(object)`. If no taxa grouping variables are present, expressions
#'   are evaluated once globally; otherwise they are evaluated *per subgroup*
#'   on a subsetted object, and results are stitched back in the original order.
#' - **No string rewriting:** expressions are evaluated using a `data` mask so the
#'   symbols `netw` and `comm` resolve directly to the intended slots, avoiding
#'   fragile text substitution.
#' - **Selected links:** if `select_link()` has been used earlier, only those edges
#'   are present during the computation (via a temporary subgraph), then the original
#'   network is restored.
#' - **Result length:** each expression must return either a scalar or a vector of
#'   length equal to the number of taxa in scope (whole object or subgroup).
#'
#' @return The modified `omic` object, with new columns appended to `taxa(object)`.
#'
#' @export
#' @importFrom dplyr mutate select pull filter
#' @importFrom tidyr unite
#' @importFrom rlang enquos get_expr expr_text eval_tidy current_env :=
#' @importFrom cli cli_abort cli_warn
#' @importFrom tidyselect all_of
#' @name mutate_netw
#' @aliases mutate_netw,omic-method mutate_netw,omics-method
setGeneric("mutate_netw", function(object, ...) standardGeneric("mutate_netw"))

setMethod("mutate_netw", "omic", function(object, ...) {
  
  # 0) Basic availability checks ------------------------------------------------
  if (miss_slot(object, "netw") && miss_slot(object, "comm")) {
    cli::cli_abort("No network and communities available in the {.cls omic} object.")
  }
  
  # 1) Capture user expressions and validate reserved keywords -----------------
  quosures <- rlang::enquos(...)
  if (length(quosures) == 0L) cli::cli_abort("{.fn mutate_netw} requires at least one expression in `...`.")
  check_reserved_keywords(quosures)
  
  # 2) Read taxa grouping from group_omic() state ------------------------------
  #    We assume taxa_vars(object) returns the taxa-level grouping variables
  taxa_groups <- setdiff(get_group_omic(object), meta_vars(object))
  if (length(taxa_groups) == 0) taxa_groups <- NULL
  
  # 3) Temporarily filter edges if user has selected links ---------------------
  selected_links <- get_selected_links(object)
  if (!is.null(selected_links)) {
    netw0 <- netw(object)
    netw(object) <- igraph::subgraph_from_edges(
      graph = netw(object),
      eids  = selected_links,
      delete.vertices = FALSE
    )
  }
  
  # 4) Determine names for the new output columns ------------------------------
  expressions_text <- purrr::map_chr(quosures, ~ rlang::expr_text(rlang::get_expr(.x)))
  quosures_names  <- rlang::names2(quosures)  
  quosures_names  <- ifelse(quosures_names == "", expressions_text, quosures_names)
  
  # 5) Evaluate expressions (no string replacement; use a data mask) -----------
  if (is.null(taxa_groups)) {
    # Ungrouped case: evaluate with a mask that binds `netw` and `comm`
    mask <- list(
      netw = netw(object),
      comm = comm(object)
    )
    
    for (i in seq_along(quosures)) {

      val <- rlang::eval_tidy(quosures[[i]], data = mask, env = rlang::current_env())

      # Recycle scalar or validate vector length
      if (length(val) == 1L) {
        val <- rep(val, ntaxa(object))
      } else if (length(val) != ntaxa(object)) {
        cli::cli_abort(
          "Expression '{quosures_names[i]}' must return length 1 or {ntaxa(object)}, not {length(val)}."
        )
      }
      
      taxa(object) <- taxa(object, .fmt = "tbl") |>
        dplyr::mutate(!!quosures_names[i] := val)
    }
    
  } else {
    # Grouped case: evaluate per subgroup on a subset object --------------------
    
    # Build subgroup keys from taxa metadata
    subgroups <- taxa(object, .fmt = "tbl") |>
      dplyr::select(tidyselect::all_of(taxa_groups)) |>
      tidyr::unite("_internal_") |>
      dplyr::pull("_internal_")
    
    unique_keys <- unique(subgroups)
    
    for (i in seq_along(quosures)) {
      # Store per-taxa results; allow per-group scalars or vectors
      result <- vector(length = ntaxa(object))
      
      for (key in unique_keys) {
        idx_key <- which(subgroups == key)
        object_subset <- object[, idx_key]
        
        # Per-group mask: bind netw/comm of the subset
        mask_subset <- list(
          netw = netw(object_subset),
          comm = comm(object_subset)
        )
        
        val_key <- rlang::eval_tidy(quosures[[i]], data = mask_subset, env = rlang::current_env())
        
        # Recycle scalar or validate vector length for this group
        if (length(val_key) == 1L) {
          result[idx_key] <- rep(val_key, length(idx_key))
        } else if (length(val_key) == length(idx_key)) {
          result[idx_key] <- val_key
        } else {
          cli::cli_abort(
            "Expression '{quosures_names[i]}' returned length {length(val_key)} for a group of size {length(idx_key)}."
          )
        }
      }
      
      taxa(object) <- taxa(object, .fmt = "tbl") |>
        dplyr::mutate(!!quosures_names[i] := result)
    }
  }
  
  # 6) Restore the original network if we temporarily filtered edges -----------
  if (!is.null(selected_links)) netw(object) <- netw0
  
  # 7) Validate and return ------------------------------------------------------
  validObject(object)
  object
})


setMethod("mutate_netw", "omics", function(object, ...) {
  
  # 0) Basic availability checks ------------------------------------------------
  if (miss_slot(object, "netw", "any") && miss_slot(object, "comm", "any")) {
    cli::cli_abort("No network and communities available in at least one {.cls omic} object.")
  }
  
  # 1) Capture user expressions and validate reserved keywords -----------------
  quosures <- rlang::enquos(...)
  if (length(quosures) == 0L) cli::cli_abort("{.fn mutate_netw} requires at least one expression in `...`.")
  check_reserved_keywords(quosures)
  
  # 2) Read taxa grouping from group_omic() state ------------------------------
  #    We assume taxa_vars(object) returns the taxa-level grouping variables
  taxa_groups <- setdiff(get_group_omic(object), meta_vars(object))
  if (length(taxa_groups) == 0) taxa_groups <- NULL
  
  # 3) Temporarily filter edges if user has selected links ---------------------
  selected_links <- get_selected_links(object)
  if (all(all(purrr::map_lgl(selected_links, \(x) !is.null(x))))) {
    netw0 <- netw(object)
    for(mg in 1:length(object)){
      netw(object[[mg]]) <- igraph::subgraph_from_edges(graph = netw(object[[mg]]),
                                                        eids = get_selected_links(object[[mg]]),
                                                        delete.vertices = FALSE)
    }
  }
  
  # 4) Determine names for the new output columns ------------------------------
  expressions_text <- purrr::map_chr(quosures, ~ rlang::expr_text(rlang::get_expr(.x)))
  quosures_names  <- rlang::names2(quosures)  
  quosures_names  <- ifelse(quosures_names == "", expressions_text, quosures_names)
  
  # 5) Evaluate expressions (use a data mask) ----------------------------------
  if (is.null(taxa_groups)) {
    # Caso non raggruppato: per ciascun mg creiamo un mask con netw/comm
    for (i in seq_along(quosures)) {
      for (mg in seq_along(object)) {
        mask <- list(
          netw = netw(object[[mg]]),
          comm = comm(object[[mg]])
        )
        val <- rlang::eval_tidy(quosures[[i]], data = mask, env = rlang::current_env())
        taxa(object[[mg]]) <- taxa(object[[mg]], .fmt = "tbl") %>%
          dplyr::mutate(!!quosures_names[i] := val)
      }
    }
    
  } else {
    # Caso raggruppato: ricaviamo i sottogruppi e valutiamo su subset per mg ----
    subgroups <- taxa(object, .collapse = TRUE) %>%
      dplyr::select(tidyselect::all_of(c("omic", taxa_groups))) %>%
      tidyr::unite("_internal_", remove = FALSE) %>%
      split(.data[["omic"]])
    subgroups <- sapply(subgroups, \(x) dplyr::pull(x, "_internal_"), USE.NAMES = TRUE, simplify = FALSE)
    
    for (mg in seq_along(object)) {
      unique_keys <- unique(subgroups[[mg]])
      for (i in seq_along(quosures)) {
        result <- vector(length = ntaxa(object[[mg]]))
        for (key in unique_keys) {
          idx_key <- which(subgroups[[mg]] == key)
          object_subset <- object[[mg]][, idx_key]
          mask_subset <- list(
            netw = netw(object_subset),
            comm = comm(object_subset)
          )
          result_key <- rlang::eval_tidy(quosures[[i]], data = mask_subset, env = rlang::current_env())
          result[idx_key] <- result_key
        }
        taxa(object[[mg]]) <- taxa(object[[mg]], .collapse = TRUE) %>%
          dplyr::mutate(!!quosures_names[i] := result)
      }
    }
  }
  
  if (all(all(purrr::map_lgl(selected_links, \(x) !is.null(x))))) {
    for(mg in seq_along(object)) netw(object[[mg]]) <- netw0[[mg]]
  }
  validObject(object)
  return(object)
  
})