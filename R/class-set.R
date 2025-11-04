#' Set slots of an `omic` object using tidy-style expressions
#'
#' @description
#' `set_omic()` lets you compute and assign one or more slots of an `omic`
#' object using rlang-style expressions that can reference existing slots by
#' their bare names (e.g., `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`).
#' Expressions are evaluated **sequentially** in the order provided, so later
#' expressions can depend on values assigned by earlier ones.
#'
#' @param object An object of class `omic`.
#' @param ... Named expressions, where each name must be one of the allowed
#'   slot names: `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`.
#'
#' @return The modified `omic` object with updated slots.
#'
#' @section How name resolution works:
#' Within each expression, the symbols `abun`, `rela`, `norm`, `meta`, `taxa`,
#' `netw`, and `comm` resolve to the **current** working values (initially taken
#' from the object, then updated as each expression is assigned).
#'
#' @export
setGeneric("set_omic", function(object, ...) standardGeneric("set_omic"))

#' @rdname set_omic
#' @export
setMethod("set_omic", signature(object = "omic"), function(object, ...) {
  # Capture user expressions as quosures, preserving their environment
  dots <- rlang::enquos(..., .named = TRUE)
  
  # Allowed target slots
  allowed <- c("abun", "rela", "norm", "meta", "taxa", "netw", "comm")
  
  # 1) Basic validations on `...`
  if (length(dots) == 0L) {
    cli::cli_abort("No expressions supplied to {.fn set_omic}.")
  }
  
  nm <- names(dots)
  if (is.null(nm) || any(nm == "")) {
    cli::cli_abort("All arguments to {.fn set_omic} must be named.")
  }
  unknown <- setdiff(nm, allowed)
  if (length(unknown)) {
    cli::cli_abort(
      c(
        "Some targets are not valid {'.omic'} slots:",
        "x" = paste0(unknown, collapse = ", "),
        "i" = "Valid targets are: {paste(allowed, collapse = ', ')}"
      )
    )
  }
  
  # 2) Build the initial data mask from current object values
  mask <- list(
    abun = abun(object),
    rela = rela(object),
    norm = norm(object),
    meta = meta(object),
    taxa = taxa(object),
    netw = netw(object),
    comm = comm(object)
  )
  
  # 3) Evaluate each expression with rlang::eval_tidy in the mask
  for (target in nm) {
    quo <- dots[[target]]
    
    value <- tryCatch(
      rlang::eval_tidy(quo, data = mask, env = rlang::get_env(quo)),
      error = function(e) {
        cli::cli_abort(c(
          "Failed while computing {.field {target}} in {.fn set_omic}.",
          "x" = conditionMessage(e),
          "i" = "Within expressions you can reference: {paste(allowed, collapse = ', ')}"
        ))
      }
    )
    
    # Assign via your replacement accessors and refresh the mask
    object <- do.call(paste0(target, "<-"), list(object, value))
    mask[[target]] <- value
  }
  
  # 4) Return the updated object (validity is handled inside your setters)
  object
})


#' @rdname set_omic
#' @export
setMethod("set_omic", signature(object = "omics"), function(object, ...) {

  dots <- rlang::enquos(..., .named = TRUE)
  allowed <- c("abun", "rela", "norm", "meta", "taxa", "netw", "comm")
  
  if (length(dots) == 0L) {
    cli::cli_abort("No expressions supplied to {.fn set_omic}.")
  }
  
  nm <- names(dots)
  if (is.null(nm) || any(nm == "")) {
    cli::cli_abort("All arguments to {.fn set_omic} must be named.")
  }
  unknown <- setdiff(nm, allowed)
  if (length(unknown)) {
    cli::cli_abort(
      c(
        "Some targets are not valid {'.omic'} slots:",
        "x" = paste0(unknown, collapse = ", "),
        "i" = "Valid targets are: {paste(allowed, collapse = ', ')}"
      )
    )
  }
  
  # Apply set_omic() to each element of omics
  n <- length(object)
  if (n == 0L) return(object)
  
  for (i in seq_len(n)) {
    el_name <- names(object)[i]
    label <- if (!is.null(el_name) && nzchar(el_name)) {
      paste0("element '", el_name, "' (index ", i, ")")
    } else {
      paste0("element at index ", i)
    }
    
    om <- object[[i]]
    if (!methods::is(om, "omic")) {
      cli::cli_abort(c(
        "All elements of an {.cls omics} object must be {.cls omic}.",
        "x" = paste0("Found element of class {.cls ", class(om)[1], "} at index ", i, ".")
      ))
    }
    
    # Excecute set_omic on the current element with original quosures
    object[[i]] <- tryCatch(
      rlang::inject(set_omic(om, !!!dots)),
      error = function(e) {
        cli::cli_abort(c(
          "Failed while running {.fn set_omic} on {.val {label}}.",
          "x" = conditionMessage(e)
        ))
      }
    )
  }
  
  object
})

