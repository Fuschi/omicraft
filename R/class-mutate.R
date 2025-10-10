#' Modify and Augment `omic` Objects by Transforming the `meta` Slot
#'
#' This function dynamically manipulates the `sample` slot within `omic` or `omics` objects,
#' applying user-defined transformations. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `omic` or `omics` object.
#'        The function targets the `sample` slot, which contains metadata for each sample.
#' @param ... Dynamic expressions or functions to be applied to the data.
#'        These expressions can manipulate both abundance data (e.g., 'abun', 'rela', 'norm') and
#'        metadata within the sample slot. This allows for a comprehensive data transformation
#'        experience that supports all standard and custom `tidyverse` manipulation techniques.
#'
#' @details The function is designed to integrate seamlessly with the `tidyverse`, allowing users
#'          to utilize familiar and potent data manipulation verbs such as `mutate`, `filter`.
#'          It supports using any `tidyverse`-compatible expressions, including conditional operations,
#'          summarizations, and complex transformations involving both abundance and metadata fields.
#'          This flexibility makes it particularly useful for ecological and biological data analysis,
#'          where combining different data types and conditions is common.
#'
#'          ### Keywords in `omic` and `omics`:
#'          - **abun, rela, norm**: Slots within `omic` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **sample_id**: An essential identifier used to uniquely reference individual samples within an `omic` object. 
#'          - **omic**: Used exclusively within `omics` objects to differentiate between multiple `omic` objects 
#'            contained in the list.
#'            
#' @return Returns the `omic` or `omics` object with updated `sample` slots reflecting the applied transformations.
#'         All other structures within the object remain unchanged, ensuring that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_meta,omic-method mutate_meta,omics-method
#' @importFrom dplyr mutate group_by ungroup distinct relocate arrange
#' @importFrom tidyr expand_grid 
#' @importFrom tidyselect any_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_meta", function(object, ...) {standardGeneric("mutate_meta")})

setMethod("mutate_meta", "omic", function(object, ...) {
  
  # 1) Quick check for samples
  if (miss_sample(object)) {
    stop("Error: No samples available in the 'omic' object.")
  }
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for meta & abundance
  meta_tbl   <- meta(object, "tbl")
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 5) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    
    # Determine grouping columns
    meta_groups <- setdiff(get_group_omic(object), taxa_vars(object))
    if(length(meta_groups) == 0 && is_abun_expr) meta_groups <- "sample_id"
    
    # Apply mutate logic
    if (is_abun_expr) {
      # Abundance-related expression => join abundance first
      meta_tbl <- long_abun %>%
        dplyr::left_join(meta_tbl, by = "sample_id") %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct()
    } else {
      # Non-abundance expression => mutate in-place
      meta_tbl <- meta_tbl %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
  }
  
  # --------------------------------------------------------------------------
  # 6) Overwrite meta slot and return
  # --------------------------------------------------------------------------
  meta(object) <- meta_tbl
  object
})