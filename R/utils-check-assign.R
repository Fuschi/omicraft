# Internal function to check the assigned elements for the setter methods of omics
#------------------------------------------------------------------------------#

#' Check list in assign methods for omics
#'
#' @description
#' Internal validator ensuring that a named list supplied to an `omics`
#' setter has the correct length and name set (order may differ).
#'
#' @param object An `omics` object.
#' @param value  A named list intended to be assigned to the `omics`.
#'
#' @return `invisible(TRUE)` if all checks pass; otherwise throws an error.
#' @keywords internal
#' @export
is_list_omics_assign <- function(object, value) {
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  
  # 1) Type check on object ----------------------------------------------------
  if (!inherits(object, "omics")) {
    cli::cli_abort("{.arg {objectName}} must be an {.cls omics} object.")
  }
  
  # 2) Basic type/length checks on value ---------------------------------------
  if (!is.list(value)) {
    cli::cli_abort("{.arg {valueName}} must be a list.")
  }
  
  obj_len <- length(object@omics)
  val_len <- length(value)
  if (obj_len != val_len) {
    cli::cli_abort(
      "Lengths must match: {obj_len} in {.arg {objectName}}@omics vs {val_len} in {.arg {valueName}}."
    )
  }
  
  # 3) Names: must exist, non-empty, no duplicates -----------------------------
  val_names <- names(value)
  if (is.null(val_names) || any(is.na(val_names)) || any(val_names == "")) {
    cli::cli_abort("{.arg {valueName}} must be a named list with non-empty names.")
  }
  
  dup_idx <- anyDuplicated(val_names)
  if (dup_idx) {
    dups <- unique(val_names[duplicated(val_names)])
    dups_fmt <- paste0("{", paste(dups, collapse = "}, {"), "}")
    cli::cli_abort("{.arg {valueName}} cannot contain duplicated names: {dups_fmt}.")
  }
  
  # 4) Name set equality (order can differ) ------------------------------------
  obj_names <- names(object@omics)
  if (!is.null(obj_names)) {
    miss_in_value  <- setdiff(obj_names, val_names)
    extra_in_value <- setdiff(val_names, obj_names)
    
    if (length(miss_in_value) || length(extra_in_value)) {
      miss_fmt  <- if (length(miss_in_value)) paste(miss_in_value, collapse = ", ") else "none"
      extra_fmt <- if (length(extra_in_value)) paste(extra_in_value, collapse = ", ") else "none"
      cli::cli_abort(c(
        "Names of {.arg {valueName}} must match names of {.arg {objectName}}@omics (order can differ).",
        "x" = "Missing in {.arg {valueName}}: {miss_fmt}",
        "x" = "Extra in {.arg {valueName}}: {extra_fmt}"
      ))
    }
  }
  
  invisible(TRUE)
}



#' Check Tibble in Assign Methods for omics
#'
#' @description
#' Internal validator for a tibble/data frame replacing all sample/taxa
#' metadata across an `omics` object. Requires `omic` and `sample_id`/`taxa_id`
#' columns and forbids row names.
#'
#' @param object An `omics` object.
#' @param value  A tibble/data frame with columns `omic` and either
#'   `sample_id` (if `sample_or_taxa = "sample"`) or `taxa_id` (if `"taxa"`).
#' @param sample_or_taxa Character, `"sample"` or `"taxa"`.
#'
#' @return `invisible(TRUE)` if all checks pass; otherwise throws an error.
#' @keywords internal
#' @importFrom cli cli_abort
#' @importFrom tibble has_rownames
#' @export
is_assign_omics_tbl <- function(object, value, sample_or_taxa) {
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  sample_or_taxa <- match.arg(sample_or_taxa, c("sample", "taxa"))
  
  # 1) Object type -------------------------------------------------------------
  if (!inherits(object, "omics")) {
    cli::cli_abort("{.arg {objectName}} must be an {.cls omics} object.")
  }
  
  # 2) Value type --------------------------------------------------------------
  if (!(inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg {valueName}} must be a tibble or data frame.")
  }
  
  # 3) Required columns --------------------------------------------------------
  required_cols <- if (sample_or_taxa == "sample") c("omic", "sample_id") else c("omic", "taxa_id")
  value_cols <- colnames(value)
  missing_cols <- setdiff(required_cols, value_cols)
  if (length(missing_cols) > 0L) {
    cli::cli_abort(c(
      "Missing required columns in {.arg {valueName}}.",
      "x" = "Required: {paste(required_cols, collapse = ', ')}",
      "x" = "Missing: {paste(missing_cols, collapse = ', ')}"
    ))
  }
  
  # 4) Rownames must NOT be set ------------------------------------------------
  if (tibble::has_rownames(value)) {
    cli::cli_abort("Row names must not be set in {.arg {valueName}}.")
  }
  
  # 5) 'omic' values must match names(object) ----------------------------------
  unique_omics  <- unique(value$omic)
  object_omics  <- names(object)
  missing_omics <- setdiff(unique_omics, object_omics)
  if (length(missing_omics) > 0L) {
    cli::cli_abort(c(
      "Some {.val omic} values in {.arg {valueName}} are not present in {.arg {objectName}}.",
      "x" = "Not found: {paste(missing_omics, collapse = ', ')}",
      "i" = "Allowed values: {paste(object_omics, collapse = ', ')}"
    ))
  }
  
  invisible(TRUE)
}