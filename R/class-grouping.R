#' @include class-omic.R class-omics.R class-base-methods.R
NULL

#==============================================================================#
# Grouping for `omic` / `omics`
#==============================================================================#

#' Grouping helpers for `omic` / `omics`
#'
#' Set or get grouping variables on `omic` / `omics` objects, with semantics
#' similar to `dplyr::group_by()` / `dplyr::ungroup()`.
#'
#' @details
#' The grouping state is stored as an attribute named "omic_groups" on the
#' object (for both `omic` and `omics`).
#'
#' Semantics:
#' * `group_omic(x, a, b)` replaces the current groups with `c("a","b")`.
#' * `group_omic(x, a, .add = TRUE)` adds `"a"` to the existing groups.
#' * `group_omic(x)` with no grouping vars returns `x` unchanged (keep current).
#' * `group_omic(x, NULL)` clears the groups (like `group_by(NULL)`).
#' * `ungroup_omic(x)` clears all groups.
#' * `ungroup_omic(x, a)` removes `"a"` from the current groups (if present).
#'
#' Validation:
#' For `omic`, valid fields come from `meta_taxa_vars(object)`.
#' For `omics`, valid fields come from `meta_taxa_vars(object, .fmt = "unique")`.
#'
#' @param object An `omic` or `omics` object.
#' @param ... One or more variable names (bare or strings). May also include a
#' single `NULL` to drop all groups. In `ungroup_omic()`, names in `...` are
#' removed from the current grouping vector.
#' @param .add Logical; if `TRUE`, add to existing groups; if `FALSE`, replace.
#'
#' @return
#' - `group_omic()` returns the modified object.
#' - `ungroup_omic()` returns the modified object.
#' - `get_group_omic()` returns a character vector (possibly length-0).
#'
#' @name group-omic
NULL

# -----------------------------------------------------------------------------#
# Getters
# -----------------------------------------------------------------------------#

#' @rdname group-omic
#' @export
setGeneric("get_group_omic", function(object) standardGeneric("get_group_omic"))

#' @rdname group-omic
#' @export
setMethod("get_group_omic", "omic", function(object) {return(attr(object, "omic_groups"))})

#' @rdname group-omic
#' @export
setMethod("get_group_omic", "omics", function(object) {return(attr(object, "omic_groups"))})

#' @importFrom rlang enquos as_name quo_is_null
NULL

# -----------------------------------------------------------------------------#
# Setters (group_omic)
# -----------------------------------------------------------------------------#

#' @noRd
.group_omic <- function(object, ..., .add = FALSE){
  qs <- rlang::enquos(..., .ignore_empty = "all")
  
  # No arguments or NULL remove all groups
  if (length(qs) == 0L || (length(qs) == 1L && rlang::quo_is_null(qs[[1L]]))) {
    attr(object, "omic_groups") <- character(0)
    return(object)
  }
  
  # Convert quosures in strings
  vars <- unique(vapply(qs, rlang::as_name, character(1)))
  
  # Validate the arguments
  available <- meta_taxa_vars(object, .fmt = "unique")
  missing_vars <- setdiff(vars, available)
  if (length(missing_vars)) {
    cli::cli_abort(c(
      "x" = "Unknown grouping variable{?s}: {.val {missing_vars}}.",
      "v" = "Available variables are: {.val {available}}."
    ))
  }
  
  if (isTRUE(.add)) {
    current <- get_group_omic(object)
    attr(object, "omic_groups") <- unique(c(current, vars))
  } else {
    attr(object, "omic_groups") <- vars
  }
  object
}

#' @rdname group-omic
#' @export
setGeneric("group_omic", function(object, ..., .add = FALSE) standardGeneric("group_omic"))

#' @rdname group-omic
#' @export
setMethod("group_omic", signature(object = "omic"),
          function(object, ..., .add = FALSE) .group_omic(object, ..., .add = .add))

#' @rdname group-omic
#' @export
setMethod("group_omic", signature(object = "omics"),
          function(object, ..., .add = FALSE) .group_omic(object, ..., .add = .add))

# -----------------------------------------------------------------------------#
# Ungroup
# -----------------------------------------------------------------------------#

.ungroup_omic <- function(object, ...){
  qs <- rlang::enquos(..., .ignore_empty = "all")
  if (length(qs) == 0L) {
    attr(object, "omic_groups") <- character(0)   # ungroup all
    return(object)
  }
  drop <- unique(vapply(qs, rlang::as_name, character(1)))
  keep <- setdiff(get_group_omic(object), drop)
  attr(object, "omic_groups") <- keep
  object
}

#' @rdname group-omic
#' @export
setGeneric("ungroup_omic", function(object, ...) standardGeneric("ungroup_omic"))

#' @rdname group-omic
#' @export
setMethod("ungroup_omic", signature(object = "omic"), function(object, ...) {
  .ungroup_omic(object, ...)
})

#' @rdname group-omic
#' @export
setMethod("ungroup_omic", signature(object = "omics"), function(object, ...) {
  .ungroup_omic(object, ...)
})
