#' @include class-omic.R class-omics.R class-base-methods.R
NULL

# PULL INFO SAMPLE
#------------------------------------------------------------------------------#
#' Pull Sample Information from omic or omics Objects
#'
#' @description
#' Retrieves specific sample information from the `meta` data frame within `omic` objects,
#' or each `omic` object within an `omics`. This function simplifies direct access to specific columns of interest
#' using dynamic column name handling.
#'
#' @param object An `omic` or `omics` object.
#' @param var The name or position of the column to retrieve from the `meta` data frame. 
#'        If -1 (default), the last column of the data frame is returned. You can specify the column name unquoted due to non-standard evaluation.
#' @return For a single `omic` object, a vector containing the data from the specified column 
#'         of the `meta` data frame is returned. For an `omics` object, a list of such vectors
#'         is returned, each corresponding to one `omic` object in the list.
#'
#' @details
#' This function supports dynamic evaluation of the column name using `rlang` for unquoted names.
#'
#' @export
#' @importFrom dplyr pull
#' @importFrom rlang ensym
#' @importFrom purrr map
#' @name pull_meta
#' @aliases pull_meta,omic-method pull_meta,omics-method
setGeneric("pull_meta", function(object, var = -1) standardGeneric("pull_meta"))

setMethod("pull_meta", signature(object = "omic"), function(object, var = -1) {
  if(miss_sample(object)) {stop("Error: No sample available.")}
  var <- rlang::ensym(var)
  dplyr::pull(meta(object, .fmt = "tbl"), var)
})

setMethod("pull_meta", signature(object = "omics"), function(object, var = -1) {
  if(miss_sample(object, "any")) {stop("Error: No sample available in any of the omic objects.")}
  var <- rlang::ensym(var)
  meta(object, .collapse = TRUE) %>% 
    dplyr::select(omic, var) %>% 
    dplyr::group_split(.data$omic) %>%
    purrr::map(\(x) pull(x, !!var))
  
})


# PULL TAXA
#------------------------------------------------------------------------------#
#' Pull Taxa Information from omic or omics Objects
#'
#' @description
#' Retrieves specific taxonomic information from the `taxa` data frame within `omic` objects,
#' or each `omic` object within an `omics`. This function simplifies direct access to specific columns of interest
#' using dynamic column name handling.
#'
#' @param object An `omic` or `omics` object.
#' @param var The name or position of the column to retrieve from the `taxa` data frame. 
#'        If -1 (default), the last column of the data frame is returned. You can specify the column name unquoted due to non-standard evaluation.
#' @return For a single `omic` object, a vector containing the data from the specified column 
#'         of the `taxa` data frame is returned. For an `omics` object, a list of such vectors
#'         is returned, each corresponding to one `omic` object in the list.
#'
#' @details
#' This function supports dynamic evaluation of the column name using `rlang` for unquoted names, allowing more
#' flexible and intuitive usage within data manipulation workflows.
#'
#' @export
#' @importFrom dplyr pull
#' @importFrom rlang ensym
#' @importFrom purrr map
#' @name pull_taxa
#' @aliases pull_taxa,omic-method pull_taxa,omics-method
setGeneric("pull_taxa", function(object, var = -1) standardGeneric("pull_taxa"))

setMethod("pull_taxa", signature(object = "omic"), function(object, var = -1) {
  if(miss_taxa(object)) {stop("Error: No taxa available.")}
  var <- rlang::ensym(var)
  dplyr::pull(taxa(object, .fmt = "tbl"), var)
})

setMethod("pull_taxa", signature(object = "omics"), function(object, var = -1) {
  if(miss_taxa(object, "any")) {stop("Error: No taxa available in any of the omic objects.")}
  var <- rlang::ensym(var)
  taxa(object, .collapse = TRUE) %>% 
    dplyr::select(omic, var) %>% 
    dplyr::group_split(.data$omic) %>%
    purrr::map(\(x) pull(x, !!var))
  
})