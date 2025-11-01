#' Summarize `omic` and `omics` Objects in Long Format
#'
#' Converts `omic` or `omics` objects into a longer format including abundances 
#' (absolute, relative, and normalized) and metadata into a single structured tibble. 
#' This transformation facilitates the use of tidyverse tools for
#' subsequent data manipulation and analysis.
#'
#' @param object An `omic` or `omics` object.
#'        For `omic`, this method consolidates sample and taxa information along with 
#'        abundance data and associated metadata into a long-format tibble. For `omics`,
#'        it applies the transformation to each `omic` object within the list and combines the results
#'        into a single tibble, incorporating an additional column named `omic` that identifies the source
#'        `omic` object.
#'
#' @details
#' The `omic_longer` method is tailored to prepare `omic` data for easy integration with tidyverse 
#' workflows by:
#' - Transforming data into a long format based on `sample_id` and `taxa_id`, making it more accessible for
#'   various tidyverse functions such as `mutate`, `filter`, and `summarize`.
#' - Aggregating both abundance metrics (e.g., `abun`, `rela`, `norm`) and metadata associated (e.g meta, taxa, comm_id) 
#'   with samples and taxa, which allows comprehensive data analysis within a single data frame structure.
#' - In the context of `omics` objects, besides the standard data transformation applied to each `omic`,
#'   an additional column using the reserved keyword `omic` is included to keep track of the original `omic` 
#'   object each row of data pertains to, facilitating analysis across multiple datasets.
#'
#' @return A tibble in a long format that consolidates all the relevant data fields under sample and taxon identifiers.
#'         For `omics` objects, the tibble includes an extra column `omic` that indicates the source `omic` object.
#'
#' @export
#' @aliases omic_longer,omic-method omic_longer,omics-method
setGeneric("omic_longer", function(object) standardGeneric("omic_longer"))


setMethod("omic_longer", "omic", function(object){
  
  # empty
  if(miss_sample(object) || miss_taxa(object)) cli::cli_abort("Error: No sample or taxa available.")
  
  
  long_omic <- tidyr::expand_grid(sample_id = sample_id(object),
                                   taxa_id = taxa_id(object))
  if(has_slot(object, "abun")) {
    long_omic <- long_omic %>% 
      dplyr::left_join(
        abun(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "abun"),
        dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(has_slot(object, "rela")) {
    long_omic <- long_omic %>% 
      dplyr::left_join(
        rela(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "rela"),
        dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(has_slot(object, "norm")) {
    long_omic <- long_omic %>% 
      dplyr::left_join(
        norm(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "norm"),
        dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(has_slot(object, "meta")){
    long_omic <- long_omic %>% 
      dplyr::left_join(
        meta(object, .fmt="tbl"), by = "sample_id")
  }
  
  if(has_slot(object, "taxa")){
    long_omic <- long_omic %>% 
      dplyr::left_join(
        taxa(object, .fmt="tbl"), by = "taxa_id")
  }
  
  return(long_omic)
  
})

setMethod("omic_longer", "omics", function(object) {
  
  # empty
  if( miss_sample(object, "any") || miss_taxa(object, "any") ) cli::cli_abort("Error: No sample or taxa available in at least one of the omic objects.")
  
  purrr::map(object, omic_longer) %>%
    purrr::imap(\(x,y){
      x <- dplyr::mutate(x, omic = y, .before = 1)
    }) %>%
    purrr::list_rbind() %>%
    return()
  
})