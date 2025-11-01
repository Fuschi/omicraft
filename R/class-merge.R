#' Merge Multiple `omic` or `omics` Objects
#'
#' @description
#' Merges multiple `omic` or `omics` objects into a single `omic` or `omics` object by merging their respective slots according to user-specified functions. This method provides flexibility to define how each slot is combined, which is essential because a default merging behavior for such complex data structures is not feasible.
#'
#' @param ... `omic` objects, `omics` objects, or lists of such objects to be merged.
#' @param abun A function to merge the 'abun' (abundance) slots from the provided `omic` objects, or `NULL` if no merging is required for this slot.
#' @param rela A function to merge the 'rela' (relative abundance) slots from the provided `omic` objects, or `NULL` if no merging is required for this slot.
#' @param norm A function to merge the 'norm' (normalized abundance) slots from the provided `omic` objects, or `NULL` if no merging is required for this slot.
#' @param meta A function to merge the 'meta' (metadata) slots from the provided `omic` objects, or `NULL` if no merging is required for this slot.
#' @param taxa A function to merge the 'taxa' (taxonomic data) slots from the provided `omic` objects, or `NULL` if no merging is required for this slot.
#' @param netw A function to merge the 'netw' (network) slots from the provided `omic` objects, or `NULL` if no merging is required for this slot.
#' @param comm A function to merge the 'comm' (community structure) slots from the provided `omic` objects, or `NULL` if no merging is required for this slot.
#'
#' @return A new `omic` or `omics` object containing the merged data from the input objects. Each slot of the returned object is the result of the merging function applied to that slot, or remains `NULL` if no function was specified for that slot.
#'
#' @details
#' The function allows specifying custom merging strategies for each slot:
#' - **Abundance data (`abun`)**: Functions like `dplyr::bind_cols` or `dplyr::bind_rows` might be appropriate.
#' - **Metadata (`meta`)**: Merging can be handled by functions such as `dplyr::full_join` or `dplyr::bind_rows`.
#' - **Network structures (`netw`)**: These might be combined using `igraph::union` or similar graph-specific functions.
#' This approach ensures flexibility in handling the complex and varied data stored within `omic` objects.
#'
#' The method checks for consistency among `omics` objects. If merging `omics` objects, it first ensures that all have identical names in the 'omics' slot. Mismatched names will cause the function to stop and issue an error.
#'
#' @importFrom methods new
#' @importFrom purrr map reduce map_lgl list_flatten
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#'   # Assuming mg1 and mg2 are omic objects
#'   merged_omic <- merge_omics(mg1, mg2, abun = rbind)
#'   # For merging omics objects with custom functions for each slot
#'   merged_list <- merge_omics(list1, list2, 
#'                               meta = function(x) dplyr::bind_rows(x),
#'                               netw = igraph::graph.union)
#' }
merge_omics <- function(..., 
                         abun = NULL, rela = NULL, norm = NULL,
                         meta = NULL, taxa = NULL,
                         netw = NULL, comm = NULL) {
  
  # Flatten input to handle lists of omic objects or omics
  omics <- purrr::list_flatten(list(...))
  
  if(length(omics) <= 1) return(omics)
  
  # Check for non-omic objects and stop execution if found
  are_omic     <- all(purrr::map_lgl(omics, \(x) inherits(x, "omic")))
  are_omics <- all(purrr::map_lgl(omics, \(x) inherits(x, "omics")))
  if (isFALSE(are_omic) & isFALSE(are_omics)) {
    stop("Error: All elements must be 'omic' or 'omics'.")
  }
  
  # Define a helper function to check if all items are equal
  all_equal <- function(items) {
    all(sapply(2:length(items), function(i) identical(items[[1]], items[[i]])))
  }
  
  # MGNET
  #--------------------------------------------#
  if(are_omic){
    new_omic <- new("omic")
    
    # Merge each slot using the provided functions or leave them as NULL
    if(!is.null(abun)){
      new_omic@abun <- purrr::map(omics, \(x) omicraft::abun(x, .fmt = "df")) %>%
        purrr::reduce(abun) %>%
        as.matrix()
    } else if(all_equal(purrr::map(omics, \(x) omicraft::abun(x)))){
      new_omic@abun <- omicraft::abun(omics[[1]])
    }
    
    if(!is.null(rela)){
      new_omic@rela <- purrr::map(omics, \(x) omicraft::rela(x, .fmt = "df")) %>%
        purrr::reduce(rela) %>%
        as.matrix()
    } else if(all_equal(purrr::map(omics, \(x) omicraft::rela(x)))){
      new_omic@rela <- omicraft::rela(omics[[1]])
    }
    
    if(!is.null(norm)){
      new_omic@norm <- purrr::map(omics, \(x) omicraft::norm(x, .fmt = "df")) %>%
        purrr::reduce(norm) %>%
        as.matrix()
    } else if(all_equal(purrr::map(omics, \(x) omicraft::norm(x)))){
      new_omic@norm <- omicraft::norm(omics[[1]])
    }
    
    if(!is.null(meta)){
      new_omic@meta <- purrr::map(omics, \(x) omicraft::meta(x)) %>%
        purrr::reduce(meta)
    } else if(all_equal(purrr::map(omics, \(x) omicraft::meta(x)))){
      new_omic@meta <- omicraft::meta(omics[[1]])
    }
    
    if(!is.null(taxa)){
      new_omic@taxa <- purrr::map(omics, \(x) omicraft::taxa(x)) %>%
        purrr::reduce(taxa)
    } else if(all_equal(purrr::map(omics, \(x) omicraft::taxa(x)))){
      new_omic@taxa <- omicraft::taxa(omics[[1]])
    }
    
    if(!is.null(netw)){
      new_omic@netw <- purrr::map(omics, \(x) omicraft::netw(x)) %>%
        purrr::reduce(netw)
    } else if(all_equal(purrr::map(omics, \(x) omicraft::netw(x)))){
      new_omic@netw <- omicraft::netw(omics[[1]])
    }
    
    if(!is.null(comm)){
      new_omic@comm <- purrr::map(omics, \(x) omicraft::comm(x)) %>%
        purrr::reduce(comm)
    } else if(all_equal(purrr::map(omics, \(x) omicraft::comm(x)))){
      new_omic@comm <- omicraft::comm(omics[[1]])
    }
    
    validObject(new_omic)
    return(new_omic)
  }
  
  
  # MGNETLIST
  #--------------------------------------------#
  if(are_omics){
    
    names_list <- purrr::map(omics, names)
    
    # Check if all lists of names are equal
    are_names_equal <- purrr::reduce(names_list, \(a,b) all(a%in%b) & all(b%in%a))
    
    # If names are not equal, stop the execution and report the inconsistency
    if (!are_names_equal) {
      stop("Error: The omics objects have different names in the 'omics' slot.")
    }
    
    unique_names <- unique(unlist(names_list))
    merged_list <- purrr::map(unique_names, \(name) {
      
      omics_name <- purrr::map(omics, \(x){x[[name]]})
      return(merge_omics(omics_name,
                          abun = abun, rela = rela, norm = norm,
                          meta = meta, taxa = taxa,
                          netw = netw, comm = comm))
    }) %>% 
      stats::setNames(unique_names) %>% 
      omics()
    return(merged_list)
  }
  
}


#' Collapse Specified Combinations within an omics
#'
#' @description
#' Merges specific combinations of `omic` objects within an `omics` according to user-defined rules,
#' allowing for targeted merging of subsets within the list. This function facilitates focused analyses
#' by enabling the aggregation of related data into meaningful groups.
#'
#' @param object An `omics` containing multiple `omic` objects.
#' @param by A list specifying the combinations of `omic` object names to be merged. Each element in
#'           the list should be a character vector containing the names of the `omic` objects to merge.
#'           If NULL, all `omic` objects in the list are merged into a single `omic` object.
#' @param abun A function to merge the 'abun' (abundance) slots from the selected `omic` objects, or
#'             NULL if no merging is required for this slot.
#' @param rela A function to merge the 'rela' (relative abundance) slots, or NULL if not required.
#' @param norm A function to merge the 'norm' (normalized abundance) slots, or NULL if not required.
#' @param meta A function to merge the 'meta' (metadata) slots, or NULL if not required.
#' @param taxa A function to merge the 'taxa' (taxonomic data) slots, or NULL if not required.
#' @param netw A function to merge the 'netw' (network) slots, or NULL if not required.
#' @param comm A function to merge the 'comm' (community structure) slots, or NULL if not required.
#'
#' @return An `omics` containing the merged `omic` objects. Each entry in the returned `omics`
#'         corresponds to a merged group as defined in the `by` parameter. The names of the entries
#'         reflect the groups specified in `by` or are derived by concatenating the names of the merged
#'         `omic` objects with dashes if `by` is unnamed.
#'
#' @details
#' The function checks for the validity of names specified in the `by` parameter against the names
#' of `omic` objects in the `omics`. It ensures that all specified groups contain valid `omic`
#' object names and that each group's name is unique. This function is useful for scenarios where
#' specific subsets of data within an `omics` need to be aggregated based on certain criteria
#' or experimental conditions.
#'
#' If `by` is NULL, a default merge of all `omic` objects in the list is performed, resulting in a
#' single `omic` object.
#'
#' Example usage:
#' \dontrun{
#'   # Assuming list1 contains several omic objects named 'A', 'B', 'C', etc.
#'   result <- omic_collapse(list1, by = list(Group1 = c("A", "B"), Group2 = c("C")),
#'                            meta = function(x) dplyr::bind_rows(x),
#'                            netw = igraph::graph.union)
#' }
#'
#' @importFrom stats setNames
#' @name omic_collapse
#' @aliases omic_collapse,omics-method
#' @export
setGeneric("omic_collapse", function(object, by = NULL,
                                      abun = NULL, rela = NULL, norm = NULL, 
                                      meta = NULL, taxa = NULL, 
                                      netw = NULL, comm = NULL) {
  standardGeneric("omic_collapse")
})

setMethod("omic_collapse", "omics", function(object, by = NULL,
                                                  abun = NULL, rela = NULL, norm = NULL, 
                                                  meta = NULL, taxa = NULL, 
                                                  netw = NULL, comm = NULL) {
  
  if(is.null(by)){
    return(merge_omics(as.list(object), 
                        abun = abun, rela = rela, norm = norm,
                        meta = meta, taxa = taxa, 
                        netw = netw, comm = comm))
  } 
  
  if (!is.list(by)) {
    stop("'by' must be a list of character vectors specifying combinations of omic names.")
  }
  
  # Check the validity of elements in 'by'
  if (!all(sapply(by, \(x) is.character(x) && all(x %in% names(object))))) {
    stop("All elements in 'by' must be character vectors containing valid names of omic objects within the omics.")
  }
  
  # Prepare names for the output based on 'by' or create descriptive names
  names(by) <- ifelse(nzchar(names(by)), names(by), sapply(by, \(x) paste(x, collapse = "-")))
  
  # Merge specified omic combinations
  results <- stats::setNames(vector("list", length(by)), names(by))
  for (i in seq_along(by)) {
    sub_object <- omics(object)[by[[i]]]
    results[[i]] <- merge_omics(sub_object, 
                                 abun = abun, rela = rela, norm = norm,
                                 meta = meta, taxa = taxa, 
                                 netw = netw, comm = comm)
  }
  
  return(omics(results))
})
