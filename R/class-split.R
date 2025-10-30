# SPLIT OMIC SAMPLE
#------------------------------------------------------------------------------#
#' Split an omic Object into an omics Based on Specified Columns
#'
#' @description
#' This method splits an `omic` object into subsets based on unique combinations
#' of values in specified columns from the `meta` data frame. Each subset
#' corresponds to an `omic` object in the resulting `omics`, facilitating
#' analyses specific to each unique combination of metadata values.
#'
#' @param object An `omic` object.
#' @param ... Columns in the `meta` data frame used to define subsets.
#'        You can specify columns by name, using any of the selection methods supported by `dplyr::select()`.
#' @return An `omics` object containing `omic` objects for each unique combination
#'         of specified column values.
#'
#' @details
#' The `split_meta` function leverages `dplyr` functionality to identify unique combinations
#' of specified column values in the `meta` data frame. It then creates a new `omic`
#' object for each combination, allowing for tailored analysis per group.
#'
#' This function is particularly useful in exploratory data analysis and pre-processing stages
#' where data need to be examined or analyzed based on specific grouping variables.
#'
#' @note
#' This function depends on `dplyr` for data manipulation. Ensure that `dplyr` is installed
#' and loaded in your R session.
#'
#' @export
#' @aliases split_meta,omic-method
#' @importFrom dplyr select distinct group_by group_split semi_join
#' @importFrom rlang enquos
setGeneric("split_meta", function(object, ...) {
  standardGeneric("split_meta")
})

setMethod("split_meta", "omic", function(object, ...) {
  
  if(miss_sample(object)) stop("Error: No sample available.")
  
  # Capture the column names as quosures
  split_cols <- enquos(...)
  
  # Check if split_cols is empty
  if (length(split_cols) == 0) {
    stop("No columns specified for splitting. Please specify one or more columns.")
  }
  
  # Get the unique combinations from sample for the specified columns
  meta <- meta(object, "tbl")
  groups <- meta %>%
    dplyr::select(!!!split_cols) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!!split_cols) %>%
    dplyr::group_split()
  
  # Initialize an empty list to store the subsetted omic objects
  omic_list <- omics()
  
  # For each group, filter the sample based on the unique combination and create a new omic object
  for (i in seq_along(groups)) {
    # Filter sample based on the group
    group_df <- groups[[i]]
    filtered_sample <- dplyr::semi_join(meta, group_df, by = names(group_df))
    
    # Filter the original omic object to create a subset based on the filtered sample
    subsetted_omic <- object[filtered_sample$sample_id, , drop = FALSE]
    
    # Add the subsetted omic object to the list with a meaningful name
    group_name <- apply(group_df, 1, function(row) {paste(row, collapse = "-")})
    omic_list[[group_name]] <- subsetted_omic
  }
  
  return(omic_list)
})


# SPLIT OMIC TAXA
#------------------------------------------------------------------------------#
#' Split an omic Object into an omics Based on Specified Columns
#'
#' @description
#' This method splits an `omic` object into subsets based on unique combinations
#' of values in specified columns from the `taxa` data frame. Each subset
#' corresponds to an `omic` object in the resulting `omics`, facilitating
#' analyses specific to each unique combination of metadata values.
#'
#' @param object An `omic` object.
#' @param ... Columns in the `taxa` data frame used to define subsets.
#'        You can specify columns by name, using any of the selection methods supported by `dplyr::select()`.
#' @return An `omics` object containing `omic` objects for each unique combination
#'         of specified column values.
#'
#' @details
#' The `split_taxa` function leverages `dplyr` functionality to identify unique combinations
#' of specified column values in the `taxa` data frame. It then creates a new `omic`
#' object for each combination, allowing for tailored analysis per group.
#'
#' This function is particularly useful in exploratory data analysis and pre-processing stages
#' where data need to be examined or analyzed based on specific grouping variables.
#'
#' @note
#' This function depends on `dplyr` for data manipulation. Ensure that `dplyr` is installed
#' and loaded in your R session.
#'
#' @export
#' @aliases split_taxa,omic-method
#' @importFrom dplyr select distinct group_by group_split semi_join
#' @importFrom rlang enquos
setGeneric("split_taxa", function(object, ...) {
  standardGeneric("split_taxa")
})

setMethod("split_taxa", "omic", function(object, ...) {
  
  if(miss_taxa(object)) stop("Error: No taxa available.")
  
  # Capture the column names as quosures
  split_cols <- enquos(...)
  
  # Check if split_cols is empty
  if (length(split_cols) == 0) {
    stop("No columns specified for splitting. Please specify one or more columns.")
  }
  
  # Get the unique combinations from taxa for the specified columns
  taxa <- taxa(object, "tbl")
  groups <- taxa %>%
    dplyr::select(!!!split_cols) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!!split_cols) %>%
    dplyr::group_split()
  
  # Initialize an empty list to store the subsetted omic objects
  omic_list <- omics()
  
  # For each group, filter the taxa based on the unique combination and create a new omic object
  for (i in seq_along(groups)) {
    # Filter taxa based on the group
    group_df <- groups[[i]]
    filtered_taxa <- dplyr::semi_join(taxa, group_df, by = names(group_df))
    
    # Filter the original omic object to create a subset based on the filtered taxa
    subsetted_omic <- object[, filtered_taxa$taxa_id, drop = FALSE]
    
    # Add the subsetted omic object to the list with a meaningful name
    group_name <- apply(group_df, 1, function(row) {paste(row, collapse = "-")})
    omic_list[[group_name]] <- subsetted_omic
  }
  
  return(omic_list)
})