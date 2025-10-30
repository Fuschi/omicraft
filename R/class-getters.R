#'@include class-omic.R class-omics.R
NULL

# ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Abundance Data 
#'
#' @description
#' Retrieves and aggregates abundance data for an `omic` or `omics` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain abundances at different taxonomic ranks or other categorical variables.
#'
#' @param object An `omic` or `omics` object containing abundance and informational data.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `omic` objects, the row names of the abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original abundance data without any aggregation.
#' @param .fun A function to specify how the abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `omic` object, returns abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `omics` object, returns a list of such formatted data from each `omic`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name abun
#' @aliases abun,omic-method abun,omics-method
setGeneric("abun", function(object, .fmt = "mat", .var = NULL, .fun = sum) standardGeneric("abun"))

setMethod("abun", "omic", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  # Checks
  if(!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if(!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handling empty abundance data
  if(length(object@abun) == 0) {
    switch(.fmt,
           mat = matrix(nrow = 0, ncol = 0),
           df = data.frame(),
           tbl = tibble::tibble())
  } else if(is.null(.var)) {
    # Return data in specified format without aggregation
    switch(.fmt,
           mat = object@abun,
           df = as.data.frame(object@abun),
           tbl = tibble::as_tibble(object@abun, rownames = "sample_id"))
  } else {
    # Ensure .var is a valid column name in taxa
    if(!(.var %in% taxa_vars(object))) {
      stop(".var must be present in the taxa columns. Available choices are: ", 
           paste(colnames(object@taxa), collapse = ", "))
    }
    
    # Prepare data for aggregation
    data_frame <- object@abun %>% 
      tibble::as_tibble(rownames = "sample_id") %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "abun")
    
    taxonomic_info <- taxa(object, .fmt = "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonomic_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
      dplyr::summarise(abun = .fun(abun), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "abun") %>%
      dplyr::arrange(match(sample_id, sample_id(object)))
    
    result_matrix <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
    
    switch(.fmt,
           mat = {return(result_matrix)},
           df  = {return(as.data.frame(result_matrix))},
           tbl = {return(tibble::as_tibble(result_matrix, rownames="sample_id"))}
    )
  }
})

setMethod("abun", "omics", function(object,  .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object@omics, function(x) abun(object = x, .var = .var, 
                                                   .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# RELATIVE ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Relative Abundance Data Based on Taxonomic Information
#'
#' @description
#' Retrieves and aggregates relative abundance data for an `omic` or `omics` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original relative abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain relative abundances at different taxonomic ranks or other categorical variables.
#'
#' @param object An `omic` or `omics` object containing relative abundance and informational data.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `omic` objects, the row names of the relative abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original relative abundance data without any aggregation.
#' @param .fun A function to specify how the relative abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `omic` object, returns relative abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `omics` object, returns a list of such formatted data from each `omic`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name rela
#' @aliases rela,omic-method rela,omics-method
setGeneric("rela", function(object, .fmt = "mat", .var = NULL, .fun = sum) standardGeneric("rela"))

setMethod("rela", "omic", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  # Checks
  if(!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if(!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handling empty rela data
  if(length(object@rela) == 0) {
    switch(.fmt,
           mat = matrix(nrow = 0, ncol = 0),
           df = data.frame(),
           tbl = tibble::tibble())
  } else if(is.null(.var)) {
    # Return data in specified format without aggregation
    switch(.fmt,
           mat = object@rela,
           df = as.data.frame(object@rela),
           tbl = tibble::as_tibble(object@rela, rownames = "sample_id"))
  } else {
    # Ensure .var is a valid column name in taxa
    if(!(.var %in% colnames(object@taxa))) {
      stop(".var must be present in the taxa columns. Available choices are: ", 
           paste(colnames(object@taxa), collapse = ", "))
    }
    
    # Prepare data for aggregation
    data_frame <- object@rela %>% 
      tibble::as_tibble(rownames = "sample_id") %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "rela")
    
    taxonomic_info <- taxa(object, "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonomic_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
      dplyr::summarise(rela = .fun(rela), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "rela") %>%
      dplyr::arrange(match(sample_id, sample_id(object)))
    
    result_matrix <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
    
    switch(.fmt,
           mat = {return(result_matrix)},
           df  = {return(as.data.frame(result_matrix))},
           tbl = {return(tibble::as_tibble(result_matrix, rownames="sample_id"))}
    )
  }
})

setMethod("rela", "omics", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object@omics, function(x) rela(object = x, .var = .var, 
                                                   .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# NORMALIZED ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Normalized Abundance Data Based on Taxonomic Information
#'
#' @description
#' Retrieves and aggregates normalized abundance data for an `omic` or `omics` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original normalized abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain normalized abundances at different taxonomic ranks or other categorical variables.
#'
#' @param object An `omic` or `omics` object containing normalized abundance and informational data.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `omic` objects, the row names of the normalized abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original normalized abundance data without any aggregation.
#' @param .fun A function to specify how the relative abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `omic` object, returns normalized abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `omics` object, returns a list of such formatted data from each `omic`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name norm
#' @aliases norm,omic-method norm,omics-method
setGeneric("norm", function(object, .fmt = "mat", .var = NULL, .fun = sum) standardGeneric("norm"))

setMethod("norm", "omic", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  # Checks
  if(!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if(!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handling empty norm data
  if(length(object@norm) == 0) {
    switch(.fmt,
           mat = matrix(nrow = 0, ncol = 0),
           df = data.frame(),
           tbl = tibble::tibble())
  } else if(is.null(.var)) {
    # Return data in specified format without aggregation
    switch(.fmt,
           mat = object@norm,
           df = as.data.frame(object@norm),
           tbl = tibble::as_tibble(object@norm, rownames = "sample_id"))
  } else {
    # Ensure .var is a valid column name in taxa
    if(!(.var %in% colnames(object@taxa))) {
      stop(".var must be present in the taxa columns. Available choices are: ", 
           paste(colnames(object@taxa), collapse = ", "))
    }
    
    # Prepare data for aggregation
    data_frame <- object@norm %>% 
      tibble::as_tibble(rownames = "sample_id") %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "norm")
    
    taxonomic_info <- taxa(object, "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonomic_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
      dplyr::summarise(norm = .fun(norm), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "norm") %>%
      dplyr::arrange(match(sample_id, sample_id(object)))
    
    result_matrix <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
    
    switch(.fmt,
           mat = {return(result_matrix)},
           df  = {return(as.data.frame(result_matrix))},
           tbl = {return(tibble::as_tibble(result_matrix, rownames="sample_id"))}
    )
  }
})

setMethod("norm", "omics", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object, function(x) norm(object = x, .var = .var, .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# META
#------------------------------------------------------------------------------#
#' Get Sample Metadata Information
#'
#' Retrieves the sample information stored in the `meta` slot of an `omic` object
#' or for each `omic` object within an `omics`, with the option to format the output as
#' a `data.frame`, `tibble`.
#'
#' @param object An `omic` or `omics` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `omic` objects, the row names of
#'          the abundance matrix are converted into a new column named `sample_id`, ensuring alignment
#'          with the reserved keyword in `omic-class`.
#' @param .collapse Logical, only for `omics`. If TRUE, ignore `.fmt` and return
#'   a single tibble with all metas row-bound together and an `omic` column  
#'   indicating the source object. Default: FALSE.
#'          
#' @return The content of the `meta` slot for `omic` or a list of such contents for `omics`.
#' 
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @export
#' @name meta
#' @aliases meta,omic-method meta,omics-method
setGeneric("meta", function(object, .fmt = "df", .collapse = FALSE) standardGeneric("meta"))

setMethod("meta", "omic", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if( length(object@meta) == 0 ){
    
    switch(.fmt,
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
    
  } else {
    
    switch(.fmt,
           df  = {return(object@meta)},
           tbl = {return(tibble::as_tibble(object@meta, rownames="sample_id"))})
    
  }
})

setMethod("meta", "omics", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if (.collapse) {
    meta_collapsed <- object %>% 
      purrr::map(\(x) meta(x, .fmt = "tbl")) %>% 
      purrr::imap(\(x, y) tibble::add_column(x, omic = y, .before = 1)) %>% 
      purrr::list_rbind()
      return(meta_collapsed)
  }
  
  if(.fmt == "df") {
    
    if(length(object) == 0) return(data.frame())
    return(purrr::map(object, function(x) meta(x, .fmt = "df")))
    
  } else if(.fmt == "tbl") {
    
    if(length(object) == 0) return(tibble::tibble())
    return(purrr::map(object, function(x) meta(x, .fmt = "tbl")))
    
  }
  
})


# TAXA
#------------------------------------------------------------------------------#
#' Get Taxa Metadata Information
#'
#' Retrieves the taxa information stored in the `taxa` slot of an `omic` object
#' or for each `omic` object within an `omics`, with the option to format the output as
#' a `data.frame`, `tibble`, or a combined `tibble` for multiple omic objects.
#'
#' @param object An `omic` or `omics` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `omic` objects, the the row names of
#'          the abundance matrix are converted into a new column named `taxa_id`, ensuring alignment
#'          with the reserved keyword in `omic-class`.
#' @param .collapse Logical, only for `omics`. If TRUE, ignore `.fmt` and return
#'   a single tibble with all taxas row-bound together and an `omic` column  
#'   indicating the source object. Default: FALSE.
#'          
#' @return The content of the `taxa` slot for `omic` or a list of such contents for `omics`.
#' 
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @export
#' @name taxa
#' @aliases taxa,omic-method taxa,omics-method
setGeneric("taxa", function(object, .fmt = "df", .collapse = FALSE) standardGeneric("taxa"))

setMethod("taxa", "omic", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if(length(object@taxa) == 0 && length(object@comm) == 0){
    switch(.fmt,
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
  }
  
  if(length(object@taxa) != 0 && length(object@comm) == 0){
    switch(.fmt,
           df  = {return(object@taxa)},
           tbl = {return(tibble::as_tibble(object@taxa, rownames="taxa_id"))})
  }
  
  if(length(object@taxa) == 0 && length(object@comm) != 0){
    switch(.fmt,
           df  = {return(data.frame(taxa_id = taxa_id(object),
                                    comm_id = comm_id(object)))},
           tbl = {return(comm_id(object, .fmt = "tbl"))})
  }
  
  if(length(object@taxa) != 0 && length(object@comm) != 0){
    
    result <- comm_id(object, .fmt = "tbl") %>%
      left_join(tibble::rownames_to_column(object@taxa, "taxa_id"),  
                by = "taxa_id")
    
    switch(.fmt,
           df  = {return(tibble::column_to_rownames(result, "taxa_id"))},
           tbl = {return(result)})
    
  }
  
  
})

setMethod("taxa", "omics", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if (.collapse) {
    taxa_collapsed <- object %>% 
      purrr::map(\(x) taxa(x, .fmt = "tbl")) %>% 
      purrr::imap(\(x, y) tibble::add_column(x, omic = y, .before = 1)) %>% 
      purrr::list_rbind()
    return(taxa_collapsed)
  }
  
  if(.fmt == "df") {
    
    if(length(object) == 0) return(data.frame())
    return(sapply(object, taxa, .fmt = "df", simplify = FALSE, USE.NAMES = TRUE))
    
  } else if(.fmt == "tbl") {
    
    if(length(object) == 0) return(tibble::tibble())
    return(sapply(object, taxa, .fmt = "tbl", simplify = FALSE, USE.NAMES = TRUE))
    
  } 
  
})


#' @title Retrieve Network Graph from omic Objects
#'
#' @description
#' Retrieves the network graph stored in the \code{netw} slot of an \code{omic} object,
#' or from each \code{omic} object within an \code{omics}. The network graph is
#' an \code{igraph} object representing interactions among different taxa.
#'
#' If the object has a non-empty \code{taxa} slot, all taxa information is automatically
#' attached as vertex attributes in the returned graph, making node-level metadata 
#' immediately accessible for further analysis or visualization.
#'
#' @param object An \code{omic} or \code{omics} object.
#'
#' @return 
#' \itemize{
#'   \item For an \code{omic}, a single \code{igraph} object is returned, potentially 
#'         with vertex attributes derived from the \code{taxa} slot.
#'   \item For an \code{omics}, a \code{list} of \code{igraph} objects is returned
#'         (one per contained \code{omic}), each enriched with node-level attributes
#'         if available.
#' }
#'
#' @importFrom igraph vertex_attr set_edge_attr is_weighted E membership
#' @aliases netw,omic-method netw,omics-method
#' @export
setGeneric("netw", function(object) standardGeneric("netw"))

setMethod("netw", "omic", function(object) {
  
  g <- object@netw
  if (length(g) == 0) return(g)
  
  # If taxa slot is non-empty, set each column as a vertex attribute
  if (has_metataxa(object)) {
    metataxa <- taxa(object)
    for (vertex_attr in colnames(metataxa)) {
      igraph::vertex_attr(g, vertex_attr) <- object@taxa[[vertex_attr]]
    }
  }
  
  # If the links are selected get the network with the filtered links
  if(are_selected_links(object)){
    g <- igraph::subgraph_from_edges(graph = g,
                                     eids = get_selected_links(object), 
                                     delete.vertices = FALSE)
  }
  
  g
})

setMethod("netw", "omics", function(object) {
  sapply(object@omics, function(x) netw(x), simplify = FALSE, USE.NAMES = TRUE)
})


# COMMUNITY
#------------------------------------------------------------------------------#
#' Get Community Detection Results
#'
#' Retrieves the community detection results stored in the `comm` slot of an `omic` object
#' or each `omic` object within an `omics`. These results are typically derived from network
#' analysis methods and indicate the grouping or clustering of taxa into communities based on their
#' network interactions.
#'
#' @param object An `omic` or `omics` object.
#' @return The community detection result object stored in the `comm` slot for an `omic` object,
#'         or a list of such objects for an `omics` object, each representing the community
#'         detection results of a contained `omic` object.
#' @export
#' @name comm
#' @aliases comm,omic-method comm,omics-method
setGeneric("comm", function(object) standardGeneric("comm"))

setMethod("comm", "omic", function(object) {
  object@comm
})

setMethod("comm", "omics", function(object) {
  sapply(object@omics, function(x) x@comm, simplify = FALSE, USE.NAMES = TRUE)
})