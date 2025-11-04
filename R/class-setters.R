#'@include class-omic.R class-omics.R
NULL

# ABUNDANCE
#------------------------------------------------------------------------------#
#' Set Abundance Data
#'
#' This setter function allows you to update the abundance data for `omic` objects
#' and each `omic` object within a `omics`. The abundance data must be a numeric matrix
#' for `omic` objects. For `omics` objects, the abundance data should be a named list
#' of numeric matrices corresponding to each `omic` object within the list.
#'
#' @param object An `omic` or `omics` object.
#' @param value The new abundance data to be set.
#' @return The `omic` or `omics` object with updated abundance data.
#' @export
#' @importFrom methods validObject
#' @name abun<-
#' @aliases abun<-,omic-method abun<-,omics-method
setGeneric("abun<-", function(object, value) standardGeneric("abun<-"))

setMethod("abun<-", c("omic","ANY"), function(object, value){
  
  if (!(inherits(value, "matrix") || inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg value} must be a tibble or data frame.")
  }
  
  if ("sample_id" %in% colnames(value) && !is.null(rownames(value))) {
    cli::cli_abort("{.arg value} cannot have both a {.val sample_id} column and row names.")
  }
  
  if (!"sample_id" %in% colnames(value) && is.null(rownames(value))) {
    cli::cli_abort("{.arg value} must have either a {.val sample_id} column or row names.")
  }
  
  if ("sample_id" %in% colnames(value)) {
    if (any(duplicated(value$sample_id))) {
      cli::cli_abort("Column {.var sample_id} contains duplicated values in {.arg value}.")
    }
    rownames(value) <- value[,"sample_id"]
    value <- value[,colnames(value) != "sample_id", drop = F]
  }
  
  if (has_sample(object)) {
    if (!all(rownames(value) %in% sample_id(object))) {
      cli::cli_abort("Provided {.var sample_id} values do not match those in the object.")
    }
    if (!identical(rownames(value), sample_id(object))) {
      value <- value[sample_id(object), , drop = FALSE]
    }
  }
  
  object@abun <- as.matrix(value)
  validObject(object)
  object
})

setMethod("abun<-", c("omics","ANY"), function(object, value){
  is_list_omics_assign(object, value)
  for(i in names(object)) { abun(object@omics[[i]]) <- value[[i]] }
  validObject(object)
  object
})


# rela
#------------------------------------------------------------------------------#
#' Set rela Data
#'
#' This setter function allows you to update the rela data for `omic` objects
#' and each `omic` object within a `omics`. The rela data must be a numeric matrix
#' for `omic` objects. For `omics` objects, the rela data should be a named list
#' of numeric matrices corresponding to each `omic` object within the list.
#'
#' @param object An `omic` or `omics` object.
#' @param value The new rela data to be set.
#' @return The `omic` or `omics` object with updated abundance data.
#' @export
#' @importFrom methods validObject
#' @name rela<-
#' @aliases rela<-,omic-method rela<-,omics-method
setGeneric("rela<-", function(object, value) standardGeneric("rela<-"))

setMethod("rela<-", c("omic","ANY"), function(object, value){
  
  if (!(inherits(value, "matrix") || inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg value} must be a tibble or data frame.")
  }
  
  if ("sample_id" %in% colnames(value) && !is.null(rownames(value))) {
    cli::cli_abort("{.arg value} cannot have both a {.val sample_id} column and row names.")
  }
  
  if (!"sample_id" %in% colnames(value) && is.null(rownames(value))) {
    cli::cli_abort("{.arg value} must have either a {.val sample_id} column or row names.")
  }
  
  if ("sample_id" %in% colnames(value)) {
    if (any(duplicated(value$sample_id))) {
      cli::cli_abort("Column {.var sample_id} contains duplicated values in {.arg value}.")
    }
    rownames(value) <- value[,"sample_id"]
    value <- value[,colnames(value) != "sample_id", drop = F]
  }
  
  if (has_sample(object)) {
    if (!all(rownames(value) %in% sample_id(object))) {
      cli::cli_abort("Provided {.var sample_id} values do not match those in the object.")
    }
    if (!identical(rownames(value), sample_id(object))) {
      value <- value[sample_id(object), , drop = FALSE]
    }
  }
  
  object@rela <- as.matrix(value)
  validObject(object)
  object
})

setMethod("rela<-", c("omics","ANY"), function(object, value){
  is_list_omics_assign(object, value)
  for(i in names(object)) { rela(object@omics[[i]]) <- value[[i]] }
  validObject(object)
  object
})


# norm
#------------------------------------------------------------------------------#
#' Set Log-Transformed Abundance Data
#'
#' This setter function allows you to update the log-transformed abundance data for `omic` objects
#' and each `omic` object within a `omics`. The log-transformed abundance data must be a numeric matrix
#' for `omic` objects. For `omics` objects, the data should be a named list of numeric matrices
#' corresponding to each `omic` object within the list.
#'
#' @param object An `omic` or `omics` object.
#' @param value The new log-transformed abundance data to be set.
#' @return The `omic` or `omics` object with updated log-transformed abundance data.
#' @export
#' @importFrom methods validObject
#' @name norm<-
#' @aliases norm<-,omic-method norm<-,omics-method
setGeneric("norm<-", function(object, value) standardGeneric("norm<-"))

setMethod("norm<-", "omic", function(object, value) {
  
  if (!(inherits(value, "matrix") || inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg value} must be a tibble or data frame.")
  }
  
  if ("sample_id" %in% colnames(value) && !is.null(rownames(value))) {
    cli::cli_abort("{.arg value} cannot have both a {.val sample_id} column and row names.")
  }
  
  if (!"sample_id" %in% colnames(value) && is.null(rownames(value))) {
    cli::cli_abort("{.arg value} must have either a {.val sample_id} column or row names.")
  }
  
  if ("sample_id" %in% colnames(value)) {
    if (any(duplicated(value$sample_id))) {
      cli::cli_abort("Column {.var sample_id} contains duplicated values in {.arg value}.")
    }
    rownames(value) <- value[,"sample_id"]
    value <- value[,colnames(value) != "sample_id", drop = F]
  }
  
  if (has_sample(object)) {
    if (!all(rownames(value) %in% sample_id(object))) {
      cli::cli_abort("Provided {.var sample_id} values do not match those in the object.")
    }
    if (!identical(rownames(value), sample_id(object))) {
      value <- value[sample_id(object), , drop = FALSE]
    }
  }
  
  object@norm <- as.matrix(value)
  validObject(object)
  object
})

setMethod("norm<-", "omics", function(object, value) {
  is_list_omics_assign(object, value)
  for(i in names(object)) { norm(object@omics[[i]]) <- value[[i]] }
  validObject(object)
  object
})



# META
#------------------------------------------------------------------------------#
#' Update Sample Metadata for omic and omics Objects
#'
#' This function sets the sample metadata for `omic` objects and also updates
#' each `omic` object within a `omics`. The metadata for `omic` objects must be
#' provided as a dataframe. For `omics` objects, the metadata should be either
#' a named list of dataframes corresponding to each `omic` object within the list,
#' or a single dataframe containing a 'omic' and 'sample_id' columns to split the 
#' metadata accordingly.
#'
#' @param object An `omic` or `omics` object to be updated.
#' @param value The new sample metadata to set, either a dataframe or a named list
#' of dataframes as per the object type.
#' @return The updated `omic` or `omics` object with new sample metadata.
#' @export
#' @importFrom methods validObject
#' @importFrom tibble column_to_rownames has_rownames
#' @name meta<-
#' @aliases meta<-,omic-method meta<-,omics-method
setGeneric("meta<-", function(object, value) standardGeneric("meta<-"))

setMethod("meta<-", c("omic", "ANY"), function(object, value) {
  
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  
  if (!(inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg {valueName}} must be a tibble or data frame.")
  }
  
  if ("sample_id" %in% colnames(value) && tibble::has_rownames(value)) {
    cli::cli_abort("{.arg {valueName}} cannot have both a {.val sample_id} column and row names.")
  }
  
  if (!"sample_id" %in% colnames(value) && !tibble::has_rownames(value)) {
    cli::cli_abort("{.arg {valueName}} must have either a {.val sample_id} column or row names.")
  }
  
  if ("sample_id" %in% colnames(value)) {
    if (any(duplicated(value$sample_id))) {
      cli::cli_abort("Column {.var sample_id} contains duplicated values in {.arg value}.")
    }
    value <- tibble::column_to_rownames(value, "sample_id")
  }
  
  if (has_sample(object)) {
    if (!all(rownames(value) %in% sample_id(object))) {
      cli::cli_abort("Provided {.var sample_id} values do not match those in the object.")
    }
    if (!identical(rownames(value), sample_id(object))) {
      value <- value[sample_id(object), , drop = FALSE]
    }
  }
  
  object@meta <- as.data.frame(value)
  validObject(object)
  object
})

setMethod("meta<-", c("omics","ANY"), function(object, value){
  
  if(class(value)[[1]] == "list"){
    
    is_list_omics_assign(object, value)
    for(i in names(object)) meta(object[[i]]) <- value[[i]] 
    
  } else if(is.data.frame(value)){
    
    is_assign_omics_tbl(object, value, "sample")
    splitted_value <- split(value, value$omic)
    splitted_value <- lapply(splitted_value, \(x){
      x$omic <- NULL
      return(x)
    })
    for(i in names(object)) meta(object[[i]]) <- splitted_value[[i]]
    
  } else {
    
    valueName <- deparse(substitute(value))
    cli::cli_abort(
      "{.arg value} must be either a named list of data frames or a single data frame with columns {.val omic} and {.val sample_id}.")
    
  }
  
  validObject(object)
  object
})


# TAXA
#------------------------------------------------------------------------------#
#' Update Taxa Metadata for omic and omics Objects
#'
#' This function sets the taxa metadata for `omic` objects and also updates
#' each `omic` object within a `omics`. The metadata for `omic` objects must be
#' provided as a dataframe. For `omics` objects, the metadata should be either
#' a named list of dataframes corresponding to each `omic` object within the list,
#' or a single dataframe containing a 'omic' and 'taxa_id' columns to split the 
#' metadata accordingly.
#'
#' @param object An `omic` or `omics` object to be updated.
#' @param value The new taxa metadata to set, either a dataframe or a named list
#' of dataframes as per the object type.
#' @return The updated `omic` or `omics` object with new taxa metadata.
#' @export
#' @importFrom methods validObject
#' @importFrom tibble column_to_rownames
#' @name taxa<-
#' @aliases taxa<-,omic-method taxa<-,omics-method
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))

setMethod("taxa<-", c("omic", "ANY"), function(object, value) {
  
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  
  if (!(inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg {valueName}} must be a tibble or data frame.")
  }
  
  if ("taxa_id" %in% colnames(value) && tibble::has_rownames(value)) {
    cli::cli_abort("{.arg {valueName}} cannot have both a {.val taxa_id} column and row names.")
  }
  
  if (!"taxa_id" %in% colnames(value) && !tibble::has_rownames(value)) {
    cli::cli_abort("{.arg {valueName}} must have either a {.val taxa_id} column or row names.")
  }
  
  if ("taxa_id" %in% colnames(value)) {
    if (any(duplicated(value$taxa_id))) {
      cli::cli_abort("Column {.var taxa_id} contains duplicated values in {.arg value}.")
    }
    value <- tibble::column_to_rownames(value, "taxa_id")
  }
  
  if (has_sample(object)) {  # (se qui intendevi has_taxa(object), lascio nota sotto)
    if (!all(rownames(value) %in% taxa_id(object))) {
      cli::cli_abort("Provided {.var taxa_id} values do not match those in the object.")
    }
    if (!identical(rownames(value), taxa_id(object))) {
      value <- value[taxa_id(object), , drop = FALSE]
    }
  }
  
  if ("comm_id" %in% colnames(value)) {
    if (!all(value$comm_id == comm_id(object))) {
      cli::cli_abort("Community memberships in column {.var comm_id} differ from those stored in the object.")
    }
    value$comm_id <- NULL
  }
  
  object@taxa <- as.data.frame(value)
  validObject(object)
  object
})

setMethod("taxa<-", c("omics","ANY"), function(object, value){
  
  if(class(value)[[1]] == "list"){
    
    is_list_omics_assign(object, value)
    for(i in names(object)) taxa(object[[i]]) <- value[[i]] 
    
  } else if(is.data.frame(value)){
    
    is_assign_omics_tbl(object, value, "taxa")
    splitted_value <- split(value, value$omic)
    splitted_value <- lapply(splitted_value, \(x){
      x$omic <- NULL
      return(x)
    })
    for(i in names(object)) taxa(object[[i]]) <- splitted_value[[i]]
    
  } else {
    
    valueName <- deparse(substitute(value))
    cli::cli_abort("{.arg value} must be either a named list of data frames or a single data frame with columns {.val omic} and {.val taxa_id}.")

  }
  
  validObject(object)
  object
})


# NETWORK
#------------------------------------------------------------------------------#
#' Set Network Data
#'
#' This setter function updates the network data for `omic` objects and
#' each `omic` object within a `omics`. The network data must be an `igraph` object
#' for `omic` objects. For `omics` objects, the network data should be a named list of `igraph` objects
#' corresponding to each `omic` object within the list.
#'
#' @param object An `omic` or `omics` object.
#' @param value The new network data to be set.
#' @return The `omic` or `omics` object with updated network data.
#' @export
#' @importFrom methods validObject
#' @name netw<-
#' @aliases netw<-,omic-method netw<-,omics-method
setGeneric("netw<-", function(object, value) standardGeneric("netw<-"))

setMethod("netw<-", "omic", function(object, value) {
  
  # --- Ensure link_id exists if edges are present ----------------------------#
  if (igraph::ecount(value) > 0L && is.null(igraph::edge_attr(value, "link_id"))) {
    if (igraph::ecount(value) > 0L && is.null(igraph::edge_attr(value, "link_id"))) {
      igraph::E(value)$link_id <- paste(
        igraph::as_data_frame(value)[["from"]],
        igraph::as_data_frame(value)[["to"]],
        seq_len(igraph::ecount(value)),
        sep = "--")
    }
  }
  
  object@netw <- value
  validObject(object)
  object
})

setMethod("netw<-", "omics", function(object, value) {
  is_list_omics_assign(object, value)
  for(i in names(object)) { netw(object@omics[[i]]) <- value[[i]] }
  validObject(object)
  object
})

# LINK (SET)
#------------------------------------------------------------------------------#
#' Update Network Links (Edges) for omic and omics Objects
#'
#' This function sets the network links (edges) for `omic` objects and also
#' updates each `omic` inside an `omics`. For `omic`, the value must be a
#' data frame / tibble with at least columns `from` and `to`. For `omics`,
#' the value can be either:
#' - a **named list** of data frames (one per contained `omic`, names matching `names(object)`), or
#' - a **single** data frame with a routing column `omic` (plus `from` and `to`)
#'   which is split and dispatched to each element.
#'
#' The existing vertex set and graph directedness are preserved from the
#' current network. Edge endpoints must be a subset of the existing vertex
#' names; otherwise an error is raised.
#'
#' @param object An `omic` or `omics` object to be updated.
#' @param value  The new edge table(s): a data frame/tibble for `omic`, or
#'   a named list of data frames / a single data frame with `omic` for `omics`.
#'
#' @return The updated `omic` or `omics` object with new network links.
#' @export
#' @name link<-
#' @aliases link<-,omic-method link<-,omics-method
setGeneric("link<-", function(object, value) standardGeneric("link<-"))

# -- helpers ------------------------------------------------------------------#
.is_edge_df <- function(x) {
  (is.data.frame(x) || inherits(x, "tbl_df")) &&
    all(c("from", "to") %in% colnames(x))
}

# Rebuild graph preserving vertices + directedness; enforce endpoints subset
.rebuild_graph_with_edges <- function(object, edge_df) {
  if (miss_netw(object)) {
    cli::cli_abort("No existing network available to derive vertices/directedness from.")
  }
  old_g  <- netw(object, selected = FALSE)
  vtab   <- igraph::as_data_frame(old_g, what = "vertices")
  is_dir <- igraph::is_directed(old_g)
  
  if (!("name" %in% colnames(vtab))) {
    cli::cli_abort("Vertex table must contain a {.val name} column.")
  }
  endpoints <- unique(c(edge_df$from, edge_df$to))
  unknown   <- setdiff(endpoints, vtab$name)
  if (length(unknown)) {
    cli::cli_abort(c(
      "x" = "Some edge endpoints are not present among existing vertex names.",
      "i" = "Unknown ids: {toString(utils::head(unknown, 10))}{if (length(unknown) > 10) ' …' else ''}"
    ))
  }
  
  igraph::graph_from_data_frame(d = edge_df, directed = is_dir, vertices = vtab)
}


setMethod("link<-", c("omic", "ANY"), function(object, value) {
  
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  
  if (!.is_edge_df(value)) {
    cli::cli_abort("{.arg {valueName}} must be a tibble or data frame with columns {.val from} and {.val to}.")
  }
  if ("omic" %in% colnames(value)) {
    cli::cli_abort("{.arg {valueName}} cannot include a routing column {.val omic} when assigning to a single {.cls omic}.")
  }
  
  new_g <- .rebuild_graph_with_edges(object, value)
  object@netw <- new_g
  
  validObject(object)
  object
})


setMethod("link<-", c("omics","ANY"), function(object, value){
  
  if (is.data.frame(value) || inherits(value, "tbl_df")) {
    # single tibble/data.frame with `omic` routing column
    is_assign_omics_link_tbl(object, value)
    splitted_value <- split(value, value$omic)
    splitted_value <- lapply(splitted_value, \(x) { x$omic <- NULL; x })
    
    for (i in names(object)) link(object[[i]]) <- splitted_value[[i]]
    validObject(object)
    return(object)
    
  } else if (is.list(value)) {
    # named list path (like taxa<-)
    is_list_omics_assign(object, value)
    for (i in names(object)) link(object[[i]]) <- value[[i]]
    validObject(object)
    return(object)
    
  } else {
    cli::cli_abort("{.arg value} must be either a named list or a data.frame/tibble.")
  }
})



# COMMUNITY
#------------------------------------------------------------------------------#
#' Set Community Detection Results
#'
#' This setter function updates the community detection results for `omic` objects and
#' each `omic` object within a `omics`. The community data must be an object of class `communities`
#' for `omic` objects. For `omics` objects, the community data should be a named list of `communities` objects
#' corresponding to each `omic` object within the list.
#'
#' @param object An `omic` or `omics` object.
#' @param value The new community detection results to be set, which should be compatible
#' with the structure expected by the `omic` object's community slot. For individual `omic`
#' objects, this is typically an object of class `communities` as returned by community detection
#' functions in the `igraph` package. For `omics` objects, provide a named list where each element
#' is a `communities` object corresponding to the respective `omic` object within the list.
#'
#' @return The `omic` or `omics` object with updated community detection results.
#'
#' @export
#' @importFrom methods validObject
#' @name comm<-
#' @aliases comm<-,omic-method comm<-,omics-method
setGeneric("comm<-", function(object, value) standardGeneric("comm<-"))

setMethod("comm<-", "omic", function(object, value) {
  object@comm <- value
  validObject(object)
  object
})

setMethod("comm<-", "omics", function(object, value) {
  is_list_omics_assign(object, value)
  for(i in names(object)) { object@omics[[i]]@comm <- value[[i]] }
  validObject(object)
  object
})