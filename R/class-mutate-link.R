#'@include class-omic.R class-omics.R class-links.R
NULL

#' Mutate Link Data for an omic object
#'
#' @description
#' `mutate_link,omic-method` provides a data-frame-like interface for mutating 
#' edge-level attributes in an igraph network associated with an `omic` object. 
#'
#' It can utilize columns from the "from" and "to" nodes via local helpers:
#'
#' - `from(var)`: Pull the column `var` from the node on the "from" side.
#' - `to(var)`: Pull the column `var` from the node on the "to" side.
#' - `either(expr)`: For each edge, returns `TRUE` if \emph{either} side meets `expr`.
#' - `both(expr)`: For each edge, returns `TRUE` if \emph{both} sides meet `expr`.
#' - `neither(expr)`: For each edge, returns `TRUE` if \emph{neither} side meets `expr`.
#' - `one(expr)`: Returns \code{TRUE} if \emph{exactly one} side satisfies `expr`.
#'
#' After constructing (or modifying) columns, the new columns are attached as
#' edge attributes in the igraph object. The `omic` object is then returned
#' with an updated network.
#'
#' @param object An `omic` object containing an igraph network (retrieved by
#'   `netw(object)`) and taxon data.
#' @param ... `dplyr`-style mutate expressions. Inside these expressions, you 
#'   can call the local helpers described above.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the link grouping
#'        from the object at the end.
#' @param .deselect Logical, default `FALSE`. If `TRUE`, remove the link selection
#'        from the object at the end.
#'
#' @return The original `omic` object, updated with new or modified edge attributes.
#'
#' @rdname mutate_link
#' @aliases mutate_link,omic-method mutate_link,omics-method
#' @importFrom dplyr mutate group_by
#' @importFrom igraph set_edge_attr E
#' @importFrom rlang .data
#' @export
setGeneric("mutate_link", function(object, ...) {standardGeneric("mutate_link")})

setMethod("mutate_link", "omic", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  # 1) Prepare the edges, expressions, and local helpers via .link_prepare()
  setup <- .link_prepare(object, ...)
  g     <- setup$graph
  edges <- setup$edges
  quos  <- setup$quos
  
  if (is_link_grouped(object)) {
    edges <- edges %>%
      dplyr::mutate(`_internal_` = get_grouped_link(object)) %>%
      dplyr::group_by(.data[["_internal_"]])
  }
  
  # 3) Perform the mutate on (potentially) grouped edges
  edges <- dplyr::mutate(edges, !!!quos)
  
  # 4) Attach any new columns as edge attributes
  edge_cols <- setdiff(names(edges), c("from", "to", "_internal_"))
  g0 <- netw(object, selected = FALSE)
  for (col_name in edge_cols) {
    g0 <- igraph::set_edge_attr(
      g0,
      name  = col_name,
      index = igraph::E(g0),
      value = edges[[col_name]]
    )
  }
  
  # 5) Update the igraph object in 'object' and return
  netw(object) <- g0
  
  if(isTRUE(.ungroup)) object <- ungroup_link(object)
  if(isTRUE(.deselect)) object <- deselect_link(object)
  object
})


setMethod("mutate_link", "omics", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  # 1) Prepare the edges, expressions, and local helpers via .link_prepare()
  setup <- .link_prepare(object, ...)
  g <- setup$graph
  edges <- setup$edges
  quos  <- setup$quos
  
  if(is_link_grouped(object)) {
    edges <- edges %>%
      dplyr::mutate(`_internal_` = get_grouped_link(object)) %>%
      dplyr::group_by(.data[["_internal_"]])
  }
  
  # 3) Perform the mutate on (potentially) grouped edges
  edges <- dplyr::mutate(edges, !!!quos)
  edges_names <- names(edges)
  
  edges <- split(edges, edges$omic)
  edges <- sapply(edges, \(x){
    x$omic <- NULL
    return(x)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # 4) Attach any new columns as edge attributes
  edge_cols <- setdiff(edges_names, c("from", "to", "_internal_"))
  g0 <- netw(object, selected = FALSE)
  for (nm in names(object)){
    for (col_name in edge_cols) {
      g0[[nm]] <- igraph::set_edge_attr(
        g0[[nm]],
        name  = col_name,
        index = igraph::E(g0[[nm]]),
        value = edges[[nm]][[col_name]]
      )
    }
    # 5) Update the igraph object in 'object' and return
    netw(object[[nm]]) <- g0[[nm]]
  }
  
  if(isTRUE(.ungroup)) object <- ungroup_link(object)
  if(isTRUE(.deselect)) object <- deselect_link(object)
  object
})

  