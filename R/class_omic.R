setOldClass("igraph")
setOldClass("communities")
################################################################################
# CLASS OMIC
################################################################################
#' Class omic
#'
#' An S4 class designed for organizing and managing omics data in a tidy and modular fashion.
#' The class stores key components such as abundance matrices (raw, relative, normalized),
#' sample metadata, feature annotations, network structures, and community assignments.
#'
#' This class is designed to integrate omics data — abundances, metadata, annotations, networks —
#' with structural consistency checks. It performs no modeling or normalization and focuses on
#' tidy-style manipulation and coherent data export.
#'
#' @slot abun A numeric matrix of raw abundance values (samples x features).
#' @slot rela A numeric matrix of relative abundances.
#' @slot norm A numeric matrix of normalized abundances.
#' @slot meta A data.frame with sample metadata.
#' @slot taxa A data.frame with feature (e.g. taxon, gene) metadata.
#' @slot netw An igraph object representing a network among features or samples.
#' @slot comm A communities object (from igraph) representing community structure.
#'
#' @section Reserved Keywords:
#' The `omic` class reserves the following keywords for internal use and method compatibility:
#' `sample_id`, `taxa_id`, `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`, `omic`.
#' Avoid using these as column names in your metadata or abundance matrices.
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods setClass
#'
#' @name omic-class
#' @rdname omic-class
#' @exportClass omic
setClass("omic",
         slots = c(
           abun = "ANY",
           rela = "ANY",
           norm = "ANY",
           meta = "ANY",
           taxa = "ANY",
           netw = "ANY",
           comm = "ANY"
         ),
         prototype = prototype(
           abun = matrix(numeric(0), nrow = 0, ncol = 0),
           rela = matrix(numeric(0), nrow = 0, ncol = 0),
           norm = matrix(numeric(0), nrow = 0, ncol = 0),
           meta = data.frame(),
           taxa = data.frame(),
           netw = igraph::make_empty_graph(0, directed = F),
           comm = igraph::cluster_fast_greedy(
             igraph::make_empty_graph(0, directed = F)
           )
         )
)


