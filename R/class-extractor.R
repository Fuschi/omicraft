#' @include class-omic.R
NULL

# OMIC EXTRACTOR
#------------------------------------------------------------------------------#
#' Subset omic Object
#'
#' Subsets an `omic` object based on sample (i) and taxa (j) indices.
#'
#' @name subset-omic
#' @rdname subset-omic
#' @aliases [,omic-method \S4method{[}{omic}
#'
#' @usage \S4method{[}{omic}(x, i, j, ..., drop = FALSE)
#'
#' @param x An `omic` object to be subsetted.
#' @param i Indices or a logical vector selecting **samples**. If missing, all samples are included.
#' @param j Indices or a logical vector selecting **taxa**. If missing, all taxa are included.
#' @param ... Ignored. Present only for compatibility with the S4 `[` signature.
#' @param drop Ignored. Always treated as `FALSE`; passing `drop = TRUE` has no effect.
#'
#' @details
#' Subsetting by samples filters abundance matrices and sample metadata; subsetting by taxa
#' filters taxa-related data and, if present, updates dependent structures accordingly.
#' Subsetting by samples removes the `network` and `community` slots; a warning is issued.
#'
#' @return A new `omic` object containing only the requested subset.
#' @seealso \code{\link[base]{Extract}}
#' @importFrom igraph make_empty_graph cluster_fast_greedy 
#' @importFrom methods new 
#' @importFrom igraph graph_from_adjacency_matrix as_adjacency_matrix V
#' @importFrom cli cli_abort cli_warn
#' @export
setMethod(f="[", signature="omic", definition = function(x, i, j, ..., drop = FALSE)  {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Helper function to check uniqueness (only for numeric or character indices)
  check_unique <- function(indices, label) {
    if (any(duplicated(indices))) {
      cli::cli_abort("'{.var {label}}' contains duplicate values.")
    }
  }
  
  # Validate and process 'i' indices (samples)
  if (!missing(i)) {
    if (is.numeric(i)) {
      if (any(i < 1 | i > nsample(x))) {
        cli::cli_abort("{.var i} (samples) indices out of range (1..{nsample(x)}).")
      }
      check_unique(i, "i (samples)")  # Check uniqueness
    } else if (is.character(i)) {
      i <- match(i, sample_id(x))
      if (any(is.na(i))) {
        cli::cli_abort("Some sample names in {.var i} not found in {.fn sample_id}().")
      }
      check_unique(i, "i (samples)")  # Check uniqueness
    } else if (is.logical(i)) {
      if (length(i) != nsample(x)) {
        cli::cli_abort("Logical {.var i} must have length equal to the number of samples ({nsample(x)}).")
      }
      i <- which(i)  # Convert logical to numeric
    } else {
      cli::cli_abort("Invalid {.var i} (samples). Must be numeric, character, or logical.")
    }
  } else {
    i <- seq_len(nsample(x))  # Default: all samples
  }
  
  # Validate and process 'j' indices (taxa)
  if (!missing(j)) {
    if (is.numeric(j)) {
      if (any(j < 1 | j > ntaxa(x))) {
        cli::cli_abort("{.var j} (taxa) indices out of range (1..{ntaxa(x)}).")
      }
      check_unique(j, "j (taxa)")  # Check uniqueness
    } else if (is.character(j)) {
      j <- match(j, taxa_id(x))
      if (any(is.na(j))) {
        cli::cli_abort("Some taxa names in {.var j} not found in {.fn taxa_id}().")
      }
      check_unique(j, "j (taxa)")  # Check uniqueness
    } else if (is.logical(j)) {
      if (length(j) != ntaxa(x)) {
        cli::cli_abort("Logical {.var j} must have length equal to the number of taxa ({ntaxa(x)}).")
      }
      j <- which(j)  # Convert logical to numeric
    } else {
      cli::cli_abort("Invalid {.var j} (taxa). Must be numeric, character, or logical.")
    }
  } else {
    j <- seq_len(ntaxa(x))  # Default: all taxa
  }
  
  # Determine if full sample set is used to retain network and community
  full_i <- all(seq_len(nsample(x)) %in% i) || all(sample_id(x) %in% i)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # Handle empty indices cases
  if ( length(i) == 0 && length(j) == 0 ){
    
    return(new("omic")) # Return a completely new, empty omic object
    
  } else if ( length(i) == 0 && length(j) > 0 ){
    
    # Case where only taxa are selected
    new_omic <- x
    new_omic@abun <- matrix(nrow=0, ncol=0)
    new_omic@rela <- matrix(nrow=0, ncol=0)
    new_omic@norm <- matrix(nrow=0, ncol=0)
    new_omic@meta <- data.frame()
    new_omic@taxa <- x@taxa[j,,drop = F]
    new_omic@netw <- igraph::make_empty_graph(0)
    new_omic@comm = igraph::cluster_fast_greedy(igraph::make_empty_graph(0,directed=F))
    validObject(new_omic)
    return(new_omic)
    
  } else if ( length(i) > 0 && length(j) == 0 ){
    
    # Case where only samples are selected
    new_omic <- x
    new_omic@abun <- matrix(numeric(0), nrow=0, ncol=0)
    new_omic@rela <- matrix(numeric(0), nrow=0, ncol=0)
    new_omic@norm <- matrix(numeric(0), nrow=0, ncol=0)
    new_omic@meta <- x@meta[i,,drop = F]
    new_omic@taxa <- data.frame()
    new_omic@netw <- igraph::make_empty_graph(0)
    new_omic@comm = igraph::cluster_fast_greedy(igraph::make_empty_graph(0,directed=F))
    validObject(new_omic)
    return(new_omic)
    
  } else {
    
    # Standard case where both i and j are present
    # Initialize new data subsets
    if(length(x@abun) != 0) abun.new <- x@abun[i,j,drop=F] else abun.new <- x@abun
    if(length(x@rela) != 0) rela.new <- x@rela[i,j,drop=F] else rela.new <- x@rela
    if(length(x@norm) != 0) norm.new <- x@norm[i,j,drop=F] else norm.new <- x@norm
    if(length(x@meta) != 0) meta.new <- x@meta[i, ,drop=F] else meta.new <- x@meta
    if(length(x@taxa) != 0) taxa.new <- x@taxa[j, ,drop=F] else taxa.new <- x@taxa
    
    if(length(x@netw) == 0){
      
      # When network is missing only abundances and metadata are returned
      return(omic(abun = abun.new,
                  rela = rela.new,
                  norm = norm.new,
                  meta = meta.new,
                  taxa = taxa.new))
      
    } else if (length(x@netw)!=0 & full_i) {
      
      # Obtain the adjacency matrix
      if(igraph::is_weighted(x@netw)){
        adj <- igraph::as_adjacency_matrix(x@netw, attr = "weight", sparse = FALSE)
      } else {
        adj <- igraph::as_adjacency_matrix(x@netw, sparse = FALSE)
      }
      
      #Reorder vertices
      adj_reorder <- as.matrix(adj[j,j,drop=F])
      rownames(adj_reorder) <- colnames(adj_reorder) <- V(x@netw)$name[j]
      
      # Extract edge attributes as a data frame
      edge_df <- igraph::as_data_frame(x@netw, what = "edges")
      
      # Reconstruct the graph from the reordered adjacency matrix
      netw.new <- igraph::graph_from_adjacency_matrix(adj_reorder, 
                                                      mode = "undirected", 
                                                      weighted = TRUE, 
                                                      diag = FALSE)
      
      # Match reordered edges with the original edge attributes
      new_edge_df <- igraph::as_data_frame(netw.new, what = "edges")
      original_names <- paste(edge_df$from, edge_df$to, sep = "_")
      reordered_names <- paste(new_edge_df$from, new_edge_df$to, sep = "_")
      
      # Transfer edge attributes
      for (col in colnames(edge_df)) {
        if (!col %in% c("from", "to", "weight")) {  # Skip 'from', 'to', and 'weight'
          new_edge_df[[col]] <- edge_df[[col]][match(reordered_names, original_names)]
        }
      }
      
      for (col in colnames(new_edge_df)) {
        if (!col %in% c("from", "to", "weight")) {
          netw.new <- igraph::set_edge_attr(netw.new, name = col, value = new_edge_df[[col]])
        }
      }
      
      
      if(length(x@comm) != 0){
        comm.new <- x@comm
        if(is.character(j)) j <- which(taxa_id(x)%in%j)
        comm.new$membership <- x@comm$membership[j]
        comm.new$vcount <- length(comm.new$membership)
        comm.new$modularity <- NA
      } else {
        comm.new <- x@comm
      }
      
      return(omic(abun = abun.new,
                  rela = rela.new,
                  norm = norm.new,
                  meta = meta.new,
                  taxa = taxa.new,
                  netw = netw.new,
                  comm = comm.new))
      
    } else if (length(x@netw)!=0 & !full_i) {
      # SUBCASE NETWORK PRESENT AND I PRESENT
      
      cli::cli_warn("Subsetting by samples removes {.field netw} and {.field comm} slots.")
      return(omic(abun = abun.new,
                  rela = rela.new,
                  norm = norm.new,
                  meta = meta.new,
                  taxa = taxa.new))
      
    } else {
      
      cli::cli_abort("Why are you here? Only for suffering? (cit. Kaz).")
      
    }
    
  }
})