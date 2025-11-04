setMethod("show", "omic", function(object) {
  cat("==== omic Object Summary ====\n")
  
  # -- blocchi disponibili
  A <- abun(object); R <- rela(object); N <- norm(object)
  present_blocks <- c(if (!is.null(A) && length(A)) "abun",
                      if (!is.null(R) && length(R)) "rela",
                      if (!is.null(N) && length(N)) "norm")
  present_blocks <- present_blocks[!is.na(present_blocks)]
  
  # -- dimensioni
  ns <- suppressWarnings(tryCatch(nsample(object), error = function(e) NA_integer_))
  nt <- suppressWarnings(tryCatch(ntaxa(object),   error = function(e) NA_integer_))
  cat("General Info:\n")
  cat(sprintf("  Samples: %s\n", ifelse(is.na(ns), "NA", ns)))
  cat(sprintf("  Taxa:    %s\n", ifelse(is.na(nt), "NA", nt)))
  
  # -- statistiche blocco (priorità: abun > rela > norm)
  if (length(present_blocks)) {
    X <- if ("abun" %in% present_blocks) A else if ("rela" %in% present_blocks) R else N
    total_el <- as.double(nrow(X)) * as.double(ncol(X))
    nnz <- suppressWarnings(tryCatch(sum(X != 0, na.rm = TRUE), error = function(e) NA_real_))
    zero_pct <- if (is.na(nnz) || total_el == 0) NA_real_ else 100 * (1 - (as.double(nnz) / total_el))
    na_cnt <- suppressWarnings(tryCatch(sum(is.na(X)), error = function(e) NA_integer_))
    cat(sprintf("  Active block: %s  |  Zeros: ~%s  |  NAs: %s\n",
                paste(present_blocks, collapse = ","),
                ifelse(is.na(zero_pct), "NA%", sprintf("%.2f%%", zero_pct)),
                format(na_cnt, big.mark = ",")))
  } else {
    cat("  No abundance-like data available (abun/rela/norm all empty).\n")
  }
  
  # -- meta/taxa
  M <- meta(object)
  if (!is.null(M) && ncol(M) > 0) {
    cols <- utils::head(colnames(M), 6)
    cat(sprintf("  Sample meta cols: %s%s\n", toString(cols), if (ncol(M) > 6) ", ..." else ""))
  } else cat("  Sample meta: none.\n")
  
  Tt <- taxa(object)
  if (!is.null(Tt) && ncol(Tt) > 0) {
    cols <- utils::head(colnames(Tt), 6)
    cat(sprintf("  Taxa meta cols:   %s%s\n", toString(cols), if (ncol(Tt) > 6) ", ..." else ""))
  } else cat("  Taxa meta: none.\n")
  
  # -- network e link selezionati
  rawG <- suppressWarnings(tryCatch(object@netw, error = function(e) NULL))
  raw_nodes <- if (!is.null(rawG) && inherits(rawG, "igraph")) igraph::vcount(rawG) else 0L
  raw_edges <- if (!is.null(rawG) && inherits(rawG, "igraph")) igraph::ecount(rawG) else 0L
  
  sel_active <- suppressWarnings(tryCatch(isTRUE(are_selected_links(object)), error = function(e) FALSE))
  sel_ids    <- suppressWarnings(tryCatch(if (sel_active) get_selected_links(object) else integer(0),
                                          error = function(e) integer(0)))
  selE    <- length(sel_ids)
  sel_pct <- if (sel_active && raw_edges > 0) 100 * selE / raw_edges else NA_real_
  
  G <- netw(object)
  if (!is.null(G) && inherits(G, "igraph") && igraph::vcount(G) > 0) {
    dens_shown <- suppressWarnings(tryCatch(100 * igraph::edge_density(G), error = function(e) NA_real_))
    cat(sprintf("  Network (shown): %d nodes, %d edges, density ~%s\n",
                igraph::vcount(G), igraph::ecount(G),
                ifelse(is.na(dens_shown), "NA%", sprintf("%.2f%%", dens_shown))))
    if (sel_active) {
      cat(sprintf("    Links selected: YES (%d of %d, ~%.2f%%). Showing selected subgraph.\n",
                  selE, raw_edges, sel_pct))
    } else if (raw_edges > 0L) {
      cat("    Links selected: NO (showing full graph).\n")
    }
  } else {
    if (raw_nodes > 0L) {
      cat(sprintf("  Network: empty after filtering. Raw graph had %d nodes, %d edges.\n",
                  raw_nodes, raw_edges))
      if (sel_active) {
        cat(sprintf("    Links selected: YES (%d of %d, ~%.2f%%). Selection produced empty view.\n",
                    selE, raw_edges, sel_pct))
      }
    } else {
      cat("  Network: none.\n")
    }
  }
  
  # -- communities
  Cc <- comm(object)
  if (!is.null(Cc)) {
    sizes <- suppressWarnings(tryCatch(igraph::sizes(Cc), error = function(e) NULL))
    if (length(sizes)) {
      sizes <- as.integer(sizes)
      n_iso <- sum(sizes == 1L)
      n_comm_gt1 <- sum(sizes > 1L)
      iso_pct <- if (!is.null(G) && igraph::vcount(G) > 0) 100 * (n_iso / igraph::vcount(G)) else NA_real_
      if (n_comm_gt1 > 0) {
        top_sz <- utils::head(sort(sizes[sizes > 1L], decreasing = TRUE), 6)
        cat(sprintf("  Communities (>1): %d  |  Isolated: %s (~%s)\n",
                    n_comm_gt1, n_iso,
                    ifelse(is.na(iso_pct), "NA%", sprintf("%.2f%%", iso_pct))))
        if (length(top_sz)) {
          cat(sprintf("  Top community sizes: %s%s\n",
                      toString(top_sz),
                      if (sum(sizes > 1L) > length(top_sz)) ", ..." else ""))
        }
      } else {
        cat("  Communities: all nodes isolated.\n")
      }
    } else cat("  Communities: empty.\n")
  } else cat("  Communities: none.\n")
  
  cat("==== End of omic Object Summary ====\n")
})

setMethod("show", "omics", function(object) {
  cat("==== omics Object Summary ====\n")
  
  k  <- length(object)
  nm <- names(object); has_names <- !is.null(nm) && length(nm) == k
  cat(sprintf("Elements: %d\n", k))
  
  if (k == 0L) {
    cat("Empty omics: no omic objects.\n")
    cat("==== End of omics Object Summary ====\n")
    return(invisible(NULL))
  }
  
  cat("Preview (up to 5 elements):\n")
  for (i in seq_len(min(k, 5L))) {
    tag <- if (has_names && nzchar(nm[i])) sprintf("'%s'", nm[i]) else sprintf("#%d", i)
    om  <- object[[i]]
    
    if (!methods::is(om, "omic")) {
      cat(sprintf("  - %s: [not an 'omic' instance]\n", tag))
      next
    }
    
    # -- dimensioni
    ns <- suppressWarnings(tryCatch(nsample(om), error = function(e) NA_integer_))
    nt <- suppressWarnings(tryCatch(ntaxa(om),   error = function(e) NA_integer_))
    
    # -- blocchi
    A <- abun(om); R <- rela(om); N <- norm(om)
    present_blocks <- c(if (!is.null(A) && length(A)) "abun",
                        if (!is.null(R) && length(R)) "rela",
                        if (!is.null(N) && length(N)) "norm")
    present_blocks <- present_blocks[!is.na(present_blocks)]
    
    X <- if ("abun" %in% present_blocks) A else if ("rela" %in% present_blocks) R else if ("norm" %in% present_blocks) N else NULL
    zero_pct <- NA_real_; na_cnt <- NA_integer_
    if (!is.null(X) && length(X)) {
      total_el <- as.double(nrow(X)) * as.double(ncol(X))
      nnz <- suppressWarnings(tryCatch(sum(X != 0, na.rm = TRUE), error = function(e) NA_real_))
      zero_pct <- if (is.na(nnz) || total_el == 0) NA_real_ else 100 * (1 - (as.double(nnz) / total_el))
      na_cnt <- suppressWarnings(tryCatch(sum(is.na(X)), error = function(e) NA_integer_))
    }
    
    # -- network
    rawG <- suppressWarnings(tryCatch(om@netw, error = function(e) NULL))
    raw_edges <- if (!is.null(rawG) && inherits(rawG, "igraph")) igraph::ecount(rawG) else 0L
    
    sel_active <- suppressWarnings(tryCatch(isTRUE(are_selected_links(om)), error = function(e) FALSE))
    sel_ids    <- suppressWarnings(tryCatch(if (sel_active) get_selected_links(om) else integer(0),
                                            error = function(e) integer(0)))
    selE    <- length(sel_ids)
    sel_pct <- if (sel_active && raw_edges > 0) 100 * selE / raw_edges else NA_real_
    
    G <- netw(om)
    shown_nodes <- if (!is.null(G) && inherits(G, "igraph")) igraph::vcount(G) else 0L
    shown_edges <- if (!is.null(G) && inherits(G, "igraph")) igraph::ecount(G) else 0L
    dens_shown  <- if (shown_nodes > 0) suppressWarnings(tryCatch(100 * igraph::edge_density(G), error = function(e) NA_real_)) else NA_real_
    
    # -- communities
    Cc <- comm(om)
    n_comm_gt1 <- NA_integer_; iso_pct <- NA_real_
    if (!is.null(Cc)) {
      sz <- suppressWarnings(tryCatch(igraph::sizes(Cc), error = function(e) NULL))
      if (!is.null(sz)) {
        sz <- as.integer(sz)
        n_comm_gt1 <- sum(sz > 1L)
        if (shown_nodes > 0) iso_pct <- 100 * (sum(sz == 1L) / shown_nodes)
      }
    }
    
    cat(sprintf("  - %s\n", tag))
    cat(sprintf("      dims: samples=%s, taxa=%s\n",
                ifelse(is.na(ns), "NA", ns),
                ifelse(is.na(nt), "NA", nt)))
    cat(sprintf("      data: active=%s%s\n",
                if (length(present_blocks)) paste(present_blocks, collapse = ",") else "none",
                if (!is.na(zero_pct) || !is.na(na_cnt))
                  sprintf(" | zeros~%s%s",
                          ifelse(is.na(zero_pct), "NA%",
                                 sprintf("%.2f%%", zero_pct)),
                          if (!is.na(na_cnt)) sprintf(", NAs=%s", format(na_cnt, big.mark=",")) else "")
                else ""))
    
    link_tag <- if (sel_active && raw_edges > 0) {
      sprintf("selected %d/%d (~%.2f%%)", selE, raw_edges, sel_pct)
    } else if (raw_edges > 0) {
      "all"
    } else {
      "none"
    }
    cat(sprintf("      netw: nodes=%d, edges=%d, dens~%s | links=%s\n",
                shown_nodes, shown_edges,
                ifelse(is.na(dens_shown), "NA%", sprintf("%.2f%%", dens_shown)),
                link_tag))
    cat(sprintf("      comm: groups>1=%s, iso~%s\n",
                ifelse(is.na(n_comm_gt1), "NA", n_comm_gt1),
                ifelse(is.na(iso_pct), "NA%", sprintf("%.2f%%", iso_pct))))
  }
  
  if (k > 5L) cat(sprintf("  ... and %d more.\n", k - 5L))
  cat("==== End of omics Object Summary ====\n")
})

