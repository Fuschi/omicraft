################################################################################
##                             SHARED SLOT CHECKS                             ##
################################################################################

#------------------------------------------------------------------------------#
#' Assert: unique row and column names
#' @keywords internal
.assert_named_and_uniqueness <- function(x, field) {
  if (is.null(rownames(x))) {
    cli::cli_abort("The {.field {field}} must have row names.",
                   class = "omic_validators_error")
  }
  if (anyDuplicated(rownames(x)) > 0) {
    cli::cli_abort("Duplicate row names in {.field {field}}.",
                   class = "omic_validators_error")
  }
  if (is.null(colnames(x))) {
    cli::cli_abort("The {.field {field}} must have column names.",
                   class = "omic_validators_error")
  }
  if (anyDuplicated(colnames(x)) > 0) {
    cli::cli_abort("Duplicate column names in {.field {field}}.",
                   class = "omic_validators_error")
  }
}

#------------------------------------------------------------------------------#
#' Assert: no reserved keywords in column names
#' @keywords internal
.OMICS_RESERVED_COLNAMES <- c("sample_id", "taxa_id", "comm_id", "link_id", 
                              "abun", "rel", "norm", "omic", "_internal_")
.assert_no_reserved_names <- function(x, field) {
  conflicted <- intersect(colnames(x), .OMICS_RESERVED_COLNAMES)
  if (length(conflicted) > 0) {
    cli::cli_abort(
      c(
        "Column name{?s} in {.field {field}} must not match reserved keyword{?s}.",
        "x" = "Found: {.val {conflicted}}"
      ),
      class = "omic_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: only values >=0
#' @keywords internal
.assert_positive_values <- function(x, field) {
  if (any(x < 0, na.rm = TRUE)) {
    cli::cli_abort("Negative values found in {.field {field}}; all values must be >= 0.",
                   class = "omic_validators_error")
  }
}

################################################################################
##                             SHARED RECIPROCAL CHECKS                       ##
################################################################################

#------------------------------------------------------------------------------#
#' Assert: reciprocal names of rows and columns must be identical
#' @keywords internal
.assert_reciprocal_dimensions <- function(x1, field1, x2, field2) {
  if (nrow(x1) != nrow(x2)) {
    cli::cli_abort(
      c(
        "Row count mismatch between {.field {field1}} and {.field {field2}}.",
        "i" = "Rows represent samples. Ensure both have the same number of samples.",
        "x" = "{.field {field1}} has {.val {nrow(x1)}} rows.",
        "x" = "{.field {field2}} has {.val {nrow(x2)}} rows."
      ),
      class = "omic_validators_error"
    )
  }
  if (ncol(x1) != ncol(x2)) {
    cli::cli_abort(
      c(
        "Column count mismatch between {.field {field1}} and {.field {field2}}.",
        "i" = "Columns represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {ncol(x1)}} columns.",
        "x" = "{.field {field2}} has {.val {ncol(x2)}} columns."
      ),
      class = "omic_validators_error"
    )
  }
  if (!identical(rownames(x1), rownames(x2))) {
    cli::cli_abort(
      c(
        "Row names of {.field {field1}} and {.field {field2}} must be identical.",
        "i" = "Row names correspond to sample IDs and must match exactly in content and order."
      ),
      class = "omic_validators_error"
    )
  }
  if (!identical(colnames(x1), colnames(x2))) {
    cli::cli_abort(
      c(
        "Column names of {.field {field1}} and {.field {field2}} must be identical.",
        "i" = "Column names correspond to taxa IDs and must match exactly in content and order."
      ),
      class = "omic_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of rows
#' @keywords internal
.assert_identical_samples <- function(x1, field1, x2, field2) {
  if (nrow(x1) != nrow(x2)) {
    cli::cli_abort(
      c(
        "Row count mismatch between {.field {field1}} and {.field {field2}}.",
        "i" = "Rows represent samples. Ensure both have the same number of samples.",
        "x" = "{.field {field1}} has {.val {nrow(x1)}} rows.",
        "x" = "{.field {field2}} has {.val {nrow(x2)}} rows."
      ),
      class = "omic_validators_error"
    )
  }
  if (!identical(rownames(x1), rownames(x2))) {
    cli::cli_abort(
      c(
        "Row names of {.field {field1}} and {.field {field2}} must be identical.",
        "i" = "Row names correspond to sample IDs and must match exactly in content and order."
      ),
      class = "omic_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of cols of x1 and rows of x2
#' @keywords internal
.assert_identical_taxa <- function(x1, field1, x2, field2) {
  if (ncol(x1) != nrow(x2)) {
    cli::cli_abort(
      c(
        "Columns count in {.field {field1}} and rows in {.field {field2}} mismatch.",
        "i" = "Columns in {.field {field1}}, rows in {.field {field2}} represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {ncol(x1)}} columns.",
        "x" = "{.field {field2}} has {.val {nrow(x2)}} rows."
      ),
      class = "omic_validators_error"
    )
  }
  if (!identical(colnames(x1), rownames(x2))) {
    cli::cli_abort(
      c(
        "Column names in {.field {field1}} and row names in {.field {field2}} must match.",
        "i" = "These names represent taxa and must be identical and in the same order."
      ),
      class = "omic_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of cols of x1 and nodes in netw
#' @keywords internal
.assert_identical_abun_nodes <- function(x1, field1, netw) {
  if (ncol(x1) != igraph::vcount(netw)) {
    cli::cli_abort(
      c(
        "Columns count in {.field {field1}} and nodes in {.field netw} mismatch.",
        "i" = "Columns in {.field {field1}}, nodes in {.field netw} represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {ncol(x1)}} columns.",
        "x" = "{.field netw} has {.val {igraph::vcount(netw)}} rows."
      ),
      class = "omic_validators_error"
    )
  }
  if (!identical(colnames(x1), igraph::V(netw)$name)) {
    cli::cli_abort(
      c(
        "Column names in {.field {field1}} and nodes names in {.field netw} must match.",
        "i" = "These names represent taxa and must be identical and in the same order."
      ),
      class = "omic_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of rows of x1 and nodes in netw
#' @keywords internal
.assert_identical_taxa_nodes <- function(x1, field1, netw) {
  if (nrow(x1) != igraph::vcount(netw)) {
    cli::cli_abort(
      c(
        "Rows count in {.field {field1}} and nodes in {.field netw} mismatch.",
        "i" = "Rows in {.field {field1}}, nodes in {.field netw} represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {nrow(x1)}} rows.",
        "x" = "{.field netw} has {.val {igraph::vcount(netw)}} rows."
      ),
      class = "omic_validators_error"
    )
  }
  if (!identical(rownames(x1), igraph::V(netw)$name)) {
    cli::cli_abort(
      c(
        "Row names in {.field {field1}} and node names in {.field netw} must match.",
        "i" = "These names represent taxa and must be identical and in the same order."
      ),
      class = "omic_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: identical zero position in the two matrices
#' @keywords internal
.assert_abun_zero_structure <- function(x1, field1, x2, field2) {
  if (!identical(x1 == 0, x2 == 0)) {
    cli::cli_abort(
      c(
        "Inconsistent zero patterns between {.field {field1}} and {.field {field2}}.",
        "i" = "Zeros must appear in the same positions in both matrices."
      ),
      class = "omic_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: consistency between network and community
#' @keywords internal
.assert_netw_comm_consistency <- function(netw, comm) {
  if (length(netw) == 0 && length(comm) > 0) {
    cli::cli_abort("The {.field comm} slot cannot be present without a {.field netw} slot.",
                   class = "omic_validators_error")
  }
  if (length(netw) > 0 && length(comm) > 0) {
    if (length(comm$membership) != igraph::vcount(netw)) {
      cli::cli_abort("The number of membership elements in {.field comm} must match the number of vertices in {.field netw}.",
                     class = "omic_validators_error")
    }
  }
}

################################################################################
##                                SLOT CHECKS                                 ##
################################################################################

#------------------------------------------------------------------------------#
#' Validate abun slot
#' @keywords internal
.validate_abun <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("The {.field abun} slot must be a numeric matrix (samples x taxa).",
                   class = "omic_validators_error")
  }
  .assert_named_and_uniqueness(x, "abun")
  .assert_no_reserved_names(x, "abun")
  .assert_positive_values(x, "abun")
}

#------------------------------------------------------------------------------#
#' Validate rela slot
#' @keywords internal
.validate_rela <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("The {.field rela} slot must be a numeric matrix (samples x taxa).",
                   class = "omic_validators_error")
  }
  .assert_named_and_uniqueness(x, "rela")
  .assert_no_reserved_names(x, "rela")
  .assert_positive_values(x, "rela")
}

#------------------------------------------------------------------------------#
#' Validate norm slot
#' @keywords internal
.validate_norm <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("The {.field norm} slot must be a numeric matrix (samples x taxa).",
                   class = "omic_validators_error")
  }
  .assert_named_and_uniqueness(x, "norm")
  .assert_no_reserved_names(x, "norm")
}

#------------------------------------------------------------------------------#
#' Validate meta slot
#' @keywords internal
.validate_meta <- function(x) {
  if (!is.data.frame(x)) {
    cli::cli_abort("The {.field meta} slot must be a {.cls data.frame} with sample metadata (samples x features).",
                   class = "omic_validators_error")
  }
  .assert_named_and_uniqueness(x, "meta")
  .assert_no_reserved_names(x, "meta")
}

#------------------------------------------------------------------------------#
#' Validate taxa slot
#' @keywords internal
.validate_taxa <- function(x) {
  if (!is.data.frame(x)) {
    cli::cli_abort("The {.field taxa} slot must be a {.cls data.frame} with taxa metadata (taxa x features).",
                   class = "omic_validators_error")
  }
  .assert_named_and_uniqueness(x, "taxa")
  .assert_no_reserved_names(x, "taxa")
}

#------------------------------------------------------------------------------#
#' Validate netw slot
#' @keywords internal
.validate_netw <- function(x) {
  if (!igraph::is_igraph(x) || !igraph::is_named(x)) {
    cli::cli_abort("The {.field netw} slot must be a {.cls igraph} with names vertices.",
                   class = "omic_validators_error")
  }
}

#------------------------------------------------------------------------------#
#' Validate comm slot
#' @keywords internal
.validate_comm <- function(x) {
  if (!inherits(x, "communities")) {
    cli::cli_abort("The {.field comm} slot must be a {.cls communities} object.",
                   class = "omic_validators_error")
  }
}

################################################################################
##                             CHECK RECIPROCAL                               ##
################################################################################
#' Validate comm slot
#' @keywords internal
.validate_reciprocal <- function(obj) {
  
  if (length(obj@abun) > 0){
    if (length(obj@rela) > 0) {.assert_reciprocal_dimensions(obj@abun, "abun", obj@rela, "rela")
                               .assert_abun_zero_structure(obj@abun, "abun", obj@rela, "rela")}
    if(length(obj@norm) > 0) .assert_reciprocal_dimensions(obj@abun, "abun", obj@rela, "norm")
    if(length(obj@meta) > 0) .assert_identical_samples(obj@abun, "abun", obj@meta, "meta")
    if(length(obj@taxa) > 0) .assert_identical_taxa(obj@abun, "abun", obj@taxa, "taxa")
    if(length(obj@netw) > 0) .assert_identical_abun_nodes(obj@abun, "abun", obj@netw)
  } 
  
  if (length(obj@rela) > 0){
    if(length(obj@norm) > 0) .assert_reciprocal_dimensions(obj@rela, "rela", obj@norm, "norm")
    if(length(obj@meta) > 0) .assert_identical_samples(obj@rela, "rela", obj@meta, "meta")
    if(length(obj@taxa) > 0) .assert_identical_taxa(obj@rela, "rela", obj@taxa, "taxa")
    if(length(obj@netw) > 0) .assert_identical_abun_nodes(obj@rela, "rela", obj@netw)
  } 
  
  if (length(obj@norm) > 0){
    if(length(obj@meta) > 0) .assert_identical_samples(obj@norm, "norm", obj@meta, "meta")
    if(length(obj@taxa) > 0) .assert_identical_taxa(obj@norm, "norm", obj@taxa, "taxa")
    if(length(obj@netw) > 0) .assert_identical_abun_nodes(obj@norm, "norm", obj@netw)
  } 
  
  if (length(obj@taxa) > 0){
    if(length(obj@netw) > 0) .assert_identical_taxa_nodes(obj@taxa, "taxa", obj@netw)
  } 
  
  .assert_netw_comm_consistency(obj@netw, obj@comm)
  
}

################################################################################
##                                 VALIDATOR                                  ##
################################################################################

#------------------------------------------------------------------------------#
#' Validity method for omic class
#' @name omic-validator
#' @keywords internal
setValidity("omic", function(object) {
  tryCatch({
    if (length(object@abun) > 0) .validate_abun(object@abun)
    if (length(object@rela) > 0) .validate_rela(object@rela)
    if (length(object@norm) > 0) .validate_norm(object@norm)
    if (length(object@meta) > 0) .validate_meta(object@meta)
    if (length(object@taxa) > 0) .validate_taxa(object@taxa)
    if (length(object@netw) > 0) .validate_netw(object@netw)
    if (length(object@comm) > 0) .validate_comm(object@comm)
    .validate_reciprocal(object)
    TRUE
  }, omic_validators_error = function(e) {
    conditionMessage(e)
  })
})
#------------------------------------------------------------------------------#
# End of validity definition


