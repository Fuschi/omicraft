# Define omics class
#------------------------------------------------------------------------------#
#' Class omics
#'
#' @description
#' `omics` is an S4 class designed to hold and manage a collection of `omic` objects,
#' providing a structured way to handle multiple NGS datasets. It ensures that all
#' elements within the list adhere to the `omic` class structure, making it useful for
#' analyses involving multiple datasets. This class allows for consistent and
#' type-safe manipulation of multiple `omic` objects.
#'
#' @slot omics A list where each element is an instance of the `omic` class.
#' If the list is not empty, each `omic` object within the list must be named to facilitate
#' easy access and management.
#'
#' @section Validity Checks:
#' The class includes validity checks to ensure that, if not empty, all elements within the
#' `omics` slot are instances of the `omic` class and are named. These checks help maintain
#' integrity and usability, preventing errors during downstream analyses.
#'
#' @note
#' For clarity, note that the class is named `omics` and it contains a slot named `omics` as well.
#'
#' @seealso \code{\link[=omic]{omic}} for details on the `omic` class and its functionalities.
#'
#' @name omics-class
#' @rdname omics-class
#' @exportClass omics
setClass(
  "omics",
  slots = c(omics = "list"),
  prototype = prototype(omics = list()),
  validity = function(object) {
    errors <- character()
    omic_errors <- character()
    
    if (length(object@omics) != 0) {
      if (!all(sapply(object@omics, inherits, "omic"))) {
        errors <- c(errors, "All elements of `omics` must be instances of the `omic` class.")
      }
      if (is.null(names(object@omics)) || any(!nzchar(names(object@omics)))) {
        errors <- c(errors, "All elements in `omics` must be named.")
      }
      if (length(unique(names(object@omics))) != length(object@omics)) {
        errors <- c(errors, "All element names must be unique.")
      }
      
      for (i in seq_along(object@omics)) {
        omic_obj <- object@omics[[i]]
        omic_name <- names(object@omics)[i]
        
        omic_validity <- tryCatch(
          {
            validObject(omic_obj, test = TRUE)
            TRUE
          },
          error = function(e) {
            e$message
          }
        )
        
        if (omic_validity != TRUE) {
          if (is.null(omic_name) || !nzchar(omic_name)) omic_name <- "NA"
          omic_errors <- c(omic_errors, sprintf("==== '%s' error report ====", omic_name), omic_validity)
        }
      }
      
      all_errors <- character()
      if (length(errors) > 0) {
        errors <- lapply(errors, function(x) paste("- ", x, sep = ""))
        all_errors <- c("\n==== omics error report ====", paste(errors, collapse = "\n"))
      }
      if (length(omic_errors) > 0) {
        all_errors <- c(all_errors, omic_errors)
      }
      if (length(all_errors) > 0) {
        return(paste(all_errors, collapse = "\n"))
      }
    }
    TRUE
  }
)

# Define constructor for omics
#------------------------------------------------------------------------------#
#' Create an omics Object
#'
#' @description
#' Constructs an `omics` object to encapsulate and manage a collection of `omic` objects.
#' This function ensures that all provided `omic` objects are validated and organized into a
#' structured list, facilitating easy access to multiple datasets. The constructor
#' supports the initialization of an empty list, individual `omic` objects, or a list of `omic` objects.
#'
#' @param ... Optional `omic` objects or a list of `omic` objects to be included in the `omics` object.
#'
#' @return An `omics` object containing any provided `omic` objects. If not empty,
#' each `omic` object must be named.
#'
#' @export
#' @name omics
omics <- function(...) {
  omics <- list(...)
  
  if (length(omics) == 1 && is.list(omics[[1]]) && all(sapply(omics[[1]], inherits, "omic"))) {
    omics <- omics[[1]]
  }
  
  if (length(omics) != 0 && (is.null(names(omics)) || any(!nzchar(names(omics))))) {
    cli::cli_abort("All omic objects must be named.")
  }
  
  new("omics", omics = omics)
}

# Convert an omics Object to a List
#------------------------------------------------------------------------------#
#' Convert an omics Object to a List
#'
#' @param x An `omics` object to be converted to a list.
#' @param ... Not used.
#' @return A list of `omic` objects.
#' @export
as.list.omics <- function(x, ...) {
  x@omics
}

#' @export
setAs("omics", "list", function(from) from@omics)

# Access or Assign Elements in an omics Object
#------------------------------------------------------------------------------#
#' @title Access or Assign Elements in an omics Object
#' @name omics-access
#' @rdname omics-access
#' @aliases $,omics-method $<-,omics-method
#' @export
setMethod("$", "omics", function(x, name) {
  if (!name %in% names(x@omics)) {
    cli::cli_abort("No {.cls omic} object named {.val {name}} found in the {.cls omics} object.")
  }
  x@omics[[name]]
})


#' @rdname omics-access
#' @export
setMethod("$<-", "omics", function(x, name, value) {
  if (!inherits(value, "omic")) {
    cli::cli_abort("Assigned value must be an {.cls omic} object.")
  }
  x@omics[[name]] <- value
  validObject(x)
  x
})
