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
    if (length(object@omics) == 0) return(TRUE)
    
    if (!all(sapply(object@omics, inherits, "omic"))) {
      cli::cli_abort("All elements of {.cls omics} must be instances of the {.cls omic} class.")
    }
    
    if (is.null(names(object@omics)) || any(!nzchar(names(object@omics)))) {
      cli::cli_abort("All elements of {.cls omics} must be named.")
    }
    
    if (length(unique(names(object@omics))) != length(object@omics)) {
      cli::cli_abort("All element names in {.cls omics} must be unique.")
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
#' @importFrom methods new
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
  
  methods::new("omics", omics = omics)
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


#' Access or assign elements in an omics object
#'
#' @description
#' These methods allow access to or assignment of `omic` objects within an `omics` object using the `$` operator.
#'
#' @param x An object of class `omics`.
#' @param name The name of the `omic` object to access or assign within the `omics` list.
#' @param value An object of class `omic` to assign to the given name (used in `$<-` only).
#'
#' @return For `$`, returns the `omic` object with the specified name.
#' For `$<-`, returns the updated `omics` object with the new or modified entry.
#' 
#' @importFrom methods validObject
#' @aliases $,omics-method $<-,omics-method
#' @rdname omics-access
#' @export
setMethod("$<-", "omics", function(x, name, value) {
  if (!inherits(value, "omic")) {
    cli::cli_abort("Assigned value must be an {.cls omic} object.")
  }
  x@omics[[name]] <- value
  methods::validObject(x)
  x
})
