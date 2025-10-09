#' @include class-omic.R
NULL

#' Class `omics` – a collection of `omic` objects
#'
#' @description
#' `omics` is an S4 class to hold and manage a collection of `omic` objects,
#' enabling multi-dataset analyses with a consistent, type-safe interface.
#'
#' @slot omics A list whose elements are instances of class `omic`.
#' If not empty, the list must be named and names must be unique.
#'
#' @section Validity checks:
#' If not empty:
#' - all elements must inherit from class `omic`;
#' - all elements must have non-empty names;
#' - names must be unique.
#'
#' @seealso [omic()] for details about the single-omic class used as elements.
#'
#' @name omics-class
#' @rdname omics-class
#' @exportClass omics
#' @importFrom methods setClass setValidity
setClass(
  "omics",
  slots = c(omics = "list"),
  prototype = prototype(omics = list()),
  validity = function(object) {
    if (length(object@omics) == 0L) return(TRUE)

    # all elements are omic
    if (!all(vapply(object@omics, function(z) inherits(z, "omic"), logical(1)))) {
      cli::cli_abort("All elements of {.cls omics} must be instances of the {.cls omic} class.")
    }
    # names exist and are non-empty
    nm <- names(object@omics)
    if (is.null(nm) || any(!nzchar(nm))) {
      cli::cli_abort("All elements of {.cls omics} must be named with non-empty strings.")
    }
    # names are unique
    if (anyDuplicated(nm) > 0L) {
      cli::cli_abort("All element names in {.cls omics} must be unique.")
    }
    TRUE
  }
)

#' Create an `omics` object
#'
#' @description
#' Build an `omics` object to encapsulate a collection of `omic` objects.
#' Accepts multiple `omic` arguments or a single *list* of `omic`.
#'
#' @param ... `omic` objects, or a *list* of `omic` objects.
#'
#' @return An object of class `omics`. If not empty, elements must be named with
#' unique names.
#'
#' @importFrom methods new
#' @export
#' @name omics
omics <- function(...) {
  x <- list(...)

  # single list argument
  if (length(x) == 1L && is.list(x[[1L]])) {
    x <- x[[1L]]
  }

  if (length(x) > 0L) {
    # type check
    if (!all(vapply(x, function(z) inherits(z, "omic"), logical(1)))) {
      cli::cli_abort("All inputs must be {.cls omic} objects or a list of {.cls omic}.")
    }
    # names present, non-empty
    nm <- names(x)
    if (is.null(nm) || any(!nzchar(nm))) {
      cli::cli_abort("All omic objects must be named with non-empty strings.")
    }
    # unique names
    if (anyDuplicated(nm) > 0L) {
      cli::cli_abort("Names must be unique.")
    }
  }

  methods::new("omics", omics = x)
}

#------------------------------------------------------------------------------#
# `omics` General methods 
#------------------------------------------------------------------------------#

#' General methods for `omics` objects
#'
#' Basic utilities to interact with `omics` objects: coerce to list, get length,
#' get and set names.
#'
#' @section Methods:
#' \describe{
#'   \item{`as.list(x)`}{Convert an `omics` object to a plain list of `omic` objects.}
#'   \item{`as(x, 'list')`}{Coerce the `omics` object to a list.}
#'   \item{`length(x)`}{Return the number of contained `omic` objects.}
#'   \item{`names(x)`}{Return element names.}
#'   \item{`names(x) <- value`}{Set element names.}
#' }
#'
#' @param x An object of class `omics`.
#' @param value A character vector for the replacement form `names<-`.
#' @param ... Not used.
#'
#' @return
#' - `as.list()` / `as(., 'list')`: a list of `omic` objects. \cr
#' - `length()`: an integer with the number of elements. \cr
#' - `names()`: a character vector with element names. \cr
#' - `names<-()`: the updated `omics` object.
#'
#' @aliases
#' as.list.omics
#' coerce,omics,list-method
#' length,omics-method
#' names,omics-method
#' names<-,omics,character-method
#'
#' @importFrom methods setAs setMethod setReplaceMethod validObject
#' @name omics-methods
#' @rdname omics-methods
NULL

# --- as.list (S3) ---
#' @rdname omics-methods
#' @export
as.list.omics <- function(x, ...) {
  x@omics
}

# --- Coercion to list (S4 setAs) ---
# no roxigen
methods::setAs("omics", "list", function(from) from@omics)

# --- length (S4) ---
#' @rdname omics-methods
#' @export
setMethod("length", "omics", function(x) {
  length(x@omics)
})

# --- names (getter S4) ---
#' @rdname omics-methods
#' @export
setMethod("names", "omics", function(x) {
  names(x@omics)
})

# --- names<- (replacement S4) ---
#' @rdname omics-methods
#' @export
setReplaceMethod(
  "names",
  signature(x = "omics", value = "character"),
  function(x, value) {
    names(x@omics) <- value
    methods::validObject(x)
    x
  }
)


#------------------------------------------------------------------------------#
# List-like access for `omics`: $, $<-, [[, [[<-
#------------------------------------------------------------------------------#

#' Access or assign `omic` elements in an `omics` object
#'
#' List-like accessors using `$`, `$<-`, `[[`, and `[[<-`.
#'
#' @param x An `omics` object.
#' @param name Length-1 character: the element name (for `$`).
#' @param i Index or name (for `[[` / `[[<-`).
#' @param j Ignored (signature compatibility).
#' @param value An `omic` object to assign.
#' @param ... Not used.
#'
#' @return
#' - `$` / `[[`: the selected `omic` object. \cr
#' - `$<-` / `[[<-`: the updated `omics` object.
#'
#' @aliases $,omics-method
#' $<-,omics,omic-method
#' [[,omics-method
#' [[<-,omics,ANY,ANY,omic-method
#'
#' @name omics-access
#' @rdname omics-access
NULL

# --- $ (getter S4) ---
#' @rdname omics-access
#' @export
setMethod("$", "omics", function(x, name) {
  if (!name %in% names(x@omics)) {
    cli::cli_abort("No {.cls omic} object named {.val {name}} found in the {.cls omics} object.")
  }
  x@omics[[name]]
})

# --- $<- (replacement S4) ---
#' @rdname omics-access
#' @export
methods::setReplaceMethod(
  "$",
  signature(x = "omics", value = "omic"),
  function(x, name, value) {
    if (!inherits(value, "omic")) {
      cli::cli_abort("Assigned value must be an {.cls omic} object.")
    }
    x@omics[[name]] <- value
    methods::validObject(x)
    x
  }
)

# --- [[ (getter S4) ---
#' @rdname omics-access
#' @export
setMethod("[[", "omics", function(x, i, j, ...) {
  x@omics[[i]]
})

# --- [[<- (replacement S4) ---
#' @rdname omics-access
#' @export
methods::setReplaceMethod(
  "[[",
  signature(x = "omics", i = "ANY", j = "ANY", value = "omic"),
  function(x, i, j, value) {
    if (!inherits(value, "omic")) {
      cli::cli_abort("Assigned value must be an {.cls omic} object.")
    }
    x@omics[[i]] <- value
    methods::validObject(x)
    x
  }
)