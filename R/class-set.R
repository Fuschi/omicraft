#' Set slots of an `omic` object using tidy-style expressions
#'
#' @description
#' `set_omic()` lets you compute and assign one or more slots of an `omic`
#' object using rlang-style expressions that can reference existing slots by
#' their bare names (e.g., `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`).
#' Expressions are evaluated **sequentially** in the order provided, so later
#' expressions can depend on values assigned by earlier ones.
#'
#' @param object An object of class `omic`.
#' @param ... Named expressions, where each name must be one of the allowed
#' slot names: `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`.
#'
#' @return The modified `omic` object with updated slots.
#'
#' @section How name resolution works:
#' Within each expression, the symbols `abun`, `rela`, `norm`, `meta`, `taxa`,
#' `netw`, and `comm` resolve to the **current** working values (initially taken
#' from the object, then updated as each expression is assigned).
#' 
#' @export
setGeneric("set_omic", function(object, ...) standardGeneric("set_omic"))


#' @rdname set_omic
#' @export
setGeneric("set_omic", function(object, ...) standardGeneric("set_omic"))


#' @rdname set_omic
#' @export
setMethod("set_omic", signature(object = "omic"), function(object, ...) {
  
  dots <- rlang::enquos(..., .named = TRUE)
  
  allowed <- c("abun", "rela", "norm", "meta", "taxa", "netw", "comm")
  if (length(dots) == 0L) {
    cli::cli_abort("No expressions supplied to {.fn set_omic}.")
  }
  nm <- names(dots)
  if (is.null(nm) || any(nm == "")) {
    cli::cli_abort("All arguments to {.fn set_omic} must be named.")
  }
  unknown <- setdiff(nm, allowed)
  if (length(unknown)) {
    cli::cli_abort(c(
      "Some targets are not valid {'.omic'} slots:",
      "x" = paste0(unknown, collapse = ", "),
      "i" = "Valid targets are: {paste(allowed, collapse = ', ')}"
    ))
  }
  
  # mask iniziale dai valori correnti
  mask <- list(
    abun = abun(object),
    rela = rela(object),
    norm = norm(object),
    meta = meta(object),
    taxa = taxa(object),
    netw = netw(object),
    comm = comm(object)
  )
  
  # valuta e assegna in sequenza
  for (target in nm) {
    quo <- dots[[target]]
    
    value <- tryCatch(
      rlang::eval_tidy(quo, data = mask, env = rlang::get_env(quo)),
      error = function(e) {
        cli::cli_abort(c(
          "Failed while computing {.field {target}} in {.fn set_omic}.",
          "x" = conditionMessage(e),
          "i" = "Within expressions you can reference: {paste(allowed, collapse = ', ')}"
        ))
      }
    )
    
    # usa i tuoi replacement accessors: abun(object) <- value, ecc.
    setter <- tryCatch(match.fun(paste0(target, "<-")), error = function(e) NULL)
    if (is.null(setter)) {
      cli::cli_abort(c(
        "Replacement accessor not found for {.field {target}}.",
        "i" = "Expected a setter function named {target}<(object, value)."
      ))
    }
    object <- setter(object, value)   # <- assegna allo slot
    mask[[target]] <- value           # <- aggiorna la data mask
  }
  
  # opzionale: validazione finale
  tryCatch(methods::validObject(object), error = function(e) {
    cli::cli_abort(c(
      "The updated object failed validity after {.fn set_omic}.",
      "x" = conditionMessage(e)
    ))
  })
  
  object
})