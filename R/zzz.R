#' EndTB: A package for simulate, estimate, and project Tuberculosis models
#'
#' @section The model: Based on WHO's model with adaptations.
#'
#' @docType package
#' @name EndTB
#' @useDynLib EndTB, .registration=TRUE
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to EndTB")
}

.onUnload <- function(libpath) {
  library.dynam.unload("EndTB", libpath)
}