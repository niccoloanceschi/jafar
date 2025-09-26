
#' @useDynLib jafar, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppParallel
NULL

.onLoad <- function(libname, pkgname) {
  # Force load RcppParallel before our DLL
  requireNamespace("RcppParallel", quietly = TRUE)
}
