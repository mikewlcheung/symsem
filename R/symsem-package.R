#' A collection of functions for symbolic computation using 'caracas'
#' package for structural equation models and other statistical analyses.
#' Among its features is the ability to calculate the model-implied
#' covariance (and correlation) matrix and the sampling covariance matrix of
#' variable functions using the delta method.
#' @name symSEM-package
#' @note As 'caracas' uses 'SymPy" in the backend. Reserved words in SymP,
#' such as"lambda" and "I" are converted to some random strings first.
#' These random strings are converted back to R.
#' @import caracas
#' @importFrom metaSEM Diag vec2symMat as.symMatrix
#' @importFrom OpenMx vech vechs
NULL
