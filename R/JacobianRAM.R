#' Compute a Jacobian Matrix of the Implied Covariance/Correlation Matrix based on a RAM model.
#'
#' It computes a symbolic Jacobian matrix of the model-implied covariance (or correlation) matrix in
#' SEM using the RAM specification.
#'
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}
#' @param vars A vector of characters of the random variables. If the random
#' variables are not listed in `vars`, they are treated as constants. If `vars`
#' is missing, all names in `RAM` are treated as random variables.
#' @param corr Whether the model implied matrix is covariance (default) or
#' correlation structure.
#' @return A Jacobian matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @export
#' @examples
#'
#' #### A mediation model
#' model1 <- "y ~ c*x + b*m
#'            m ~ a*x
#'            ## Means
#'            y ~ b0*1
#'            m ~ m0*1
#'            x ~ x0*1"
#'
#' RAM1 <- metaSEM::lavaan2RAM(model1)
#'
#' ## Model-implied covariance matrix and mean structure
#' JacobianRAM(RAM1, corr=FALSE)
#'
#' ## Model-implied correlation matrix
#' JacobianRAM(RAM1, corr=TRUE)
#'
#' #### A CFA model
#' model2 <- "f =~ x1 + x2 + x3 + x4#'
#'            ## Mean
#'            f ~ fmean*1"
#'
#' RAM2 <- metaSEM::lavaan2RAM(model2)
#'
#' ## Model-implied covariance matrix
#' JacobianRAM(RAM2, corr=FALSE)
#'
#' ## Model-implied correlation matrix
#' JacobianRAM(RAM2, corr=TRUE)
#'
JacobianRAM <- function(RAM, vars, corr=FALSE) {
  RAM.old <- RAM
  
  ## Convert reserved words to random strings
  RAM$A <- .convert(RAM$A, type="sympy2random")
  RAM$S <- .convert(RAM$S, type="sympy2random")
  RAM$M <- .convert(RAM$M, type="sympy2random")
  dim(RAM$M) <- dim(RAM.old$M)
  dimnames(RAM$M) <- dimnames(RAM.old$M)

  ## Keep the random strings as sympy reserved words are not allowed in
  ## Jacobian in Sympy
  x <- impliedS(RAM, corr=corr, convert=FALSE)
   
  vecM <- x$Mu
  Sigma <- x$Sigma
  labels.mu <- paste0("Mu_", colnames(vecM))
  
  p <- ncol(Sigma)
  labels.sigma <- outer(seq_len(p), seq_len(p), function(x, y) paste0("Cov", x, "_", y))
  
  ## Vector of correlation matrix only
  if (corr) {        
    vecS <- matrix(OpenMx::vechs(Sigma), ncol=1)
    labels <- OpenMx::vechs(labels.sigma)
  } else {
    ## Vector of covariance matrix and then vector of means     
    vecS <- matrix(c(OpenMx::vech(Sigma), vecM), ncol=1)
    labels <- c(OpenMx::vech(labels.sigma), labels.mu)
  }
  
  varlist <- sort(all.vars(parse(text=vecS)))    
  vecS <- caracas::as_sym(vecS)
  
  ## Ensure the column order of the Jmatrix is sorted according to the parameters
  if (missing(vars)) {
    vars <- sort(varlist)
  } else {
    vars <- .convert(sort(vars), type="sympy2random")
    if (any(!vars %in% varlist)) {
      stop("Some of \"vars\" do not agree with those names in \"x\".\n")
    }
  }
  
  ## Jmatrix <- caracas::jacobian(vecS, vars)
  ## It cannot handle many vars, say 15
  ## Error in py_call_impl(callable, call_args$unnamed, call_args$named) : 
  ## ValueError: too many values to unpack (expected 2)
  
  ## Adhoc solution to split vars to batches
  vars.batch <- split(x=vars, 
                      f=rep(1:length(vars), each=10, length.out = length(vars)) )
  J.batch <- lapply(vars.batch, function(x) caracas::jacobian(vecS, x))
  Jmatrix <- do.call(cbind, J.batch)
  Jmatrix <- caracas::as_character_matrix(Jmatrix)
  Jmatrix <- .convert(Jmatrix, type="random2sympy")
  dimnames(Jmatrix) <- list(labels, .convert(vars, type="random2sympy"))   
  Jmatrix
}

