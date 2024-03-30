## fn <- "log(a)^b+c"
## fn <- c("x^2 + log(y)*z", "x+y+sin(z)")

#' Compute the Variance-Covariance Matrix of Functions using the first-order
#' Delta Method
#'
#' It computes the variance-covariance matrix of functions using the
#' first-order delta method.
#'
#'
#' @param fn A function in character strings or a vector of functions.
#' @param Covvars Variance-covariance matrix of the variables. Users must
#' ensure the order of variables is the same as that in \code{vars}; 
#' Otherwise, the results are likely incorrect. If it is not specified,
#' they are automatically generated.
#' @param vars A vector of characters of the random variables. If the random
#' variables are not listed in `vars`, they are treated as constants. If `vars`
#' is missing, all names in `RAM` are treated as random variables.
#' @param Var.name Name of the variances.
#' @param Cov.name Name of the covariances.
#' @param simplify Attempt to simplify the output.
#' @return Variance-covariance matrix of the functions.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @export
#' @examples
#' 
#' #### Fisher-z-transformation
#' fn  <- "0.5*log((1+r)/(1-r))"
#'
#' ## Sampling variance of r
#' Covvars <- "(1-r^2)^2/n"
#'
#' deltamethod(fn=fn, Covvars=Covvars, vars="r")
#' ## $fn
#' ##     [,1]
#' ## fn1 "0.5*log((r+1)/(1-r))"
#'
#' ## $Covfn
#' ##     fn1
#' ## fn1 "1/n"
#'
#' ## $vars
#' ## [1] "r"
#'
#' ## $Covvars
#' ##   r
#' ## r "(1-r^2)^2/n"
#'
#' ## $Jmatrix
#' ##     r
#' ## fn1 "(0.5*(1-r+r+1)*(1-r))/((1-r)^2*(r+1))"
#'
#' #### Raw mean difference: y_treatment - y_control
#' fn <- "yt - yc"
#'
#' ## Sampling covariance matrix
#' ## S2p: pooled variance
#' ## nt: n_treatment
#' ## nc: n_control
#' Covvars <- matrix(c("S2p/nt", 0,
#'                     0, "S2p/nc"),
#'                   ncol=2, nrow=2)
#'
#' deltamethod(fn=fn, Covvars=Covvars, vars=c("yt", "yc"))
#' ## $fn
#' ##     [,1]
#' ## fn1 "yt-yc"
#'
#' ## $Covfn
#' ##     fn1
#' ## fn1 "(S2p*nt+S2p*nc)/(nt*nc)"
#'
#' ## $vars
#' ## [1] "yt" "yc"
#'
#' ## $Covvars
#' ##    yt       yc
#' ## yt "S2p/nt" "0"
#' ## yc "0"      "S2p/nc"
#'
#' ## $Jmatrix
#' ##     yt  yc
#' ## fn1 "1" "-1"
#'
#' #### log(odds)
#' fn <- "log(p/(1-p))"
#'
#' ## Sampling variance of p
#' Covvars <- "p*(1-p)/n"
#'
#' ## Though it is correct, the simplification does not work well.
#' deltamethod(fn=fn, Covvars=Covvars, vars="p")
#' ## $fn
#' ##     [,1]
#' ## fn1 "log(p/(1-p))"
#'
#' ## $Covfn
#' ##     fn1
#' ## fn1 "(3*p^2-p^3-3*p+1)/((p^4-4*p^3+6*p^2-4*p+1)*p*n)"
#'
#' ## $vars
#' ## [1] "p"
#'
#' ## $Covvars
#' ##   p
#' ## p "(p*(1-p))/n"
#'
#' ## $Jmatrix
#' ##     p
#' ## fn1 "((1-p+p)*(1-p))/((1-p)^2*p)"
deltamethod <- function(fn, Covvars, vars, Var.name="V", Cov.name="C",
                        simplify=TRUE) {
    
  ## Univariate or multivariate
  fn.p <- length(fn)
  
  ## function names
  fn.names <- paste0("fn", seq_len(fn.p))
  
  ## convert it to a column vector
  fn <- matrix(fn, ncol=1)
  rownames(fn) <- fn.names

  fn.S <- as_sym(fn)
  
  ## Get the variable names
  varlist <- sort(unique(all.names(parse(text=fn), functions=FALSE)))
  if (missing(vars)) {
    vars <- varlist
  } else {
    if (any(!vars %in% varlist)) {
      stop("Some of \"vars\" do not agree with the names in \"fn\".\n")
    }
  }
  
  ##
  if (missing(Covvars)) {
    ## Variance covariance matrix of x
    Covvars <- outer(vars, vars, function(y, z) paste0(Cov.name, y, z))
    metaSEM::Diag(Covvars) <- paste0(Var.name, vars)
    ## Make it symmetric
    Covvars <- metaSEM::vec2symMat(OpenMx::vech(Covvars))
    Covvars.S <- caracas::as_sym(Covvars)
  } else {
    Covvars.S <- caracas::as_sym(Covvars)
  }

  ## Jacobian matrix
  Jmatrix <- caracas::jacobian(fn.S, vars)

  Covfn <- Jmatrix %*% Covvars.S %*% t(Jmatrix)
  
  if (simplify) {Covfn <- caracas::simplify(Covfn)}
  
  Covfn <- caracas::as_character_matrix(Covfn)
  dimnames(Covfn) <- list(fn.names, fn.names)
  Covvars <- as.matrix(Covvars)
  dimnames(Covvars) <- list(vars, vars)

  Jmatrix <- caracas::as_character_matrix(Jmatrix)
  dimnames(Jmatrix) <- list(fn.names, vars)

  list(fn=fn, Covfn=Covfn, vars=vars, Covvars=Covvars, Jmatrix=Jmatrix)
}

