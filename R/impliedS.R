#' Compute a Symbolic Model-Implied Covariance/Correlation Matrix
#'
#' It computes a symbolic model-implied covariance (or correlation) matrix in
#' SEM using the RAM specification inputs. 
#'
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}
#' @param corr Whether the model implied matrix is covariance (default) or
#' correlation structure.
#' @param replace.constraints Whether to replace the parameters with the constraints
#' in the \code{mxalgebras} slot. Suppose the formula is \code{para1==para2+para3},
#' \code{para1} will be replaced by \code{para2+para3} if this argument is \code{TRUE}.
#' @param convert Whether to convert random strings back to parameters. For internal
#' use only. Users unlikely need to use this argument. 
#' @return A list of object with class \code{implieS}. It stores the A, S, and F
#' matrices and the model implied covariance (or correlation) matrix and the
#' vector of the means.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @export
#' @examples
#' \dontrun{
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
#' impliedS(RAM1, corr=FALSE)
#'
#' ## Model-implied correlation matrix
#' impliedS(RAM1, corr=TRUE)
#'
#' #### A CFA model
#' model2 <- "f =~ x1 + x2 + x3 + x4
#'            ## Mean
#'            f ~ fmean*1"
#'
#' RAM2 <- metaSEM::lavaan2RAM(model2)
#'
#' ## Model-implied covariance matrix
#' impliedS(RAM2, corr=FALSE)
#'
#' ## Model-implied correlation matrix
#' impliedS(RAM2, corr=TRUE)
#' }
impliedS <- function(RAM, corr=FALSE, replace.constraints=FALSE, convert=TRUE) {
  ## Remove starting values and "*" in RAM
  ## When there is no label in lavaan syntax, the starting values are
  ## considered as fixed values.
  RAM1 <- metaSEM::as.symMatrix(RAM)

  if (replace.constraints & !is.null(RAM1$mxalgebras)) {
    x <- RAM1$mxalgebras

    ## Get all constraints
    constraints <- x[grep("^constraint", names(x))]
    
    if (length(constraints)!=0) {
      ## Replace the constraints
      for (i in seq_len(length(constraints))) {
        str <- as.character(constraints[[i]]$formula)

        RAM1$A <- gsub(str[2], str[3], RAM1$A, fixed=TRUE)
        RAM1$S <- gsub(str[2], str[3], RAM1$S, fixed=TRUE)
        RAM1$M <- gsub(str[2], str[3], RAM1$M, fixed=TRUE)
      }
    }
  }
  
  ## No. of observed and latent variables
  p <- nrow(RAM1$S)
  A <- caracas::as_sym(.convert(RAM1$A, type="sympy2random"))
  S_new <- .convert(RAM1$S, type="sympy2random")
  S <- caracas::as_sym(S_new)
  F <- caracas::as_sym(RAM1$F)

  ## Inverse(I-A)
  invIA <- caracas::inv(caracas::as_sym(diag(1L,p)) - A)
  ## Model implied covariance matrix with observed and latent variables
  SigmaAll <- invIA %*% S %*% t(invIA)
    
  if (corr) {
    ## When it is NA after as.numeric, it is a parameter.
    ## Index of error variances in the diagonals of S including latent variables
    index <- suppressWarnings(which(is.na(as.numeric(diag(S_new)))))
    Slabels <- diag(S_new)[index]

    ## representation in R characters
    Sigma.R <- as_character_matrix(SigmaAll)
        
    ## Assuming there are at least one free error variance in S
    ## including latent variables
    for (i in seq_along(index)) {
      j <- index[i]
      ## Sigma.R[i, i]: 1 + ErrVar - (sigma)
      ## It will become 1 on the diagonals the after simplification.
      ## For off-diagonal elements, the ErrVar are replaced by terms without ErrVar.
      ## For example, Slabels[1]="mWITHm", Sigma.R[2,2]="a^2+mWITHm".
      ## After replacing "mWITHm" with "(1+mWITHm-(a^2+mWITHm))",
      ## it becomes "a^2+(1+mWITHm-(a^2+mWITHm))".
      ## It becomes 1 on the diagonal after simplification.
      Sigma.R <- gsub(Slabels[i],
                      paste0("(1+", Slabels[i], "-(", Sigma.R[j, j], "))"),
                      Sigma.R)
    }
    SigmaAll <- caracas::as_sym(Sigma.R)

    ## Handle cases when the diagonals are numbers but not 1, e.g., var=4.
    ## SD: inverse of sds
    if (any(as_character_matrix(diag(SigmaAll)) != 1)) {
      SDinv <- caracas::as_diag(1/sqrt(diag(SigmaAll)))
      SigmaAll <- SDinv %*% SigmaAll %*% SDinv
    }
    
    ## Means are zeros as it is a correlation matrix
    Mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))

    ## Corr=FALSE
  } else {
    ## M is a zero matrix
    if (all(RAM1$M==0)) {
      Mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
    } else {
      M <- caracas::as_sym(.convert(RAM1$M, type="sympy2random"))
      Mu <- M %*% t(F %*% invIA)
      Mu <- caracas::as_character_matrix(Mu)
    }
  }

  SigmaAll <- caracas::simplify(SigmaAll)

  ## Convert it into R matrix of string
  SigmaAll <- caracas::as_character_matrix(SigmaAll)
    
  if (convert) {
    SigmaAll <- .convert(SigmaAll, type="random2sympy")
    Mu <- .convert(Mu, type="random2sympy")
  }

  if (!is.matrix(SigmaAll)) SigmaAll <- as.matrix(SigmaAll)
  dimnames(SigmaAll) <- list(colnames(RAM$F), colnames(RAM$F))
    
  index <- colSums(RAM$F)==1
  SigmaObs <- SigmaAll[index, index]

  if (!is.matrix(SigmaObs)) SigmaObs <- as.matrix(SigmaObs)
  dimnames(SigmaObs) <- list(rownames(RAM$F), rownames(RAM$F))
    
  Mu <- matrix(Mu, nrow=1, dimnames=list("1", colnames(RAM$M)[colSums(RAM$F)==1]))
    
  if (suppressWarnings(all(!is.na(as.numeric(RAM1$A))))) {
    Amatrix <- matrix(as.numeric(RAM1$A), nrow=nrow(RAM1$A), ncol=ncol(RAM1$A),
                      dimnames = dimnames(RAM1$A))
  } else {
    Amatrix <- RAM1$A
  }
    
  if (suppressWarnings(all(!is.na(as.numeric(RAM1$S))))) {
    Smatrix <- matrix(as.numeric(RAM1$S), nrow=nrow(RAM1$S), ncol=ncol(RAM1$S),
                      dimnames = dimnames(RAM1$S))
  } else {
    Smatrix <- RAM1$S
  }

  out <- list(A=Amatrix, S=Smatrix, F=RAM$F, M=RAM1$M, Sigma=SigmaObs,
              SigmaAll=SigmaAll, Mu=Mu, corr=corr)
  class(out) <- "impliedS"
  out
}

#' @export
print.impliedS <- function(x, ...) {
  if (!is.element("impliedS", class(x)))
    stop("\"x\" must be an object of class \"impliedS\".")
  cat("Correlation matrix:", x$corr) 
  cat("\n\nAmatrix:\n")
  print(x$A)
  cat("\nSmatrix:\n")
  print(x$S)
  cat("\nFmatrix:\n")
  print(x$F)
  cat("\nMmatrix:\n")
  print(x$M)
  cat("\nModel implied covariance matrix (Sigma):\n")
  print(x$Sigma)
  cat("\nModel implied mean vector (Mu):\n")
  print(x$Mu)
}

## Alternative version.
## It solves system of equations of error variances.
## Twice longer than the original version.
# impliedS2 <- function(RAM, corr=FALSE, simplify=TRUE) {
# 
#   ## Ad-hoc replacement of reserved words in sympy
#   adhoc_replace <- function(x) {
#     original <- c("1i", "exp(1)")
#     replace <- c("I", "E")
# 
#     for (i in seq_along(original)) {
#       x <- gsub(original[i], replace[i], x, fixed=TRUE)
#     }
#     x
#   }
# 
#   ## Remove starting values and "*" in RAM
#   RAM1 <- metaSEM::as.symMatrix(RAM)
# 
#   ## No. of observed and latent variables
#   p <- nrow(RAM1$S)
#   A <- caracas::as_sym(.convert_reserved(RAM1$A, type="sympy2random"))
#   S <- caracas::as_sym(.convert_reserved(RAM1$S, type="sympy2random"))
#   F <- caracas::as_sym(RAM1$F)
# 
#   ## Inverse(I-A)
#   invIA <- caracas::inv(caracas::as_sym(diag(1L,p)) - A)
#   ## Model implied covariance matrix with latent variables
#   SigmaAll <- invIA %*% S %*% t(invIA)
# 
#   if (corr) {
#     ## F*S*t(F) rather than S
#     ## Index of error variances in the diagonals of S
#     index <- suppressWarnings(which(is.na(as.numeric(diag(RAM1$S)))))
#     Slabels <- diag(RAM1$S)[index]
#     ## Solve for error variances
#     sol <- caracas::solve_sys(lhs=diag(SigmaAll)[index],
#                               rhs=as_sym(matrix(1, nrow=length(index))),
#                               vars=Slabels)
#     ## Replace the error variances with the solutions
#     SigmaAll <- caracas::subs(SigmaAll, sol[[1]])
#     SigmaAll <- caracas::simplify(SigmaAll)
# 
#     ## Means are zeros as it is a correlation matrix
#     mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
# 
#     ## Corr=FALSE
#   } else {
#     ## M is a zero matrix
#     if (all(RAM1$M==0)) {
#       mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
#     } else {
#       M <- caracas::as_sym(RAM1$M)
#       mu <- M %*% t(F %*% invIA)
#       mu <- caracas::as_character_matrix(mu)
#     }
#   }
# 
#   ## Extract the parts with observed variables
#   Sigma <- F %*% SigmaAll %*% t(F)
# 
#   if (simplify) {Sigma <- caracas::simplify(Sigma)}
# 
#   ## Convert it into R matrix of string
#   Sigma <- caracas::as_character_matrix(Sigma)
#   Sigma <- .convert_reserved(Sigma, type="random2sympy")
#   dimnames(Sigma) <- list(rownames(RAM$F), rownames(RAM$F))
#   dimnames(mu) <- list("1", colnames(RAM$M)[colSums(RAM$F)==1])
# 
#   list(Sigma=Sigma, mu=mu, corr=corr)
# }
# 
# library(testthat)
# expect_equal(impliedS(RAM1, corr=FALSE), impliedS2(RAM1, corr=FALSE))
# expect_equal(impliedS(RAM1, corr=TRUE), impliedS2(RAM1, corr=TRUE))
# expect_equal(impliedS(RAM2, corr=FALSE), impliedS2(RAM2, corr=FALSE))
# expect_equal(impliedS(RAM2, corr=TRUE), impliedS2(RAM2, corr=TRUE))
# 
# library(microbenchmark)        ## Ad-hoc solution to replace reserved words in sympy

# microbenchmark(old=impliedS(RAM1, corr=TRUE), new=impliedS2(RAM1, corr=TRUE), times=50)
# microbenchmark(old=impliedS(RAM2, corr=TRUE), new=impliedS2(RAM2, corr=TRUE), times=50)




