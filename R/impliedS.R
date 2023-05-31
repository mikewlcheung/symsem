#' Compute a Symbolic Model-Implied Covariance/Correlation Matrix
#'
#' It computes a symbolic model-implied covariance (or correlation) matrix in
#' SEM using the RAM inputs.
#'
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}
#' @param corr Whether the model implied matrix is covariance (default) or
#' correlation structure.
#' @param simplify Attempt to simplify the output.
#' @return The model implied covariance (or correlation) matrix and means vector.
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
#' model2 <- "f =~ x1 + x2 + x3 + x4#'
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
impliedS <- function(RAM, corr=FALSE, simplify=TRUE) {
    ## Remove starting values and "*" in RAM
    RAM1 <- metaSEM::as.symMatrix(RAM)

    ## No. of observed and latent variables
    p <- nrow(RAM1$S)

    A <- caracas::as_sym(RAM1$A)
    S <- caracas::as_sym(RAM1$S)
    F <- caracas::as_sym(RAM1$F)

    ## Inverse(I-A)
    invIA <- inv(caracas::as_sym(diag(1L,p)) - A)
    ## Model implied covariance matrix with latent variables
    SigmaLat <- invIA %*% S %*% t(invIA)

    if (corr) {
        ## ## An alternative approach of standardization.
        ## ## Error variances are included as parameters.
        ## SD <- Ryacas::diag(SigmaLat)
        ## SD <- 1/(sqrt(SD))
        ## SD <- Ryacas::y_fn(SD, "DiagonalMatrix")
        ## SigmaLat <- SD * SigmaLat * SD

        ## F*S*t(F) rather than S
        ## Index of error variances in the diagonals of S
        index <- suppressWarnings(which(is.na(as.numeric(diag(RAM1$S)))))
        Slabels <- diag(RAM1$S)[index]

        ## representation in R characters
        Sigma.R <- as_character_matrix(SigmaLat)

        ## Assuming there are at least one free error variance in S including latent variables
        for (i in seq_along(index)) {
            j <- index[i]
            ## Sigma.R[i, i]: 1 + ErrVar - (sigma)
            ## It will become 1 on the diagonals the after simplification.
            ## For off-diagonal elements, the ErrVar are replaced by terms without ErrVar.
            Sigma.R <- gsub(Slabels[i],
                            paste0("(1+", Slabels[i], "-(", Sigma.R[i, i], "))"),
                            Sigma.R)
        }

        SigmaLat <- caracas::as_sym(Sigma.R)

        ## Means are zeros as it is a correlation matrix
        mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))

    ## Corr=FALSE
    } else {
        ## M is a zero matrix
        if (all(RAM1$M==0)) {
            mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
        } else {
            M <- caracas::as_sym(RAM1$M)
            mu <- M %*% t(F %*% invIA)
            mu <- caracas::as_character_matrix(mu)
        }
    }

    ## Extract the parts with observed variables
    Sigma <- F %*% SigmaLat %*% t(F)

    if (simplify) {Sigma <- caracas::simplify(Sigma)}

    ## Convert it into R matrix of string
    Sigma <- caracas::as_character_matrix(Sigma)
    dimnames(Sigma) <- list(rownames(RAM$F), rownames(RAM$F))
    dimnames(mu) <- list("1", colnames(RAM$M)[colSums(RAM$F)==1])

    list(Sigma=Sigma, mu=mu, corr=corr)
}
