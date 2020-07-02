impliedS <- function(RAM, corr=FALSE) {

    ## Remove starting values and "*" in RAM
    RAM1 <- metaSEM::as.symMatrix(RAM)
           
    ## No. of observed and latent variables
    p <- nrow(RAM1$S) 

    ## Pass them to Ryacas
    A <- sym(RAM1$A)
    S <- sym(RAM1$S)
    F <- sym(RAM1$F)

    ## F %*% Inverse(I-A)
    FinvIA <- F * solve(diag(1,p) - A)
    ## Model implied covariance matrix
    Sigma <- FinvIA * S * t(FinvIA)

    if (corr) {

        ## Index of error variances in the diagonals of S
        index <- suppressWarnings(which(is.na(as.numeric(diag(RAM1$S)))))
        Slabels <- diag(RAM1$S)[index]

        if (any(sapply(Slabels, function(x) grepl("_", x, fixed=TRUE)))) {
            warning("The character \"_\" is found in the variable labels. The results are likely incorrect. Please change the variable labels and rerun it.\n")
        }
            
        ## Assuming there are at least one free error variance in S
        for (i in seq_along(index)) {
            j <- index[i]
            ## sigma: diagonal including Error variance in S. It should be constrained as 1.
            sigma <- paste0("Ryacas::yac_str(Sigma[",j,",",j,"])")
            sigma <- eval(parse(text=sigma))
            ## 1 + Err - (sigma): It will be 1 after simplification
            Sigma$yacas_cmd <- gsub(Slabels[i],
                                    paste0("(1+", Slabels[i], "-(", sigma, "))"),
                                    Sigma$yacas_cmd)
        }

        Sigma <- Ryacas::simplify(Sigma)
        mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
        
    } else {

        ## M is a zero matrix
        if (all(RAM1$M==0)) {
            mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
        } else {
            M <- sym(RAM1$M)
            mu <- M * t(FinvIA) 
            mu <- as.matrix(mu)
        }
    }

    ## Convert it into R matrix of string
    Sigma <- as.matrix(Sigma)
    dimnames(Sigma) <- list(rownames(RAM$F), rownames(RAM$F))
    dimnames(mu) <- list("1", colnames(RAM$M)[colSums(RAM$F)==1])

    list(Sigma=Sigma, mu=mu, corr=corr)
}
