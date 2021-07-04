## TODO: it is slow as gsub has to be applied multiple times to fix " _" to "_".
impliedS <- function(RAM, corr=FALSE) {

    ## Remove starting values and "*" in RAM
    RAM1 <- metaSEM::as.symMatrix(RAM)
           
    ## No. of observed and latent variables
    p <- nrow(RAM1$S) 

    ## Pass them to Ryacas
    A <- sym(RAM1$A)
    S <- sym(RAM1$S)
    F <- sym(RAM1$F)

    ## Inverse(I-A)
    invIA <- solve(diag(1,p) - A)
    ## Model implied covariance matrix with latent variables
    SigmaLat <- invIA * S * t(invIA)    
    SigmaLat$yacas_cmd <- gsub("\\s", "", SigmaLat$yacas_cmd)
    
    if (corr) {

        ## ## An alternative approach of standaridzation.
        ## ## Error variances are included as parameters.
        ## SD <- Ryacas::diag(SigmaLat)
        ## SD <- 1/(sqrt(SD))
        ## SD <- Ryacas::y_fn(SD, "DiagonalMatrix")
        ## SigmaLat <- SD * SigmaLat * SD
     
        ## F*S*t(F) rather than S
        ## Index of error variances in the diagonals of S
        index <- suppressWarnings(which(is.na(as.numeric(diag(RAM1$S)))))
        Slabels <- diag(RAM1$S)[index]

        ## Assuming there are at least one free error variance in S including latent variables
        for (i in seq_along(index)) {
            j <- index[i]
            ## sigmalat: diagonals including Error variance in S with observed and latent variables. It should be constrained as 1.
            sigmalat <- paste0("Ryacas::yac_str(SigmaLat[",j,",",j,"])")
            sigmalat <- eval(parse(text=sigmalat))
            sigmalat <- gsub("\\s", "", sigmalat)
            ## 1 + ErrVar - (sigma): It will become 1 on the diagonals the after simplification.
            ## For off-diagonal elements, the ErrVar are replaced by terms without ErrVar.
            SigmaLat$yacas_cmd <- gsub(Slabels[i],
                                       paste0("(1+", Slabels[i], "-(", sigmalat, "))"),
                                       SigmaLat$yacas_cmd)
        }
        
        SigmaLat <- Ryacas::simplify(SigmaLat)
        SigmaLat$yacas_cmd <- gsub("\\s", "", SigmaLat$yacas_cmd)
        mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
        
    } else {

        ## M is a zero matrix
        if (all(RAM1$M==0)) {
            mu <- matrix(0, nrow=1, ncol=nrow(RAM$F))
        } else {
            M <- sym(RAM1$M)
            M$yacas_cmd <- gsub("\\s", "", M$yacas_cmd)
            mu <- M * t(F*invIA)
            mu$yacas_cmd <- gsub("\\s", "", mu$yacas_cmd)
        }
    }

    Sigma <- F*SigmaLat*t(F)
    Sigma$yacas_cmd <- gsub("\\s", "", Sigma$yacas_cmd)
    
    ## Convert it into R matrix of string
    Sigma <- as.matrix(Sigma)
    dimnames(Sigma) <- list(rownames(RAM$F), rownames(RAM$F))
    mu <- as.matrix(mu)
    dimnames(mu) <- list("1", colnames(RAM$M)[colSums(RAM$F)==1])
        
    list(Sigma=Sigma, mu=mu, corr=corr)
}
