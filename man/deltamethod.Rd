\name{deltamethod}
\alias{deltamethod}
\title{Compute the Variance-Covariance Matrix of Functions using the Delta Method
}
\description{It computes the variance-covariance matrix of functions using the delta method.
}
\usage{
deltamethod(fn, Covvars, vars, Var.name="V", Cov.name="C", simplify=TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{A function in character strings or a vector of functions.}
  \item{Covvars}{Variance-covariance matrix of the variables.}
  \item{vars}{A vector of characters of random variables.}
  \item{Var.name}{Name of the variances.}
  \item{Cov.name}{Name of the covariances.}
  \item{simplify}{Attempts to simpify the expression}
}  

\value{Variance-covariance matrix of the functions. 
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
## Fisher-z-transformation
fn  <- "0.5*log((1+r)/(1-r))"

## Sampling variance of r
Covvars <- "(1-r^2)^2/n"

deltamethod(fn=fn, Covvars=Covvars, vars="r")
## $fn
##     [,1]                  
## fn1 "0.5*log((r+1)/(1-r))"

## $Covfn
##     fn1  
## fn1 "1/n"

## $vars
## [1] "r"

## $Covvars
##   r            
## r "(1-r^2)^2/n"

## $Jmatrix
##     r                                      
## fn1 "(0.5*(1-r+r+1)*(1-r))/((1-r)^2*(r+1))"

## Raw mean difference: y_treatment - y_control
fn <- "yt - yc"

## Sampling covariance matrix
## S2p: pooled variance
## nt: n_treatment
## nc: n_control
Covvars <- matrix(c("S2p/nt", 0,
                    0, "S2p/nc"),
                  ncol=2, nrow=2)

deltamethod(fn=fn, Covvars=Covvars, vars=c("yt", "yc"))
## $fn
##     [,1]   
## fn1 "yt-yc"

## $Covfn
##     fn1                      
## fn1 "(S2p*nt+S2p*nc)/(nt*nc)"

## $vars
## [1] "yt" "yc"

## $Covvars
##    yt       yc      
## yt "S2p/nt" "0"     
## yc "0"      "S2p/nc"

## $Jmatrix
##     yt  yc  
## fn1 "1" "-1"

## log(odds)
fn <- "log(p/(1-p))"

## Sampling variance of p
Covvars <- "p*(1-p)/n"

## Though it is correct, the simplification does not work well.
deltamethod(fn=fn, Covvars=Covvars, vars="p")
## $fn
##     [,1]          
## fn1 "log(p/(1-p))"

## $Covfn
##     fn1                                              
## fn1 "(3*p^2-p^3-3*p+1)/((p^4-4*p^3+6*p^2-4*p+1)*p*n)"

## $vars
## [1] "p"

## $Covvars
##   p            
## p "(p*(1-p))/n"

## $Jmatrix
##     p                            
## fn1 "((1-p+p)*(1-p))/((1-p)^2*p)"
}

