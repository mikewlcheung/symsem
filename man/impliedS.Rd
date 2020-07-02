\name{impliedS}
\alias{impliedS}
\title{Compute a Symbolic Model Implied Covariance/Correlation Matrix
}
\description{It computes a symbolic model implied covariance (or
  correlation) matrices based on the RAM inputs.
}
\usage{
impliedS(RAM, corr=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{RAM}{A RAM object including a list of matrices of the model
	returned from \code{\link[metaSEM]{lavaan2RAM}}}
\item{corr}{Whether the model implied matrix is covariance (default) or correlation.}
}  

\value{The model implied covariance (or correlation) matrix. 
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
## A mediation model
model1 <- "y ~ c*x + b*m
           m ~ a*x
           ## Means
           y ~ b0*1
           m ~ m0*1
           x ~ x0*1"
    
RAM1 <- metaSEM::lavaan2RAM(model1)

## Implied covariance matrix and mean structure
impliedS(RAM1, corr=FALSE)

## Implied correlation matrix
impliedS(RAM1, corr=TRUE)

## A CFA model
model2 <- "f =~ x1 + x2 + x3 + x4"
    
RAM2 <- metaSEM::lavaan2RAM(model2)

## Implied covariance matrix
impliedS(RAM2, corr=FALSE)

## Implied correlation matrix
impliedS(RAM2, corr=TRUE)
}

