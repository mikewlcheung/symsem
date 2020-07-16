\name{as.matrix}
\alias{as.matrix}
\alias{as.matrix.yac_symbol}
\title{Convert a Ryacas Symbolic Matrix to a Character Matrix
}
\description{It converts a Ryacas Symbolic Matrix to a Character Matrix
  in R.
}
\usage{
\method{as.matrix}{yac_symbol}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A symbolic matrix of class "yac_symbol".}
  \item{\dots}{not used}
}
\value{A character matrix.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\seealso{\code{\link[symSEM]{sym}}
}

\examples{
## Create a character matrix
A1 <- matrix(c(1, 2, 3, "a", "b", "c"), ncol=2, nrow=3)
#      [,1] [,2]
# [1,] "1"  "a" 
# [2,] "2"  "b" 
# [3,] "3"  "c"

## Convert it to an Ryacas matrix for symbolic calculations
A2 <- sym(A1)
# {{1, a},
#  {2, b},
#  {3, c}} 

## Convert it back to a character matrix
as.matrix(A2)
#      [,1] [,2]
# [1,] "1"  "a" 
# [2,] "2"  "b" 
# [3,] "3"  "c" 
}

