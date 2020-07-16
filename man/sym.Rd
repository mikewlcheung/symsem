\name{sym}
\alias{sym}
\title{Convert a Character Expression to Ryacas Symbol or matrix
}
\description{It converts a Character Expression or Matrix in R to Ryacas
  symbolic expression or matrix.
}
\details{It converts mathematic operators from R to mathematic operators
  in yacas and then uses \code{\link[Ryacas]{ysym}} to convert the R character expression or
  matrix to Ryacas symbolic expression or matrix. The Ryacas expression
  or matrix can be used for further symbolic calculations.
 }
\usage{
sym(x)
}
\arguments{
  \item{x}{A character expression or matrix.}
}
\value{An expression or matrix of class "yac_symbol". 
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\seealso{ \code{\link[Ryacas]{ysym}}, \code{\link[symSEM]{as.matrix}}
}

\examples{
A1 <- "log(x) + exp(y)"

sym(A1)
# y: Ln(x)+Exp(y)

A2 <- matrix(c(1, 2, 3, "a", "sqrt(b)", "sin(c)"), ncol=2, nrow=3)
A2
#     [,1] [,2]     
# [1,] "1"  "a"      
# [2,] "2"  "sqrt(b)"
# [3,] "3"  "sin(c)" 

sym(A2)
# {{      1,       a},
#  {      2, Sqrt(b)},
#  {      3,  Sin(c)}} 
}
