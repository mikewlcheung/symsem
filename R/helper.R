## Convert math operators between R and yacas
.math.convert <- function(x, r2yacas=TRUE) {
    ## Conversion math functios from yac to R
    ## It is copied from yac-symbol.R Ryacas
    Math_transtab <- matrix( c(
        ## Need to replace asin before sin;
        ## otherwise, asin will be replaced by aSin
        ## R	    yacas
        "asin",	  	"ArcSin",
        "acos",	  	"ArcCos",
        "atan",    	"ArcTan",
        "asinh", 	"ArcSinh", 
        "acosh",    "ArcCosh", 
        "atanh",   	"ArcTanh",  
        "exp", 	  	"Exp",
        "log", 	  	"Ln",
        "sqrt", 	"Sqrt",
        "sin",		"Sin",
        "cos",		"Cos",
        "tan",		"Tan",
        ## Additional functions
        "abs",      "Abs",
        "sign",     "Sign"
    ), byrow = TRUE, ncol = 2)
    colnames(Math_transtab) <- c("R", "yacas")

    if (r2yacas) {
        j = 1; k = 2
    } else {
        j = 2; k = 1
    }

    ## Replace math functions from yac to R
    for (i in 1:nrow(Math_transtab)) {
        x <- gsub(Math_transtab[i, j], Math_transtab[i, k], x, fixed=TRUE)
    }
    x
}

## Convert yac to R
as.matrix.yac_symbol <- function(x, ...) {
    if (!inherits(x, "yac_symbol")) {
        stop("\"x\" is not a class of \"yac_symbol\".\n")
    }

    ## Remove all {} and white spaces
    z <- gsub("\\{|\\}|\\s", "", x$yacas_cmd)
    
    ## Convert yacas to R math operators
    z <- .math.convert(z, r2yacas=FALSE)

    ## Split the elements by ","
    z <- strsplit(z, ",")[[1]]

    x_dim <- dim(x)
    ## If no dimension, assume a column vector
    if (is.null(x_dim)) {
        x_dim <- c(length(z), 1)
    }
    
    out <- matrix(NA, nrow=x_dim[1], ncol=x_dim[2])
    k <- 1
    for (i in seq_len(x_dim[1]))
        for (j in seq_len(x_dim[2])) {
            out[i, j] <- z[k]
            k <- k+1
        }
    out
}

## Convert R scaloar, vectors, and matrices to Ryacs objects
sym <- function(x) {
    ## Need to convert the math operators for some math operators
    ## e.g., log(x) to Ln(x)
    x <- .math.convert(x, r2yacas=TRUE)
    Ryacas::ysym(x)
}
