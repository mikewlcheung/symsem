.convert <- function(x, type=c("sympy2random", "random2sympy")) {
  x.old <- x
  
  ## Conversion the reserved words in sympy with some random strings
  words <- matrix( c(
    ## random	    sympy
    "a4dbepz",	  	"I",
    "ajljael",	  	"lambda",
    "akr0dsf",    	"E",
    "data_",        "data.",
    "jeljsxe",      "oo"
  ), byrow = TRUE, ncol = 2)
  colnames(words) <- c("random", "sympy")
  
  type <- match.arg(type)
  switch (type,
          sympy2random = for (i in 1:nrow(words)) {
            x <- gsub(words[i, 2], words[i, 1], x, fixed=TRUE)},
          random2sympy = for (i in 1:nrow(words)) {
            x <- gsub(words[i, 1], words[i, 2], x, fixed=TRUE)}
  )
  
  ## Try to convert these matrices from string to numeric, if possible
  if (suppressWarnings(all(!is.na(as.numeric(x))))) {
    x <- as.numeric(x)
  }
  
  if (is.matrix(x.old)) {
    x <- matrix(x, nrow=nrow(x.old), ncol=ncol(x.old),
                dimnames=dimnames(x.old))
  } else if (is.matrix(x)) {
    names(x) <- names(x.old)
  }
  x
}
