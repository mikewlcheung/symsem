## fn <- "log(a)^b+c"
## fn <- c("x^2 + log(y)*z", "x+y+sin(z)")
## Note. Initiate symbols first does not work as there are problems in creating 
## Better: replace R -> yac, create character matrices, ysym convert

## TODO: ensure dimensions of Covvars and vars are matched
deltamethod <- function(fn, Covvars, vars, Var.name="V", Cov.name="C", simplify=TRUE) {
    ## Univariate or multivariate
    fn.p <- length(fn)
    
    ## function names
    fn.names <- paste0("fn", seq_len(fn.p))

    ## Convert it into a column vector if multivariate
    if (fn.p>1) {
        fn <- matrix(fn, ncol=1)
    }
    
    ## Get the variable names
    varlist <- sort(unique(all.names(parse(text=fn), functions=FALSE)))
    if (missing(vars)) {
        vars <- varlist
    } else {
        if (any(!vars %in% varlist)) {
            stop("Some of \"vars\" do not agree with the names in \"fn\".\n")
        }        
    }

    ## Need to convert all variables into symbols first;
    ## otherwise, some functions, e.g., Ln, will not be converted from R to yac.
    #for (i in seq_along(varlist)) {
    #    eval(parse(text=paste0(varlist[i], " <- sym('", varlist[i],"')")))
    #}    
    #fn <- eval(parse(text=fn))
    fn <- sym(fn)

    ##
    if (missing(Covvars)) {
    
        ## Variance covariance matrix of x
        Covvars <- outer(vars, vars, function(y, z) paste0(Cov.name, y, z))
        metaSEM::Diag(Covvars) <- paste0(Var.name, vars)
        ## Make it symmetric
        Covvars <- metaSEM::vec2symMat(OpenMx::vech(Covvars))
        Covvars <- sym(Covvars)
    } else {
        Covvars <- sym(Covvars)
    }
 
    ## Jacobian matrix
    ## Univariate
    if (fn.p==1) {
        J <- deriv(fn, vars)
        ## Possibly bugs in Ryacas?
        ## Need to clean up as too many {}
        J <- sym(as.matrix(J))
        J <- t(J)
    } else {
    ## Multivariate
        J <- Ryacas::Jacobian(fn, vars)
        J <- sym(as.matrix(J))
        # J <- Ryacas::simplify(J)
    }

    fn <- as.matrix(fn)
    rownames(fn) <- fn.names
    
    Covfn <- J * Covvars * t(J)

    if (simplify) {
        Covfn <- Ryacas::simplify(Covfn)
        Jmatrix <- Ryacas::simplify(J)
    }

    Covfn <- as.matrix(Covfn)
    dimnames(Covfn) <- list(fn.names, fn.names)
    
    Covvars <- as.matrix(Covvars)
    dimnames(Covvars) <- list(vars, vars)

    Jmatrix <- as.matrix(J)
    dimnames(Jmatrix) <- list(fn.names, vars)
    
    list(fn=fn, Covfn=Covfn, vars=vars, Covvars=Covvars, Jmatrix=Jmatrix)
}
