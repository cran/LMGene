#Computes two sets of p-values per gene or probe
#using gene-by-gene ANOVA with individual gene MSE using 
#both the gene-specific MSE and the posterior mean MSE for
#each term in the ANOVA. P-values are not adjusted for multiple #testing.
genediff <- function (eS, model = NULL, method = c("MLE", "MOM", "MOMlog"), 
verbose = TRUE){
    method <- match.arg(method)
    if (class(eS) != "ExpressionSet"){
	stop("'eS' must be an object of class 'ExpressionSet'")
    }
    if (!is.element(method, c("MLE", "MOM", "MOMlog"))){
	stop("'method' must be one of 'MLE', 'MOM', 'MOMlog'")
    }
    

    mat1 <- as.matrix(exprs(eS))
    if (is.null(model)) {
        model <- ""
        vars <- varLabels(eS)
        for (i in 1:length(vars)) {
            model <- paste(model, vars[i], ifelse(i < length(vars), 
                "+", ""), sep = "")
        }
    }
    model2 <- paste("y ~", model)
    mat2 <- as.matrix(mat1)
    p <- dim(mat2)[1]
    n <- dim(mat2)[2]
    owaov <- function(y) {
        formobj <- as.formula(model2)
        tmp <- row.names(anova(lm(formobj, data = pData(eS))))
        return(tmp)
    }
    effnames <- owaov(mat2[1, ])
    numeff <- length(effnames)
    tmp1 <- rowaov(eS, model)
    if (is.null(tmp1)) {
        return(NULL)
    }
    msmat <- tmp1[1:numeff, ]
    dfmat <- tmp1[(numeff + 1):(2 * numeff), ]
    numf <- numeff - 1
    nu <- median(dfmat[numeff, ])
    dfvec <- dfmat[numeff, ]
    msevec <- msmat[numeff, ]
    mn <- mean(msevec)
    v <- var(msevec)
    alpha.MOM <- (mn^2 - 2 * mn^2/nu + 2 * v)/(v - 2 * mn^2/nu)
    beta.MOM <- 1/(mn * (alpha.MOM - 1))
    if (method == "MOM"){
	alpha <- alpha.MOM
    	eta <- alpha.MOM * beta.MOM
    }
    if (method == "MLE"){
	if (alpha.MOM > 0 && beta.MOM > 0){
		startvals <- c(alpha.MOM, beta.MOM)
	} else{
		startvals <- c(1,1)
	}
	parests <- fitScaledF(msevec, nu, startvals)
	alpha <- parests[1]
	beta <- parests[2]
	eta <- alpha*beta
    }
    if (method == "MOMlog"){
	parests <- fitFDist(msevec, nu)
	s20 <- parests$scale
	df0 <- parests$df2
	alpha <- df0/2
	eta <- 1/s20
    }
    prior <- 1/eta
    df <- 2 * alpha
    if (df < 0 || df == Inf){
	adjmsevec <- rep(1/eta, length(msevec))
	adjdfvec <- rep(Inf, length(dfvec))
	warning('negative or infinite estimate for prior d.f.,
		posterior estimates set to inverse of prior mean,
			exercise caution in using posterior estimates')
    } else{
    	adjdfvec <- dfvec + df
    	adjmsevec <- (msevec * dfvec + prior * df)/adjdfvec
    }
    pmat1 <- matrix(numeric(numf * p), ncol = numf)
    pmat2 <- pmat1
    for (i in 1:numf) {
        pmat1[, i] <- 1 - pf(msmat[i, ]/msevec, dfmat[i, ], dfvec)
        pmat2[, i] <- 1 - pf(msmat[i, ]/adjmsevec, dfmat[i, ], 
            adjdfvec)
    }
    colnames(pmat1) <- effnames[-numeff]
    colnames(pmat2) <- effnames[-numeff]
    if (verbose == TRUE){
    	cat("Prior d.f. = ", df, "\n")
    	cat("Prior mean reciprocal precision = ", prior, "\n")
    }
    return(list(Gene.Specific = pmat1, Posterior = pmat2))
}
