tranest2 <-
function (eS, starting = FALSE, lambda = 1000, alpha = 0, gradtol = 0.001, 
    lowessnorm=FALSE, method=1, model=NULL, SD = FALSE, rank = TRUE, model.based = TRUE, rep.arrays = NULL) 
{
    control <- NULL
    starttime <- Sys.time()
    mat1 <- as.matrix(exprs(eS))
    if (starting) {
        lamstart <- log(lambda)
        alphastart <- alpha
    }
    else {
        lamstart <- log(median(abs(mat1))^2)
        alphastart <- quantile(abs(as.vector(mat1)), 0.1, names = FALSE)
    }
    typsize <- c(lamstart, alphastart)
    startvar <- c(lamstart, alphastart)
    names(startvar) <- NULL
    stepmax <- 3

    #Calculate I minus hat matrix 
    if (!(SD == TRUE && model.based == FALSE)){
    	mod <- GetLMObj(eS, model)
    	X <- mod$x
    	df <- mod$df * dim(mat1)[1]
    	U <- svd(X)$u
    	H <- U %*% t(U)
    	n <- dim(H)[1]
    	R <- diag(rep(1,n)) - H
    } else{
	R <- NULL
    }

    if (SD == TRUE){
	  startvar <- c(exp(lamstart), alphastart)
	  control <- list(reltol = gradtol)
        stability <- function(x){
	  	lambda <- x[1]
	  	alpha <- x[2]
	  	tmp <- arrayGlogSDStability(eS,lambda,alpha,lowessnorm,rank,model.based,R,rep.arrays)
        	return(tmp)
        }
        objfun <- stability
    } else{
    	 like <- function(x) {
        	lnlam <- x[1]
        	alpha <- x[2]
        	lam <- exp(lnlam)
        	tmp <- msecalc(eS, lam, alpha, lowessnorm, R)
        	tmpv <- tmp[1] / df
        	attr(tmpv, "gradient") <- c(tmp[2] * lam, tmp[3]) / df
        	return(tmpv)
    	  }
	  objfun <- like
    }

    if (method==1) {
       opt <- nlm(objfun, startvar, stepmax=4, typsize=typsize,
         check.analyticals=FALSE, gradtol=gradtol)
	if (SD ==TRUE){
		return(c(opt$estimate[1], opt$estimate[2]))
	} else{
       	return(c(exp(opt$estimate[1]), opt$estimate[2]))
	}
    } else {
       if (method==2) {optype <- 'Nelder-Mead'}
       if (method==3) {optype <- 'BFGS'}
       if (method==4) {optype <- 'CG'}
       opt <- optim (startvar, objfun, method=optype, control = control) 
	 if (SD == TRUE){
		return(c(opt$par[1], opt$par[2]))
	 } else{
       	return(c(exp(opt$par[1]), opt$par[2]))
	}
    }
}
