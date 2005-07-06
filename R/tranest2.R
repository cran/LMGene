"tranest2" <-
function (eS, starting = FALSE, lambda = 1000, alpha = 0, gradtol = 0.001, 
    lowessnorm=FALSE, method=1, model=NULL) 
{
    starttime <- Sys.time()
    mat1 <- as.matrix(exprs(eS))
    like <- function(x) {
        lnlam <- x[1]
        alpha <- x[2]
        lam <- exp(lnlam)
        tmp <- msecalc(eS, lam, alpha, lowessnorm, R)
        tmpv <- tmp[1] / df
        attr(tmpv, "gradient") <- c(tmp[2] * lam, tmp[3]) / df
        return(tmpv)
    }
    if (starting) {
        lamstart <- log(lambda)
        alphastart <- alpha
    }
    else {
        lamstart <- log(median(abs(mat1))^2)
        alphastart <- quantile(abs(as.vector(mat1)), 0.1)
    }
    typsize <- c(lamstart, alphastart)
    startvar <- c(lamstart, alphastart)
    names(startvar) <- NULL
    stepmax <- 3

    #Calculate residuals matrix R (from "hat matrix")
    mod <- GetLMObj(eS, model)
    X <- mod$x
    df <- mod$df * dim(mat1)[1]
    U <- svd(X)$u
    H <- U %*% t(U)
    n <- dim(H)[1]
    R <- diag(rep(1,n)) - H

    if (method==1) {
       opt <- nlm(like, startvar, stepmax=4, typsize=typsize,
         check.analyticals=FALSE, gradtol=gradtol)
       return(c(exp(opt$estimate[1]), opt$estimate[2]))
    } else {
       if (method==2) {optype <- 'Nelder-Mead'}
       if (method==3) {optype <- 'BFGS'}
       if (method==4) {optype <- 'CG'}
       opt <- optim (startvar, like, method=optype) 
       return(c(exp(opt$par[1]), opt$par[2]))
    }
}
