tranestmult <-
function (eS, starting = FALSE, lambda = 1000, alpha = 0, gradtol = 0.001,
    lowessnorm=FALSE, method=1, max_iter=200, model=NULL)
{
    starttime <- Sys.time()
    mat1 <- as.matrix(exprs(eS))
    like <- function(x) {
        lnlam <- x[1]
        alpha <- x[2:length(x)]
        lam <- exp(lnlam)

        tmp <- msecalcmult(eS, lam, alpha, lowessnorm, R, grads=grads)
        tmpv <- tmp[1] / df
        attr(tmpv, "gradient") <- c(tmp[2] * lam, tmp[3:length(tmp)]) / df

        return(tmpv)
    }
    k <- dim(mat1)[2]
    p <- dim(mat1)[1]
    if (starting) {
        lamstart <- log(lambda)
        if (length(alpha)==1) {alphastart <- rep(alpha, k)}
        else {alphastart <- alpha}
    }
    else {
        lamstart <- log(median(abs(mat1))^2)
        lamstart <- log(1000)
        alphastart <- rep(0, k)
        for (i in 1:k) {
          #alphastart[i] = quantile(abs(as.vector(mat1[,i])), 0.1)
          alphastart[i] = min(abs(as.vector(mat1[,i])))
        }
    }
    typsize <- c(lamstart, alphastart)
    startvar <- c(lamstart, alphastart)
    names(startvar) <- NULL
    stepmax <- 3

    #Calculate residuals matrix R (from "hat matrix")
    mod <- GetLMObj(eS, model)
    X <- mod$x
    df <- p * mod$df

    U <- svd(X)$u
    H <- U %*% t(U)
    n <- dim(H)[1]
    R <- diag(rep(1,n)) - H

    if (method==1) {
       grads <- TRUE
       opt <- nlm(like, startvar, stepmax=4, typsize=typsize,
         check.analyticals=FALSE, gradtol=gradtol, steptol=1e-8, iterlim=max_iter)
       return(list(lambda=exp(opt$estimate[1]), alpha=opt$estimate[2:length(opt$estimate)]))
    }
    else {
       if (method==2) {grads <- FALSE; optype <- 'Nelder-Mead'}
       if (method==3) {grads <- TRUE; optype <- 'BFGS'}
       if (method==4) {grads <- TRUE; optype <- 'CG'}
       if (method==5) {grads <- FALSE; optype <- 'SANN'}
       opt <- optim(startvar, like, method=optype, control=list(maxit=max_iter))
       print (opt)
       return(list(lambda=exp(opt$par[1]), alpha=opt$par[2:length(opt$par)]))
    }
}
