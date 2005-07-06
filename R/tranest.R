"tranest" <-
function (eS, ngenes = -1, starting = FALSE, lambda = 1000, alpha = 0,
    gradtol = 1e-3, lowessnorm = FALSE, method=1, mult=FALSE, model=NULL)
{
    mat1 <- as.matrix(exprs(eS))
    n <- dim(mat1)[2]
    p <- dim(mat1)[1]
    if (p > 100000) {ngenes <- 50000}
    if ((ngenes < p) & (ngenes > 0))
        mat2 <- mat1[sample(p, ngenes), ]
    else mat2 <- mat1
    eS2 <- new("ExpressionSet", exprs = mat2, phenoData = phenoData(eS))

    if (mult==FALSE) {
      tranpar <- tranest2(eS2, starting, lambda, alpha, gradtol,
        lowessnorm, method, model)
      return(list(lambda = (tranpar[1]), alpha = tranpar[2:length(tranpar)]))
    }
    else {
      return(tranestmult(eS2, starting, lambda, alpha, gradtol,
        lowessnorm, method, 200, model))
    }
}
