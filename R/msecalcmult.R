"msecalcmult" <-
function (eS, lam, alpha, lowessnorm=FALSE, R, grads=TRUE) 
{
    lambda <- FALSE
    mat1 <- as.matrix(exprs(eS))
    n <- dim(mat1)[2]
    p <- dim(mat1)[1]
    y <- mat1
    r <- dim(y)[1]
    onevec <- matrix(1, nrow=r, ncol=1)

    #This does what jggrad does, but for msecalcmult
    #it has to be done internally. The reason is because
    #making all of the jga matrices would blow the RAM.
    #Instead, we calculate the alpha gradients sequentially.
    ya <- y - onevec%*%alpha
    g <- glog(ya, lam)
    z <- sqrt(ya^2 + lam)
    za <- -ya/z
    zl <- 1/(2 * z)
    ga.cat <- -1/z #Is sum of ga.i's, one per column
    gl <- 2 * z * (ya + z)
    gl <- 1/gl
    J <- exp(mean(log(z)))
    k <- length(alpha) #Will be used several times
    Ja <- J * colMeans(za/z) / k #Is a vector now, not a scalar
    Jl <- J * mean(zl/z)
    jg <- g * J
    jgl <- gl * J + g * Jl

    if (lowessnorm) {
        mat2l <- lnorm(jgl)
        #mat2a <- lnorm(mat2[, (2 * n + 1):(3 * n)])
        mat2 <- lnorm(jg)
    }
    else {
        mat2l <- norm(jgl)
        #mat2a <- norm(mat2[, (2 * n + 1):(3 * n)])
        mat2 <- norm(jg)
    }

    yres <- R %*% t(mat2)
    SSE <- sum(yres^2)

    #Rescaling test
    #SSE <- SSE / lam
    ##
    if (grads) {
      ylres <- R %*% t(mat2l)
      SSEl <- 2 * sum(yres * ylres)

      SSEa <- vector(length=k)
      for (i in 1:k) {
        ga.i <- matrix(0, nrow=r, ncol=k)
        ga.i[,i] <- ga.cat[,i]
        jga.i <- ga.i * J + g * Ja[i]
        if (lowessnorm) {mat2a <- lnorm(jga.i)}
        else {mat2a <- norm(jga.i)}
        yares.i <- R %*% t(mat2a)
        SSEa[i] <- 2 * sum(yres * yares.i)
      }

    #Scaling test
    #SSEl <- SSEl / lam - SSE / lam
    #SSEa <- SSEa / lam
    ##
    return(c(SSE, SSEl, SSEa))
  }
  else { return(SSE) }
}
