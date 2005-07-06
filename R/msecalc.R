"msecalc" <- 
function (eS, lam, alpha, lowessnorm, R) 
{
    starttime <- Sys.time()
    mat1 <- as.matrix(eS@exprs)
 #   for (i in 1:length(eS@phenoData@varLabels)) {
 #       assign(paste("x", i, sep = ""), as.factor(eS@phenoData@pData[, 
 #           i]))
 #   }
 #   fchar = ""
 #   for (i in 1:length(eS@phenoData@varLabels)) {
 #       fchar = paste(fchar, paste("x", i, sep = ""), ifelse(i < 
 #           length(eS@phenoData@varLabels), "+", ""), sep = "")
 #   }
    mat2 <- as.matrix(mat1)

    n <- dim(mat2)[2]
    p <- dim(mat2)[1]
    mat2 <- jggrad2(mat2, lam, alpha)
    if (lowessnorm) {
        mat2l <- lnorm(mat2[, (n + 1):(2 * n)])
        mat2a <- lnorm(mat2[, (2 * n + 1):(3 * n)])
        mat2 <- lnorm(mat2[, 1:n])
    }
    else {
        mat2l <- norm(mat2[, (n + 1):(2 * n)])
        mat2a <- norm(mat2[, (2 * n + 1):(3 * n)])
        mat2 <- norm(mat2[, 1:n])
    }

    yres = R %*% t(mat2)
    ylres = R %*% t(mat2l)
    yares = R %*% t(mat2a)
    r1 <- sum(yres^2)
    r2 <- sum(2 * yres * ylres)
    r3 <- sum(2 * yres * yares)

    return(c(r1, r2, r3))
    #return(r1)
}
