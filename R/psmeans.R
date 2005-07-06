"psmeans" <-
function (eS, ind)
{
  r <- dim(exprs(eS))[1]
  c <- dim(exprs(eS))[2]

  I <- max(ind)
  outmat <- matrix(0, nrow=I, ncol=c)
  for (i in 1:I) {
    k <- sum(ind==i)
    v <- rep(1/k, k)
    outmat[i,] <- v %*% exprs(eS)[ind==i,]
  }

  meaneS <- neweS(outmat, pData(eS))
  return(meaneS)
}
