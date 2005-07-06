"psmeans" <-
function (eS, ind)
{
r <- dim(eS@exprs)[1]
c <- dim(eS@exprs)[2]

I = max(ind)
outmat <- matrix(0, nrow=I, ncol=c)
for (i in 1:I) {
  k = sum(ind==i)
  v = rep(1/k, k)
  outmat[i,] <- v %*% eS@exprs[ind==i,]
}

meaneS <- neweS (outmat, eS@phenoData@pData)
return(meaneS)
}
