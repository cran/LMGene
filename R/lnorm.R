"lnorm" <-
function(mat1,span=.1)
{
  mat2 <- as.matrix(mat1)
  p <- dim(mat2)[1]
  n <- dim(mat2)[2]
  rmeans <- apply(mat2,1,mean)
  rranks <- rank(rmeans,ties.method="first")
  matsort <- mat2[order(rranks),]  
  r0 <- 1:p
  lcol <- function(x)
  {
    lx <- lowess(r0,x,f=span)$y
  }
  lmeans <- apply(matsort,2,lcol)
  lgrand <- apply(lmeans,1,mean)
  lgrand <- matrix(rep(lgrand,n),byrow=FALSE,ncol=n)
  matnorm0 <- matsort-lmeans+lgrand
  matnorm1 <- matnorm0[rranks,]
  return(matnorm1)
}
