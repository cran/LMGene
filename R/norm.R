"norm" <-
function(mat1)
{
  mat2 <- as.matrix(mat1)
  p <- dim(mat2)[1]
  n <- dim(mat2)[2]
  cmean <- apply(mat2,2,mean)
  cmean <- cmean - mean(cmean)
  mnmat <- matrix(rep(cmean,p),byrow=TRUE,ncol=n)
  return(mat2-mnmat)
}
