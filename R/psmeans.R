#Converts an ExpressionSet object with one row of expression 
#data per probeset into an ExpressionSet object with one row per 
#probe.  
psmeans <-
function (eS, ind)
{
  if (!is.element(class(eS), c("ExpressionSet", "AffyBatch"))){
 	stop("'eS' must be an object of class 'ExpressionSet' or 'AffyBatch'")
  }
  if (length(ind) != dim(exprs(eS))[1]){
	stop("'ind' must have length equal to the number of rows in 'exprs(eS)'")
  }
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
