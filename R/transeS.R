transeS <-
function (eS, lambda, alpha) 
{
  if (!is.element(class(eS),c("ExpressionSet","AffyBatch"))){
	stop("'eS' must be an object of class 'ExpressionSet' or 'AffyBatch'")
  }
  mat <- exprs(eS)
  if (length(alpha)==1) {
    mat.cor <- mat-alpha
  } else {
    r <- dim(mat)[1]
    onevec <- matrix(1, nrow=r, ncol=1)
    mat.cor <- mat - onevec %*% alpha
  }

  mat.trans <- glog(mat.cor, lambda)
  eS.trans <- eS
  exprs(eS.trans) <- mat.trans
  return(eS.trans)
}
