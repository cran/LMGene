#Like lnorm, but applies to and returns an ExpressionSet object #instead of a matrix
lnormeS <-
function (eS, span=0.1)
{
  if (!is.element(class(eS), c("ExpressionSet", "AffyBatch"))){
 	stop("'eS' must be an object of class 'ExpressionSet' or 'AffyBatch'")
  }
  normat <- lnorm(exprs(eS), span)
  normed.eS <- eS
  exprs(normed.eS) <- normat
  return(normed.eS)
}
