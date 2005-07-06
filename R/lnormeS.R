"lnormeS" <-
function (eS, span=0.1)
{
  normat <- lnorm(exprs(eS), span)
  normed.eS <- neweS(normat, pData(eS))
  return(normed.eS)
}
