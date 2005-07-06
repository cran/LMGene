"lnormeS" <-
function (eS, span=0.1)
{
normat <- lnorm (eS@exprs, span)
normed.eS <- neweS (normat, eS@phenoData@pData)
return(normed.eS)
}
