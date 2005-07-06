"normeS" <-
function (eS)
{
normat <- norm (eS@exprs)
normed.eS <- neweS (normat, eS@phenoData@pData)
return(normed.eS)
}
