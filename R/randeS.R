"randeS" <-
function (eS, n) 
{
output = eS
M = eS@exprs
rows = dim(M)[1]
q = sample(1:rows, n)
output@exprs = M[q,]
return(output)
}
