"msa" <-
function(v)
{
  tmp <- sum(v^2)
  if (tmp==0)tmp <- 1
  return(v/sum(v))
  
}
