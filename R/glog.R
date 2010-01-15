#This function transforms the input values by the generalized #log function.
glog <-
function(y,lambda)
{
  yt <- log(y+sqrt(y^2+lambda))
  return(yt)
}
