"arrplotd" <-
function(eS,cs="")
{
#  if (high==-1) high=dim(slides)[1]
#  tmp <- rank(apply(slides,1,mean))
#  ind <- (tmp >= low) & (tmp <= high)
#  slides <- slides[ind,]

  tmp <- arrplot1(eS)
  dplot(tmp,cs)
#  return(tmp)
}

