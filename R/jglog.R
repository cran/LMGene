"jglog" <-
function(y,lambda)
{
  z <- sqrt(y^2+lambda)
  gmn <- exp(mean(log(z)))
  y1 <- glog(y,lambda)
  y1 <- y1*gmn
  return(y1)
}

