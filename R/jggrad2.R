"jggrad2" <-
function(y,lambda,alpha)
{
  ya <- y-alpha
  g <- glog(ya,lambda)
  z <- sqrt(ya^2+lambda)
  za <- -ya/z
  zl <- 1/(2*z)
  ga <- -1/z
  gl <- 2*z*(ya+z)
  gl <- 1/gl
  J <- exp(mean(log(z)))
  Ja <- J*mean(za/z)
  Jl <- J*mean(zl/z)
  jg <- g*J
  jga <- ga*J+g*Ja
  jgl <- gl*J+g*Jl
  return(cbind(jg,jgl,jga))
}

