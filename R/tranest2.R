"tranest2" <-
function(eS,starting=FALSE,lambda=1000,alpha=0,gradtol=1e-3,lowessnorm)
{
  starttime <- Sys.time()
  mat1 <- as.matrix(eS@exprs)
  
  like <- function(x)
  {
    lnlam <- x[1]
    alpha <- x[2]
    lam <- exp(lnlam)
print(c("lam lnlam alpha", format(lam, nsmall=2), format(lnlam, nsmall=2), format(alpha, nsmall=2)))
    tmp <- msecalc(eS,lam,alpha,lowessnorm)
    tmpv <- tmp[1]

    attr(tmpv,"gradient") <- c(tmp[2]*lam,tmp[3])
#    attr(tmpv,"gradient") <- c(tmp[2],tmp[3])

print(c("tmp", format(tmp[1:3],nsmall=2), lam))
print(c("tmpv",format(tmpv, nsmall=2), format(attr(tmpv,"gradient"), nsmall=2)) )
    return(tmpv)
  }

  if (starting)
  {
    lamstart <- log(lambda)
    alphastart <- alpha
  }
  else{
    lamstart <- log(median(abs(mat1))^2)
    alphastart <- quantile(abs(as.vector(mat1)),.1)
  }
  
  typsize <- c(lamstart,alphastart)
  startvar <- c(lamstart,alphastart)
  stepmax <- 3
  
  #opt <- nlm(like,startvar,stepmax=4,typsize=typsize,check.analyticals=FALSE,gradtol=gradtol, print.level = 2)
  #print(difftime(Sys.time(),starttime))
  #print(opt)
  #return(opt$estimate)
  
  ##
  return(beams(eS, startvar))
  ##
  
  
}

