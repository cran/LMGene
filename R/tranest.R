"tranest" <-
function(eS,ngenes=-1,starting=FALSE,lambda=1000,alpha=0,gradtol=1e-3,lowessnorm=FALSE)
{
  mat1 <- as.matrix(eS@exprs)
  
  n <- dim(mat1)[2]
  p <- dim(mat1)[1]
  
  if((ngenes < p)&(ngenes>0)) mat2=mat1[sample(p,ngenes),]
  else mat2=mat1
  
  eS2=new("exprSet", exprs=mat2, phenoData=eS@phenoData)
  tranpar <- tranest2(eS2,starting,lambda,alpha,gradtol,lowessnorm)
  #return(list(lambda=exp(tranpar[1]),alpha=tranpar[2]))
  return(list(lambda=(tranpar[1]),alpha=tranpar[2]))
}

