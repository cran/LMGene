"msecalc" <-
function(eS,lam,alpha,lowessnorm)
#
# computes the mean square error and gradient for the global anova
#   Assumes a fixed effects model
# mat1 is a p by n matrix of expression values
# vlist is a list of the form list(v1=v1,v2=v2,..) of variables, 
#   each of length n, the number of chips
# fchar is the model for each gene in the form of a character 
#   string with variables from vlist. An example is "v1 + v2"
#
# returns the mse and gradient
#
# create variables
#
{ 
  starttime <- Sys.time()

  mat1 <- as.matrix(eS@exprs)
  
  # model information 
  for(i in 1:length(eS@phenoData@varLabels)){
  assign(paste('x', i, sep=''),as.factor(eS@phenoData@pData[,i]))
  }
  
  fchar=''
  for(i in 1:length(eS@phenoData@varLabels)){
  fchar=paste(fchar, paste('x', i, sep=''), ifelse(i<length(eS@phenoData@varLabels), '+', ''), sep='')
  }

  mat2 <- as.matrix(mat1)
  n <- dim(mat2)[2]
  p <- dim(mat2)[1]
  mat2 <- jggrad2(mat2,lam,alpha)
  if (lowessnorm)
  {
    mat2l <- lnorm(mat2[,(n+1):(2*n)])
    mat2a <- lnorm(mat2[,(2*n+1):(3*n)])
    mat2 <- lnorm(mat2[,1:n])
  }
  else
  {
    mat2l <- norm(mat2[,(n+1):(2*n)])
    mat2a <- norm(mat2[,(2*n+1):(3*n)])
    mat2 <- norm(mat2[,1:n])
  }
  #nvar <- length(vlist)
  #for (ivar in 1:nvar)
  #{
  #  assign(names(vlist)[ivar],vlist[[ivar]])
  #}


  fchar1 <- paste("y ~",fchar)
  y <- t(mat2)
  formobj <- as.formula(fchar1)
  lm1 <- lm(formobj)
  y <- t(mat2l)
  lm2 <- lm(formobj)
  y <- t(mat2a)
  lm3 <- lm(formobj)
  df <- p*lm1$df
  r1 <- sum( (lm1$resid)^2)/df
  r2 <- sum( 2*(lm1$resid)*(lm2$resid))/df
  r3 <- sum( 2*(lm1$resid)*(lm3$resid))/df
#  print(difftime(Sys.time(),starttime))
  return(c(r1,r2,r3))
}

