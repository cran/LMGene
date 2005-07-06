"genediff" <-
function(eS, model=NULL)
#
# computes two vectors of p-values per gene or probe
#   using gene-by-gene anova with individual gene MSE using 
#   both the gene-specific MSE and the posterior mean MSE for
#   each term in the anova. Assumes a fixed effects model and
#   the correct denominator for all comparisons is the MSE
# mat1 is a p by n matrix of expression values
# vlist is a list of the form list(v1=v1,v2=v2,..) of variables, 
#   each of length n, the number of chips
# model is the model for each gene in the form of a character 
#   string with variables from vlist. An example is "v1 + v2"
#
# create variables
#
{ 
  mat1 <- as.matrix(eS@exprs)
  
  if (is.null(model)) {
    model=''
    vars <- eS@phenoData@varLabels
    for (i in 1:length(vars)){
      model=paste(model, vars[i], ifelse(i<length(vars), '+', ''), sep='')
    }
  } 

  model2 <- paste("y ~",model)
  mat2 <- as.matrix(mat1)
  p <- dim(mat2)[1]
  n <- dim(mat2)[2]
#
# retrieve effect names
#
  owaov <- function(y)
  {
    formobj <- as.formula(model2)
    tmp <- row.names(anova(lm(formobj, data=eS@phenoData@pData)))
    return(tmp)
  }
  effnames <- owaov(mat2[1,])
#
# Perform ANOVA's and retrieve mean squares and df's
#
  numeff <- length(effnames)
  tmp1 <- rowaov(eS, model)
  #Check if overfit:
  if (is.null(tmp1)) {return(NULL)}
  #Otherwise proceed
  msmat <- tmp1[1:numeff,]
  dfmat <- tmp1[(numeff+1):(2*numeff),]
  numf <- numeff-1
  nu <- median(dfmat[numeff,])
  dfvec <- dfmat[numeff,]
  msevec <- msmat[numeff,]
#
# compute hyperparameters for prior
#

  mn <- mean(msevec)
  v <- var(msevec)
  alpha <- (mn^2-2*mn^2/nu +2*v)/(v-2*mn^2/nu)
  beta <- 1/(mn*(alpha-1))
  eta <- alpha*beta
  prior <- 1/eta
  df <- 2*alpha
#print(c(mn,v,alpha,beta,eta,prior,df))
#
# compute adjusted mean square errors
#
  adjdfvec <- dfvec+df
  adjmsevec <- (msevec*dfvec+prior*df)/adjdfvec
#
# compute F statistics and p-values
#
  pmat1 <- matrix(numeric(numf*p),ncol=numf)
  pmat2 <- pmat1
  for (i in 1:numf)
  { 
    pmat1[,i] <- 1-pf(msmat[i,]/msevec,dfmat[i,],dfvec) 
    pmat2[,i] <- 1-pf(msmat[i,]/adjmsevec,dfmat[i,],adjdfvec) 
  }
  colnames(pmat1) <- effnames[-numeff]
  colnames(pmat2) <- effnames[-numeff]
return(list("Gene.Specific"=pmat1,"Posterior"=pmat2))

}
