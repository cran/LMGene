"rowaov" <-
function(eS)
#
# computes the mean squares and degrees of freedom for
#   gene-by-gene anovas. Assumes a fixed effects model and
#   the correct denominator for all comparisons is the MSE
# mat1 is a p by n matrix of expression values
# vlist is a list of the form list(v1=v1,v2=v2,..) of variables, 
#   each of length n, the number of chips
# fchar is the model for each gene in the form of a character 
#   string with variables from vlist. An example is "v1 + v2"
#
# returns matrix with p rows and as many twice as many columns 
#   as lines in the anova table
#
# create variables
#
{ 
  mat1 <- as.matrix(eS@exprs)
  
  # model information 
  for(i in 1:length(eS@phenoData@varLabels)){
  assign(paste('x', i, sep=''),as.factor(eS@phenoData@pData[,i]))
  }
  
  fchar=''
  for(i in 1:length(eS@phenoData@varLabels)){
  fchar=paste(fchar, paste('x', i, sep=''), ifelse(i<length(eS@phenoData@varLabels), '+', ''), sep='')
  }

  fchar2 <- paste("y ~",fchar)
  mat2 <- as.matrix(mat1)
  n <- dim(mat2)[2]
  p <- dim(mat2)[1]
#
# run regression and anovas
#

  y <- t(mat2)
   formobj <- as.formula(fchar2)
  tmp <- lm(formobj)
  for (i in 1:p)
  {
#    y <- mat2[i,]
#    tmp2 <- lm(formobj)
    tmp2 <- mlm2lm(tmp,i)
    tmp3 <- anova(tmp2)
    tmp4 <- c(tmp3$Mean,tmp3$Df)
     if ( i == 1)
    {
      resmat <- tmp4
    }
    else
    {
      resmat <- cbind(resmat,tmp4)
    }
  }
  return(resmat)
}

