"arrplot1" <-
function(eS)
{
  slides <- as.matrix(eS@exprs)
  
  # model information 
  for(i in 1:length(eS@phenoData@varLabels)){
  assign(paste('x', i, sep=''),as.factor(eS@phenoData@pData[,i]))
  }
  
  fchar=''
  for(i in 1:length(eS@phenoData@varLabels)){
  fchar=paste(fchar, paste('x', i, sep=''), ifelse(i<length(eS@phenoData@varLabels), '+', ''), sep='')
  }

  fchar2 <- paste("y ~",fchar)
  y <- t(slides)
  formobj <- as.formula(fchar2)

  n <- dim(slides)[2]
  p <- dim(slides)[1]
  tmp <- lm(formobj)
  for (i in 1:p)
  {
    tmp2 <- mlm2lm(tmp,i)
    tmp3 <- anova(tmp2)$Mean
    tmp4 <- msa(tmp3)
    if ( i == 1)
    {
      resmat <- tmp4
    }
    else
    {
      resmat <- cbind(resmat,tmp4)
    }
  }
  
  rownames(resmat)=c( names(eS@phenoData@varLabels) ,"Error")

  return(resmat)
}

