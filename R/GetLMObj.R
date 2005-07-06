"GetLMObj" <-
function (eS, model=NULL) 
{
 mat1 <- as.matrix(eS@exprs)
  
  # model information 

  if (is.null(model)) {
    model=''
    vars <- eS@phenoData@varLabels
    for (i in 1:length(vars)){
      model=paste(model, vars[i], ifelse(i<length(vars), '+', ''), sep='')
    }
  } 


  model2 <- paste("y ~",model)
  y <- mat1[1,]
  formobj <- as.formula(model2)
  mod <- lm(formobj, x=TRUE, data=eS@phenoData@pData)
  return(mod)
}
