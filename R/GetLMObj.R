"GetLMObj" <-
function (eS, model=NULL) 
{
  mat1 <- as.matrix(exprs(eS))

  # model information 

  if (is.null(model)) {
    model <- ''
    vars <- varLabels(eS)
    for (i in 1:length(vars)){
      model <- paste(model, vars[i], ifelse(i<length(vars), '+', ''), sep='')
    }
  }

  model2 <- paste("y ~",model)
  y <- mat1[1,]
  formobj <- as.formula(model2)
  mod <- lm(formobj, x=TRUE, data=pData(eS))
  return(mod)
}
