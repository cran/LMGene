"rowaov" <-
function (eS, model=NULL) 
{
    mat1 <- as.matrix(exprs(eS))
    for (i in 1:length(varLabels(eS))) {
        assign(paste("x", i, sep = ""), pData(eS)[, i])
    }
    if (is.null(model)) {
      model <- ""
      for (i in 1:length(varLabels(eS))) {
        model <- paste(model, paste("x", i, sep = ""), ifelse(i < 
            length(varLabels(eS)), "+", ""), sep = "")
      }
    }
    model2 <- paste("y ~", model)
    mat2 <- as.matrix(mat1)
    n <- dim(mat2)[2]
    p <- dim(mat2)[1]
    y <- t(mat2)
    formobj <- as.formula(model2)
    tmp <- lm(formobj, x=TRUE, data=pData(eS))
    if (tmp$df <= 0) {
      print("Error: model is overfit. Try a simpler model.")
      return(NULL)
    }
    for (i in 1:p) {
      tmp2 <- mlm2lm(tmp, i)
      tmp3 <- anova(tmp2)
      tmp4 <- c(tmp3$Mean, tmp3$Df)
      if (i == 1) {
          resmat <- tmp4
      }
      else {
          resmat <- cbind(resmat, tmp4)
      }
    }
    return(resmat)
}
