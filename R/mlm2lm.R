"mlm2lm" <-
function(lmobj,i)
{
  lmobj2 <- lmobj
  lmobj2$coefficients <- (lmobj$coefficients)[,i]
  lmobj2$residuals <- (lmobj$residuals)[,i]
  lmobj2$effects <- (lmobj$effects)[,i]
  lmobj2$fitted.values <- (lmobj$fitted.values)[,i]
  class(lmobj2) <- "lm"
  return(lmobj2)
}

