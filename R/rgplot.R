"rgplot" <-
function(eS, norm="lnorm",red=TRUE,green=TRUE,title="",span=.1)
# plots the individual dye/slide distributions either with or without 
# normalization

{
  slides <- as.matrix(eS@exprs)
  factor=as.factor(eS@phenoData@pData[,1])
  

  drawden<-function(slides){
    if(factor[1]==levels(factor)[1])
      plot(density(slides[,1]), main=title, ylab="Density", col=2)
    else
      plot(density(slides[,1]), main=title, ylab="Density", col=3)
    for(i in 2:length(factor)){
      if(factor[i]==levels(factor)[1]){
        lines(density(slides[,i]),col=2)
      }
      else{
        lines(density(slides[,i]),col=3)
      }    
    }
  }
  
  slides<-log(slides)
  if (norm=="norm")  {
    slides <- norm(slides)
  }
  else if (norm=="lnorm")   {
    slides <- lnorm(slides,span=span)
  }
  drawden(slides)

}

