"dplot" <-
function(scoef, cs ='')
{
  csum <- apply(scoef,2,sum)
  ind <- !is.na(csum)
  scoef <- scoef[,ind]
  cut <- .1
  
  ymax=0
  for(i in 1:dim(scoef)[1]){
    assign(paste('d', i, sep=''),density(scoef[i,]))
    ymax=max(ymax,get(paste('d', i, sep=''))$y)
  }

  ts <- paste("Smoothed Histogram",cs,sep=" ")
  
  plot(d1,ylim=c(0,ymax), xlab="Relative Mean Square",  main=ts, lwd=2, col=1)
  for(i in 2:dim(scoef)[1]){
    lines(get(paste('d', i, sep='')),col=i,lwd=2)
  }

  legend(.6,ymax*4/5,col=1:dim(scoef)[1],legend=rownames(scoef),  lty=1,lwd=2)

}

