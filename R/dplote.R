"dplote" <-
function(scoef,cs='')
{
  csum <- apply(scoef,2,sum)
  ind <- !is.na(csum)
  scoef <- scoef[,ind]
  for(i in 1:dim(scoef)[1]){
    assign(paste('d', i, sep=''),scoef[i,])
  }
  ts <- paste("Cumulative Frequency Distribution",cs,sep=" ")

  plot.stepfun(ecdf(d1),do.points=FALSE, main=ts,lwd=2,xlim=c(0,1), xlab="Relative Mean Square", ylab="Cumulative Density")
  for(i in 2:dim(scoef)[1]){
    plot.stepfun(ecdf(get(paste('d', i, sep=''))),do.points=FALSE,lwd=2,col.hor=i,col.vert=i,add=TRUE)
  }

  legend(.6,.4,col=1:dim(scoef)[1],legend=rownames(scoef), lty=1,lwd=2)
}

