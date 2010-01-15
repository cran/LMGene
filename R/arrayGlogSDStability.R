# transform data by glog and then compute standard deviation #across arrays  
arrayGlogSDStability <- function(eS,lambda,alpha,lowessnorm,rank,model.based,R,rep.arrays)
{
  # Code by Shiquan Wu, 2009
  # Modified to add switch for using rank vs mean
  # Modified to add switch for using lowess normalization vs. #simple normalization
  # Modified to use sd of residuals from linear model in addition to replicate sd
  # Blythe Durbin-Johnson, Dec. 2009  
  # transform data by glog, followed by lnorm
  exprs.glog <- transeS(eS, lambda, alpha)
  if (lowessnorm == TRUE){
  	exprs.glog.norm <- lnormeS(exprs.glog)
  } else{
	normat <- norm(exprs(exprs.glog))
      exprs.glog.norm <- neweS(normat, pData(eS))
  }
  
  exprs.glog.norm <- exprs(exprs.glog.norm)
  if (model.based == TRUE){
      ngenes <- dim(exprs.glog.norm)[1]
      mean.unsorted <- apply(exprs.glog.norm, 1, mean)
      mean.order <- order(mean.unsorted)
      mean <- mean.unsorted[mean.order]
      mean.rank <- c(1:ngenes)
	resids <- exprs.glog.norm%*%R
      resid.sd <- apply(resids, 1, sd)
	resid.sd.sorted <- resid.sd[mean.order]
      ydata <- resid.sd.sorted
  } else{  
	pooled.sd <- function(x, rep.arrays){
		k <- length(rep.arrays)
		s2 <- vector(length = k)
		n <- vector(length = k)
		for (i in 1:k){
			arrays <- rep.arrays[[i]]
			xi <- x[arrays]
			xi <- xi[!is.na(xi)]
			s2[i] <- var(xi)
			n[i] <- length(xi)
		}
		s2p <- sum(s2*(n-1))/sum(n-1)
		sp <- sqrt(s2p)
		return(sp)
	}
	sd.unsorted <- apply(exprs.glog.norm, 1, pooled.sd, rep.arrays = rep.arrays)
	mean.unsorted <- apply(exprs.glog.norm[,unlist(rep.arrays)], 1, mean)
      mean.order <- order(mean.unsorted)
	mean <- mean.unsorted[mean.order]
	standard.deviation <- sd.unsorted[mean.order]

	ngenes <- dim(exprs.glog.norm)[1]
	mean.rank <- 1:ngenes

	ydata <- standard.deviation
  }

  # stability
  if (rank == TRUE){
  	sd.vs.mean.regression <- lm(ydata ~ mean.rank)
  } else {
	sd.vs.mean.regression <- lm(ydata ~ mean)
  }
  coef.regression <- coef(sd.vs.mean.regression)
  score <- ngenes*abs(coef.regression[2]) 
  return(score)
}