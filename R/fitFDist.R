fitFDist <- function(x,df1) {
#	Moment estimation of the parameters of a scaled F-distribution
#	The first degrees of freedom is given
#	Gordon Smyth
#	8 Sept 2002.  Last revised 6 April 2006.
#     
#     From limma package in Bioconductor
#
#	Remove missing or infinite values and zero degrees of freedom
	o <- is.finite(x) & is.finite(df1) & (x >= 0) & (df1 > 0)
	if(any(!o)) {
		x <- x[o]
		df1 <- df1[o]
	}
	n <- length(x)
	if(n==0) return(list(scale=NA,df2=NA))

#	Avoid exactly zero values
	m <- median(x)
	if(m==0) {
		warning("More than half of residual variances are exactly zero: eBayes unreliable")
		m <- 1
	} else {
		if(any(x==0)) warning("Zero sample variances detected, have been offset",call.=FALSE)
	}
	x <- pmax(x, 1e-5 * median(x))

#	Better to work on with log(F)
	z <- log(x)
	e <- z-digamma(df1/2)+log(df1/2)
	emean <- mean(e)
	evar <- mean(n/(n-1)*(e-emean)^2-trigamma(df1/2))
	if(evar > 0) {
		df2 <- 2*trigammaInverse(evar)
		s20 <- exp(emean+digamma(df2/2)-log(df2/2))
	} else {
		df2 <- Inf
		s20 <- exp(emean)
	}
	list(scale=s20,df2=df2)
}
