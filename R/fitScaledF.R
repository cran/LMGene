fitScaledF <- function(vars, nu, startvals){	
# Fits scaled F distribution by maximum likelihood, as in
# Simon and Wright (2003)
# Input variables: vars = Vector of length p of MSEs from linear #model fit to each gene
#			 nu = Scalar degrees of freedom from linear #model fit to each gene, assumed common 
#			 startvals = Initial estimates of alpha and beta
# Value: Vector of length 2 containing MLEs of alpha and beta
	like <- function(x){
		lna <- x[1]
		a <- exp(lna)
		lnb <- x[2]
		b <- exp(lnb)
		eta <- a*b
		loglik <- sum(log(eta*df(eta*vars, df1 = nu, df2 = 2*a)))
		return(-loglik)
	}

	a.start <- startvals[1]
	b.start <- startvals[2]
	startvar <- c(log(a.start), log(b.start))
	typsize <- startvar

	opt <- optim(startvar, like, method = "BFGS")
	return(exp(opt$par))
}