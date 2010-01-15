# Plots row standard deviation against row mean or rank of row mean
# Input data must be a matrix or data frame of expression values, or a
#    ExpressionSet object
plotMeanSD <- function(indata, by.rank = TRUE, line = FALSE, ymax = NULL){
	if (is.element(class(indata), c("matrix", "data.frame"))){
		mat <- as.matrix(indata)
	} else{
		if (is.element(class(indata), c("ExpressionSet", "AffyBatch"))){
			mat <- exprs(indata)
		} else{
			stop("class of 'indata' must be one of 'matrix',
				'data.frame','ExpressionSet', or 'AffyBatch'") 
		}
	}
	# remove NaN's and infinite values from mat
	noNanInf <- as.logical(prod(is.finite(mat)))
	if (noNanInf == FALSE){
		tmp <- is.finite(mat)
		test.vector <- apply(tmp, 1, prod)
		mat2 <- mat[test.vector ==1,]
		warning('rows with NaN and infinite values stripped from input data')
	} else{
		mat2 <- mat
	}

	msd.row<- function(data)
	{
 		mean.values <- apply(data,1,mean) 
 		standard.deviation<- sqrt(apply(data,1,var))   
		return(cbind(mean.values,standard.deviation))  
	}
	msd <- msd.row(mat2)

	mean <- msd[,1]
	sort.msd <- msd[order(mean),]
	index <- dim(mat2)[1]

	rank.mean <- c(1:index)
	standard.deviation <- msd[,2]
	sort.standard.deviation <- sort.msd[,2]
	
	if (by.rank == TRUE){
		xdata <- rank.mean
		ydata <- sort.standard.deviation
		xlabel <- "Rank of Mean"
		mainlabel <- "Rank of Mean vs. Standard Deviation"
	} else{
		xdata <- mean
		ydata <- standard.deviation
		xlabel <- "Mean"
		mainlabel <- "Mean vs. Standard Deviation"
	}
	if (is.null(ymax)){
		plot.ylim <- NULL
	} else{
		plot.ylim <- c(0, ymax)
	} 
	plot(ydata~xdata,main=mainlabel,xlab=xlabel,ylab="Standard Deviation",ylim=plot.ylim)
	if (line == TRUE){
		lines(lowess(y = ydata, x = xdata), col = "blue")
	}
} 

