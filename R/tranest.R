#Estimates parameters for the glog transformation, by maximum 
#likelihood or by minimizing the stability score.
tranest <- 
function (eS, ngenes = -1, starting = FALSE, lambda = 1000, alpha = 0,
    gradtol = 1e-3, lowessnorm = FALSE, method=1, mult=FALSE, model=NULL, 
	SD = FALSE, rank = TRUE, model.based = TRUE, rep.arrays = NULL)
{
    if (class(eS) != "ExpressionSet"){
	stop("'eS' must be an object of class 'ExpressionSet'")
    }
    if (SD == TRUE && mult == TRUE){
	warning('estimation of vector alpha not implemented
		   for minimum stability score method, defaulting
               to scalar alpha')
	mult <- FALSE
    }
    if (SD == FALSE && model.based == FALSE){
	warning("'SD = FALSE', ignoring input variable 'model.based'")
    }
    if (SD == FALSE && !is.null(rep.arrays)){
	warning("'SD = FALSE', ignoring input variable 'rep.arrays'")
    }
    if (SD == TRUE && model.based == TRUE && !is.null(rep.arrays)){
	warning("'model.based = TRUE', ignoring input variable 'rep.arrays'")
    }
    if (model.based == FALSE && !is.null(model)){
	warning("'model.based = FALSE', ignoring input variable 'model'")
    }
    if (model.based == FALSE && is.null(rep.arrays)){
	stop("if 'model.based = FALSE', 'rep.arrays' must be specified")
    }
    if (!is.element(method,1:5)){
	stop("'method' must be an integer between 1 and 5")
    } 
    if (method == 5 && mult == FALSE){
	stop("if 'mult = FALSE', 'method' must be between 1 and 4")
    }
    if (length(alpha) > 1 && mult == FALSE){
	stop("vector alpha requires 'mult = TRUE'")
    }
    if (lambda <= 0 || length(lambda) > 1){
	stop("'lambda' must be a positive scalar")
    }
    if (length(rep.arrays) > 0){
	for (i in 1:length(rep.arrays)){
		tmp <- length(rep.arrays[[i]])
            if (tmp < 2){
			stop("each element of 'rep.arrays' must have length at least 2")
		}
	}
    }
    if (!identical(unlist(rep.arrays), unique(unlist(rep.arrays)))){
		stop("elements of 'rep.arrays' may not overlap")
    }

    mat1 <- as.matrix(exprs(eS))
    n <- dim(mat1)[2]
    p <- dim(mat1)[1]

    if (length(alpha) > 1 && length(alpha)!= n){
	stop("vector alpha must have length equal to number of arrays")
    }

    if (p > 100000) {ngenes <- 50000}
    if ((ngenes < p) & (ngenes > 0))
        mat2 <- mat1[sample(p, ngenes), ]
    else mat2 <- mat1
    eS2 <- new("ExpressionSet", exprs = mat2, phenoData = phenoData(eS))

    if (mult==FALSE) {
      tranpar <- tranest2(eS2, starting, lambda, alpha, gradtol,
        lowessnorm, method, model, SD, rank, model.based, rep.arrays)
      return(list(lambda = (tranpar[1]), alpha = tranpar[2:length(tranpar)]))
    }
    else {
      return(tranestmult(eS2, starting, lambda, alpha, gradtol,
        lowessnorm, method, 200, model))
    }
}
