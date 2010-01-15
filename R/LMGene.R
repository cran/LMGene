#LMGene calls function genediff to calculate the unadjusted #gene-specific and posterior p-values of all genes
#and then calls function pvadjust to calculate the FDR-adjusted #p-values of all genes.
#Significant genes for each factor in model (based on either the #gene-specific or posterior FDR-adjusted p-values) are output.
LMGene <-
function (eS, model=NULL, level = 0.05, posterior = FALSE, method = c("MLE", "MOM", "MOMlog")) 
{
    method <- match.arg(method)
    if (class(eS) != "ExpressionSet"){
    	stop("'eS' must be an object of class 'ExpressionSet'")
    }
    if (level < 0 || level > 1){
    	stop("'level' must be between 0 and 1")
    }
    if (!is.element(method, c("MLE", "MOM", "MOMlog"))){
	stop("'method' must be one of 'MLE', 'MOM', 'MOMlog'")
    }

    pvlist <- genediff(eS, model, method, verbose = FALSE)
    #Check for overfitting
    if (is.null(pvlist)) {return(NULL)}
    #Otherwise proceed
    apvlist <- pvadjust(pvlist)
    numeff <- ncol(apvlist$Posterior.FDR)
    for (effnum in 1:numeff) {
        tmp <- rowlist(exprs(eS), effnum, apvlist, level, posterior)
        if (effnum == 1) {
            if (length(tmp) > 1) {
                lmres <- list(tmp = tmp)
            }
            else if (tmp != -1) {
                lmres <- list(tmp = tmp)
            }
            else {
                lmres <- list(tmp = "No significant genes")
            }
        }
        else {
            if (length(tmp) > 1) {
                lmres <- c(lmres, list(tmp = tmp))
            }
            else if (tmp != -1) {
                lmres <- c(lmres, list(tmp = tmp))
            }
            else {
                lmres <- c(lmres, list(tmp = "No significant genes"))
            }
        }
        effname <- colnames(apvlist$Posterior.FDR)[effnum]
        names(lmres)[effnum] <- effname
    }
    return(lmres)
}
