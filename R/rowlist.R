"rowlist" <-
function(genemat,effnum,apvlist,level,posterior=TRUE)
#
# genemat is an n-by-p matrix of expression values 
# effnum is the coumn number for the effect of interest
# apvlist is a matrix of p-values from pvadjust or genediff
# the routine returns a list of genes whose FDR p-value is 
#   less than level using either individual gene or posterior 
#   MSE's. This is gene names if rownames(genemat) is not null,
#   and gene numbers otherwise.
#
{
  if(posterior)
  {
    ind <- apvlist$Posterior.FDR[,effnum] < level
  }
  else
  {
    ind <- apvlist$Gene.Specific.FDR[,effnum] < level
  }
  numsig <- sum(ind)
  if (is.null(rownames(genemat)))
  {
    p <- dim(genemat)[1]
    if (numsig > 0)
    {
      return((1:p)[ind])
    }
    else
    {
      return(-1)
    }
  }
  else
  {
    if (numsig > 0)
    {
      return(rownames(genemat)[ind])
    }
    else
    {
      return(-1)
    }
  }
}

