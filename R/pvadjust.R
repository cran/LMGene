pvadjust <-
function(pvlist)
#
# pvlist is the output from genediff containing p-values from
#   gene-specific MSE's and posterior MSE's. This routine
#   adds FDR adjusted p-values using the multtest routine
#   rawp2adjp
#
{
  library(multtest)
  pv1 <- pvlist$Gene.Specific
  pv2 <- pvlist$Posterior
  nump <- dim(pv1)[2]
  for (i  in 1:nump)
  {
    ap <- mt.rawp2adjp(pv1[,i],"BH")
    pv1[,i] <- ap$adjp[order(ap$index),2]
    ap <- mt.rawp2adjp(pv2[,i],"BH")
    pv2[,i] <- ap$adjp[order(ap$index),2]
  }
  pvlist2 <- c(pvlist,list("Gene.Specific.FDR"=pv1,"Posterior.FDR"=pv2))
  return(pvlist2)
}
