"LMGene" <-
function(eS,level=.05)
#
#
#
{
  pvlist <- genediff(eS)
  apvlist <- pvadjust(pvlist)
  numeff <- ncol(apvlist$Posterior.FDR)
  for (effnum in 1:numeff)
  {
    tmp <- rowlist(eS@exprs,effnum,apvlist,level)
    if (effnum==1){
      if (length(tmp)>1  ){
        lmres <- list(tmp=tmp)
      }else if(tmp != -1){
        lmres <- list(tmp=tmp)
      }else{
        lmres <- list(tmp="No significant genes")
      }
    }else {
      if ( length(tmp)>1 ){
        lmres <- c(lmres,list(tmp=tmp))
      }else if(tmp != -1){
        lmres <- c(lmres,list(tmp=tmp))
      }else{
        lmres <- c(lmres,list(tmp="No significant genes"))
      }
    }
    effname <- colnames(apvlist$Posterior.FDR)[effnum]
    names(lmres)[effnum] <- effname
  }
  return(lmres)
}

