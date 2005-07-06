"neweS" <-
function(mat, vlist, vlabel=as.list(names(vlist))){
names(vlabel)=names(vlist)

#Must add appropriate names to the variables in vlist for R 2.3 compatibility.
for (i in 1:length(vlist)) {
  names(vlist[[i]]) <- colnames(mat)
}
###

pdata <- new("phenoData", pData=as.data.frame(vlist), varLabels=vlabel)
eset <- new("exprSet", exprs=as.matrix(mat), phenoData=pdata)
return(eset)
}
