"neweS" <-
function(mat, vlist, vlabel=as.list(names(vlist))){
names(vlabel)=names(vlist)
pdata <- new("phenoData", pData=as.data.frame(vlist), varLabels=vlabel)
eset <- new("exprSet", exprs=as.matrix(mat), phenoData=pdata)
return(eset)
}

