"neweS" <-
function(mat, vlist, annotat="", vlabel=as.list(names(vlist)))
{
  for (i in 1:length(vlist)) {
    names(vlist[[i]]) <- colnames(mat)
  }

  metdat <-data.frame(labelDescription=names(vlist))
  pdata <- new("AnnotatedDataFrame", data=as.data.frame(vlist),varMetadat=metdat)
  eset<-new("ExpressionSet",exprs=as.matrix(mat),phenoData=pdata,annotation=annotat)

  return(eset)
}


