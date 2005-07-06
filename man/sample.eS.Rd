\name{sample.eS}
\alias{sample.eS}
\docType{data}
\title{ Sample array data for LMGene}
\description{
  'exprSet' class data 
}
\usage{data(sample.eS)}
\format{
  The format is:
Formal class 'exprSet' [package "Biobase"] with 6 slots
  ..@ exprs   
  ..@ se.exprs 
  ..@ phenoData
  ..@ description
  ..@ annotation
  ..@ notes 
}
\details{
  identical with 'neweS(sample.mat,vlist)'
}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

Smpdt<-neweS(sample.mat,vlist)

data(sample.eS)
identical(sample.eS, Smpdt)
}
\keyword{datasets}
