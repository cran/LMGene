\name{sample.eS}
\alias{sample.eS}
\docType{data}
\title{Sample array data for LMGene}
\description{
  Sample 'ExpressionSet' class data.
}
\usage{data(sample.eS)}
\format{
  Formal class 'ExpressionSet' [package "Biobase"].
}
\details{
  identical with 'neweS(sample.mat, vlist)'
}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

Smpdt <- neweS(sample.mat,vlist)

data(sample.eS)
identical(sample.eS, Smpdt)
}
\keyword{datasets}
