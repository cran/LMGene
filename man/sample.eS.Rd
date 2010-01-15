\name{sample.eS}
\alias{sample.eS}
\docType{data}
\title{Sample array data for LMGene}
\description{
  Sample \code{ExpressionSet} class data.
}
\usage{data(sample.eS)}
\format{
  Formal class \code{ExpressionSet} [package \code{Biobase}].
}
\details{
  Identical with \code{neweS(sample.mat, vlist)}, up to metadata
}
\examples{
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

Smpdt <- neweS(sample.mat,vlist)

data(sample.eS)
identical(exprs(sample.eS), exprs(Smpdt))
identical(pData(sample.eS), pData(Smpdt))
}
\keyword{datasets}
