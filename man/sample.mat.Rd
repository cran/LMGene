\name{sample.mat}
\alias{sample.mat}
\docType{data}
\title{ Sample array data for LMGene package}
\description{
  A \code{matrix} of array data
}
\usage{data(sample.mat)}
\format{
  A \code{matrix} measuring 613 rows (probes) by 32 columns (samples).
}
\examples{
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

Smpdt<-neweS(sample.mat,vlist)

data(sample.eS)
identical(exprs(sample.eS), exprs(Smpdt))
identical(pData(sample.eS), pData(Smpdt))
}
\keyword{datasets}
