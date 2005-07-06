\name{sample.mat}
\alias{sample.mat}
\docType{data}
\title{ Sample array data for LMGene package}
\description{
  A matrix of array data
}
\usage{data(sample.mat)}
\format{
  A data frame measuring 613 probes on the 32 samples.
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
