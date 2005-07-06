% --- Source file: sample.eS.Rd ---
\name{sample.eS}
\alias{sample.eS}
\docType{data}
\title{ Sample array data for LMGene}
\description{
  'ExpressionSet' class data 
}
\usage{data(sample.eS)}
\format{
  The format is:
Formal class 'ExpressionSet' (See package "Biobase" for details) 
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

