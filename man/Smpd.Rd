\name{Smpd}
\alias{Smpd}
\docType{data}
\title{ Sample array data for LMGene package}
\description{
  A matrix of array data
}
\usage{data(Smpd)}
\format{
  A data frame with 500 observations on the following 16 variables.
  \describe{
    \item{X1}{a numeric vector}
    \item{X2}{a numeric vector}
    \item{X3}{a numeric vector}
    \item{X4}{a numeric vector}
    \item{X5}{a numeric vector}
    \item{X6}{a numeric vector}
    \item{X7}{a numeric vector}
    \item{X8}{a numeric vector}
    \item{X9}{a numeric vector}
    \item{X10}{a numeric vector}
    \item{X11}{a numeric vector}
    \item{X12}{a numeric vector}
    \item{X13}{a numeric vector}
    \item{X14}{a numeric vector}
    \item{X15}{a numeric vector}
    \item{X16}{a numeric vector}
  }
}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd)
data(vlist)

Smpdt<-neweS(Smpd,vlist)

data(Smpd0)
identical(Smpd0, Smpdt)
}
\keyword{datasets}
