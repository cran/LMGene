\name{neweS}
\alias{neweS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coercing to exprSet function 
 }
\description{
  This function converts a matrix data and its experimental data into an object of 'exprSet' class.
}
\usage{
neweS(mat, vlist, vlabel = as.list(names(vlist)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{ A matrix data to be converted }
  \item{vlist}{ A list which contains several factors representing the experiment description }
  \item{vlabel}{ A list of labels for the variables represented by the columns of \code{pData} of the 'exprSet' object to be made }
}
\details{
  Must load Biobas package first for converting. \cr
  \code{vlist} contains all the considered factors in which level values of each element represent the corresponding column of \code{mat}.
}
\value{
  \item{eset }{The converted object of 'exprSet' class}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{as}} }
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
\keyword{ manip }% at least one, from doc/KEYWORDS

