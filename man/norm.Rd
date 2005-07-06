\name{norm}
\alias{norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Additive normalization function
 }
\description{
  This function normalizes the matrix in additive way.
}
\usage{
norm(mat1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat1}{ A matrix data to be normalized }
}
\value{
  \item{matnorm }{Normalized matrix}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{lnorm}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
LoggedSmpd<-norm(log(sample.mat))

}
\keyword{ math }% at least one, from doc/KEYWORDS

