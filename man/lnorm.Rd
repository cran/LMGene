\name{lnorm}
\alias{lnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lowess normalization function }
\description{
  Lowess normalization function
}
\usage{
lnorm(mat1, span = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat1}{ A matrix data to be normalized }
  \item{span}{ A parameter for lowess }
}
\details{
  \code{mat1} must be a \code{n}by\code{p} matrix, where \code{n} is the number of genes and \code{p} is the number of expression levels for each gene.
}
\value{
  \item{matnorm1 }{Normalized matrix}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{  \code{\link{norm}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
LoggedSmpd<-lnorm(log(sample.mat))

}
\keyword{ math }% at least one, from doc/KEYWORDS

