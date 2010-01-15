\name{lnorm}
\alias{lnorm}
\title{ Lowess normalization function }
\description{
  Lowess normalization function
}
\usage{
lnorm(mat1, span = 0.1)
}
\arguments{
  \item{mat1}{ A data matrix to be normalized }
  \item{span}{ Lowess smoother span.  Larger values give more smoothness.}
}
\details{
  \code{mat1} must be a \code{p} by \code{n} matrix, where \code{p} is the number of genes and \code{n} is the number of arrays or samples}
\value{
  \item{matnorm1}{Normalized matrix}
}
\references{ \url{http://dmrocke.ucdavis.edu} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{\code{\link{lnormeS}}, \code{\link{norm}}}
\examples{
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

raw.eS <- neweS(sample.mat, vlist)

# glog transform data
trans.eS <- transeS(raw.eS, lambda = 727, alpha = 56)

# normalize
normed.exprs <- lnorm(exprs(trans.eS))
}
\keyword{ smooth}

