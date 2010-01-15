\name{norm}
\alias{norm}
\title{Normalization function
 }
\description{
  This function normalizes a matrix by subtracting the column (sample) mean from each element and adding the grand mean.
}
\usage{
norm(mat1)
}
\arguments{
  \item{mat1}{ A \code{matrix} to be normalized }
}
\value{
  \item{matnorm}{Normalized \code{matrix}}
}
\references{ \url{http://dmrocke.ucdavis.edu} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{lnorm}}}
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
normed.exprs <- norm(exprs(trans.eS))
}


