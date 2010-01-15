\name{neweS}
\alias{neweS}
\title{Coerce a matrix to class ExpressionSet}
\description{
  This function converts a data matrix into an \code{ExpressionSet} object.
}
\usage{
neweS(mat, vlist, vlabel = as.list(names(vlist)))
}
\arguments{
  \item{mat}{A data matrix to be converted.}
  \item{vlist}{A list, each component of which describes a factor in the experimental design.}
  \item{vlabel}{A list of labels for each component of \code{vlist}.}
}
\details{
Each element of a component of \code{vlist} corresponds to a column of \code{mat}.
See \code{\link{vlist}} for an example.
}
\value{
  \item{eset}{An \code{ExpressionSet} object.}
}
\references{\url{http://dmrocke.ucdavis.edu}}
\author{David Rocke and Geun-Cheol Lee}
\seealso{\code{\link{vlist}}}
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
\keyword{manip}
