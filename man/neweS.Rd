\name{neweS}
\alias{neweS}
\title{Coercing to an ExpressionSet code}
\description{
  This function converts a matrix data and its experimental data into an object of 'ExpressionSet' class.
}
\usage{
neweS(mat, vlist, vlabel = as.list(names(vlist)))
}
\arguments{
  \item{mat}{A matrix data to be converted.}
  \item{vlist}{A list which contains several factors representing the experiment description.}
  \item{vlabel}{A list of labels for the variables represented by the columns of \code{pData} of the 'ExpressionSet' object to be made.}
}
\details{
  Must load Biobase package first for converting. \cr
  \code{vlist} contains all the considered factors in which level values of each element represent the corresponding column of \code{mat}.
}
\value{
  \item{eset }{The converted object of 'ExpressionSet' class.}
}
\references{\url{http://www.idav.ucdavis.edu/~dmrocke/}}
\author{David Rocke and Geun-Cheol Lee}
\seealso{\code{\link[methods]{as}}}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

Smpdt <- neweS(sample.mat,vlist)

data(sample.eS)
identical(sample.eS, Smpdt)
}
\keyword{manip}
