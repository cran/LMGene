\name{arrplote}
\alias{arrplote}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for plotting a cumulative frequency distribution }
\description{
  This function calls \code{\link{arrplot1}} first to calculate the relative mean square values of the considered factors,
  then calls \code{\link{dplote}} to plot the obtained relative mean square values in a cumulative frequency distribution
}
\usage{
arrplote(eS, cs = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array Data. Must be a exprSet type }
  \item{cs}{ Addable caption in main title of the graph }
}
\value{
  The input argument, eS, must be exprSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into exprSet.
  Please see \code{\link{neweS}} in more detail.
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{arrplot1}}, \code{\link{dplote}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd0)

#
arrplote(Smpd0)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

