\name{dplote}
\alias{dplote}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting cumulative frequency distribution
 }
\description{
  This function is a sub-function of \code{arrplote} and draws a cumlative frequency function.
}
\usage{
dplote(scoef, cs = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scoef}{ A matrix containing the relative mean square information of all factors for all genes }
  \item{cs}{ Addable caption in main title of the graph }
}
\details{
  \code{scoef} can be made by the function \code{\link{arrplot1}}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{arrplot1}}, \code{\link{arrplote}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd0)

#
RMS<-arrplot1(Smpd0)
dplote(RMS, "Example")

}
\keyword{ hplot }% at least one, from doc/KEYWORDS

