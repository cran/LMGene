\name{arrplot1}
\alias{arrplot1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate relative mean squares }
\description{
  Calculate relative mean squares of all considered factors for all genes. 
  The obtained relative mean square values are used for plotting the smoothed histogram and the cumulative frequency distributioin.
}
\usage{
arrplot1(eS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data. must be a exprSet type }
}
\details{
  The input argument, eS, must be exprSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into exprSet.
  Please see \code{\link{neweS}} in more detail.
}
\value{
  \item{resmat }{Matrix containing the relative mean square values of all factors considered including Error for all genes}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{arrplotd}}, \code{\link{arrplote}}, \code{\link{mlm2lm}}, \code{\link{msa}} }
\examples{

#library
library(Biobase)
library(LMGene)

#data
data(Smpd0)

#
RMS<-arrplot1(Smpd0)
RMS[,1:5]
}

\keyword{ dplot }% at least one, from doc/KEYWORDS

