\name{beams}
\alias{beams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Beam search rountine for finding the optimal parameters
 }
\description{
  This routine finds optimal values of the two parameters, alpha and lambda, for the glog transformation function
}
\usage{
beams(eS, startvar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data. must be a exprSet type  }
  \item{startvar}{ Given (initial) value of log(alpha) and lambda }
}
\details{
  The input argument, eS, must be exprSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into exprSet.
  Please see \code{\link{neweS}} in more detail.
  \cr
  eS is used only as an argument when calling \code{msecalc} in the function.
}
\value{
  \item{bestsol }{A vector containing best alpha and lambda value}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{msecalc}}, \code{\link{tranest2}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd0)

#
mat1 <- as.matrix(Smpd0@exprs)
lamstart <- log(median(abs(mat1))^2)
alphastart <- quantile(abs(as.vector(mat1)), 0.1)
startvar <- c(lamstart, alphastart)

beams(Smpd0, startvar)
}
\keyword{ optimize }% at least one, from doc/KEYWORDS

