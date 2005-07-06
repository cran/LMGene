\name{tranest2}
\alias{tranest2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Glog transformation parameter estimation function 2
 }
\description{
  A sub-function of \code{tranest} which search the best parameters for glog transformation.
}
\usage{
tranest2(eS, starting = FALSE, lambda = 1000, alpha = 0, gradtol = 0.001, lowessnorm, method=1, model=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data. must be exprSet type. }
  \item{starting}{ TRUE, if the given initial parameter values are used }
  \item{lambda}{ Initial parameter value for lambda }
  \item{alpha}{ Initial parameter value for alpha }
  \item{gradtol}{ a positive scalar giving the tolerance at which the scaled
          gradient is considered close enough to zero to terminate the
          algorithm }
  \item{lowessnorm}{ TRUE, if lowess method is going to be used }
  \item{method}{ Set optimization method; default is modified Gauss-Newton (nlm). See \code{\link{tranest}}.}
  \item{model}{ Model in terms of vlist which is compared to transformed expression data. See \code{\link{tranest}}.}
}
\details{
  The input argument, eS, must be exprSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into exprSet.
  Please see \code{\link{neweS}} in more detail.
  'model' is an optional character string, constructed like the right-hand side of a formula for lm. It specifies which of the variables in the exprSet will be used in the
model and whether interaction terms will be included. If model=NULL, it uses all variables
from the exprSet without interactions. Be careful of using interaction terms with factors: this often leads to overfitting, which will yield an error.

}
\value{
  \item{tranpar }{A numeric vector containing the best parameter for 'lambda' and 'alpha'}
}
\references{ 
B. Durbin and D.M. Rocke, (2003) Estimation of Transformation Parameters for Microarray Data,  Bioinformatics, 19, 1360-1367.

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{jggrad2}}, \code{\link{tranest2}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.eS)

tranpar <- tranest2(sample.eS, lambda= 500, alpha=50)
tranpar

}
\keyword{ math }% at least one, from doc/KEYWORDS

