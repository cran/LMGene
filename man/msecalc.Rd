\name{msecalc}
\alias{msecalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MSE calculation function
 }
\description{
  Computes the mean square error and gradient for the global anova
}
\usage{
msecalc(eS, lam, alpha, lowessnorm, R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data. must be exprSet type. }
  \item{lam}{ A parameter for glog transformation }
  \item{alpha}{ A parameter for glog transformation }
  \item{lowessnorm}{ TRUE, if lowess method is going to be used }
  \item{R}{ The residual matrix, i.e., identity minus the hat matrix}
}
\details{
  The input argument, eS, must be exprSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into exprSet.
  Please see \code{\link{neweS}} in more detail.
}
\value{
  \item{msev }{A vector which contains MSE and gradient of two parameters}
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

lmod <- GetLMObj(sample.eS)
X <- lmod$x
U <- svd(X)$u
H = crossprod(t(U), t(U))
n = dim(H)[1]
R = diag(rep(1,n)) - H

msecalc(sample.eS,500,50, FALSE, R)

}
\keyword{ math }% at least one, from doc/KEYWORDS

