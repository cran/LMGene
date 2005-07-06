\name{msecalc}
\alias{msecalc}
\title{MSE calculation function}
\description{
  Computes the mean square error and gradient for the global ANOVA.
}
\usage{
msecalc(eS, lam, alpha, lowessnorm, R)
}
\arguments{
  \item{eS}{Array data. must be an \code{ExpressionSet} object.}
  \item{lam}{A parameter for glog transformation.}
  \item{alpha}{A parameter for glog transformation.}
  \item{lowessnorm}{TRUE, if lowess method is going to be used.}
  \item{R}{The residual matrix, i.e., identity minus the hat matrix.}
}
\details{
  The argument \code{eS} must be an \code{ExpressionSet} object from the Biobase package. 
  If you have a data in a \code{matrix} and information about the considered factors, then you
  can use \code{\link{neweS}} to convert the data into an \code{ExpressionSet} object. Please see
  \code{\link{neweS}} in more detail.
}
\value{
  \item{msev }{A vector which contains MSE and gradient of two parameters.}
}
\references{ 
B. Durbin and D.M. Rocke, (2003) Estimation of Transformation Parameters for Microarray Data,  Bioinformatics, 19, 1360-1367.

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{David Rocke and Geun-Cheol Lee}
\seealso{\code{\link{jggrad2}}, \code{\link{tranest2}}}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.eS)

lmod <- GetLMObj(sample.eS)
X <- lmod$x
U <- svd(X)$u
H <- crossprod(t(U), t(U))
n <- dim(H)[1]
R <- diag(rep(1,n)) - H

msecalc(sample.eS,500,50, FALSE, R)
}
\keyword{math}
