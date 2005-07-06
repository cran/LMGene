\name{tranestmult}
\alias{tranestmult}
\title{Glog transformation parameter estimation function for multiple parameters}
\description{
  A sub-function of \code{tranest} which searches the best parameters for glog transformation.
}
\usage{
tranestmult (eS, starting = FALSE, lambda = 1000, alpha = 0, gradtol = 0.001, lowessnorm=FALSE, method=1, max_iter=200, model=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{Array data. must be an \code{ExpressionSet} object.}
  \item{starting}{TRUE, if the given initial parameter values are used.}
  \item{lambda}{Initial parameter value for lambda.}
  \item{alpha}{Initial parameter value for alpha.}
  \item{gradtol}{a positive scalar giving the tolerance at which the scaled
          gradient is considered close enough to zero to terminate the
          algorithm.}
  \item{lowessnorm}{TRUE, if lowess method is going to be used.}
  \item{method}{Set optimization method; default is modified Gauss-Newton (nlm). See \code{\link{tranest}}.}
  \item{max_iter}{Max. number of iterations of \code{\link{nlm}} to use in optimization.}
  \item{model}{Model in terms of vlist which is compared to transformed expression data. See \code{\link{tranest}}.}
}
\details{
  This is primarily an internal function. The normal way of calling it would 
  be to call \code{\link{tranest}} with the option mult=TRUE.

  The argument \code{eS} must be an \code{ExpressionSet} object from the Biobase package. 
  If you have a data in a \code{matrix} and information about the considered factors, then you
  can use \code{\link{neweS}} to convert the data into an \code{ExpressionSet} object. Please see
  \code{\link{neweS}} in more detail.

  The \code{model} argument is an optional character string, constructed like the right-hand
  side of a formula for lm. It specifies which of the variables in the \code{ExpressionSet} will
  be used in the model and whether interaction terms will be included. If \code{model=NULL},
  it uses all variables from the \code{ExpressionSet} without interactions. Be careful of using
  interaction terms with factors; this often leads to overfitting, which will yield an error.
}
\value{
  \item{tranpar}{A list (not a vector) containing the best parameter for 'lambda' and the best vector for 'alpha'.}
}
\references{ 
B. Durbin and D.M. Rocke, (2003) Estimation of Transformation Parameters for Microarray Data,  Bioinformatics, 19, 1360-1367.

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{David Rocke and Geun-Cheol Lee}
\seealso{\code{\link{tranest}}, \code{\link{tranest2}}}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.eS)

tranpar <- tranestmult(sample.eS, lambda= 500, alpha=50)
tranpar
}
\keyword{math}
