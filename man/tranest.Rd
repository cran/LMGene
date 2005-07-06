\name{tranest}
\alias{tranest}
\title{Glog transformation parameter estimation function}
\description{
  Finds the best parameters for glog transformation.
}
\usage{
tranest(eS, ngenes = -1, starting = FALSE, lambda = 1000, alpha = 0, gradtol = 0.001, lowessnorm = FALSE, method=1, mult=FALSE, model=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{Array data. must be an \code{ExpressionSet} object.}
  \item{ngenes}{Number of genes that is going to be used for the parameter estimation.}
  \item{starting}{TRUE, if the given initial parameter values are used.}
  \item{lambda}{Initial parameter value for lambda.}
  \item{alpha}{Initial parameter value for alpha.}
  \item{gradtol}{a positive scalar giving the tolerance at which the scaled
          gradient is considered close enough to zero to terminate the
          algorithm.}
  \item{lowessnorm}{TRUE, if lowess method is going to be used.}
  \item{method}{Determines optimization method. Default is 1, which corresponds to a Newton-type method (see \code{\link{nlm}}). Method 2 is based on the Nelder-Mead method (see \code{\link{optim}}).}
  \item{mult}{If true, tranest will use a vector alpha with one entry per sample. Default is false (same alpha for every sample).}
  \item{model}{Specifies model to be used. Default is to use all variables from eS without interactions. See details.}

}
\details{
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
  \item{tranpar}{A list containing the best parameter for 'lambda' and 'alpha'.}
}
\references{ 
B. Durbin and D.M. Rocke, (2003) Estimation of Transformation Parameters for Microarray Data,  Bioinformatics, 19, 1360-1367.

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{ David Rocke, Geun-Cheol Lee and John Tillinghast}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.eS)

tranpar <- tranest(sample.eS, 100)
tranpar
tranpar <- tranest(sample.eS, mult=TRUE)
tranpar
}
\keyword{math}
