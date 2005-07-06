% --- Source file: rowaov.Rd ---
\name{rowaov}
\alias{rowaov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gene by gene anova functioin
 }
\description{
  Computes the mean squares and degrees of freedom for gene-by-gene anovas.
}
\usage{
rowaov(eS, model=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data. must be a ExpressionSet type and the log-transformation and the normalization of exprs(eS) are recommended }
  \item{model}{ Model used for comparison. See details and \code{\link{LMGene}}.}

}
\details{
  The input argument, eS, must be ExpressionSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into ExpressionSet.
  Please see \code{\link{neweS}} in more detail.

  'model' is an optional character string, constructed like the right-hand side of a formula for lm. It specifies which of the variables in the ExpressionSet will be used in the
model and whether interaction terms will be included. If model=NULL, it uses all variables
from the ExpressionSet without interactions. Be careful of using interaction terms with factors: this often leads to overfitting, which will yield an error.

}
\value{
  \item{resmat }{ A matrix of MSE and DF of all factors for all genes}
}
\references{ 
David M. Rocke (2004), Design and analysis of experiments with high throughput biological assay data, Seminars in Cell & Developmental Biology, 15, 703-713. 

\url{http://dmrocke.ucdavis.edu/} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{genediff}}, \code{\link{mlm2lm}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)
LoggedSmpd0<-neweS(lnorm(log(sample.mat)),vlist)

resmat <- rowaov(LoggedSmpd0)
resmat[,1:3]
}
\keyword{ models }% at least one, from doc/KEYWORDS

