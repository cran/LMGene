\name{LMGene}
\alias{LMGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LMGene main function }
\description{
  LMGene calls function \code{\link{genediff}} to calculate the raw p-values of all genes
  and then calls function \code{\link{pvadjust}} to calculate the adjusted p-values of all genes.
  Finally, calls function \code{\link{rowlist}} to list the names of genes that are selected as significant under the specified significance level.
  
}
\usage{
LMGene(eS, model=NULL, level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data. must be exprSet type and the log-transformation and the normalization of exprSet@exprs are recommended. }
  \item{model}{Specifies model to be used. Default is to use all variables from eS without interactions. See details.}
  \item{level}{ Significance level }
}
\details{
  The input argument, eS, must be exprSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into exprSet.
  Please see \code{\link{neweS}} in more detail.
  \cr
  \code{level} indicates False Discovery Rate, e.g. level=0.05 means 5% FDR.

  'model' is an optional character string, constructed like the right-hand side of a formula for lm. It specifies which of the variables in the exprSet will be used in the
model and whether interaction terms will be included. If model=NULL, it uses all variables
from the exprSet without interactions. Be careful of using interaction terms with factors: this often leads to overfitting, which will yield an error.
}
\value{
  \item{lmres }{ A list which contains significant gene names for each considered factor}
}
\references{ 
David M. Rocke (2004), Design and analysis of experiments with high throughput biological assay data, Seminars in Cell & Developmental Biology, 15, 703-713. 

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
 }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{genediff}}, \code{\link{pvadjust}}, \code{\link{rowlist}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)
LoggedSample<-neweS(lnorm(log(sample.mat)),vlist)

siggeneslist<-LMGene(LoggedSample, 'patient + dose', 0.01)

}
\keyword{ models }% at least one, from doc/KEYWORDS

