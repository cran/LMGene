\name{rowaov}
\alias{rowaov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gene by gene anova functioin
 }
\description{
  Computes the mean squares and degrees of freedom for gene-by-gene anovas.
}
\usage{
rowaov(eS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data. must be a exprSet type and the log-transformation and the normalization of exprSet@exprs are recommended }
}
\details{
  The input argument, eS, must be exprSet type from Biobase package. 
  If you have a matrix data and information about the considered factors,
  then you can use \code{\link{neweS}} to conver the data into exprSet.
  Please see \code{\link{neweS}} in more detail.
}
\value{
  \item{resmat }{ A matrix of MSE and DF of all factors for all genes}
}
\references{ 
David M. Rocke (2004), Design and analysis of experiments with high throughput biological assay data, Seminars in Cell & Developmental Biology, 15, 703-713. 

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{genediff}}, \code{\link{mlm2lm}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd)
data(vlist)
LoggedSmpd0<-neweS(lnorm(log(Smpd)),vlist)

resmat <- rowaov(LoggedSmpd0)
resmat[,1:3]
}
\keyword{ models }% at least one, from doc/KEYWORDS

