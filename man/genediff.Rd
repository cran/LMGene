\name{genediff}
\alias{genediff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Raw p-value calculation function
 }
\description{
 Computes two vectors of p-values per gene or probe
   using gene-by-gene anova with individual gene MSE using 
   both the gene-specific MSE and the posterior mean MSE for
   each term in the anova. 
   \cr
   Assumes a fixed effects model and the correct denominator for all comparisons is the MSE
   
}
\usage{
genediff(eS)
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
  \item{pvlist }{a list containing two sets of p-values obtained by gene specific MSE and the posterior MSE methods}
}
\references{ 
David M. Rocke (2004), Design and analysis of experiments with high throughput biological assay data, Seminars in Cell & Developmental Biology, 15, 703-713. 

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{LMGene}}, \code{\link{rowaov}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd)
data(vlist)
LoggedSmpd0<-neweS(lnorm(log(Smpd)),vlist)

pvlist<-genediff(LoggedSmpd0)
pvlist$Posterior[1:5,]

}
\keyword{ models }% at least one, from doc/KEYWORDS

