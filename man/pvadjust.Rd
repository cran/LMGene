\name{pvadjust}
\alias{pvadjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ P-value adjusting function
 }
\description{
  This function converts the given raw p-values into the FDR adjusted p-values using R package 'multtest'.
}
\usage{
pvadjust(pvlist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvlist}{ A list containing raw p-values }
}
\details{
 \code{pvlist} is the output from \code{genediff} containing p-values from
   gene-specific MSE's and posterior MSE's.
}
\value{
  \item{pvlist2 }{A list with the raw p-values and the newly computed FDR adjusted p-values}
}
\references{ 
David M. Rocke (2004), Design and analysis of experiments with high throughput biological assay data, Seminars in Cell & Developmental Biology, 15, 703-713. 

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{genediff}}}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)
LoggedSmpd0<-neweS(lnorm(log(sample.mat)),vlist)

pvlist<-genediff(LoggedSmpd0)
pvlist$Posterior[1:5,]

apvlist<-pvadjust(pvlist)
names(apvlist)
apvlist$Posterior.FDR[1:5,]


}
\keyword{ models }% at least one, from doc/KEYWORDS

