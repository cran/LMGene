\name{rowlist}
\alias{rowlist}
\title{Gene name listing function}
\description{
  This function makes significant gene list for a specified factor, where genes are selected as significant by the given p-values and significance level.
}
\usage{
rowlist(genemat, effnum, apvlist, level, posterior = TRUE)
}
\arguments{
  \item{genemat}{A matrix data of array.}
  \item{effnum}{Factor number.}
  \item{apvlist}{A vector with FDR adjusted p-value.}
  \item{level}{Significance level.}
  \item{posterior}{TRUE, if adjusted p-values are to be computed with Posterior method.}
}
\details{
  \code{genemat} is an n-by-p matrix of expression values. 
  \code{effnum} is the column number for the effect of interest.
  \code{apvlist} is a matrix of p-values from pvadjust or genediff
  the routine returns a list of genes whose FDR p-value is less
  than level using either individual gene or posterior MSE's.
  This function returns gene names if \code{rownames(genemat)} is
  not \code{NULL}, or gene numbers otherwise.
  \cr
  \code{level} indicates False Discovery Rate. e.g.) level 0.05 means 5% FDR.
}
\value{
  \item{genelist }{A vector containing gene names if \code{rownames(genemat)} is
  not \code{NULL}, or gene numbers otherwise.}
}
\references{ 
David M. Rocke (2004), Design and analysis of experiments with high throughput biological assay data, Seminars in Cell & Developmental Biology, 15, 703-713. 

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{David Rocke and Geun-Cheol Lee}
\seealso{\code{\link{LMGene}}, \code{\link{rowaov}}}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)
LoggedSmpd0<-neweS(lnorm(log(sample.mat)), vlist)

pvlist <- genediff(LoggedSmpd0)
apvlist <- pvadjust(pvlist)

genelist <- rowlist(exprs(LoggedSmpd0), 2, apvlist, 0.01)
genelist
}
\keyword{manip}
