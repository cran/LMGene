\name{psmeans}
\alias{psmeans}
\title{Function to take means of probesets.}
\description{
Converts an \code{ExpressionSet} or \code{AffyBatch} object with one row of expression data per probeset into an \code{ExpressionSet} 
or \code{AffyBatch} object with one row per probe.  

}
\usage{
psmeans(eS, ind)
}
\arguments{
  \item{eS}{An \code{ExpressionSet} or \code{AffyBatch} object}
  \item{ind}{ A vector used to indicate which probes go into which probesets.}

}
\details{
Each entry of \code{ind} corresponds to one probe and tells the number of the probeset it belongs to.  See \code{\link{tranestAffyProbeLevel}}  
and \code{\link{sample.ind}} for examples.

}
\value{
Returns an \code{ExpressionSet} or \code{AffyBatch} object with the expression matrix rows corresponding to probesets instead of individual
probes.  Elements of the returned \code{ExpressionSet} or \code{AffyBatch} object are means over each probeset.
}
\author{ John Tillinghast }
\seealso{\code{\link{tranestAffyProbeLevel}}, \code{\link{sample.ind}}}
\examples{
library(LMGene)
library(Biobase)

data(sample.eS)
data(sample.ind)

# glog transform data
trs.eS <- transeS (sample.eS, 667, 65) 

# lowess normalize
ntrs.eS <- lnormeS(trs.eS) 

# take means over probesets
genesample.eS<- psmeans (ntrs.eS, sample.ind)
}

\keyword{manip}

