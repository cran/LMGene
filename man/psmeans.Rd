\name{psmeans}
\alias{psmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to take means of probesets. }
\description{
This is used to estimate expression levels of genes based on the measurements for the relevant probes.

}
\usage{
psmeans(eS, ind)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ A transformed, normalized expression set. }
  \item{ind}{ A vector used to indicate which probes go into which probesets. }

}
\details{

The vector ind has form like c(1,1,1,2,2,2,2,3,3,4,4,4,...)
Each entry corresponds to one probe and tells the number of the probeset it belongs to.

}
\value{
Returns an expression set with the same vlist as eS, but the matrix rows now correspond to probesets instead of individual probes.
}
\references{ }
\author{ John Tillinghast }
\note{ 
}
\seealso{ }
\examples{

data(sample.eS)
data(sample.ind)
transeS (sample.eS, 667, 65) -> trs.eS
lnormeS(trs.eS) -> ntrs.eS
psmeans (ntrs.eS, sample.ind) -> genesample.eS

}

\keyword{ math }% at least one, from doc/KEYWORDS

