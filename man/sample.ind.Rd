\name{sample.ind}
\alias{sample.ind}
\docType{data}
\title{ Sample probeset index vector}
\description{
  Vector indicating which probeset each probe belongs to 
}
\usage{data(sample.ind)}
\format{
  A vector of integers, e.g., c(1,1,1,2,2,3,3,3,4,4,...).
  Length is of course equal to the number of probes (rows) in sample.mat.

}
\details{
  
}
\examples{
     data(sample.eS)
     data(sample.ind)
     transeS (sample.eS, 667, 65) -> trs.eS
     lnormeS(trs.eS) -> ntrs.eS
     psmeans (ntrs.eS, sample.ind) -> genesample.eS
}
\keyword{datasets}
