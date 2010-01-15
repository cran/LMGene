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
  Length is equal to the number of probes (rows) in sample.mat.

}
\examples{
     data(sample.eS)
     data(sample.ind)
     trs.eS <- transeS (sample.eS, 667, 65)
     ntrs.eS <- lnormeS(trs.eS) 
     genesample.eS <- psmeans (ntrs.eS, sample.ind) 
}
\keyword{datasets}
