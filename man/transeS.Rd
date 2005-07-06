\name{transeS}
\alias{transeS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to apply the glog transform to an expression set. Returns the transformed expression set (not normalized). }
\description{
  For each element in the array of expression data, this applies the glog transform
y -> glog (y-alpha, lambda).
  If alpha is a vector, it must have one entry per sample, and transeS will use the appropriate entry from the vector.
}
\usage{
transeS(eS, lambda, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ An unprocessed expression set. }
  \item{lambda}{ The parameter lambda to be used in the glog transform (Durbin and Rocke 2003). }
  \item{alpha}{ The alpha parameter(s) for the glog transform. May be a single number used for all samples, or a vector with one entry per sample.}
}
\details{

}
\value{
Returns an expression set with the same vlist as eS, but the matrix is now glog-transformed. That matrix can be normalized with \code{\link{norm}} or \code{\link{lnorm}}.
}
\references{ }
\author{ John Tillinghast }
\note{ 
}
\seealso{ }
\examples{

data(sample.eS)
transeS (sample.eS, 667, 65) -> trsample.eS
}

\keyword{ math }% at least one, from doc/KEYWORDS

