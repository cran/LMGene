\name{randeS}
\alias{randeS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to pick a random sample of probes from an expression set. }
\description{
Takes an expression set eS and a number of probes n. Returns a random subset of eS with n probes.

}
\usage{
randeS(eS, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ An expression set. }
  \item{n}{ A positive integer smaller than the number of probes in eS. }

}
\details{


}
\value{
Returns an expression set with the same vlist as eS, but only n probes chosen at random.
}
\references{ }
\author{ John Tillinghast }
\note{ 
}
\seealso{ }
\examples{

data(sample.eS)
randeS (sample.eS, 100) -> sample.100.eS

}

\keyword{ math }% at least one, from doc/KEYWORDS

