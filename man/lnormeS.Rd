\name{lnormeS}
\alias{lnormeS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to apply lowessnorm to a transformed expression set. Returns the normalized expression set. }
\description{
Basically the same as \code{\link{lnorm}}, but it applies to, and returns, expression sets instead of matrices.

}
\usage{
lnormeS(eS, span=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ A transformed expression set. }
  \item{span}{A parameter for lowess. }
}
\details{

}
\value{
Returns an expression set with the same vlist as eS, but the matrix has been normalized by \code{\link{lnorm}}.
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/}  }
\author{ John Tillinghast }
\note{ 
}
\seealso{ \code{\link{lnorm}}, \code{\link{norm}}}
\examples{

data(sample.eS)
transeS (sample.eS, 667, 65) -> trsample.eS
lnormeS (trsample.eS) -> normtrsample.eS
}

\keyword{ math }% at least one, from doc/KEYWORDS

