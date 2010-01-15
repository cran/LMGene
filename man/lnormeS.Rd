\name{lnormeS}
\alias{lnormeS}
\title{Function to apply lowess normalization to an expression set.}
\description{
Like \code{\link{lnorm}}, but applies to and returns an \code{ExpressionSet} or \code{AffyBatch} object instead of a \code{matrix}.
}
\usage{
lnormeS(eS, span=0.1)
}
\arguments{
  \item{eS}{ An \code{ExpressionSet} or \code{AffyBatch} object}
  \item{span}{Smoothing parameter for lowess.  Larger values correspond to more smoothness.}
}
\value{
Returns an \code{ExpressionSet} with \code{exprs(eS)} normalized by \code{\link{lnorm}}.
}
\references{ \url{http://dmrocke.ucdavis.edu}  }
\author{ John Tillinghast, Blythe Durbin-Johnson}
\seealso{ \code{\link{lnorm}}, \code{\link{norm}}}
\examples{
library(LMGene)
library(Biobase)

data(sample.eS)

# glog transform expression set
trsample.eS <- transeS (sample.eS, 667, 65)

# normalize expression set
normtrsample.eS <- lnormeS (trsample.eS) 
}

\keyword{ smooth }

