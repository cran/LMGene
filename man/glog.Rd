\name{glog}
\alias{glog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized log transformation function
 }
\description{
  This function transforms the input values with the generalized log function.
}
\usage{
glog(y, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ A matrix data }
  \item{lambda}{ Parameter that should be determined }
}
\details{
  Usually, matrix \code{y} is a modified matrix from an original matrix, after deducting parameter alpha.
  Lambda is one of the parameters that should be determined when using the glog function,
  and these parameters are decided by using the function \code{\link{tranest}}
}
\value{
  \item{yt }{ A matrix containing a transformmed values by glog}
}
\references{ 
Durbin, B.P., Hardin, J.S., Hawkins, D.M., and Rocke, D.M. (2002) A variance-stabilizing transformation for gene-expression microarray
data, Bioinformatics, 18, S105--S110.

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{tranest}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd)
Smpd[1:5,1:4]

GloggedSmpd<-glog(Smpd-50,500)
GloggedSmpd[1:5,1:4]

}
\keyword{ math }% at least one, from doc/KEYWORDS

