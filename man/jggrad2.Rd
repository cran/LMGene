\name{jggrad2}
\alias{jggrad2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generating Jacobian-corrected data
 }
\description{
  This function returns a Jacobian-corrected data with the given parameters lambda and alpha.
}
\usage{
jggrad2(y, lambda, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ A matrix data containing array information }
  \item{lambda}{ A parameter for glog transformation }
  \item{alpha}{ A parameter for glog transformation }
}
\details{
  The input arguments here would be rarely dealt by users directly. 
}
\value{
  \item{data_matrix }{ A matrix containing Jacobian-corrected data, gradient data by lambda and gradient data by alpha}
}
\references{ 
B. Durbin and D.M. Rocke, (2003) Estimation of Transformation Parameters for Microarray Data,  Bioinformatics, 19, 1360-1367.

\url{http://www.idav.ucdavis.edu/~dmrocke/} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{msecalc}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd)
dim(Smpd)

JCSmpd<-jggrad2(Smpd, 500, 50)
dim(JCSmpd)

}
\keyword{ math }% at least one, from doc/KEYWORDS

