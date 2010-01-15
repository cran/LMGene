\name{glog}
\alias{glog}
\title{ Generalized log transformation function
 }
\description{
  This function transforms the input values by the generalized log function.
}
\usage{
glog(y, lambda)
}
\arguments{
  \item{y}{ A data matrix}
  \item{lambda}{Transformation parameter}
}
\details{
  The glog transformation of a variable y is defined as \code{log(y + sqrt(y^2 + lambda))}.  Using \code{lambda = 0} corresponds to the log 
transformation, up to a scale factor of 2.  (Other, equivalent expressions exist for the glog transformation.  See Durbin et al. (2002) and 
Huber et al. (2002) for futher details.)

  The input matrix \code{y} may be modified prior to transformation by subtracting a constant or vector ("\code{alpha}").
  The parameters \code{lambda} and \code{alpha} may be estimated from \code{\link{tranest}}.  
}
\value{
  \item{yt }{ A matrix of glog-transformed values}
}
\references{ 
Durbin, B.P., Hardin, J.S., Hawkins, D.M., and Rocke, D.M. (2002) A variance-stabilizing transformation for gene-expression microarray
data, \emph{Bioinformatics}, \bold{18}, S105--S110.

Huber, W., Von Heydebreck, A., Sueltmann, H., Poustka, A., and Vingron, M. (2002) Variance 
stabilization applied to microarray data calibration and to the quantification of differential
expression, \emph{Bioinformatics}, \bold{18}, S96--S104.

\url{http://dmrocke.ucdavis.edu} 
}
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{tranest}}, \code{\link{transeS}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
sample.mat[1:5,1:4]

GloggedSmpd<-glog(sample.mat-50,500)
GloggedSmpd[1:5,1:4]

}
\keyword{ math }

