\name{transeS}
\alias{transeS}
\title{ Function to apply the glog transform to an expression set.}
\description{
  For each element in the array of expression data, this function applies the glog transform
y -> glog (y-alpha, lambda).
  If alpha is a vector, it must have one element for each column in \code{exprs(eS)}.
}
\usage{
transeS(eS, lambda, alpha)
}
\arguments{
  \item{eS}{ An \code{ExpressionSet} or \code{AffyBatch} object}
  \item{lambda}{ The parameter lambda to be used in the glog transform.}
  \item{alpha}{ The alpha parameter(s) for the glog transform. May be a single number used for all samples, or a vector with one entry per sample.}
}

\details{
The glog transformation of a variable y is defined as \code{log(y + sqrt(y^2 + lambda))}.  Using \code{lambda = 0} corresponds to the log 
transformation, up to a scale factor of 2.  (Other, equivalent expressions exist for the glog transformation.  See Durbin et al. (2002) and 
Huber et al. (2002) for futher details.)

\code{transeS} subtracts a (scalar or vector) parameter \code{alpha} prior to application of the glog transformation, resulting in the expression
\code{log(y - alpha + sqrt((y - alpha)^2 + lambda))}.

The parameters \code{lambda} and \code{alpha} may be estimated using \code{\link{tranest}}.
}

\value{
Returns an \code{ExpressionSet} or \code{AffyBatch} object with the expression matrix glog-transformed. 
}
\references{Durbin, B.P., Hardin, J.S., Hawkins, D.M., and Rocke, D.M. (2002) A variance-stabilizing transformation for gene-expression microarray
data, \emph{Bioinformatics}, \bold{18}, S105--S110.

Huber, W., Von Heydebreck, A., Sueltmann, H., Poustka, A., and Vingron, M. (2002) Variance 
stabilization applied to microarray data calibration and to the quantification of differential
expression, \emph{Bioinformatics}, \bold{18}, S96--S104.

\url{http://dmrocke.ucdavis.edu} }
\author{ John Tillinghast }
\seealso{\code{\link{glog}}, \code{\link{tranest}}}
\examples{
library(LMGene)
library(Biobase)

data(sample.eS)
trsample.eS <- transeS (sample.eS, 667, 65) 
}

\keyword{ math }

