\name{GetLMObj}
\alias{GetLMObj}
\title{Function to get a simple lm object for a regression on the relevant model.}
\description{
  Internal to routines. Primarily used to get the X matrix corresponding to the model given (or the default model for the eS). Typically this is used to find residuals efficiently.
}
\usage{
GetLMObj(eS, model=NULL)
}
\arguments{
  \item{eS}{An unprocessed \code{ExpressionSet} object.}
  \item{model}{Model used in the regression. Uses only variables from \code{pData(eS)}.}
}
\details{
}
\value{
Returns an lm object than corresponds to regressing one probe from the eS on the model specified (or the default model). See \code{\link{lm}}.}
\references{ }
\author{John Tillinghast}
\note{
}
\seealso{ }
\examples{
data(sample.eS)
lmod <- GetLMObj (sample.eS)
X <- lmod$x
}
\keyword{math}
