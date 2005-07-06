\name{GetLMObj}
\alias{GetLMObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to get a simple lm object for a regression on the relevant model. }
\description{
  Internal to routines. Primarily used to get the X matrix corresponding to the model given (or the default model for the eS). Typically this is used to find residuals efficiently.
}
\usage{
GetLMObj(eS, model=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ An unprocessed expression set. }
  \item{model}{ Model used in the regression. Uses only variables from eS@phenoData@pData. }
}
\details{

}
\value{
Returns an lm object than corresponds to regressing one probe from the eS on the model specified (or the default model). See \code{\link{lm}}.}
\references{ }
\author{ John Tillinghast }
\note{ 
}
\seealso{ }
\examples{

data(sample.eS)
lmod <- GetLMObj (sample.eS)
X <- lmod$x
}

\keyword{ math }% at least one, from doc/KEYWORDS

