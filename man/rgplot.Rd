\name{rgplot}
\alias{rgplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting density function
 }
\description{
  This function plotting density distribution of each channel of the given array data.
}
\usage{
rgplot(eS, norm = "lnorm", red = TRUE, green = TRUE, title = "", span = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eS}{ Array data, must be 'exprSet' class }
  \item{norm}{ Data normalization method }
  \item{red}{ TRUE, if Cy5 channels are to be plotted }
  \item{green}{ TRUE, if Cy3 channels are to be plotted }
  \item{title}{ Title of the graph }
  \item{span}{ A parameter for lowess normalization method }
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{arrplot1}}, \code{\link{arrplote}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd0)

#
rgplot(Smpd0)

}
\keyword{ hplot }% at least one, from doc/KEYWORDS
