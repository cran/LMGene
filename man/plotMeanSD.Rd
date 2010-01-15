\name{plotMeanSD}
\alias{plotMeanSD}
\title{
Plotting function for gene means and standard deviations
}
\description{
Plots the row standard deviation of a matrix of expression data against the row mean, or the rank of the row mean.
}
\usage{
plotMeanSD(indata, by.rank = TRUE, line = FALSE, ymax = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indata}{
An object of class \code{matrix}, \code{data.frame}, \code{ExpressionSet}, or \code{AffyBatch}
}
  \item{by.rank}{
If \code{TRUE}, the row standard deviations are plotted against the ranks of the row means.  Otherwise, the row standard deviations are plotted
against the row means themselves.
}
  \item{line}{
If \code{TRUE}, a lowess smoother line is drawn on the plot.
}
  \item{ymax}{
The upper limit for the plot y-axis.  If missing, axis limits are generated automatically by \code{plot}.
}
}
\details{
Generates a scatter plot of the row standard deviations of a matrix of expression data against the row means or ranks of the row means.
}
\value{
NULL
}
\author{
Rachel Chen and Blythe Durbin-Johnson
}
\examples{
library(LMGene)
library(Biobase)

data(sample.eS)
# transform data
trans.eS <- transeS(sample.eS, lambda = 727, alpha = 56)

# plot SD against rank of mean
plotMeanSD(trans.eS, line = TRUE)
plotMeanSD(sample.eS, line = TRUE, ymax = 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}

