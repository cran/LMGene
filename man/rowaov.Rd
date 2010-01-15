\name{rowaov}
\alias{rowaov}
\title{Gene by gene ANOVA function}
\description{
  Computes the mean squares and degrees of freedom for gene-by-gene ANOVAs.
}
\usage{
rowaov(eS, model=NULL)
}
\arguments{
  \item{eS}{An \code{ExpressionSet} object.  Any transformation and normalization of \code{eS} should be
done prior to use in \code{rowaov}.}
  \item{model}{Model used for comparison. See details and \code{\link{LMGene}}.}
}
\details{
 If you have data in a \code{matrix} and information about experimental design factors, then you
  can use \code{\link{neweS}} to convert the data into an \code{ExpressionSet} object. Please see
  \code{\link{neweS}} for more detail.

  The \code{model} argument is an optional character string, constructed like the right-hand
  side of a formula for \code{lm}. It specifies which of the variables in the \code{ExpressionSet} will
  be used in the model and whether interaction terms will be included. If \code{model=NULL},
  it uses all variables from the \code{ExpressionSet} without interactions. Be careful of using
  interaction terms with factors; this often leads to overfitting, which will yield an error.
}
\value{
  \item{resmat}{A matrix of MSEs and degrees of freedom for all model factors and all genes.  The first rows of \code{resmat} contain MSE's
for each effect in \code{model}, ending with the residual MSE.  The remaining rows contain degrees of freedom for each effect in the model,
ending with the residual d.f.  Each column corresponds to a gene.}
}
\references{ 
David M. Rocke (2004), Design and analysis of experiments with high throughput biological assay data, 
\emph{Seminars in Cell & Developmental Biology}, \bold{15}, 703--713. 

\url{http://dmrocke.ucdavis.edu} 
}
\author{David Rocke and Geun-Cheol Lee}
\seealso{\code{\link{genediff}}, \code{\link{LMGene}}}
\examples{
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

raw.eS <- neweS(sample.mat, vlist)

# glog transform data
trans.eS <- transeS(raw.eS, lambda = 727, alpha = 56)

# Perform gene-by-gene anova
resmat <- rowaov(trans.eS)
resmat[,1:3]
}
\keyword{ models }
