\name{LMGene}
\alias{LMGene}
\title{LMGene main function}
\description{
LMGene calls function \code{\link{genediff}} to calculate the unadjusted gene-specific and posterior p-values of all genes
and then calculates the FDR-adjusted p-values of all genes.
Significant genes for each factor in \code{model} (based on either the gene-specific or posterior FDR-adjusted p-values) are output.

}
\usage{
LMGene(eS, model = NULL, level = 0.05, posterior = FALSE, 
method = c("MLE", "MOM", "MOMlog"))
}
\arguments{
  \item{eS}{An \code{ExpressionSet} object.  Any transformation and normalization of \code{exprs(eS)} should be conducted prior 
to use in \code{LMGene}.}
  \item{model}{Specifies model to be used. Default is to use all variables from eS without interactions. See details.}
  \item{level}{Significance level}
  \item{posterior}{If \code{TRUE}, the posterior FDR-adjusted p-values are used in listing significant genes for each factor.  Default
is to use gene-specific FDR-adjusted p-values.}
  \item{method}{Method by which the posterior p-values are calculated.  Default is \code{"MLE"}.}
}
\details{
  If you have data in a \code{matrix} and information about experimental design factors, then you
  can use \code{\link{neweS}} to convert the data into an \code{ExpressionSet} object. Please see
  \code{\link{neweS}} for more detail.

The \code{level} argument indicates the False Discovery Rate, e.g. level=0.05 means a 5 percent FDR.

  The \code{model} argument is an optional character string, constructed like the right-hand
  side of a formula for \code{lm}. It specifies which of the variables in the \code{ExpressionSet} will
  be used in the model and whether interaction terms will be included. If \code{model=NULL},
  it uses all variables from the \code{ExpressionSet} without interactions. Be careful of using
  interaction terms with factors; this often leads to overfitting, which will yield an error.

See \code{\link{genediff}} for details of \code{method}.
}
\value{
\item{lmres}{A list with one component for each factor in \code{model}. Each component consists of a character vector with one element
per significant gene.  If no genes are significant for a given factor, the component for that factor is set to \code{"No significant genes"}.}  
}
\references{
Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical and powerful
approach to multiple testing, \emph{Journal of the Royal Statistical Society, Series B}, \bold{57},
289--300.

David M. Rocke (2004) Design and analysis of experiments with high throughput biological assay data, 
\emph{Seminars in Cell & Developmental Biology}, \bold{15}, 703--713. 

\url{http://dmrocke.ucdavis.edu}
}
\author{David Rocke and Geun-Cheol Lee}
\seealso{\code{\link{genediff}}, \code{\link{neweS}}}
\examples{
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

raw.eS <- neweS(sample.mat, vlist)

# glog transform data
trans.eS <- transeS(raw.eS, lambda = 727, alpha = 56)

# Identify significant genes, using an FDR of 1 percent
LMGene(trans.eS, level = 0.01)
}

