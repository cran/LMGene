\name{tranest}
\alias{tranest}
\title{Glog transformation parameter estimation function}
\description{
Estimates parameters for the glog transformation, by maximum likelihood or by minimizing
the stability score.
}
\usage{
tranest(eS, ngenes = -1, starting = FALSE, lambda = 1000, alpha = 0,
    gradtol = 1e-3, lowessnorm = FALSE, method=1, mult=FALSE, model=NULL, 
	SD = FALSE, rank = TRUE, model.based = TRUE, rep.arrays = NULL)
}
\arguments{
  \item{eS}{
An \code{ExpressionSet} object
}
  \item{ngenes}{Number of genes to be used in parameter estimation.  Default is to use all genes
unless there are more than 100,000, in which case a subset of 50,000 genes is selected at random.}
  \item{starting}{If \code{TRUE}, user-specified starting values for \code{lambda} and \code{alpha} are input to 
the optimization routine}
  \item{lambda}{Starting value for parameter \code{lambda}. Ignored unless \code{starting = TRUE}}
  \item{alpha}{Starting value for parameter \code{alpha}. Ignored unless \code{starting = TRUE}}
  \item{gradtol}{A positive scalar giving the tolerance at which the scaled 
gradient is considered close enough to zero to terminate the algorithm
}
  \item{lowessnorm}{If \code{TRUE}, lowess normalization (using \code{\link{lnorm}}) is used in calculating 
the likelihood.}
  \item{method}{Determines optimization method. Default is 1, 
which corresponds to a Newton-type method (see \code{nlm} and details.)}
  \item{mult}{If \code{TRUE}, tranest will use a vector alpha with one (possibly different) entry per sample. 
Default is to use same alpha for every sample.  \code{SD} and \code{mult} may not both be \code{TRUE}.}
  \item{model}{Specifies model to be used. Default is to use all variables from eS without interactions. See details.}
  \item{SD}{If \code{TRUE}, transformation parameters are estimated by minimizing the stability score rather than by maximum likelihood.  
	See details.}
  \item{rank}{If \code{TRUE}, the stability score is calculated by regressing the replicate standard deviations
on the ranks of the gene/row means (rather than on the means themselves).  Ignored unless \code{SD = TRUE}}
  \item{model.based}{If \code{TRUE}, the stability score is calculated using the standard deviations of residuals from the linear
    model in \code{model}.  Ignored unless \code{SD = TRUE}}
  \item{rep.arrays}{List of sets of replicate arrays.  Each element of \code{rep.arrays} should be a vector with entries
corresponding to arrays (columns) in \code{exprs(eS)} conducted under the same experimental conditions, i.e., with identical
rows in \code{pData(eS)}. Ignored unless \code{SD = TRUE} and \code{model.based = FALSE}}
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

  The default estimation method is maximum likelihood.  The likelihood is derived by assuming that there exist values for \code{lambda}
  and \code{alpha} such that the residuals from the linear model in \code{model}, fit to glog-transformed data using those values 
  for \code{lambda} and \code{alpha}, follow a normal distribution.  See Durbin and Rocke (2003) for details.

  If \code{SD = TRUE}, \code{lambda} and \code{alpha} are estimated by minimizing the stability score rather than by maximum likelihood.
  The stability score is defined as the absolute value of the slope coefficient from the regression of the replicate/residual 
  standard deviation on the gene/row means, or on the rank of the gene/row means.  If \code{model.based = TRUE}, the stability
  score is calculated using the standard deviation of residuals from the linear model in \code{model}.  Otherwise, the stability score is 
  calculated using the pooled standard deviation over sets of replicates in \code{rep.arrays}. See Wu and Rocke (2009) for details. 

Optimization methods in \code{method} are as follows:
\describe{
\item{1 = }{Newton-type method, using \code{nlm}}
\item{2 = }{Nelder-Mead, using \code{optim}}
\item{3 = }{BFGS, using \code{optim}}
\item{4 = }{Conjugate gradients, using \code{optim}}
\item{5 = }{Simulated annealing, using \code{optim} (may only be used when \code{mult = TRUE})}
}
}
\value{
A list with components: 
\item{lambda}{Estimate of transformation parameter lambda}
\item{alpha}{Estimate of transformation parameter alpha}
}
\references{
Durbin, B.P and Rocke, D.M. (2003) Estimation of Transformation Parameters for Microarray Data,  
\emph{Bioinformatics}, \bold{19}, 1360--1367.

Wu, S. and Rocke, D.M. (2009) Analysis of Illumina BeadArray data using variance stabilizing transformations.

\url{http://dmrocke.ucdavis.edu} 

}
\author{David Rocke, Geun-Cheol Lee, John Tillinghast, Blythe Durbin-Johnson, and Shiquan Wu}
\seealso{
\code{\link{tranestAffyProbeLevel}}, \code{\link{lnorm}}, \code{\link{glog}}}
\examples{
library(Biobase)
library(LMGene)

#data
data(sample.eS)

tranpar <- tranest(sample.eS, 100)
tranpar
tranpar <- tranest(sample.eS, mult=TRUE)
tranpar
}
\keyword{math}

