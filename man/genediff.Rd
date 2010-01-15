\name{genediff}
\alias{genediff}
\title{
Gene-by-gene and posterior p-value calculation function.
}
\description{
  Computes two sets of p-values per gene or probe
  via gene-by-gene ANOVA, using 
both the gene-specific MSE and the posterior MSE for
  each term in the ANOVA. P-values are not adjusted for multiple testing.
  \cr
  Assumes a fixed effects model and that the correct denominator for all comparisons is the MSE.
}
\usage{
genediff(eS, model = NULL, method = c("MLE", "MOM", "MOMlog"), 
verbose = TRUE)
}
\arguments{
  \item{eS}{An \code{ExpressionSet} object.  Any transformation and normalization of \code{exprs(eS)} should be conducted prior to use in \code{genediff}.}
  \item{model}{Model used for comparison; see details and \code{\link{LMGene}}.}
  \item{method}{Method by which posterior p-values are calculated.  Default \code{"MLE"}.}
  \item{verbose}{If \code{TRUE}, the prior degrees of freedom and mean reciprocal precision are printed.  See details.}
}
\details{  The argument \code{eS} must be an \code{ExpressionSet} object from the Biobase package. 
  If you have data in a \code{matrix} and information about experimental design factors, then you
  can use \code{\link{neweS}} to convert the data into an \code{ExpressionSet} object. Please see
  \code{\link{neweS}} for more detail.

  The \code{model} argument is an optional character string, constructed like the right-hand
  side of a formula for \code{lm}. It specifies which of the variables in the \code{ExpressionSet} will
  be used in the model and whether interaction terms will be included. If \code{model=NULL},
  it uses all variables from the \code{ExpressionSet} without interactions. Be careful of using
  interaction terms with factors; this often leads to overfitting, which will yield an error.

The \code{method} argument specifies how the adjusted MSE and degrees of freedom should be calculated for use in computation of the 
posterior p-values:

\describe{
\item{"MLE"}{Default.  Calculate adjusted MSE and degrees of freedom by maximum likelihood estimation, as described in Wright and Simon (2003).}
\item{"MOM"}{Calculate adjusted MSE and degrees of freedom by method of moments, as described in Rocke (2003).}
\item{"MOMlog"}{Calculate adjusted MSE and degrees of freedom by method of moments on log scale, as described in Smyth (2004).  
 Uses functions \code{fitFdist} and \code{trigammainverse} from the package \code{limma}.  Note that the method of Smyth (2004) is  
 used here to calculate the posterior MSE, but not to directly calculate the posterior p-values.}
}

All three methods assume that the gene-specific MSE's follow a gamma distribution with mean tau. (NB: Notation and parameterization vary somewhat 
between each of the source papers.) The mean of the gamma distribution, tau, is modeled with an inverse gamma prior with hyperparameters alpha
and beta.   Empirical Bayes methods are used to estimate the prior hyperparameters, either by maximum likelihood, method of moments, 
or method of moments on the log scale.  The "posterior MSE" is the posterior mean of the variances given the observed gene-specific MSE's.  

If \code{verbose = TRUE}, the function prints the estimated prior degrees of freedom, which equals twice the prior shape parameter alpha, 
and the estimated prior mean reciprocal precision, or 1/(alpha*beta).  

All p-values are calculated from fixed-effects ANOVA F statistics, using either the gene-specific MSE or the posterior MSE as the denominator. 
}
\value{
A list with components:
\item{Gene.Specific}{A matrix of p-values calculated using the gene-specific MSE, with one row for each gene/probe and one column for each factor}
\item{Posterior}{A matrix of p-values calculated using the posterior MSE, with one row for each gene/probe and one column for each factor}
}
\references{
Rocke, D.M. (2003) Heterogeneity of variance in gene expression microarray data. \url{http://dmrocke.ucdavis.edu/papers/empbayes2.pdf}

Smyth, G.K (2004) Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. 
\emph{Statistical Applications in Genetics and Molecular Biology} \bold{3}, Article 3. 
\url{http://www.bepress.com/sagmb/vol3/iss1/art3/}

Wright, G.W. and Simon, R.M. (2003) A random variance model for detection of differential gene expression in small microarray experiments.
\emph{Bioinformatics} \bold{19}, 2448--2455.

\url{http://dmrocke.ucdavis.edu}
}
\author{David Rocke, Geun-Cheol Lee, and Blythe Durbin-Johnson}

\seealso{\code{\link{LMGene}}, \code{\link{rowaov}}, \code{\link{neweS}}}
\examples{
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

raw.eS <- neweS(sample.mat, vlist)

# glog transform data
trans.eS <- transeS(raw.eS, lambda = 727, alpha = 56)

# calculate p-values
pvlist <- genediff(trans.eS)
pvlist$Posterior[1:5,]
}

\keyword{models}

