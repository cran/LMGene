\name{msa}
\alias{msa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Relative mean square calculation function
 }
\description{
  Calculate the relative mean squre values
}
\usage{
msa(v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{ A vetor containing mean square values of all the factors}
}
\value{
  \item{rv }{relative mean squre values for all factors.}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ }
\examples{
#library
library(Biobase)
library(LMGene)

#data
#data
data(sample.eS)
Smpd0 = sample.eS
# model information 
for(i in 1:length(Smpd0@phenoData@varLabels)){
  assign(paste('x', i, sep=''),as.factor(Smpd0@phenoData@pData[,i]))
}
  
fchar=''
for(i in 1:length(Smpd0@phenoData@varLabels)){
  fchar=paste(fchar, paste('x', i, sep=''), ifelse(i<length(Smpd0@phenoData@varLabels), '+', ''), sep='')
}
fchar2 <- paste("y ~",fchar)
#
# run regression and anovas
y <- t(as.matrix(Smpd0@exprs))
formobj <- as.formula(fchar2)
tmp <- lm(formobj)
tmp2 <- mlm2lm(tmp,i)
tmp3 <- anova(tmp2)$Mean
tmp4 <- msa(tmp3)
rbind(tmp3, tmp4)

}
\keyword{ math }% at least one, from doc/KEYWORDS

