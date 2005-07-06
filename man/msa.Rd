\name{msa}
\alias{msa}
\title{Relative mean square calculation function}
\description{
  Calculate the relative mean square values.
}
\usage{
msa(v)
}
\arguments{
  \item{v}{A vector containing mean square values of all the factors.}
}
\value{
  \item{rv}{relative mean square values for all factors.}
}
\references{\url{http://www.idav.ucdavis.edu/~dmrocke/}}
\author{David Rocke and Geun-Cheol Lee}
\seealso{ }
\examples{
#library
library(Biobase)
library(LMGene)

#data
#data
data(sample.eS)
Smpd0 <- sample.eS
# model information 
for(i in 1:length(varLabels(Smpd0))){
  assign(paste('x', i, sep=''), as.factor(pData(Smpd0)[,i]))
}
  
fchar <- ''
for(i in 1:length(varLabels(Smpd0))){
  fchar <- paste(fchar, paste('x', i, sep=''), ifelse(i<length(varLabels(Smpd0)), '+', ''), sep='')
}
fchar2 <- paste("y ~", fchar)

# run regression and ANOVAs
y <- t(as.matrix(exprs(Smpd0)))
formobj <- as.formula(fchar2)
tmp <- lm(formobj)
tmp2 <- mlm2lm(tmp, i)
tmp3 <- anova(tmp2)$Mean
tmp4 <- msa(tmp3)
rbind(tmp3, tmp4)
}
\keyword{math}
