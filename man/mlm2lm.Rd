\name{mlm2lm}
\alias{mlm2lm}
\title{Linear Model converting function}
\description{
  This function rule out the specified 'lm' class data out of the given 'c("mlm", "lm")' class data.
}
\usage{
mlm2lm(lmobj, i)
}
\arguments{
  \item{lmobj}{An object of class 'c("mlm", "lm")'.}
  \item{i}{A specific number that indicates a 'lm' in \code{lmobj}.}
}
\details{
  In case of multiple response from 'lm' function, this function can used.
}
\value{
  \item{lmobj2}{Selected 'lm' class data.}
}
\references{\url{http://www.idav.ucdavis.edu/~dmrocke/}}
\author{David Rocke and Geun-Cheol Lee}
\seealso{ }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.eS)
Smpd0 <- sample.eS
# model information 
for(i in 1:length(varLabels(Smpd0))) {
  assign(paste('x', i, sep=''), as.factor(pData(Smpd0)[,i]))
}

fchar <- ''
for(i in 1:length(varLabels(Smpd0))) {
  fchar <- paste(fchar, paste('x', i, sep=''), ifelse(i<length(varLabels(Smpd0)), '+', ''), sep='')
}
fchar2 <- paste("y ~",fchar)

# run regression and ANOVAs
y <- t(as.matrix(exprs(Smpd0)))
formobj <- as.formula(fchar2)
tmp <- lm(formobj)
class(tmp)

tmp2 <- mlm2lm(tmp, i)
class(tmp2)
}
\keyword{models}
