\name{mlm2lm}
\alias{mlm2lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model converting function
 }
\description{
  This function rule out the specified 'lm' class data out of the given 'c("mlm", "lm")' class data.
}
\usage{
mlm2lm(lmobj, i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmobj}{ An object of class 'c("mlm", "lm")' }
  \item{i}{ A specific number that indicates a 'lm' in \code{lmobj}  }
}
\details{
  In case of multiple response from 'lm' function, this function can used.
}
\value{
  \item{lmobj2 }{Selected 'lm' class data}
}
\references{ \url{http://www.idav.ucdavis.edu/~dmrocke/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link{arrplot1}}, \code{\link{arrplotd}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd0)

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
class(tmp)

tmp2 <- mlm2lm(tmp,i)
class(tmp2)


}
\keyword{ models }% at least one, from doc/KEYWORDS

