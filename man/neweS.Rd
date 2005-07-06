% --- Source file: neweS.Rd ---
\name{neweS}
\alias{neweS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coercing to ExpressionSet function 
 }
\description{
  This function converts a matrix data and its experimental data into an object of 'ExpressionSet' class.
}
\usage{
neweS(mat, vlist, annotat, vlabel = as.list(names(vlist)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{ A matrix data to be converted whose column names must be unique }
  \item{vlist}{ A list which contains several factors representing the experiment description }
  \item{annotat}{ An optional character string identifying the type of chips used in the experiment }
  \item{vlabel}{ A list of labels for the variables represented by the columns of \code{pData} of the 'ExpressionSet' object to be made }
}
\details{
  Must load Biobase package first for converting. \cr
  \code{vlist} contains all the considered factors in which level values of each element represent the corresponding column of \code{mat}.
}
\value{
  \item{eset }{The converted object of 'ExpressionSet' class}
}
\references{ \url{http://dmrocke.ucdavis.edu/} }
\author{ David Rocke and Geun-Cheol Lee }
\seealso{ \code{\link[methods]{as}} }
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(sample.mat)
data(vlist)

Smpdt<-neweS(sample.mat,vlist)

data(sample.eS)
identical(sample.eS, Smpdt)

}
\keyword{ manip }% at least one, from doc/KEYWORDS

