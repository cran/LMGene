\name{Smpd0}
\alias{Smpd0}
\docType{data}
\title{ Sample array data for LMGene}
\description{
  'exprSet' class data 
}
\usage{data(Smpd0)}
\format{
  The format is:
Formal class 'exprSet' [package "Biobase"] with 6 slots
  ..@ exprs      : int [1:500, 1:16] 474 2156 166 246 344 269 249 127 136 147 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "570" "27" "42" "155" ...
  .. .. ..$ : chr [1:16] "X1" "X2" "X3" "X4" ...
  ..@ se.exprs   : logi[0 , 0 ] 
  ..@ phenoData  :Formal class 'phenoData' [package "Biobase"] with 6 slots
  .. .. ..@ : Factor w/ 2 levels "R","G": 1 2 1 2 1 2 1 2 1 2 ...
  .. .. ..@ : Factor w/ 8 levels "1","2","3","4",..: 1 1 2 2 3 3 4 4 5 5 ...
  .. .. ..@ : Factor w/ 2 levels "P","S": 1 2 2 1 1 2 2 1 1 2 ...
  .. .. ..@ : chr "dye"
  .. .. ..@ : chr "slide"
  .. .. ..@ : chr "treat"
  .. .. ..- attr(*, "dim")= int [1:2] 3 2
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:3] "dye" "slide" "treat"
  .. .. .. ..$ : chr [1:2] "pData" "varLabels"
  ..@ description:Formal class 'MIAME' [package "Biobase"] with 11 slots
  .. .. ..@ name          : chr ""
  .. .. ..@ lab           : chr ""
  .. .. ..@ contact       : chr ""
  .. .. ..@ title         : chr ""
  .. .. ..@ abstract      : chr ""
  .. .. ..@ url           : chr ""
  .. .. ..@ samples       : list()
  .. .. ..@ hybridizations: list()
  .. .. ..@ normControls  : list()
  .. .. ..@ preprocessing : list()
  .. .. ..@ other         : list()
  ..@ annotation : chr ""
  ..@ notes      : chr ""
}
\details{
  identical with 'neweS(Smpd,vlist)'
}
\examples{
#library
library(Biobase)
library(LMGene)

#data
data(Smpd)
data(vlist)

Smpdt<-neweS(Smpd,vlist)

data(Smpd0)
identical(Smpd0, Smpdt)
}
\keyword{datasets}
