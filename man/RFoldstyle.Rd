\name{RFoldstyle}
\alias{RFoldstyle}
\title{RFoldstyle}
\description{
 This functions is written only for package writers who have based
 their code on RandomFields version 2.

 It avoids warnings if the old style is used, and sets
 \code{spConform = FALSE}.
 
 }
\usage{
RFoldstyle(old=TRUE) 
}
\arguments{
 \item{old}{logical}
}
\value{
 \code{NULL}
}

\note{
  see \sQuote{doc/version2.pdf} for partial documentation of the former
  functions. Note that the documentation is not maintained.
}

\author{
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
  RFoldstyle()
\dontshow{RFoldstyle(FALSE)}
\dontshow{FinalizeExample()}
}

\keyword{spatial}


