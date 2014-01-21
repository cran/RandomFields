\name{RFoldstyle}
\alias{RFoldstyle}
\title{RFoldstyle}
\description{
 This functions is written only for packages writers who have based
 their code on RandomFields version 2.

 It avoids warnings if the old style is used.
 
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
 \url{http://ms.math.uni-mannheim.de}
}

\examples{
RFoptions(seed=0)
  RFoldstyle()
\dontshow{RFoldstyle(FALSE)}
\dontshow{RFoptions(seed=NA)}
}

\keyword{spatial}


