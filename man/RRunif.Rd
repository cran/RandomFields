\name{RRunif}
\alias{RRunif}
\title{Random scaling used with balls}
\description{
 The model refers to the d-dimensional univariate distribution on a rectangular
 window.
  }
 \usage{
RRunif(min, max) 
}

\arguments{
 \item{min,max}{lower and upper corner of a rectangular window}
}

\value{
 \command{\link{RRunif}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RRdistr}},
 \command{\link{RRgauss}},
 \command{\link{RRspherical}},
 }

 
\examples{
RFoptions(seed=0)
## uniform distribution on [0,1] x [-2, -1]
RFrdistr(RRunif(c(0, -2), c(2, -1)), n=5, dim=2)
RFpdistr(RRunif(c(0, -2), c(2, -1)), q=c(1, -1.5), dim=2)
RFddistr(RRunif(c(0, -2), c(2, -1)), x=c(1, -1.5), dim=2)
\dontshow{RFoptions(seed=NA)}
}

 
\keyword{spatial}
\keyword{models}




