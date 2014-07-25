\name{RRunif}
\alias{RRunif}
\title{Random scaling used with balls}
\description{
 The model refers to the d-dimensional univariate distribution on a rectangular
 window.
  }
 \usage{
RRunif(min, max, normed) 
}

\arguments{
  \item{min,max}{lower and upper corner of a rectangular window}
  \item{normed}{logical with default value \code{TRUE}.

    Advanced. If \code{FALSE} then the identicator function
    for the window is not normed to get a probability
    distribution. Nonetheless random drawing from the distribution
    still works.
  }
}

\details{
  In the one-dimensional case it has the same effect as
  \code{\link{RRdistr}(\link[=runif]{unif}(min=min, max=max, log=log))}
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
 \command{\link{RRspheric}},
 }

 
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
## uniform distribution on [0,1] x [-2, -1]
RFrdistr(RRunif(c(0, -2), c(2, -1)), n=5, dim=2)
RFpdistr(RRunif(c(0, -2), c(2, -1)), q=c(1, -1.5), dim=2)
RFddistr(RRunif(c(0, -2), c(2, -1)), x=c(1, -1.5), dim=2)
\dontshow{FinalizeExample()}
}

 
\keyword{spatial}
\keyword{models}




