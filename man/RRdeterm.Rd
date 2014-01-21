\name{RRdeterm}
\alias{RRdeterm}
\title{Random scaling used with balls}
\description{
 \command{RRdeterm} refers to the determinstic distribution.
 }
 \usage{
RRdeterm(mean) 
}

\arguments{
 \item{mean}{the deterministic value}
}

\value{
 \command{\link{RRdeterm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\details{
  It has the same effect as
  \code{\link{RRdistr}(\link{norm}(mu=mu, sd=sd, log=log))}
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RRdistr}},
\command{\link{RRgauss}},
 }

 

\examples{
RFoptions(seed=0)
x <- seq(-2, 2, 0.001)
p <- RFpdistr(RRdeterm(mean=1), q=x)
plot(x, p, type="l")
\dontshow{RFoptions(seed=NA)}
}

 
\keyword{spatial}
\keyword{models}




