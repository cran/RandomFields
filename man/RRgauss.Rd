\name{RRgauss}
\alias{RRgauss}
\title{Random scaling used with balls}
\description{
  \command{RRgauss} defines
  the d-dimensional vector of independent Gaussian random
 variables. 
 }
 \usage{
RRgauss(mu, sd, log) 
}

\arguments{
 \item{mu, sd, log}{see \link[stats]{Normal}. Here the components can be
 vectors, leading to multivariate distibution with independent components}
}

\value{
 \command{\link{RRgauss}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

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
 \command{\link{RRunif}},
 }

\examples{
RFoptions(seed=0)
r <- RFrdistr(RRgauss(mu=c(1,5)), n=1000, dim=2)
plot(r[1,], r[2, ])
\dontshow{RFoptions(seed=NA)}
}

 
\keyword{spatial}
\keyword{models}




