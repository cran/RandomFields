\name{RRloc}
\alias{RRloc}
\title{Random scaling used with balls}
\description{

 \command{RRloc} modifies location and scale of a distribution.
 
 }
 \usage{
RRloc(phi, mu, scale) 
}

\arguments{
  \item{phi}{distribution \command{\link{RMmodel}}.}
  \item{mu}{location shift}
 \item{scale}{scale modification}
}

\value{
 \command{\link{RRloc}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

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
set.seed(0)
hist(RFrdistr(RRspherical(balldim=2), n=1000), 50)
hist(RFrdistr(RRloc(mu=3, RRspherical(balldim=2)), n=1000), 50)
}

\keyword{spatial}
\keyword{models}




