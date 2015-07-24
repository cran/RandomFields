\name{RMball}
\alias{RMball}
\title{RMball}
\description{
  \command{RMball} referes to the indicator function of a ball
 with radius 1.
}
\usage{
RMball(var, scale, Aniso, proj)
}

\arguments{
 \item{var, scale, Aniso, proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}

\seealso{
  \command{\link{RMpolygon}}
  \command{\link{RMspheric}},
 \command{\link{RFsimulate}},
 \command{\link{RMmodel}}.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\keyword{spatial}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
\dontshow{\dontrun{
model <- RPcoins(shape=RMball())
x <- seq(0, 10, 0.02)
plot(RFsimulate(model, x=x))
}}
\dontshow{FinalizeExample()}
}