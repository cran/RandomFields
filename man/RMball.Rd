\name{RMball}
\alias{RMball}
\title{RMball}
\description{
  \command{RMball} referes to the indicator function of a ball
 with radius 1.
}
\usage{
RMball()
}
\seealso{
  \command{\link{RMpolygon}}
  \command{\link{RMspheric}},
 \command{\link{RFsimulate}},
 \command{\link{RMmodel}}.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}
\examples{
RFoptions(seed=0)
\dontshow{\dontrun{
model <- RPcoins(shape=RMball())
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(RFsimulate(model, x=x))
}}
\dontshow{RFoptions(seed=NA)}
}