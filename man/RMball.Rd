\name{RMball}
\alias{RMball}
\title{RMball}
\description{
  \command{RMball} refers to the indicator function of a ball
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
  \command{\link{RMpolygon}},
  \command{\link{RMspheric}},
 \command{\link{RFsimulate}},
 \command{\link{RMmodel}}.
}

\me

\keyword{spatial}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- seq(0,10,len=100)
model <- RMball(var=2,scale=1.5)
plot(model)
z <- RFsimulate(RPpoisson(model),x=x,y=x,intensity=0.1)
plot(z)

\dontshow{\dontrun{
model <- RPcoins(shape=RMball())
x <- seq(0, 10, 0.02)
plot(RFsimulate(model, x=x))
}}
\dontshow{FinalizeExample()}}