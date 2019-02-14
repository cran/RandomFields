\name{RMsign}
\alias{RMsign}
\alias{RRsign}
\title{Random sign}
\description{
  \command{RMsign} defines a random sign. It can be used as part of the
  model definition of a Poisson field. 

 }
\usage{
RMsign(phi, p)
}

\arguments{
 \item{phi}{shape function of class \code{\link[=RMmodel-class]{RMmodel}}}
 \item{p}{probability of keeping the sign}
}

\details{
  \command{RMsign} changes the sign of the shape function \code{phi}
  with probability \code{1-p} and keeps it otherwise.
}
 
\value{
 \command{\link{RMsign}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}

\note{
  Random univariate or multivariate objects usually start with \code{RR}, not with \code{RM}.
  This is an exception here, as it operates on shape functions.
}

\me
\seealso{
  \command{\link{RMmodel}},
  \link{RR}.

 }

\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RPpoisson(RMsign(RMtent(), p=0.8))
x <- seq(0, 10, 0.02)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}}