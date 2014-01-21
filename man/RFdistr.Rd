\name{RFdistr}
\alias{RFdistr}
\alias{RFddistr}
\alias{RFpdistr}
\alias{RFqdistr}
\alias{RFrdistr}
\title{
  Evaluating distribution families
}
\description{
  Through \command{\link{RRdistr}} distribution families can the
  passed to \pkg{RandomFields} to create random parameters
}
\usage{
RFddistr(model, x, dim=1, ...)
RFpdistr(model, q, dim=1, ...)
RFqdistr(model, p, dim=1, ...)
RFrdistr(model, n, dim=1, ...)
RFdistr(model, x, q, p, n, dim=1, ...)
}
\arguments{
  \item{model}{
    an \code{RRmodel}}
  \item{x}{the location where the density is evaluated}
  \item{q}{the location there the probability function is evaluated}
  \item{p}{the value where the quantile function is evaluated}
  \item{n}{the number of random values to be drawn}
  \item{dim}{the dimension of the vector to be drawn}
  \item{\dots}{for advanced use:
    further options and control parameters for the simulation
    that are passed to and processed by \command{\link{RFoptions}}}
}
\details{
  \command{RFdistr} is the generic function for the 4 functions
  belonging to a distribution. 
}
\value{
  as described in the arguments
}
\author{
Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de
}
}

\seealso{
  \command{\link{RRgauss}}, \link{RR}
}
\examples{
RFoptions(seed=0)

Print(RFoptions())
## a very toy example to understand the use
model <- RRdistr(norm())
v <- 0.5
Print(RFdistr(model=model, x=v), dnorm(x=v)) 
Print(RFdistr(model=model, q=v), pnorm(q=v))
Print(RFdistr(model=model, p=v), qnorm(p=v))

n <- 10
r <- RFdistr(model=model, n=n, seed=0)
set.seed(0); Print(r, rnorm(n=n))


## note that a conditional covariance function given the
## random parameters is given here:
\dontshow{RFoptions(seed=NA)}
model <- RMgauss(scale=exp())
for (i in 1:3) {
  readline(paste("Simulation no.", i, ": press return", sep=""))
  plot(model)
  plot(RFsimulate(model, x=seq(0,10,0.1)))
}
\dontshow{RFoptions(seed=NA)}
}
\keyword{spatial}
