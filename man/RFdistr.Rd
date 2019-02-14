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
  Through \command{\link{RRdistr}} distribution families can be passed to \pkg{RandomFields} to create distributions available
  in the \command{\link{RMmodel}} definitions.
}
\usage{
RFddistr(model, x, params, dim=1, ...)
RFpdistr(model, q, params, dim=1, ...)
RFqdistr(model, p, params, dim=1, ...)
RFrdistr(model, n, params, dim=1, ...)
RFdistr(model, x, q, p, n, params, dim=1, ...)
}
\arguments{
  \item{model,params}{an \code{\link{RRmodel}}.}
  \item{x}{the location where the density is evaluated}
  \item{q}{the location where the probability function is evaluated}
  \item{p}{the value where the quantile function is evaluated}
  \item{n}{the number of random values to be drawn}
  \item{dim}{the dimension of the vector to be drawn}
  \item{\dots}{for advanced use:
    further options and control arguments for the simulation
    that are passed to and processed by \command{\link{RFoptions}}}
}
\details{
  \command{RFdistr} is the generic function for the 4 functions
  belonging to a distribution. 
}
\value{
  as described in the arguments
}
\me

\seealso{
  \command{\link{RRgauss}}, \link{RR}
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

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
model <- RMgauss(scale=exp())
for (i in 1:3) {
  RFoptions(seed = i + 10)
  readline(paste("Model no.", i, ": press return", sep=""))
  plot(model)
  readline(paste("Simulation no.", i, ": press return", sep=""))  
  plot(RFsimulate(model, x=seq(0,10,0.1)))
}
\dontshow{FinalizeExample()}}
\keyword{spatial}
