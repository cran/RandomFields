\name{ExtremalGaussian}
\alias{RPschlather}
\alias{extremal Gaussian}
\alias{extremal Gaussian process}

\title{Extremal Gaussian process}
\description{ 
  \command{RPschlather} defines an extremal Gaussian process.
}

\usage{
RPschlather(phi, tcf, xi, mu, s)
}

\arguments{
 \item{phi}{an \command{\link{RMmodel}}, see Details.}
 \item{tcf}{an \command{\link{RMmodel}} specifying the
   extremal correlation function; either \code{phi} or \code{tcf} must
   be given.}
 \item{xi,mu,s}{the extreme value index, the location parameter and the
   scale parameter, respectively, of the generalized extreme value
   distribution. See Details.
 }
}

\details{
  \GEV


 The argument \code{phi} can be any random field for
 which the expectation of the positive part is known at the origin.

 It simulates an Extremal Gaussian process \eqn{Z} (also
 called \dQuote{Schlather model}), which is defined by
 \deqn{Z(x) = \max_{i=1}^\infty X_i \max(0, Y_i(x)),
 }{Z(x) = max_{i=1, 2, ...} X_i * max(0, Y_i(x)),}
 where the \eqn{X_i} are the points of a Poisson point process on the
 positive real half-axis with intensity \eqn{c x^{-2} dx}{c/x^2 dx},
 \eqn{Y_i \sim Y}{Y_i ~ Y}
 are iid stationary Gaussian processes with a covariance function
 given by \code{phi}, and \eqn{c} is chosen such
 that \eqn{Z} has standard Frechet margins. \code{phi} must
 represent a stationary covariance model.
 }
 
 \note{Advanced options
   are \code{maxpoints} and \code{max_gauss}, see
   \command{\link{RFoptions}}.}

\me

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RPgauss}},
 \command{\link{maxstable}},
 \command{\link{maxstableAdvanced}}.
}

\keyword{spatial}


\examples{\dontshow{StartExample()}
RFoptions(seed=0, xi=0)
## seed=0: *ANY* simulation will have the random seed 0; set
##         RFoptions(seed=NA) to make them all random again

## xi=0: any simulated max-stable random field has extreme value index 0
x <- seq(0, 2,0.01)

## standard use of RPschlather (i.e. a standardized Gaussian field)
model <- RMgauss()
z1 <- RFsimulate(RPschlather(model), x)
plot(z1, type="l")

## the following refers to the generalized use of RPschlather, where
## any random field can be used. Note that 'z1' and 'z2' have the same
## margins and the same .Random.seed (and the same simulation method),
## hence the same values
model <- RPgauss(RMgauss(var=2))
z2 <- RFsimulate(RPschlather(model), x)
plot(z2, type="l")
all.equal(z1, z2) # true

## Note that the following definition is incorrect
try(RFsimulate(model=RPschlather(RMgauss(var=2)), x=x))


## check whether the marginal distribution (Gumbel) is indeed correct:
model <- RMgauss()
z <- RFsimulate(RPschlather(model, xi=0), x, n=100)
plot(z)
hist(unlist(z@data), 50, freq=FALSE)
curve(exp(-x) * exp(-exp(-x)), from=-3, to=8, add=TRUE) 


\dontshow{FinalizeExample()}}

