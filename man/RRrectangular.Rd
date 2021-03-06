\name{RRrectangular}
\alias{RRrectangular}
\title{Random scaling used with balls}
\description{
 Approximates an isotropic decreasing density function
 by a density function that is isotropic with respect to the \eqn{l_1} norm.
 
 }
 \usage{
RRrectangular(phi, safety, minsteplen, maxsteps, parts, maxit,
             innermin, outermax, mcmc_n, normed, approx, onesided) 
}

\arguments{
  \item{phi}{a shape function; it is the user's responsibility that it
    is non-negative. See Details.}
  \item{safety, minsteplen, maxsteps, parts, maxit, innermin, outermax, mcmc_n}{
    Technical arguments to run an algorithm to simulate from this
    distribution. See \command{\link{RFoptions}} for the default values.
  }
  \item{normed}{logical. If \code{FALSE} then the norming constant
    \eqn{c} in the Details is set to \eqn{1}.
    This affects the values of the density function, the
    probability distribution and the quantile function, but not
    the simulation of random variables.
  }
  \item{approx}{logical.
    Default is \bold{\code{TRUE}}. If \code{TRUE}
    the isotropic distribution with respect to the \eqn{l_1} norm
    is returned. If \code{FALSE} then the exact isotropic distribution
    with respect to the \eqn{l_2} norm is simulated.
    Neither the density function, nor the probability distribution, nor
    the quantile function will be available if \code{approx=TRUE}.
  }
  \item{onesided}{logical.
    Only used for univariate distributions.
    If \code{TRUE} then the density is assumed to be non-negative only
    on the positive real axis. Otherwise the density is assumed to be
    symmetric.
    }
}

\details{
  This model defines an isotropic density function $f$ with respect to the
  \eqn{l_1} norm, i.e. \eqn{f(x) = c \phi(\|x\|_{l_1})}
  with some function \eqn{\phi}.
  Here, \eqn{c} is a norming constant so that the integral of \eqn{f}
  equals one.

  In case that \eqn{\phi} is monotonically decreasing then rejection sampling
  is used, else MCMC.

  The function \eqn{\phi} might have a polynomial pole at the origin
  and asymptotical decreasing of the form \eqn{x^\beta
    exp(-x^\delta)}.
}

\value{
 \command{\link{RRrectangular}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}


\me
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RRdistr}},
\command{\link{RRgauss}}.
 }

\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
# simulation of Gaussian variables (in a not very straightforward way):
distr <- RRrectangular(RMgauss(), approx=FALSE)
z <- RFrdistr(distr, n=1000000)
hist(z, 200, freq=!TRUE)
x <- seq(-10, 10, 0.1)
lines(x, dnorm(x, sd=sqrt(0.5)))


#creation of random variables whose density is proportional
# to the spherical model:
distr <- RRrectangular(RMspheric(), approx=FALSE)
z <- RFrdistr(distr, n=1000000)
hist(z, 200, freq=!TRUE)
\dontshow{StartExample(reduced=FALSE, save.seed=FALSE)}
x <- seq(-10, 10, 0.01)
lines(x, 4/3 * RFcov(RMspheric(), x))
\dontshow{FinalizeExample()}}
 



