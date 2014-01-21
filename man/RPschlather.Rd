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
  The parameter \code{xi} is always a number, i.e. \eqn{\xi} is constant in
 space. In contrast, \eqn{\mu} and \eqn{s} might be constant
 numerical value or given a \code{\link{RMmodel}}, in particular by a
 \code{\link{RMtrend}} model. 


 The parameter \code{phi} can be any random field for
 which the expectation of the positive part is known at the origin.

 It simulates Extremal Gaussian process \eqn{Z} (also
 called \dQuote{Schlather model}), which is defined by
 \deqn{Z(x) = \max_{i=1}^\infty X_i \max(0, Y_i(x)),
 }{Z(x) = max_{i=1, 2, ...} X_i * max(0, Y_i(x)),}
 where the \eqn{X_i} are the points of a Poisson point process on the
 positive real half-axis with intensity \eqn{c x^{-2} dx}{c/x^2 dx},
 \eqn{Y_i \sim Y}{Y_i ~ Y}
 are iid stationary Gaussian processes with a covariance function
 given by \code{model}, and \eqn{c} is chosen such
 that \eqn{Z} has standard Frechet margins. \code{model} must
 represent a stationary covariance model.
 }
 
 \note{Advanced options
   are \code{maxpoints} and \code{max_gauss}, see
   \command{\link{RFoptions}}.}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}


\examples{
RFoptions(seed=0)
model <- RMgauss()
x <- seq(0, 10, if (interactive()) 0.1 else 1)
z <- RFsimulate(RPschlather(model, xi=0, s=0.1), x)
plot(z)

\dontshow{if (.C("isAuthor", a=integer(1))$a) {

model <- RMgauss()
x <- seq(0,10, 0.02)
z <- RFsimulate(RPschlather(model, xi=0), x,
                n=if (interactive()) 100 else 1)
plot(z)
hist(unlist(z@data), 50, freq=FALSE)
curve(exp(-x) * exp(-exp(-x)), from=-3, to=8, add=TRUE)

## for some more sophisticated models see 'maxstableAdvanced'

}}
\dontshow{RFoptions(seed=NA)}

}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RPgauss}},
 \command{\link{maxstable}},
 \command{\link{maxstableAdvanced}}
}

\keyword{spatial}


