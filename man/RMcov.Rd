\name{RMcov}
\alias{RMcov}
\alias{RMcov}
\alias{RMCOV_X}
\title{Non-stationary covariance model corresponding to a variogram model}
\description{
  This function generalizes the well-known non-stationary covariance
  function \eqn{2\min\{x,y\}} of the Brownian motion with variogram
  \eqn{\gamma(x,y) = |x-y|}, \eqn{x,y\ge 0}
  to arbitrary variogram models any spatial processes of any dimension
  and multivariability.

  Furthermore, the
  standard condition for the Brownian motion \eqn{W} is that 
  variance equals \eqn{0} at the origin,
  i.e., \eqn{W(x) =^d Z(x) -Z(0)} for any zero mean Gaussian process
  \eqn{Z} with variogram \eqn{\gamma(x,y) = |x-y|} is replaced by
  \eqn{W(x) = Z(x) -\sum_{i=1}^n a_i Z(x_i)} with \eqn{\sum_{i=1}^n a_i
    = 1}.

  For a given variogram \eqn{\gamma}, \eqn{a_i} and \eqn{x_i}, the model
  equals
  \eqn{C(x, y) = \sum_{i=1}^n a_i (\gamma(x, x_i) + \gamma(x_i, y)) -
    \gamma(x, y) - \sum_{i=1}^n \sum_{j=1}^n a_i a_j \gamma(x_i, y_i)
  }
}
\usage{
RMcov(gamma, x, y=NULL, z=NULL, T=NULL, grid, a,
       var, scale, Aniso, proj, raw, norm)
}
\arguments{
  \item{gamma}{a variogram model. Possibly multivariate.}
  \item{x,y,z,T,grid}{
    The usual arguments as in \command{\link{RFsimulate}} to define the
    locations where the covariates are given.
    Additional \code{x} might be set to one of the values
    \code{"origin"}, \code{"center"}, \code{"extremals"}, or \code{"all"}.
    If \code{x} is not given, \code{x} is set to \code{"origin"}.
  }
  \item{a}{vector of weights. The length of \code{a} must equal the
    number of points given by \code{x}, \code{y}, \code{z} and \code{T}.
    The values of \code{a} must sum up to \eqn{1}.
    If \code{a} is not given, equals weights are used.
  }
  \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
  \command{\link{RMmodel}}. If not passed, the above covariance function
  remains unmodified.} 
  \item{raw}{\argRaw}
  \item{norm}{\argNorm}
}
%\details{}
\value{
 \command{\link{RMcov}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
bm <- RMfbm(alpha=1)
plot(bm)

x <- seq(0, 6, if (interactive()) 0.125 else 3)
plot(RFsimulate(bm, x))

## standardizing with the random variable at the origin
z1 <- RFsimulate(RMcov(bm), x)
plot(z1)
z1 <- as.vector(z1)
zero <- which(abs(x) == 0)
stopifnot(abs(z1[zero]) < 1e-13)

## standardizing with the random variable at the center of the interval
z2 <- RFsimulate(RMcov(bm, "center"), x)
plot(z2)
z2 <- as.vector(z2)
stopifnot(abs(z2[(length(z2) + 1) / 2]) < 1e-13)


## standardizing with the random variables at the end points of the interval
z3 <- RFsimulate(RMcov(bm, "extremals"), x)
plot(z3)
z3 <- as.vector(z3)
stopifnot(abs(z3[1] + z3[length(z3)]) < 1e-13)

\dontshow{FinalizeExample()}
}