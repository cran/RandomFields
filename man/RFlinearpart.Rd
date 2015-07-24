\name{RFlinearpart}
\alias{RFlinearpart}
\title{Linear part of \command{\link{RMmodel}}}
\description{
 \command{\link{RFlinearpart}} returns the linear part of a model
}
\usage{
RFlinearpart(model, x, y = NULL, z = NULL, T = NULL, grid,
                data, distances, dim, set=0, ...)

}
\arguments{
 \item{model}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 the covariance or variogram model, which is to be evaluated}
 \item{x}{vector or \eqn{(n \times \code{dim})}{(n x
 \code{dim})}-matrix, where \eqn{n} is the number of points at
 which the covariance function is to be evaluated;
 in particular,
 if the model is isotropic or \code{dim=1} then \code{x}
 is a vector. \code{x}}
 \item{y}{second vector or matrix for non-stationary covariance
 functions}
 \item{z}{z-component of point if xyzT-specification of points is used}
 \item{T}{T-component of point if xyzT-specification of points is used}
 \item{grid}{boolean; whether xyzT specify a grid}
 \item{data}{vector or matrix of values measured at \code{coord};
 If a matrix is given then the columns are interpreted as independent
 realisations.\cr
 If also a time component is given, then in the data the indices for
 the spatial components run the fastest.
 
 If an \code{m}-variate model is used, then each realisation is given as
 \code{m} consecutive columns of \code{data}.
 }
 \item{distances}{vector;
   the lower triangular part of the distance matrix column-wise;
 equivalently the upper triangular part of the distance matrix row-wise;
 either \code{x} or \code{distances} must be missing}
 \item{dim}{dimension of the coordinate space in which the model is
   applied; only necesary for given \code{distances}}
 \item{set}{integer. See section Value for details.}
 \item{...}{for advanced
   further options and control arguments for the simulation
   that are passed to and processed by \command{\link{RFoptions}}
 }

}


\note{
  In the linear part of the model specification the parameters
  that are NA must be the first model part. I.e.
  \code{NA * sin(R.p(new="isotropic")) + NA + R.p(new="isotropic")}
  is OK, but not
  \code{sin(R.p(new="isotropic")) * NA + NA + R.p(new="isotropic")}

}

\value{  
  \command{\link{RFlinearpart}} returns a list
  of three components, \code{Y}, \code{X}, \code{vdim} returning
  the deterministic trend, the design matrix, and the multivariability,
  respectively.
  If \code{set} is positive, \code{Y} and \code{X} contain
  the values for the \code{set}-th set of coordinates.
  Else, \code{Y} and \code{X} are both lists containing
  the values for all the sets.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{	
\link{Bayesian},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFlikelihood}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


x <- seq(0, pi, len=10)
trend <- 2 * sin(R.p(new="isotropic")) + 3
model <- RMexp(var=2, scale=1) + trend
print(RFlinearpart(model, x=x))  ## only a deterministic part

trend <- NA * sin(R.p(new="isotropic")) + NA + R.p(new="isotropic") / pi
model <- RMexp(var=NA, scale=NA) + trend
print(RFlinearpart(model, x=x))


\dontshow{FinalizeExample()}

}
\keyword{spatial}






