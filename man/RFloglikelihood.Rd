\name{RFloglikelihood}
\alias{RFloglikelihood}
\alias{RFlikelihood}
\title{Likelihood and estimation of linear models}
\description{
 \command{\link{RFloglikelihood}} returns the log likelihood for Gaussian
 random fields. In case NAs are given that refer to linear modeling, the
 ML of the linear model is returned.
}
\usage{
RFlikelihood(model, x, y = NULL, z = NULL, T = NULL, grid = NULL,
                data, distances, dim, likelihood,
                estimate_variance =NA, ...) 

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
 \item{likelihood}{ not programmed yet. Character.
   choice of kind of likehood ("full", "composite", etc.),
   see also \code{likelihood} for \command{\link{RFfit}}
   in \command{\link{RFoptions}}.
 }
 %\item{log}{logical. If \code{TRUE} the loglikelihood is returned.
% }
 \item{estimate_variance}{logical or \code{NA}. See Details.
 } 
 \item{...}{for advanced
   further options and control arguments for the simulation
   that are passed to and processed by \command{\link{RFoptions}}
 }

}
\details{
  The function calculates the likelihood for data of a Gassian process
  with given covariance structure.
  The covariance structure may not have \code{NA} values in the
  parameters except for a global variance. In this case the variance
  is returned that maximizes the likelihood.
  Additional to the covariance structure the model may include a
  trend. The latter may contain unknown linear parameters.
  In this case again, the unknown parameters are estimated, and returned.
}
\value{
  \command{\link{RFloglikelihood}} returns a list
  containing the likelihood, the log likelihood, and
  the global variance (if estimated -- see details).
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{	
\link{Bayesian},
 \command{\link{RMmodel}},
 \command{\link{RFfit}},
 \command{\link{RFsimulate}},
 \command{\link{RFlinearpart}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


require("mvtnorm")

pts <- 5
repet <- 3
model <- RMexp()
x <- runif(n=pts, min=-1, max=1)
y <- runif(n=pts, min=-1, max=1)
data <- as.matrix(RFsimulate(model, x=x, y=y, n=repet, spC = FALSE))
print(cbind(x, y, data))
print(unix.time(likeli <- RFlikelihood(model, x, y, data=data)))
str(likeli, digits=8)

L <- 0
C <- RFcovmatrix(model, x, y)
for (i in 1:ncol(data)) {
  print(unix.time(dn <- dmvnorm(data[,i], mean=rep(0, nrow(data)),
sigma=C, log=TRUE)))
  L <- L + dn
}
print(L)
stopifnot(all.equal(likeli$log, L))



%--------------------------------------------------------------


pts <- 5
repet <- 1
trend <- 2 * sin(R.p(new="isotropic")) + 3
#trend <- RMtrend(mean=0)
model <- 2 * RMexp() + trend
x <- seq(0, pi, len=10)
data <- as.matrix(RFsimulate(model, x=x, n=repet, spC = FALSE))
print(cbind(x, y, data))

print(unix.time(likeli <- RFlikelihood(model, x, data=data)))
str(likeli, digits=8)

L <- 0
tr <- RFfctn(trend, x=x, spC = FALSE)
C <- RFcovmatrix(model, x)
for (i in 1:ncol(data)) {
  print(unix.time(dn <- dmvnorm(data[,i], mean=tr, sigma=C, log=TRUE)))
  L <- L + dn
}
print(L)

stopifnot(all.equal(likeli$log, L))



%--------------------------------------------------------------



pts <- c(4, 5)
repet <- c(2, 3)
trend <- 2 * sin(R.p(new="isotropic")) + 3
model <- 2 * RMexp() + trend
x <- y <- data <- list()
for (i in 1:length(pts)) {
  x[[i]] <- list(x = runif(n=pts[i], min=-1, max=1),
                 y = runif(n=pts[i], min=-1, max=1))
  data[[i]] <- as.matrix(RFsimulate(model, x=x[[i]]$x, y=x[[i]]$y,
                                    n=repet[i], spC = FALSE))
}

print(unix.time(likeli <- RFlikelihood(model, x, data=data)))
str(likeli, digits=8)

L <- 0
for (p in 1:length(pts)) {
  tr <- RFfctn(trend, x=x[[p]]$x, y=x[[p]]$y,spC = FALSE)
  C <- RFcovmatrix(model, x=x[[p]]$x, y=x[[p]]$y)
  for (i in 1:ncol(data[[p]])) {
    print(unix.time(dn <- dmvnorm(data[[p]][,i], mean=tr, sigma=C,
log=TRUE)))
    L <- L + dn
  }
}
print(L)


stopifnot(all.equal(likeli$log, L))



\dontshow{FinalizeExample()}

}
\keyword{spatial}






