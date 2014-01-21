\name{RMnugget}
\alias{RMnugget}
\title{Nugget Effect Covariance Model}
\description{
  \command{\link{RMnugget}} is a multivariate stationary isotropic covariance model
 called \dQuote{nugget effect}.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given for \eqn{i,j}
 in \eqn{1,...,}\code{vdim} by
 \deqn{C_{ij}(r) = \delta_ij 1_{0}(r),}
 where \eqn{\delta_{ij}=1} if \eqn{i=j} and \eqn{\delta_{ij}=0} otherwise.
}
\usage{
RMnugget(tol, vdim, var, scale, Aniso, proj)
}
\arguments{
 \item{tol}{Only for advanced users. See \command{\link{RPnugget}}.}
 \item{vdim}{Only for advanced users. See \command{\link{RPnugget}}.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 Note that the parameter \code{scale} does not affect the covariance
 model; \code{Aniso} has an effect in case of zonal anisotropy.

 The nugget effect belongs to Gaussian white noise and is often used for
 modeling measurement errors.

 The 
 locations at a distance less than or equal to \code{nugget.tol}
 are considered as being identical. This strategy applies to
 the simulation method and the covariance function itself.
 Hence, the covariance function is only positive definite
 if \code{nugget.tol=0.0}. However, if the anisotropy matrix
 does not have full rank and \code{nugget.tol=0.0} then,
 the simulations are likely to be odd. 
 The value of \code{nugget.tol}
 should be of order \eqn{10^{-15}}{1e-15}.
}
\value{
 \command{\link{RMnugget}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
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




\examples{
RFoptions(seed=0)
model <- RMnugget(Aniso=matrix(1, nr=2, nc=2))
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(RFsimulate(model, x=x, y=x, tol=1e-10))
\dontshow{RFoptions(seed=NA)}
}