\name{RMnugget}
\alias{RMnugget}
\alias{RM_NUGGET}
\title{Nugget Effect Covariance Model}
\description{
  \command{\link{RMnugget}} is a multivariate stationary isotropic covariance model
 called \dQuote{nugget effect}.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given for \eqn{i,j}
 in \eqn{1,...,}\code{vdim} by
 \deqn{C_{ij}(r) = \delta_{ij} 1_{0}(r),}
 where \eqn{\delta_{ij}=1} if \eqn{i=j} and \eqn{\delta_{ij}=0} otherwise.
}
\usage{
RMnugget(tol, vdim, var, Aniso, proj)
}
\arguments{
 \item{tol}{Only for advanced users. See \command{\link{RPnugget}}.}
 \item{vdim}{Must be set only for multivariate models (advanced).} 
 \item{var}{optional argument; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
 \item{Aniso,proj}{(zonal modelling and repeated
   measurements(advanced)); see \command{\link{RPnugget}} for details.
 }
}
\details{
   The nugget effect belongs to Gaussian white noise and is used for
   modeling measurement errors or to model spatial \sQuote{nuggets}.
 }
 \value{
 \command{\link{RMnugget}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}


\me
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RPnugget}} (advanced users).
}

\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- y <- 1:2
xy <- as.matrix(expand.grid(x, y)) ## we get 4 locations

## Standard use of the nugget effect
model <- RMnugget(var = 100)
RFcovmatrix(model, x=xy)
as.vector(RFsimulate(model, x=x, y=x, tol=1e-10))

\dontshow{FinalizeExample()}}
