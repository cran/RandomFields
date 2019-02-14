\name{RMderiv}
\alias{RMderiv}
\alias{derivative}
\alias{gradient}
\title{Gradient of a field}
\description{
 \command{\link{RMderiv}}
 is a multivariate covariance model which
 models a field and its gradient.

 For an isotropic covariance model \eqn{varphi}, the covariance \eqn{C} given by
 \command{RMderiv} equals
 \deqn{C_{11}(x,y) = \varphi(\| x - y\|)}
 \deqn{C_{j1}(x,y) = -C_{1j}(x,y) = \partial \varphi(\|x - y\|) /
   \partial x}
 \deqn{C_{i,j}(x,y) = \partial^2 \varphi(\|x - y\|) /
   \partial x \partial y}
 for \eqn{i,j = 2,\ldots, d} where \eqn{d} is the dimension of the field.
 }
\usage{
RMderiv(phi, which, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a univariate stationary covariance model (in 2 or 3 dimensions).}
 \item{which}{vector of integers. If not given all components are
   returned; otherwise the selected components are returned.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}

\value{
 \command{\link{RMderiv}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Matheron
 }
}

\me
\seealso{
 \command{\link{RMcurlfree}},
 \command{\link{RMdivfree}},
 \command{\link{RMvector}}
}

\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMderiv(RMgauss(), scale=4)
plot(model, dim=2)

x.seq <- y.seq <- seq(-10, 10, 0.4)
simulated <- RFsimulate(model=model, x=x.seq, y=y.seq)

plot(simulated)
\dontshow{FinalizeExample()}}
