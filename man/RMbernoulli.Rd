\name{RMbernoulli}
\alias{RMbernoulli}
\title{Covariance Model for binary field based on a Gaussian field}
\description{
  \command{RMbernoulli} gives
  the centered \bold{correlation} function of a binary field,
  obtained by thresholding a Gaussian field.
}
\usage{
RMbernoulli(phi, threshold, correlation, centred, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{threshold}{real valued threshold, see
   \command{\link{RPbernoulli}}.
   Currently, only \command{threshold=0.0} is possible. %to do

   Default: 0.
 }
 \item{correlation}{logical. If \code{FALSE} the corresponding
   covariance function is returned.

   Default: \code{TRUE}.
 }
 \item{centred}{logical. If \code{FALSE} the uncentred covariance is
   returned.

   Default: \code{TRUE}.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model yields the covariance function of the field
 that is returned by \command{\link{RPbernoulli}}.
}
\note{
\bold{Previous to version 3.0.33 the covariance function was returned,
  not the correlation function.}
}
\value{
 \command{\link{RMbernoulli}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 Ballani, Schlather
}

\me
\seealso{
 \command{\link{RPbernoulli}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}}.
}


\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

threshold <- 0
x <- seq(0, 5, 0.02)
GaussModel <- RMgneiting()

n <- 1000
z <- RFsimulate(RPbernoulli(GaussModel, threshold=threshold), x=x, n=n)
plot(z)

model <- RMbernoulli(RMgauss(), threshold=threshold, correlation=FALSE)
plot(model, xlim=c(0,5))
z1 <- as.matrix(z)
estim.cov <- apply(z1, 1, function(x) cov(x, z1[1,]))
points(coordinates(z), estim.cov, col="red")

\dontshow{FinalizeExample()}}