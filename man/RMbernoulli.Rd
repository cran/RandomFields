\name{RMbernoulli}
\alias{RMbernoulli}
\title{Covariance Model for binary field based on a Gaussian field}
\description{
  \command{RMbernoulli} gives
  the centered covariance function of a binary field.
}
\usage{
RMbernoulli(phi, threshold, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{threshold}{real valued threshold, see \command{\link{RPbernoulli}}.
 }
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model yields the covariance function of the field
 that is returned by \command{\link{RPbernoulli}}
}
\value{
 \command{\link{RMbernoulli}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 Ballani, Schlather
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de}}
\seealso{
 \command{\link{RPbernoulli}}
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0)

threshold <- 0
x <- seq(0, 5, if (interactive()) 0.02 else 1)
GaussModel <- RMgneiting()
z <- RFsimulate(RPbernoulli(GaussModel, threshold=threshold),
                x=x, n=2)
plot(z)

z <- RFsimulate(RPbernoulli(GaussModel, threshold=threshold),
                x=x, n=if (interactive()) 10000 else 2, spConform=FALSE)
estim.cov <- apply(z, 1, function(x) cov(x, z[1,]))
plot(x, estim.cov)
lines(x, RFcov(x=x, RMbernoulli(RMgauss(), threshold=threshold)))
\dontshow{RFoptions(seed=NA)}
}