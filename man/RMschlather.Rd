\name{RMschlather}
\alias{RMschlather}
\title{Covariance Model for binary field based on Gaussian field}
\description{
  \command{RMschlather} gives
  the tail correlation function of the extremal Gaussian
 process, i.e.

 \deqn{C(h) = 1 - \sqrt{ (1-\phi(h)/\phi(0)) / 2 }}

 where \eqn{\phi} is the covariance of a stationary Gaussian field.
 
}
\usage{
RMschlather(phi, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model yields the tail correlation function of the field
 that is returned by \command{\link{RPschlather}}
}
\value{
 \command{\link{RMschlather}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}}
\seealso{
 \command{\link{RPschlather}}
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

#Example 1
x <- if (interactive()) 0:1000 else 0:10
model <- RMgauss(scale=50)
z <- RFsimulate(RPschlather(model), x)
plot(z)

y  <- if (interactive()) 0:100 else 0:10
cov <- RFcov(model, y)
ecf <- RFcov(RMschlather(model), y)
matplot(y, cbind(cov, ecf), type="l")


# Example 2
x <- if (interactive()) 0:1000 else 0:10
model <- RMschlather(model)
z <- RFsimulate(RPschlather(model), x)
plot(z)

y  <- if (interactive()) 0:100 else 0:10
cov <- RFcov(model, y)
ecf <- RFcov(RMschlather(model), y)
matplot(y, cbind(cov, ecf), type="l")
\dontshow{FinalizeExample()}

}

\keyword{spatial}
\keyword{models}
