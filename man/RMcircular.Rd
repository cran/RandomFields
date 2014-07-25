\name{RMcircular}
\alias{RMcircular}
\title{Circular Covariance Model}
\description{
 \command{\link{RMcircular}} is a stationary isotropic covariance model
 which is only valid for dimensions \eqn{d \le 2 }{d \le 2 }.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = 1 - 2/\pi (r \sqrt(1-r^2) + arcsin(r)) 1_{[0,1]}.}{C(r) = 1 - 2/\pi (r \sqrt(1-h^2) + asin(r)) 1_{[0,1]}.}
 
}
\usage{
RMcircular(var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{The model is only valid for dimensions \eqn{d \le 2 }{d \le 2}.
 It is a covariance function with compact support
 (cf. Chiles, J.-P. and Delfiner, P. (1999), p. 82).
}
\value{
 \command{\link{RMcircular}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 }
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
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMcircular()
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}