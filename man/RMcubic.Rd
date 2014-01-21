\name{RMcubic}
\alias{RMcubic}
\title{Cubic Covariance Model}
\description{
 \command{\link{RMcubic}} is a stationary isotropic covariance model
 which is only valid for dimensions \eqn{d \le 3 }{d \le 3 }.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = (1 -7 r^(2) + 8.75 r^3 - 3.5 r^5 + 0.75 r^7) 1_{[0,1]}(r).}{C(r) = (1 -7 r^(2) + 8.75 r^3 - 3.5 r^5 + 0.75 r^7) 1_{[0,1]}(r).}
 
}
\usage{
RMcubic(var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{The model is only valid for dimensions \eqn{d \le 3 }{d \le 3 }.
 It is a 2 times differentiable covariance function with compact support
 (cf. Chiles, J.-P. and Delfiner, P. (1999), p. 84).
}
\value{
 \command{\link{RMcubic}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

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
RFoptions(seed=0)
model <- RMcubic()
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{RFoptions(seed=NA)}
}