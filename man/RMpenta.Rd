\name{RMpenta}
\alias{RMpenta}
\title{Penta Covariance Model}
\description{
 \command{\link{RMpenta}} is a stationary isotropic covariance model, which is only valid for dimensions 
 \eqn{d \le 3}{d \le 3}.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = (1 - \frac{22}{3}r^{2} + 33 r^{4} - \frac{77}{2} r^{5} + \frac{33}{2} r^{7} - \frac{11}{2} r^{9} + \frac{5}{6}r^{11}) 1_{[0,1]}(r) . }{C(r)=(1 - 22/3 r^{2} + 33 r^{4} - 77/2 r^{5} + 33/2 r^{7} - 11/2 r^{9} + 5/6 r^{11}) 1_{[0,1]}(r).}
 
}
\usage{
RMpenta(var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
The model is only valid for dimensions \eqn{d \le 3}{d \le 3}.

It has a 4 times differentiable covariance function with compact
support (cf. Chiles, J.-P. and Delfiner, P. (1999), p. 84).

}
\value{
 \command{\link{RMpenta}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}
\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 } 
}

\me
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMpenta()
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}}