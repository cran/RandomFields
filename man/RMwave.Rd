\name{RMwave}
\alias{RMwave}
\alias{RMcardinalsine}
\title{Wave Covariance Model / Cardinal Sine}
\description{
 \command{\link{RMwave}} is a stationary isotropic covariance model, which is valid only for dimension 
 \eqn{d \le 3}{d \le 3}.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = sin(r)/r 1_{r>0} + 1_{r=0} .}{C(r)=sin(r)/r 1_{r>0} + 1_{r=0} .}
 It is a special case of \command{\link{RMbessel}}.
}
\usage{
RMwave(var, scale, Aniso, proj)
RMcardinalsine(var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
The model is valid only for dimension \eqn{d \le 3}{d \le 3}. It is a
special case of \command{\link{RMbessel}} for \eqn{\nu = 0.5}{\nu = 0.5}.

This covariance models a hole effect (cf. Chiles, J.-P. and Delfiner, P. (1999), p. 92).
}
\value{
 \command{\link{RMwave}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{ \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMbessel}}
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMwave(scale=0.1)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}