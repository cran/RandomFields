\name{RMspheric}
\alias{RMspheric}
\title{Spherical Covariance Model}
\description{
 \command{\link{RMhyperbolic}} is a stationary isotropic covariance model
 which is only valid up to dimension 3.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{C(r) = \left(1 - \frac{3}{2} r + \frac{1}{2} r^3\right)
 1_{[0,1]}(r)}{C(r) = (1 - 1.5 r + 0.5 r^3) 1_{[0,1]}(r).}
}
\usage{
RMspheric(var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This covariance model is valid only for dimensions less than or equal
 to 3.

 The covariance function has a finite range.
}
\value{
 \command{\link{RMspheric}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.
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
set.seed(0)
model <- RMspheric()
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
}