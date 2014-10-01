\name{RMgauss}
\alias{RMgauss}
\title{Gaussian Covariance Model}
\description{
 \command{\link{RMgauss}} is a stationary isotropic covariance model.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{C(r) = e^{-r^2}}{C(r)=e^{-r^2}.}
}
\usage{
RMgauss(var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model is called Gaussian because of the functional similarity of
 the spectral density of a process with that covariance function to the
 Gaussian probability density function.
 
 The Gaussian model has an infinitely differentiable covariance
 function. This smoothness is artificial. Furthermore, this often leads to
 singular matrices and therefore numerically instable procedures
 (cf. Stein, M. L. (1999), p. 29).
% See \command{\link{RMgneiting}} for an alternative model that does not
% have the disadvantages of the Gaussian model.

 The Gaussian model is included in the symmetric stable class (see
 \command{\link{RMstable}}) for the choice \eqn{\alpha = 2}{alpha = 2}.
}
\value{
 \command{\link{RMgauss}} returns an object of class
 \code{\link[=RMmodel-class]{RMmodel}}
}

\note{
  The use of \command{RMgauss} is questionable from both a theoretical
  (analytical paths) and a practical point of view (e.g., speed of
  algorithms).
  Instead, \command{\link{RMgneiting}} should be used.
}

\references{
 Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.

 Stein, M. L. (1999) \emph{Interpolation of Spatial Data.} New York: Springer-Verlag 
 }

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMstable}} and \command{\link{RMmatern}} for generalisations;
 \cr
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.

 Do not mix up with \command{\link{RPgauss}} or \command{\link{RRgauss}}.
}


\keyword{spatial}
\keyword{models}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMgauss(scale=0.4)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
lines(RMgauss(), col="red")
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}
