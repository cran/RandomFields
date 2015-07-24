\name{RMdampedcos}
\alias{RMdampedcos}
\title{Exponentially Damped Cosine}
\description{
 \command{\link{RMdampedcos}} is a stationary isotropic covariance model.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = exp(-\lambda r) \cos(r).}{C(r) = exp(-\lambda r) cos(r).}
 }
\usage{
RMdampedcos(lambda, var, scale, Aniso, proj)
}
\arguments{
 \item{lambda}{numeric. The range depends on the dimension of the random
   field (see details)} 
 \item{var, scale, Aniso, proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{The model is valid for any dimension \eqn{d}{d}. However, depending on the dimension of
  the random field the following bound for the arguments
  \eqn{\lambda}{\lambda} has to be respected: 

\deqn{\lambda \ge 1/{\tan(\pi/(2d))}.}{\lambda \ge 1/{tan(\pi/(2d))}.}

This covariance models a hole effect (cf. Chiles, J.-P. and Delfiner,
P. (1999), p. 92).

For \eqn{\lambda = 0} we obtain the covariance function
\deqn{C(r)=\cos(r)}{C(r)=cos(r)} which is only valid for \eqn{d=1}{d=1}
and corresponds to \command{\link{RMbessel}} for
\eqn{\nu=-0.5}{\nu=-0.5}, there.

}
\value{
 \command{\link{RMdampedcos}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.

 \item Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010)
 \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMbessel}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
model <- RMdampedcos(lambda=0.3, scale=0.1)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}