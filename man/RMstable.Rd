\name{RMstable}
\alias{RMstable}
\alias{RMpoweredexp}
\alias{RMpoweredexponential}
\title{Stable Family / Powered Exponential Model}
\description{
 \command{\link{RMstable}} is a stationary isotropic covariance model
 belonging to the so called stable family.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{C(r) = e^{-r^\alpha}}{C(r)=e^{-r^\alpha}}
 where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2}.
}
\usage{
RMstable(alpha, var, scale, Aniso, proj)
RMpoweredexp(alpha, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value; should be in the interval (0,2]
 to provide a valid covariance function for a random field of any
 dimension. 
}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 The parameter \eqn{\alpha}{\alpha} determines the fractal dimension
 \eqn{D}{D} of the Gaussian sample paths:
 \deqn{ D = d + 1 - \frac{\alpha}{2}}{D = d + 1 - \alpha/2}
 where \eqn{d}{d} is the dimension of the random field.
 For \eqn{\alpha < 2}{\alpha < 2} the Gaussian sample paths are not
 differentiable (cf. Gelfand et al., 2010, p. 25).

 Each covariance function of the stable family is a normal scale mixture.
 
 The stable family includes the exponential model (see
 \command{\link{RMexp}}) \eqn{\alpha = 1}{\alpha = 1} and the Gaussian
 model (see \command{\link{RMgauss}}) for \eqn{\alpha = 2}{\alpha = 2}.

 The model is called stable, because in the 1-dimensional case the
 covariance is the characteristic function of a stable random variable
 (cf. Chiles, J.-P. and Delfiner, P. (1999), p. 90).
}
\value{
 \command{\link{RMstable}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 
 \item Diggle, P. J., Tawn, J. A. and Moyeed, R. A. (1998) Model-based
 geostatistics (with discussion). \emph{Applied Statistics} \bold{47},
 299--350.
 
 \item Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp, P. (eds.)
 (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.

 \item Strokorb, K., Ballani, F., and Schlather, M. (2014)
 In Preparation.
 }

}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMexp}},
 \command{\link{RMgauss}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMstable(alpha=1.9, scale=0.4)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}