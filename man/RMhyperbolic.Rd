\name{RMhyperbolic}
\alias{RMhyperbolic}
\title{Generalized Hyperbolic Covariance Model}
\description{
 \command{\link{RMhyperbolic}} is a stationary isotropic covariance model
 called \dQuote{generalized hyperbolic}.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{C(r) = \frac{(\delta^2+r^2)^{\nu/2}
 K_\nu(\xi(\delta^2+r^2)^{1/2})}{\delta^\nu K_\nu(\xi
 \delta)}}{C(r) = \delta^(-\nu) (K_\nu(\nu \delta))^{-1}
 (\delta^2+r^2)^{\nu/2} K_\nu(\xi(\delta^2+r^2)^{1/2})}
 where \eqn{K_{\nu}}{K_\nu} denotes the modifies Bessel function of
 second kind. 
}
\usage{
RMhyperbolic(nu, lambda, delta, var, scale, Aniso, proj)
}
\arguments{
 \item{nu, lambda, delta}{numerical values; should either satisfy\cr
 \eqn{\delta \ge 0}{\delta \ge 0}, \eqn{\lambda > 0}{\lambda > 0}
 and \eqn{\nu > 0}{\nu > 0}, or\cr
 \eqn{\delta > 0}{\delta > 0}, \eqn{\lambda > 0}{\lambda > 0} and
 \eqn{\nu = 0}{\nu = 0}, or\cr
 \eqn{\delta > 0}{\delta > 0}, \eqn{\lambda \ge 0}{\lambda \ge 0}
 and \eqn{\nu < 0}{\nu < 0}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This class is over-parametrized, i.e. it can be reparametrized by
 replacing the three parameters \eqn{\lambda}{\lambda},
 \eqn{\delta}{\delta} and scale by two other parameters. This means
 that the representation is not unique.

 Each generalized hyperbolic covariance function is a normal scale
 mixture.

 The model contains some other classes as special cases;
 for \eqn{\lambda = 0}{\lambda = 0} we get Cauchy covariance function
 (see \command{\link{RMcauchy}}) with \eqn{\gamma =
 -\frac{\nu}2}{\gamma = -\nu/2} and scale=\eqn{\delta}{\delta};
 the choice \eqn{\delta = 0}{\delta = 0} yields a covariance model of type
 \command{\link{RMwhittle}} with smoothness parameter \eqn{\nu}{\nu}
 and scale parameter \eqn{\lambda^{-1}}{1/\lambda}. 
}
\value{
 \command{\link{RMhyperbolic}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Shkarofsky, I.P. (1968) Generalized turbulence space-correlation and
 wave-number spectrum-function pairs. \emph{Can. J. Phys.} \bold{46},
 2133-2153.
 
 \item Barndorff-Nielsen, O. (1978) Hyperbolic distributions and
 distributions on hyperbolae. \emph{Scand. J. Statist.} \bold{5}, 151-157.
 
 \item Gneiting, T. (1997). Normal scale mixtures and dual
 probability densities. \emph{J. Stat. Comput. Simul.} \bold{59}, 375-384.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMcauchy}},
 \command{\link{RMwhittle}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
model <- RMhyperbolic(nu=1, lambda=2, delta=0.2)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}