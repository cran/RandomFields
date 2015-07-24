\name{RMcauchy}
\alias{RMcauchy}
\title{Cauchy Family Covariance Model}
\description{
 \command{\link{RMcauchy}} is a stationary isotropic covariance model
 belonging to the Cauchy family.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = (1 + r^2)^(-\gamma)}{C(r) = (1 + r^2)^(-\gamma)}
 where \eqn{\gamma > 0}{\gamma > 0}.
 See also \command{\link{RMgencauchy}}.
}
\usage{
RMcauchy(gamma, var, scale, Aniso, proj)
}
\arguments{
 \item{gamma}{a numerical value; should be positive
 to provide a valid covariance function for a random field of any dimension.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 The paramater \eqn{\gamma}{\gamma} determines the asymptotic power law. The smaller \eqn{\gamma}{\gamma}, the longer the long-range dependence. The covariance function is very regular near the origin, because its Taylor expansion only contains even terms and reaches its sill slowly.
 
 Each covariance function of the Cauchy family is a normal scale mixture.

 The generalized Cauchy Family (see \command{\link{RMgencauchy}})
 includes this family for the choice \eqn{\alpha = 2}{\alpha = 2} and
 \eqn{\beta = 2 \gamma}{\beta = 2 \gamma}.
 The generalized Hyperbolic Family (see \command{\link{RMhyperbolic}})
 includes this family for the choice \eqn{\xi = 0}{\xi = 0} and
 \eqn{\gamma = -\nu/2}{\gamma = -\nu/2}; in this case scale=\eqn{\delta}{\delta}.
 
 
}
\value{
 \command{\link{RMcauchy}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Gneiting, T. and Schlather, M. (2004)
 Stochastic models which separate fractal dimension and Hurst effect.
 \emph{SIAM review} \bold{46}, 269--282.
 
 \item Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMcauchytbm}},
 \command{\link{RMgencauchy}},
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
model <- RMcauchy(gamma=1)
x <- seq(0, 10, 0.02)
plot(model, xlim=c(-3, 3))
plot(RFsimulate(model, x=x, n=4))
\dontshow{FinalizeExample()}
}