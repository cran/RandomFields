\name{RMgengneiting}
\alias{RMgengneiting}
\alias{RMwendland}
\title{Gneiting-Wendland Covariance Models}
\description{
 \command{\link{RMgengneiting}} is a stationary isotropic covariance model family whose elements 
 are specified by the two parameters \eqn{\kappa}{\kappa} and \eqn{\mu}{\mu} with \eqn{n}{n} a non-negative integer and 
 \eqn{\mu \ge \frac{d}{2}}{\mu \ge d/2} with \eqn{d}{d} denoting the dimension of the random field
 (the models can be used for any dimension). 
 A corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points. For the case \eqn{\kappa = 0} the Gneiting-Wendland model
 equals the Askey model \command{\link{RMaskey}},
 \deqn{C(r) = (1-r)^\beta 1_{[0,1]}(r),\qquad\beta = \mu +1/2 = \mu +
 2\kappa + 1/2.}{
 C(r) = (1-r)^\beta 1_{[0,1]}(r), \beta = \mu + 1/2 = \mu + 2\kappa + 1/2.}
 For \eqn{\kappa = 1}
 the Gneiting model is given by
 \deqn{C(r) = \left(1+\beta r \right)(1-r)^{\beta} 1_{[0,1]}(r),
 \qquad \beta = \mu +2\kappa+1/2.}{ C(r) = (1+\beta r)(1-r)^\beta
 1_{[0,1]}(r), \beta = \mu + 2\kappa + 1/2.}
 If \eqn{\kappa = 2}
 \deqn{C(r) = \left(1 + \beta r + \frac{\beta^{2} -
 1}{3}r^{2} \right)(1-r)^{\beta} 1_{[0,1]}(r), \qquad
 \beta = \mu+2\kappa+1/2.}{ C(r) = (1 + \beta r + (\beta^2 - 1) r^2 / 3) (1 - r)^\beta
 1_{[0,1]}(r), \beta = \mu + 2\kappa + 1/2.} 
 In the case \eqn{\kappa = 3}
 \deqn{ C(r) = \left( 1 + \beta r + \frac{(2\beta^{2}-3)}{5} r^{2}+
 \frac{(\beta^2 - 4)\beta}{15} r^{3} \right)(1-r)^\beta 1_{[0,1]}(r),
 \qquad \beta = \mu+2\kappa + 1/2.}{ C(r) = (1 + \beta r + (2 \beta^2 - 3
 )r^2 / 5 + (\beta^2 - 4) \beta
 r^3 / 15)(1-r)^\beta 1_{[0,1]}(r), \beta = \mu + 2\kappa + 1/2.}
 A special case of this model is \command{\link{RMgneiting}}.
 \eqn{\ell}
}
\usage{
RMgengneiting(kappa, mu, var, scale, Aniso, proj)
}
\arguments{
 \item{kappa}{\eqn{0,\ldots,3}{0,...,3}}; it chooses between the three
 different covariance models above
 \item{mu}{\code{mu} has to be greater than or equal to
 \eqn{\frac{d}{2}}{d / 2} where \eqn{d} is the dimension of the
 random field.} 
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This isotropic family of covariance functions is valid for any dimension of the random field. 

A special case of this family is \command{\link{RMgneiting}} (with \eqn{s = 1}{s = 1} there) for the choice \eqn{\kappa = 3, \mu = 3/2}{\kappa = 3, \mu = 3/2}.

}
\value{
 \command{\link{RMgengneiting}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Gneiting, T. (1999)
 Correlation functions for atmospherical data analysis.
 \emph{Q. J. Roy. Meteor. Soc} Part A \bold{125}, 2449-2464.

 \item Wendland, H. (2005) \emph{Scattered Data Approximation.}
 Cambridge Monogr. Appl. Comput. Math.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMaskey}},
 \command{\link{RMbigneiting}},
 \command{\link{RMgneiting}},
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
model <- RMgengneiting(kappa=1, mu=1.5)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))


## same models:
model2 <- RMgengneiting(kappa=3, mu=1.5, scale= 1 / 0.301187465825)
plot(RMgneiting(), model2=model2, type=c("p", "l"), pch=20)


\dontshow{FinalizeExample()}

}
