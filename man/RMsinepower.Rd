\name{RMsinepower}
\alias{RMsinepower}
\alias{sine power function}
\title{The Sinepower Covariance Model on the Sphere}
\description{
  \command{\link{RMsinepower}} is a isotropic covariance model. The
  corresponding covariance function, the sine power function of
  Soubeyrand, Enjalbert and Sache, only depends on the angle \eqn{\theta \in [0,\pi]}{0 \le \theta \le \pi} between two points on the sphere and is given by
 \deqn{\psi(\theta) = 1 - ( sin\frac{\theta}{2} )^{\alpha}}{\psi(\theta) = 1 - ( sin(\theta/2) )^{\alpha},}
 where \eqn{\alpha\in (0,2]}{0 < \alpha \le 2}. 
}

\usage{
RMsinepower(alpha, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value in \eqn{(0,2]}}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 For the sine power function of Soubeyrand, Enjalbert and Sache, see
 Gneiting, T. (2013) equation (17). For a more general form see \command{\link{RMchoquet}}.
}
\value{
 \command{\link{RMsinepower}} returns an object of class \command{\link[=RMmodel-class]{RMmodel}}
}
\references{
  Gneiting, T. (2013)
  \emph{Strictly and non-strictly positive definite functions on
 spheres} \emph{Bernoulli}, \bold{19}(4), 1327-1349.
}

\author{Christoph Berreth, \email{cberreth@mail.uni-mannheim.de},
  Martin Schlather
}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{spherical models}},
 \command{\link{RMchoquet}}
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}

RFoptions(coord_system="sphere")
model <- RMsinepower(alpha=1.7)
plot(model, dim=2)

## the following two pictures are the same
x <- seq(0, 0.4, 0.01)
z1 <- RFsimulate(model, x=x, y=x)
plot(z1)

x2 <- x * 180 / pi
z2 <- RFsimulate(model, x=x2, y=x2, coord_system="earth")
plot(z2)

stopifnot(all.equal(as.array(z1), as.array(z2)))

RFoptions(coord_system="auto")
\dontshow{FinalizeExample()}
}