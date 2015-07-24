\name{RMmultiquad}
\alias{RMmultiquad}
\alias{Poisson spline}
\alias{Inverse multiquadric}
\alias{multiquadric family}
\title{The Multiquadric Family Covariance Model on th Sphere}
\description{
  \command{\link{RMmultiquad}} is a isotropic covariance model. The
  corresponding covariance function, the multiquadric family, only
  depends on the angle \eqn{\theta \in [0,\pi]}{0 \le \theta \le \pi}
  between two points on the sphere and is given by 
 \deqn{\psi(\theta) = (1 - \delta)^{2*\tau} / (1 + delta^2 - 2*\delta*cos(\theta))^{\tau}}{\psi(\theta) = (1 - \delta)^{2*\tau} / (1 + delta^2 - 2*\delta*cos(\theta))^{\tau},}
 where \eqn{\delta \in (0,1)}{0 < \delta < 1} and \eqn{\tau > 0}{\tau > 0}. 
}

\usage{
RMmultiquad(delta, tau, var, scale, Aniso, proj)
}
\arguments{
  \item{delta}{a numerical value in \eqn{(0,1)}}
  \item{tau}{a numerical value greater than \eqn{0}}
  \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
    \command{\link{RMmodel}}. If not passed, the above
    covariance function remains unmodified.}
}
\details{
 Special cases (cf. Gneiting, T. (2013), p.1333) are known for fixed parameter
 \eqn{\tau=0.5}{\tau=0.5} which leads to the covariance function called
 'inverse multiquadric'\deqn{\psi(\theta) = (1 - \delta) / \sqrt( 1 + delta^2 - 2*\delta*cos(\theta) )}{\psi(\theta) = (1 - \delta) / \sqrt( 1 + delta^2 - 2*\delta*cos(\theta) )} and for fixed parameter \eqn{\tau=1.5}{\tau=1.5}
 which gives the covariance function called 'Poisson spline' 
 \deqn{\psi(\theta) = (1 - \delta)^{3} / (1 + delta^2 -
   2*\delta*cos(\theta))^{1.5}}{\psi(\theta) = (1 - \delta)^{3} / (1 +
   delta^2 - 2*\delta*cos(\theta))^{1.5}.}
 For a more general form, see \command{\link{RMchoquet}}.
}
\value{
 \command{\link{RMmultiquad}} returns an object of class \command{\link[=RMmodel-class]{RMmodel}}
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
 \command{\link{RMchoquet}},
 \command{\link{spherical models}}

}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
RFoptions(coord_system="sphere")
model <- RMmultiquad(delta=0.5, tau=1)
plot(model, dim=2)

## the following two pictures are the same
%nicht mehr als 0.12 da inversen der cov matrizen sich dann stark unterscheiden
x <- seq(0, 0.12, 0.01)
z1 <- RFsimulate(model, x=x, y=x)
plot(z1)

x2 <- x * 180 / pi
z2 <- RFsimulate(model, x=x2, y=x2, coord_system="earth")
plot(z2)

stopifnot(all.equal(as.array(z1), as.array(z2)))

RFoptions(coord_system="auto")
\dontshow{FinalizeExample()}

\dontshow{FinalizeExample()}
}