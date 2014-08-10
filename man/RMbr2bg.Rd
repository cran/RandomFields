\name{RMbr2bg}
\alias{RMbr2bg}
\title{Transformation from Brown-Resnick to Bernoulli}
\description{
  This function can be used to model a max-stable process
  based on the a binary field,  with the same extremal correlation
  function as a Brown-Resnick process
  \deqn{
    C_{bg}(h) = \cos(\pi (2\Phi(\sqrt{\gamma(h) / 2}) -1) )
  }
  Here, \eqn{\Phi} is the standard normal distribution
  function, and \eqn{\gamma} is a \bold{semi-}variogram with sill
  \deqn{4(erf^{-1}(1/2))^2 = 2 * { \Phi^{-1}( 3 / 4  ) }^2 =
    1.819746 / 2 = 0.9098728}
}
\usage{
RMbr2bg(phi, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\details{
  \command{\link{RMbr2bg}} \cr
  binary random field \command{\link{RPbernoulli}}
  simulated with \code{\link{RMbr2bg}(\link{RMmodel}())} has
  a uncentered covariance function that equals
  \enumerate{
    \item 
    the tail correlation function of
    the max-stable process constructed with this binary random field
    \item
    the tail correlation function of Brown-Resnick process with
    variogram \command{\link{RMmodel}}.
  }
 
  Note that the reference paper is based on the notion of the
  (genuine) variogram, whereas the package \pkg{RandomFields}
  is based on the notion of semi-variogram. So formulae
  differ by factor 2.
}
\references{
 Strokorb, K. (2012) PhD thesis. Goettingen.
}
\seealso{
  \link{maxstableAdvanced},
  \command{\link{RMbr2eg}},
  \command{\link{RMmodel}},
  \command{\link{RMm2r}},
  \command{\link{RPbernoulli}},
  \command{\link{RPbrownresnick}},
  \command{\link{RPschlather}},
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\keyword{spatial}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMexp(var=1.62 / 2) 
step <- if (interactive()) 0.05 else 2
y <- seq(0, 10, step)
z <- RFsimulate(RPschlather(RMbr2eg(model)), y, y)
plot(z)

\dontshow{FinalizeExample()}
}