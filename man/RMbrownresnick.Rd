\name{RMbrownresnick}
\alias{RMbrownresnick}
\title{Tail correlation function of the Brown-Resnick process}
\description{
  \command{RMbrownresnick} defines
  the tail correlation function of the Brown-Resnick process.
  \deqn{C(h) = 2 - 2\Phi(\sqrt{\gamma(h) / 2})}
  where \eqn{Phi} is the standard normal distribution function.
}
\usage{
RMbrownresnick(phi, var, scale, Aniso, proj);
}
\arguments{
 \item{phi}{variogram of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\details{
 For a given \command{\link{RMmodel}} the function
 \code{\link{RMbrownresnick}(\link{RMmodel}())} 'returns' the tail correlation
 function of a Brown-Resnick process with variogram \command{\link{RMmodel}}.
}
\references{
  Kabluchko, Z., Schlather, M. & de Haan, L (2009)
  Stationary max-stable random fields associated to negative definite functions
  \emph{Ann. Probab.} \bold{37}, 2042-2065.
}
\seealso{
 \command{\link{RFsimulate}},
 \command{\link{RMstrokorb}},
 \command{\link{RMmodel}}.
}
\examples{
set.seed(0)
  #plot covariance model of type RMbrownresnick
  RMmodel <- RMfbm(alpha=1.5, scale=0.2)
  plot(RMbrownresnick(RMmodel))

  #simulate and plot corresponding Gaussian random field
  x <- if (interactive()) seq(-5, 5, 0.25) else seq(-5, 5, 5)
  z <- RFsimulate(RMbrownresnick(RMmodel), x=x, y=x, grid=TRUE, print=6)
  plot(z)
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}
