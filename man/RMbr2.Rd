\name{RMbr2}
\alias{RMbr2bg}
\alias{RMbr2eg}
\title{Transformation among tail correlation functions}
\description{
 These functions can be used to model different max-stable processes
 with the same tail correlation functions.

 \deqn{
   C(h) = \cos(\pi (1 - 2\Phi(\sqrt{\gamma(h) / 8}) ) )
 }
 and
  \deqn{
    C(h) =1 - 2 (1 - 2 \Phi(\sqrt{\gamma(h) / 8}) )^2 
    ,
 }
 Here, \eqn{\Phi} is the standard normal distribution
 function, and \eqn{\gamma} is a variogram with sill
 \eqn{ 8(erf^{-1}(1/sqrt 2))^2 = 1.8197} and
 \eqn{8(erf^{-1}(1/2))^2 = 4.425098}, respectively.

}
\usage{
RMbr2bg(phi, var, scale, Aniso, proj)
RMbr2eg(phi, var, scale, Aniso, proj)

}
\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\details{
 \itemize{
 \item \command{\link{RMbr2bg}} \cr
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
 
 \item \command{\link{RMbr2eg}} \cr
 The extremal Gaussian model \command{\link{RPschlather}}
 simulated with \code{\link{RMbr2eg}(\link{RMmodel}())} has
 tail correlation function that equals
 the tail correlation function of Brown-Resnick process with
 variogram \command{\link{RMmodel}}.
 } 
}
\references{
 Strokorb, K. (2012) PhD thesis.
}
\seealso{
 \command{\link{RPbernoulli}},
 \command{\link{RPbrownresnick}},
 \command{\link{RPschlather}},
 \command{\link{RMmodel}}.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}
\examples{
set.seed(0)
## see 'maxstableAdvanced' for some examples
}