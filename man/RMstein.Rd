\name{RMstein}
\alias{RMstein}
\title{Stein nonseparable space-time model}
\description{
 \command{\link{RMstein}} is a univariate stationary covariance model
 whose corresponding covariance function only depends on the difference
 \eqn{h}{h} between two points and is given by 
 \deqn{C(h, t) = W_{\nu}(y) - ( < h, z > t)/((\nu - 1)(2\nu + d))
 * W_{\nu-1}(y)}{ C(h, t) = W_\nu(y) - < h, z > t W_{\nu-1}(y) / [
 (\nu-1) (2\nu + d + 1) ]}

 Here \eqn{W_\nu} is the covariance of the
 \command{\link{RMwhittle}} model with
 smoothness parameter \eqn{\nu};
 \eqn{y=\|(h,t)\|}{y = ||(h,t)||} is the norm of the vector
 \eqn{(h,t)}{(h,t)},
 \eqn{d}{d} is the dimension of the space on which the random field is considered.
}

\usage{
RMstein(nu, z, var, scale, Aniso, proj)
}
\arguments{
 \item{nu}{numerical value; greater than 1; smoothness parameter of the
 RMwhittle model}
 \item{z}{a vector; the norm of \eqn{z}{z} must be less or equal to 1.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{See Stein (2005).
}
\value{
 \command{\link{RMstein}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Stein, M.L. (2005) Space-time covariance functions.
 \emph{J. Amer. Statist. Assoc.} \bold{100}, 310-321. Equation (8).
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMstein(nu=1.5, z=0.9)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(RFsimulate(model, x=x, y=x))
\dontshow{FinalizeExample()}
}