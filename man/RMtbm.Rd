\name{RMtbm}
\alias{RMtbm}
\title{Turning Bands Method}
\description{
 \command{\link{RMtbm}} is a univariate stationary isotropic covariance
 model in dimension \code{reduceddim} which depends on a univariate stationary
 isotropic covariance \eqn{\phi}{phi} in a bigger dimension \code{fulldim}.
 For formulas for the covariance function see details.
}
\usage{
RMtbm(phi, fulldim, reduceddim, layers, var, scale, Aniso, proj)
}
\arguments{ 
 \item{phi, fulldim, reduceddim, layers}{see \command{\link{RPtbm}}.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 \command{\link{RMtbm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\details{
 The turning bands method stems from the 1:1 correspondence between the
 isotropic covariance functions of different dimensions. See Gneiting
 (1999).

 The standard case \code{reduceddim=1} and \code{fulldim=3}.
 If only one of the parameters are given, then the difference of two
 parameters equals 2.
 
 For \code{d == n + 2}, where \code{n=reduceddim} and 
 \code{d==fulldim} the original dimension, we have

 \deqn{
 C(r) = \phi(r) + r \phi'(r) / n
 }{
 C(r) = phi(r) + r phi'(r) / n
 }
 which, for \code{n=1} reduced to the standard TBM operator
 \deqn{
 C(r) =\frac {d}{d r} r \phi(r)
 }{
 C(r) = d/dr [ r phi(r) ]
 }
 
 For \code{d == 2 && n == 1} we have
 
 \deqn{
 C(r) = \frac{d}{dr}\int_0^r \frac{u\phi(u)}{\sqrt{r^2 - u^2}} d u
 }{
 C(r) = d/dr int_0^r [ r phi(r) ] / [ sqrt{r^2 - u^2} ] d u
 }
 
}
\references{
 \itemize{
 \item Gneiting, T. (1999)
 On the derivatives of radial positive definite function.
 \emph{J. Math. Anal. Appl}, \bold{236}, 86-99

 \item
 Matheron, G. (1973).
 The intrinsic random functions and their applications.
 \emph{Adv . Appl. Probab.}, \bold{5}, 439-468.
 }
}
\seealso{
 \command{\link{RPtbm}},
 \command{\link{RFsimulate}}.
}
\examples{
RFoptions(seed=0)
x <- seq(0,25, if (interactive()) 0.02 else 5)
model <- RPtbm(RMspheric())
z <- RFsimulate(model, x, x)
plot(z)
\dontshow{RFoptions(seed=NA)}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}
\keyword{models}
