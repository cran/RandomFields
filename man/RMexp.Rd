\name{RMexp}
\alias{RMexp}
\title{Exponential Covariance Model}
\description{
 \command{\link{RMexp}} is a stationary isotropic covariance model whose
 corresponding covariance function only depends 
 on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = e^(-r).}{C(r) = exp(-r).}
}
\usage{
RMexp(var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model is a special case of the Whittle covariance model (see \command{\link{RMwhittle}}) if \eqn{\nu=\frac{1}{2}}{\nu=0.5} and of the symmetric stable family (see \command{\link{RMstable}}) if \eqn{\nu = 1}{\nu=1}. Moreover,
 it is the continuous-time analog of the first order autoregressive time
 series covariance structure.

 The exponential covariance function is a normal scale mixture.
}
\value{
 \command{\link{RMexp}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{ Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.
 }

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{
 \command{\link{RMwhittle}},
 \command{\link{RMstable}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMexp()
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}