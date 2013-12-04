\name{RMintexp}
\alias{RMintexp}
\title{Integral exponential operator}
\description{
 \command{\link{RMintexp}} is a univariate stationary covariance model
 depending on a univariate variogram model \eqn{\phi}{phi}.
 The corresponding covariance function only depends on the difference
 \eqn{h}{h} between two points and is given by
 \deqn{C(h)=(1 - exp(-\phi(h)))/\phi(h)}{C(h)=(1 - exp(-phi(h)))/phi(h)}}
 \usage{
RMintexp(phi, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a variogram \command{\link{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
%\details{}
\value{
 \command{\link{RMintexp}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Schlather, M. (2012)
 Construction of covariance functions and unconditional simulation of
 random fields.
 \emph{Lecture Notes in Statistics, Proceedings}, \bold{207}, 25--54.
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
set.seed(0)
model <- RMintexp(RMfbm(alpha=1.5, scale=0.2))
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
}