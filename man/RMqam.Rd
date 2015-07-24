\name{RMqam}
\alias{RMqam}
\alias{RMqam}
\title{Quasi-arithmetic mean}
\description{
 \command{\link{RMqam}} is a univariate stationary covariance model depending
 on a submodel \eqn{\phi}{phi} such that
 \eqn{\psi(\cdot) := \phi(\sqrt(\cdot))}{psi( . ) := phi(sqrt( . ))}
 is completely monotone, and depending on further stationary
 covariance models \eqn{C_i}. The covariance is given by
 
 \deqn{C(h) = \phi(\sqrt(\sum_i \theta_i (\phi^{-1}(C_i(h)))^2))}{C(h)
 = phi(sqrt(sum_i theta_i (phi^{-1}(C_i(h)))^2))}
}

\usage{
RMqam(phi, C1, C2, C3, C4, C5, C6, C7, C8, C9, theta, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a valid covariance \command{\link{RMmodel}} that is a normal
   scale mixture. See, for instance,
   \cr
   \code{\link{RFgetModelNames}(monotone="normal mixture")}
   }
 \item{C1, C2, C3, C4, C5, C6, C7, C8, C9}{optional further univariate
 stationary \command{\link{RMmodel}}.}
 \item{theta}{a vector with positive entries}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any \command{\link{RMmodel}}. If not passed, the above covariance function remains unmodified.}
}
\details{
 Note that \eqn{\psi(\cdot) := \phi(\sqrt(\cdot))}{psi( . ) :=
 phi(sqrt( . ))} is completely monotone if and only if
 \eqn{\phi}{phi} is a valid covariance function for all dimensions,
 e.g. \command{\link{RMstable}}, \command{\link{RMgauss}}, \command{\link{RMexponential}}.
 
 Warning: \code{RandomFields} cannot check whether the combination
 of \eqn{\phi}{phi} and \eqn{C_i} is valid.
}
\value{
 \command{\link{RMqam}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Porcu, E., Mateu, J. & Christakos, G. (2007) Quasi-arithmetic
 means of covariance functions with potential applications to
 space-time data. Submitted to Journal of Multivariate Analysis.
 }
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmqam}},
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
model <- RMqam(phi=RMgauss(), RMexp(), RMgauss(),
               theta=c(0.3, 0.7), scale=0.5)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}
