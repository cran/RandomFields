\name{RMmqam}
\alias{RMmqam}
\alias{RMmqam}
\title{multivariate quasi-arithmetic mean}
\description{
 \command{\link{RMqam}} is a multivariate stationary covariance model depending
 on a submodel \eqn{\phi}{phi} such that
 \eqn{\psi(\cdot) := \phi(\sqrt(\cdot))}{psi( . ) := phi(sqrt( . ))}
 is completely monotone, and depending on further stationary
 covariance models \eqn{C_i}. The covariance is given by

 \deqn{C_{ij}(h) = \phi(\sqrt(\theta_i (\phi^{-1}(C_i(h)))^2 + \theta_j (\phi^{-1}(C_j(h)))^2 ))}
 where \eqn{\phi} is a completely monotone function, \eqn{C_i} are
 suitable covariance functions and \eqn{\theta_i\ge0} such that
 \eqn{\sum_i \theta_i=1}.
 
}
\usage{
RMmqam(phi, C1, C2, C3, C4, C5, C6, C7, C8, C9, theta, var, scale, Aniso, proj)
}
\arguments{
  \item{phi}{a valid covariance \command{\link{RMmodel}}
    that is a normal scale mixture,
    cf. \code{\link{RFgetModelNames}(monotone="normal mixture")}
  }
 \item{C1, C2, C3, C4, C5, C6, C7, C8, C9}{optional further stationary \command{\link{RMmodel}}}
 \item{theta}{is a vector of values in \eqn{[0,1]}, summing up to \eqn{1}.}
% \item{rho}{a matrix with positive entries}
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
 \command{\link{RMmqam}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Porcu, E., Mateu, J. & Christakos, G. (2009) Quasi-arithmetic
 means of covariance functions with potential applications to
 space-time data. \emph{Journal of Multivariate Analysis}, \bold{100},
 1830--1844.
 }
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMqam}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{\dontrun{ ## irgendo ein Fehler
model <- RMmqam(phi=RMgauss(), RMexp(), RMgauss(), theta=c(0.7, 0.3))
x <- seq(0, 4, if (interactive()) 0.3 else 1) 
plot(model)
plot(RFsimulate(model, x=x, print=10))
}}
\dontshow{FinalizeExample()}
}