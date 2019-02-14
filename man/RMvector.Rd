\name{RMvector}
\alias{RMvector}
\title{Vector Covariance Model}
\description{
 \command{\link{RMvector}} is a multivariate covariance model which depends on
 a univariate covariance model that is stationary in the first \eqn{Dspace}
 coordinates \eqn{h}{h} and where the covariance function phi(h,t)
 is twice differentiable in the first component \eqn{h}{h}. 
 
 The corresponding matrix-valued covariance function C of the model
 only depends on the difference \eqn{h}{h} between two points in the
 first component.
 It is given by
 \deqn{
 C(h,t)=( -0.5 * (a + 1) \Delta + a \nabla \nabla^T ) C_0(h, t)
 } 
 where the operator is applied to the first component \eqn{h}{h} only. 
}
\usage{
RMvector(phi, a, Dspace, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{an \command{\link{RMmodel}}; has two components \eqn{h}{h} (2- or 3-dimensional and stationary) and \eqn{t}{t} (arbitrary dimension).}
 \item{a}{a numerical value; should be in the interval \eqn{[-1,1]}.}
 \item{Dspace}{an integer; either 2 or 3; the first \eqn{Dspace}
 coordinates give the first component \eqn{h}{h}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 \eqn{C_0} is either a spatio-temporal model (then \eqn{t} is the
 time component) or it is an isotropic model. Then, the first \eqn{Dspace}
 coordinates are considered as \eqn{h} coordinates and the remaining
 ones as \eqn{t} coordinates. By default, \eqn{Dspace} equals the
 dimension of the field (and there is no \eqn{t} component).
 If \eqn{a=-1} then the field is curl free; if \eqn{a=1} then the field is
 divergence free.
}
\value{
 \command{\link{RMvector}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Scheuerer, M. and Schlather, M. (2012)
 Covariance Models for Divergence-Free and Curl-Free Random Vector Fields.
 \emph{Stochastic Models} \bold{28:3}.
 }
}

\me
\seealso{
 \command{\link{RMcurlfree}},
 \command{\link{RMdivfree}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMvector(RMgauss(), scale=0.3)
x <- seq(0, 10, 0.4)
plot(RFsimulate(model, x=x, y=x, z=0), select.variables=list(1:2))
\dontshow{FinalizeExample()}}