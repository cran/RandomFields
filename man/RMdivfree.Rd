\name{RMdivfree}
\alias{RMdivfree}
\title{Divfree Covariance Model}
\description{
 \command{\link{RMdivfree}}
 is a multivariate covariance model which depends on
 a univariate stationary covariance model where the covariance function phi(h)
 is twice differentiable. 
 
 The corresponding matrix-valued covariance function C of the model
 only depends on the difference \eqn{h}{h} between two points and it is
 given by the following components,
 \itemize{
   \item the potential
   \item the vector field given by
   \deqn{
     C(h)=( - \Delta E + \nabla \nabla^T ) C_0(h)
   }
   \item the curl field
 }
 }
\usage{
RMdivfree(phi, which, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a univariate stationary covariance model (in 2 or 3 dimensions).}
 \item{which}{vector of integers. If not given all components are
   returned; otherwise the selected components are returned.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 The model returns the potential field in the first component, the
 corresponding divfree field and the field of curl strength in the last component.

 See also the models \command{\link{RMcurlfree}} and \command{\link{RMvector}}.
}
\value{
 \command{\link{RMdivfree}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Scheuerer, M. and Schlather, M. (2012)
 Covariance Models for Divergence-Free and Curl-Free Random Vector Fields.
 \emph{Stochastic Models} \bold{28:3}.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMcurlfree}},
 \command{\link{RMvector}},
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
model <- RMdivfree(RMgauss(), scale=4)
plot(model, dim=2)

x.seq <- y.seq <- seq(-10, 10, 0.2)
simulated <- RFsimulate(model=model, x=x.seq, y=y.seq)

plot(simulated)
plot(simulated, select.variables=1)
plot(simulated, select.variables=2:3)
plot(simulated, select.variables=list(2:3))
plot(simulated, select.variables=list(1, 2:3, 4))
plot(simulated, select.variables=list(1, c(1, 2:3), 4))

\dontshow{FinalizeExample()}
}
