\name{RMpower}
\alias{RMpower}
\alias{RMpower}
\title{Power operator for Variograms and Covariance functions}
\description{
 \command{\link{RMpower}} yields a variogram or covariance model 
 from a given variogram or covariance model.
 The variogram \eqn{\gamma}{gamma} of the model is given by
 \deqn{ \gamma = \phi^\alpha }{gamma=phi^alpha}
 if \eqn{\phi}{phi} is a variogram model.
 The covariance \eqn{C}{C} of the model is given by
 \deqn{ C(h) = \phi(0)-(\phi(0)-\phi(h))^\alpha }{C(h) = phi(0)-(phi(0)-phi(h))^alpha}
 if \eqn{\phi}{phi} is a covariance model.
}
\usage{
RMpower(phi, alpha, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a valid \command{\link{RMmodel}}; either a variogram model or a
 covariance model}
 \item{alpha}{a numerical value in the interval [0,1]}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any \command{\link{RMmodel}}. If not passed, the above covariance function remains unmodified.}
}
\details{
 If \eqn{\gamma} is a variogram, then \eqn{\gamma^\alpha} is a valid
 variogram for \eqn{\alpha} in the interval [0,1].
}
\value{
 \command{\link{RMpower}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 Schlather, M. (2012) Construction of covariance functions and
 unconditional simulation of random fields. In
 Porcu, E., Montero, J. M., Schlather, M. 
 \emph{Advances and Challenges in Space-time Modelling of Natural Events},
 Springer, New York.
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
model <- RMpower(RMgauss(), alpha=0.5)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}