\name{RMgneitingdiff}
\alias{RMgneitingdiff}
\title{Gneiting Covariance Model Used as Tapering Function}
\description{
 \command{\link{RMgneitingdiff}} is a stationary isotropic covariance model
 which is only valid up to dimension 3.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between 
 two points and is given by
 \deqn{C(h) = C_0(h / t) W_\nu(h / s)}
 where \eqn{C_0} is Gneiting's model \command{\link{RMgneiting}}
 and \eqn{W_\nu} is the Whittle model\command{\link{RMwhittle}}.
}
\usage{
RMgneitingdiff(nu, taper.scale, scale, var, Aniso, proj)
}
\arguments{
 \item{nu}{see \command{\link{RMwhittle}} }
 \item{taper.scale}{is the paramter \eqn{t} in the above formula}
 \item{scale}{is the paramter \eqn{s} in the above formula}
 \item{var,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 The model allows to a certain degree the smooth modelling of the
 differentiability of a covariance function with compact support.
}

\value{
 \command{\link{RMgneitingdiff}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
 }
\references{
 \itemize{
 \item Gneiting, T. (1999)
 Correlation functions for atmospherical data analysis.
 \emph{Q. J. Roy. Meteor. Soc} Part A \bold{125}, 2449-2464.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMbigneiting}},
 \command{\link{RMgneiting}},
 \command{\link{RMgengneiting}},
 \command{\link{RMgauss}},
 \command{\link{RMmodel}},
 \command{\link{RMwhittle}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
model <- RMgneitingdiff(nu=2, taper.scale=1, scale=0.2)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}