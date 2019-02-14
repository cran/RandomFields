\name{RMscale}
\alias{RMscale}
\title{Scale model for arbitrary areas of scales}
\description{
  Let \eqn{s_x} the scaling at location \eqn{x} and
  \eqn{p} a bijective penalizing function for (different) scales.
  Then covariance function is given by
  \deqn{C(x,y) = \phi(\|x-y\| + |p(s_x) - p(s_y)|)}
}
\usage{
RMscale(phi, scaling, penalty, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{isotropic submodel}
 \item{scaling}{model that gives the non-stationary scaling \eqn{s_x}}
 \item{penalty}{bijective function \eqn{p} applied to the \code{scaling}}
  \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 \command{\link{RMscale}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}

\references{
  \itemize{
    \item Bonat, W.H. , Ribeiro, P. Jr. and Schlather, M. (2019)
    Modelling non-stationarity in scale. In preparation.
  }
}


\me

\seealso{
  \command{\link{RMSadvanced}},
  \command{\link{RMblend}},
  \command{\link{RMbubble}}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- seq(0,1, 0.01)
scale <- RMcovariate(x=c(0,1), y=c(1,0),#2 areas separated by the 1st bisector
                     grid=FALSE, data=c(1, 3))

model <- RMscale(RMexp(), scaling = scale, penalty=RMid() / 2)
plot(z <- RFsimulate(model, x, x))


\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
