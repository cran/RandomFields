\name{RMschlather}
\alias{RMschlather}
\title{Covariance Model for binary field based on Gaussian field}
\description{
  \command{RMschlather} gives
  the tail correlation function of the extremal Gaussian
 process, i.e.

 \deqn{C(h) = 1 - \sqrt{ (1-\phi(h)/\phi(0)) / 2 }}

 where \eqn{\phi} is the covariance of a stationary Gaussian field.
 
}
\usage{
RMschlather(phi, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model yields the tail correlation function of the field
 that is returned by \command{\link{RPschlather}}
}
\value{
 \command{\link{RMschlather}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}}
\seealso{
 \command{\link{RPschlather}}
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## This examples considers an extremal Gaussian random field
## with Gneiting's correlation function.

## first consider the covriance model and its corresponding tail
## corrlation function
model <- RMgneiting()
plot(model, model.tail.corr.fct=RMschlather(model),  xlim=c(0, 5))


## the extremal Gaussian field with the above underlying
## correlation function that has the above tail correlation function tcf
x <- seq(0, 10, if (interactive()) 0.1 else 3)
z <- RFsimulate(RPschlather(model), x)
plot(z)

## Note that in RFsimulate R-P-schlather was called, not R-M-schlather.
## The following lines give a Gaussian random field with corrlation
## function equal to the above tail correlation function.
z <- RFsimulate(RMschlather(model), x)
plot(z)


\dontshow{FinalizeExample()}

}

\keyword{spatial}
\keyword{models}
