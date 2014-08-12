\name{RMqexp}
\alias{RMqexp}
\title{Variant of the exponential model}
\description{
 The covariance function is

 \deqn{C(x)= ( 2 e^{-x} - \alpha e^{-2x} ) / ( 2 - \alpha )}
 
 }
 \usage{
RMqexp(alpha, var, scale, Aniso, proj)
}

\arguments{
 \item{alpha}{value in \eqn{[0,1]}}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}

\value{
 \command{\link{RMqexp}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
\itemize{
 \item ?
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 }


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMqexp(alpha=0.95, scale=0.2)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}