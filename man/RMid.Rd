\name{RMid}
\alias{RMid}
\title{Identical Model}
\description{
  \command{RMid} is the identical operator
  for objects of class \code{\link[=RMmodel-class]{RMmodel}}
 }
 \usage{
RMid(phi, vdim, var, scale, Aniso, proj)
}

\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{vdim}{for internal purposes}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}

 
\value{
 \command{\link{RMid}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMexp()
x <- 0:10
z <- RFsimulate(model, x)

model2 <- RMid(model)
z2 <- RFsimulate(model, x)
sum(abs(as.vector(z)- as.vector(z2))) == 0 # TRUE
\dontshow{FinalizeExample()}

}


\keyword{spatial}
\keyword{models}

