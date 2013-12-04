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
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
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
set.seed(0)
model <- RMexp()
x <- 0:10
set.seed(0)
z <- RFsimulate(model, x, spC=FALSE)
model2 <- RMid(model)
set.seed(0)
z2 <- RFsimulate(model, x, spC=FALSE)
sum(abs(z-z2)) == 0 # TRUE

}


\keyword{spatial}
\keyword{models}

