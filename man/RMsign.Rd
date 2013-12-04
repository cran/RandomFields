\name{RMsign}
\alias{RMsign}
\title{Random sign}
\description{
  \command{RMsign} defined a random sign. It can be used as part of the
  model definition of a Poisson field. 

 }
\usage{
RMsign(phi, p)
}

\arguments{
 \item{phi}{shape function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{p}{probability of keeping the sign}
}

\details{
  \command{RMsign} changes the sign of the shape function \code{phi}
  with probability \code{1-p} and keeps it otherwise.
}
 
\value{
 \command{\link{RMsign}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RMmatrix}},
 }


\keyword{spatial}
\keyword{models}




\examples{
set.seed(0)
model <- RPpoisson(RMsign(RMtent(), p=0.8))
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(RFsimulate(model, x=x))
}