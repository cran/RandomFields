\name{RRspherical}
\alias{RRspherical}
\title{Random scaling used with balls}
\description{
  This model describes the \bold{radius} of a \eqn{d} ball obtained by the
  intersection of a \code{balldim}-dimensional ball
  with \bold{diameter} 1 with a \eqn{d}-dimensional hyperplane that has uniform
  distance from the center.
 }
 \usage{
RRspherical(spacedim, balldim) 
}

\arguments{
  \item{spacedim}{dimension \eqn{d} of the hyperplane; defaults to 1.}
 \item{balldim}{the dimension of the ball}
}

\value{
 \command{\link{RRspherical}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
  \command{\link{RMmodel}},
  \command{\link{RMball}},
}

\examples{
RFoptions(seed=0)
hist(RFrdistr(RRspherical(balldim=2), n=1000), 50)
\dontshow{RFoptions(seed=NA)}
}

\keyword{spatial}
\keyword{models}




