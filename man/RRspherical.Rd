\name{RRspheric}
\alias{RRspheric}
\title{Random scaling used with balls}
\description{
  This model delivers the distribution
  of the \bold{radius} of a ball 
  obtained by the
  intersection of a \code{balldim}-dimensional ball
  with \bold{diameter} R by a \eqn{spacedim}-dimensional
  hyperplane that has uniform
  distance from the center.
 }
 \usage{
RRspheric(spacedim, balldim, R) 
}

\arguments{
  \item{spacedim}{dimension of the hyperplane; defaults to 1.}
  \item{balldim}{the dimension of the ball}
  \item{R}{radius. Default: 1}
}

\value{
 \command{\link{RRspheric}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
  \command{\link{RMmodel}},
  \command{\link{RMball}},
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
hist(RFrdistr(RRspheric(balldim=2), n=1000), 50)
\dontshow{FinalizeExample()}
}

\keyword{spatial}
\keyword{models}




