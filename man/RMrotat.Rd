\name{RMrotat}
\alias{RMrotat}
\alias{RMrotation}
\title{Rotation matrices}
\description{

  \command{RMrotat} and \command{RMrotation} are
  auxiliary space-time functions that create some rotation
 \deqn{f(h, t) = s (\cos(\phi t) h_1 + \sin(\phi t) h_2) / \|h\|}
 and
 \deqn{f(h, t) = (\cos(\phi t) h_1 + \sin(\phi t) h_2,
 - \sin(\phi t) h_1 + \cos(\phi t) h_2, t)}
respectively
}

 \usage{
RMrotat(speed, phi) 
RMrotation(phi)
}

\arguments{
 \item{speed}{real value \eqn{s} }
 \item{phi}{angle}
}

\details{
 \command{\link{RMrotat}} and \command{\link{RMrotation}} are
 space-time models for 
 two-dimensial space.
}

 
\value{
 \command{\link{RMrotat}} and \command{\link{RMrotation}} return an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{  
  \command{\link{RMmodel}},
  \link{S10}
 }


\keyword{spatial}
\keyword{models}


\examples{
# see S10
\dontshow{FinalizeExample()}
}
