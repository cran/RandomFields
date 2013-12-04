\name{RMcoord}
\alias{RMcoord}
\title{Coordinates in the mixed model representation}
\description{
 This function is used to define, in the extremely rare case,
 coordinates that differ from the original coordinates to define
 a covariance matrix for a random effect model

 \code{~ RMmodel() + Z @ RMcoord(coord=X, RMmodel2())}
 
}

\usage{
RMcoord(C0, coord, dist)
}

\arguments{
 \item{C0}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{coord, dist}{either coordinates or a
 the lower matrix of a distance matrix can be passed
 }
}

\value{
 \command{\link{RMcoord}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
set.seed(0)
 ## For examples see the help page of 'RFformula' ##
}