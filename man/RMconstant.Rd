\name{RMconstant}
\alias{RMconstant}
\alias{RMconstant}
\title{Covariance Matrix Constant in Space}
\description{
  \command{\link{RMconstant}} defines a spatially constant covariance function
}
\usage{
RMconstant(M, var)
}
\arguments{
 \item{M}{a numerical matrix defining the user-defined covariance for a
 random field; The matrix should be positive definite, symmetric and
 its dimension should be equal to the length of observation or
 simulation vector.}
\item{var}{ variance }
 }

\value{
 \command{\link{RMconstant}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
  
 \command{\link{RMfixcov}},
\command{\link{RMmodel}},
}


\keyword{spatial}
\keyword{models}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## to do

\dontshow{FinalizeExample()}
}