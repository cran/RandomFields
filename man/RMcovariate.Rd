\name{RMcovariate}
\alias{RMcovariate}
\alias{RMcovariateIntern}
\title{Model for covariates}
\description{
 The model makes covariates available
}
\usage{
RMcovariate(c, x, y=NULL, z=NULL, T=NULL, grid, factor, var)
RMcovariateIntern(norm, c, x, T, grid, factor)
}
\arguments{
  \item{c}{vector or matrix of data
  }
  \item{x,y,z,T,grid}{optional.
    The usual arguments as in \command{\link{RFsimulate}} to define the
    locations where the covariates are given
  }
  \item{factor}{vector or matrix of numerical values.
    The length of the vector must match the given number of locations.
  }
  \item{var}{variance, i.e. factor multiplied to the data, which can be
    estimated through ML
  }
  \item{norm}{optional model that gives the norm between locations}
  
  %\item{t}{takes the role \code{T} in the internal function}
}
\details{
  The functions interpolates (nearest neighbour) between the values.

  \command{RMcovariate_intern} should not be used by the user.
}
\value{
 \command{\link{RMcovariate}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
% 785 getnset

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

z <- 0.2 + (1:10)
RFfctn(RMcovariate(z, 1:10), c(2, 2.1, 2.5, 3))

RFfctn(RMcovariate(z), 1:10)

\dontshow{FinalizeExample()}
}

\keyword{spatial}
\keyword{models}
