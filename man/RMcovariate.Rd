\name{RMcovariate}
\alias{RMcovariate}
\title{Model for covariates}
\description{
 The model makes covariates available.
}
\usage{
RMcovariate(c, x, y=NULL, z=NULL, T=NULL, grid, var, scale,
                        Aniso, proj,  raw, norm, addNA)
}
\arguments{
  \item{scale, Aniso, proj, var}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}

  \item{c}{vector or matrix of data
  }
  \item{x,y,z,T,grid}{optional.
    The usual arguments as in \command{\link{RFsimulate}} to define the
    locations where the covariates are given
  }
%  \item{factor}{vector or matrix of numerical values.
%    The length of the vector must match the given number of locations.
%  }
% 
%  \item{var}{variance, i.e. factor multiplied to the data, which can be
%    estimated through ML
%  }
  \item{raw}{
    logical. If \code{FALSE} then the data are interpolated. This
    approach is always save, but might be slow.

    If \code{TRUE} then the data may be accessed when covariance
    matrices are calculated. No rescaling or anisotropy definition
    is allowed in combination with the model. The use is dangerous,
    but fast.

    Default: FALSE
  }
  \item{norm}{optional model that gives the norm between locations}
  \item{addNA}{
    If \code{addNA} is \code{TRUE}, then an additional
    (linear) factor is estimated in an estimation framework.
    This parameter must be set in particular when \command{RMcovariate}
    passes several covariates.
  }
}
\note{
  \itemize{
    \item   \code{c}, \code{x}, 
    also accept lists of data. However, its use is not in an advanced
    stage yet.
   }
}


\details{
  The functions interpolates (nearest neighbour) between the values.
}
\value{
 \command{\link{RMcovariate}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMfixcov}},
 \command{\link{RMmodel}},
 \command{\link{RMtrend}}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

z <- 0.2 + (1:10)
RFfctn(RMcovariate(z), 1:10)
RFfctn(RMcovariate(z, 1:10), c(2, 2.1, 2.5, 3))


\dontshow{FinalizeExample()}
}

\keyword{spatial}
\keyword{models}
