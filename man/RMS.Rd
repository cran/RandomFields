\name{RMS}
\alias{RMS}
\title{Scaling operator}
\description{
 \command{\link{RMS}} is an operator that modifies the variance and the
 coordinates or distances of a submodel \eqn{\phi} by 
 \deqn{C(h) = v * \phi(A*h/s).}

 \bold{Most users will never call \command{\link{RMS}} directly, see
   Details.} However, the following describes the arguments
 \code{var}, \code{scale}, \code{Aniso}, \code{proj} that are common to
 nearly all models. See \code{\link{RMSadvanced}} for advanced use of
   these arguments.
}
\usage{
RMS(phi, var, scale, Aniso, proj, anisoT)
}
\arguments{
  \item{phi}{submodel}
  \item{var}{is the optional variance parameter \eqn{v}.
  }
  \item{scale}{scaling parameter \eqn{s} which is positive.
  }
  \item{Aniso}{matrix or \code{\link{RMmodel}}.
    The optional anisotropy matrix \eqn{A}, multiplied from
    the right by a distance vector \eqn{x}, i.e. \eqn{Ax}.
  } 
  \item{proj}{is the optional projection vector which defines a diagonal
    matrix of zeros and ones and \code{proj} gives the 
    positions of the ones (integer values between 1 and the dimension of
    \eqn{x}). It also allows for the values \code{'space'} and
    \code{'time'} in case of space-time modelling.
  }
  \item{anisoT}{the transpose of the anisotropy matrix \eqn{B},
    multiplied from the left by a distance vector \eqn{x}, i.e.
    \eqn{x^\top B}.
  }
}
\value{
 \command{\link{RMS}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\details{
 The call in the usage section is equivalent to 
 \code{phi(..., var, scale, anisoT, Aniso, proj)}, where \code{phi} has
 to be replaced by a valid \command{\link{RMmodel}}.

 Most users will never call \command{\link{RMS}} directly.
}
\note{At most one of the arguments \code{Aniso}, \code{anisoT} and \code{proj} may be given at the same time.
}

\me
\seealso{
  \command{\link{RMSadvanced}},
  \command{\link{RMmodel}},

  \command{\link{RMprod}} for an alternative way to define
  an arbitrary, location dependent variance. There, the standard
  deviation is given so that \command{\link{RMprod}} might be used
  even in the multivariate case.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model1 <- RMS(RMexp(), scale=2)
model2 <- RMexp(scale=2)
x <- seq(0, 10, 0.02)
print(all(RFcov(model1, x) == RFcov(model2, x))) # TRUE
\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
