\name{RMSadvanced}
\alias{RMSadvanced}
\title{Scaling operator -- comments for advanced applications}
\description{
  Here advances uses are given for the arguments
  \code{var}, \code{scale}, \code{Aniso}, \code{proj} that are available
  to most of the models
}
\section{Usage}{
RMS(phi, var, scale, Aniso, proj, anisoT)
}
\section{Arguments}{
  \describe{
    \item{phi}{submodel}
    \item{var}{Instead of a constant it can be 
      also an arbitrary non-negative function, see \code{\link{R.}}
      and \command{\link{RMuser}} for defining arbitrary functions.
    }
    \item{scale}{instead of a positive constant  it can be an arbitrary,
      positive 
      deterministic function. In case of the latter, the scale should be
      given by one of the functions \command{\link{RMbubble}} or
      \command{\link{RMscale}}. In case none of them are given,
      \command{\link{RMscale}} is assumed with scale penality
      \eqn{\|s(x) - s(y)\|^2} for the square of the norm.
      
      The scale can be also a random variable
      in case of Bayesian modelling.
    }
    \item{Aniso}{matrix or \code{\link{RMmodel}}.
      Instead of a matrix, \code{Aniso} can be an arbitrary, vector-valued
      function .
    } 
    \item{proj}{is the optional projection vector which defines a diagonal
      matrix of zeros and ones and \code{proj} gives the 
      positions of the ones (integer values between 1 and the dimension of
      \eqn{x}). It also allows for the values \code{'space'} and
      \code{'time'} in case of space-time modelling.
    }
    
    \item{anisoT}{the transpose of the anisotropy matrix \eqn{B},
      multiplied from the left by a distance vector \eqn{x}, i.e. \eqn{x^\top B}.
    }
  }
}

\section{Details}{
  See the reference for Gneitings nsst model used for modelling scales.
  See also the example below.
}

\section{Value}{
 \command{\link{RMSadvanced}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}

\references{
  \itemize{
    \item Bonat, W.H. , Ribeiro, P. Jr. and Schlather, M. (2019)
    Modelling non-stationarity in scale. In preparation.
  }
}


\me

\seealso{
  \command{\link{RMS}},
  \command{\link{RMblend}} for a different approach on modelling
  different scales 
}

\examples{\dontshow{StartExample(FALSE)}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- seq(0,1, if (interactive()) 0.01 else 0.5)
d <- sqrt(rowSums(as.matrix(expand.grid(x-0.5, x-0.5))^2))
d <- matrix(d < 0.25, nc=length(x))
image(d)

scale <- RMcovariate(data=as.double(d) * 2 + 0.5, raw=TRUE)

S <- RMexp(scale = scale)
plot(zS <- RFsimulate(S, x, x))
CS <- RFcovmatrix(S, x, x)


\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
