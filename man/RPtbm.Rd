\name{Tbm}
\alias{Tbm}
\alias{RPtbm}
\title{Turning Bands method}
\description{
 The Turning Bands method is a simulation method for stationary, isotropic
 random fields in any dimension and defined on arbitrary points or
 arbitrary grids. It performs a multidimensional simulation
 by superposing lower-dimensional fields. In fact, the Turning Bands
 method is called with the Turning Bands model, see  
 \command{\link{RMtbm}}.
 \cr
 For details see \command{\link{RMtbm}}.
}

\usage{
RPtbm(phi, boxcox, fulldim, reduceddim, layers, lines,
      linessimufactor, linesimustep, center, points)
}

\arguments{
  \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
    specifies the covariance function to be simulated;
    a univariate stationary isotropic covariance model
    (see \code{RFgetModelNames(type="positive definite",
      domain="single variable", isotropy="isotropic",  vdim=1)})
    which is valid in dimension \code{fulldim}.
  }
\item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  see \command{\link{RFboxcox}} for Details.
 }
  % \item{loggauss}{see \command{\link{RPgauss}}.}
  \item{fulldim}{a positive integer. The dimension of the space of the
    random field to be simulated}
  
  \item{reduceddim}{a positive integer; less than \code{fulldim}. 
    The dimension of the auxiliary hyperplane (most frequently a line,
    i.e. \code{reduceddim=1} used in the simulation.
  }
  
  \item{layers}{a boolean value; for space time model. If \code{TRUE}
    then the turning layers are 
    used whenever a time component is given.
    If \code{NA} the turning layers are used only when the
    traditional TBM is not applicable.
    If \code{FALSE} then turning layers may never be used.
    
    Default: \code{TRUE}.
  }
  
  \item{lines}{
    Number of lines used.
    Default: \code{60}.
  }
  
  \item{linessimufactor}{ \code{linessimufactor} or
    \code{linesimustep} must be non-negative; if
    \code{linesimustep}
    is positive then \code{linesimufactor} is ignored.
    If both
    arguments are naught then \code{points} is used (and must be
    positive).
    The grid on the line is \code{linesimufactor}-times
    finer than the smallest distance. 
    See also \code{linesimustep}.
    
    Default: \code{2.0}.
  }
  
  \item{linesimustep}{
    If \code{linesimustep} is positive the grid on the line has lag
    \code{linesimustep}. 
    See also \code{linesimufactor}.
    
    Default: \code{0.0}.
  }
  
  \item{center}{Scalar or vector.
    If not \code{NA}, the \code{center} is used as the center of
    the turning bands for \code{fulldim}.
    Otherwise the center is determined
    automatically such that the line length is minimal.
    See also \code{points} and the examples below.
    
    Default: \code{NA}.
  }
  
  \item{points}{integer. If greater than 0,
    \code{points} gives the number of points simulated on the TBM
    line, hence 
    must be greater than the minimal number of points given by
    the size of the simulated field and the two paramters
    \code{linesimufactor} and \code{linesimustep}.
    If \code{points} is not positive the number of points is
    determined automatically.
    The use of \code{center} and \code{points} is highlighted
    in an example below.
    
    Default: \code{0}.
  } 
}

 
\note{Both the precision and the simulation time
  depend heavily on \code{linesimustep} and
  \code{linesimufactor}.
  For covariance models with larger values of the scale parameter,
  \code{linesimufactor=2} is too small.
}

 

\details{
  \itemize{
    \item 2-dimensional case\cr
    It is generally difficult to use the turning bands method
    (\command{RPtbm}) directly in the 2-dimensional space.
    Instead, 2-dimensional random fields are frequently obtained
    by simulating a 3-dimensional random field (using
    \command{RPtbm}) and taking a 2-dimensional cross-section.
    See also the arguments \code{fulldim} and \code{reduceddim}.
 
    \item 4-dimensional case\cr
    The turning layers can be used for the simulations with a (formal)
    time component. It works for all isotropic models, 
    some special models such as \command{\link{RMnsst}}, and
    multiplicate models with that separate the time component. 
  }
}


\value{
 \code{RPtbm} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
  Turning bands
  \itemize{
    \item Lantuejoul, C. (2002)
    \emph{Geostatistical Simulation: Models and Algorithms.}
    Springer.
    
    \item
    Matheron, G. (1973).
    The intrinsic random functions and their applications.
    \emph{Adv. Appl. Probab.}, \bold{5}, 439-468.

    \item
    Strokorb, K., Ballani, F., and  Schlather, M. (2014)
    Tail correlation functions of max-stable processes: Construction principles, recovery and diversity of some mixing max-stable processes with identical TCF.
    \emph{Extremes}, \bold{} Submitted.
  }

  Turning layers
  \itemize{
    \item Schlather, M. (2011) Construction of covariance functions and
    unconditional simulation of random fields. In Porcu, E., Montero, J.M.
    and Schlather, M., \emph{Space-Time Processes and Challenges Related
      to Environmental Problems.} New York: Springer.
  }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
## isotropic example that forces the use of the turning bands method
model <- RPtbm(RMstable(s=1, alpha=1.8))
x <- seq(-3, 3, 0.1)
z <- RFsimulate(model=model, x=x, y=x)
plot(z)

## anisotropic example that forces the use of the turning bands method
model <- RPtbm(RMexp(Aniso=matrix(nc=2, rep(1,4))))
z <- RFsimulate(model=model, x=x, y=x)
plot(z)

## isotropic example that uses the turning layers method
model <- RMgneiting(orig=FALSE, scale=0.4)
x <- seq(0, 10, 0.1)
z <- RFsimulate(model, x=x, y=x, z=x, T=c(1,1,5))
plot(z, MARGIN.slices=4, MARGIN.movie=3)


\dontshow{FinalizeExample()}
}

\seealso{
  \command{\link{RMtbm}},
  \link{RP},
 \command{\link{RPspectral}}
}

\keyword{methods}
