\name{CondSimu}
\alias{CondSimu}
\title{Conditional Simulation}
\description{
  the function returns conditional simulations of a random field
}
\usage{
CondSimu(method, x, y=NULL, z=NULL, grid, gridtriple=FALSE, n=1,
         model, param, given, data, tol=1e-05)
}
\arguments{
  \item{method}{Assumption for the random field which corresponds to
    the respective kriging method;
    currently only \code{"S"} (simple
    kriging) and \code{"O"} (ordinary kriging) possible.}
  \item{x}{\eqn{(n \times d)}{(n x d)} matrix or vector of \code{x} coordinates; coordinates of
    \eqn{n} points to be kriged}
  \item{y}{vector of \code{y} coordinates}
  \item{z}{vector of \code{z} coordinates}
  \item{grid}{logical; determines whether the vectors \code{x},
    \code{y}, and \code{z} should be
    interpreted as a grid definition, see Details.}
  \item{gridtriple}{logical. Only relevant if \code{grid==TRUE}.
    If \code{gridtriple==TRUE}
    then \code{x}, \code{y}, and \code{z} are of the
    form \code{c(start,end,step)}; if
    \code{gridtriple==FALSE} then \code{x}, \code{y}, and \code{z}
    must be vectors of ascending values
  }
  \item{n}{number of realisations to generate}
  \item{model}{string; covariance model, see \code{\link{CovarianceFct}}, or
    type \code{\link{PrintModelList}()} to get all options}
  \item{param}{parameter vector:
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order; further parameters are to be added in case of a
    parametrised class of covariance functions,
    see \code{\link{CovarianceFct}};
    the value of \code{mean} must be finite
    in the case of simple kriging, and is ignored otherwise}
  \item{given}{matrix or vector of locations where data are available;
    note that it is not possible to give the points in form of a grid
    definition}
  \item{data}{the values measured}
  \item{tol}{considered only if \code{grid=TRUE};
    tolerated distances of a given point to the nearest grid point to
  be regarded as being zero; see Details}
}
\details{
  The same way as \code{GaussRF} the function
  \code{CondSimu} allows for simulating on grids or arbitrary
  locations. However simulation on a grid is sometimes performed
  as if the points were at arbitrary locations, what may
  imply a great reduction in speed! This happens when the \code{given}
  locations do not ly on the specified grid, since in an intermediate
  step simulation has to be performed simultaneously on both the grid
  defined by \code{x}, \code{y}, \code{z}, and the locations
  of \code{given}.\cr

  Comments on specific parameters
  \itemize{
    \item \code{grid==FALSE} : the vectors \code{x}, \code{y},
    and \code{z} are interpreted as vectors of coordinates
    \item \code{(grid==TRUE) && (gridtriple==FALSE)} : the vectors
    \code{x}, \code{y}, and \code{z}
    are increasing sequences with identical lags for each sequence. 
    A corresponding
    grid is created (as given by \code{expand.grid}). 
    \item \code{(grid==TRUE) && (gridtriple==FALSE)} : the vectors
    \code{x}, \code{y}, and \code{z}
    are triples of the form (start,end,step) defining a grid
    (as given by \code{expand.grid(seq(x$start,x$end,x$step),
      seq(y$start,y$end,y$step),
      seq(z$start,z$end,z$step))})
  }
}
value{
  The returned object depends on the parameters \code{n} and
  \code{grid}:\cr
    \code{n==1}:\cr
    * \code{grid==FALSE}.  A vector of simulated values is
    returned (independent of the dimension of the random field)\cr
    * \code{grid==TRUE}.  An array of the dimension of the
    random field is returned.\cr
    
    \code{n>1}:\cr
    * \code{grid==FALSE}.  A matrix is returned.  The columns
    contain the repetitions.\cr
    * \code{grid==TRUE}.  An array of dimension
    \eqn{d+1}{d+1}, where \eqn{d}{d} is the dimension of
    the random field as given by \code{x}, \code{y}, and \code{z},
    is returned.  The last dimension contains the repetitions.
}
\references{
 Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.

 Cressie, N.A.C. (1993)
 \emph{Statistics for Spatial Data.}
 New York: Wiley.
 
 Goovaerts, P. (1997) \emph{Geostatistics for Natural Resources
   Evaluation.} New York: Oxford University Press.
 
 Wackernagel, H. (1998) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 2nd edition.  
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
%\note{}


\seealso{
  \code{\link{CovarianceFct}},
  \code{\link{GaussRF}},
  \code{\link{Kriging}}
  \code{\link{RandomFields}},
}

\examples{
## creating random variables first
## here, a grid is chosen, but any arbitrary points for which
## data are given are fine. Indeed if the data are given on a
## grid, the grid has to be expanded before calling `CondSimu',
## see below.
## However, locations where values are to be simulated,
## should be given in form of a grid definition whenever
## possible 
param <- c(0, 1, 0, 1)
model <- "exponential"
RFparameters(PracticalRange=FALSE)$null
p <- 1:7
data <- GaussRF(x=p, y=p, grid=TRUE, model=model, param=param)

# another grid, where values are to be simulated
step <- 0.25 # or 0.3
x <-  seq(0, 7, step)

# standardisation of the output
lim <- range( c(x, p) )
zlim <- c(-2.6, 2.6)
colour <- rainbow(100)

## visualise generated spatial data
image(p, p, data, xlim=lim, ylim=lim, zlim=zlim, col=colour)

#conditional simulation
method <- "O" ## random field assumption corresponding to
           ## those of ordinary kriging
cz <- CondSimu(method=method, x, x,  grid=TRUE,
               model=model, param=param,
               given=expand.grid(p,p),## data are given on a grid
                                      ## expand the grid first
               data=data)

image(x, x, cz, col=colour, xlim=lim, ylim=lim, zlim=zlim)
}
\keyword{spatial}
