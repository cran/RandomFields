\name{Kriging}
\alias{Kriging}
\title{Kriging methods}
\description{
  The function allows for different methods of kriging.
}
\usage{
Kriging(krige.method, x, y=NULL, z=NULL, T=NULL, grid,
        gridtriple=FALSE, model, param, given, data, trend, pch=".",
        return.variance=FALSE, internal=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{krige.method}{kriging method; currently only 'S' (simple
    kriging) and 'O' (ordinary kriging) implemented.}
  \item{x}{\eqn{(n \times d)}{(n x d)} matrix or vector of \code{x}
    coordinates; coordinates of \eqn{n} points to be kriged}
  \item{y}{vector of \code{y} coordinates.}
  \item{z}{vector of \code{z} coordinates.}
  \item{T}{vector in grid triple form for the time coordinates.}
  \item{grid}{logical; determines whether the vectors \code{x},
    \code{y}, and \code{z} should be
    interpreted as a grid definition, see Details.}
  \item{gridtriple}{logical.  Only relevant if \code{grid==TRUE}. 
    If \code{gridtriple==TRUE}
    then \code{x}, \code{y}, and \code{z} are of the
    form \code{c(start,end,step)}; if
    \code{gridtriple==FALSE} then \code{x}, \code{y}, and \code{z}
    must be vectors of ascending values.
  }
  \item{model}{string; covariance model, see \code{\link{CovarianceFct}}, or
    type \code{\link{PrintModelList}()} to get all options.}
  \item{param}{parameter vector:
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order.  Further parameters are to be added in case of a
    parametrised class of covariance functions, see
    \link{CovarianceFct}. 
    The value of \code{mean} must be finite
    in the case of simple kriging, and is ignored otherwise.}
  \item{given}{matrix or vector of points where data are available.}
  \item{data}{the data values given at \code{given}; it might be a
    vector or a matrix. If a matrix is given multivariate data are
    assumed which are kriged \emph{separately}.}
  \item{trend}{not programmed yet (will be used in case of universal kriging)}
  \item{pch}{Kriging procedures are quite time consuming in general. 
    The character \code{pch} is printed after roughly
    each 80th part of calculation.}
  \item{return.variance}{logical. If \code{FALSE} the kriged field is
    returned. If \code{TRUE} a list of two elements, \code{estim} and
    \code{var}, i.e. the kriged field and the kriging variances,
    is returned.}
  \item{internal}{
    \code{FALSE}. \code{internal} should not be set to \code{TRUE}
    by the user.
    (In case \code{internal=TRUE},
    various consistency checks for the input variables
    are not performed. Further, \code{grid} must be \code{FALSE},
    and \code{model} must be given in the output format of
    \code{\link{PrepareModel}}.)
  }
}
\details{
  \itemize{
    \item \code{grid==FALSE} : the vectors \code{x}, \code{y},
    and \code{z} are interpreted as vectors of coordinates
    \item \code{(grid==TRUE) && (gridtriple==FALSE)} : the vectors
    \code{x}, \code{y}, and \code{z}
    are increasing sequences with identical lags for each sequence. 
    A corresponding
    grid is created (as given by \code{expand.grid}). 
    \item \code{(grid==TRUE) && (gridtriple==TRUE)} : the vectors
    \code{x}, \code{y}, and \code{z}
    are triples of the form (start,end,step) defining a grid
    (as given by \code{expand.grid(seq(x$start,x$end,x$step),
      seq(y$start,y$end,y$step),
      seq(z$start,z$end,z$step))})
  }
}
\value{
  If \code{variance.return=FALSE} \code{Kriging} returns a vector or matrix
  of kriged values corresponding to the
  specification of \code{x}, \code{y}, \code{z}, and
  \code{grid}, and \code{data}.\cr
  
 \code{data}: a vector or matrix with \emph{one} column\cr
    * \code{grid==FALSE}.  A vector of simulated values is
    returned (independent of the dimension of the random field)\cr
    * \code{grid==TRUE}.  An array of the dimension of the
    random field is returned (according to the specification
    of \code{x}, \code{y}, and \code{z}).\cr
    
  \code{data}: a matrix with \emph{at least two} columns\cr
    * \code{grid==FALSE}.  A matrix with the \code{ncol(data)} columns
    is returned.\cr
    * \code{grid==TRUE}.  An array of dimension
    \eqn{d+1}{d+1}, where \eqn{d}{d} is the dimension of
    the random field, is returned (according to the specification
    of \code{x}, \code{y}, and \code{z}).  The last
    dimension contains the repetitions.

    If \code{variance.return=TRUE} a list of two elements, \code{estim} and
    \code{var}, i.e. the kriged field and the kriging variances,
    is returned. The format of \code{estim} is the same as described
    above.
    The format of \code{var} is accordingly.
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
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
%\note{}

\seealso{
  \code{\link{CondSimu}},
  \code{\link{CovarianceFct}},
  \code{\link{EmpiricalVariogram}},
  \code{\link{RandomFields}},
}

\examples{
## creating random variables first
## here, a grid is chosen, but does not matter
step <- 0.25 
x <-  seq(0,7,step)
param <- c(0,1,0,1)
model <- "exponential"
RFparameters(PracticalRange=FALSE)
p <- 1:7
points <- as.matrix(expand.grid(p,p))
data <- GaussRF(points, grid=FALSE, model=model, param=param)

## visualise generated spatial data
zlim <- c(-2.6,2.6)
colour <- rainbow(100)
image(p, p, xlim=range(x), ylim=range(x),
      matrix(data,ncol=length(p)),
      col=colour,zlim=zlim)

## now: kriging
krige.method <- "O" ## ordinary kriging
z <-  Kriging(krige.method=krige.method,
              x=x, y=x, grid=TRUE,
              model=model, param=param,
              given=points, data=data)
image(x,x,z,col=colour,zlim=zlim)
}
\keyword{spatial}%-- one or more ...


