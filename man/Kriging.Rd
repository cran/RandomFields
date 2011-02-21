\name{Kriging}
\alias{Kriging}
\title{Kriging methods}
\description{
  The function allows for different methods of kriging.
}
\usage{
Kriging(krige.method, x, y=NULL, z=NULL, T=NULL, grid,
        gridtriple=FALSE, model, param, given, data, trend=NULL, pch=".",
        return.variance=FALSE, allowdistanceZero = FALSE, cholesky=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{krige.method}{kriging method; currently only 'S' (simple
    kriging), 'O' (ordinary kriging), 'U' (universal kriging) and 'I'
    (intrinsic kriging) implemented.}
  \item{x}{\eqn{(n \times d)}{(n x d)} matrix or vector of \code{x}
    coordinates; coordinates of \eqn{n} points to be kriged}
  \item{y}{vector of \code{y} coordinates.}
  \item{z}{vector of \code{z} coordinates.}
  \item{T}{vector in grid triple form for the time coordinates.}
  \item{grid}{logical; determines whether the vectors \code{x},
    \code{y}, and \code{z} should be
    interpreted as a grid definition, see Details.}
  \item{gridtriple}{logical.  Only relevant if \code{grid=TRUE}. 
    If \code{gridtriple=TRUE}
    then \code{x}, \code{y}, and \code{z} are of the
    form \code{c(start,end,step)}; if
    \code{gridtriple=FALSE} then \code{x}, \code{y}, and \code{z}
    must be vectors of ascending values.
  }
  \item{model}{string; covariance model, see \command{\link{CovarianceFct}}, or
    type \command{\link{PrintModelList}}\code{()} to get all options.}
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
  \item{trend}{only used for universal and intrinsic kriging. In case of
  universal kriging \code{trend} is a non-negative integer (monomials
  up to order k as trend functions), a list of functions or a formula (the
  summands are the trend functions); you have the choice of using either
  x, y, z or X1, X2, X3,... as spatial variables;
  in case of intrinsic kriging trend should be a nonnegative integer which
  is the order of the underlying model.
  }
  \item{pch}{Kriging procedures are quite time consuming in general. 
    The character \code{pch} is printed after roughly
    each 80th part of calculation.}
  \item{return.variance}{logical. If \code{FALSE} the kriged field is
    returned. If \code{TRUE} a list of two elements, \code{estim} and
    \code{var}, i.e. the kriged field and the kriging variances,
    is returned.}
  \item{allowdistanceZero}{if \code{TRUE} then
    identical locations are slightly scattered}
  \item{cholesky}{
    if \code{TRUE} cholesky decomposition is used instead of LU.
  }
}
\details{
  \itemize{
    \item \code{grid=FALSE} : the vectors \code{x}, \code{y},
    and \code{z} are interpreted as vectors of coordinates
    \item \code{(grid=TRUE) && (gridtriple=FALSE)} : the vectors
    \code{x}, \code{y}, and \code{z}
    are increasing sequences with identical lags for each sequence. 
    A corresponding
    grid is created (as given by \code{expand.grid}). 
    \item \code{(grid=TRUE) && (gridtriple=TRUE)} : the vectors
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
    * \code{grid=FALSE}.  A vector of simulated values is
    returned (independent of the dimension of the random field)\cr
    * \code{grid=TRUE}.  An array of the dimension of the
    random field is returned (according to the specification
    of \code{x}, \code{y}, and \code{z}).\cr
    
  \code{data}: a matrix with \emph{at least two} columns\cr
    * \code{grid=FALSE}.  A matrix with the \code{ncol(data)} columns
    is returned.\cr
    * \code{grid=TRUE}.  An array of dimension
    \eqn{d+1}{d+1}, where \eqn{d}{d} is the dimension of
    the random field, is returned (according to the specification
    of \code{x}, \code{y}, and \code{z}).  The last
    dimension contains the realisations.

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
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
%\note{}

\seealso{
  \command{\link{CondSimu}},
  \command{\link{Covariance}},
  \command{\link{CovarianceFct}},
  \command{\link{EmpiricalVariogram}},
  \code{\link{RandomFields}},
}

\examples{
###Example 1: Ordinary Kriging
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



%###Example 2: Universal Kriging (trend as formula)
%##creating random variables
%model <- "gauss"
%step <- pi/30
%#step <- pi/30 is nicer, but also time consuming
%param <- c(1,0,3)
%given <- expand.grid(x=seq(0,4*pi,pi/2), y=seq(0,4*pi, pi/2))
%data <- 0.5*(jitter(sin(rowSums(given)),factor=3)+1)
%x <- seq(0,4*pi,step)
%
%#Universal Kriging with trend as formula
%z <- Kriging("U", x=x, y=x, grid=TRUE, model=model, param=param, given=given,
% data=data, trend= z ~ sin(X1+X2) + 1, return.variance=FALSE)
%
%zlim <- c(min(z),max(z))
%colour <- rainbow(1000)
%image(x,x,z,col=colour,zlim=zlim)



%###Example 3: Universal Kriging (trend as list of functions)
%##creating random variables
%model <- "cubic"
%step <- pi/3
%#step <- pi/30 is nicer, but also time consuming
%param <- c(1,0,2)
%given <- expand.grid(x=seq(-2*pi,2*pi,pi/2), y=seq(-2*pi,2*pi, pi/2))
%data <- jitter(given[,1]^2 + given[,2]^2,factor=15)
%x <- seq(-2*pi,2*pi,step)

%#Universal Kriging with trend as list of functions
%z <- Kriging("U", x=x, y=x, grid=TRUE, model=model, param=param, given=given,
% data=data, trend=list(function(x,y) x^2 + y^2, function(x,y) (x^2 + y^2)^0.5, function(x,y) 1), return.variance=FALSE)
%
%zlim <- c(min(z),max(z))
%colour <- rainbow(1000)
%image(x,x,z,col=colour,zlim=zlim)



%###Example 4: Intrinsic Kriging
%model <- list("$", var=1, scale=4, list("fractalB", a=1))
%step=0.5
%given=expand.grid(x=seq(-3,3,0.25), y=seq(-3,3,0.25))
%data <- jitter((given[,1]+given[,2])^2, factor=10)
%x <- seq(-5,5,step)
%z <- Kriging("I", x=x, y=x, grid=TRUE, model=model, given=given, data=data, trend=2, return.variance=TRUE)
%Print(z)
%zlim <- c(min(z$estim),max(z$estim))
%colour <- rainbow(1000)
%image(x,x,z$estim,col=colour,zlim=zlim)
}
\keyword{spatial}%-- one or more ...
