\name{EmpiricalVariogram}
\alias{EmpiricalVariogram}
\title{Empirical (Semi-)Variogram}
\description{
  \code{EmpiricalVariogram} calculates the empirical (semi-)variogram
  of a random field realisation
}
\usage{
EmpiricalVariogram(x, y=NULL, z=NULL, data, grid, bin, gridtriple=FALSE)}
\arguments{
  \item{x}{vector of coordinates}
  \item{y}{vector of coordinates}
  \item{z}{vector of coordinates}
  \item{data}{vector or matrix of data}
  \item{grid}{logical; if \code{TRUE} then
    \code{x}, \code{y}, and \code{z} define a grid; otherwise
    \code{x}, \code{y}, and \code{z} are interpreted as points}
  \item{bin}{vector of ascending values giving the bin boundaries}
  \item{gridtriple}{logical. Only relevant if \code{grid==TRUE}.
    If \code{gridtriple==TRUE}
    then \code{x}, \code{y}, and \code{z} are of the
    form \code{c(start,end,step)}; if
    \code{gridtriple==FALSE} then \code{x}, \code{y}, and \code{z}
    must be vectors of ascending values}
}
\value{
  The function returns
  \code{list(centers,emp.vario)} where \code{centers} are the central
  points of the bins and \code{emp.vario} gives the empirical variogram.
  Both elements are
  vectors of length \code{(length(bin)-1)}.
}
\details{
%  A version of the algorithm is custom-tailored to large grids with 
%  constant grid length in each direction), and invoked in such cases.  

  Comments on specific parameters:
  \itemize{
    \item \code{data}: the number of values must match the number of
    points (given by \code{x}, \code{y}, \code{z}, \code{grid}, and
    \code{gridtriple}).  That is, it must equal the number of points or be
    a multiple of it.  In case the number of data equals \eqn{n}{n} times the
    number of points, the data are interpreted as \eqn{n}{n} independent
    realisations for the given set of points. 
    \item \code{(grid==FALSE)}: the vectors \code{x}, \code{y}, and
    \code{z}, are interpreted as
    vectors of coordinates
    \item \code{(grid==TRUE) && (gridtriple==FALSE)}: the vectors \
    code{x}, \code{y}, and
    \code{z}
    are increasing sequences with identical lags for each sequence.
    A corresponding
    grid is created (as given by \code{expand.grid}).
    \item \code{(grid==TRUE) && (gridtriple==FALSE)}: the vectors
    \code{x}, \code{y}, and \code{z}
    are triples of the form (start,end,step) defining a grid (as given by
    \code{expand.grid(seq(x$start,x$end,x$step),
      seq(y$start,y$end,y$step),
      seq(z$start,z$end,z$step))})
    \item 
    The bins are left open, right closed intervals, i.e.,
    \eqn{(b_i,b_{i+1}]}{(\code{b[i]},\code{bin[i+1]}]} for
    \eqn{i=1,\ldots,}{i=1,...,}\code{length(bin)}\eqn{-1}{-1}. 
    Hence, to include zero, \code{bin[1]} must be negative.
  }
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}
}
\seealso{\code{\link{GaussRF}} and \code{\link{RandomFields}}}
\examples{
  #############################################################
  ## this example checks whether a certain simulation method ##
  ## works well for a specified covariance model and         ##
  ## a configuration of points                               ##
  #############################################################
  x <- seq(0, 10, 0.5)
  y <- seq(0, 10, 0.5)
  grid <- TRUE
  gridtriple <- FALSE   ## see help("GaussRF")
  model <- "wh"         ## whittlematern
  alpha <- 2
  mean <- 1
  variance <- 10
  nugget <- 5
  scale <- 2
  method <- "TBM3"
  bins <- seq(0, 5, 0.001)
  repetition <- 20 ## by far too small to get reliable results!!
                   ## It should be of order 500,
                   ## but then it will take some time
                   ## to do the simulations
  param <- c(mean, variance, nugget, scale, alpha)
  f <- GaussRF(x=x, y=y, grid=grid, gridtriple=gridtriple,
                  model=model, param=param, meth=method,
                  n=repetition)
  binned <- EmpiricalVariogram(x=x, y=y, data=f,
                 grid=grid, gridtriple=gridtriple, bin=bins)
  truevariogram  <- Variogram(binned$c, model, param)
  matplot(binned$c, cbind(truevariogram,binned$e), pch=c("*","e"))
  ##black curve gives the theoretical values
%  #############################################################
%  ## this example shows under which conditions the algorithm ##
%  ## still works reasonably fast (about 15 sec. at 500 MHz), ##
%  ## using all information available (and not only a random  ##
%  ## subset of the points)                                   ##
%  #############################################################
%  grid <- TRUE; 
%  y <- x <- seq(0,50,  0.1) ## large grid (of size 500 x 500)
%  bins <-   seq(0, 2,0.001) ## 2 << 50  -- this is the key point
%
%  model <- "sph"         ## spherical
%  mean <- 1
%  variance <- 10
%  nugget <- 5
%  scale <- 1.5
%  param <- c(mean, variance, nugget, scale)
%  f <- GaussRF(x=x, y=y, grid=grid, model=model, param=param)
%  binned <-
%      EmpiricalVariogram(x=x, y=y, data=f, grid=grid, bin=bins)
%  truevariogram  <- Variogram(binned$c, model, param)
%  matplot(midbin, cbind(truevariogram,binned$e), pch=c("*","e")) 
%
}
\keyword{spatial}



