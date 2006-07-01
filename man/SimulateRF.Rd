\name{SimulateRF}
\alias{SimulateRF}
\alias{InitSimulateRF}
\alias{DoSimulateRF}
\title{Simulation of Random Fields}
\description{
  \code{DoSimulateRF} performs an already initialised simulation.
  
  \code{InitSimulateRF} internal function;
  use \command{\link{InitGaussRF}} and \command{\link{InitMaxStableRF}}, instead.


% and some exotic methods, like hyperplane tessellations. 
}
\usage{
DoSimulateRF(n=1, register=0, paired=FALSE)

InitSimulateRF(x, y=NULL, z=NULL, T=NULL, grid, model, param, trend,
               method=NULL, register=0, gridtriple=FALSE, distribution=NA)
}
\arguments{
  \item{x}{matrix of coordinates, or vector of x coordinates}
  \item{y}{vector of y coordinates}
  \item{z}{vector of z coordinates}
  \item{T}{time instances}
  \item{grid}{logical; determines whether the vectors \code{x},
    \code{y}, and \code{z} should be
    interpreted as a grid definition, see Details.}
  \item{model}{string; covariance or variogram model,
    see \command{\link{CovarianceFct}}, or
    type \command{\link{PrintModelList}}\code{()} to get all options}
  \item{param}{vector or list. 
    \code{param=c(mean, variance, nugget, scale, ...)}, 
    \code{param=list(c(variance, scale,
      ...), ..., c(variance,scale,...))},
    \code{param=matrix(...)}, or
    \code{param=list(list(variance, anisotropy, kappa),...,
      list(variance, anisotropy, kappa))};
      the parameters must be given
      in this order; further parameters are to be added in case of a
      parametrised class of models, see \command{\link{CovarianceFct}}}
    \item{trend}{Not programmed yet.
      trend surface: number (mean),  vector of length
    \eqn{d+1} (linear trend \eqn{a_0+a_1 x_1 + \ldots + a_d x_d}{
      a_0 +a_1 x_1 + ... + a_d x_d}), or function}
  \item{method}{\code{NULL} or string; Method used for simulating,
    see \command{\link{RFMethods}}, or
    type \command{\link{PrintMethodList}}\code{()} to get all options}
  \item{register}{0:9; place where intermediate calculations are stored;
    the numbers are aliases for 10 internal registers}
  \item{gridtriple}{logical;  if \code{gridtriple=FALSE} ascending
    sequences for the parameters 
    \code{x}, \code{y}, and \code{z} are
    expected; if \code{gridtriple=TRUE} triples of form
    \code{c(start,end,step)} 
    expected; this parameter is used only
    if \code{grid=TRUE}}
  \item{distribution}{marginal distribution:\cr
    'Gauss', 'Poisson', or 'MaxStable'}
  \item{n}{number of realisations to generate; if \code{paired=TRUE}
    then \code{n} must be even.}
  \item{paired}{
    logical. \code{paired} may be \code{TRUE} only for the simulation of
    Gaussian random fields.
    If \code{TRUE} then every second simulation is obtained by
    only changing the signs of the standard Gaussian random variables, the
    simulation is based on (\dQuote{antithetic pairs}).
  }
}
\value{
  \code{InitSimulateRF} returns 0 if no error has occurred during the
  initialisation process, and a positive value if failed.\cr

  \code{DoSimulateRF} returns \code{NULL}
  if an error has occurred; otherwise the returned object
  depends on the parameters \code{n} and \code{grid}:\cr
    \code{n=1}:\cr
    * \code{grid=FALSE}.  A vector of simulated values is
    returned (independent of the dimension of the random field)\cr
    * \code{grid=TRUE}.  An array of the dimension of the
    random field is returned.\cr
    
    \code{n>1}:\cr
    * \code{grid=FALSE}.  A matrix is returned.  The columns
    contain the realisations.\cr
    * \code{grid=TRUE}.  An array of dimension
    \eqn{d+1}{d+1}, where \eqn{d}{d} is the dimension of
    the random field, is returned.  The last
    dimension contains the realisations.    
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\seealso{
  \command{\link{GaussRF}}, \command{\link{MaxStableRF}}, \code{\link{RandomFields}}
}
\keyword{spatial}



