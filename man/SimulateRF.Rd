\name{SimulateRF}
\alias{SimulateRF}
\alias{InitSimulateRF}
\alias{DoSimulateRF}
\title{Simulation of Random Fields}
\description{
  \code{DoSimulateRF} performs an already initialized simulation.
  
  \code{InitSimulateRF} internal function;
  use \code{\link{InitGaussRF}} and \code{\link{InitMaxStableRF}}, instead.

  \code{SimulateRF} internal function;
  use \code{\link{GaussRF}} and \code{\link{MaxStableRF}}, instead.

% and some exotic methods, like hyperplane tessellations. 
}
\usage{
DoSimulateRF(n=1, register=0)

InitSimulateRF(x, y=NULL, z=NULL, grid, model, param,
               method=NULL, register=0, gridtriple=FALSE,
               distribution=NA)

SimulateRF(x, y=NULL, z=NULL, grid, model, param, method=NULL,
          n=1, register=0, gridtriple=FALSE, distribution=NA)
}
\arguments{
  \item{x}{matrix of coordinates, or vector of x coordinates}
  \item{y}{vector of y coordinates}
  \item{z}{vector of z coordinates}
  \item{grid}{logical; determines whether the vectors \code{x},
    \code{y}, and \code{z} should be
    interpreted as a grid definition, see Details.}
  \item{model}{string; covariance or variogram model,
    see \code{\link{CovarianceFct}}, or
    type \code{\link{PrintModelList}()} to get all options}
  \item{param}{parameter vector:
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order; further parameters are to be added in case of a
    parametrised class of models, see \code{\link{CovarianceFct}}}
  \item{method}{\code{NULL} or string; Method used for simulating,
    see \code{\link{RFMethods}}, or
    type \code{\link{PrintMethodList}()} to get all options}
  \item{n}{number of realisations to generate}
  \item{register}{0:9; place where intermediate calculations are stored;
    the numbers are aliases for 10 internal registers}
  \item{gridtriple}{logical;  if \code{gridtriple==FALSE} ascending
    sequences for the parameters 
    \code{x}, \code{y}, and \code{z} are
    expected; if \code{gridtriple==TRUE} triples of form
    \code{c(start,end,step)} 
    expected; this parameter is used only
    if \code{grid==TRUE}}
  \item{distribution}{marginal distribution:\cr
    \code{"Gauss"}, \code{"Poisson"}, or \code{"MaxStable"}.

  }
}
\value{
  \code{InitSimulateRF} returns 0 if no error has occured during the
  initialisation process, and a positive value
  if failed.\cr

  \code{SimulateRF} and \code{DoSimulateRF} return \code{NULL}
  if an error has occured; otherwise the returned object
  depends on the parameters \code{n} and \code{grid}:\cr
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
    the random field, is returned.  The last
    dimension contains the repetitions.    
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
\seealso{
  \code{\link{GaussRF}}, \code{\link{MaxStableRF}}, \code{\link{RandomFields}}
}
\keyword{spatial}



