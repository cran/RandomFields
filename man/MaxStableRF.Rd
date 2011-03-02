\name{MaxStableRF}
\alias{MaxStableRF}
\alias{InitMaxStableRF}
\title{Max-Stable Random Fields}
\description{
  These functions simulate stationary and isotropic max-stable
  random fields with unit Frechet margins.
}
\usage{
MaxStableRF(x, y=NULL, z=NULL, grid, model, param, maxstable,
            method=NULL, n=1, register=0, gridtriple=FALSE,...)

InitMaxStableRF(x, y=NULL, z=NULL, grid, model, param, maxstable,
               method=NULL, register=0, gridtriple=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{matrix of coordinates, or vector of x coordinates}
  \item{y}{vector of y coordinates}
  \item{z}{vector of z coordinates}
  \item{grid}{logical; determines whether the vectors \code{x},
    \code{y}, and \code{z} should be
    interpreted as a grid definition, see Details.}
  \item{model}{string; see \command{\link{CovarianceFct}}, or
    type \command{\link{PrintModelList}}\code{()} to get all options;
    interpretation depends on the value of \code{maxstable},
    see Details.}
  \item{param}{parameter vector:
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order; further parameters are to be added in case of a
    parametrised class of covariance functions,
    see \command{\link{CovarianceFct}}, or be given in one of the extended
    forms, see Details}
  \item{maxstable}{string. Either 'extremalGauss' or
    'BooleanFunction'; see Details.}
  \item{method}{\code{NULL} or string; method used for simulating,
    see \command{\link{RFMethods}}, or
    type \command{\link{PrintMethodList}}\code{()} to get all options;
    interpretation depends on the value of \code{maxstable}.}
  \item{n}{number of realisations to generate}
  \item{register}{0:9; place where intermediate calculations are stored;
    the numbers are aliases for 10 internal registers}
  \item{gridtriple}{logical;  if \code{gridtriple=FALSE} ascending
    sequences for the parameters 
    \code{x}, \code{y}, and \code{z} are
    expected; if \code{gridtriple=TRUE} triples of form
    \code{c(start,end,step)} 
    expected; this parameter is used only
    if \code{grid=TRUE}}
  \item{...}{\command{\link{RFparameters}} that are locally used only.}
}
\details{
  There are two different kinds of models for max-stable processes
  implemented:
  \itemize{
    \item \code{maxstable="extremalGauss"}\cr
    Gaussian random fields are multiplied by independent
    random factors,
    and the maximum is taken. The random factors are such that
    the resulting random field has unit
    Frechet margins; the specification of the random factor
    is uniquely given by the specification of the random
    field. The parameter vector \code{param}, the \code{model},
    and the \code{method} are interpreted
    in the same way as for Gaussian random fields, see
    \command{\link{GaussRF}}.
    
    \item \code{maxstable="BooleanFunction"}\cr
    Deterministic or random, upper semi-continuous
    \eqn{L_1}{L1}-functions are randomly centred and multiplied by
    suitable, independent random factors; the pointwise maximum over all
    these functions yields a max-stable random field.
    The simulation technique is related to the random coin
    method for Gaussian random field simulation,
    see \command{\link{RFMethods}}. Hence, only
    models that are suitable for the random coin method
    are suitable for this technique, see \command{\link{PrintModelList}}\code{()}
    for a complete list of suitable covariance models.\cr
    The only value allowed for \code{method} is 'max.MPP' (and
    \code{NULL}),
    see \command{\link{PrintMethodList}()}. In the parameter list
    \code{param} the first two entries, namely \code{mean} and
    \code{variance}, are ignored. If the nugget is positive,
    for each point an additional independent unit Frechet variable
    with scale parameter
    \code{nugget} is involved when building the maximum
    over all functions.

    The model may be defined alternatively in one of the two extended
    ways as introduced in \code{CovarianceFct} and \code{GaussRF}.
    However only a single model may be given! The model may be
    anisotropic.
  }
} \value{
  \code{InitMaxStableRF} returns 0 if no error has occurred, and
  a positive value if failed.\cr

  \code{MaxStableRF} and \command{\link{DoSimulateRF}} return \code{NULL}
  if an error has occurred; otherwise the returned object
  depends on the parameters:\cr
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
\references{
   Schlather, M. (2002) Models for stationary max-stable
   random fields. \emph{Extremes} \bold{5}, 33-44.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\seealso{
  \command{\link{CovarianceFct}},
  \link{sophisticated},
  \command{\link{GaussRF}},
  \code{\link{RandomFields}},
  \command{\link{RFMethods}},
  \command{\link{RFparameters}},
  \command{\link{DoSimulateRF}},
.
}
\examples{
n <- 30 ## nicer, but time consuming if n <- 100
x <- y <- 1:n
ms0 <- MaxStableRF(x, y, grid=TRUE, model="exponen",
                 param=c(0,1,0,40), maxstable="extr",
                 CE.force = TRUE)
image(x,y,ms0)



% par(mfcol=c(2,2)); for (i in list(ms1,ms2,ms3,ms4)) image(i)

}
\keyword{spatial}%-- one or more ...
