\name{GetRegisterInfo}
\alias{GetRegisterInfo}
\alias{GetModelInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal information}
\description{
  The function returns internal information about the simulation
  of a random field
}
\usage{
GetRegisterInfo(register=0, ignore.active = FALSE, max.elements=10^6)
GetModelInfo(register=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{register}{-2, -1, 0:9;
    place where intermediate calculations are stored;
    the numbers \code{0:9} are aliases for 10 internal registers, see also
    \command{\link{GaussRF}}}
  \item{ignore.active}{logical. If \code{FALSE} and the register has
    non-active flag (because of an error or it is deleted) then a list
    is returned that contains only the element \code{active=FALSE}.
    Otherwise, the full list is always returned.
  }
  \item{max.elements}{integer; since \command{GetRegisterInfo}
    might be a dump of simulation that needs a huge amount of
    memory and since all entries are copied, the maximal amount
    of available memory might be easily exceeded.  Therefore,
    only the size of the critical parts are returned and not the
    vector or matrix itself, if the number of elements exceeds
    \code{max.elements}. 
  }
}
\details{
  \code{GetRegisterInfo(register, ignore.active=TRUE)} is
  useful for debugging and specialists' need to control
  the algorithm, see the examples in \code{\link{RFparameters}}
  and \code{\link{GaussRF}}.
  
  If \command{\link{RFparameters}}\code{()$Storage=FALSE} then values of
  the internal registers are not kept if \command{\link{GaussRF}} or
  \command{\link{DoSimulateRF}} has been called.
  Hence \code{GetRegisterInfo} cannot
  provide any information.

  \code{GetModelInfo} returns the partial information on the internal
  model structure. It allows for \code{register==-1} returning the
  model structure for the last use of \command{\link{CovarianceFct}}
  or \command{\link{Variogram}} or similar commands.
  \code{register==-2} is for internal use only.
}
\note{
  Put \code{Storing=TRUE}, see \command{\link{RFparameters}}
  if you like to have more
  internal information in case of an expected failure of an initialisation of
  a random field simulation.
}
\value{
  List of internal information is returned.
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}

\seealso{\command{\link{GaussRF}}}
\examples{
# lets regard at a non-initialised register first:
DeleteAllRegisters()
str(GetRegisterInfo(0))

# now, a register that has successfully initialised
try(GaussRF(1:4, grid=TRUE, model="exp", param=c(1,2,3,4), Storing=TRUE))
str(GetRegisterInfo(0))

# finally, a register that could not be successfully initialised
try(GaussRF(runif(4), grid=FALSE, model="exp", param=c(1,2,3,4),
            me="ci", Storing=TRUE))
str(GetRegisterInfo(0, TRUE))
}
\keyword{spatial}
