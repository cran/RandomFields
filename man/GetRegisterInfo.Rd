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
GetRegisterInfo(register=0, ignore.active = FALSE)
GetModelInfo(register=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{register}{0:9; place where intermediate calculations are stored;
    the numbers are aliases for 10 internal registers, see also
    \command{\link{GaussRF}}}
  \item{ignore.active}{logical. If \code{FALSE} and the register has
    non-active flag (because of an error or it is deleted) then a list
    is returned that contains only the element \code{active=FALSE}.
    Otherwise, the full list is always returned.
  }
}
\details{
  \code{GetRegisterInfo(register, ignore.active=TRUE)} can be
  useful if an error of simulation has been occured.
  
  If \command{\link{RFparameters}}\code{()$Storage=FALSE} then values of
  the internal registers are not kept. Hence \code{GetRegisterInfo} cannot
  provide information.

  \code{GetModelInfo} returns the partial information on the internal
  model structure. It allows for \code{register==-1} returning the
  model structure for the last use of \command{CovarianceFct}
  or \command{Variogram} or similar commands.
  \code{register==-2} is for internal use only.
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
