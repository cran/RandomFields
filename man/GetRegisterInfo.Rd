\name{GetRegisterInfo}
\alias{GetRegisterInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal information}
\description{
  The function returns internal information about the simulation
  of a random field
}
\usage{
GetRegisterInfo(register, ignore.active = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{register}{0:9; place where intermediate calculations are stored;
    the numbers are aliases for 10 internal registers, see also
    \code{\link{GaussRF}}} 
  \item{ignore.active}{logical. If \code{FALSE} and the register has
    non-active flag (because of an error or it is deleted) then a list
    is returned that contains only the element \code{active=FALSE}.
    Otherwise, the full list is always returned.
  }
}
\details{
  \code{GetRegisterInfo(register, ignore.active=TRUE)} can be
  useful if an error of simulation has been occured.
  
  If \code{\link{RFparameters}()$Storage=FALSE} then values of
  the internal registers are not kept. Hence \code{GetRegisterInfo} cannot
  provide information.
}
\value{
  List of internal information is returned.
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}

\seealso{\code{\link{GaussRF}}}
\examples{
str(GetRegisterInfo(0))
try(GaussRF(1:4, grid=TRUE, model="exp", param=c(1,2,3,4)))
str(GetRegisterInfo(0))

try(GaussRF(runif(4), grid=FALSE, model="exp", param=c(1,2,3,4), me="ci"))
str(GetRegisterInfo(0, TRUE))
}
\keyword{spatial}
