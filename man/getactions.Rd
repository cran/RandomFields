\name{getactions}
\alias{getactions}
\alias{getactionlist}
\title{Get input behaviour}
\description{
  The functions return values stored by \code{\link{useraction}},
  \code{\link{Locator}} and \code{\link{Readline}}
}
\usage{
  getactions()
  getactionlist()
}
\value{
  \command{getactionlist} returns a list of stored values
  created by \code{\link{Locator}} and \code{\link{Readline}}
  
  \command{getactions} returns a list of all parameters that influence
  the behaviour of \code{\link{Locator}} and \code{\link{Readline}}.
}
\seealso{
  \code{\link{Locator}},
  \code{\link{Readline}},
  \code{\link{useraction}},
  %\code{\link{userinput}},
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\examples{
## see useraction
}
\keyword{utilities}

