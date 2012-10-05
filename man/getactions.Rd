\name{getactions}
\alias{getactions}
\alias{getactionlist}
\title{Get input behaviour}
\description{
  The functions return values stored by \command{\link{useraction}},
  \command{\link{Locator}} and \command{\link{Readline}}
}
\usage{
  getactions()
  getactionlist()
}
\value{
  \command{getactionlist} returns a list of stored values
  created by \command{\link{Locator}} and \command{\link{Readline}}
  
  \command{getactions} returns a list of all parameters that influence
  the behaviour of \command{\link{Locator}} and \command{\link{Readline}}.
}
\seealso{
  \command{\link{Locator}},
  \command{\link{Readline}},
  \command{\link{useraction}},
  %\command{\link{userinput}},
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de}}
\examples{
## see useraction
}
\keyword{utilities}

