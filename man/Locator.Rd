\name{Locator}
\alias{Locator}
\title{Graphical Input}
\description{
  Reads the position of the graphics cursor when the (first) mouse
  button is pressed, or replays from storage.
}
\usage{
Locator(n, type="n", info=NULL, ...)
}
\arguments{
  \item{n}{the maximum number of points to locate; the value of \code{n}
    is ignored if }
  \item{type}{One of '"n"', '"p"', '"l"' or '"o"'. If '"p"' or '"o"' the
    points are plotted; if '"l"'  or '"o"' they are joined by lines.}
  \item{info}{arbitrary object; tracing information that is useful
    in case the user likes to edit the stored sequence of inputs;
    the value of \code{info} is considered only if \command{\link{useraction}}
    has been called by \code{action="start.register"}
    or \code{action="continue.register"} 
}
  \item{...}{additional graphics parameters}
}
\value{
    A list containing 'x' and 'y' components which are the coordinates
     of the identified points in the user coordinate system, i.e., the
     one specified by 'par("usr")'.
}
\details{
  The behaviour of \command{Locator} depends on the
  the value of \code{action} set by \command{\link{useraction}},
  see there for more information.

  See also \command{\link[graphics]{locator}} for further information on
  the parameters.
}
\seealso{
  \command{\link{getactions}},
  \command{\link{Readline}},
  \command{\link{useraction}},
  \command{\link{userinput}}
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\examples{
## see useraction
}
\keyword{iplot}
\keyword{utilities}

