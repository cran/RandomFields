\name{Readline}
\alias{Readline}
\title{Read a Line}
\description{
  \command{Readline} reads a line from the terminal or from storage
}
\usage{
   Readline(prompt="", info=NULL)
}
\arguments{
  \item{prompt}{the string printed when prompting the user for input. Should
    usually end with a space \code{" "}}
  \item{info}{arbitrary object; tracing information that is useful
    in case the user likes to edit the stored sequence of inputs;
    the value of \code{info} is considered only if \command{\link{useraction}}
    has been called by \code{action="start.register"}
    or \code{action="continue.register"} 
   }
}
\value{
 A character vector of length one.
}
\details{
  The behaviour of \command{Readline} depends on the
  the value of \code{action} set by \command{\link{useraction}},
  see there for more information.
  
  The prompt string will be truncated to a maximum allowed length,
  normally 256 chars (but can be changed in the source code), since the
  function is based on \command{\link[base]{readline}}.
}
\seealso{
  \command{\link{getactions}},
  \command{\link{Locator}},
  \command{\link[base]{readline}}
  \command{\link{useraction}},
  \command{\link{userinput}},
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\examples{
## see useraction
}
\keyword{utilities}

