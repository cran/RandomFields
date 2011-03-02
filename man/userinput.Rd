\name{userinput}
\alias{userinput}
\title{General input function}
\description{
  The function provides a general functionality to read and
  reply user input by \command{\link[graphics]{locator}} and
  \command{\link[base]{readline}} 
}
\usage{
userinput(fct, info=NULL, prompt, n,  type="n", pch=par()$pch,
          cex=par()$cex, ...)
}
\arguments{
  \item{fct}{either \sQuote{locator} or \sQuote{readline}}
  \item{info}{arbitrary object; tracing information that is useful
    in case the user likes to edit the stored sequence of inputs}
  \item{prompt}{prompt parameter of \command{\link[base]{readline}}}
  \item{n}{parameter \code{n} of \command{\link[graphics]{locator}}}
  \item{type}{parameter \code{type} of \command{\link[graphics]{locator}}}
  \item{pch}{character type used by \command{\link[graphics]{locator}}}
  \item{cex}{character size used by \command{\link[graphics]{locator}}}
  \item{...}{additional graphics parameters used by
    \command{\link[graphics]{locator}}}
}
\value{
  If \code{fct="readline"} then a character vector of length one is
  returned. 

  If \code{fct="locator"} then the function returns a
    list containing 'x' and 'y' components which are the coordinates
     of the identified points in the user coordinate system, i.e., the
     one specified by 'par("usr")'.
  
}
\details{
  The behaviour of \command{userinput} depends on the
  the value of \code{action} set by \command{\link{useraction}},
  see there for more information.  
}
\seealso{
  \command{\link{getactions}},
  \command{\link{Locator}},
  \command{\link{Readline}},
  \command{\link{useraction}}
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\keyword{iplot}
\keyword{utilities}
\keyword{internal}
