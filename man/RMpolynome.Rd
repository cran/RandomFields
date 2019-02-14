\name{RMpolynome}
\alias{RMpolynome}
\title{Creating polynomial models}
\description{
  Polynomials, mainly used in trend models, can be created easily with
  this function.
}

\usage{
RMpolynome(degree, dim, value=NA, coordnames = c("x", "y", "z", "T"),
           proj=1:4)
}

\arguments{
  \item{degree}{degree of the polynome
  }
  \item{dim}{number of variables in the polynome
  }
  \item{value}{ values of the coefficients. See Details.
  }
  \item{coordnames}{the names of the variables
  }
  \item{proj}{the projection to certain dimensions
  }
}

\details{
  If the length of \code{value} is smaller than the number of monomials,
  the remaining terms are filled with \code{NA}s. If the length is
  larger, the vector is cut.
}

\value{
 \command{\link{RMpolynome}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}


\me
\seealso{
  \command{\link{RMtrend}},
 \command{\link{RFfit}}.
}

\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
 ## For examples see the help page of 'RFformula' ##

RMpolynome(1, 1)
RMpolynome(1, 2)
RMpolynome(2, 1)
RMpolynome(2, 2)
RMpolynome(3, 3)

\dontshow{FinalizeExample()}}