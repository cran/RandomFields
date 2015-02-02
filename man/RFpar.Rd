\name{RFpar}
\alias{RFpar}
\title{Graphical parameters for plots}
\description{
  This functions sets globally graphical parameters for plots
  of \code{RMmodels}, simulations and estimations
}
\usage{
RFpar(...)
}
\arguments{
  \item{...}{see \command{\link[graphics]{par}}
  }
}
\value{
  If \command{RFpar} is called without arguments, the current list is
  returned.

  If \command{RFpar} is called with \code{NULL} only, the current list is
  deleted.

  Otherwise the arguments are stored for global use in RandomFields.
 
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\seealso{
  \link{plot-method}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFpar(col="red")
plot(RMexp())

\dontshow{FinalizeExample()}
}

\keyword{spatial}
\keyword{models}
