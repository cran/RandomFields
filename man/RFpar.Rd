\name{RFpar}
\alias{RFpar}
\title{Graphical parameters for plots}
\description{
  This function sets globally graphical parameters for plots
  of \code{RMmodels}, simulations and estimations.
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

\me

\seealso{
  \link{plot-method}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFpar(col="red")
plot(RMexp())

\dontshow{FinalizeExample(); RFpar(col=NULL)}}

\keyword{spatial}
\keyword{models}
