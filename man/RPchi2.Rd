\name{RPchi2}
\alias{RPchi2}

\title{Simulation of a Chi2 Random Fields}
\description{
  \command{RPchi2} defines a chi2 fields.
}

\usage{

RPchi2(phi, f)

}

\arguments{
 \item{phi}{the \command{\link{RMmodel}}. If a model for the
 distribution is not specified, \command{\link{RPgauss}} is used as
 default and a covariance model is expected.}
 \item{f}{integer. Degree of freedom.}
}

\value{
 The function returns an object of class \code{\link[=RMmodel]{RMmodel}}
}

 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}


\seealso{
  \command{\link{Auxiliary RMmodels}},
  \command{\link{RP}},
  \command{\link{RPgauss}}
}


\keyword{spatial}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RPchi2(RMexp(), f=2)
x <- seq(0, 10, if (interactive()) 0.01 else 2)
z <- RFsimulate(model=model, x, x, n=4)
plot(z)

\dontshow{FinalizeExample()}
}
