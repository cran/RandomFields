\name{RPchi2}
\alias{RPchi2}

\title{Simulation of Chi2 Random Fields}
\description{
  \command{RPchi2} defines a chi2 field.
}

\usage{

RPchi2(phi, boxcox, f)

}

\arguments{
 \item{phi}{the \command{\link{RMmodel}}. If a model for the
 distribution is not specified, \command{\link{RPgauss}} is used as
 default and a covariance model is expected.}
\item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  See \command{\link{RFboxcox}} for details.
 }
 \item{f}{integer. Degree of freedom.}
}

\value{
 The function returns an object of class \code{\link[=RMmodel]{RMmodel}}.
}

\me

\seealso{
  \command{\link{Auxiliary RMmodels}},
  \command{\link{RP}},
  \command{\link{RPgauss}}.
}


\keyword{spatial}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RPchi2(RMexp(), f=2)
x <- seq(0, 10, 0.1)
z <- RFsimulate(model=model, x, x, n=4)
plot(z)

\dontshow{FinalizeExample()}}
