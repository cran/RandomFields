\name{RMconstant}
\alias{RMconstant}
\alias{RMconstant}
\title{Covariance Matrix Constant in Space}
\description{
  \command{\link{RMconstant}} defines a spatially constant covariance function.
}
\usage{
RMconstant(M, var)
}
\arguments{
 \item{M}{a numerical matrix defining the user-defined covariance for a
 random field; the matrix should be positive definite, symmetric and
 its dimension should be equal to the length of observation or
 simulation vector.}
\item{var}{ variance }
 }

\value{
 \command{\link{RMconstant}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}

\me

\seealso{
 \command{\link{RMfixcov}},
\command{\link{RMmodel}}.
}


\keyword{spatial}
\keyword{models}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMconstant(diag(2),var=3)
plot(model)
x <- seq(0,10,length=100)
z <- RFsimulate(model=model,x=x)

\dontshow{FinalizeExample()}}