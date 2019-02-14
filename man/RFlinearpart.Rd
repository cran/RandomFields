\name{RFlinearpart}
\alias{RFlinearpart}
\title{Linear part of \command{\link{RMmodel}}}
\description{
 \command{\link{RFlinearpart}} returns the linear part of a model
}
\usage{
RFlinearpart(model, x, y = NULL, z = NULL, T = NULL, grid=NULL,
                data, params, distances, dim, set=0, ...)

}
\arguments{
 \item{model,params}{\argModel}
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{distances,dim}{\argDistances}
 \item{data}{\argData}
 \item{set}{integer. See section Value for details.}
 \item{...}{\argDots}
}


\note{
  In the linear part of the model specification the parameters
  that are NA must be the first model part. I.e.
  \code{NA * sin(R.p(new="isotropic")) + NA + R.p(new="isotropic")}
  is OK, but not
  \code{sin(R.p(new="isotropic")) * NA + NA + R.p(new="isotropic")}

}

\value{  
  \command{\link{RFlinearpart}} returns a list
  of three components, \code{Y}, \code{X}, \code{vdim} returning
  the deterministic trend, the design matrix, and the multivariability,
  respectively.
  If \code{set} is positive, \code{Y} and \code{X} contain
  the values for the \code{set}-th set of coordinates.
  Else, \code{Y} and \code{X} are both lists containing
  the values for all the sets.
}

\me
\seealso{	
\link{Bayesian},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFlikelihood}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


x <- seq(0, pi, len=10)
trend <- 2 * sin(R.p(new="isotropic")) + 3
model <- RMexp(var=2, scale=1) + trend
print(RFlinearpart(model, x=x))  ## only a deterministic part

trend <- NA * sin(R.p(new="isotropic")) + NA + R.p(new="isotropic") / pi
model <- RMexp(var=NA, scale=NA) + trend
print(RFlinearpart(model, x=x))

\dontshow{FinalizeExample()}}

\keyword{spatial}






