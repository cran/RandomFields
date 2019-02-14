\name{RMfixcov}
\alias{RMfixcov}
\title{Fixed Covariance Matrix}
\description{
  \command{\link{RMfixcov}} is a user-defined covariance according to
  the given covariance matrix.

  It extends to the space through a Voronoi tessellation.
}
\usage{
RMfixcov(M, x, y=NULL, z=NULL, T=NULL, grid, var, proj, raw, norm)
}
\arguments{
 \item{M}{a numerical matrix defining the user-defined covariance for a
   random field; the matrix should be positive definite, symmetric and
   its dimension should be equal to the length of observation or
   simulation vector.}
 \item{x,y,z,T,grid}{optional.
   The usual arguments as in \command{\link{RFsimulate}} to define the
   locations where the covariates are given.
 }
 \item{var,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.}
 
 \item{raw}{\argRaw}
  \item{norm}{\argNorm}

 % \item{vdim}{an integer value; defining the response dimension.}
}

\note{
  Starting with version 3.0.64, the former argument \code{element}
  is replaced by the \code{general} option \code{set} in
  \command{\link{RFoptions}}. 
}
\details{
  The covariances passed are implemented for the given locations.
  Within any Voronoi cell (around a given location) the correlation is
  assumed to be one.

 In particular, it is used in \command{\link{RFfit}} to define neighbour or network structure in the data.
}
\value{
 \command{\link{RMfixcov}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Ober, U., Ayroles, J.F., Stone, E.A., Richards, S., Zhu, D., Gibbs, R.A., Stricker, C., Gianola, D., Schlather, M., Mackay, T.F.C., Simianer, H. (2012): \emph{Using Whole Genome Sequence Data to Predict Quantitative Trait Phenotypes in Drosophila melanogaster}. PLoS Genet 8(5): e1002685.
 }
}
\me
\seealso{
 \command{\link{RMcovariate}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RMuser}}.
}


\keyword{spatial}
\keyword{models}
\examples{\dontshow{StartExample(reduce=FALSE)}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


## Example 1 showing that the covariance structure is correctly implemented
n <- 10
C <- matrix(runif(n^2), nc=n)
(C <- C \%*\% t(C))
RFcovmatrix(RMfixcov(C), 1:n)


## Example 2 showing that the covariance structure is interpolated
RFcovmatrix(RMfixcov(C, 1:n), c(2, 2.1, 2.5, 3))


## Example 3 showing the use in a separable space-time model
model <- RMfixcov(C, 1:n, proj="space") * RMexp(s=40, proj="time")
(z <- RFsimulate(model, x = seq(0,12, 0.5), T=1:100))
plot(z)

\dontshow{FinalizeExample()}}