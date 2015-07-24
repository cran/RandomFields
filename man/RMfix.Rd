\name{RMfixcov}
\alias{RMfixcov}
\title{Fixed Covariance Matrix}
\description{
  \command{\link{RMfixcov}} is a user-defined covariance according to
  the given covariance matrix.

  It extends to the space through a Voronoi tesselation.
}
\usage{
RMfixcov(M, x, y=NULL, z=NULL, T=NULL, grid, var, scale, Aniso, proj, 
raw, norm)
}
\arguments{
 \item{scale, Aniso, proj,var}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.}
 
 \item{M}{a numerical matrix defining the user-defined covariance for a
   random field; The matrix should be positive definite, symmetric and
   its dimension should be equal to the length of observation or
   simulation vector.}
 \item{x,y,z,T,grid}{optional.
   The usual arguments as in \command{\link{RFsimulate}} to define the
   locations where the covariates are given
 }
 \item{raw}{
    logical. If \code{FALSE} then the data are interpolated. This
    approach is always save, but might be slow.

    If \code{TRUE} then the data may be accessed when covariance
    matrices are calculated. No rescaling or anisotropy definition
    is allowed in combination with the model. The use is dangerous,
    but fast.

    Default: FALSE (outside mixed models)
  }
  \item{norm}{optional model that gives the norm between locations}

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
 \command{\link{RMfixcov}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Ober, U., Ayroles, J.F., Stone, E.A., Richards, S., Zhu, D., Gibbs, R.A., Stricker, C., Gianola, D., Schlather, M., Mackay, T.F.C., Simianer, H. (2012): \emph{Using Whole Genome Sequence Data to Predict Quantitative Trait Phenotypes in Drosophila melanogaster}. PLoS Genet 8(5): e1002685.
 }
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMcovariate}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RMuser}}
}


\keyword{spatial}
\keyword{models}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


## Example 1 showing that the covariance structure is correctly implemented
n <- 10
z <- matrix(runif(n^2), nc=n)
(z <- z \%*\% t(z))
RFcovmatrix(RMfixcov(z), 1:n)


## Example 2 showing that the covariance structure is interpolated
RFcovmatrix(RMfixcov(z, 1:n), c(2, 2.1, 2.5, 3))


## Example 3 showing the use in a separable space-time model
model <- RMfixcov(z, 1:n, proj="space") * RMexp(s=40, proj="time")
(z <- RFsimulate(model, x = seq(0,12, 0.5), T=1:100))
plot(z)


\dontshow{FinalizeExample()}
}