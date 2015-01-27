\name{RMfixcov}
\alias{RMfixcov}
\alias{RMfixcov}
\title{Fixed Covariance Matrix}
\description{
  \command{\link{RMfixcov}} is a user-defined covariance according to
  the given covariance matrix, multiplied with \code{var} 
}
\usage{
RMfixcov(M, vdim, element, var)
}
\arguments{
 \item{M}{a numerical matrix defining the user-defined covariance for a
 random field; The matrix should be positive definite, symmetric and
 its dimension should be equal to the length of observation or
 simulation vector.}
 \item{vdim}{an integer value; defining the response dimension (for details see \command{\link{RFgetModelNames}}).}
 \item{element}{integer. If a list of several matrices is given then
 \code{element} chooses the one that should currently be
 used. Default: 0}
 \item{var}{variance, i.e. factor multiplied to the data, which can be
    estimated through ML
  }
}
\details{
 Note, this is not a covariance model for arbitrary locations, the
 covariance is fixed and user-defined by a given matrix corresponding
 to the locations. 
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
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}
\examples{

## to do

\dontshow{FinalizeExample()}
}