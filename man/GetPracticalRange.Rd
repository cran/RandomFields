\name{GetPracticalRange}
\alias{GetPracticalRange}
\alias{PracticalRange}
\title{Determination of the practical range}
\description{
  the function returns the practical range of a covariance model,
  i.e. the distance for which the model with standard parameters
  (\code{variance=1}, \code{nugget=0}, \code{scale=1}) values 0.05
}
\usage{
GetPracticalRange(model, param, dim=1)
}
\arguments{
  \item{model}{covariance model in the new list not;
    the model must be primitive, i.e. may not call other models.
    See \command{\link{Covariance}} and \command{\link{CovarianceFct}}
    for details.}
  \item{param}{optional vector of parameters for the model,
  see \command{\link{CovarianceFct}}}
  \item{dim}{space-time dimension; except for some few anisotropic
    models, this parameter does not have an effect}
}
\details{
  In case the practical range has to be determined numerically, only a
  rough approximation is provided.

  The function can only be applied to isotropic models
}
\value{
  real number (the practical range)
}
\references{
  Goovaerts, P. (1997) \emph{Geostatistics for Natural Resources
    Evaluation.} New York: Oxford University Press.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}

\seealso{  \command{\link{CovarianceFct}},
  \code{\link{RandomFields}},
  \command{\link{RFparameters}},
}

\examples{
GetPracticalRange(list("exponential"))
GetPracticalRange(list("whittle", nu=2))
try(GetPracticalRange(list("bessel", nu=-0.5))) # leads to an error
}
\keyword{spatial}
