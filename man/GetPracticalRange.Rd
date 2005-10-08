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
GetPracticalRange(model, kappas=NULL, dim=1)
}
\arguments{
  \item{model}{string; covariance model, see \command{\link{CovarianceFct}}, or
    type \command{\link{PrintModelList}}\code{()} to get all options.}
  \item{kappas}{additional parameters of a parametrised model, see
    \command{\link{CovarianceFct}}.}
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
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}

\seealso{  \command{\link{CovarianceFct}},
  \code{\link{RandomFields}},
  \command{\link{RFparameters}},
}

\examples{
GetPracticalRange("exponential")
GetPracticalRange("whittlematern",kappa=2)
GetPracticalRange("gengneiting",kappa=c(a=3,b=5))
}
\keyword{spatial}
