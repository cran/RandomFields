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
GetPracticalRange(model, kappas=NULL)
}
\arguments{
  \item{model}{string; covariance model, see \code{\link{CovarianceFct}}, or
    type \code{\link{PrintModelList}()} to get all options.}
  \item{kappas}{additional parameters of a parametrised model, see
    \code{\link{CovarianceFct}}.}
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
\author{Martin Schlather, \email{martin.schlather@cu.lu}
  \url{http://www.cu.lu/~schlathe}}

\seealso{  \code{\link{CovarianceFct}},
  \code{\link{RandomFields}},
  \code{\link{RFparameters}},
}

\examples{
GetPracticalRange("exponential")
GetPracticalRange("whittlematern",kappa=2)
GetPracticalRange("gengneiting",kappa=c(a=3,b=5))
}
\keyword{spatial}
