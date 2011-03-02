\name{parameter.range}
\alias{parameter.range}

\title{Range of model specific parameters}
\description{
  \code{parameter.range} returns the range of the parameters specific to
  a parametrised variogram model or a covariance function

  --- NOT PROGRAMMED YET --
}
\usage{
parameter.range(model, param, dim=1)
}
\arguments{
  \item{model}{the variogram model or covariance function, see \code{CovarianceFct}}
  \item{param}{see \code{CovarianceFct}}
  \item{dim}{the dimension of the random field}
}
\value{
  If the model is not valid or not specified for the given dimension,
  \code{NaN} is returned.
  
  If the model does not have additional parameters, \code{NULL} is
  returned. Otherwise a list is returned.
  
  If the parameter space is simple, a list of the following elements
  is returned
  \item{theoretical}{The theoretical ranges of the parameters. Note that
  open and closed intervals are not distinguished here.}
  \item{practical}{The ranges of the parameters usually not exceeded in
    practical applications.}
  Both elements are matrices with two rows and the number of columns
  being the number of parameters.
  \cr\cr
  If the parameter space is complex, the parameter space is devided into
  rectangular areas and \code{theoretical} and \code{practical} are
  lists where each element is a matrix as above.

}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}

\seealso{\command{\link{CovarianceFct}}}

\examples{
  parameter.range("exponential", param=NA, dim=1) # NULL
  try(parameter.range("power", param=NA, dim=1))
  try(parameter.range("power", param=c(NA, NA, NA, NA, 2), dim=2))
  Print(parameter.range("gengneiting", param=NA, dim=2))
  Print(parameter.range("nsst", param=NA,  dim=1))
  Print(parameter.range("hyper",param, dim=1))
}

\keyword{spatial}



