\name{parameter.range}
\alias{parameter.range}

\title{Range of model specific parameters}
\description{
  \code{parameter.range} returns the range of the parameters specific to
  a parametrised variogram model or a covariance function
}
\usage{
parameter.range(model, dim)
}
\arguments{
  \item{model}{the variogram model or covariance function, see \code{CovarianceFct}}
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

\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}

\seealso{\command{\link{CovarianceFct}}}

\examples{
  parameter.range("exponential", 1) # NULL
  str(parameter.range("power", 1))
  str(parameter.range("power", 2))
  str(parameter.range("gengneiting", 2))
  str(parameter.range("nsst",1))
  str(parameter.range("nsst2",1))
  str(parameter.range("hyper",1))
}

\keyword{spatial}



