\name{CheckAndComplete}
\alias{CheckAndComplete}
\title{Validity of a covariance model}
\description{
  \command{CheckAndComplete} checks the validity of a covariance model;
  further a matrix of internal parameters is returned.
}
\usage{
CheckAndComplete(model, param, dim)
}
\arguments{
  \item{model}{definition of a covariance function,
    see \code{\link{CovarianceFct}} for examples}
  \item{param}{see model}
  \item{dim}{dimension for the space}
}

\value{
  A list of two parameters is returned:

  \code{error}
    a value different from zero signifies an error;
    to get a readable error information put
    \code{RFparameters()$PrintLevel}
    to a higher level and rerun.

    \code{param}
    (currently for internal purposes only).
    matrix of the internal parameter list; columns correspond to
    the simple models and can be extracted using
    \command{\link{paramextract}}.
  }
}
\author{Martin Schlather, \email{schlath@hsu-hh.de},
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}

\seealso{\command{\link{GaussRF}},  \command{\link{paramextract}},
  \code{\link{RandomFields}}, \code{\link{RFMethods}}.}

\examples{
 ## for an example see \link{RFMethods}
}
\keyword{spatial}





