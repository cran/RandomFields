\name{RMfixed}
\alias{RMfixed}
\title{Fixed Effect Model}

%\usage{
%}

\description{Expressions of the form \code{X@RMfixed(beta)} can be used
 within a formula of the type
 \deqn{response ~ fixed effects + random effects + error term}
 that specifies the Linear Mixed Model.

 Important Remark: \code{RMfixed} is NOT a function although
 the parantheses notation is used to specify the vector of coefficients.

 The matrix \eqn{X} is the design matrix  and \eqn{\beta} is a vector of coefficients

 Note that a fixed effect of the form \eqn{X} is interpreted as
 \code{X@RMfixed(beta=NA)} by default (and \eqn{\beta} is estimated
 provided that the formula is used in \command{\link{RFfit}}).
 Note that the \code{1} in an expression \code{1@RMfixed(beta)} is
 interpreted as the identity matrix. 
}

%\arguments{
%}

%\value{
%}

%\references{Chiles, J. P., Delfiner, P. (1999) \emph{Geostatistics:
% Modelling Spatial Uncertainty.} New York: John Wiley & Sons.
%}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFformula}},
 \command{\link{RFsimulate}},
}

\examples{
 ## For examples see the help page of 'RFformula' ##
\dontshow{RFoptions(seed=NA)}
}
		 
\keyword{models}


