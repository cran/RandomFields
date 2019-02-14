\name{RMid}
\alias{RMid}
\title{Identical Model}
\description{
  \command{RMid} is the identical function \eqn{f(x) = x}
  where \eqn{x} is a vector of coordinates and \eqn{f(x)} is a model value.
 }
 \usage{
RMid()%var, scale, Aniso, proj)
}

%\arguments{
% \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
%   \command{\link{RMmodel}}. If not passed, the above
%   covariance function remains unmodified.}
%}

 
\value{
 \command{\link{RMid}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}


\me

\seealso{
  \command{\link{RMmodel}}.
  \command{RMid},
  \command{RMtrafo},
  \command{RMprod}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## C(x,y) = < x, y >
RFcov(RMprod(RMid()), as.matrix(1:10), as.matrix(1:10), grid=FALSE)
\dontshow{FinalizeExample()}}


\keyword{spatial}
\keyword{models}

