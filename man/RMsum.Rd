\name{RMsum}
\alias{RMsum}
\title{Plain scalar product}
\description{
 \command{\link{RMsum}} is given by
 \deqn{C(x,y) =  \phi(x) + \phi(y) }{C(x,y) =  \phi(x) + \phi(y).}
 It is a negative definite function although not a variogram.
}
\usage{
RMsum(phi, var, scale, Aniso, proj)
}
\arguments{
  \item{phi}{any function of class \code{\link[=RMmodel-class]{RMmodel}}
  }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
%\details{}
\value{
 \command{\link{RMsum}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\note{
  Do not mix up this model with  \command{\link{RMplus}}.
}
%\references{
%  Wendland, H. \emph{Scattered Data Approximation} (2005) Cambridge: Cambridge
%  University Press.  
%}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RMplus}}
 \command{\link{RMprod}}
}


\keyword{spatial}
\keyword{models}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

\dontshow{\dontrun{
RFfctn(RMsum(RMexp()), 1:10, grid=FALSE)
RFfctn(RMsum(RMexp()), 1:10, 1:10, grid=FALSE)
RFfctn(RMsum(RMexp()), as.matrix(1:10), as.matrix(1:10), grid=FALSE)
}}


\dontshow{FinalizeExample()}
}