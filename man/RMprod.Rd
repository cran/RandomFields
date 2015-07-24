\name{RMprod}
\alias{RMprod}
\title{Plain scalar product}
\description{
 \command{\link{RMprod}} is a non-stationary covariance model given by
 \deqn{C(x,y) = \langle \phi(x), \phi(y)\rangle }{C(x,y) = \langle \phi(x), \phi(y)\rangle.}
}
\usage{
RMprod(phi, var, scale, Aniso, proj)
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
 \command{\link{RMprod}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\details{
  This model defines 
}

\note{
  Do not mix up this model with  \command{\link{RMmult}}.

  See also \command{RMS} for a simple, alternative method to set
  an arbitrary, i.e. location dependent, univariate variance.  
 }
\references{
  Wendland, H. \emph{Scattered Data Approximation} (2005) Cambridge: Cambridge
  University Press.  
 }

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{
  \command{\link{RMsum}},
  \command{\link{RMmodel}},
  \command{\link{RMmult}}
}


\keyword{spatial}
\keyword{models}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFcov(RMprod(RMexp()), as.matrix(1:10), as.matrix(1:10), grid=FALSE)

\dontshow{FinalizeExample()}
}