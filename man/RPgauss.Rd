
\name{RPgauss}
\alias{RPgauss}

\title{Simulation of Gaussian Random Fields}
\description{
 This function is used to specify a Gaussian random field that
 is to be estimated or simulated.
 Return an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}

\usage{

RPgauss(phi, stationary_only)


}

\arguments{
 \item{phi}{the \command{\link{RMmodel}}.}
% \item{loggauss}{logical. If \code{TRUE} then a log-Gaussian random
% field is returned. Default is \code{FALSE}.
% }
 \item{stationary_only}{Logical or NA. Used for the automatic
 choice of methods. 
 \itemize{
 \item \code{TRUE}: the simulation of non-stationary random
 fields is refused. In particular, the intrinsic
 embedding method is excluded and
 the simulation of Brownian motion is rejected.
 \item \code{FALSE}: intrinsic embedding is always allowed,
 actually it's the first one considered in the automatic
 selection algorithm.
 \item \code{NA}: the simulation of the Brownian motion allowed,
 but intrinsic embedding is not used for stationary random fields.
 }
 Default: \code{NA}
 }

}

\value{
 The function returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\note{
  In most cases, \code{RPgauss} need not be given explicitely as
  Gaussian random fields are assumed as default.

  \command{\link{RPgauss}} may not find the fastest
 method neither the most precise one. It just finds any method
 among the available methods. (However it guesses what is a good
 choice.)
 See \command{\link{RFgetMethodNames}} for further information.
 Note that some of the methods do not work for all covariance 
 or variogram models, see \code{\link{RFgetModelNames}(intern=FALSE)}

By default, all Gaussian random fields have zero mean. 
 Simulating with trend can be done by including \command{\link{RMtrend}}
 in the model. 
 
\command{\link{RPgauss}} allows to simulate different classes of random fields,
 controlled by the wrapping model:

 If the submodel is a pure covariance or
 variogram model, i.e. of class \code{\link[=RMmodel-class]{RMmodel}}, a
 corresponding centered Gaussian field is simulated.
 Not only stationary fields but also non-stationary and anisotropic
 models can be used, e.g. zonal anisotropy, geometrical
 anisotropy, separable models, non-separable space-time smodels,
 multiplicative or nested models;
 see \command{\link{RMmodel}} for a list of all available models. 

 
}

 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMexp()
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x, seed=0))
plot(RFsimulate(RPgauss(model), x=x, seed=0), col=2) ## the same
\dontshow{FinalizeExample()}
}


\seealso{
  \link{RP},
  \link{Gaussian},
 \command{\link{RMmodel}},
 \command{\link{RFoptions}},
 \command{\link{RPbrownresnick}},
 \command{\link{RPchi2}},
 \command{\link{RPschlather}}.

 Do not mix up with \command{\link{RMgauss}} or \command{\link{RRgauss}}.

}




\keyword{spatial}


