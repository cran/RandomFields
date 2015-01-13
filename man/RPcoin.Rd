\name{Coins}
\alias{Coins}
\alias{RPcoins}
\alias{Average}
\alias{RPaverage}
\title{Random coin method}
\description{
 The random coin method (or dilution method) is simulation method for
 stationary Gaussian random fields. It is based on the following procedure:
 For a stationary Poisson point process on \eqn{{\bf R}^d}{R^d}
 consider the random field
 \deqn{Y(y) = \sum_{x\in X} f(y-x)}{Y(y) = \sum_{x\in X} f(y-x)}
 for a function \eqn{f}{f}. The covariance of \eqn{Y}{Y} is
 proportional to the convolution
 \deqn{C(h) = \int f(x)f(x+h) dx }{C(h) = \int f(x)f(x+h) dx}
 If the intensity of the Poisson point process increases, the
 random field \eqn{Y}{Y} is approaches a Gaussian random field
 with covariance function \eqn{C}{C}.
}

\usage{
RPcoins(phi, shape, intensity, method) 

RPaverage(phi, shape, intensity, method) 

}

\arguments{
  \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
    specifies the covariance function of the Poisson process;
    either \code{phi} or \code{shape} must be given.
  }
  \item{shape}{object of class \code{\link[=RMmodel-class]{RMmodel}};
    specifies the function which is attached to the Poisson points;
    note that this is not the covariance function of the simulated
    random field.}
  \item{intensity}{positive number, intensity of the underlying Poisson
    point process.
  }
  \item{method}{integer.
    Default is the value \code{0} which addresses the current standard
    procedure. There might be further methods implemented mainly for
    internal purposes.
  }
}
 

 

%\details{}



\value{
 \command{\link{RPcoins}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
 \itemize{

 \item Lantuejoul, C. (2002)
 \emph{Geostatistical Simulation: Models and Algorithms.}
 Springer.
 }}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
  \command{\link{RFgetMethodNames}}
  \link{RP},
 \command{\link{RPhyperplane}},
 \command{\link{RPspectral}},
 \command{\link{RPtbm}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

\dontshow{\dontrun{
x <- seq(0,25, if (interactive()) 0.02 else 5)
model <- RPcoins(RMspheric())
z <- RFsimulate(model, x, x, spConform=FALSE) # takes 20 sec
Print(sd(as.vector(z)), mean(z))
image(z)


### Gaussian field approximates by Poisson fields
x <- seq(0,10, if (interactive()) 0.02 else 5)
for (intensity in c(1, 10, 100)) {
 z <- RFsimulate(x=x, model=RPcoins(RMspheric(), intensity = intensity))
 plot(z)
}
}}

\dontshow{FinalizeExample()}
}


\keyword{methods}
