\name{RMnatsc}
\alias{RMnatsc}
\title{Natural scale}
\description{
 \command{\link{RMnatsc}} is a stationary isotropic covariance model that
 depends on a stationary isotropic covariance model \eqn{\phi}{phi}.
 The covariance is given by 
 \deqn{ C(h) = \phi(h / s) }{ C(h) = phi(h / s) },
 where the argument \code{s} is chosen by \command{\link{RMnatsc}}
 such that the practical range or the mathematical range, if finite) is 1.
}
\usage{
RMnatsc(phi, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a stationary isotropic covariance \command{\link{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 For internal use only.
}
\value{
 \command{\link{RMnatsc}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
%\references{}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMnatsc(RMexp())
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(RMexp(), model=model)
RFcov(model, 1)
\dontshow{FinalizeExample()}
}