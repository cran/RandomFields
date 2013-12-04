\name{RMnatsc}
\alias{RMnatsc}
\title{Natural scale}
\description{
 \command{\link{RMnatsc}} is a stationary isotropic covariance model that
 depends on a stationary isotropic covariance model \eqn{\phi}{phi}.
 The covariance is given by 
 \deqn{ C(h) = \phi(h / s) }{ C(h) = phi(h / s) },
 where the parameter \code{s} is chosen by \command{\link{RMnatsc}}
 such that the practical range or the mathematical range, if finite) is 1.
}
\usage{
RMnatsc(phi, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a stationary isotropic covariance \command{\link{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
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
set.seed(0)
model <- RMnatsc(RMexp())
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1)) 
print(RFcov(model, 1)) 
plot(RFsimulate(model, x=x))
}