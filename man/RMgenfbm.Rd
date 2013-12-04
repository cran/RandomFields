\name{RMgenfbm}
\alias{RMgenfbm}

\title{Generalized Fractal Brownian Motion Variogram Model}
\description{
  \command{\link{RMgenfbm}} is an intrinsically stationary isotropic
 variogram model.
 The corresponding centered semi-variogram only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{\gamma(r) = (r^{\alpha}+1)^{\delta}-1}{\gamma(r)=(r^{\alpha}+1)^{\delta}-1}
 where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2} and \eqn{0< \delta <
 1}.\cr
 By now, the model is implemented for dimensions up to 3.\cr
 See also \command{\link{RMfbm}}.

}
\usage{
RMgenfbm(alpha, delta, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value; should be in the interval (0,2].}
 \item{delta}{a numerical value; should be in the interval (0,1).}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 variogram remains unmodified.}
}
\details{
 Here the variogram of \command{\link{RMfbm}} is slightly modified by
 the transformation \eqn{(\gamma+1)^\delta-1} on variograms \eqn{\gamma}.
 The model is implemented up to dimension \eqn{d=3}. 
}
\value{
 \command{\link{RMgenfbm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Stein, M.L. (2002) Fast and exact simulation of fractional
 Brownian surfaces. \emph{J. Comput. Graph. Statist.} \bold{11},
 587--599
 % \item Martin's Toledo-Chapter: Construction of covariance functions
 % and unconditional simulation of random fields, Application to variograms
 }
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMfbm}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
set.seed(0)
model <- RMgenfbm(alpha=1, delta=0.5)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
}