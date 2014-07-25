\name{RMgenfbm}
\alias{RMgenfbm}

\title{Generalized Fractal Brownian Motion Variogram Model}
\description{
  \command{\link{RMgenfbm}} is an intrinsically stationary isotropic
 variogram model.
 The corresponding centered semi-variogram only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{\gamma(r) = (r^{\alpha}+1)^{\beta/\alpha}-1}{\gamma(r)=(r^{\alpha}+1)^{\beta/\alpha}-1}
 where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2} and \eqn{\beta \in (0,2]}.\cr
 See also \command{\link{RMfbm}}.
}
\usage{
RMgenfbm(alpha, beta, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value; should be in the interval (0,2].}
 \item{beta}{a numerical value; should be in the interval (0,2].}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 variogram remains unmodified.}
}
\details{
 Here the variogram of \command{\link{RMfbm}} is modified by
 the transformation \eqn{(\gamma+1)^{\delta/-1}} on variograms
 \eqn{\gamma}
 for \eqn{delta \in (0,1]}. This original modification allows for
 further generalization, cf. \command{\link{RMbcw}}.
}
\value{
 \command{\link{RMgenfbm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Gneiting, T. (2002) Nonseparable, stationary covariance
 functions for space-time data, \emph{JASA} \bold{97}, 590-600.
 
 \item
 Schlather, M. (2010)
 On some covariance models based on normal scale mixtures.
 \emph{Bernoulli}, \bold{16}, 780-797.

 
 % \item Martin's Toledo-Chapter: Construction of covariance functions
 % and unconditional simulation of random fields, Application to variograms
 }
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
  \command{\link{RMbcw}}
 \command{\link{RMfbm}},
 \command{\link{RMmodel}},
 \command{\link{RMflatpower}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMgenfbm(alpha=1, beta=0.5)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}