\name{RMbcw}
\alias{RMbcw}

\title{Model bridging stationary and intrinsically stationary processes}
\description{
  \command{\link{RMbcw}} is a variogram model
  that bridges between some intrinsically stationary isotropic processes
  and some stationary ones. It reunifies the
  \link{RMgenfbm}, \link{RMdewijsian} and \link{RMgencauchy}.

 The corresponding centered semi-variogram only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{\gamma(r) =
   \frac{(r^{\alpha}+1)^{\beta/alpha}-1}{2^{\beta/alpha} -1}}{
   \gamma(r)=[(r^{\alpha}+1)^{\beta/alpha}-1] / (2^{\beta/alpha}-1)}
 where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2} and \eqn{beta \le 2}{\beta <= 2}.\cr
}
\usage{
RMbcw(alpha, beta, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value; should be in the interval (0,2].}
 \item{beta}{a numerical value; should be in the interval (-infty,2].}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 variogram remains unmodified.}
}
\details{
  For \eqn{betaa >0}, \eqn{beta<0}, \eqn{beta=0}
  we have the generalised fractal Brownian motion \command{\link{RMgenfbm}},
  the generalised Cauchy model \command{\link{RMgencauchy}},
  and the de Wisjian model \command{\link{RMdewijsian}}, respectively.

  Hence its two arguments \code{alpha} and \code{beta}
  allow for modelling the smoothness and a wide range of tail behaviour,
  respectively. 
  
}
\value{
 \command{\link{RMbcw}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Schilling, R.L., Song, R., Vondracek, Z. (2010) Berstein
 functions. Theory and Applications. De Gruyter
 
 \item Schlather, M (2014) A parameteric variogram model bridging
 between stationary 
and intrinsically stationary processes. \emph{Submitted} \bold{}, --
 % \item Martin's Toledo-Chapter: Construction of covariance functions
 % and unconditional simulation of random fields, Application to variograms
 }
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMgenfbm}},
 \command{\link{RMgencauchy}},
 \command{\link{RMdewijsian}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMbcw(alpha=1, beta=0.5)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}