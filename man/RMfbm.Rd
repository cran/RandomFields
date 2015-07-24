\name{RMfbm}
\alias{RMfbm}
\title{Variogram Model of Fractal Brownian Motion}
\description{
 \command{\link{RMfbm}} is an intrinsically stationary isotropic variogram
 model. The corresponding centered semi-variogram only depends on the
 distance \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{\gamma(r) = r^\alpha}{\gamma(r) = r^\alpha}
 where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2}.\cr
 By now, the model is implemented for dimensions up to 3.\cr
 For a generalized model see also \command{\link{RMgenfbm}}.
}
\usage{
RMfbm(alpha, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{numeric in \eqn{(0,2]}; refers to the fractal dimension of the
 process}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 variogram remains unmodified.}
}
\details{
The variogram is unbounded and belongs to a non-stationary process with
stationary increments. For \eqn{\alpha=1}{\alpha=1} and \code{scale=2}
we get a variogram corresponding to a standard Brownian Motion.

For \eqn{\alpha \in (0,2)}{0 < \alpha < 2} the quantity \eqn{H =
 \frac{\alpha} 2}{H=\alpha/2} is called Hurst index and determines
the fractal dimension \eqn{D} of the corresponding Gaussian sample paths
\deqn{D = d + 1 - H}
where \eqn{d}{d} is the dimension of the random field (see Chiles and
 Delfiner, 1999, p. 89).

}
\value{
 \command{\link{RMfbm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Chiles, J.-P. and P. Delfiner (1999) \emph{Geostatistics. Modeling
 Spatial Uncertainty.} New York, Chichester: John Wiley & Sons.

 \item Stein, M.L. (2002) Fast and exact simulation of fractional
 Brownian surfaces. \emph{J. Comput. Graph. Statist.} \bold{11}, 587--599
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMgenfbm}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
model <- RMfbm(alpha=1)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}