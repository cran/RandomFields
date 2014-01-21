\name{RMdewijsian}
\alias{RMdewijsian}
\title{Generalized DeWijsian Variogram Model}
\description{
  \command{\link{RMdewijsian}} is an intrinsically stationary isotropic variogram model.
 The corresponding centered semi-variogram only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{\gamma(r) = \log(r^{\alpha}+1)}{\gamma(r)=log(r^{\alpha}+1)}
 where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2}.
}
\usage{
RMdewijsian(alpha, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value; should be in the interval (0,2]
 to provide a valid variogram for a random field of any dimension.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 variogram remains unmodified.}
}
\details{
 The parameter \eqn{\alpha}{\alpha} must satisfy \eqn{\alpha \in (0,2]}{0
 < \alpha \le 2}.
 
 Originally, the logarithmic model \eqn{\gamma(r) = \log(r)} was named
 after de Wijs and reflects a principle of similarity (cf. Chiles,
 J.-P. and Delfiner, P. (1999), p. 90).
 But note that \eqn{\gamma(r) = \log(r)} is not a valid variogram
 (\eqn{\gamma(0)} does not vanish)
 whereas \eqn{\gamma(r) = \log(r^{\alpha}+1)} is valid
 (cf. Wackernagel, H. (2003), p. 336).
}

\value{
 \command{\link{RMdewijsian}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.

 \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 3nd edition.
 % \item Martin's Toledo-Chapter: Construction of covariance functions
 % and unconditional simulation of random fields, Example 7
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0)
model <- RMdewijsian(alpha=1)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{RFoptions(seed=NA)}
}