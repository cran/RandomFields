\name{RMdewijsian}
\alias{RMdewijsian}
\title{Modified DeWijsian Variogram Model}
\description{
  The modified \command{RMdewijsian} model
  is an intrinsically stationary isotropic variogram model.
 The corresponding centered semi-variogram only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{\gamma(r) = \log(r^{\alpha}+1)}{\gamma(r)=log(r^{\alpha}+1)}
 where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2}.
}
\usage{
RMdewijsian(alpha, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value; in the interval (0,2].}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 variogram remains unmodified.}
}
\details{
 
 Originally, the logarithmic model \eqn{\gamma(r) = \log(r)} was named
 after de Wijs and reflects a principle of similarity (cf. Chiles,
 J.-P. and Delfiner, P. (1999), p. 90).
 But note that \eqn{\gamma(r) = \log(r)} is not a valid variogram
 (\eqn{\gamma(0)} does not vanish) and can only be understood as a
 characteristic of a generalised random field.

 The modified \code{RMdewijsian} model
 \eqn{\gamma(r) = \log(r^{\alpha}+1)} is a valid variogram model
 (cf. Wackernagel, H. (2003), p. 336).
}

\value{
 \command{\link{RMdewijsian}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}

\note{Note that the (non-modified) de Wijsian model equals
\eqn{\gamma(r) = \log(r)}.}

\references{
 \itemize{
% \item Chiles, J.-P. and Delfiner, P. (1999)
% \emph{Geostatistics. Modeling Spatial Uncertainty.}
% New York: Wiley.

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
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
model <- RMdewijsian(alpha=1)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}