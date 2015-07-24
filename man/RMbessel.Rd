\name{RMbessel}
\alias{RMbessel}
\alias{RMjbessel}
\title{Bessel Family Covariance Model}
\description{
 \command{\link{RMbessel}} is a stationary isotropic covariance model
 belonging to the Bessel family. 
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = 2^\nu \Gamma(\nu+1) r^{-\nu} J_\nu(r)}{C(r) = 2^\nu \Gamma(\nu+1) r^{-\nu} J_\nu(r)}
 where \eqn{\nu \ge \frac{d-2}2}{\nu \ge (d-2)/2},
 \eqn{\Gamma} denotes the gamma function and
 \eqn{J_\nu}{J_\nu} is a Bessel function of first kind.
}
\usage{
RMbessel(nu, var, scale, Aniso, proj)
}
\arguments{
 \item{nu}{a numerical value; should be equal to or greater than
 \eqn{\frac{d-2}2}{(d-2)/2} to provide a valid
 covariance function for a random field of dimension \eqn{d}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This covariance models a hole effect (cf. Chiles, J.-P. and Delfiner,
 P. (1999), p. 92, cf. Gelfand et al. (2010), p. 26).
 
 An important case is \eqn{\nu=-0.5}{\nu=-0.5} 
 which gives the covariance function
 \deqn{C(r)=\cos(r)}{C(r)=cos(r)}
 and which is only valid for \eqn{d=1}{d=1}. This equals \command{\link{RMdampedcos}} for \eqn{\lambda = 0}, there.

 A second important case is \eqn{\nu=0.5}{\nu=0.5} with covariance function
 \deqn{C(r)=\sin(r)/r}{C(r)=sin(r)/r} and which is valid for \eqn{d \le 3}{d \le 3}. 
 This coincides with \command{\link{RMwave}}. 
 
 Note that all valid continuous stationary isotropic covariance
 functions for \eqn{d}{d}-dimensional random fields
 can be written as scale mixtures of a Bessel type
 covariance function with \eqn{\nu=\frac{d-2}2}{\nu=(d-2)/2}
 (cf. Gelfand et al., 2010, pp. 21--22).
}
\value{
 \command{\link{RMbessel}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
  \itemize{
    \item Chiles, J.-P. and Delfiner, P. (1999)
    \emph{Geostatistics. Modeling Spatial Uncertainty.}
    New York: Wiley.
    
    \item Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
    P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
    Boca Raton: Chapman & Hall/CRL.
    
    \item
    \url{http://homepage.tudelft.nl/11r49/documents/wi4006/bessel.pdf}
  }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMdampedcos}},
 \command{\link{RMwave}},
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
model <- RMbessel(nu=1, scale=0.1)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}