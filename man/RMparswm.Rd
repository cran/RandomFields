\name{RMparswm}
\alias{RMparswm}
\alias{RMparswmX}
\title{Parsimonious Multivariate Whittle Matern Model}
\description{
 \command{\link{RMparswm}} is a multivariate stationary isotropic
 covariance model 
 whose corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between
 two points and is given for \eqn{i,j \in \{1,2\}}{i,j = 1,2} by
 \deqn{C_{ij}(r)= c_{ij} W_{\nu_{ij}}(r).}{C_{ij}(r)=c_{ij} W_{\nu_{ij}}(r).}
 Here \eqn{W_\nu} is the covariance of the
 \command{\link{RMwhittle}} model.

 \command{RMparswmX} ist defined as
 \deqn{\rho_{ij} C_{ij}(r)}
 where \eqn{\rho_{ij}} is any covariance matrix.
}

\usage{
RMparswm(nudiag, var, scale, Aniso, proj)
RMparswmX(nudiag, rho, var, scale, Aniso, proj)
}
\arguments{
 \item{nudiag}{a vector of arbitrary length of positive values; each entry
 positive; the vector \eqn{(\nu_{11},\nu_{22},...)}.
 The offdiagonal elements \eqn{\nu_{ij}} are calculated as
 \eqn{0.5 (\nu_{ii} + \nu_{jj})}.}
 \item{rho}{any positive definite \eqn{m \times m}{m x m}
 matrix; 
 here \eqn{m} equals \code{length(nudiag)}
 For the calculation of \eqn{c_{ij}} see Details.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.
 }
}
\details{
 In the equation above we have
 \deqn{c_{ij} = \rho_{ij} \sqrt{G_{ij}}
 }{
 c_{ij} = \rho_{ ij} \sqrt G_{ij}}
 and
 \deqn{G_{ij} = \frac{\Gamma(\nu_{11} + d/2) \Gamma(\nu_{22} + d/2)
 \Gamma(\nu_{12})^2}{\Gamma(\nu_{11}) \Gamma(\nu_{22})
 \Gamma(\nu_{12}+d/2)^2}
 }{
 G_{ij} = \Gamma(\nu_{11} + d/2) \Gamma(\nu_{22} + d/2) \Gamma(\nu_{12}) /
 (\Gamma(\nu_{11}) \Gamma(\nu_{22}) \Gamma(\nu_{12}+d/2))^2)}
 where \eqn{\Gamma} is the Gamma function and \eqn{d} is the dimension
 of the space.

 Note that the definition of \command{RMparswmX} is
 \code{RMschur(M=rho, RMparswm(nudiag, var, scale, Aniso, proj))} .
}


\value{
 \command{\link{RMparswm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Gneiting, T., Kleiber, W., Schlather, M. (2010)
 Matern covariance functions for multivariate random fields
 \emph{JASA}
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMbiwm}},
 \command{\link{RMwhittle}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\examples{

RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}

rho <- matrix(nc=3, c(1, 0.5, 0.2, 0.5, 1, 0.6, 0.2, 0.6, 1))
model <- RMparswmX(nudiag=c(1.3, 0.7, 2), rho=rho)
x.seq <- y.seq <- seq(-10, 10, 0.1)
z <- RFsimulate(model = model, x=x.seq, y=y.seq)
plot(z)

\dontshow{FinalizeExample()}
} % examples

\keyword{spatial}
\keyword{models}



