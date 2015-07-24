\name{RMnonstwm}
\alias{RMnonstwm}
\title{RMnonstwm}
\description{
  The non-stationary Whittle-Matern model
  \eqn{C} is given by 
 
 \deqn{C(x, y)=\Gamma(\mu) \Gamma(\nu(x))^{-1/2} \Gamma(\nu(y))^{-1/2}
 W_{\mu} (|x-y|)}{C(x, y)=\Gamma(\mu) \Gamma(\nu(x))^{-1/2}
 \Gamma(\nu(y))^{-1/2} 
 W_{\mu} (|x-y|)}
 
 where \eqn{\mu = [\nu(x) + \nu(y)]/2}, and
 \eqn{\nu} must a positive function.

 \eqn{W_{\mu}} is the
 covariance of the \command{\link{RMwhittle}} model or the
 \command{\link{RMmatern}} model and
}

\details{
  The non-stationary Whittle-Matern models are obtained by
  the respective stationary model, replacing the 
  real-valued argument for \code{nu}  by a non-negative function.
}

\note{
  It cannot be checked whether \code{nu} only takes positive number. So
  the responsibility is completely left to the user.
}

\references{
 \itemize{
 \item Stein, M. (2005)
 Nonstationary Spatial Covariance Functions.
 Tech. Rep., 2005
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMwhittle}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
x <- seq(-1.2, 1.2, len=50)
model <- RMwhittle(nu=RMgauss())
%plot(model, dim=2, xlim=c(-1,1))
z <- RFsimulate(model=model, x, x, n=4)
plot(z)
\dontshow{FinalizeExample()}
}


\keyword{spatial}
\keyword{models}


