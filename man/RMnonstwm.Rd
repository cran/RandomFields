\name{RMnonstwm}
\alias{RMnonstwm}
\title{RMnonstwm}
\description{
 \command{\link{RMnonstwm}} is a covariance model whose corresponding
 covariance \eqn{C}{C} is given by 
 
 \deqn{C(x, y)=\Gamma(\mu) \Gamma(\nu(x))^{-1/2} \Gamma(\nu(y))^{-1/2}
 W_{\mu} (|x-y|)}{C(x, y)=\Gamma(\mu) \Gamma(\nu(x))^{-1/2}
 \Gamma(\nu(y))^{-1/2} 
 W_{\mu} (|x-y|)}
 
 where \eqn{\mu = [\nu(x) + \nu(y)]/2} and \eqn{W_{\mu}} is the
 covariance of the \command{\link{RMwhittle}} model and
 \eqn{\nu} is a positive function or a constant.

}
\usage{
RMnonstwm(nu, var, scale, Aniso, proj)
}
\arguments{
 \item{nu}{positive argument defining the smoothness of the random
 field, or it is an arbitrary \command{\link[=RMmodelsAuxiliary]{shape function}}.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
%\details{}

\value{
 \command{\link{RMnonstwm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
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
\dontshow{\dontrun{
x <- seq(-1.2, 1.2, len=if (interactive()) 70 else 3)
model <- RMnonstwm(nu=RMgauss())
%plot(model, dim=2, xlim=c(-1,1))
z <- RFsimulate(model=model, x, x, n=4)
plot(z)
}}
\dontshow{FinalizeExample()}
}


\keyword{spatial}
\keyword{models}


