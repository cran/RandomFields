\name{Independent Variables}
\alias{Nugget}
\alias{RPnugget}
\title{Method to simulate the Nugget effect}



\description{
 Method to simulate the Nugget effect.
}

\usage{
RPnugget(phi, loggauss, tol, vdim)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 specifies the covariance model to be simulated. The only possible
 model for \code{phi} is \command{\link{RMnugget}}.}
 \item{loggauss}{optional parameters;
 same meaning as for 
 \command{\link{RPgauss}}.}
 \item{tol}{
   points at a distance less than or equal to \code{nugget.tol}
 are considered as being identical. This strategy applies to
 the simulation method and the covariance function itself.
 Hence, the covariance function is only positive definite
 if \code{nugget.tol=0.0}. However, if the anisotropy matrix
 does not have full rank and \code{nugget.tol=0.0} then,
 the simulations are likely to be odd. 
 The value of \code{nugget.tol}
 should be of order \eqn{10^{-15}}{1e-15}.
 
 Default: \code{0.0}
 }
 \item{vdim}{positive integer; the model is treated
   \code{vdim}-variate, \code{vdim=1} (default) corresponds to a
   univariate random field.
   Mostly, the value of \code{vdim} is set automatically.

   Default is that it takes the value of the submodel \code{phi}}
}
 

\details{ 
 This method only allows \command{\link{RMnugget}} as a submodel.

 The method also allows for zonal nugget effects. Only there the
 parameter \code{tol} becomes important.
 For the zonal nugget effect, the anisotropy matrix \code{Aniso}
 should be given in \command{\link{RMnugget}}. There, only the kernal of the
 matrix is important.
}

\value{
 \code{RPnugget} returns an object of class
 \code{\link[=RMmodel-class]{RMmodel}} 

}

\references{
 \itemize{
 \item
 Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University. 
 }
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\examples{
RFoptions(seed=0)
model <- RMnugget()
z <- RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4)
plot(z)

model <- RPnugget(RMnugget(var=0.01, Aniso=matrix(nc=2, rep(1,4))))
z <- RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4)
plot(z)
\dontshow{RFoptions(seed=NA)}
}

\seealso{
  \link{RP},
 \command{\link{RPcoins}},
 \command{\link{RPhyperplane}},
 \command{\link{RPspectral}},
 \command{\link{RPtbm}}.
}

\keyword{methods}
