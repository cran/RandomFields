\name{RMmult}
\alias{RMmult}
\alias{*}
\alias{RM_MULT}
  
\title{Multiplication of Random Field Models}
\description{
 \command{\link{RMmult}} is a multivariate covariance model which depends on
 up to 10 submodels \eqn{C_0, C_1, ..., C_9}.
 In general, realizations of the created \command{\link{RMmodel}} are pointwise
 products of independent realizations of the submodels.
 
 In particular, if all submodels are given through a covariance
 function, the resulting model is defined through its covariance
 function, which is the product of the submodels' covariances.
}
\usage{
RMmult(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, var, scale, Aniso, proj)
}
\arguments{
 \item{C0}{an \command{\link{RMmodel}}.}
 \item{C1,C2,C3,C4,C5,C6,C7,C8,C9}{optional; each an \command{\link{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 model remains unmodified.}
}
\details{
 \command{\link{RMmodel}}s can also be multiplied via the
 \code{*}-operator, e.g. C0 * C1. 

 The global arguments \code{scale,Aniso,proj} of \command{\link{RMmult}}
 are multiplied to the corresponding argument of the submodels
 (from the right side). E.g.,
\cr
 \code{RMmult(Aniso=A1, RMexp(Aniso=A2), RMspheric(Aniso=A3))}
 \cr
 equals
 \cr
 \code{RMexp(Aniso=A2 \%*\% A1) * RMspheric(Aniso=A3 \%*\% A1)}
 
 In case that all submodels are given through a covariance function,
 the global argument \code{var} of \command{\link{RMmult}}
 is multiplied to the product covariance of \command{\link{RMmult}}.
}

\value{
 \command{\link{RMmult}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}


\me
\seealso{	
 \command{\link{RMplus}},	
 \command{\link{RMmodel}},
 \command{\link{RMprod}},
 \command{\link{RFsimulate}}, 
 \command{\link{RFfit}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

# separable, multiplicative model
model <- RMgauss(proj=1) * RMexp(proj=2, scale=5)
z <- RFsimulate(model=model, 0:10, 0:10, n=4)
plot(z)

\dontshow{FinalizeExample()}}


\keyword{spatial}
\keyword{models}


