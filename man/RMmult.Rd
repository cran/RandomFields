\name{RMmult}
\alias{RMmult}
\title{Multiplication of Random Field Models}
\description{
 \command{\link{RMmult}} is a multivariate covariance model which depends on
 up to 10 submodels \eqn{C_0, C_1, ..., C_10}.
 In general, realizations of the created \command{\link{RMmodel}} are pointwise
 product of independent realizations of the submodels.
 
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
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 model remains unmodified.}
}
\details{
 The global parameters \code{scale,Aniso,proj} of \command{\link{RMmult}}
 are multiplied to the corresponding parameters of the submodels
 (from the right side).
 
 In case that all submodels are given through a covariance function,
 the global parameter \code{var} of \command{\link{RMmult}}
 is multiplied to the product covariance of \command{\link{RMmult}}.
}
\value{
 \command{\link{RMmult}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{	
 \command{\link{RMplus}},	
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{
set.seed(0)

# separable, multiplicative model
model <- RMmult(RMgauss(proj=1), RMexp(proj=2, scale=5))
z <- RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4)
plot(z)


}


\keyword{spatial}
\keyword{models}


