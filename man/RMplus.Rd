\name{RMplus}
\alias{RMplus}
\alias{+}
\title{Addition of Random Field Models}
\description{
 \command{\link{RMplus}} is an additive covariance model which depends on
 up to 10 submodels \eqn{C_0, C_1, ..., C_10}.
 In general, realizations of the created \command{\link{RMmodel}} are pointwise
 sums of independent realizations of the submodels.
 
 In particular, if all submodels are given through a covariance
 function, the resulting model is defined through its covariance
 function, which is the sum of the submodels' covariances.
 Analogously, if all submodels are given through a variogram.
}
\usage{
RMplus(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, var, scale, Aniso, proj)
}
\arguments{
 \item{C0}{an \command{\link{RMmodel}}.}
 \item{C1,C2,C3,C4,C5,C6,C7,C8,C9}{optional; each an \command{\link{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 model remains unmodified.}
}
\details{
 \command{\link{RMmodel}}s can also be summed up via the
 \code{+}-operator, e.g.: C0 + C1 

 The global parameters \code{var,scale,Aniso,proj} of \command{\link{RMplus}}
 are multiplied to the corresponding parameters of the submodels
 (from the right side).
}
\value{
 \command{\link{RMplus}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmult}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}
\examples{
RFoptions(seed=0)
model <- RMgauss() + RMnugget(var=0.01)
z <- RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4)
plot(z)


# gives the same model as 
model <- RMplus(RMgauss(), RMnugget(var=0.01))
z <- RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4)
plot(z)
\dontshow{RFoptions(seed=NA)}
}

\keyword{spatial}
\keyword{models}
