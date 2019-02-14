\name{RMplus}
\alias{RMplus}
\alias{+}
\alias{RM_PLUS}
\title{Addition of Random Field Models}
\description{
 \command{\link{RMplus}} is an additive covariance model which depends on
 up to 10 submodels \eqn{C_0, C_1, ..., C_9}.
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
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 model remains unmodified.}
}
\details{
 \command{\link{RMmodel}}s can also be summed up via the
 \code{+}-operator, e.g. C0 + C1. 

 The global arguments \code{var,scale,Aniso,proj} of \command{\link{RMplus}}
 are multiplied to the corresponding arguments of the submodels
 (from the right side).
}
\value{
 \command{\link{RMplus}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}


\me
\seealso{
 \command{\link{RMmult}},
 \command{\link{RMmodel}},
 \command{\link{RMsum}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMplus(RMgauss(), RMnugget(var=0.1))
model2<- RMgauss() + RMnugget(var=0.1)
plot(model, "model.+"=model2, type=c("p", "l"), pch=20, xlim=c(0,3)) # the same

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
