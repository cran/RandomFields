\name{RMmatrix}
\alias{RMmatrix}
\alias{coregionalisation}
\alias{coregionalization}
%\alias{diag}
%\alias{diag,RMmodel,missing,missing,missing-method}
\title{Matrix operator}
\description{
 \command{\link{RMmatrix}} is a multivariate covariance model
 depending on one multivariate covariance model, or
 one or several univariate covariance models \eqn{C0,\ldots}{C0,\ldots}.
 The corresponding covariance function is given by
 \deqn{ C(h) = M \phi(h) M^t}{ C(h) = M phi(h) M^t}
 if a multivariate case is given. Otherwise it returns a matrix
 whose diagonal elements are filled with the univarate model(s)
 \code{C0}, \code{C1}, etc, and the
 offdiagonals are all zero.
}
\usage{
RMmatrix(C0, C1,  C2, C3, C4, C5, C6, C7, C8, C9, M, vdim,
         var, scale, Aniso, proj)
%\S4method{diag}{RMmodel,missing,missing,missing}(x,nrow,ncol,names)
}
\arguments{
  \item{C0}{a k-variate covariance \command{\link{RMmodel}}
    or a univariate model or a list of models joined by
    \command{\link{c}}
  }
  \item{C1,C2,C3,C4,C5,C6,C7,C8,C9}{optional univariate models}
  \item{M}{a k times k matrix, which is multiplied from left and right to
    the given model; \eqn{M} may depend on the location, hence it is then a
    matrix-valued function and \eqn{C} will be non-stationary with
    \deqn{C(x, y) = M(x) \phi(x, y) M(y)^t}{ C(x, y) = M(x) phi(x, y) M(y)^t}}
 \item{vdim}{positive integer. This argument should be given if and only
    if a multivariate model is created from a single univariate model and
    \code{M} is not given. (In fact, if \code{M} is given, \code{vdim} must
    equal the number of columns of \code{M})}
  \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
    \command{\link{RMmodel}}. If not passed, the above
    covariance function remains unmodified.}
  % \item{x}{an \link{R.c} list of \link{RMmodel}s}
%  \item{nrow,ncol,names}{ignored arguments}
}
\note{
  \itemize{
%    \item
%    The more intuitive code
%    \code{diag(c(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9))}
%    is identical to \code{RMmatrix(C0, C1, C2, C3, C4, C5, C6, C7, C8,
%      C9)}. Here \code{C0},..., \code{C9} must be univariate models.
    
    \item\command{RMmatrix} also allows variogram models are arguments.
  }
}
\value{\command{\link{RMmatrix}} returns an object of class
  \code{\link[=RMmodel-class]{RMmodel}}.}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}

RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


\dontrun{
## first example: bivariate Linear Model of Coregionalisation
x <- y <- seq(0, 10, 0.2)

model1 <- RMmatrix(M = c(0.9, 0.43), RMwhittle(nu = 0.3)) + 
  RMmatrix(M = c(0.6, 0.8), RMwhittle(nu = 2))
plot(model1)
simu1 <- RFsimulate(RPdirect(model1), x, y)
plot(simu1)


## second, equivalent way of defining the above model
model2 <- RMmatrix(M = matrix(ncol=2, c(0.9, 0.43, 0.6, 0.8)),
                  c(RMwhittle(nu = 0.3), RMwhittle(nu = 2)))
simu2 <- RFsimulate(RPdirect(model2), x, y)
stopifnot(all.equal(as.array(simu1), as.array(simu2)))


## third, equivalent way of defining the above model
model3 <- RMmatrix(M = matrix(ncol=2, c(0.9, 0.43, 0.6, 0.8)),
                   RMwhittle(nu = 0.3), RMwhittle(nu = 2))
simu3 <- RFsimulate(RPdirect(model3), x, y)
stopifnot(all(as.array(simu3) == as.array(simu2)))
}


## second example: bivariate, independent fractional Brownian motion
## on the real axis
x <- seq(0, 10, 0.1) % to do: naechste Zeile RMmatrix durch diag ersetzen
modelB <- RMmatrix(c(RMfbm(alpha=0.5), RMfbm(alpha=1.5))) ## see the Note above
print(modelB)
simuB <- RFsimulate(modelB, x)
plot(simuB)


## third example: bivariate non-stationary field with exponential correlation
## function. The variance of the two components is given by the
## variogram of fractional Brownian motions.
## Note that the two components have correlation 1.
x <- seq(0, 10, 0.1)
modelC <- RMmatrix(RMexp(), M=c(RMfbm(alpha=0.5), RMfbm(alpha=1.5))) 
print(modelC)
simuC <- RFsimulate(modelC, x, x, print=1)
#print(as.vector(simuC))
plot(simuC)

\dontshow{FinalizeExample()}
}