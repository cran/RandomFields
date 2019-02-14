\name{Obsolete Functions Version 3}
\alias{RFempiricalvariogram}
\alias{RFempiricalcovariance}
\alias{RFempiricalmadogram}
\alias{RMstrokorbMono}
\alias{RMstrokorbBall}
\alias{RMstrokorbPoly}
\title{Obsolete functions Version 3}
\description{
  Some functions of RandomFields Version 3 have been replaced by more
  powerful functions
}
\usage{
RFempiricalvariogram(...)
RFempiricalcovariance(...)
RFempiricalmadogram(...)
}
\arguments{
 \item{...}{See for the recent functions given in the Details}
}

\details{
  \describe{
    \item{RFempiricalvariogram}{see \command{\link{RFvariogram}}}
    \item{RFempiricalcovariance}{see \command{\link{RFcov}}}
    \item{RFempiricalmadogram}{see \command{\link{RFmadogram}}}
 %  \item{RFempirical} see \command{\link{RF}}    
 % \item{RFempirical} see \command{\link{RF}}
 \item{Strokorb's M3/M4 functions}{are called
 \command{\link{RMm2r}}, \command{\link{RMm3b}},
  \command{\link{RMmps}}}
} 
   
}

\value{
  see the respective recent function
}
\seealso{
  \command{\link{RFcov}},
  \command{\link{RFcovmatrix}},
  \command{\link{RFvariogram}}.
  \command{\link{RFpseudovariogram}},
  \command{\link{RFmadogram}},
  \command{\link{RFpseudomadogram}}
}

\me


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

%\dontshow{RFoptions(internal.warn_newstyle=FALSE,
%internal.warn_oldstyle=FALSE)}



## no examples given, as command is obsolete
\dontshow{\dontrun{ # ok
x <- seq(0, 5, 0.05) 
model <- "exp"
param <- c(0, 1, 0, 1)

plot(x, Covariance(x, model=model, param=param))
plot(x, Variogram(x, model=model, param=param))

z <- GaussRF(x, model=model, param=param)
plot(x, z, type="l")

z <- GaussRF(x=x, y=x, model=model, param=param, grid=TRUE)
image(x, x, z)
EmpiricalVariogram(x, x, data=z)

g1 <- runif(5) * 2
g2 <- runif(5) * 2
z <- GaussRF(g1, g2, grid=FALSE, model=model, param=param)
k <- Kriging("S", x, x, given=cbind(g1, g2),
        data=z, model=model, param=param)
Print(range(z), range(k))
col <- rainbow(20)
zlim <- range(z, k)
image(x, x, k, col=col, zlim=zlim)
points(g1, g2, pch=20,cex=2)
points(g1, g2, pch=16, col=col[pmax(1, (z-zlim[1]) / diff(zlim) * 20 )])
#%   farben stimmen in etwa.

estparam <- rep(NA, 4)
v <- fitvario(x=g1, y=g2, grid=FALSE, model=model, param=estparam, data=z)
Print(v$ml, v$ml$ml.value)
}}
\dontshow{FinalizeExample()}}

\keyword{spatial}
