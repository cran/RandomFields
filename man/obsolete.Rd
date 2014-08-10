\name{Obsolete Functions}
\alias{CondSimu}
\alias{Covariance}
\alias{CovarianceFct}
\alias{CovMatrix}
\alias{DeleteRegister}
\alias{DeleteAllRegisters}
\alias{DoSimulateRF}
\alias{EmpiricalVariogram}
\alias{fitvario}
\alias{fractal.dim}
\alias{GaussRF}
\alias{hurst}
\alias{InitSimulateRF}
\alias{InitGaussRF}
\alias{InitMaxStableRF}
\alias{Kriging}
\alias{MaxStableRF}
\title{Obsolete functions}
\alias{RFparameters}
\alias{Variogram}
\description{
 This part gives the obsolete functions of RandomFields Version 2
}
\usage{
Covariance(x, y = NULL, model, param = NULL, dim = if
 (!missing(Distances)) { if (is.matrix(x)) ncol(x) else 1},
 Distances, fctcall = c("Cov", "Variogram", "CovMatrix"))
CovarianceFct(x, y = NULL, model, param = NULL, dim = if
 (!missing(Distances)) { if (is.matrix(x)) ncol(x) else 1},
 Distances, fctcall = c("Cov", "Variogram", "CovMatrix"))
CovMatrix(x, y = NULL, model, param = NULL, dim = if
 (!missing(Distances)) { if (is.matrix(x)) ncol(x) else 1}, Distances)
DeleteAllRegisters()
DeleteRegister(nr=0)
DoSimulateRF(n = 1, register = 0, paired=FALSE, trend=NULL) 
InitSimulateRF(x, y = NULL, z = NULL, T=NULL, grid=!missing(gridtriple),
 model, param, trend, method = NULL, register = 0,
 gridtriple, distribution=NA)
InitGaussRF(x, y = NULL, z = NULL, T=NULL, grid=!missing(gridtriple),
 model, param, trend=NULL, method = NULL, register = 0, gridtriple) 
GaussRF(x, y = NULL, z = NULL, T=NULL, grid=!missing(gridtriple), model,
 param, trend=NULL, method = NULL, n = 1, register = 0, gridtriple,
 paired=FALSE, PrintLevel=1, Storing=TRUE, ...) 
Variogram(x, model, param = NULL, dim = if (!missing(Distances))
 { if (is.matrix(x)) ncol(x) else 1}, Distances)
InitMaxStableRF(x, y = NULL, z = NULL, grid, model, param, maxstable,
 method = NULL, register = 0, gridtriple = FALSE)
MaxStableRF(x, y = NULL, z = NULL, grid, model, param, maxstable,
 method = NULL, n = 1, register = 0, gridtriple = FALSE, ...)
EmpiricalVariogram(x, y = NULL, z = NULL, T=NULL, data, grid, bin,
 gridtriple = FALSE, phi, theta, deltaT)
Kriging(krige.method, x, y=NULL, z=NULL, T=NULL, grid, gridtriple=FALSE,
 model, param, given, data, trend=NULL,pch=".", return.variance=FALSE,
 allowdistanceZero = FALSE, cholesky=FALSE) 
CondSimu(krige.method, x, y=NULL, z=NULL, T=NULL, grid, gridtriple=FALSE,
 model, param, method=NULL, given, data, trend=NULL, n=1, register=0, 
 err.model=NULL, err.param=NULL, err.method=NULL, err.register=1, 
 tol=1E-5, pch=".", paired=FALSE, na.rm=FALSE) 
RFparameters(...)
hurst(x, y = NULL, z = NULL, data, gridtriple = FALSE, sort=TRUE,
 block.sequ = unique(round(exp(seq(log(min(3000, dim[1] / 5)),
 log(dim[1]), len=min(100, dim[1]))))),
 fft.m = c(1, min(1000, (fft.len - 1) / 10)),
 fft.max.length = Inf, 
 method=c("dfa", "fft", "var"), mode=c("plot", "interactive"),
 pch=16, cex=0.2, cex.main=0.85,
 PrintLevel=RFoptions()$general$printlevel,height=3.5, ...)
fractal.dim(x, y = NULL, z = NULL, data, grid=TRUE, gridtriple = FALSE,
 bin, vario.n=5, sort=TRUE, fft.m = c(65, 86), fft.max.length=Inf,
 fft.max.regr=150000, fft.shift = 50, method=c("variogram", "fft"),
 mode=c("plot", "interactive"), pch=16, cex=0.2, cex.main=0.85,
 PrintLevel = RFoptions()$general$printlevel, height=3.5, ...)
fitvario(x, y=NULL, z=NULL, T=NULL, data, model, param, lower=NULL,
 upper=NULL, sill=NA, grid=!missing(gridtriple), gridtriple=FALSE, ...)
}
\arguments{
 \item{x, y, model, param, dim, Distances, fctcall, n, register,
 paired, trend, z, T, grid, method, gridtriple, distribution,
 PrintLevel, Storing, ..., maxstable, data, bin, 
	phi, theta, deltaT, krige.method, pch, return.variance,
 allowdistanceZero, cholesky, given,
	err.model, err.param, err.method, err.register, 
	tol, na.rm, sort,
	block.sequ, fft.m, fft.max.length, mode,
	cex, cex.main, height, vario.n, fft.max.regr, fft.shift,
	lower, upper, sill, nr}{
	as the functions are obsolete, all these arguments
	are not documented anymore.
 }
}
\value{
 not given anymore.
}
\seealso{
  The functions that should be used instead are, for instance,
 \command{\link{RFcov}}, \command{\link{RFcovmatrix}},
 \command{\link{RFvariogram}}, \command{\link{RFsimulate}},
 \command{\link{RFinterpolate}},
 \command{\link{RFempiricalvariogram}}, \command{\link{RFfit}},
 \command{\link{RFoptions}}, \command{\link{RFhurst}},
 \command{\link{RFfractaldim}} 
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

%\dontshow{RFoptions(internal.warn_newstyle=FALSE, internal.warn_oldstyle=FALSE)}

## no examples given, as command is obsolete
\dontshow{\dontrun{

x <- seq(0, 5, if (interactive()) 0.05 else 2)
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
\dontshow{FinalizeExample()}
}

\keyword{spatial}
