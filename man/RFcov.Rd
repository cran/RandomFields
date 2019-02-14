\name{RFcov}
\alias{RFcov}
\title{(Cross-)Covariance function}
\description{
 Calculates both the empirical and the theoretical (cross-)covariance function.
}
\usage{
RFcov(model, x, y = NULL, z = NULL, T=NULL, grid, params, distances, dim, ...,
      data, bin=NULL, phi=NULL, theta = NULL, deltaT = NULL, vdim=NULL)
}
\arguments{
  \item{model,params}{\argModel }
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{data}{\argData}
 \item{bin}{\argBin}
 \item{phi}{\argPhi} 
 \item{theta}{\argTheta} 
 \item{deltaT}{\argDeltaT}
 \item{distances,dim}{\argDistances}
 \item{vdim}{\argVdim}
 \item{...}{\argDots}
}
\details{ \command{\link{RFcov}} computes the empirical
 cross-covariance function for given (multivariate) spatial data. 

 The empirical
 (cross-)covariance function of two random fields \eqn{X}{X} and \eqn{Y}{Y} is given by 
 \deqn{\gamma(r):=\frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
   (X(t_{i})Y(t_{j})) - m_{X} m_{Y}}{\gamma(r):=1/N(r)
   \sum_{(t_{i},t_{j})|t_{i,j}=r} (X(t_{i})Y(t_{j})) - m_{X} m_{Y}}
where \eqn{t_{i,j}:=t_{i}-t_{j}}{t_{i,j}:=t_{i}-t_{j}}, \eqn{N(r)}{N(r)}
denotes the number of pairs of data points with distancevector
\eqn{t_{i,j}=r}{t_{i,j}=r} and where
\eqn{m_{X} := \frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
  X_{t_{i}}}{m_{X} := \frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
  X_{t_{i}}} and \eqn{m_{Y} := \frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
  Y_{t_{i}}}{m_{Y} := 1/N(r) \sum_{(t_{i},t_{j})|t_{i,j}=r}
  Y_{t_{i}}} denotes the mean of data points with distancevector
\eqn{t_{i,j}=r}{t_{i,j}=r}.

 The spatial coordinates \code{x}, \code{y}, \code{z}
   should be vectors. For random fields of
 spatial dimension \eqn{d > 3} write all vectors as columns of matrix x. In
 this case do neither use y, nor z and write the columns in
 \code{gridtriple} notation.

 If the data is spatially located on a grid a fast algorithm based on
 the fast Fourier transformed (fft) will be used.
 As advanced option the calculation method can also be changed for grid
 data (see \command{\link{RFoptions}}.)
 
 It is also possible to use \command{\link{RFcov}} to calculate
 the pseudocovariance function (see \command{\link{RFoptions}}).

}
\value{
 \command{\link{RFcov}} returns objects of class
 \command{\link[=RFempVariog-class]{RFempVariog}}. 
}
\references{
 Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.

 Stein, M. L. (1999) \emph{Interpolation of Spatial Data.}
 New York: Springer-Verlag 
 }


 \author{Jonas Auel; Sebastian Engelke; Johannes Martini; \martin}
 
\seealso{
  \command{\link{RFvariogram}},
  \command{\link{RFmadogram}},
 \command{\link{RMstable}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

n <- 1 ## use n <- 2 for better results

## isotropic model
model <- RMexp()
x <- seq(0, 10, 0.02)
z <- RFsimulate(model, x=x, n=n)
emp.vario <- RFcov(data=z)
plot(emp.vario, model=model)


## anisotropic model
model <- RMexp(Aniso=cbind(c(2,1), c(1,1)))
x <- seq(0, 10, 0.05)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFcov(data=z, phi=4)
plot(emp.vario, model=model)


## space-time model
model <- RMnsst(phi=RMexp(), psi=RMfbm(alpha=1), delta=2)
x <- seq(0, 10, 0.05)
T <- c(0, 0.1, 100)
z <- RFsimulate(x=x, T=T, model=model, n=n)
emp.vario <- RFcov(data=z, deltaT=c(10, 1))
plot(emp.vario, model=model, nmax.T=3)


## multivariate model
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
x <- seq(0, 20, 0.1)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFcov(data=z)
plot(emp.vario, model=model)


## multivariate and anisotropic model
model <- RMbiwm(A=matrix(c(1,1,1,2), nc=2),
                nudiag=c(0.5,2), s=c(3, 1, 2), c=c(1, 0, 1))
x <- seq(0, 20, 0.1)
dta <- RFsimulate(model, x, x, n=n)
ev <- RFcov(data=dta, phi=4)
plot(ev, model=model, boundaries=FALSE)


\dontshow{\dontrun{
######################################################
# NOTE: distinguish the different uses of x and y
x <- c(1,2,1)
y <- c(4,5,6)
# coordinate space 1-dimensional, evaluated at 3 points:
RFcov(model=model, x=as.matrix(x), y=as.matrix(y))
# coordinate space is 3-dimensional, evaluated at a pair of points
RFcov(model=model, x=t(x), y=t(y))
}}

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
