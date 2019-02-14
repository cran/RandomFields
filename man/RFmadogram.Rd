\name{RFmadogram}
\alias{RFmadogram}
\title{Empirical (Cross-)Madogram}
\description{
 Calculates the empirical (cross-)madogram. The empirical
 (cross-)madogram of two random fields \eqn{X}{X} and \eqn{Y}{Y} is given by 
 \deqn{\gamma(r):=\frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r} |(X(t_{i})-X(t_{j}))||(Y(t_{i})-Y(t_{j}))|}{\gamma(r):=1/N(r) \sum_{(t_{i},t_{j})|t_{i,j}=r} |(X(t_{i})-X(t_{j}))||(Y(t_{i})-Y(t_{j}))|}
where \eqn{t_{i,j}:=t_{i}-t_{j}}{t_{i,j}:=t_{i}-t_{j}}, and where \eqn{N(r)}{N(r)} denotes the number of pairs of data points with distancevector 
\eqn{t_{i,j}=r}{t_{i,j}=r}. 

}
\usage{
RFmadogram(model, x, y=NULL, z=NULL, T=NULL, grid, params, distances,
           dim, ..., data, bin=NULL, phi=NULL, theta = NULL,
           deltaT = NULL, vdim=NULL)
}
\arguments{
 \item{model,params}{\argModel }
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{distances,dim}{\argDistances}
 \item{...}{\argDots}
 \item{data}{\argData}
 \item{bin}{\argBin}
 \item{phi}{\argPhi} 
 \item{theta}{\argTheta} 
 \item{deltaT}{\argDeltaT}
 \item{vdim}{\argVdim}
}
\details{ \command{\link{RFmadogram}} computes the empirical
 cross-madogram for given (multivariate) spatial data. 

 
 The spatial coordinates \code{x}, \code{y}, \code{z}
   should be vectors. For random fields of
 spatial dimension \eqn{d > 3} write all vectors as columns of matrix x. In
 this case do neither use y, nor z and write the columns in
 \code{gridtriple} notation.

 If the data is spatially located on a grid a fast algorithm based on
 the fast Fourier transformed (fft) will be used.
 As advanced option the calculation method can also be changed for grid
 data (see \command{\link{RFoptions}}.)
 
 It is also possible to use \command{\link{RFmadogram}} to calculate
 the pseudomadogram (see \command{\link{RFoptions}}).

}
\value{
 \command{\link{RFmadogram}} returns objects of class
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
 \command{\link{RMstable}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RFcov}},
 \command{\link{RFpseudomadogram}}.
\command{\link{RFvariogram}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

n <- 1 ## use n <- 2 for better results

## isotropic model
model <- RMexp()
x <- seq(0, 10, 0.02)
z <- RFsimulate(model, x=x, n=n)
emp.vario <- RFmadogram(data=z)
plot(emp.vario)


## anisotropic model
model <- RMexp(Aniso=cbind(c(2,1), c(1,1)))
x <- seq(0, 10, 0.05)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFmadogram(data=z, phi=4)
plot(emp.vario)


## space-time model
model <- RMnsst(phi=RMexp(), psi=RMfbm(alpha=1), delta=2)
x <- seq(0, 10, 0.05)
T <- c(0, 0.1, 100)
z <- RFsimulate(x=x, T=T, model=model, n=n)
emp.vario <- RFmadogram(data=z, deltaT=c(10, 1))
plot(emp.vario, nmax.T=3)


## multivariate model
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
x <- seq(0, 20, 0.1)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFmadogram(data=z)
plot(emp.vario)


## multivariate and anisotropic model
model <- RMbiwm(A=matrix(c(1,1,1,2), nc=2),
                nudiag=c(0.5,2), s=c(3, 1, 2), c=c(1, 0, 1))
x <- seq(0, 20, 0.1)
dta <- RFsimulate(model, x, x, n=n)
ev <- RFmadogram(data=dta, phi=4)
plot(ev, boundaries=FALSE)

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
