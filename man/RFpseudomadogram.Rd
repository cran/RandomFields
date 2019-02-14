\name{RFpseudomadogram}
\alias{RFpseudomadogram}
\title{Empirical  Pseudomadogram}
\description{
 Calculates the empirical pseudomadogram. The empirical
 pseudomadogram of two random fields \eqn{X}{X} and \eqn{Y}{Y} is given by 
 \deqn{\gamma(r):=\frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r} |(X(t_{i})-X(t_{j}))||(Y(t_{i})-Y(t_{j}))|}{\gamma(r):=1/N(r) \sum_{(t_{i},t_{j})|t_{i,j}=r} |(X(t_{i})-X(t_{j}))||(Y(t_{i})-Y(t_{j}))|}
where \eqn{t_{i,j}:=t_{i}-t_{j}}{t_{i,j}:=t_{i}-t_{j}}, and where \eqn{N(r)}{N(r)} denotes the number of pairs of data points with distancevector 
\eqn{t_{i,j}=r}{t_{i,j}=r}. 

}
\usage{
RFpseudomadogram(model, x, y=NULL, z=NULL, T=NULL, grid, params, distances,
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
\details{ \command{\link{RFpseudomadogram}} computes the empirical
 pseudomadogram for given (multivariate) spatial data. 

 
 The spatial coordinates \code{x}, \code{y}, \code{z}
   should be vectors. For random fields of
 spatial dimension \eqn{d > 3} write all vectors as columns of matrix x. In
 this case do neither use y, nor z and write the columns in
 \code{gridtriple} notation.

 If the data is spatially located on a grid a fast algorithm based on
 the fast Fourier transformed (fft) will be used.
 As advanced option the calculation method can also be changed for grid
 data (see \command{\link{RFoptions}}.)
 
 It is also possible to use \command{\link{RFpseudomadogram}} to calculate
 the pseudomadogram (see \command{\link{RFoptions}}).

}
\value{
 \command{\link{RFpseudomadogram}} returns objects of class
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
 \command{\link{RFmadogram}}.
\command{\link{RFvariogram}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
n <- 2
x <- seq(0, 20, 0.1)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFpseudomadogram(data=z)
plot(emp.vario)

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
