\name{RFempiricalvariogram}
\alias{RFempiricalvariogram}
\title{Empirical (Cross-)Variogram}
\description{
 Calculates the empirical (cross-)variogram. The empirical
 (cross-)variogram of two random fields \eqn{X}{X} and \eqn{Y}{Y} is given by 
 \deqn{\gamma(r):=\frac{1}{2N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r} (X(t_{i})-X(t_{j}))(Y(t_{i})-Y(t_{j}))}{\gamma(r):=1/2N(r) \sum_{(t_{i},t_{j})|t_{i,j}=r} (X(t_{i})-X(t_{j}))(Y(t_{i})-Y(t_{j}))}
where \eqn{t_{i,j}:=t_{i}-t_{j}}{t_{i,j}:=t_{i}-t_{j}}, and where \eqn{N(r)}{N(r)} denotes the number of pairs of data points with distancevector 
\eqn{t_{i,j}=r}{t_{i,j}=r}. 

}
\usage{
RFempiricalvariogram(x, y = NULL, z = NULL, T = NULL, data, grid,
 bin, phi, theta, deltaT, distances, vdim, ...)
}
\arguments{
 \item{x}{matrix of coordinates, or vector of x coordinates, or object
 of class \code{\link[sp:GridTopology-class]{GridTopology}}  or
   \code{\link[raster]{raster}}. 
 If matrix, \code{ncol(x)} is the dimension
 of the index space. Matrix notation is required in case of more than 3 spatial
 dimensions; in this case, if \code{grid=FALSE}, \code{x_ij} is the
 i-th coordinate in the j-th dimension. Otherwise, if
 \code{grid=TRUE}, the columns of \code{x} are interpreted as
 gridtriples (see \code{grid}). if of class \code{\link[sp:GridTopology-class]{GridTopology}} , \code{x} is interpreted as grid definition and \code{grid}
 is automatically set to \code{TRUE}. Coordinates are not required if
 the data is an object of class
 \command{\link[=RFsp-class]{RFsp}}, as these objects already
 contain its coordinates
 }
 \item{y}{optional vector of y coordinates, ignored if \code{x}
 is a matrix }
 \item{z}{optional vector of z coordinates, ignored if \code{x}
 is a matrix}
 \item{T}{optional vector of time coordinates, 
 \code{T} must always be an equidistant vector or given in a
 gridtriple format (see \code{grid}); for each component
 of \code{T}, the random field is simulated at all location points;
 the argument \code{T} is in an experimental stage. % to do
 }
 \item{grid}{logical; determines whether the vectors \code{x},
 \code{y}, and \code{z} or the columns of \code{x} should be
 interpreted as a grid definition (see Details). If \code{grid=TRUE},
 either \code{x}, \code{y}, and \code{z} must
 be equidistant vectors in ascending order or the columns of \code{x}
 must be given in the gridtriple format
 \code{c(from, stepsize, len)} (see Details);
 Not required if data is of class \code{\link[=RFsp-class]{RFsp}}
 }
 \item{data}{matrix, data.frame or object of class
 \command{\link[=RFsp-class]{RFsp}};
 
 } 
 \item{bin}{a vector giving the borders of the bins; If not specified
 an array describing the empirical (pseudo-)(cross-) variogram 
 in every direction is returned.}
 \item{phi}{an integer defining the number of sectors one half of the
 X/Y plane shall be devided into. If not specified, either an array is
 returned (if bin missing) or isotropy is assumed (if bin specified)} 
 \item{theta}{an integer defining the number of sectors one half of the
 X/Z plane shall be devided into. Use only for dimension \eqn{d=3}{d=3}
 if phi is already specified} 
 \item{deltaT}{vector of length 2, specifying the temporal bins. 
   The internal bin vector becomes \code{seq(from=0, to=deltaT[1], by=deltaT[2])}}
 \item{distances}{object of class \code{\link[stats]{dist}} representing
 the upper trianguar part of the matrix of Euclidean distances
 between the points at which the field is to be simulated; only
 applicable for stationary and isotropic models; if not \code{NULL},
 \code{dim} must be given and \code{x}, \code{y}, \code{z} and
 \code{T} must be missing or \code{NULL}. 
 }
 \item{vdim}{the number of variables of a multivariate data set.   
   If not given and \code{data} is an \code{RFsp} object
   created by \pkg{RandomFields}, the information there is
   taken from there. Otherwise \code{vdim} is assumed to be one.

   NOTE: still the argument \code{vdim} is an experimental stage. %to do
 }
 \item{...}{further options and control arguments for the simulation
 that are passed to and processed by \command{\link{RFoptions}}.}
}
\details{ \command{\link{RFempiricalvariogram}} computes the empirical
 cross-variogram for given (multivariate) spatial data. 

 
 The spatial coordinates \code{x}, \code{y}, \code{z}
   should be vectors. For random fields of
 spatial dimension \eqn{d > 3} write all vectors as colums of matrix x. In
 this case do neither use y, nor z and write the colums in
 \code{gridtriple} notation.

 If the data is spatially located on a grid a fast algorithm based on
 the fast Fourier transformed (fft) will be used.
 As advanced option the calculation method can also be changed for grid
 data (see \command{\link{RFoptions}}.)
 
 It is also possible to use \command{\link{RFempiricalvariogram}} to calulate
 the pseudovariogram (see \command{\link{RFoptions}}).

}
\value{
 \command{\link{RFempiricalvariogram}} returns objects of class
 \command{\link[=RFempVariog-class]{RFempVariog}}. 
}
\references{
 Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.

 Stein, M. L. (1999) \emph{Interpolation of Spatial Data.}
 New York: Springer-Verlag 
 }


\author{
Sebastian Engelke, \email{sebastian.engelke@unil.ch}
%\url{https://hec.unil.ch}

Johannes Martini, \email{jmartin2@uni-goettingen.de}

Martin Schlather, \email{schlather@math.uni-mannheim.de}
\url{http://ms.math.uni-mannheim.de/de/publications/software}
}
 
\seealso{
 \command{\link{RMstable}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

n <- 1 ## use n <- 2 for better results

## isotropic model
model <- RMexp()
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
z <- RFsimulate(model, x=x, n=n)
emp.vario <- RFempiricalvariogram(data=z)
plot(emp.vario, model=model)


## anisotropic model
model <- RMexp(Aniso=cbind(c(2,1), c(1,1)))
x <- seq(0, 10, if (interactive()) 0.05 else 1)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFempiricalvariogram(data=z, phi=4)
plot(emp.vario, model=model)


## space-time model
model <- RMexp()
x <- seq(0, 10, if (interactive()) 0.05 else 1)
T <- c(0, if (interactive()) 0.1 else 1, if (interactive()) 100 else 10)
z <- RFsimulate(x=x, T=T, model=model, n=n)
emp.vario <- RFempiricalvariogram(data=z, deltaT=c(10, 1))
plot(emp.vario, model=model, nmax.T=3)


## multivariate model
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
x <- seq(0, 20, if (interactive()) 0.1 else 2)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFempiricalvariogram(data=z)
plot(emp.vario, model=model)


## multivariate and anisotropic model
model <- RMbiwm(A=matrix(c(1,1,1,2), nc=2),
                nudiag=c(0.5,2), s=c(3, 1, 2), c=c(1, 0, 1))
x <- seq(0, 20, if (interactive()) 0.1 else 2)
data <- RFsimulate(model, x, x, n=n)
ev <- RFempiricalvariogram(data=data, phi=4)
plot(ev, model=model, boundaries=FALSE)


\dontshow{FinalizeExample()}
}

\keyword{spatial}
\keyword{models}
