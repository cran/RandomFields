\name{Tbm}
\alias{Tbm}
\alias{RPtbm}
\title{Turning Bands method}
\description{
 The Turning Bands method is a simulation method for stationary, isotropic
 random fields in any dimension and defined on arbitrary points or
 arbitrary grids. It performs a multidimensional simulation
 by superposing lower-dimensional fields. In fact, the Turning Bands
 method is called with the Turning Bands model, see  
 \command{\link{RMtbm}}.
 \cr
 For details see \command{\link{RMtbm}}.
}

\usage{
RPtbm(phi, fulldim, reduceddim, layers, lines,
      linessimufactor, linesimustep, center, points)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
   specifies the covariance function to be simulated;
   a univariate stationary isotropic covariance model
   (see \code{RFgetModelNames(type="positive definite",
     domain="single variable", isotropy="isotropy",  vdim=1)})
 which is valid in dimension \code{fulldim}.
 }
% \item{loggauss}{see \command{\link{RPgauss}}.}
 \item{fulldim}{a positive integer. The dimension of the space of the
 random field to be simulated}

 \item{reduceddim}{a positive integer; less than \code{fulldim}. 
   The dimension of the auxiliary hyperplane (most frequently a line,
   i.e. \code{reduceddim=1} used in the simulation.
}

 \item{layers}{a boolean value; for space time model. If \code{TRUE} then the turning layers are
 used whenever a time component is given.
 If \code{NA} the turning layers are used only when the
 traditional TBM is not applicable.
 If \code{FALSE} then turning layers may never be used.
 
 Default: \code{TRUE}.
 }

 \item{lines}{
 Number of lines used.
 Default: \code{60}.
 }
 
 \item{linessimufactor}{ \code{linessimufactor} or
 \code{linesimustep} must be non-negative; if
 \code{linesimustep}
 is positive then \code{linesimufactor} is ignored.
 If both
 arguments are naught then \code{points} is used (and must be
 positive).
 The grid on the line is \code{linesimufactor}-times
 finer than the smallest distance. 
 See also \code{linesimustep}.
 
 Default: \code{2.0}.
 }
 \item{linesimustep}{
 If \code{linesimustep} is positive the grid on the line has lag
 \code{linesimustep}. 
 See also \code{linesimufactor}.
 
 Default: \code{0.0}.
 }
 \item{center}{Scalar or vector.
 If not \code{NA}, the \code{center} is used as the center of
 the turning bands for \code{fulldim}.
 Otherwise the center is determined
 automatically such that the line length is minimal.
 See also \code{points} and the examples below.
 
 Default: \code{NA}.
 }
 
 \item{points}{integer. If greater than 0,
 \code{points} gives the number of points simulated on the TBM
 line, hence 
 must be greater than the minimal number of points given by
 the size of the simulated field and the two paramters
 \code{linesimufactor} and \code{linesimustep}.
 If \code{points} is not positive the number of points is
 determined automatically.
 The use of \code{center} and \code{points} is highlighted
 in an example below.
 
 Default: \code{0}.
 }

 
 
 
 }

 \details{
\command{RPtbm}(Turning bands methods; turning layers).\cr
 It is generally difficult to use the turning bands method
 (\command{RPtbm}) directly
 in the 2-dimensional space.
 Instead, 2-dimensional random fields are frequently obtained
 by simulating a 3-dimensional random field (using
 \command{RPtbm}) and taking a 2-dimensional cross-section.
 3-dimensional \command{RPtbm} allows multiplicative models;
 in case of anisotropy the anisotropy matrices must be multiples
 of the first matrix or the anisotropy matrix consists of a time
 component only (i.e. all
 components are zero except the very last one).\cr
 \command{RPtbm} allows for arbitrary points, and
 arbitrary grids
 (arbitrary number of points in each direction, arbitrary grid length
 for each direction).

 \bold{Note:} Both the precision and the simulation time
 depend heavily on \code{linesimustep} and
 \code{linesimufactor}.
 For covariance models with larger values of the scale parameter,
 \code{linesimufactor=2} is too small.

 The turning layers are used for the simulations with time component.
 Here, 
 if the model is a
 multiplicative covariance function then the
 product may contain matrices with pure time component. All
 the other matrices must be equal up to a factor and the temporal
 part of the anisotropy matrix (right column) may contain only
 zeros, except the very last entry. 
 
}



\value{
 \code{RPtbm} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
 \itemize{

 \item Lantuejoul, C. (2002)
 \emph{Geostatistical Simulation: Models and Algorithms.}
 Springer.

 \item
 Matheron, G. (1973).
 The intrinsic random functions and their applications.
 \emph{Adv. Appl. Probab.}, \bold{5}, 439-468.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RPtbm(RMstable(s=1, alpha=1.8))
x <- seq(-3,3,0.1)
z <- RFsimulate(model=model, x=x, y=x)
plot(z)

model <- RPtbm(RMexp(Aniso=matrix(nc=2, rep(1,4))))
z <- RFsimulate(model=model, x=x, y=x)
plot(z)
\dontshow{FinalizeExample()}
}

\seealso{
  \command{\link{RMtbm}},
  \link{RP},
 \command{\link{RPspectral}}
}

\keyword{methods}
