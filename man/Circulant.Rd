\name{Circulant Embedding}
\alias{Circulant}
\alias{RPcirculant}
\alias{Cutoff}
\alias{RPcutoff}
\alias{Intrinsic} 
\alias{RPintrinsic} 
\title{Circulant Embedding methods}
\description{

 Circulant embedding is a fast simulation method for stationary
 (possibly anisotropic) Gaussian random
 fields on regular grids based on Fourier transformations. It is garantueed to be an exact method
 for covariance functions with finite support, e.g. the spherical
 model. The method is admissible for any dimension apart from memory
 restictions. \cr
 The simulation is performed on a torus which represents the bended
 grid. To remove wrong dependencies occuring at different borders of the
 grid which would be close on the torus, the simulation area is
 multiplied by a natural number.
 There is also a multivariate version of circulant embedding.

 Cut-off embedding is a fast simulation method for stationary,
 isotropic Gaussian random fields on square lattices based on
 the standard \command{\link{RPcirculant}} method,
 so that exact simulation is garantueed
 for further covariance models, e.g. the \command{\link{RMwhittle}} model.

 In fact, the circulant embedding is called with the cutoff
 hypermodel, see \command{\link{RMcutoff}}.
 \code{Cutoff} halfens the maximum number of
 elements models used to define the covariance function of interest
 (from 10 to 5).

 Here multiplicative models are not allowed (yet).
 \cr
 For details see \command{\link{RMcutoff}}.


 Intrinsic embedding is a fast simulation method for intrinsically stationary,
 isotropic Gaussian random fields on square lattices based on
 the standard \command{\link{RPcirculant}} method,
 for further \emph{variogram} models, e.g. \command{\link{RMfbm}}.

 Note that the simulated random field is always \emph{non-stationary}.
 In fact, the circulant embedding is called with the Intrinsic
 hypermodel, see \command{\link{RMintrinsic}}.

 Here multiplicative models are not allowed (yet).
 \cr
 For details see \command{\link{RMintrinsic}}.

 
}

\usage{
RPcirculant(phi, loggauss, force, mmin, strategy,
 maxmem, tolIm, tolRe, trials, useprimes, dependent,
 approx_step, approx_maxgrid)

RPcutoff(phi, loggauss,  force, mmin, strategy,
 maxmem, tolIm, tolRe, trials, useprimes, dependent,
 approx_step, approx_maxgrid, diameter, a) 
 

RPintrinsic(phi, loggauss, force, mmin, strategy,
 maxmem, tolIm, tolRe, trials, useprimes, dependent,
 approx_step, approx_maxgrid, diameter, rawR) 
}

\arguments{
 \item{phi}{See \command{\link{RPgauss}}}
 \item{loggauss}{See \command{\link{RPgauss}}}
 \item{force}{
   Logical. Circulant embedding does not work
   if the constructed circulant
   matrix has negative eigenvalues. Sometimes it is convenient
   to replace all the negative eigenvalues by zero
   (\code{force=TRUE}) after \code{trials} number of trials. 
   Default: \code{FALSE}
 }
 \item{mmin}{Scalar or vector, integer if positive.
   \code{CE.mmin} determines the initial size of the circulant
   matrix. If \code{CE.mmin=0} the minimal starting size is
   determined automatically according to the
   dimensions of the grid.
   If \code{CE.mmin>0} then the absolute starting size is given.
   If \code{CE.mmin<0} then the automatically determined
   matrix size is multiplied
   by \eqn{|\code{CE.mmin}|}; here \code{CE.mmin} must be smaller
   than -1; the 
   value -1 takes over the minimal starting size.\cr
   Note: in any cases, the initial size might be increased according
   to \code{CE.useprimes}.\cr
   Default: \code{0}
 } 
 \item{strategy}{Logical. \code{0}, if the circulant
 matrix has negative eigenvalues then the
 size in each direction is doubled; \cr \code{1}: the size is
 enhanced only in 
 one direction, namely that one where the covariance function has the
 largest value at the end point of the grid --- note that
 the default value of \code{trials} is probably too small
 in that case.

 In some cases \code{strategy=0} works better, in other cases
 \code{strategy=1}. Just try. Clearly, if the field is
 isotropic and a square grid should be simulated, then
 \code{strategy=0} is the better choice.
 
 Default: \code{0}
 }
 
 \item{maxmem}{Integer.
   maximal number of entries in a row of the circulant matrix. 
   The total amount of memory needed for the internal calculations
   is about 16 (=2 * sizeof(double))
   times as large as \code{maxmem}
   if \code{\link{RFoptions}}\code{()$Storing=FALSE},
   and 32 (=4 * sizeof(double)) time as large if \code{Storing=TRUE}.
   
   The value of \code{maxmem} must be at least twice as large as
   the number of points to be simulated. In some case even much larger.
   
   Note that \code{maxmem} can be used to control the automatic
   choice of the simulation algorithm. Namely, in case of huge
   circulant matrices, other simulation
   methods (TBM) are faster and might be preferred by the user.
   
   Default: \code{4096^2 = 16777216}
 }
 
 \item{tolIm}{
   If the modulus of the imaginary part is less than
   \code{tolIm} then the eigenvalue is always considered as
   real (independently of the value of \code{force}).
   
   Default: \code{1E-3}
 }
 
 \item{tolRe}{
   Eigenvalues between \code{tolRe} and 0 are always considered as
   0 and set 0 (independently of the value of \code{force}).

   Default: \code{-1E-7}
 }

 \item{trials}{Integer.
   A larger circulant matrix is likely to make more eigenvalues
   non-negative. If at least one of the thresholds \code{tolRe} and
   \code{tolIm} are missed then the matrix size is doubled
   according to \code{strategy},
   and the matrix is checked again. This procedure is repeated
   up to \code{trials - 1} times. If there are still negative
   eigenvalues, the simulation method fails if \code{force=FALSE}. 
   
   Default: \code{3}
 }
 
 \item{useprimes}{Logical. If \code{FALSE}
   the columns of the circulant matrix
   have length \eqn{2^k} for some \eqn{k}. Otherwise the algorithm
   tries to find a nicely factorizable number close to the size of the
   given matrix.
   
   Default: \code{TRUE}
 }
 \item{dependent}{Logical. If \code{FALSE}
 then independent random fields are created. If \code{TRUE}
 then at least 4 non-overlapping rectangles are taken out of the
 the expanded grid defined by the circulant matrix.
 These simulations are dependent.
 See below for an example.
 See \code{trials} for some
 more information on the circulant matrix.

 Default: \code{FALSE}
 }
 \item{approx_step}{Real value.
   It gives the grid size of the approximating grid in case
   circulant embedding is used although the points do not ly on a grid.
   
   If \code{NA} then \code{approx_step} is chosen such that
   \code{approx_maxgrid} is nearly reached.   
 }
 \item{approx_maxgrid}{
   It defaults to \code{maxmem}.
 }

 \item{diameter}{See \command{\link{RMcutoff}} or \command{\link{RMintrinsic}}}
 \item{a}{See \command{\link{RMcutoff}}}
 \item{rawR}{See \command{\link{RMintrinsic}}}
 }
 

 

\details{
 Here, the algorithms by Dietrich and Newsam are implemented.
}



\value{
 an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
 Circulant Embedding
 \itemize{

 \item Chan, G. and Wood, A.T.A. (1997)
 An Algorithm for Simulating Stationary Gaussian Random Fields.
 \emph{Journal of the Royal Statistical Society: Series C (Applied
 Statistics)}, \bold{46}, 171--181.

 \item Dietrich, C. R. and G. N. Newsam (1993)
 A fast and exact method for multidimensional gaussian
 stochastic simulations.
 \emph{Water Resour. Res.} \bold{29(8)}, 2861--2869.


 \item Dietrich, C. R. and G. N. Newsam (1996)
 A fast and exact method for multidimensional Gaussian stochastic
 simulations: Extension to realizations conditioned on direct and
 indirect measurements 
 \emph{Water Resour. Res.} \bold{32(6)}, 1643--1652.

 
 \item Dietrich, C. R. and Newsam, G. N. (1997)
 Fast and Exact Simulation of Stationary Gaussian Processes through
 Circulant Embedding of the Covariance Matrix.
 \emph{SIAM J. Sci. Comput.} \bold{18}, 1088--1107.

 \item Wood, A. T. A. and Chan, G. (1994)
 Simulation of Stationary Gaussian Processes in \eqn{[0, 1]^d}.
 \emph{Journal of Computational and Graphical Statistics}
 \bold{3}, 409--432. 
 
 }

 Cutoff and Intrinsic
 \itemize{
 \item Gneiting, T., Sevecikova, H, Percival, D.B., Schlather M.,
 Jiang Y. (2006) Fast and Exact Simulation of Large {G}aussian
 Lattice Systems in {$R^2$}: Exploring the Limits.
 \emph{J. Comput. Graph. Stat.} \bold{15}, 483--501.
 \item Stein, M.L. (2002) Fast and exact simulation of fractional
 Brownian surfaces. \emph{J. Comput. Graph. Statist.} \bold{11}, 587--599 
 }

}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\seealso{
  \link{RP},  \command{\link{RPtbm}}
}

\examples{
RFoptions(seed=0)

model <- RMstable(s=1, alpha=1.8)
x <- seq(-3,3,0.1)
z <- RFsimulate(model=RPcirculant(model), x=x, y=x, grid=TRUE, n=1)
plot(z)

model <- RMexp(var=10, s=10)
z <- RFsimulate(model=RPcirculant(model), 1:10)
plot(z)

model <- RMfbm(Aniso=diag(c(1,2)), alpha=1.5)
z <- RFsimulate(model, x=1:10, y=1:10, grid=TRUE)
plot(z)

\dontshow{RFoptions(seed=NA)}

}

\keyword{methods}
