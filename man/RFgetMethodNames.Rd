\name{RFgetMethodNames}
\alias{RFgetMethodNames}
\title{Simulation Techniques}
\description{
 \command{\link{RFgetMethodNames}} prints and returns a list of currently
 implemented methods for simulating Gaussian random fields
}
\usage{
RFgetMethodNames(show=TRUE)
}
\arguments{
  \item{show}{logical. Whether the methods are shown.
  }
}
\value{
  an invisible string vector of the Gaussian methods.
}
\details{
 By default, \command{\link{RFsimulate}} automatically chooses an appropriate
 method for simulation. The method can also be set explicitly by the
 user via \command{\link{RFoptions}}, in particular by passing
 \code{gauss.method=_a valid method string_} as an additional argument
 to \command{\link{RFsimulate}} or by globally changing the options via
 \code{\link{RFoptions}(gauss.method=_a valid method string_)}.
 The following methods are available:
 \itemize{
 \item (random spatial) Averages\cr
 -- details soon
 
 \item Boolean functions.\cr
 See marked point processes.

 \item \code{circulant embedding}. \cr
 Introduced by Dietrich & Newsam (1993) and Wood and Chan (1994). 

 Circulant embedding is a fast simulation method based on
 Fourier transformations. It is garantueed to be an exact method
 for covariance functions with finite support, e.g. the spherical
 model.

 See also \code{cutoff embedding} and \code{intrinsic embedding} for
 variants of the method.
 
 
 \item \code{cutoff embedding}. \cr
 Modified circulant embedding method so that exact simulation is garantueed
 for further covariance models, e.g. the whittle matern model.
 In fact, the circulant embedding is called with the cutoff
 hypermodel, see \command{\link{RMmodel}}, and \eqn{A=B} there.
 \code{cutoff embedding} halfens the maximum number of
 elements models used to define the covariance function of interest
 (from 10 to 5).

 Here multiplicative models are not allowed (yet).
 
 
 \item \code{direct matrix decomposition}.\cr
 This method is based on the well-known method for simulating
 any multivariate Gaussian distribution, using the square root of the
 covariance matrix. The method is pretty slow and limited to
 about 8000 points, i.e. a 20x20x20 grid in three dimensions. 
 This implementation can use the Cholesky decomposition and
 the singular value decomposition. 
 It allows for arbitrary points and arbitrary grids.

 \item \code{hyperplane method}.\cr
 The method is based on a tessellation of the space by
 hyperplanes. Each cell takes a spatially constant value
 of an i.i.d. random variables. The superposition of several
 such random fields yields approximatively a Gaussian random field.
 
 \item \code{intrinsic embedding}. \cr
 Modified circulant embedding so that exact simulation is garantueed
 for further \emph{variogram} models, e.g. the fractal brownian one.
 Note that the simulated random field is always \emph{non-stationary}.
 In fact, the circulant embedding is called with the Stein
 hypermodel, see \command{\link{RMmodel}}, and \eqn{A=B} there.

 Here multiplicative models are not allowed (yet).

 \item Marked point processes.\cr
 Some methods are based on marked point process
 \eqn{\Pi=\bigcup [x_i,m_i]}{P = ([x_1,m_1], [x_2,m_2], ...)}
 where the marks \eqn{m_i}{m_i}
 are deterministic or i.i.d. random functions on \eqn{R^d}{R^d}.
 \itemize{
 \item \code{add.MPP} (Random coins).\cr
 Here the functions are elements
 of the intersection \eqn{L_1 \cap L_2}{(L1 cap L2)}
 of the Hilbert spaces \eqn{L_1}{L1} and \eqn{L_2}{L2}.
 A random field Z is obtained by adding the marks:
 \deqn{ Z(\cdot) = \sum_{[x_i,m_i] \in \Pi} m_i(\cdot - x_i)}{
	Z(.) = sum_i m_i( . - x_i)}
 In this package, only stationary Poisson point fields
 are allowed
 as underlying unmarked point processes.
 Thus, if the marks \eqn{m_i}{m_i}
 are all indicator functions, we obtain
 a Poisson random field. If the intensity of the Poisson
 process is high we obtain an approximative Gaussian random
 field by the central limit theorem - this is the
 \code{add.mpp} method.
 
 \item \code{max.MPP} (Boolean functions).\cr
 If the random functions are multiplied by suitable,
 independent random values, and then the maximum is
 taken, a max-stable random field with unit Frechet margins
 is obtained - this is the \code{max.mpp}
 method.
 }
 
% \item \code{Markov} Gaussian Random Fields.\cr
% The method is due to Havard Rue and uses the property that the
% precision matrix (i.e. the inverse of the covariance matrix)
% of a Markov Gaussian random field has a small band width.
%
% \bold{Note:} This method is only available on Intel Linux systems
% when installed as follows:\cr
% 
% 1. install GMRFLib by Havard Rue from \url{www.math.ntnu.no/~hrue/GMRFsim}
% 
% 2. install all external libraries given in GMRFLib/EXTLIBS/
% (it may happen that precompiled libraries included in your distribution
% work better)
%
% 3. install RandomFields with option --enable-GMRF and 
% --with-GMRF-lib=LIB_PATH and --with-GMRF-EXT-lib=LIB_PATH 
% for the path ot GMRFLib and the EXTLIBS, respectively.
% E.g., \sQuote{R CMD INSTALL RandomFields_* --configure-args="--enable-GMRF --with-GMRF-EXT-lib=/home/schlather/local/lib"}

 \item \code{nugget}.\cr
 The method allows for generating a random field of 
 independent Gaussian random variables. In the isotropic case
 and if the simple notation of a model (with \code{model} and
 \code{param})
 is used, this method is called automatically if the nugget
 effect is positive except the method \code{"circulant embedding"}
 or \code{"direct"} have been explicitely.

 The method has been extended to zonal anisotropies, see
 also argument \code{nugget.tol} in \command{\link{RFoptions}}.

 \item \code{particular} method\cr
 -- details missing --
 
 \item Random coins.\cr
 See marked point processes.

 \item \code{sequential}
 This method is programmed for spatio-temporal models
 where the field is modelled sequentially in the time direction
 conditioned on the previous \eqn{k} instances.
 For \eqn{k=5} the method has its limits for about 1000 spatial
 points. It is an approximative method. The larger \eqn{k} the
 better.
 It also works for certain grids where the last dimension should
 contain the highest number of grid points.
 
 \item \code{spectral TBM} (Spectral turning bands).\cr
 The principle of \code{spectral TBM}
 does not differ from the other
 turning bands methods. However, line simulations are performed by a
 spectral technique (Mantoglou and Wilson, 1982).

 The standard method
 allows for the simulation of 2-dimensional random
 fields defined on arbitrary points or arbitrary grids. Here
 realisation is given as the cosine with random
 amplitude and random phase. 

 

 \item \code{TBM2}, \code{TBM3} (Turning bands methods; turning layers).\cr
 It is generally difficult to use the turning bands method
 (\code{TBM2}) directly
 in the 2-dimensional space.
 Instead, 2-dimensional random fields are frequently obtained
 by simulating a 3-dimensional random field (using
 \code{TBM3}) and taking a 2-dimensional cross-section.
 TBM3 allows for multiplicative models; in case of anisotropy the
 anisotropy matrices must be multiples of the first matrix or the
 anisotropy matrix consists of a time component only (i.e. all
 components are zero except the very last one).\cr
 \code{TBM2} and \code{TBM3} allow for arbitrary points, and
 arbitrary grids
 (arbitrary number of points in each direction, arbitrary grid length
 for each direction).

 \bold{Note:} Both the precision and the simulation time
 depend heavily on \code{TBM*.linesimustep} and
 \code{TBM*.linesimufactor}
 that can be set by \command{\link{RFoptions}}.
 For covariance models with larger values of the scale parameter,
 \code{TBM*.linesimufactor=2} is too small.

 The turning layers are used for the simulations with time component.
 Here, 
 if the model is a
 multiplicative covariance function then the
 product may contain matrices with pure time component. All
 the other matrices must be equal up to a factor and the temporal
 part of the anisotropy matrix (right column) may contain only
 zeros, except the very last entry. 
 }
}
\note{Most methods possess additional arguments,
 see \code{\link{RFoptions}()}
 that control the precision of the result. The default arguments 
 are chosen such that the simulations are fine for many models
 and their parameters. 
 The example in \code{\link{RFempiricalvariogram}()}
 shows a way of checking the precision.
}
\references{
 Gneiting, T. and Schlather, M. (2004)
 Statistical modeling with covariance functions.
 \emph{In preparation.}

 Lantuejoul, Ch. (2002) \emph{Geostatistical simulation.}
 \bold{New York:} Springer.
 
 Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University. 

 Original work:
 \itemize{
 \item Circulant embedding:
 
 Chan, G. and Wood, A.T.A. (1997)
 An algorithm for simulating stationary Gaussian random fields.
 \emph{J. R. Stat. Soc., Ser. C} \bold{46}, 171-181.
 
 Dietrich, C.R. and Newsam, G.N. (1993)
 A fast and exact method for multidimensional Gaussian
 stochastic simulations.
 \emph{Water Resour. Res.} \bold{29}, 2861-2869. 
 
 Dietrich, C.R. and Newsam, G.N. (1996)
 A fast and exact method for multidimensional {G}aussian stochastic
 simulations: Extensions to realizations conditioned on direct and
 indirect measurement 
 \emph{Water Resour. Res.} \bold{32}, 1643-1652. 

 Wood, A.T.A. and Chan, G. (1994)
 Simulation of stationary Gaussian processes in \eqn{[0,1]^d}{[0,1]^d}
 \emph{J. Comput. Graph. Stat.} \bold{3}, 409-432.

 The code used in \cite{RandomFields} is based on
 Dietrich and Newsam (1996).
 
 \item Intrinsic embedding and Cutoff embedding:
 
 Stein, M.L. (2002)
 Fast and exact simulation of fractional Brownian surfaces.
 \emph{J. Comput. Graph. Statist.} \bold{11}, 587--599.
 
 
 Gneiting, T., Sevcikova, H., Percival, D.B., Schlather, M. and
 Jiang, Y. (2005)
 Fast and Exact Simulation of Large Gaussian Lattice Systems in
 \eqn{R^2}: Exploring the Limits
 \emph{J. Comput. Graph. Statist.} Submitted.

 \item Markov Gaussian Random Field:
 
 Rue, H. (2001) Fast sampling of Gaussian Markov random fields.
 \emph{J. R. Statist. Soc., Ser. B}, \bold{63} (2), 325-338.
 
 Rue, H., Held, L. (2005) \emph{Gaussian Markov Random Fields:
 Theory and Applications.}
 Monographs on Statistics and Applied Probability, no \bold{104},
 Chapman \& Hall.


 \item Turning bands method (TBM), turning layers:
 
 Dietrich, C.R. (1995) A simple and efficient space domain implementation
 of the turning bands method. \emph{Water Resour. Res.} \bold{31},
 147-156.
 
 Mantoglou, A. and Wilson, J.L. (1982) The turning bands method for
 simulation of random fields using line generation by a spectral
 method. \emph{Water. Resour. Res.} \bold{18}, 1379-1394.

 Matheron, G. (1973)
 The intrinsic random functions and their applications.
 \emph{Adv. Appl. Probab.} \bold{5}, 439-468.

 Schlather, M. (2004)
 Turning layers: A space-time extension of turning bands.
 \emph{Submitted}

 \item Random coins:
 
 Matheron, G. (1967) \emph{Elements pour une Theorie des Milieux
 Poreux}. Paris: Masson.
 }
 }
 \section{Automatic selection algorithm}{
 --- details coming soon ---
 }
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \link[=RandomFields-package]{RandomFields}.
}
\keyword{spatial}




\examples{
RFgetMethodNames()
\dontshow{FinalizeExample()}
}


