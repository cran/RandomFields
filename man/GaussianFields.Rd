\name{GaussianFields}
\alias{Gaussian}
\title{Methods for Gaussian Random Fields}
\description{
  Here, all the methods (models) for simulating
  Gaussian random fields are listed
}
\section{Implemented models}{
  \tabular{ll}{
 \command{\link{RPcirculant}} \tab simulation by circulant embedding \cr
 \command{\link{RPcutoff}} \tab simulation by a variant of circulant embedding \cr
 \command{\link{RPcoins}} \tab simulation by random coin / shot noise
 \cr
 \command{\link{RPdirect}} \tab through the square root of the covariance matrix\cr
 \command{\link{RPgauss}} \tab generic model that chooses automatically
 among the specific methods
 \cr
 \command{\link{RPhyperplane}} \tab simulation by hyperplane tessellation \cr
 \command{\link{RPintrinsic}} \tab simulation by a variant of circulant
 embedding \cr 
 \command{\link{RPnugget}} \tab simulation of (anisotropic) nugget effects \cr
 \command{\link{RPsequential}} \tab sequential method \cr
 \command{\link{RPspecific}} \tab model specific methods (very advanced)\cr
 \command{\link{RPspectral}} \tab spectral method \cr
 \command{\link{RPtbm}} \tab turning bands \cr
}
}


\section{Computing demand for simulations}{
  Assume at \eqn{n} locations in \eqn{d} dimensions a \eqn{v}-variate
  field has to be simulated.
  Let \deqn{f(n, d) = 2^d  n \log(n)}{f(n, d) = 2^d * n * log(n)}
  The following table gives in particular the time and memory needed for
  the specific simulation method.
\tabular{lllllll}{
  \tab grid \tab \eqn{v} \tab \eqn{d} \tab time \tab memory \tab comments\cr
  \command{\link{RPcirculant}}
  \tab yes \tab any \tab \eqn{\le 13}{<=13} \tab\eqn{O(v^3f(n, d))} \tab
  \eqn{O(v^2f(n, d))}\tab \cr
  \tab no \tab any \tab \eqn{\le 13}{<=13} \tab \eqn{O(v^3 f(k, d))} \tab \eqn{O(v^2f(k, d))}
  \tab \eqn{k \sim }{k ~ }\command{\link[=RFoptions]{approx_step}}\eqn{{}^{-d}}{^{-d}}\cr
 \command{\link{RPcutoff}} \tab \tab \tab \tab \tab \tab see
  RPcirculant above
 \cr
 \command{\link{RPcoins}} \tab yes \tab \eqn{1} \tab \eqn{\le 4}{<=4} \tab \eqn{O(k
  n)}{O(k * n)} \tab\eqn{ O(n) }\tab
 \eqn{k \sim}{k ~ }\eqn{(lattice spacing)^{-d}} \cr
 \tab no \tab \eqn{1} \tab \eqn{\le 4}{<=4} \tab \eqn{O(k n)}{O(k * n)} \tab
  \eqn{O(n)} \tab \eqn{k} depends on the geometry
 \cr
 \command{\link{RPdirect}}
 \tab any\tab any\tab any \tab\eqn{O(1)..O(v^2 n^2)}{ O(v^2 * n^2)}\tab\eqn{ O(v^2
   n^2)}{ O(v^2 * n^2)}\tab effort to investigate the covariance matrix, if
 \code{\link[=RFoptions]{matrix_methods}} is not specified (default)\cr
 \tab \tab \tab  \tab\eqn{O(v n)}{ O(v * n)}\tab\eqn{ O(v
  n)}{ O(v * n)}\tab covariance matrix is diagonal \cr
 \tab \tab \tab  \tab see \pkg{\link[spam:SPAM]{spam}} \tab \eqn{O(z + v
  n)}{O(z + v * n)}\tab covariance matrix is sparse
 matrix with \eqn{z} non-zeros\cr
 \tab\tab \tab  \tab \eqn{O(v^3 n^3)}{O(v^3 * n^3)}\tab
  \eqn{O(v^2 n^2)}{O(v^2*n^2)}\tab arbitrary covariance matrix (preparation) \cr
  \tab\tab \tab  \tab \eqn{O(v^2 n^2)}{O(v^2*n^2)}\tab
  \eqn{O(v^2 n^2)}{O(v^2*n^2)}\tab arbitrary covariance matrix (simulation)\cr
 \command{\link{RPgauss}} \tab any \tab any \tab any \tab \eqn{O(1)
  \ldots O(v^3n^3)}{O(1)..O(v^3*n^3)} \tab
 \eqn{O(1)\ldots O(n^2)}{O(1)..O(n^2)}\tab \bold{only} the selection
 process; \eqn{O(1)} if first method tried is successful
 \cr
 \command{\link{RPhyperplane}} \tab any \tab \eqn{1} \tab \eqn{2} \tab
 \eqn{O(n / s^d)} \tab \eqn{O(n / s^d)}\tab
 \eqn{s = }\code{\link[=RMmodels]{scale}}\cr
 \command{\link{RPintrinsic}} \tab  \tab \tab \tab \tab \tab see
  RPcirculant above\cr 
 \command{\link{RPnugget}} \tab any  \tab any\tab any\tab \eqn{O(v n)}
  \tab \eqn{O(v n)} \tab\cr
 \command{\link{RPsequential}} \tab any\tab \eqn{1} \tab any \tab
 \eqn{O(S^3 b^3)}{O(S^3 * b^3)}
 \tab
 \eqn{O(S^2 b^2)}{O(S^2*b^2)}
 \tab
 \eqn{n=ST}{n = S * T};
 \eqn{S} and \eqn{T} the number of spatial and temporal locations,
  respectively; \eqn{b = }\code{\link[=RPsequential]{back_steps}} (preparation) \cr
 \tab \tab  \tab  \tab
 \eqn{O(n S b^2)}{O(n * S * b^2)}
 \tab
 \eqn{O(S^2 b^2) + O(n)}{O(n)} 
 \tab (simulation) \cr
 \command{\link{RPspectral}} \tab any \tab \eqn{1} \tab \eqn{\le 2}{<=2}
  \tab\eqn{O(C(d) n)}{O(C(d) * n)}
 \tab \eqn{ O(n) }\tab  \eqn{C(d)} : large constant increasing in \eqn{d} \cr
 \command{\link{RPtbm}} \tab any \tab \eqn{1} \tab \eqn{\le 4}{<=4} \tab
  \eqn{ O(C(d) (n + L) }{ O(C(d) * (n + L))} \tab\eqn{ O(n + L)}  \tab\eqn{C(d)} :
  large constant increasing in \eqn{d}; \eqn{L} is the effort needed to
  simulate on a line (or plane)\cr
 \command{\link{RPspecific}}  \tab \tab \tab \tab \tab \tab \bold{only}
  the specific part \cr
 * * \command{\link{RMplus}}\tab any \tab any \tab any \tab O(v n) \tab
  O(v n) \tab\cr
 * * \command{\link{RMS}}\tab any \tab any \tab any  \tab O(1) \tab O(v n) \tab\cr
 * * \command{\link{RMmult}}\tab any \tab any \tab any  \tab O(v n) \tab
  O(v n) \tab\cr
}
}



\section{Computing demand for interpolation}{
  Assume \eqn{v}-variate data are given at 
  \eqn{n} locations in \eqn{d} dimensions.
  To interpolate at \eqn{k} locations RandomFields needs
\tabular{lllllll}{
  grid \tab \eqn{v} \tab \eqn{d} \tab time \tab memory \tab comments\cr
  any\tab any\tab any \tab \eqn{O(1)..O(v^2 n^2)}{ O(v^2 * n^2)}\tab\eqn{ O(v^2
   n^2)}{ O(v^2 * n^2)}\tab effort to investigate the covariance matrix, if
 \code{\link[=RFoptions]{matrix_methods}} is not specified (default)
 \cr
  \tab \tab  \tab\eqn{O(v ^2 n k)}{ O(v^2 * n k)}\tab\eqn{ O(v
    (n + k))}{ O(v * (n + k))}\tab covariance matrix is diagonal
  \cr
  \tab \tab  \tab see \pkg{\link[spam:SPAM]{spam}}+ O(v^2nk) \tab \eqn{O(z + v
  (n + k))}{O(z + v * (n + k))}\tab covariance matrix is sparse
 matrix with \eqn{z} non-zeros\cr
 \tab \tab  \tab \eqn{O(v^3 n^3 + v^2nk)}{O(v^3*n^3 + v^2*n*k)}\tab
 \eqn{O(v^2 n^2 + v*k)}{O(v^2*n^2 + v*k)}\tab arbitrary covariance matrix
}
}



\section{Computing demand for conditional simulation}{
  Assume \eqn{v}-variate data are given at 
  \eqn{n} locations \eqn{x_1,\ldots, x_n}{x_1,...,x_n} in \eqn{d} dimensions.
  To conditionally simulated at \eqn{k} location
  \eqn{y_1,\ldots, y_k}{y_1,...,y_k}, the
  computing demand equals the
  sum of the demand for interpolating and
  the demand for simulating on the \eqn{k+n} location.
  (Grid algorithms for simulating will apply if the \eqn{k} locations
  \eqn{y_1,\ldots, y_k}{y_1,...,y_k}
  are defined by a grid and the \eqn{n} locations
  \eqn{x_1,\ldots, x_n}{x_1,...,x_n} are a subset of
  \eqn{y_1,\ldots, y_k}{y_1,...,y_k}, a situation typical in image analysis.)
}
 


\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 % \item Gneiting, T. and Schlather, M. (2004)
 % Statistical modeling with covariance functions.
 % \emph{In preparation.}
 \item Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University.
 
 \item Schlather, M. (2010)
 On some covariance models based on normal scale mixtures.
 \emph{Bernoulli}, \bold{16}, 780-797.
 
 \item Schlather, M. (2011) Construction of covariance functions and
 unconditional simulation of random fields. In Porcu, E., Montero, J.M.
 and Schlather, M., \emph{Space-Time Processes and Challenges Related
 to Environmental Problems.} New York: Springer.
 % \item Schlather, M. (2002) Models for stationary max-stable
 % random fields. \emph{Extremes} \bold{5}, 33-44.
 \item Yaglom, A.M. (1987) \emph{Correlation Theory of Stationary and
 Related Random Functions I, Basic Results.}
 New York: Springer.
 \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 3nd edition.
 }
}
\seealso{
   \link{RP},
  \command{\link{Other models}},
  \command{\link{RMmodel}},
  \command{\link{RFgetMethodNames}}
  \command{\link{RFsimulateAdvanced}}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
set.seed(1)
x <- runif(90, 0, 500)
z <- RFsimulate(RMspheric(), x)
z <- RFsimulate(RMspheric(), x, max_variab=10000)
\dontshow{FinalizeExample()}

}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\keyword{spatial}

