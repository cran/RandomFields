\name{GaussianFields}
\alias{Gaussian}
\title{Methods for Gaussian Random Fields}
\description{
  Hier all the methods (models) for simulating
  Gaussian random fields are listed
}
\section{Implemented models}{
  \tabular{ll}{
 \command{\link{RPcirculant}} \tab simulation by circulant embedding \cr
 \command{\link{RPcutoff}} \tab simulation by a variant of circulant embedding \cr
 \command{\link{RPcoins}} \tab simulation by random coin / shot noise
 \cr
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
  \command{\link{RFsimulateAdvanced}}
}
\examples{
RFoptions(seed=0)

RFoptions(printlevel=4)
x <- seq(0, 10, 0.01)
z <- RFsimulate(RMexp(), x) # circulant embedding is chosen
\dontshow{RFoptions(seed=NA)}

}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}

\keyword{spatial}

