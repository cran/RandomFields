\name{RFMethods}
\alias{RFMethods}
\alias{PrintMethodList}
\alias{GetMethodNames}
\title{Simulation Techniques}
\description{
  \code{PrintMethodList} prints the list of currently implemented methods for
  simulating random fields

  \code{GetMethodNames} returns a list of currently implemented methods
}
\usage{
PrintMethodList()

GetMethodNames()
}
\details{
  \itemize{
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
    hypermodel, see \command{\link{CovarianceFct}}, and \eqn{A=B} there.
    \code{cutoff embedding} halfens the maximum number of
    elements models used to define the covariance function of interest
    (from 10 to 5).

    Here multiplicative models are not allowed (yet).
  
     
    \item \code{direct matrix decomposition}.\cr
    This method is based on the well-known method for simulating
    any multivariate Gaussian distribution, using the square root of the
    covariance matrix.  The method is pretty slow and limited to
    about 1000 points, i.e. a 10x10x10 grid in three dimensions. 
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
    hypermodel, see \command{\link{CovarianceFct}}, and \eqn{A=B} there.

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
	Z(.) = sum_i  m_i( . - x_i)}
      In this package, only stationary Poisson point fields
      are allowed
      as underlying unmarked point processes.
      Thus, if the marks \eqn{m_i}{m_i}
      are all indicator functions, we obtain
      a Poisson random field. If the intensity of the Poisson
      process is high we obtain an approximate Gaussian random
      field by the central limit theorem - this is the
      \code{add.mpp} method.
      
      \item \code{max.MPP} (Boolean functions).\cr
      If the random functions are multiplied by suitable,
      independent random values, and then the maximum is
      taken, a max-stable random field with unit Frechet margins
      is obtained - this is the \code{max.mpp}
      method.
    }
    
    \item \code{nugget}.\cr
    One may specify this method (and "nugget" as covariance
    function) to generate a random field of
    independent Gaussian random variables.  However, any other
    method and any covariance function, called with zero
    variance, generates also such a random field (without loss
    of speed). 
    This method exists mainly for reasons of internal
    implementation.

    \item Random coins.\cr
    See marked point processes.
    
    \item \code{spectral TBM} (Spectral turning bands).\cr
    The principle of \code{spectral TBM}
    does not differ from the other
    turning bands methods.  However, line simulations are performed by a
    spectral technique (Mantoglou and Wilson, 1982); a
    realisation is given as the cosine with random
    amplitude and random phase.   
    The implementation allows the simulation of 2-dimensional random
    fields defined on arbitrary points or arbitrary grids.

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

    The turning layers are used for the simulations with time component.
    Here, 
    if the model is a
    multiplicative covariance function then the
    product may contain matrices with pure time component.  All
    the other matrices must be equal up to a factor and the temporal
    part of the anisotropy matrix (right column) may contain only
    zeros, except the very last entry.
    
  }
}
\note{Most methods possess additional parameters,
  see \command{\link{RFparameters}}\code{()}
  that control the precision of the result.  The default parameters
  are chosen such that the simulations are fine for many models
  and their parameters. 
  The example in \command{\link{EmpiricalVariogram}}\code{()}
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
        
    Wood, A.T.A. and Chan, G. (1994)
    Simulation of stationary Gaussian processes in \eqn{[0,1]^d}{[0,1]^d}
    \emph{J. Comput. Graph. Stat.} \bold{3}, 409-432.
    
    \item Intrinsic embedding and  Cutoff embedding:
    
    Stein, M.L. (2002)
    Fast and exact simulation of fractional Brownian surfaces.
    \emph{J. Comput. Graph. Statist.} \bold{11}, 587--599.
    
  
    Gneiting, T., Sevcikova, H., Percival, D.B., Schlather, M. and
    Jiang, Y. (2005)
    Fast and Exact Simulation of Large Gaussian Lattice Systems in
    \eqn{R^2}: Exploring the Limits
    \emph{J. Comput. Graph. Statist.} Submitted.


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
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}
  
  Yindeng Jiang \email{jiangyindeng@gmail.com} (circulant embedding
  methods \sQuote{cutoff} and \sQuote{intrinsic})
}
\seealso{\command{\link{GaussRF}}, \command{\link{MaxStableRF}},
  \command{\link{PrintModelList}},
  \code{\link{RandomFields}}.}

\examples{
 PrintMethodList()


############################################################
##                                                        ##
##                     Figure 1 in                        ##
##           Gneiting, T., Sevcikova, H.,                 ##
##   Percival, D.B., Schlather, M. and Jiang, Y. (2005)   ##
##                                                        ##
############################################################
## the example below shows for which parameter combinations
## of the stable model, the cutoff method works
## note the method 'cutoff' itself has an automatic
## selection algorithm implemented and is called with
## the stable model only.
\dontrun{

stabletest <- function(alpha, theta, size=512) {
  RFparameters(CE.trials=1, CE.tolIm = 1e-8, CE.tolRe=0, CE.force = FALSE,
               CE.useprimes=TRUE, CE.strategy=0, Print=2)
  model <- list(list(model="cutoff", var=1, kappa=c(1, theta, 1)),
                "(",
                list(model="stable", var=1, scale=1, kappa=alpha)
                )
  check <- CheckAndComplete(model=model, dim=2)
  if (check$error != 0) {
    cat("error", check$error, "\n")
    return(NA)
  }
  r <- paramextract(check$param[ , 1])$cutoffr
#  x <- c(r/size, r,  r /  size)
  x <- c(0, r, r / (size - 1)) * theta
  RFparameters(CE.mmin=rep(2 * size, 2))
  return(if (InitGaussRF(x, x, grid=TRUE, gridtriple=TRUE,
                         model=model, meth="circulant") == 0)
         r else NA)
}

alphas <- seq(1.52, 2.0, 0.01) # 0.02
thetas <- seq(0.05, 3.5, 0.05) #0.05( 2x)
m <- matrix(NA, nrow=length(thetas), ncol=length(alphas))
for (it in 1:length(thetas)) {
  theta <- thetas[it]
  for (ia in 1:length(alphas)) {
    alpha <- alphas[ia]
    cat("alpha=", alpha, "theta=", theta,"\n")
    print(unix.time(m[it, ia] <- stabletest(alpha=alpha, theta=theta)))
    if (is.na(m[it, ia])) break
  }
  if (any(is.finite(m))) image(thetas, alphas, m, col=rainbow(100))
}

} % dont run

}
\keyword{spatial}



