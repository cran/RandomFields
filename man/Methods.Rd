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
    The way the method is implemented in this package
    it only allows for quadratic grids with
    common grid lengths.
    
    \item \code{direct matrix decomposition}.\cr
    This method is based on the well-known method for simulating
    any multivariate Gaussian distribution, using the square root of the
    covariance matrix.  The method is pretty slow and limited to
    about 1000 points, i.e. a 10x10x10 grid in three dimensions. 
    This implementation can use the Cholesky decomposition and
    the singular value decomposition. 
    It allows for arbitrary points and arbitrary grids.

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
      field by the central limit theorem -- this is the
      \code{add.mpp} method.
      
      \item \code{max.MPP} (Boolean functions).\cr
      If the random functions are multiplied by suitable,
      independent random values, and then the maximum is
      taken, a max-stable random field with unit Frechet margins
      is obtained -- this is the \code{max.mpp}
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

    \item \code{TBM2}, \code{TBM3} (Turning bands methods).\cr
    It is generally difficult to use the turning bands method
    (\code{TBM2}) directly
    in the 2-dimensional space.
    Instead, 2-dimensional random fields are frequently obtained
    by simulating a 3-dimensional random field (using
    \code{TBM3}) and taking a 2-dimensional cross-section.\cr
    \code{TBM2} and \code{TBM3} allow for arbitrary points, and
    arbitrary grids
    (arbitrary number of points in each direction, arbitrary grid length
    for each direction)
  }
}
\note{Most methods possess additional parameters,
  see \code{\link{RFparameters}()}
  that control the precision of the result.  The default parameters
  are chosen such that the simulations are fine for many models
  and their parameters. 
  The example in \code{\link{EmpiricalVariogram}()}
  shows a way of checking the precision.
}
\references{
  Gneiting, T. and Schlather, M. (2001)
  Statistical modeling with covariance functions.
  \emph{In preparation.}

  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields.}
  Technical report ST 99-10, Dept. of Maths and Statistics,
  Lancaster University. 

  Original work:
  \itemize{
    \item Circulant embedding:
    
    Chan, G. and Wood, A.T.A. (1997)
    An algorithm for simulating stationary Gaussian random fields.
    \emph{J. R. Stat. Soc., Ser. C} \bold{46}, 171--181.

    Dietrich, C.R. and Newsam, G.N. (1993)
    A fast and exact method for multidimensional Gaussian
    stochastic simulations.
    \emph{Water Resour. Res.} \bold{29}, 2861--2869.    
        
    Wood, A.T.A. and Chan, G. (1994)
    Simulation of stationary Gaussian processes in \eqn{[0,1]^d}{[0,1]^d}
    \emph{J. Comput. Graph. Stat.} \bold{3}, 409--432.

    
    \item Turning bands method (TBM):
    
    Dietrich, C.R. (1995) A simple and efficient space domain implementation
    of the turning bands method. \emph{Water Resour. Res.} \bold{31},
    147--156.
    
    Mantoglou, A. and Wilson, J.L. (1982) The turning bands method for
    simulation of random fields using line generation by a spectral
    method. \emph{Water. Resour. Res.} \bold{18}, 1379--1394.

    Matheron, G. (1973)
    The intrinsic random functions and their applications.
    \emph{Adv. Appl. Probab.} \bold{5}, 439--468.

    \item Random coins:
    
    Matheron, G. (1967) \emph{Elements pour une Theorie des Milieux
      Poreux}. Paris: Masson.
   }
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
\seealso{\code{\link{GaussRF}}, \code{\link{MaxStableRF}}, and
  \code{\link{RandomFields}}.}

\examples{
 PrintMethodList()
}
\keyword{spatial}



