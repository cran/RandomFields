\name{RandomFields}
\alias{RandomFields}
\title{Simulation and Analysis of Random Fields}
\description{
  The package \code{RandomFields} allows for simulating various kinds
  of random fields.  The emphasis is (still) on unconditional simulation
  of stationary and isotropic
  Gaussian random fields, but the package provides also basics for
  conditional simulation of Gaussian random fields.  Furthermore,
  non-Gaussian random fields, i.e.
  max-stable random fields, can also be simulated.
  
  Additionally, the package provides tools for analysing spatial data,
  namely the calculation of the empirical variogram, interactive
  fitting of a variogram model, and the MLE estimation of variogram
  parameters. A very few, basic, kriging procedures are also provided.

  Diverse extensions of the package are planned.
}
\details{
  The following random fields and related functionallities are
  provided by the package.
  \enumerate{
    \item stationary and isotropic Gaussian random fields\cr
    \itemize{
      \item \code{\link{CondSimu}} : conditional simulation
      \item \code{\link{CovarianceFct}} :  covariance functions or
      variogram models
      \item \code{\link{EmpiricalVariogram}} : empirical variogram
      \item \code{\link{GaussRF}} : simulation of Gaussian random
      fields
      \item \code{\link{Kriging}} : simple and ordinary kriging
      \item \code{\link{mleRF}} : maximum likelihood estimator for
      random field parameters
      \item \code{\link{RFparameters}} : control parameters (advanced settings)
      \item \code{\link{PrintMethodList}} : list of implemented
      simulation methods
      \item \code{\link{ShowModels}} : interactive, graphical choice of
      models
    }
    \item stationary and isotropic max-stable random fields\cr
    \itemize{
     \item \code{\link{CovarianceFct}} : covariance models for
     extremal Gaussian random fields
     \item \code{\link{MaxStableRF}} : simulation of max-stable
     random fields
     \item \code{\link{RFparameters}} : control parameters (advanced)     
    }
    
    \item stationary and isotropic Poisson random fields
    (convenient access not implemented yet)\cr
%    \itemize{
%      \item \code{%\link{
%	PoissonRF%}
%      } : simulation of Poisson random fields
%    }
  }

  Functions used for diverse simulation methods:
  \itemize{
    \item \code{\link{DeleteRegister}} : deleting internal registers
  }
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
\section{Acknowledgement}{
  Thanks to Martin Maechler, Paulo Ribeiro, and Tilmann Gneiting
  for proof-reading parts of the code and the help text of this
  package. 
  
}
\keyword{spatial}


