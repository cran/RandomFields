\name{RandomFields}
\alias{RandomFields}
\title{Simulation and Analysis of Random Fields}
\description{
  The package \code{RandomFields} allows for simulating various kinds
  of random fields.  The emphasis is (still) on unconditional simulation
  of stationary and isotropic
  Gaussian random fields. Furthermore, algorithms for  
  conditional simulation and simulation of
  max-stable random fields are provided.
  
  Additionally, the package includes tools for analysing spatial data:
  empirical variogram, interactive fitting, and MLE estimation of
  parameters.
  Basic kriging procedures are also provided.

  Version 2.0 will contain various extensions of the packages
  (anisotropic random fields, user defined nested
  models and multiplicative models). It is is supposed to be
  released the latest by the end of 2002.
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
         fields; nice examples to get familiar with the
         simulation features of the package
      \item \code{\link{Kriging}} : simple and ordinary kriging
      \item \code{\link{mleRF}} : maximum likelihood estimator for
      random field parameters
      \item \code{\link{RFparameters}} : control parameters (advanced settings)
      \item \code{\link{PrintMethodList}} : list of implemented
      simulation methods
      \item \code{\link{ShowModels}} : interactive, graphical choice of
      models
      \item \code{\link{soil}} : Soil physical and chemical data;
      the \code{example} gives a simple geostatistical analysis using
      features of the package      
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
\references{
  Schlather, M. (2001) Simulation of stationary and isotropic
    random fields. \emph{R News} \bold{1} (2), 18-20.
  
  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields.}
  Technical report ST 99-10, Dept. of Maths and Statistics,
  Lancaster University. 
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
\section{Acknowledgement}{
  Many thanks to Martin Maechler, Paulo Ribeiro, and Tilmann Gneiting
  for proof-reading parts of the code and the help text of this
  package. 

  The work has been supported by the EU TMR network ERB-FMRX-CT96-0095 on
  ``Computational and statistical methods for the analysis of spatial
  data'' in 1999,
  and by the German Federal Ministry of Research and Technology 
  (BMFT) grant PT BEO 51-0339476C during 2000--2001.
}
\keyword{spatial}


