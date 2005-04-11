\name{RandomFields}
\alias{RandomFields}
\title{Simulation and Analysis of Random Fields}
\description{
  The package \code{RandomFields} allows for simulating various kinds
  of random fields, including anisotropic
  processes. Furthermore, algorithms for  
  conditional simulation and simulation of
  max-stable random fields are provided.
  
  Additionally, the package includes tools for analysing spatial data:
  Hurst parameter, fractal dimension, empirical variogram,
  interactive fitting of parameters, LSQ and MLE estimation of parameters.
  Basic kriging procedures are also provided.

  Further extensions of the package are planned.
}
\details{
  The following random fields and related functionalities are
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
      \item \code{\link{fitvario}} : variogram/covariance function fit
      by least squares, maximum likelihood and cross validation
      techniques
      \item \code{\link{PrintMethodList}} : list of implemented
      simulation methods
      \item \code{\link{ShowModels}} : interactive, graphical choice of
      models
      \item \code{\link{soil}} : Soil physical and chemical data;
      the \code{example} gives a simple geostatistical analysis using
      features of the package      
    }
    \item stationary (and isotropic) max-stable random fields\cr
    \itemize{
     \item \code{\link{CovarianceFct}} : covariance models for
     extremal Gaussian random fields
     \item \code{\link{MaxStableRF}} : simulation of max-stable
     random fields
    }
    
    \item stationary and isotropic Poisson random fields
    (not implemented yet)\cr
%    \itemize{
%      \item \code{%\link{
%	PoissonRF%}
%      } : simulation of Poisson random fields
%    }
  }

  Functions used in diverse simulation methods:
  \itemize{
    \item \code{\link{DeleteRegister}} : deleting internal registers
    \item \code{\link{RFparameters}} : control parameters (advanced settings)
  }

  Functions of general purpose:
   \itemize{
    \item \code{\link{eval.parameters}} : provides an interactive menu 
    \item \code{\link{fractal.dim}} : estimation of the fractal dimension
    \item \code{\link{hurst}} : estimation of the Hurst parameter
    \item \code{\link{regression}} : interactive regression plot
    \item \code{\link{Locator}}, \code{\link{Readline}} :
    \code{\link[graphics]{locator}} and \code{\link[base]{readline}},
    respectively, with storage and replay functionality
  }
 
}
\note{
  Version 1.1.21 : Default value of
  \code{RFparameters()$Storing} is now \code{FALSE}.
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html};
  
  Yindeng Jiang \email{jiangyindeng@gmail.com} (circulant embedding
  methods \sQuote{cutoff} and \sQuote{intrinsic})
}
\references{
  This package is announced in:

  Schlather, M. (2001) Simulation of stationary and isotropic random
  fields. \emph{R-News} \bold{1} (2), 18-20.
  
  Schlather, M. (2004)
  Simulation of random fields and applications. \emph{In preparation.}
}
\section{Acknowledgement}{
  Many thanks to Martin Maechler, Paulo Ribeiro, and Tilmann Gneiting
  for proof-reading parts of the code and the help text for V1.0 of
  this package. 

  The coding of V1.0 has been supported by the EU TMR network
  ERB-FMRX-CT96-0095 on
  ``Computational and statistical methods for the analysis of spatial
  data'' in 1999,
  and by the German Federal Ministry of Research and Technology 
  (BMFT) grant PT BEO 51-0339476C during 2000-03.
}
\keyword{spatial}


