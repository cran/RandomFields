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

  Starting with version 2.0, it also allows for the simulation
  of random fields that are non-stationary or multivariate or
  sophisticated space-time fields.
  fitvario allows for multivariate models and mixed effect models.

  \bold{There are some changings in the definitions and in the output,
    see help("changings")}
}
%\section{Installation}{
%  If like to have Havard Rue's GMRF available,
%  please follow the instruction in \link{RFMethods}.
%}
\details{
  The following random fields and related functionalities are
  provided by the package.
  \enumerate{
    \item stationary and isotropic Gaussian random fields\cr
    \itemize{
      \item \command{\link{CondSimu}} : conditional simulation
      \item \command{\link{CovarianceFct}},
      \link{sophisticated} models:  covariance functions and
      variogram models
      \item \command{\link{EmpiricalVariogram}} : empirical variogram
      \item \command{\link{GaussRF}} : simulation of Gaussian random
        fields; nice examples to get familiar with the
        simulation features of the package;	
      \item \command{\link{Kriging}} : simple and ordinary kriging
      \item \command{\link{fitvario}} : variogram/covariance function fit
      by least squares, maximum likelihood and cross validation
      techniques
      \item \command{\link{PrintMethodList}} : list of implemented
      simulation methods
      \item \command{\link{ShowModels}} : interactive, graphical choice of
      models -- currently not available.

      Use ShowModels in version 1.3.x
      \item \command{\link{soil}} : Soil physical and chemical data;
      the \code{example} gives a simple geostatistical analysis using
      features of the package      
    }
    \item stationary (and isotropic) max-stable random fields\cr
    \itemize{
     \item \command{\link{CovarianceFct}} : covariance models for
     extremal Gaussian random fields
     \item \command{\link{MaxStableRF}} : simulation of max-stable
     random fields

     \bold{Note: Simulation algorithm for Brown-Resnick processes will
    come up soon.}
    }
    
    \item stationary and isotropic Poisson random fields
    (not implemented yet)\cr
%    \itemize{
%      \item \code{%\link{
%	PoissonRF%}
%      } : simulation of Poisson random fields
%    }
  }

  Data and example code:
  \itemize{
    \item \code{\link{soil}} : soil physical data
    \item \code{\link{papers}} : code used in the papers published by
    the author
    \item \code{\link{weather}} : UWME weather data
  }
  
  Functions used in diverse simulation methods:
  \itemize{
    \item \command{\link{DeleteRegister}} : deleting internal registers
    \item \command{\link{RFparameters}} : control parameters (advanced settings)
  }

  Functions of general purpose:
   \itemize{
%    \item \command{\link{eval.parameters}} : provides an interactive menu 
    \item \command{\link{fractal.dim}} : estimation of the fractal dimension
    \item \command{\link{hurst}} : estimation of the Hurst parameter
    \item \command{\link{Print}} : nice print function
    \item \command{\link{regression}} : interactive regression plot
    \item \command{\link{Locator}}, \command{\link{Readline}} :
    \command{\link[graphics]{locator}} and \command{\link[base]{readline}},
    respectively, with storage and replay functionality
    \item \command{\link{sleep.milli}} : sets the process into sleeping status
  }


  In the beta version, the following functionalities are currently
  not available:
  \itemize{
    \item \command{\link{ShowModels}}
    \item numerical evaluation of the covariance function in tbm2
    \item Harvard Rue's Markov fields 
  }
 
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de};
  }
\references{
  This package is announced in and can be cited by:

  Schlather, M. (2001) Simulation of stationary and isotropic random
  fields. \emph{R-News} \bold{1} (2), 18-20. 
}
\section{Acknowledgement}{
  Many thanks to 
  \itemize{
    %    \item Marco Oesting     
    \item
    Peter Menck implemented the multivariate circulant embedding
    for version 2.0. 
    \item
    Yindeng Jiang \email{jiangyindeng@gmail.com} implemented the
    circulant embedding methods \sQuote{cutoff} and \sQuote{intrinsic}
    in 2004 for the versions 1.2.
     \item
    Martin Maechler, Paulo Ribeiro, and Tilmann Gneiting
    were proof-reading parts of the code and the help text for the versions
    1.0.
  }
}
\section{Financial support}{
  \itemize{
    \item
    V1.0 has been  financially supported by
    the German Federal Ministry of Research and Technology 
    (BMFT) grant PT BEO 51-0339476C during 2000-03.
    \item
    V1.0 has been financially supported by the EU TMR network
    ERB-FMRX-CT96-0095 on
    ``Computational and statistical methods for the analysis of spatial
    data'' in 1999.
  }
}
\keyword{spatial}


