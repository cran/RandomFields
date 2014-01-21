\name{RandomFields-package}
\alias{RandomFields-package}
\docType{package}
\title{Simulation and Analysis of Random Fields}
\description{
  \bold{This version 3 is a major revision of version 2 with many
    important changings, see section \sQuote{Major Revision} below.}


  The package \code{RandomFields} offers various tools
 for
 \enumerate{
 \item{\bold{simulation}} of different kinds
 of random fields, including
 \itemize{
 \item multivariate, spatial, spatio-temporal Gaussian random fields,
 \item Poisson fields, binary fields, Chi2 fields and
 \item max-stable fields.
 }
 It can also deal with non-stationarity and anisotropy of these
 processes and conditional simulation (for Gaussian random fields, currently).
 

 \item{\bold{model estimation for regionalizd variables}} and data analysis,
 using different LSQ and MLE procedures (currently for Gaussian random
 fields only)
 \item{\bold{model estimation for (geostatistical) linear (mixed) models}}
 }
 %of parameters or fitting.
 %Hurst parameter, RFfractaldimension, empirical variogram,
 %interactive fitting of parameters, LSQ and MLE estimation of parameters.
 %Basic kriging procedures are also provided.

 % Starting with version 2.0, it also allows for the simulation
 % of random fields that are non-stationary or multivariate or
 % sophisticated space-time fields.
 % fitvario allows for multivariate models and mixed effect models.

% \bold{There are some changings in the definitions and in the output,
% see help("changings")}
}
%\section{Installation}{
% If like to have Havard Rue's GMRF available
% please follow the instruction in \link{RMMethods}.
%}


\section{Major Revision: changings}{
  The following major changings took place with respect to version 2:
  \itemize{
    \item \bold{S4 objects}
    \itemize{
      \item
      \pkg{RandomFields} is now based on S4 objects using the package
      \pkg{sp}. The functions accept both \pkg{sp} objects and simple
      objects as used in version 2.
      
      \item The functions return S4 objects
      if \code{\link[=RFoptions]{spConform=TRUE}}, otherwise simple objects
      as in version 2.

      \item \command{plot} and \command{print} recognise these S4 objects
    }
    
    \item\bold{Documentation}
    \itemize{
      \item each model has now its own man page;
      \item classes of models and functions are bundled in several pages:
      Covariance models start with \code{\link{RM}}, distribution families
      with \code{\link{RR}}, processes with \code{\link{RP}}, user
      functions with \code{\link{RF}}
      \item the man pages of several functions are split into two parts:
      
      (i) a beginners man page which includes a link to
      
      (ii) man pages for advanced users
    }
    
    \item\bold{Interfaces}
    \itemize{
      \item The interfaces become simpler, at the same time more powerful
      then the functions in version 2. E.g., \code{RFsimulate} can perform
      unconditional simulation, conditional simulation and random
      imputing.
      
      \item Only those parameters are kept in the functions that are considered
      as being absolutely necessary. All the other parameters can be
      included as \link[=RFoptions]{options}.
      
      \item \command{\link{RFgui}} is an instructive interface based on tcl/tk,
      replacing the former \code{ShowModels}
    }
    
    
    \item \bold{Inference for Gaussian random fields}    
    \itemize{
      \item \command{\link{RFfit}} has undergone a major revision. E.g.:
    
      (i) estimation random effects model with spatial covariance
      structure

      (ii) automatic estimation of 10 and more parameters in
      multivariate and/or space-time models

      \item \command{\link{RFempiricalvariogram}} is now based on an fft algorithm
      if the data are on a grid, even allowing for missing values.

      \item \command{\link{RFratiotest}} has been added.
    }
    
    \item \bold{Processes}
   \itemize{
     \item{Maxstable processes} modelling of \link[=maxstable]{maxstable
       processes} has been enhanced, including
     
     (i) the simulation of Brown-Resnick processes
     
     (ii) initial support of
     \link[=RMmodelsAdvanced]{tail correlation functions};
     \item{Further processes}
     \link[=RPchi2]{chi2 processes}, \link[=RPpoisson]{compound
       Poisson processes}, \link[=RPbernoulli]{binary processes} added.
   }
   
   \item \bold{Models}
   \itemize{
     \item the \link[=RFformula]{formula notation} for
     linear models may now be defined
     \item{Novel, user friendly definition of the covariance models}
     \item \link[=RMmodelsAdvanced]{Multivariate and vector valued random
       fields} are now fully included 
     \item The \link[=RMuser]{user}
       may now define his own functions, to some extend.
     \item The \link[=RMtrend]{trend} allows for much more flexibility
     \item \link[=RR]{Random parameters} may now included which will be
     extended to \link{Baysian} modelling in future.
   }
 }
  }


\details{
 The following functions are 
 provided by the package:
 
 \enumerate{
 \item \bold{Simulation}
 \itemize{ \item \command{\link{RFsimulate}}: Simulation
 of random fields,
 including conditional simulation. For a list of all covariance
 functions and variogram models see \command{\link{RMmodel}}. 
 
 % \item \command{\link{RMmodel}},
 
 % \link{sophisticated} models: covariance functions and
 % variogram models
		
 
 % \item \command{\link{RFsimulate}} : simulation of Gaussian random
 % fields; nice examples to get familiar with the
 % simulation features of the package;	
 
 
 % \item \command{\link{PrintMethodList}} : list of implemented
 % simulation methods
 % \item \command{\link{ShowModels}} : interactive, graphical choice of
 % models -- currently not available.
 
 % Use ShowModels in version 1.3.x
 % \item \command{\link{soil}} : Soil physical and chemical data;
 % the \code{example} gives a simple geostatistical analysis using
 % features of the package 
 }

 % \item stationary (and isotropic) max-stable random fields\cr
 % \itemize{
 % \item \command{\link{RMmodel}} : covariance models for
 % extremal Gaussian random fields
 % \item \command{\link{MaxStableRF}} : simulation of max-stable
 % random fields

 % \bold{Note: Simulation algorithm for Brown-Resnick processes will
 % come up soon.}
 % }
 
 \item \bold{ Estimation} of paramters
 \itemize{ 
 \item \command{\link{RFfit}} : general function for estimating parameters;
 \item \command{\link{RFratiotest}} : Monte Carlo likelihood ratio test
 for Gaussian random fields
 \item \command{\link{RFhurst}} : estimation of the Hurst parameter
 \item \command{\link{RFfractaldim}} : estimation of the fractal dimension 
 \item \command{\link{RFempiricalvariogram}} : calculates
 the empirical variogram 
 }
 
 \item \bold{Prediction} 
 \itemize{
 \item \command{\link{RFinterpolate}} : kriging, including imputing
 }
 
 % \item stationary and isotropic Poisson random fields
 % (not implemented yet)\cr
 % \itemize{
 % \item \code{%\link{
 %	PoissonRF%}
 % } : simulation of Poisson random fields
 % }
 
 \item \bold{Models}
 \itemize{ 
 \item For a list of covariance and variogram models --e.g. for
 \bold{geostatistical} purposes-- see \command{\link{RMmodel}}. More
 sophisticated models 
 and covariance function operators are included. 
 \item To apply the offered package procedures to \bold{mixed models}
 -- e.g. appearing in genetical data analysis-- see
 \command{\link{RFformula}}.
 \item models are evaluated by \command{\link{RFcov}},
 \command{\link{RFvariogram}} and \command{\link{RFcovmatrix}}.
 } 
 
 

 \item \bold{Data and example studies}:
 Some data sets and published code are provided to illustrate the syntax and structure of the package functions.
 \itemize{
 \item \code{\link{soil}} : soil physical data
 \item \code{\link{weather}} : UWME weather data
 \item \code{\link{papers}} : code used in the papers published by
 the author(s)
 }
 
 %Functions used in diverse simulation methods:
 % \itemize{
 % \item \command{\link{DeleteRegister}} : deleting internal registers
 % \item \command{\link{RFoptions}} : control parameters (advanced settings)
 % }
 
 % 
 % \itemize{
 % \item \command{\link{eval.parameters}} : provides an interactive menu 
 % \item \command{\link{Print}} : nice print function
 % \item \command{\link{regression}} : interactive regression plot
 % \item \command{\link{Locator}}, \command{\link{Readline}} :
 % \command{\link[graphics]{locator}} and \command{\link[base]{readline}},
 % respectively, with storage and replay functionality
 % \item \command{\link{sleep.milli}} : sets the process into sleeping status
 % }

 \item \bold{Advanced} users, package programmers
 \itemize{
 \item may use \sQuote{./configure
 --with-tcl-config=/usr/lib/tcl8.5/tclConfig.sh
 --with-tk-config=/usr/lib/tk8.5/tkConfig.sh} to configure R
\item might  decide on a large variety of parameters of the 
 simulation and estimation procedures using the function
 \command{\link{RFoptions}}
 
 
 }
 %Sophisticated options for the provided simulation and estimation procedures can be found here: \command{\link{RFoptions}}.
 % \command{\link{RFtools}}. 
}}

% In the beta version, the following functionalities are currently
% not available:
% \itemize{
% \item \command{\link{ShowModels}}
% \item numerical evaluation of the covariance function in tbm2
% \item Harvard Rue's Markov fields 
% }
 
\note{
  The following packages enable further
  features in RandomFields:
  \pkg{optimx}, \pkg{soma}, \pkg{GenSA}, \pkg{minqa}, \pkg{pso},
  \pkg{DEoptim}, \pkg{nloptr}, \pkg{RColorBrewer}, \pkg{colorspace}
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de};
}
\references{
  \itemize{
    \item
    Singleton, R.C. (1979). In \emph{Programs for Digital Signal Processing} 
    Ed.: Digital Signal Processing Committee and IEEE Acoustics,
    Speech, and Signal Processing Committe (1979)
    IEEE press.
    \item
    Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and
    Strokorb, K. (2013) 
    \R package \pkg{RandomFields}: Analysis and simulation of
    multivariate random fields and more. \emph{Submitted to JSS}
  } 
 }
\section{Contributions}{
  \itemize{
    \item Contributions to version 3.0 and following:\cr
    Felix Ballani (TU Bergakademie Freiberg; Poisson Polygons, 2013) \cr
    Daphne Boecker (Univ. Goettingen; RFgui, 2011)\cr
    Katharina Burmeister (Univ. Goettingen; testing, 2012)\cr
    Sebastian Engelke (Univ. Goettingen; RFempiricalvariogram, 2011-12)\cr
    Sebastian Gross (Univ. Goettingen; tilde formulae, 2011)\cr
    Alexander Malinowski (Univ. Mannheim; S3, S4 classes 2011-13)\cr
    Juliane Manitz (Univ. Goettingen; testing, 2012)\cr
    Johannes Martini (Univ. Goettingen; RFempiricalvariogram, 2011-12)\cr
    Marco Oesting (Univ. Mannheim; Brown-Resnick processes, Kriging, Trend,
    2011-13)\cr 
    Kirstin Strokorb (Univ. Mannheim; help pages, 2011-13)\cr
    \item Contributions to version 2.0 and following:\cr
    Peter Menck (Univ. Goettingen; multivariate circulant embedding)\cr
    R Core Team, Richard Singleton (fft.c and advice) 
    \item Contributions to version 1 and following:\cr
    Ben Pfaff, 12167 Airport Rd, DeWitt MI 48820, USA making available
    an algorithm for AVL trees (avltr*)
  }
}
\section{Thanks}{
  Patrick Brown : comments on Version 3\cr
  Paulo Riberio : comments on Version 1\cr
  Martin Maechler : advice for Version 1
}
\section{Financial support}{
 \itemize{
 \item
 V3.0 has been financially supported by the German Science Foundation
 (DFG) through the Research Training Group 1953 \sQuote{Statistical
 Modeling of Complex Systems and Processes --- Advanced Nonparametric
 Approaches} (2013-2018).
 \item
 V3.0 has been financially supported by Volkswagen Stiftung within
 the project \sQuote{WEX-MOP} (2011-2014).
 \item
 Alpha versions for V3.0 have been
 financially supported by the German Science Foundation (DFG) through the
 Research Training Groups 1644 \sQuote{Scaling problems in Statistics}
 and 1023 \sQuote{Identification in Mathematical Models} (2008-13).
 \item
 V1.0 has been financially supported by
 the German Federal Ministry of Research and Technology 
 (BMFT) grant PT BEO 51-0339476C during 2000-03.
 \item
 V1.0 has been financially supported by the EU TMR network
 ERB-FMRX-CT96-0095 on
 ``Computational and statistical methods for the analysis of spatial
 data'' in 1999.
 }
}
\examples{
# Best to start with the examples in 'RFsimulate'
\dontshow{RFoptions(seed=NA)}
}

\keyword{spatial}


