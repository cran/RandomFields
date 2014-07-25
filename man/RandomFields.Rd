\name{RandomFields-package}
\alias{RandomFields-package}
\alias{RandomFields}
\docType{package}
\title{Simulation and Analysis of Random Fields}
\description{ 
  
  The package \code{RandomFields} offers various tools for
  \enumerate{
    \item{\bold{simulation}} of different kinds
    of random fields, including
    \itemize{
      \item multivariate, spatial, spatio-temporal Gaussian random fields,
      \item Poisson fields, binary fields, Chi2 fields and
      \item max-stable fields.
    }
    It can also deal with non-stationarity and anisotropy of these
    processes and conditional simulation (for Gaussian random fields,
    currently). 
    
  
    \item{\bold{model estimation (ML) and inference (tests)
	for regionalizd variables}} and data analysis,

    \item{\bold{model estimation for (geostatistical) linear (mixed) models}}
  }
}



\details{
  The following features are provided by the package:
  
  \enumerate{
    \item \bold{Simulation}
    \itemize{
      \item \command{\link{RFsimulate}}: Simulation
      of random fields,
      including conditional simulation. For a list of all covariance
      functions and variogram models see \command{\link{RMmodel}}.
      Use \command{\link{plot}} for visualisation of the result.
      
      \item \bold{ Estimation} of parameters
      \itemize{ 
	\item \command{\link{RFfit}} : general function for estimating
	parameters; 
	\item \command{\link{RFhurst}} : estimation of the Hurst parameter
	\item \command{\link{RFfractaldim}} : estimation of the fractal
	dimension  
	\item \command{\link{RFempiricalvariogram}} : calculates
	the empirical variogram 
      }
    }
    
    
    \item \bold{Prediction (for Gaussian random fields)} 
    \itemize{
      \item \command{\link{RFinterpolate}} : kriging, including imputing
    }
    
    
    \item \bold{Inference (for Gaussian random fields)} 
    \itemize{
      \item \command{\link{RFcrossvalidate}} : cross validation
      \item \command{\link{RFratiotest}} : likelihood ratio test
      \item \command{\link[=AIC.RF_fit]{AIC}}, 
      \command{\link[=AICc.RF_fit]{AICc}},
      \command{\link[=BIC.RF_fit]{BIC}}, \command{\link[=anova.RF_fit]{anova}},
      \command{\link[=logLik.RFfit]{logLik}}
    }
   
    
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
      For a quick impression use \code{\link{plot}()model}.
    } 
    
    
    \item \bold{Data and example studies}:
    Some data sets and published code are provided to illustrate the
    syntax and structure of the package functions. 
    \itemize{
      \item \code{\link{soil}} : soil physical data
      \item \code{\link{weather}} : UWME weather data
      \item \code{\link{papers}} : code used in the papers published by
      the author(s)
    }


    \item \bold{Graphics}
    \itemize{
      \item Fitting a covariance function manually
      \command{\link{RFgui}}
      \item the generic function \command{\link[graphics]{plot}}
    }
      
  
      
    \item \bold{S3 and S4 objects}
    \itemize{      
      \item The functions return S4 objects
      if \code{\link[=RFoptions]{spConform=TRUE}}.
      This is the default.
      
      
      If \code{\link[=RFoptions]{spConform=FALSE}},
      simple objects as in version 2 are returned.
      These simple objects are frequently provided with an S3 class
      
      \item \command{\link[graphics]{plot}},
      \command{\link[base]{print}}, \command{\link[base]{summary}},
      sometimes also \command{\link[utils]{str}} recognise these S3 and S4
      objects.
      
      \item
      Further generic functions are available for fitted models,
      see \sQuote{Inference} above. 
      
      \item \bold{Note} that, in many cases, \command{print} will return
      an invisible list.  This list contains the main information of the
      S4 object in an accessible way and is in many cases the
      information obtained from \code{summary}. See examples below.
    }
      
    
    \item \bold{Advanced} users, package programmers
    \itemize{
      \item may use \sQuote{./configure
	--with-tcl-config=/usr/lib/tcl8.5/tclConfig.sh
	--with-tk-config=/usr/lib/tk8.5/tkConfig.sh} to configure R
      \item might  decide on a large variety of arguments of the 
      simulation and estimation procedures using the function
      \command{\link{RFoptions}}
    }
  }

  A list of changings from  Version 2 to Version 3 can be found
  in \link{MajorRevisions}.
}

% In the beta version, the following functionalities are currently
% not available:
% \itemize{
% \item \command{\link{ShowModels}}
% \item numerical evaluation of the covariance function in tbm2
% \item Harvard Rue's Markov fields 
% }




\seealso{
  See also  \link{RC}, \link{RF}, \link{RM}, \link{RP}, \link{RR}
}
 
 
\note{
  The following packages enable further
  features in RandomFields:
  \pkg{optimx}, \pkg{soma}, \pkg{GenSA}, \pkg{minqa}, \pkg{pso},
  \pkg{DEoptim}, \pkg{nloptr}, \pkg{RColorBrewer}, \pkg{colorspace}
}

\section{Update}{
  Current updates are available through \url{http://ms.math.uni-mannheim.de/de/publications/software}.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
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
    Felix Ballani (TU Bergakademie Freiberg; Poisson Polygons, 2014) \cr
    Daphne Boecker (Univ. Goettingen; RFgui, 2011)\cr   
    Katharina Burmeister (Univ. Goettingen; testing, 2012)\cr
    Sebastian Engelke (Univ. Goettingen; RFempiricalvariogram, 2011-12)\cr
    Sebastian Gross (Univ. Goettingen; tilde formulae, 2011)\cr
    Alexander Malinowski (Univ. Mannheim; S3, S4 classes 2011-13)\cr
    Juliane Manitz (Univ. Goettingen; testing, 2012)\cr
    Johannes Martini (Univ. Goettingen; RFempiricalvariogram,
    2011-12)\cr
    Ulrike Ober (Univ. Goettingen; help pages, testing, 2011-12)\cr
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
  Paulo Ribeiro : many comments on Version 3.0.27 and code delivered by
  means of \pkg{geoR}\cr
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
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

# Compare:
model <- RMexp(scale=2) + RMnugget(var=3)
str(model)  ## S4 object

model <- summary(RMexp(scale=2) + RMnugget(var=3))
str(model)  ## list style as in version 2 of RandomFields

\dontshow{FinalizeExample()}
}

\keyword{spatial}


