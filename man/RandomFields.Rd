\name{RandomFields-package}
\alias{RandomFields-package}
\alias{RandomFields}
\docType{package}
\title{Simulation and Analysis of Random Fields}
\description{ 
  
  The package \code{RandomFields} offers various tools for
  \enumerate{
    \item{\bold{model estimation (ML) and inference (tests)}}
    for regionalized variables and data analysis,

    \item{\bold{simulation}} of different kinds
    of random fields, including
    \itemize{
      \item multivariate, spatial, spatio-temporal, and non-stationary
      Gaussian random fields,
      \item Poisson fields, binary fields, Chi2 fields, t fields and
      \item max-stable fields.
    }
    It can also deal with non-stationarity and anisotropy of these
    processes and conditional simulation (for Gaussian random fields,
    currently). 
    
 
   % \item{\bold{model estimation}} for (geostatistical) linear (mixed) models
  }

  See \url{http://ms.math.uni-mannheim.de/de/publications/software/}
  for \bold{intermediate updates}.
  
}



\details{
  The following features are provided by the package:
  
  \enumerate{
    \item \bold{Bayesian Modelling}
    \itemize{
      \item See \link{Bayesian Modelling} for an introduction to
      hierarchical modelling 
    }

        
     \item \bold{Coordinate systems}
    \itemize{
      \item Cartesian, earth and spherical coordinates are
      recognized, see \link{coordinate systems} for details.
      \item a list of valid models is given by
      \link{spherical models}.
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
     
    \item \bold{ Estimation of parameters (for second-order random fields)}
    \itemize{ 
      \item \command{\link{RFfit}} : general function for estimating
      parameters; (for Gaussian random fields)
      \item \command{\link{RFhurst}} : estimation of the Hurst parameter
      \item \command{\link{RFfractaldim}} : estimation of the fractal
      dimension  
      \item \command{\link{RFempiricalvariogram}} : calculates
      the empirical variogram 
    }
    
   \item \bold{Graphics}
    \itemize{
      \item Fitting a covariance function manually
      \command{\link{RFgui}}
      \item the generic function \command{\link[graphics]{plot}}
      \item global graphical parameters with \command{\link{RFpar}}
    }

    \item \bold{Inference (for Gaussian random fields)} 
    \itemize{
      \item \command{\link{RFcrossvalidate}} : cross validation
      \item \command{\link{RFlikelihood}} : likelihood
     \item \command{\link{RFratiotest}} : likelihood ratio test
      \item \command{\link[=AIC.RF_fit]{AIC}}, 
      \command{\link[=AICc.RF_fit]{AICc}},
      \command{\link[=BIC.RF_fit]{BIC}}, \command{\link[=anova.RF_fit]{anova}},
      \command{\link[=logLik.RFfit]{logLik}}
    }
   
    
    \item \bold{Models}
    \itemize{
      \item For an introduction an general properties, see
      \link{RMmodels}.
      \item For an overview over classes of
      covariance and variogram models --e.g. for
      \bold{geostatistical} purposes-- see \link{RM}. More
      sophisticated models 
      and covariance function operators are included. 
      \item To apply the offered package procedures to \bold{mixed models}
      -- e.g. appearing in genetical data analysis-- see
      \command{\link{RFformula}}.
      \item definite models are evaluated by \command{\link{RFcov}},
      \command{\link{RFvariogram}} and \command{\link{RFcovmatrix}}.
      For a quick impression use \code{\link{plot}(model)}.
      \item non-definte models are evaluated by \command{\link{RFfctn}} and
      \command{\link{RFcalc}}
      \item \command{\link{RFlinearpart}} returns the linear part of a
      model
      \item \command{\link{RFboxcox}} deals explicitely with Box-Cox
      transformations. In many cases it is performed implicitely.
    }
    
    \item \bold{Prediction (for second-order random fields)} 
    \itemize{
      \item \command{\link{RFinterpolate}} : kriging, including imputing
    }
    
    \item \bold{Simulation}
    \itemize{
      \item \command{\link{RFsimulate}}: Simulation
      of random fields,
      including conditional simulation. For a list of all covariance
      functions and variogram models see \command{\link{RM}}.
      Use \command{\link{plot}} for visualisation of the result.
    }
    
    \item \bold{S3 and S4 objects}
    \itemize{      
      \item The functions return S4 objects
      based on the package \pkg{sp},
      if \code{\link[=RFoptions]{spConform=TRUE}}.
      This is the default.
      
      
      If \code{\link[=RFoptions]{spConform=FALSE}},
      simple objects as in version 2 are returned. 
      These simple objects are frequently provided with an S3 class.
      This options makes the returning procedure much faster, but
      currently does not allow for the comfortable use of
      \command{\link[=plot-method]{plot}}.
      
      \item \command{\link[graphics]{plot}},
      \command{\link[base]{print}}, \command{\link[base]{summary}},
      sometimes also \command{\link[utils]{str}} recognise these S3 and S4
      objects

      \item use \command{\link{sp2RF}} for an explicite transformation
      of \pkg{sp} objects to S4 objects of \pkg{RandomFields}.
      
      \item
      Further generic functions are available for fitted models,
      see \sQuote{Inference} above. 
      
%   \item \bold{Note} that, in many cases, \command{print} will return
%      an invisible list.  This list contains the main information of the
%      S4 object in an accessible way and is in many cases the
%      information obtained from \code{summary}. See examples below.
    }
      
    \item  \bold{Xtended} features, especially for package programmers
    \itemize{
      \item might  decide on a large variety of arguments of the 
      simulation and estimation procedures using the function
      \command{\link{RFoptions}}
      \item may use \sQuote{./configure
	--with-tcl-config=/usr/lib/tcl8.5/tclConfig.sh
	--with-tk-config=/usr/lib/tk8.5/tkConfig.sh} to configure R
    }
  }
}

\section{Changings}{
  A list of major changings from  Version 2 to Version 3 can be found
  in \link{MajorRevisions}.

  \link{Changings} lists some further changings, in particular of
  argument and argument names.
}

% In the beta version, the following functionalities are currently
% not available:
% \itemize{
% \item \command{\link{ShowModels}}
% \item numerical evaluation of the covariance function in tbm2
% \item Harvard Rue's Markov fields 
% }


\seealso{
  See also  \link{RF}, \link{RM}, \link{RP}, \link{RR}, \link{RC}, \link{R.}
}
 
 
\note{
  The following packages enable further choices for the optimizer
  (instead of \command{optim}) in RandomFields:
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
    Strokorb, K. (2015) 
    Analysis, simulation and prediction of multivariate
    random fields with package \pkg{RandomFields}. \emph{
      Journal of Statistical Software}, \bold{63} (8), 1-25,
    url          = \sQuote{http://www.jstatsoft.org/v63/i08/}
    \item
    see also the corresponding \href{../doc/multivariate_jss.pdf}{vignette}.
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
    Paulo Ribeiro (Unversidade Federal do Parana; code adopted
    from \pkg{geoR}, 2014)\cr
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


\keyword{spatial}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
# simulate some data first (Gaussian random field with exponential
# covariance; 6 realisations)
model <- RMexp()
x <- seq(0, 10, 0.1)
z <- RFsimulate(model, x, x, n=6)

## select some data from the simulated data
xy <- coordinates(z)
pts <- sample(nrow(xy), min(100, nrow(xy) / 2))
data <- matrix(nrow=nrow(xy), as.vector(z))[pts, ]
data <- cbind(xy[pts, ], data)
plot(z, data)

## re-estimate the parameter (true values are 1)
estmodel <- RMexp(var=NA, scale=NA)
(fit <- RFfit(estmodel, data=data))

## show a kriged field based on the estimated parameters
kriged <- RFinterpolate(fit, x, x, data=data)
plot(kriged, data)

\dontshow{FinalizeExample()}
}

