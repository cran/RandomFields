\name{RFfunction}
\alias{RF}
\alias{RFmodel}
\alias{RFmodels}
\alias{RFfunction}
\alias{RFfunctions}
\title{Evaluation operators (RF commands)}
\description{
  Here, all the \command{RF_name_} commands are listed.
}
\section{Models that are treated internally as operators
  on \code{\link{RMmodels}}}{
    The user's \code{\link{RMmodel}} is supplemented internally by
    operators that are tacitely assumed, e.g. \command{\link{RPgauss}}.

    Further completions of the user's model determine what should be
    done with the model, e.g. calculation of the covariance
    (\command{\link{RFcov}}).
    The following list gives those \code{RFfunctions} that have
    an internal represenation as completion to the user's model.
    \tabular{ll}{
      \bold{User's function} \tab \bold{Internal task}
      \cr
    \command{\link{RFcov}} \tab assigns to a covariance model the
    covariance values
    at given locations
    \cr
    \command{\link{RFcovmatrix}} \tab assigns to a covariance model the
    matrix of covariance values at given locations
    \cr
    %   \command{\link{RFget}} \tab not implemented yet\cr
    \command{\link{RFfctn}} \tab assigns to a model the
    value of the function at given locations. In case of a covariance
    model \command{\link{RFfctn}} is identical to \command{\link{RFcov}}.
    \cr
    \command{\link{RFdistr}} \tab generic function
    assigning to a distribution family various values of the
    distribution % (density function, distribution function, quantile,
    random sample)
    \cr
    \command{\link{RFpseudovariogram}} \tab assigns to a  model the
    values of the pseudo variogram
    at given locations\cr
    \command{\link{RFsimulate}} \tab assigns to a model
    a realisation of the corresponding random field
    \cr
    \command{\link{RFvariogram}} \tab assigns to a  model the
    values of the (cross-)variogram
    at given locations\cr
  }
}

\section{Estimation and Inference}{
 \tabular{ll}{
      \bold{User's function} \tab \bold{Description}
      \cr
     \command{\link{RFcrossvalidate}} \tab cross validation for Gaussian fields
    \cr
    \command{\link{RFempiricalvariogram}} \tab
    empirical variogram
    \cr
    \command{\link{RFfit}} \tab
    (maximum likelihood) fitting of the parameters
    \cr
    \command{\link{RFinterpolate}} \tab \code{'kriging' and 'imputing'}
    \cr
%    \command{\link{RFlikelihood}} \tab not coded yet % to do \cr
    \command{\link{RFratiotest}} \tab likelihood ratio test for Gaussian fields
    \cr
  }
}


\section{Graphics for Gaussian fields}{
  \tabular{ll}{
    \command{\link{RFgui}} \tab educational tool for\cr
    \tab * manual selection of a covariance model\cr
    \tab * manual fitting to the empirical variogram
    \cr
    \command{\link{RFfractaldim}} \tab
    determination of the fractal dimension
    \cr
    \command{\link{RFhurst}} \tab
    determination of the Hurst effect (long range dependence)
    \cr
  }
}


\section{Coordinate transformations}{
 \tabular{ll}{
   \command{\link{RFearth2cartesian}} \tab
   transformation of earth coordinates to cartesian coordinates
    \cr
    \command{\link{RFearth2dist}} \tab transformation of earth
    coordinates to Euclidean distances
  }
}

\section{Information from and to RandomFields}{
  \tabular{ll}{
    \command{\link{RFgetMethodNames}} \tab
    currently implemented list of simulation methods
    \cr
    \command{\link{RFgetModel}} \tab
    returns the model used in a \code{\link{RFfunction}}, with some
    more details
    %on the implementation
   % \cr\tab  including the default values of the arguments of the model
    \cr	  
    \command{\link{RFgetModelInfo}} \tab
    similar to \command{\link{RFgetModel}}, but with detailed
    information on the implementation
    \cr
    \command{\link{RFgetModelNames}} \tab
    lists the implemented models
    \cr
    \command{\link{RFoptions}} \tab
    options of package RandomFields
  }
}

\seealso{
   \link{RC}, \link{RM}, \link{RP},
 \link{RR},
 \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

z <- RFsimulate(model=RMexp(), 1:10)
RFgetModel(RFsimulate, show.call = TRUE)  # user's definition
RFgetModel(RFsimulate, show.call = FALSE) # main internal part

\dontshow{FinalizeExample()}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}


\keyword{spatial}

