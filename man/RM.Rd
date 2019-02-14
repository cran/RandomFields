\name{RMmodels Overview}
\alias{RM}
\title{Overview over classes of \code{RMmodels}}
\description{
  Various classes of models \code{RMxxx} are implemented in
  RandomFields, that have their own man pages. Here, an overview over
  these man pages are given.
}
\section{Man pages}{
  Beginners should start with \link{RMmodels}, then go
  for \link{RMmodelsAdvanced} if more information is needed.
 \tabular{ll}{
   \link{RMmodels} \tab general introduction and a collection of simple models
   \cr
   \link{RMmodelsAdvanced} \tab includes more advanced stationary and
   isotropic models, variogram models, non-stationary
   models and trend models
   \cr
   \link{Bayesian} \tab hierarchical models
   \cr
    \link{RMmodelsMultivariate} \tab multivariate covariance models and 
   multivariate trend models
   \cr
    \link{RMmodelsNonstationary} \tab non-stationary covariance models
   \cr
   \link{RMmodelsSpaceTime} \tab space-time covariance models
   \cr
   \link{Spherical models} \tab models based on the polar coordinate
   system, usually used in earth models
   \cr
   \link{Tail correlation functions} \tab  models related to max-stable
   random fields 
   \cr
    \link{trend modelling} \tab how to pass trend specifications
   \cr
 \link[=RFcalc]{Mathematical functions} \tab simple mathematical
   functions that are typically used to build non-stationary covariance models
   and arbitrary trends
   \cr
   \link{RMmodelsAuxiliary} \tab
   rather specialized models, most of them not having positive
   definiteness property, but used internally in certain simulation
   algorithms, for instance.
   \cr
%   \link{RMmodels} \tab
%   \cr
  }
}
\seealso{
  \link{RC}, \link{RR}, \link{RF}, \link{R.}
}

\me

\keyword{spatial}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFgetModelNames(type="positive definite", domain="single variable",
                isotropy="isotropic", operator=!FALSE) ## RMmodel.Rd

\dontshow{FinalizeExample()}}

