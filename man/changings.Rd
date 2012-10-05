\name{Changings}
\alias{changings}
\title{Changings}
\description{
  Here, only the important inconsistencies
  between Version 1.3 to 2.0 of RandomFields are listed.
  
  Documentation of the extensions will be enhanced later.
}
\details{
  \itemize{
    \item output of \code{fitvario}\cr
    \code{$vario} is defunctioned. Use, for instance, \code{$ml$model} instead of
    \code{$vario$ml}.

    \item The options for \code{mle.methods} in \command{fitvario}
    are currently restricted to \sQuote{ml}. This options 
    \sQuote{reml} and \sQuote{rml1} have not been tested yet.

    \item the \code{list} definition of sophisticated covariance models
    has become obsolete. Furthermore, the parameter \code{kappa}
    for the model parameters is not accepted anymore. Instead \code{k}
    should be used.

    See \code{help("sophisticated")} for the new definition
    of sophisticated covariance models

    \item
    Currently, \command{ShowModels} does not work, neither
    Havard Rue's Markov fields. Both should
    be available in summer 2011 again.

    \item
    Currently, the numerical evaluation of the covariance function in
    tbm2 is not implemented. This should be available in
    summer 2011 again.

  }

  Please let the author know if further inconsistencies exist.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de};
  }
\keyword{spatial}
