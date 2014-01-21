\name{RFmodel}
\alias{RFmodel}
\alias{RF}
\title{Evaluation operators}
\description{
  Here, those \command{RF_name_} commands are listed that are
  implemented as operators on the
  class of \command{\link[=RMmodel-class]{RMmodel}}.

  \bold{This is an advanced topic which refers to the internal
    representation of models}
}
\section{Implemented model}{
  \tabular{ll}{
    \command{\link{RFcov}} \tab assigns to a covariance model the
    covariance values
    at given locations \cr
    \command{\link{RFcovmatrix}} \tab assigns to a covariance model the
    cross covariance values at given locations  \cr
 %   \command{\link{RFget}} \tab not implemented yet\cr
    \command{\link{RFfctn}} \tab assigns to a covariance model the
    value of the function at given locations\cr
    \command{\link{RFdistr}} \tab generic function
    assigning to distribution family various characteristics of the
    distribution
    \cr
    \command{\link{RFpseudovariogram}} \tab assigns to a  model the
    values of the pseudo variogram
    at given locations\cr
    \command{\link{RFsimulate}} \tab simulation of random fields
    \cr
    \command{\link{RFvariogram}} \tab assigns to a  model the
    values of the (cross-)variogram
    at given locations\cr
  }
}

\note{
  Note not all the functions starting with \code{RF}
  possess a corresponding implementation as operator, e.g. not
  \command{\link{RFfit}} and \command{\link{RFinterpolate}}.
}

\seealso{
   \link{RM},
 \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
 \link{RP},
 \link{RR}
}
\examples{
RFoptions(seed=0)
z <- RFsimulate(model=RMexp(), 1:10)
str(RFgetModel(RFsimulate, show.call = TRUE))
str(RFgetModel(RFsimulate, show.call = FALSE))
\dontshow{RFoptions(seed=NA)}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}


\keyword{spatial}

