\name{PrintModelList} 
\alias{PrintModelList}
\alias{GetModelList}
\alias{GetModelNames}
\title{Information about the implemented covariance models}
\description{
  \code{PrintModelList} prints the list of currently implemented models
  including the corresponding simulation methods
}
\usage{

PrintModelList(operators=FALSE, internal=FALSE, newstyle=TRUE)

}
\arguments{
  \item{operators}{logical. Flag whether operators should be
    also considered.
  }
  \item{internal}{logical. Flag whether internal models should be
    also considered. In case of \command{PrintModelList} and
    \code{internal=2}, variants of internal models are also printed.
  }
  \item{newstyle}{
    logical. If \code{FALSE} then only the old style model names (Version 2
    and earlier) are shown. These names can still be used in the list
    definition of models, see \link{RMmodelsAdvanced}.
    If \code{TRUE} then the standard names will also be shown.
  }
}
\value{
  \code{PrintModelList} prints a table of the currently implemented covariance
  functions and the matching methods.
  \code{PrintModelList} returns \code{NULL}.
}
\details{
 See \link{RMmodel} for a description of the models and their use.
}
\note{
  From version 3.0 on, the command \code{PrintModelList()}
  is replaced by the call
  \code{\link{RFgetModelNames}(internal=FALSE)}
}



 \examples{
set.seed(0)
%  library(RandomFields)
PrintModelList()
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}
