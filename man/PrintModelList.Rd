\name{PrintModelList}
\alias{PrintModelList}
\alias{GetModelList}
\alias{GetModelNames}
\title{Information about the implemented covariance models}
\description{
  \code{PrintModelList} prints the list of currently implemented models
  including the corresponding simulation methods

  \code{GetModelList} returns a matrix of currently implemented models
  and their simulation methods

  \code{GetModelNames} returns a list of currently implemented models

}
\usage{

PrintModelList(operators=FALSE, internal=FALSE)

GetModelList(abbr=TRUE, internal=FALSE)

GetModelNames(stationary = c("any", "stationary", "variogram", 
                             "irf", "auxmatrix", "auxvector", "nonstationary",
                             "genvariogram", "gencovariance", "prev.model"),
              isotropy = c("any", "isotropic", "spaceisotropic",
                           "anisotropic", "prev.model"),
              multivariate = -9999:9999,
              operator, normalmix, finiterange, methods, internal=FALSE,
              dim = 1)
}
\arguments{
  \item{operators}{logical. Flag whether operators should be
    also considered.
  }
  \item{internal}{logical. Flag whether internal models should be
    also considered. In case of \command{PrintModelList} and
    \code{internal=2}, variants of internal models are also printed.
  }
  \item{abbr}{logical or numerical. If \code{TRUE} the names for the methods are
    abbreviated. If numerical, \code{abbr} gives the number of letters.}
  \item{stationary}{Selection criterion for the models}
  \item{isotropy}{Selection criterion for the models}
  \item{multivariate}{the dimension of the multivariate model}
  \item{operator}{logical. If missing any covariance model is returned;
    otherwise, the (non)primitive models are returned if \code{FALSE}
    (\code{TRUE});
    see \command{\link{Covariance}} for details.
  }
  \item{normalmix}{logical. If missing any covariance model is returned;
    otherwise all normal scale mixture models are returned if
    \code{TRUE}}
  \item{finiterange}{logical. If missing any covariance model is returned;
    otherwise all models with finite range are returned if
    \code{TRUE}; note that some option are excluding, e.g. no normal
    scale model can have a finite range.}
  \item{methods}{If missing any covariance model is returned;
    otherwise the models with respective attribute are returned;
    see \code{\link{PrintMethodList}} for a complete list of methods}
  \item{dim}{dimension in which the model should be valid;
    there are some models where the validity in a given dimension
    depends on the some parameter; those models are returned partially
    only.}
}
\value{
  \code{PrintModelList} prints a table of the currently implemented covariance
  functions and the matching methods.
  \code{PrintModelList} returns \code{NULL}.

  \code{GetModelNames} returns a list of implemented models
}
\details{
  \bold{See \command{\link{CovarianceFct}} and
   \link{sophisticated} models for a description of the models
    and their use}
}
\seealso{
 \link{sophisticated}, \command{\link{CovarianceFct}},
 \command{\link{GetModel}} 
 }
 \examples{
%  library(RandomFields)
PrintModelList()
GetModelList(abbr=TRUE)
GetModelNames("stationary", "spaceisotropic")
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}
}
\keyword{spatial}
