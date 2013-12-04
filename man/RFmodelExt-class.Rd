\name{RMmodelExt-class}
\docType{class}
\alias{RMmodelExt-class}
\alias{[,RMmodelExt-method}
\alias{[,RMmodelExt,ANY,ANY-method}
\alias{[,RMmodelExt,ANY,ANY,ANY-method}
\alias{[<-,RMmodelExt-method}
\alias{[<-,RMmodelExt,ANY,ANY-method}
\alias{[<-,RMmodelExt,ANY,ANY,ANY-method}

\title{Class "RMmodelExt" }
\description{ Extension of Class \code{\link[=RMmodel-class]{RMmodel}}
    which additionally contains
 the likelihood of the data w.r.t. the covariance model represented by
 the "RMmodel" part, the estimated trend of the data if it is a
 constant trend, and the residuals of the data w.r.t. the model. Objects
 of this class only occur as slots in the output of "RFfit".
}

\section{Creating Objects}{
 Objects are only ment to be created by the function
 \code{\link{RFfit}}
}

\section{Slots}{
 \describe{
 \item{\code{likelihood}:}{numeric; the likelihood of the data
 w.r.t. the covariance model} 
 \item{\code{trend}:}{numeric; the estimated mean of the data (if a
 constant mean was specified in the model)} 
 \item{\code{residuals}:}{array or of class \code{\link[=RFsp-class]{RFsp}};
 residuals of the data w.r.t. the trend model} 
 \item{\code{call}:}{see \code{\link[=RMmodel-class]{RMmodel}}.} 
 \item{\code{name}:}{see \code{\link[=RMmodel]{RMmodel}}.}
 \item{\code{submodels}:}{see \code{\link[=RMmodel]{RMmodel}}.}
 \item{\code{par.model}:}{see \code{\link[=RMmodel]{RMmodel}}.}
 \item{\code{par.general}:}{see \code{\link[=RMmodel]{RMmodel}}.}
 }
}

 
\section{Extends}{
 Class \code{"RMmodel"}, directly.
}

\section{Methods}{
 \describe{
 \item{[}{\code{signature(x = "RMmodelExt")}: enables accessing
 the slots via the \code{"["}-operator, e.g. \code{x["likelihood"]}}
 \item{[<-}{\code{signature(x = "RMmodelExt")}: enables replacing
 the slots via the \code{"["}-operator}
 
 }
}

%\section{Details}{
%}


\author{Alexander Malinowski \email{malinows@math.uni-goettingen.de}}

\seealso{
  \code{\link[=RMmodel-class]{RMmodel}}
  \command{\link{RFfit}}
}
 
\examples{
set.seed(0)
# see RFfit
}

\keyword{classes}
