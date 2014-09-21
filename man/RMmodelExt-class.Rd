\name{RMmodelFit-class}
\docType{class}
\alias{RMmodelFit-class}
\alias{RM_modelFit-class}
\alias{[,RMmodelFit-method}
\alias{[,RMmodelFit,ANY,ANY-method}
\alias{[,RMmodelFit,ANY,ANY,ANY-method}
\alias{[<-,RMmodelFit-method}
\alias{[<-,RMmodelFit,ANY,ANY-method}
\alias{[<-,RMmodelFit,ANY,ANY,ANY-method}
\alias{show,RMmodelFit-method}
\alias{print,RMmodelFit-method}
\alias{anova,RMmodelFit-method}
\alias{summary,RMmodelFit-method}
\alias{print.RMmodelFit}

\alias{print.RM_modelFit}
\alias{anova.RM_modelFit}
\alias{summary.RM_modelFit}

\title{Class "RMmodelFit" }
\description{ Extension of Class \code{\link[=RMmodel-class]{RMmodel}}
    which additionally contains
 the likelihood of the data w.r.t. the covariance model represented by
 the "RMmodel" part, the estimated trend of the data if it is a
 constant trend, and the residuals of the data w.r.t. the model. Objects
 of this class only occur as slots in the output of "RFfit".
}

%\usage{
%anova.RM_modelFit(object, ...)
%print.RM_modelFit(x, ...)
%summary.RM_modelFit(object, ..., isna.param)
%}
%\arguments{  
%  \item{object, x, ...}{
%    see the respective generic function
%  }
%  \item{isna.param}{
%    logical. Weather the vector of parameters should be considered as
%    having \code{NA}s. In this case, not the parameters, but the
%    internally used variables are reported.    
%  }
%
%}
\section{Creating Objects}{
 Objects are only ment to be created by the function
 \code{\link{RFfit}}
}

\section{Slots}{
  \describe{
    \item{\code{AIC}:}{the AIC value for the ml estimation}
    \item{\code{AICc}:}{the corrected AIC value for the ml estimation}
    \item{\code{BIC}:}{the BIC value for the ml estimation}
    \item{\code{call}:}{see \code{\link[=RMmodel-class]{RMmodel}}.} 
    \item{\code{likelihood}:}{numeric; the likelihood of the data
      w.r.t. the covariance model} 
    \item{\code{name}:}{see \code{\link[=RMmodel]{RMmodel}}.}
    \item{\code{par.model}:}{see \code{\link[=RMmodel]{RMmodel}}.}
    \item{\code{par.general}:}{see \code{\link[=RMmodel]{RMmodel}}.}
    \item{\code{param}:}{
      vector of estimated parameters
    }
    \item{\code{residuals}:}{array or of class \code{\link[=RFsp-class]{RFsp}};
      residuals of the data w.r.t. the trend model} 
    \item{\code{submodels}:}{see \code{\link[=RMmodel]{RMmodel}}.}
    \item{\code{trend}:}{numeric; the estimated mean of the data (if a
      constant mean was specified in the model)} 
    \item{\code{variab}:}{
      vector of estimated variables. Variables
      are used in the internal representation and
      can be a subset of the parameters.
    }
 }
}

 
\section{Extends}{
 Class \code{"RMmodel"}, directly.
}

\section{Methods}{
 \describe{
 \item{[}{\code{signature(x = "RMmodelFit")}: enables accessing
   the slots via the \code{"["}-operator, e.g. \code{x["likelihood"]}}
 \item{[<-}{\code{signature(x = "RMmodelFit")}: enables replacing
   the slots via the \code{"["}-operator}
  \item{show}{\code{signature(x = "RFfit")}: returns the structure
    of \code{x}}
  \item{print}{\code{signature(x = "RFfit")}: identical with
    \command{show}-method}

  \item{anova}{performs a likelihood ratio test base on a chisq approximation
  }
  \item{summary}{gives a summary}

 }
}

%\section{Details}{
%}


\author{Alexander Malinowski \email{malinows@math.uni-goettingen.de},
  Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\seealso{
  \code{\link[=RMmodel-class]{RMmodel}}
  \command{\link{RFfit}}
}
 
\keyword{classes}
\keyword{print}
\keyword{hplot}

\examples{
# see RFfit
\dontshow{FinalizeExample()}
}
