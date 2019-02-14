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

\title{Class CLASS_FIT }
\description{ Extension of class \code{\link[=RMmodel-class]{RMmodel}}
    which additionally contains
 the likelihood of the data w.r.t. the covariance model represented by
 the CLASS_CLIST part, the estimated trend of the data if it is a
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
 Objects are only meant to be created by the function
 \code{\link{RFfit}}.
}

\section{Slots}{
  \describe{
    \item{\code{model,formel}:}{See \code{\link[=RMmodel]{RMmodel}}.}
   \item{\code{variab}:}{
      vector of estimated variables. Variables
      are used in the internal representation and
      can be a subset of the parameters.
    }
  \item{\code{param}:}{
      vector of estimated parameters
    }
    \item{\code{covariate}:}{
      to do
    }
   \item{\code{globalvariance}:}{
      to do
    }
   \item{\code{hessian}:}{
      to do
    }
   \item{\code{likelihood}:}{numeric; the likelihood of the data
      w.r.t. the covariance model} 
    \item{\code{AIC}:}{the AIC value for the ml estimation}
    \item{\code{AICc}:}{the corrected AIC value for the ml estimation}
    \item{\code{BIC}:}{the BIC value for the ml estimation}
     \item{\code{residuals}:}{array or of class \code{\link[=RFsp-class]{RFsp}};
      residuals of the data w.r.t. the trend model} 

%    \item{\code{call}:}{See \code{\link[=RMmodel-class]{RMmodel}}.} 
 %    \item{\code{name}:}{See \code{\link[=RMmodel]{RMmodel}}.}
 %   \item{\code{par.model}:}{See \code{\link[=RMmodel]{RMmodel}}.}
 %   \item{\code{par.general}:}{See \code{\link[=RMmodel]{RMmodel}}.}
 %    \item{\code{trend}:}{numeric; the estimated mean of the data (if a
 %     constant mean was specified in the model)} 
  }
}

 
\section{Extends}{
 Class \code{CLASS_CLIST}, directly.
}

\section{Methods}{
 \describe{
 \item{[}{\code{signature(x = CLASS_FIT)}: enables accessing
   the slots via the \code{"["}-operator, e.g. \code{x["likelihood"]}}
 \item{[<-}{\code{signature(x = CLASS_FIT)}: enables replacing
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


\author{Alexander Malinowski; \martin}

\seealso{
  \code{\link[=RMmodel-class]{RMmodel}},
  \command{\link{RFfit}}.
}
 
\keyword{classes}
\keyword{print}
\keyword{hplot}

\examples{\dontshow{StartExample()}
# see RFfit
\dontshow{FinalizeExample()}}
