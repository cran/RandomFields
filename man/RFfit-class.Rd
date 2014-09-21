\name{RFfit-class}
\docType{class}
\alias{RFfit-class}
\alias{RF_fit-class}
\alias{show,RFfit-method}
\alias{persp,RFfit-method}
\alias{print,RFfit-method}
\alias{anova,RFfit-method}
\alias{AIC,RFfit-method}
\alias{BIC,RFfit-method}
\alias{summary,RFfit-methodt}
\alias{[,RFfit-method}
\alias{[,RFfit,ANY,ANY-method}
\alias{[,RFfit,ANY,ANY,ANY-method}
\alias{coerce,RFfit,RFempVariog-method}
\alias{print.RFfit}

\alias{AICc.RFfit}
\alias{logLik.RFfit}

\alias{print.RF_fit}
\alias{anova.RF_fit}
\alias{AIC.RF_fit}
\alias{BIC.RF_fit}
\alias{AICc.RF_fit}
\alias{summary.RF_fit}
\alias{logLik.RF_fit}
\alias{.RFfit}
\alias{.RF_fit}
%\alias{plot,RFfit-method}

\title{Class \code{RFfit}}
\description{ Class for RandomField's representation of model estimation
 results 
}

%anova.RF_fit(object, ...)
%AIC.RF_fit(object, ..., k=2, method="ml", full=TRUE)
%BIC.RF_fit(object, ..., method="ml", full=TRUE)
%summary.RF_fit(object, ...,  method="ml", full=FALSE)
%print.RF_fit(x, ...,  method="ml", full=FALSE)
%logLik.RF_fit(object, REML = FALSE, ..., method="ml")


\arguments{  
  \item{object,  ...}{
    see the respective generic function;
    The argument
  }
  \item{method}{ string.
    selects that set of parameters which has been obtained through a specific
    fitting method
  }
%  \item{k}{integer, see \command{\link[stats]{AIC}}. Unused here.  }
  \item{full}{logical.
    if \code{TRUE} submodels are reported as well (if available).
  }
}

\section{Creating Objects}{
 Objects are created by the function 
 \command{\link{RFfit}}
}

\section{Slots}{

\describe{
 \item{\code{autostart}:}{RMmodelFit; contains the estimation results
   for the method 'autostart' including a likelihood value, a constant
   trend and the residuals} 
 \item{\code{boxcox}:}{logical; whether the
   parameter of a Box Cox tranformation has been estimated
 }
 \item{\code{coord.units}:}{string giving the units of the coordinates,
   see also option \code{coordunits} of \command{\link{RFoptions}}.
 }
 \item{\code{deleted}:}{integer vector.
   positions of the parameters that has been deleted to get the set of
   variables, used in the optimazation
 }
 \item{\code{ev}:}{list; list of objects of class
   \code{\link[=RFempVariog-class]{RFempVariog}}, 
   contains the empirical variogram estimates of the data} 
 \item{\code{fixed}:}{
   list of two vectors. The fist gives the position where the
   parameters are set to zero. The second gives the position where the
   parameters are set to one.
 }
 \item{\code{internal1}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{internal2}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{internal3}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{lowerbounds}:}{RMmodel; covariance model in which each
   parameter value gives the lower bound for the respective parameter} 
 \item{\code{ml}:}{RMmodelFit; analog to slot 'autostart'
 }
 \item{\code{modelinfo}:}{ Table with information on the parameters:
   name, boundaries, type of   parameter 
 }
 \item{\code{n.covariates}:}{   number of covariates
 }
 \item{\code{n.param}:}{
   number of parameters (given by the user)
 }
 \item{\code{n.variab}:}{
   number of variables (used internally);
   \code{n.variab} is always less than or equal to \code{n.param}
 }
 \item{\code{number.of.data}:}{
   the number of data values passed to \command{\link{RFfit}} that are
   not \code{NA} or \code{NaN}
 }
 \item{\code{number.of.parameters}:}{
   total number of parameters of the model that had to be estimated
   including variances, scales, co-variables, etc.
 }
 \item{\code{p.proj}:}{vector of integers. The original position of those
   parameters that are used in the submodel
 }
 \item{\code{plain}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{report}:}{
   if not empty, it indicates that this model should be reported
   and gives a standard name of the model.
   
   Various function, e.g. \command{print.RMmodelFit} uses
   this information if their argument \code{full} equals \code{TRUE}.
  
 }
 \item{\code{self}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{sd.inv}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{sqrt.nr}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{submodels}:}{
   list. Sequence  (in some cases even nested sequence)
   of models that is used to determine an initial value in
   \command{}
 }
 \item{\code{table}:}{matrix; summary of estimation results of
   different methods} 
 \item{\code{transform}:}{function; } 
 \item{\code{true.tsdim}:}{
   time space dimension of the (original!) data,
   even for submodels that consider parts of separable models.
 }
 \item{\code{true.vdim}:}{
   multivariability of the (original!) data,
   even for submodels that consider independent models
   for the multivariate components.
 }
 \item{\code{upperbounds}:}{RMmodel; see slot 'lowerbounds'} 
 \item{\code{users.guess}:}{RMmodelFit; analog to slot 'autostart'} 
 \item{\code{ml}:}{RMmodelFit; analog to slot 'autostart'; with maximum
   likelihood method}
 \item{\code{v.proj}:}{vector of integers.
   The components selected in one of the submodels
 }
 \item{\code{variab.units}:}{string giving the units of the variables,
   see also option \code{varunits} of \command{\link{RFoptions}}.
 }
 \item{\code{x.proj}:}{
   logical or integer. If logial, it means that no
   separable model is considered there. If integer, then
   it gives the considered directions of a separable model
 }
 \item{\code{Z}:}{
   standardized list of information on the data
 }
 }
}
 
%\section{Extends}{
%}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "RFfit")}: gives a plot of the
      empirical variogram together with the fitted model, for more details see
      \command{\link{plot-method}}.
    }
    \item{show}{\code{signature(x = "RFfit")}: returns the structure
      of \code{x}
    }

    \item{persp}{code{signature(obj =
	"RFfit")}: generates \command{\link[graphics]{persp}} plots
    }
    \item{print}{\code{signature(x = "RFfit")}: identical with
      \command{show}-method, additional argument is \code{max.level}
    }
    \item{[}{\code{signature(x = "RFfit")}: enables accessing
      the slots via the \code{"["}-operator, e.g. \code{x["ml"]}
    }
    \item{as}{\code{signature(x = "RFfit")}:
      converts into other formats, only implemented for target class 
      \code{\link[=RFempVariog-class]{RFempVariog}}
    }
    \item{anova}{performs a likelihood ratio test base on a chisq approximation
    }
    \item{summary}{provides a summary}
    \item{logLik}{provides an object of class \code{"logLik"}
    }
    \item{AIC,BIC}{provides the AIC and BIC information, respectively}
  }
}

\section{Further 'methods'}{
  
  \code{AICc.RFfit(object, ..., method="ml", full=FALSE)}]
  
  \code{AICc.RF_fit(object, ..., method="ml", full=TRUE)}
}

%\section{Details}{
%}

\author{Alexander Malinowski
  \email{alexander.malinowski@math.uni-goettingen.de},
  Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\seealso{
 \code{\link{RFfit}},
 \code{\link{RFempiricalvariogram}},
 \code{\link{RMmodel-class}},
 \code{\link{RMmodelFit-class}}
 \code{\link{plot-method}}
}
 

\references{
  AICc:
  \itemize{
    \item Hurvich, C.M. and Tsai, C.-L. (1989)
    Regression and Time Series Model Selection in Small Samples
    \emph{Biometrika}, \bold{76}, 297-307.  
  }
}

\examples{
# see RFfit
\dontshow{FinalizeExample()}
}

\keyword{classes}
\keyword{print}
\keyword{hplot}
