\name{RFfit-class}
\docType{class}
\alias{RFfit-class}
\alias{show,RFfit-method}
\alias{print.RFfit}
\alias{[,RFfit-method}
\alias{[,RFfit,ANY,ANY-method}
\alias{[,RFfit,ANY,ANY,ANY-method}
\alias{coerce,RFfit,RFempVariog-method}
%\alias{plot,RFfit-method}

\title{Class \code{RFfit}}
\description{ Class for RandomField's representation of model estimation
 results 
} 

\section{Creating Objects}{
 Objects are created by the function 
 \command{\link{RFfit}}
}

\section{Slots}{
 \describe{
 \item{\code{ev}:}{list; list of objects of class
 \code{\link[=RFempVariog-class]{RFempVariog}}, 
 contains the empirical variogram estimates of the data} 
 \item{\code{table}:}{matrix; summary of estimation results of
 different methods} 
 \item{\code{lowerbounds}:}{RMmodel; covariance model in which each
 parameter value gives the lower bound for the respective parameter} 
 \item{\code{upperbounds}:}{RMmodel; see slot 'lowerbounds'} 
 \item{\code{transform}:}{function; } 
 \item{\code{autostart}:}{RMmodelExt; contains the estimation results
 for the method 'autostart' including a likelihood value, a constant
 trend and the residuals} 
 \item{\code{self}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{plain}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{sqrt.nr}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{sd.inv}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{internal1}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{internal2}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{internal3}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{users.guess}:}{RMmodelExt; analog to slot 'autostart'} 
 \item{\code{ml}:}{RMmodelExt; analog to slot 'autostart'; with maximum
   likelihood method}
 \item{\code{coord.units}:}{string giving the units of the coordinates,
   see also option \code{coord_units} of \command{\link{RFoptions}}.
 }
 \item{\code{variab.units}:}{string giving the units of the variables,
   see also option \code{variab_units} of \command{\link{RFoptions}}.
 }
 \item{\code{number.of.data}:}{
   the number of data values passed to \command{\link{RFfit}} that are
   not \code{NA} or \code{NaN}
 }
 \item{\code{number.of.parameters}:}{
   total number of parameters of the model that had to be estimated
   including variances, scales, co-variables, etc.
 }
 }
}
 
%\section{Extends}{
%}

\section{Methods}{
 \describe{
 \item{plot}{\code{signature(x = "RFfit")}: gives a plot of the
 empirical variogram together with the fitted model, for more details see
 \command{\link{plot-method}}.}
 \item{show}{\code{signature(x = "RFfit")}: returns the structure
 of \code{x}}
 \item{print}{\code{signature(x = "RFfit")}: identical with
 \command{show}-method, additional argument is \code{max.level}}
 \item{[}{\code{signature(x = "RFfit")}: enables accessing
 the slots via the \code{"["}-operator, e.g. \code{x["ml"]}}
 \item{as}{\code{signature(x = "RFfit")}:
 converts into other formats, only implemented for target class 
 \code{\link[=RFempVariog-class]{RFempVariog}} } 
 }
}

%\section{Details}{
%}

\author{Alexander Malinowski \email{alexander.malinowski@math.uni-goettingen.de}}

\seealso{
 \code{\link{RFfit}},
 \code{\link{RFempiricalvariogram}},
 \code{\link{RMmodel-class}},
 \code{\link{RMmodelExt-class}}
 \code{\link{plot-method}}
}
 
\examples{
# see RFfit
\dontshow{RFoptions(seed=NA)}
}

\keyword{classes}
