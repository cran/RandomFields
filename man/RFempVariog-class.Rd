\name{RFempVariog-class}
\docType{class}
\alias{RFempVariog-class}
\alias{coerce,RFempVariog,list-method}
\alias{show,RFempVariog-method}
\alias{persp,RFempVariog-method}
\alias{print,RFempVariog-method}
\alias{summary,RFempVariog-method}
%\alias{[,RMmodel,ANY,ANY,ANY-method}
%\alias{[<-,RMmodel,ANY,ANY,ANY-method}}
%\alias{plot,RMmodel-method}
%\alias{points.RMmodel}
%\alias{lines.RMmodel}
\alias{print.RF_empVariog}
\alias{summary.RF_empVariog}


\title{Class \code{RFempVariog}}
\description{Class for RandomField's representation of empirical variograms 
}

%\section{Creating Objects}{
% Objects are created by calling a function of class
% \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
%}

%\usage{
%summary.RF_empVariog(object, ...) 
%%print.RF_empVariog(x, ...)
%}

\arguments{  
  \item{object, x, ...}{
    see the respective generic function;
    The argument
  }
}


\section{Slots}{
 \describe{
 \item{\code{centers}:}{the bin centres of the spatial distances} 
 \item{\code{emp.vario}:}{value of the empirical variogram} 
 \item{\code{var}:}{
   the empirical (overall) variance in the data
 }
 \item{\code{sd}:}{standard deviation of the variogram cloud within each
 bin} 
 \item{\code{n.bin}:}{number of bins} 
 \item{\code{phi.centers}:}{centres of the bins with respect to the
   (first) angle (for anisotropic empricial variograms only)} 
 \item{\code{theta.centers}:}{centres of the bins with respect to the
 second angle (for anisotropic empricial variograms in 3D only)} 
 \item{\code{T}:}{the bin centres of the time axis} 
 \item{\code{vdim}:}{the multivariate dimension}
 \item{\code{coord.units}:}{string giving the units of the coordinates,
   see also option \code{coord_units} of \command{\link{RFoptions}}.
 }
 \item{\code{variab.units}:}{string giving the units of the variables,
   see also option \code{variab_units} of \command{\link{RFoptions}}.
 }
 \item{\code{call}:}{language object; the function call by which the
   object was generated}
}
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "RFempVariog")}: gives a plot of the
      empirical variogram, for more details see
      \command{\link{plot-method}}.}
    \item{as}{\code{signature(x = "RFempVariog")}:
      converts into other formats, only implemented for target class
      \code{\link[methods:list-class]{list}}.}
    \item{show}{\code{signature(x = "RFfit")}: returns the structure
      of \code{x}}
    \item{persp}{code{signature(obj =
	"RFempVario")}: generates nice \command{\link[graphics]{persp}} plots
    }
    \item{print}{\code{signature(x = "RFfit")}: identical with
      \command{show}-method}
   \item{summary}{provides a summary}
  }
}


\section{Details}{
  \command{print} returns also an invisible list that is convenient to access.
}

\author{Alexander Malinowski
  \email{alexander.malinowski@math.uni-goettingen.de}}

\seealso{
 \command{\link{RFempiricalvariogram}},
 \command{\link{plot-method}}
}
 
\examples{
#  see 'RFempiricalvariogram'
\dontshow{FinalizeExample()}
}


\keyword{classes}
\keyword{print}
\keyword{hplot}
