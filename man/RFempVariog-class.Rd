\name{RFempVariog-class}
\docType{class}
\alias{RFempVariog-class}
\alias{coerce,RFempVariog,list-method}
\alias{show,RFempVariog-method}
\alias{persp,RFempVariog-method}
\alias{print,RFempVariog-method}
\alias{summary,RFempVariog-method}
\alias{print.RF_empVariog}
\alias{summary.RF_empVariog}
\alias{plot,RFempVariog,missing-method}
\alias{RFplotEmpVariogram}

\title{Class \code{RFempVariog}}
\description{Class for RandomField's representation of empirical variograms 
}

%\section{Creating Objects}{
% Objects are created by calling a function of class
% \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
%}

\usage{
RFplotEmpVariogram(x, model = NULL, nmax.phi = NA, nmax.theta = NA,
                               nmax.T = NA,
                               plot.nbin = TRUE, plot.sd=FALSE, method = "ml",
                               variogram=TRUE,
                               boundaries = TRUE,
                               ...)
\S4method{plot}{RFempVariog,missing}(x, y, ...)
\S4method{persp}{RFempVariog}(x, ...)
}

\arguments{  
  \item{x}{object of class 
    \command{\link[=RFempVariog-class]{RFempVario}}
  }
  \item{y}{unused}
  \item{model}{object of class
    \command{\link[=RMmodel-class]{RMmodel}} or
    \code{class(x)=="RFempVario"} or \code{class(x)=="RFfit"}; a list of
    covarianve or variogram models that are to be plotted into the same
    plot as the empirical variogram (and the fitted models)}
 \item{nmax.phi}{integer; only for \code{class(x)=="RFempVario"};
    the maximal number of bins of angle phi that are to be
    plotted}
  \item{nmax.theta}{integer; only for \code{class(x)=="RFempVario"};
    the maximal number of bins of angle theta that are to be
    plotted} 
  \item{nmax.T}{integer; only for \code{class(x)=="RFempVario"};
    the maximal number of different time bins that are to be
    plotted} 
  \item{plot.nbin}{logical; only for \code{class(x)=="RFempVario"};
    indicates whether the number of pairs per bin are to be plotted}
  \item{plot.sd}{logical; only for \code{class(x)=="RFempVario"};
    indicates whether the calculated standard deviation (\code{x@sd}) is
    to be plotted (in form of arrows of length +-1*sd)}
  \item{method}{character. Currently restricted to \code{"ml"} for
    maximum-likelihood method.
    }
  \item{variogram}{logical; This argument should currently not be set
    by the user. If \code{TRUE} then the empirical variogram
    is plotted, else an estimate for the covariance function} 
  \item{boundaries}{logical; only for \code{class(x)=="RFempVario"} and
    the anisotropic case where \code{model} is given.
    As the empirical variogram is calculated on a sector of angles,
    no exact variogram curve corresponds to the mean values in this
    sector.
    If \code{boundaries=TRUE} the values of the variogram on the sector
    boundaries are plotted. If \code{FALSE} some kind of mean
    model values are plotted. Neither the boundaries may contain
    the values of emprical variogram nor does the mean values
    need to be close the emprical variogram.  
  }
  \item{...}{arguments to be passed to methods; mainly graphical
    arguments.
    %see Details.
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
 \item{\code{coordunits}:}{string giving the units of the coordinates,
   see also option \code{coordunits} of \command{\link{RFoptions}}.
 }
 \item{\code{varunits}:}{string giving the units of the variables,
   see also option \code{varunits} of \command{\link{RFoptions}}.
 }
 \item{\code{call}:}{language object; the function call by which the
   object was generated}
 \item{\code{method}:}{integer; variogram (0), covariance (2), madogram (4) }
}
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "RFempVariog")}: gives a plot of the
      empirical variogram, for more details see
      \command{\link{plot-method}}.}
    \item{plot}{\code{signature(x = "RFempVario", y = "missing")}}{Gives nice
      plots of the empirical variogram; handles binning in up to three
      space-dimensions and a time-dimension, where the empirical variogram
      is plotted along lines which are directed according to the
      angle-centers given in \code{x@phi.centers} and
      \code{x@theta.centers}; arbitrary theoretical model curves
      can be added to the plot by using the argument \code{model}.
      If no bins are given, i.e. (\code{x@bin=NULL}),
      \command{\link[graphics]{image}}-plots are generated.}
    \item{as}{\code{signature(x = "RFempVariog")}:
      converts into other formats, only implemented for target class
      \code{\link[methods:list-class]{list}}.}
    \item{show}{\code{signature(x = "RFfit")}: returns the structure
      of \code{x}}
    \item{persp}{\code{signature(obj =
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

\author{Alexander Malinowski \email{Alexander.Malinowski@web.de}}
  

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
