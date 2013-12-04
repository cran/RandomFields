\name{RFempVariog-class}
\docType{class}
\alias{RFempVariog-class}
\alias{coerce,RFempVariog,list-method}
%\alias{show,RMmodel-method}
%\alias{print.RMmodel}
%\alias{[,RMmodel,ANY,ANY,ANY-method}
%\alias{[<-,RMmodel,ANY,ANY,ANY-method}}
%\alias{plot,RMmodel-method}
%\alias{points.RMmodel}
%\alias{lines.RMmodel}

\title{Class \code{RFempVariog}}
\description{Class for RandomField's representation of empirical variograms 
}

%\section{Creating Objects}{
% Objects are created by calling a function of class
% \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
%}

\section{Slots}{
 \describe{
 \item{\code{centers}:}{the bin centres of the spatial distances} 
 \item{\code{emp.vario}:}{value of the empirical variogram} 
 \item{\code{sd}:}{standard deviation of the variogram cloud within each
 bin} 
 \item{\code{n.bin}:}{number of bins} 
 \item{\code{phi.centers}:}{centres of the bins with respect to the
   (first) angle (for anisotropic empricial variograms only)} 
 \item{\code{theta.centers}:}{centres of the bins with respect to the
 second angle (for anisotropic empricial variograms in 3D only)} 
 \item{\code{T}:}{the bin centres of the time axis} 
 \item{\code{call}:}{language object; the function call by which the
   object was generated}
 \item{\code{var}:}{
   the empirical (overall) variance in the data
 }
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
 }
}

%\section{Details}{
%}

\author{Alexander Malinowski
  \email{alexander.malinowski@math.uni-goettingen.de}}

\seealso{
 \command{\link{RFempiricalvariogram}},
 \command{\link{plot-method}}
}
 
\examples{
set.seed(0)
#  see ?RFempiricalvariogram
}



\keyword{classes}
