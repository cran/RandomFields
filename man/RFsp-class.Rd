\name{RFsp-class}
\docType{class}
\alias{RFsp-class}
\alias{RFspatialDataFrame}
\alias{RFdataFrame}
%% S4 methods
\alias{summary,RFsp-method}
\alias{variance,RFsp-method}
\alias{dimensions,RFsp-method}
\alias{dimensions,RFdataFrame-method}
\alias{dimensions,RFspatialDataFrame-method}
\alias{isGridded}
\alias{isGridded,RFsp-method}
\alias{RFspDataFrame2conventional}
\alias{RFspDataFrame2conventional,RFsp-method}
\alias{RFspDataFrame2dataArray}
\alias{RFspDataFrame2dataArray,RFsp-method}
\alias{[,RFsp-method}
\alias{[,RFsp,ANY,ANY-method}
\alias{[,RFsp,ANY,ANY,ANY-method}
\alias{[<-,RFsp-method}
\alias{[<-,RFsp,ANY,ANY-method}
\alias{[<-,RFsp,ANY,ANY,ANY-method}
\alias{RFspatialDataFrame-class}
\alias{RFspatialDataFrame}
\alias{RFdataFrame-class}
\alias{RFdataFrame}

\title{Class \code{RFsp}}

\description{ \code{"RFsp"} is a virtual class which contains the four classes
 \command{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFrame}} 
 (data on a full grid and \eqn{space-time-dimension \ge
 2}{space-time-dimension > 2}), 
 \command{\link[=RFspatialPointsDataFrame-class]{RFspatialPointsDataFrame}}
 (data not on a grid and \eqn{space-time-dimension \ge
 2}{space-time-dimension > 2}), 
 \command{\link[=RFgridDataFrame-class]{RFgridDataFrame}}
 (data on a full grid and \eqn{space-time-dimension = 1}), 
 \command{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}}
 (data not on a grid spaced and \eqn{space-time-dimension = 1})


 The first two class subclasses are summarized in
 \code{"RFspatialDataFrame"} whilst the latter two are
 summarized in \code{"RFdataFrame"}.
 
% The first two classes are summarized in the class
% \code{RFspatialDataFrame} the last two classes in the class
%\code{ RFdataFrame}.
 
}

\section{Objects from the Class}{
 are never to be generated; only derived classes can be meaningful }

%\section{Slots}{
% \describe{
% only for sub-classes
% \command{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFram%e}}
% and \command{\link[=RFspatialPointsDataFrame-class]{RFspatialPoints%DataFrame}}:
% \item{\code{bbox}:}{Object of class \code{"matrix"}; 2-column matrix
%	holding the minimum in first and maximum in second column for the
%	x-coordinate (first row), y-coordinate (second row) and optionally,
%	for points and grids only, further coordinates. The constructed 
%	Spatial object will be invalid if any bbox values are NA or
%	infinite. The column names must be \code{c("min", "max")}} 
 % \item{\code{proj4string}:}{Object of class \code{"CRS"};
%	holding a valid proj4 string, which can be used for unprojecting
%	or reprojecting coordinates; it is initialised to NA. Other strings
%	are checked for validity in the rgdal package, but attempts to
%	assign a string containing \dQuote{longlat} to data extending beyond
%	longitude [-180, 360] or lattitude [-90, 90] will be stopped. }
% }
%}
\section{Methods}{
 \describe{
% \item{bbox}{\code{signature(obj = "Spatial")}: retrieves the bbox element }
\item{summary}{\code{signature(obj = "RFsp")}: returns a summary of the
 object; uses or imitates summary method of class
 \code{\link[sp:Spatial-class]{Spatial}} from the \pkg{sp}-package}
\item{dimensions}{\code{signature(obj = "RFsp")}: retrieves the number
 of spatial or spatio-temporal dimensions spanned }
\item{isGridded}{\code{signature(obj = "RFsp")}: logical, tells whether
  the data is on a regular spatial grid }
\item{RFspDataFrame2dataArray}{\code{signature(obj = "RFsp")}:
  transforms \code{RFsp} objects to array
}
\item{RFspDataFrame2conventional}{\code{signature(obj = "RFsp")}:
  transforms \code{RFsp} objects to a list with additional information
}
\item{[}{\code{signature(obj = "RFsp")}: selects columns of the
 \code{data}-slot, while all other slots are kept unmodified}
\item{[<-}{\code{signature(obj = "RFsp")}: replaces columns of the
 \code{data}-slot, while all other slots are kept unmodified}
% \item{plot}{\code{signature(x = "Spatial", y = "missing")}: plot method
%	for spatial objects; does nothing but setting up a plotting region choosing
%	a suitable aspect if not given(see below), colouring the plot background using either a bg= argument or par("bg"), and possibly drawing axes. }
\item{variance}{\code{signature(object = "RFsp")}: returns the kriging
 variance if available}
 }
}

%\usage{
%Spatial(bbox, proj4string = CRS(as.character(NA)))
%}

%\arguments{
%\item{bbox}{a bounding box matrix}
%\item{proj4string}{a CRS object}
%}

\author{Alexander Malinowski
  \email{Alexander.Malinowski@web.de},
  Martin Schlather \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}
}


\section{Warning }{this class is not useful in itself, but the above
 mentioned classes in this package derive from it }

\seealso{
 \command{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFrame}},
 \command{\link[=RFspatialPointsDataFrame-class]{RFspatialPointsDataFrame}},
 \command{\link[=RFgridDataFrame-class]{RFgridDataFrame}},
 \command{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}},
 \command{\link{sp2RF}}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
## to do
\dontshow{FinalizeExample()}
}

\keyword{classes}
