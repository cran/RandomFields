\name{RFpointsDataFrame-class}
\docType{class}
\alias{RFpointsDataFrame-class}
\alias{RFpointsDataFrame}
%% S4 methods
\alias{RFspDataFrame2conventional,RFpointsDataFrame-method}
\alias{show,RFpointsDataFrame-method}
\alias{coerce,RFpointsDataFrame,RFgridDataFrame-method}
\alias{coerce,RFpointsDataFrame,RFgridDataFrame,ANY-method}
\alias{[,RFpointsDataFrame-method}
\alias{[,RFpointsDataFrame,ANY,ANY-method}
\alias{[,RFpointsDataFrame,ANY,ANY,ANY-method}
\alias{[<-,RFpointsDataFrame-method}
\alias{[<-,RFpointsDataFrame,ANY,ANY,ANY-method}
\alias{dimensions,RFpointsDataFrame-method}
%\alias{isGridded,RFpointsDataFrame-method}
\alias{variance,RFpointsDataFrame-method}
%% S3 methods
\alias{print.RFpointsDataFrame}
\alias{cbind.RFpointsDataFrame}
\alias{range.RFpointsDataFrame}
\alias{hist.RFpointsDataFrame}
\alias{as.matrix.RFpointsDataFrame}
\alias{as.array.RFpointsDataFrame}
\alias{as.vector.RFpointsDataFrame}
\alias{as.data.frame.RFpointsDataFrame}

\title{Class \code{RFpointsDataFrame}}
\description{ Class for attributes in one-dimensional space 
 that are not on a grid. }

\usage{
\S4method{RFspDataFrame2conventional}{RFpointsDataFrame}(obj) 
}

\arguments{
  \item{obj}{an \code{RFspatialPointsDataFrame} object}
}




\section{Creating Objects}{
 Objects can be created by using the functions
 \code{\link{RFpointsDataFrame}} or
 \code{\link{conventional2RFspDataFrame}} or by calls of the form
 \code{as(x, "RFpointsDataFrame")}, where \code{x} is of class
 \command{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}}. 
 %Ordered full grids are stored instead or unordered non-NA cells; 
}

\section{Slots}{
 \describe{
 \item{\code{data}:}{object of class \code{\link[methods:data.frame-class]{data.frame}},
 containing attribute data } 
 \item{\code{coords}:}{n-times-1 matrix of coordinates (each row is a 
 point)}
 \item{\code{.RFparams}:}{list of 2; \code{.RFparams$n} is the number
 of repetitions of the random field contained in the \code{data} slot,
 \code{.RFparams$vdim} gives the dimension of the values of the random
 field, equals 1 in most cases}
 }
}


\section{Methods}{
 \describe{
 \item{plot}{\code{signature(obj =
	"RFpointsDataFrame")}: generates nice plots
 of the random field; if \eqn{space-time-dim2}, a two-dimensional
 subspace can be selected using the argument 
 \code{MARGIN}; to get different slices in a third direction, the
 argument \code{MARGIN.slices} can be used;
 for more details see \code{\link{plot-method}} or type
 \code{method?plot("RFpointsDataFrame")}.}
 \item{show}{\code{signature(x = "RFpointsDataFrame")}:
 uses the \command{show}-method for class  \code{\link[=SpatialPointsDataFrame-class]{SpatialPointsDataFrame}}.}
 \item{print}{\code{signature(x = "RFpointsDataFrame")}:
 identical to \command{show}-method}
 \item{RFspDataFrame2conventional}{\code{signature(obj =
	"RFpointsDataFrame")}: conversion to a list of
 non-\pkg{sp}-package based objects; the \code{data}-slot
 is converted to an array of dimension \eqn{[1*(vdim>1) +
	space-time-dimension + 1*(n>1)]}}
 \item{coordinates}{\code{signature(x = "RFpointsDataFrame")}:
 returns the coordinates} 
 \item{[}{\code{signature(x = "RFpointsDataFrame")}: selects
 columns of \code{data}-slot; returns an 
 object of class \code{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}}.}
 \item{[<-}{\code{signature(x = "RFpointsDataFrame")}: replaces
 columns of \code{data}-slot; returns an 
 object of class \code{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}}.}
 \item{as}{\code{signature(x = "RFpointsDataFrame")}:
 converts into other formats, only implemented for target class
 \code{\link[=RFgridDataFrame-class]{RFgridDataFrame}} } 
 \item{cbind}{\code{signature(...)}: if arguments have identical
   topology, combine their attribute values}
 \item{range}{\code{signature(x = "RFpointsDataFrame")}: returns the range
 }
 \item{hist}{\code{signature(x = "RFpointsDataFrame")}: plots histogram
 }
 \item{as.matrix}{\code{signature(x = "RFpointsDataFrame")}: converts
          \code{data}-slot to matrix}
 \item{as.array}{\code{signature(x = "RFpointsDataFrame")}: converts
          \code{data}-slot to array}
 \item{as.vector}{\code{signature(x = "RFpointsDataFrame")}: converts
          \code{data}-slot to vector}
 \item{as.data.frame}{\code{signature(x = "RFpointsDataFrame")}:
   converts \code{data}-slot and coordinates to a data.frame}
 }
}

\section{Details}{
 Methods \command{summary} and \command{dimensions}
 are defined for the \dQuote{parent}-class
 \command{\link[=RFsp-class]{RFsp}}.
}
\author{Alexander Malinowski, \martin}


\seealso{
 \code{\link[=RFspatialPointsDataFrame-class]{RFspatialPointsDataFrame}}, which is for point
 locations in higher dimensional spaces,
 \code{\link{RFpointsDataFrame-class}} which is for one-dimensional
 locations on a grid,
 \code{\link[=RFsp-class]{RFsp}} 
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- runif(100)
f <- RFsimulate(model=RMexp(), x=x, n=3)

str(f)
str(RFspDataFrame2conventional(f))
head(coordinates(f))
str(f[2]) ## selects second column of data-slot
all.equal(f, cbind(f,f)[1:3]) ## TRUE

plot(f, nmax=2)
\dontshow{FinalizeExample()}}

\keyword{classes}
\keyword{print}
\keyword{hplot}
