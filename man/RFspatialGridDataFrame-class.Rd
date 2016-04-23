\name{RFspatialGridDataFrame-class}
\docType{class}
\alias{RFspatialGridDataFrame-class}
\alias{RFspatialGridDataFrame}
%% S4 methods
\alias{RFspDataFrame2conventional,RFspatialGridDataFrame-method}
\alias{RFspDataFrame2dataArray,RFspatialGridDataFrame-method}
\alias{show,RFspatialGridDataFrame-method}
\alias{coerce,RFspatialGridDataFrame,data.frame-method}
\alias{coerce,SpatialGridDataFrame,RFspatialGridDataFrame-method}
\alias{coerce,RFspatialGridDataFrame,RFspatialPointsDataFrame-method}
\alias{coerce,RFspatialGridDataFrame,RFspatialPointsDataFrame,ANY-method}
\alias{GridTopology2gridVectors,GridTopology-method}
\alias{GridTopology2gridVectors,matrix-method}
%\alias{coordinates,RFspatialGridDataFrame-method} only for
 %"RFpointsDataFrame" and "RFgridDataFrame" defined
\alias{[,RFspatialGridDataFrame-method}
\alias{[,RFspatialGridDataFrame,ANY,ANY-method}
\alias{[,RFspatialGridDataFrame,ANY,ANY,ANY-method}
\alias{[<-,RFspatialGridDataFrame-method}
\alias{[<-,RFspatialGridDataFrame,ANY,ANY,ANY-method}
\alias{dimensions,RFspatialGridDataFrame-method}
\alias{isGridded,RFspatialGridDataFrame-method}
\alias{variance,RFspatialGridDataFrame-method}
%% S3 methodsgr
\alias{print.RFspatialGridDataFrame}
\alias{cbind.RFspatialGridDataFrame}
\alias{range.RFspatialGridDataFrame}
\alias{hist.RFspatialGridDataFrame}
\alias{as.matrix.RFspatialGridDataFrame}
\alias{as.array.RFspatialGridDataFrame}
\alias{as.vector.RFspatialGridDataFrame}
\alias{as.data.frame.RFspatialGridDataFrame}


\title{Class "RFspatialGridDataFrame" }
\description{ Class for spatial attributes that have spatial or
  spatio-temporal locations (at least of dimension 2) 
  on a (full) regular grid. Direct extension of class
  \command{\link[sp:SpatialGridDataFrame-class]{SpatialGridDataFrame}
    from the \pkg{sp}-package. See \command{\link{sp2RF}} for an explicite
    transformation.}
}

\usage{
\S4method{RFspDataFrame2conventional}{RFspatialGridDataFrame}(obj, data.frame=FALSE) 
}

\arguments{
  \item{obj}{an \code{RFspatialGridDataFrame} object}
  \item{data.frame}{logical. If \code{TRUE} a \code{data.frame} is returned.
  }
}

\section{Creating Objects}{
 Objects can be created by using the functions
 \code{\link{RFspatialGridDataFrame}} or
 \code{\link{conventional2RFspDataFrame}} or by calls of the form
 \code{as(x, "RFspatialGridDataFrame")}, where \code{x} is of class
 \command{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFrame}}. 
 %Ordered full grids are stored instead or unordered non-NA cells; 
}

\section{Slots}{
  \describe{
    \item{\code{.RFparams}:}{list of 2; \code{.RFparams$n} is the number
      of repetitions of the random field contained in the \code{data} slot,
      \code{.RFparams$vdim} gives the dimension of the values of the random
      field, equals 1 in most cases}
    \item{\code{data}:}{object of class \code{\link{data.frame}},
      containing attribute data } 
    \item{\code{grid}:}{object of class
      \command{\link[=GridTopology-class]{GridTopology}}; grid
      parameters } 
    %\item{\code{grid.index}:}{see \link{SpatialPixels-class}; this slot
    % is of zero length for this class, as the grid is full }
    %\item{\code{coords}:}{matrix containing a subset of
    % \code{grid} }
    %numeric matrix or data.frame with coordinates
    % (each row is a 
    % point); in case of SpatialPointsDataFrame an object of class
    % SpatialPoints-class is also allowed
    % 
    % see \link{SpatialPoints}; points slot which is
    % not actually filled with all coordinates (only with min/max)
    \item{\code{bbox}:}{matrix specifying the bounding box }
    \item{\code{proj4string}:}{object of class \code{\link[sp:CRS-class]{CRS}};
      projection }
  }
}

\section{Extends}{
  Class \code{"SpatialGridDataFrame"}, directly.
  Class \code{"SpatialGrid"}, by class \code{"SpatialGridDataFrame"}.
  Class \code{"Spatial"}, by class \code{"SpatialGrid"}.
}

\section{Methods}{
 \describe{
\item{contour}{code{signature(obj =
    "RFspatialGridDataFrame")}:
  generates \command{\link[graphics]{contour}} plots
  }
 \item{plot}{\code{signature(obj =
	"RFspatialGridDataFrame")}: generates nice \command{image} plots
 of the random field; if \eqn{space-time-dim2}, a two-dimensional
 subspace can be selected using the argument 
 \code{MARGIN}; to get different slices in a third direction, the
 argument \code{MARGIN.slices} can be used;
 for more details see \code{\link{plot-method}} or type
 \code{method?plot("RFspatialGridDataFrame")}}

\item{persp}{code{signature(obj =
	"RFspatialGridDataFrame")}: generates \command{\link[graphics]{persp}} plots
}

 \item{show}{\code{signature(x = "RFspatialGridDataFrame")}:
 uses the \command{show}-method for class
 \command{\link[sp]{SpatialGridDataFrame}}.}
 \item{print}{\code{signature(x = "RFspatialGridDataFrame")}:
   identical to \command{show}-method}
 
 \item{RFspDataFrame2conventional}{\code{signature(obj =
	"RFspatialGridDataFrame")}: conversion to a list of
 non-\pkg{sp}-package based objects; the \code{data}-slot
 is converted to an array of dimension \eqn{[1*(vdim>1) +
	space-time-dimension + 1*(n>1)]}; the \code{grid}-slot is
 converted to a 3-row matrix; the grid definition of a possible
 time-dimension becomes a seperate list element}
 \item{RFspDataFrame2dataArray}{\code{signature(obj =
	"RFspatialGridDataFrame")}: conversion of the \code{data}-slot
 to an array of dimension \eqn{[space-time-dimension + 2]}, where the
 space-time-dimensions run fastest, and \eqn{vdim} and \eqn{n} are
 the last two dimensions}
 \item{coordinates}{\code{signature(x = "RFspatialGridDataFrame")}:
 calculates the coordinates from grid definition} 
 \item{[}{\code{signature(x = "RFspatialGridDataFrame")}: selects
 columns of \code{data}-slot; returns an 
 object of class \code{\link{RFspatialGridDataFrame}}.}
 \item{[<-}{\code{signature(x = "RFspatialGridDataFrame")}: replaces
 columns of \code{data}-slot; returns an 
 object of class \code{\link{RFspatialGridDataFrame}}.}
 \item{as}{\code{signature(x = "RFspatialGridDataFrame")}:
 converts into other formats, only implemented for target class
 \code{\link[=RFspatialPointsDataFrame-class]{RFspatialPointsDataFrame}} } 
 \item{cbind}{\code{signature(...)}: if arguments have identical
 topology, combine their attribute values}
 \item{range}{\code{signature(x = "RFspatialGridDataFrame")}: returns the range
 }
 \item{hist}{\code{signature(x = "RFspatialGridDataFrame")}: plots histogram
 }
 \item{as.matrix}{\code{signature(x = "RFspatialGridDataFrame")}:
 converts \code{data}-slot to matrix }
 \item{as.array}{\code{signature(x = "RFspatialGridDataFrame")}:
 converts \code{data}-slot to array }
 \item{as.vector}{\code{signature(x = "RFspatialGridDataFrame")}:
 converts \code{data}-slot to vector }
\item{as.data.frame}{\code{signature(x = "RFspatialGridDataFrame")}:
 converts \code{data}-slot and coordinates to a data.frame}
  }
}

\section{Details}{
 Note that in the \code{data}-slot, each colums is ordered according to the
 ordering of \code{coordinates(grid)}, the first dimension runs fastest
 and for all BUT the second dimension, coordinate values are in
 ascending order. In the second dimension, coordinate values run from
 high to low. Hence, when converting to conventional formats using
 \command{RFspDataFrame2conventional} or \command{RFspDataFrame2dataArray},
 the data array is re-ordered such that all dimensions are in ascending
 order. \code{as.matrix} does not perform re-ordering.

 Methods \command{summary}, \command{dimensions} and
 \command{isGridded} are defined for the \dQuote{parent}-class
 \command{\link[=RFsp-class]{RFsp}}.
}
\author{Alexander Malinowski,
  Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de/de/publications/software}}



\seealso{
 \code{\link{RFspatialPointsDataFrame-class}}, which is for point
 locations that are not on a grid,
 \code{\link{RFgridDataFrame-class}} which is for one-dimensional
 locations,
 \code{\link[=RFsp-class]{RFsp}},
 \command{\link{sp2RF}}
}


\keyword{classes}
\keyword{print}
\keyword{hplot}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
%   library(RandomFields, lib="~/TMP")

n <- 3

x <- GridTopology(cellcentre.offset=c(0, 0),
 cellsize=c(1, 0.2),
 cells.dim=c(10, 30))
f <- RFsimulate(model=RMexp(), x=x, n=n)

str(f)
str(RFspDataFrame2conventional(f))
str(RFspDataFrame2dataArray(f))
coordinates(f)[1:25,]
str(f[2]) ## selects second column of data-slot
all.equal(f, cbind(f,f)[1:3]) ## TRUE
str(as(f, "RFspatialPointsDataFrame"))

plot(f, nmax=2)

steps <- c(10, 1, 10, 10)
\dontshow{if(!interactive()){warning("modified <steps>"); steps<-c(2,1,2,2)}}%ok
x2 <- rbind(c(0, 0, 0, 0),
 c(1, 0.2, 2, 5),
 steps)
scale <- 10
\dontshow{if(!interactive()){warning("modified <scale>");scale<-0.1}}%ok
f2 <- RFsimulate(model=RMwhittle(nu=1.2, scale=scale), x=x2, n=n,
                 grid = TRUE)
plot(f2, MARGIN=c(3,4), MARGIN.slices=1, n.slices=6, nmax=2)

f.sp <- RFsimulate(model=RMexp(), x=x, n=n, seed=0)
f.old <- RFsimulate(model=RMexp(), x=x, n=n, spConform=FALSE, seed=0)
all.equal(RFspDataFrame2conventional(f.sp)$data, f.old) ## TRUE
\dontshow{FinalizeExample()}

}
