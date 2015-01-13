\name{RFspatialPointsDataFrame-class}
\docType{class}
\alias{RFspatialPointsDataFrame-class}
\alias{RFspatialPointsDataFrame}
%% S4 methods
\alias{RFspDataFrame2conventional,RFspatialPointsDataFrame-method}
\alias{show,RFspatialPointsDataFrame-method}
\alias{coerce,RFspatialPointsDataFrame,RFspatialGridDataFrame-method}
\alias{coerce,RFspatialPointsDataFrame,RFspatialGridDataFrame,ANY-method}
\alias{coerce,RFspatialPointsDataFrame,data.frame-method}
\alias{coerce,SpatialPointsDataFrame,RFspatialPointsDataFrame-method}
\alias{coerce,SpatialPointsDataFrame,RFspatialPointsDataFrame-method}
\alias{[,RFspatialPointsDataFrame-method}
\alias{[,RFspatialPointsDataFrame,ANY,ANY-method}
\alias{[,RFspatialPointsDataFrame,ANY,ANY,ANY-method}
\alias{[<-,RFspatialPointsDataFrame,ANY,ANY,ANY-method}
%% S3 methods
\alias{print.RFspatialPointsDataFrame}
\alias{cbind.RFspatialPointsDataFrame}
\alias{as.matrix.RFspatialPointsDataFrame}
\alias{as.array.RFspatialPointsDataFrame}
\alias{as.vector.RFspatialPointsDataFrame}
\alias{as.data.frame.RFspatialPointsDataFrame}

\title{Class "RFspatialPointsDataFrame" }
\description{ Class for spatial attributes that have spatial or
 spatio-temporal locations (at least of dimension 2) 
 that are not on a grid. Direct extension of class
 \code{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}}
 from the \pkg{sp}-package. See \command{\link{sp2RF}} for an explicite
 transformation.}

\section{Creating Objects}{
 Objects can be created by using the functions
 \code{\link{RFspatialPointsDataFrame}} or
 \code{\link{conventional2RFspDataFrame}} or by calls of the form
 \code{as(x, "RFspatialPointsDataFrame")}, where \code{x} is of class
 \command{\link[=RFspatialPointsDataFrame-class]{RFspatialPointsDataFrame}}. 
 %Ordered full grids are stored instead or unordered non-NA cells; 
}

\section{Slots}{
  \describe{
    \item{\code{.RFparams}:}{list of 2; \code{.RFparams$n} is the number
      of repetitions of the random field contained in the \code{data} slot,
      \code{.RFparams$vdim} gives the dimension of the values of the random
      field, equals 1 in most cases}
    \item{\code{data}:}{object of class
      \code{\link[methods:data.frame-class]{data.frame}}, 
      containing attribute data } 
    \item{\code{coords.nrs}:}{See \code{\link{SpatialPointsDataFrame}}
    }    
    \item{\code{coords}:}{matrix of coordinates (each row is a 
      point); in case of SpatialPointsDataFrame an object of class
      \code{\link[=SpatialPoints-class]{SpatialPoints}} is also allowed
      see \link{SpatialPoints}}
    \item{\code{bbox}:}{matrix specifying the bounding box }
    \item{\code{proj4string}:}{object of class \code{\link[sp:CRS-class]{CRS};
	projection }} 
  }
}

\section{Extends}{
Class
\code{\link[=SpatialPointsDataFrame-class]{SpatialPointsDataFrame}},
directly. 
Class \code{\link[=SpatialPoints-class]{SpatialPoints}}, by class
\code{\link[=SpatialPointsDataFrame]{SpatialPointsDataFrame}}. 
  Class \code{\link[=Spatial-class]{Spatial}},
    by class \code{\link[=SpatialPoints]{SpatialPoints}}.
}

\section{Methods}{
 \describe{
 \item{plot}{\code{signature(obj =
	"RFspatialPointsDataFrame")}: generates nice plots
 of the random field; if \eqn{space-time-dim2}, a two-dimensional
 subspace can be selected using the argument 
 \code{MARGIN}; to get different slices in a third direction, the
 argument \code{MARGIN.slices} can be used;
 for more details see \code{\link{plot-method}} or type
 \code{method?plot("RFspatialPointsDataFrame")}}
 \item{show}{\code{signature(x = "RFspatialPointsDataFrame")}:
 uses the \command{show}-method for class
 \command{\link[sp]{SpatialPointsDataFrame}}.}
 \item{print}{\code{signature(x = "RFspatialPointsDataFrame")}:
 identical to \command{show}-method}
 \item{RFspDataFrame2conventional}{\code{signature(obj =
	"RFspatialPointsDataFrame")}: conversion to a list of
 non-\pkg{sp}-package based objects; the \code{data}-slot
 is converted to an array of dimension \eqn{[1*(vdim>1) +
	space-time-dimension + 1*(n>1)]}}
 \item{coordinates}{\code{signature(x = "RFspatialPointsDataFrame")}:
 returns the coordinates} 
 \item{[}{\code{signature(x = "RFspatialPointsDataFrame")}: selects
 columns of \code{data}-slot; returns an 
 object of class \code{\link{RFspatialPointsDataFrame}}.}
 \item{[<-}{\code{signature(x = "RFspatialPointsDataFrame")}: replaces
 columns of \code{data}-slot; returns an 
 object of class \code{\link{RFspatialPointsDataFrame}}.}
 \item{as}{\code{signature(x = "RFspatialPointsDataFrame")}:
 converts into other formats, only implemented for target class
 \code{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFrame}} } 
 \item{cbind}{\code{signature(...)}: if arguments have identical
   topology, combine their attribute values}
 \item{as.matrix}{\code{signature(x = "RFspatialPointsDataFrame")}:
   converts \code{data}-slot to matrix }
 \item{as.array}{\code{signature(x = "RFspatialPointsDataFrame")}:
   converts \code{data}-slot to array }
 \item{as.vector}{\code{signature(x = "RFspatialPointsDataFrame")}:
   converts \code{data}-slot to vector }
 \item{as.data.frame}{\code{signature(x = "RFspatialPointsDataFrame")}:
   converts \code{data}-slot and coordinates to a data.frame}
 }
}

\section{Details}{
 Note that in the \code{data}-slot, each colums is ordered according to the
 ordering of \code{coordinates(grid)}, the first dimension runs fastest
 and for all BUT the second dimension, coordinate values are in
 ascending order. In the second dimension, coordinate values run from
 high to low. Hence, when converting to conventional formats using
 \command{\link{RFspDataFrame2conventional}} or
 \command{\link{RFspDataFrame2dataArray}}, 
 the data array is re-ordered such that all dimensions are in ascending
 order. \code{as.matrix} does not perform re-ordering.

 Methods \command{summary} and \command{dimensions}
 are defined for the \dQuote{parent}-class
 \command{\link[=RFsp-class]{RFsp}}.
}
\author{Alexander Malinowski \email{alexander.malinowski@math.uni-goettingen.de}}


\seealso{
 \code{\link{RFspatialGridDataFrame-class}}, which is for point
 locations that are on a grid,
 \code{\link{RFpointsDataFrame-class}} which is for one-dimensional
 locations,
 \code{\link[=RFsp-class]{RFsp}},
 \command{\link{sp2RF}}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- cbind(runif(50), runif(50))
f <- RFsimulate(model=RMexp(), x=x, n=3)

str(f)
str(RFspDataFrame2conventional(f))
coordinates(f)[1:25,]
str(f[2]) ## selects second column of data-slot
all.equal(f, cbind(f,f)[1:3]) ## TRUE
try(as(f, "RFspatialGridDataFrame")) # yields error

plot(f, nmax=2)

f2 <- RFsimulate(model=RMwhittle(nu=1.2, scale=10), x=cbind(x,x), n=4)
plot(f2, MARGIN=c(3,4), nmax=2)

f.sp <- RFsimulate(model=RMexp(), x=x, n=3, seed=0)
f.old <- RFsimulate(model=RMexp(), x=x, n=3, spConform=FALSE, seed=0)
all.equal(RFspDataFrame2conventional(f.sp)$data, f.old) ## TRUE

\dontshow{FinalizeExample()}
}

\keyword{classes}
\keyword{print}
\keyword{hplot}
