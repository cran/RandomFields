\name{RFgridDataFrame-class}
\docType{class}
\alias{RFgridDataFrame-class}
\alias{RFgridDataFrame}
%% S4 methods
\alias{RFspDataFrame2conventional,RFgridDataFrame-method}
\alias{RFspDataFrame2dataArray,RFgridDataFrame-method}
\alias{show,RFgridDataFrame-method}
\alias{coerce,RFgridDataFrame,RFpointsDataFrame-method}
\alias{coerce,RFgridDataFrame,RFpointsDataFrame,ANY-method}
\alias{[,RFgridDataFrame-method}
\alias{[,RFgridDataFrame,ANY,ANY-method}
\alias{[,RFgridDataFrame,ANY,ANY,ANY-method}
\alias{[<-,RFgridDataFrame,ANY,ANY,ANY-method}
%% S3 methods
\alias{print.RFgridDataFrame}
\alias{cbind.RFgridDataFrame}
\alias{range.RFgridDataFrame}
\alias{hist.RFgridDataFrame}
\alias{as.data.frame.RFgridDataFrame}
\alias{as.matrix.RFgridDataFrame}
\alias{as.array.RFgridDataFrame}
\alias{as.vector.RFgridDataFrame}
\alias{as.data.frame.RFgridDataFrame}

\title{Class \code{RFgridDataFrame} }
\description{ Class for attributes in one-dimensional space.}

\section{Creating Objects}{
 Objects can be created by using the functions
 \code{\link{RFgridDataFrame}} or
 \code{\link{conventional2RFspDataFrame}} or by calls of the form
 \code{as(x, "RFgridDataFrame")}, where \code{x} is of class
 \command{\link[=RFgridDataFrame-class]{RFgridDataFrame}}. 
 %Ordered full grids are stored instead or unordered non-NA cells; 
}

\section{Slots}{
 \describe{
   \item{\code{.RFparams}:}{list of up to 5 elements;
     \itemize{       \item{n} is the number
       of repetitions of the random field contained in the \code{data} slot
       \item \code{vdim} gives the dimension of the values of
       the random 
       field, equals 1 in most cases
       \item \code{has.variance}
       indicates whether information on the variance is available,
       \item \code{coordunits} gives the names of the
       units for the coordinates
       \item \code{varunits} gives the names of the
     units for the variables
   }
}
 \item{\code{data}:}{object of class \code{\link[methods:data.frame-class]{data.frame}},
 containing attribute data } 
 \item{\code{grid}:}{object of class \code{\link[sp:GridTopology-class]{GridTopology}}.}
 }
}


\section{Methods}{
 \describe{
 \item{plot}{\code{signature(obj =
	"RFgridDataFrame")}: generates nice plots
 of the random field; if \eqn{space-time-dim2}, a two-dimensional
 subspace can be selected using the argument 
 \code{MARGIN}; to get different slices in a third direction, the
 argument \code{MARGIN.slices} can be used;
 for more details see \code{\link{plot-method}} or type
 \code{method?plot("RFgridDataFrame")}}
 \item{show}{\code{signature(x = "RFgridDataFrame")}:
 uses the \command{show}-method for class \code{\link[=SpatialGridDataFrame-class]{SpatialGridDataFrame}}.}
 \item{print}{\code{signature(x = "RFgridDataFrame")}:
 identical to \command{show}-method}
 \item{RFspDataFrame2conventional}{\code{signature(obj =
	"RFgridDataFrame")}: conversion to a list of
 non-\pkg{sp}-package based objects; the \code{data}-slot
 is converted to an array of dimension \eqn{[1*(vdim>1) +
	space-time-dimension + 1*(n>1)]}}
 \item{coordinates}{\code{signature(x = "RFgridDataFrame")}:
 returns the coordinates} 
 \item{[}{\code{signature(x = "RFgridDataFrame")}: selects
 columns of \code{data}-slot; returns an 
 object of class \code{\link[=RFgridDataFrame-class]{RFgridDataFrame}}.}
 \item{[<-}{\code{signature(x = "RFgridDataFrame")}: replaces
 columns of \code{data}-slot; returns an 
 object of class \code{\link[=RFgridDataFrame-class]{RFgridDataFrame}}.}
 \item{as}{\code{signature(x = "RFgridDataFrame")}:
 converts into other formats, only implemented for target class
 \code{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}} } 
 \item{cbind}{\code{signature(...)}: if arguments have identical
 topology, combine their attribute values}
 \item{range}{\code{signature(x = "RFgridDataFrame")}: returns the range
 }
 \item{hist}{\code{signature(x = "RFgridDataFrame")}: plots histogram
 }
 \item{as.matrix}{\code{signature(x = "RFgridDataFrame")}: converts
          \code{data}-slot to matrix}
 \item{as.array}{\code{signature(x = "RFgridDataFrame")}: converts
          \code{data}-slot to array}
 \item{as.vector}{\code{signature(x = "RFgridDataFrame")}: converts
   \code{data}-slot to vector}
 \item{as.data.frame}{\code{signature(x = "RFgridDataFrame")}:
   converts \code{data}-slot and coordinates to a data.frame}
 }
}

\section{Details}{
 Methods \command{summary} and \command{dimensions}
 are defined for the \dQuote{parent}-class
 \command{\link[=RFsp-class]{RFsp}}.
}
\author{Alexander Malinowski \email{Alexander.Malinowski@web.de}}


\seealso{
 \code{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFrame}},
 which is for point 
 locations in higher dimensional spaces,
 \code{\link{RFpointsDataFrame-class}} which is for one-dimensional
 arbitray locations,
 \code{\link[=RFsp-class]{RFsp}} 
}


\keyword{classes}
\keyword{print}
\keyword{hplot}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- seq(0,10,length=100)
f <- RFsimulate(model=RMgauss(), x=x, n=3)

str(f)
str(RFspDataFrame2conventional(f))
head(coordinates(f))
str(f[2]) ## selects second column of data-slot
all.equal(f, cbind(f,f)[1:3]) ## TRUE

plot(f, nmax=2)
\dontshow{FinalizeExample()}

}
