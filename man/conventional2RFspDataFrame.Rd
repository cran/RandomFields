\name{conventional2RFspDataFrame}
\alias{conventional2RFspDataFrame}
\title{Coercion to class 'RFsp' objects}

\description{
 Generate an object of class \code{\link[=RFsp-class]{RFsp}} from conventional objects
}

\usage{
conventional2RFspDataFrame(data, coords=NULL, gridTopology=NULL, n=1,
vdim=1)
}

\arguments{
 \item{data}{array; of dimension \code{c(vdim, space-time-dim, n)};
 contains the values of the random fields}
 \item{coords}{matrix of coordinates}
 \item{gridTopology}{3-row-matrix or of class \code{\link[sp:GridTopology-class]{GridTopology}}; specifies
 the grid vectors; either
 \code{coords} or \code{gridTopology} must be \code{NULL}}
 \item{n}{number of iid copies of the random fields, default is 1}
 \item{vdim}{number of dimensions of the values of the random field,
 default is 1} 
}

%\details{}

\value{object of class
 \command{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFrame}},
 \command{\link[=RFspatialPointsDataFrame-class]{RFspatialPointsDataFrame}},
 \command{\link[=RFgridDataFrame-class]{RFgridDataFrame}} or
 \command{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}}
}

\author{Alexander Malinowski
  \email{alexander.malinowski@math.uni-goettingen.de}}

\examples{
set.seed(0)
x <- 1:20
z <- RFsimulate(RMexp(), x, spC=FALSE)
z2 <- conventional2RFspDataFrame(z, coord=x)
Print(z, z2)
}

\keyword{spatial}
