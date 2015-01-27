\name{RFsp2conventional}
\alias{RFspDataFrame2conventional}
\alias{RFspDataFrame2dataArray}
\title{coerce RFsp-class to conventional format}
\description{
 coerce RFsp-class to conventional format
}
\usage{
RFspDataFrame2conventional(obj)

RFspDataFrame2dataArray(obj)
}
\arguments{
 \item{obj}{object of class \code{\link[=RFsp-class]{RFsp}}.}
}
\value{
 \command{\link{RFspDataFrame2conventional}} returns a list;
 \command{\link{RFspDataFrame2dataArray}} returns an array containing the
 data-slot of \code{obj};
 
}

\author{Alexander Malinowski \email{Alexander.Malinowski@web.de}}

\seealso{
  \command{\link{RFspatialGridDataFrame}},
  \command{\link{RFspatialPointsDataFrame}}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

n <- 3
x <- GridTopology(cellcentre.offset=c(0, 0),
 cellsize=c(1, 0.2),
 cells.dim=c(10, 30))
f <- RFsimulate(model=RMexp(), x=x, n=n)
str(f)
str(RFspDataFrame2conventional(f))

\dontshow{FinalizeExample()}
}

\keyword{spatial}
\keyword{models}
