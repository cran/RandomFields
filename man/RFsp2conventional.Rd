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

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\seealso{
  \command{\link{RFspatialGridDataFrame}},
  \command{\link{RFspatialPointsDataFrame}}
}

\examples{
## for examples, see the links in section 'See Also'
\dontshow{RFoptions(seed=NA)}
}

\keyword{spatial}
\keyword{models}
