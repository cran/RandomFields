\name{RMtrafo}
\alias{RMtrafo}
\alias{RFearth2cartesian}
\alias{RFearth2dist}
\title{Transformation of coordinate systems}
\description{
  The functions transform a coordinate system into another.
  Currently, essentially only from the earth system to cartesian.
  
  \command{RMtrafo} is the internal basic function that also allows
  to reduce vectors to their norm.
}
 

\usage{
RMtrafo(isotropy)
RFearth2cartesian(coord, units=NULL)
RFearth2dist(coord, units=NULL, ...)
}

\arguments{
  \item{isotropy}{one of the values
    \code{\link{RC_ISOTROPIC}},
    \code{\link{RC_SPACEISOTROPIC}},
    %    ZERORC_SPACEISOTROPIC,
    %    VECTORRC_ISOTROPIC,
    %    SYMMETRIC,
    \code{\link{RC_CARTESIAN_COORD}},
    %   EARTH_COORD,
    %    SPHERICAL_COORD,
    %    CYLINDER_COORD
  }
  \item{coord}{matrix or vector of earth coordinates}
  \item{units}{"km" or "miles"; if not given and 
    \code{RFoptions()$general$units != ""} the latter is used.
    Otherwise \code{"km"}.
  }
  \item{...}{the optional arguments of \command{\link[stats]{dist}}    
  }
}

\details{
  The functions transform between different coordinate systems.
}

\value{
  The function \command{RMtrafo} returns a matrix, in general. For
  fixed column, 
  the results, applied to each row of the matrix, are returned.

  The function \command{RFearth2cartesian} returns a matrix in
  one-to-one correspondance with \code{coord} assuming that the
  earth is an ellipsoid.

  The function \command{RFearth2dist} calculates distances, cf.
  \command{\link[stats]{dist}}, assuming that the
  earth is an ellipsoid.
  
}
\references{
  For calculating the earth coordinates as ellipsoid:
  \itemize{
    \item \url{en.wikipedia.org/wiki/Geographic_coordinate_system}
  
    \item \url{nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html}
  }
}
  
\note{
  Important options are \code{units} and \code{coordinate_system}, see
  \command{\link{RFoptions}}.
}


 \author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 link{constants},
 \command{\link{RMangle}}
}

\examples{
%    library(RandomFields, lib="~/TMP")
data(weather)
(coord <- weather[1:5, 3:4])

(z <- RFfctn(RMtrafo(RC_CARTESIAN_COORD), coord))
(z1 <- RFearth2cartesian(coord)) ## equals t(z)
z1 - t(z) ## 0, i.e., z1 and t(z) are the same
dist(t(z))


(d <- RFearth2dist(coord)) 
d - dist(t(z)) ## 0, i.e., d and dist(t(z)) are the same


\dontshow{\dontrun{ # OK

## RFearth2dist is more precise !!!
library("fields")
(unix.time(z2 <- rdist.earth(coord, miles=FALSE)))
z2

(z2 - as.matrix(d)) / z2 * 100 ## RFearth2dist is more precise


# 
n <- 10000
unix.time(for (i in 1:n) d <- RFearth2dist (coord))) ## the same as dist(t(z))

unix.time(for (i in 1:n)z2 <- rdist.earth(xx)))


data(weather)
xx <- matrix(weather[, 3:4], ncol=2) [ 1:5, ]

xx <- weather[, 3:4][ 1:5, ]
Print(xx)
Print(yy <- RFfctn(RMtrafo(iso=RC_CARTESIAN_COORD), xx), units="miles")

(z1 <- as.matrix(dist(t(yy))))

library(fields)
(z2 <- rdist.earth(xx))

(z1-z2) / z1 * 100
z1
}}


\dontshow{FinalizeExample()}
}


\keyword{spatial}
\keyword{models}



