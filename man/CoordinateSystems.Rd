\name{Coordinate systems}
\alias{zenit}
\alias{Coordinate systems}
\alias{coordinate systems}
\alias{coordinate system}
\alias{coordinates}
\title{Coordinate systems}
\description{
  Implemented Coordinate Systems.
}
\section{Implemented coordinate systems}{
  \itemize{
    \item Cartesian coordinate system
    \item Earth coordinate systems\cr
    The earth is considered as an ellipsoid;
    The first angle takes values in \eqn{[0, 360)},
    the second angle takes values in \eqn{[-90, 90]}.
    \item Spherical coordinate systems\cr
    The earth is considered as an ellipsoid;
    The first angle takes values in \eqn{[0, 2\pi)},
    the second angle takes values in \eqn{[-\pi/2, \pi/2]}.
  }
}

\section{Transformations between the system}{
  \itemize{
    \item Earth to cartesian\cr
    The 3-dimensional resulting coordinates
    are either given in \sQuote{km} or in \sQuote{miles}.
    \item Gnomonic and orthographic projections\cr
    The 2-dimensional resulting coordinates
    are either given in \sQuote{km} or in \sQuote{miles}.
    The projection direction is given by the \code{zenit}.
    \item Earth to spherical\cr
    In this case the Earth is considered as a ball.    
  }
  Cartesian systems cannot be transformed to earth or spherical
  coordinate systems, nor a spherical system to earth coordinates.
}

\section{Options}{
  \describe{
    \item{\code{coord_system}}{character.
      One of the values \code{"auto"}, \code{"cartesian"}, \code{"earth"}
      %   \code{"polar"}, \code{"cylindric"}, \code{"spherical"},
      
      If \code{"auto"}, then the coordinates are considered as
      \code{"cartesian"} except the names of the given coordinates
      indicate a
      different system. Currently, only \code{"longitude"} and
      \code{"latidute"} (or abbreviations of them) are excepted
      as names for given coordinates
      and indicate an earth coordinate systems. See the examples below.
      
      Default: \code{"auto"}
    }

    
    \item{\code{coordidx}}{integer vector of column numbers of the variables in
      the data frame. \code{varidx} can be set alternatively to
      \code{coordnames}.
      This parameter gives the coordinate columns in a data frame
      by starting column and ending column or the sequence.
      An \code{NA} in the second component  means
      \sQuote{until the end}.    
    }
    
    \item{\code{coordnames}}{vector of characters that can be set
      alternatively to coordidx.
      This parameter gives the coordinate columns in a data frame by names.
      If it is \code{"NA"}, then, depending on the context, either
      an error message is returned or it is assumed that the first
      columns give the coordinates.    
    }
    
    \item{\code{coordunits}}{any string.     
      If \code{coordinate_system = "earth"} and longitude and latitude
      are transformed to 3d cartesian coordinates, \code{coordunits}
      determines 
      whether the radius is given in kilometers (\code{"km"}) or miles
      (\code{"miles"}).
      If empty, then \code{"km"} is chosen.
      
      Default: \code{""}
    }
    
  
    \item{\code{new_coord_system}}{
      One of the values \code{"keep"}, \code{"cartesian"}, \code{"earth"},
      \code{"plane"}.
      %   \code{"polar"}, \code{"cylindric"}, \code{"spherical"},
      
      \enumerate{
	\item \code{"keep"}\cr
	The \code{coord_system} is kept (except an explicit transformation
	is given, see \command{\link{RMtrafo}}.

	Note that some classes of models, e.g. completely monotone
	functions and compactly supported covariance models with
	range less than \eqn{\pi} are valid models on a sphere.
	In this case the models are considered as models on the sphere.
	See \link{spherical models} for lists.
      
      \item \code{"cartesian"}\cr
      If \code{coord_system} is \code{"earth"} the coordinates are transformed
      to cartesian coordinates before any model is considered.
      
      \item \code{"orthographic"}, \code{"genomic"} \cr
      If \code{coord_system} is \code{"earth"} the locations are projected
      to a plane before any
      model is considered.
    }
   
    Default: \code{"keep"}
  }

  \item{\code{new_coordunits}}{internal and should not be set by the user.
    
    Default: \code{""}
  }
  
  \item{\code{polar_coord}}{logical.
    If \code{FALSE} the spherical coordinates agree with
    the earth coordinate parametrization, except that
    radians are used for spherical coordinates instead
    of degrees for the earth coordinates.

    If \code{TRUE} the spherical coordinates signify polar coordinates.

    Default : \code{FALSE}
  }

  
  \item{\code{varidx}}{integer vector of length 2.
    \code{varidx} can be set alternatively to
    \code{varnames}.
    This parameter gives the data columns in a data frame, either
    by starting column and ending column.
    An \code{NA} in the second component  means
    \sQuote{until the end}.    
  }
 
  \item{\code{varnames}}{vector of characters that can be set
    alternatively to \code{varidx}.  This parameter gives
  the data columns in a data frame by names.

  if \code{varnames} equals \code{"NA"}
  then for keywords \sQuote{data},
    \sQuote{value} and \sQuote{variable} will be searched for keywords.
    If none of them are found, depending on the context, either
    an error message is returned or it is assumed that the last
    columns give the data.    
  }
     
  \item{\code{varunits}}{vector of characters.
    The default units of the variables.
    
    Default: \code{""}
  }


  \item{\code{xyz_notation}}{logical or \code{NA}. Used by
    \code{\link{RMuser}} only.  

    \code{NA} : automatic choice (if possible)
   
    \code{FALSE} : notation (x, y) should not be understood as kernel definition, not as xyz notation
    
    \code{TRUE}: xyz notation used

    % \code{2}:this value is only used by calls of RFcov and should not
    % be used by a user
  }

 \item{\code{zenit}}{two angles of the central
   projection direction for the gnomonic projection
   (\url{https://en.wikipedia.org/wiki/Gnomonic_projection},
    \url{https://de.wikipedia.org/wiki/Gnomonische_Projektion})
    and the orthographic projection,
  (\url{https://en.wikipedia.org/wiki/Orthographic_projection_in_cartography},
   \url{https://de.wikipedia.org/wiki/Orthografische_Azimutalprojektion}).

    If \code{any(is.na(zenit))}
    then either the value of either of the components may not be \code{NA},
    whose value will be denoted by \eqn{p}.

    If \eqn{p=1} then the mean of the locations is calculated;
    if \eqn{p=Inf} then the mean of the range is calculated.

    Default: \code{c(1, NA)}
  }
  
}

}

\references{
  Covariance models in a cartesian system
  \itemize{
    \item Schlather, M. (2011) Construction of covariance functions and
    unconditional simulation of random fields. In Porcu, E., Montero, J.M.
    and Schlather, M., \emph{Space-Time Processes and Challenges Related
      to Environmental Problems.} New York: Springer.
  }

  Covariance models on a sphere
  \itemize{
    \item Gneiting, T. (2013) Strictly and non-strictly positive definite
    functions on spheres.
    \emph{Bernoulli}, \bold{19}, 1327-1349.
  }

  Tail correlation function
  \itemize{
    \item Strokorb, K., Ballani, F., and  Schlather, M. (2014)
    Tail correlation functions of max-stable processes: Construction
    principles, recovery and diversity of some mixing max-stable
    processes with identical TCF. 
    \emph{Extremes}, \bold{} Submitted.
  }
 }


\seealso{
  \command{\link{RMtrafo}},
  \command{\link{RFearth2cartesian}},
  \command{\link{RPdirect}},
  \command{\link[=spherical models]{models valid on a sphere}},
  \command{\link{RFoptions}}.
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

z <- 1:4
x <- cbind(z, 0)
y <- cbind(0, z)
model <- RMwhittle(nu=0.5)
RFcov(model, x, y, grid=FALSE) ##  standard is (cartesian) model


## same as above, but explicit:
RFcov(model, x, y, grid=FALSE, coord_sys="cartesian") 


## model is not valid on a sphere; x,y coordinates are
## transformed from earth coordinates to spherical coordinates
RFcov(model, x, y, grid=FALSE, coord_sys="earth")


## now the scale is chosen such that the covariance
## values are comparable to those in the cartesian case
RFcov(RMS(model, s= 1 / 180 * pi), x, y, grid=FALSE,
      coord_sys="earth")


## projection onto a plane first. Then the scale is interpreted
## in the usual, i.e. cartesian, sense, i.e. the model does not
## really make sense
RFoptions(zenit = c(2.5, 2.5))
RFcov(model, x, y, grid=FALSE,
      coord_sys="earth", new_coord_sys="orthographic")


## again, here the scale is chosen to be comparable to the cartesian case
## here the (standard) units are [km]
(z1 <- RFcov(RMS(model, s= 6350 / 180 * pi), x, y, grid=FALSE,
             coord_sys="earth", new_coord_sys="orthographic"))


## as above, but in miles
(z2 <- RFcov(RMS(model, s= 6350 / 1.609344 / 180 * pi), x, y, grid=FALSE,
             coord_sys="earth", new_coord_sys="orthographic",
             new_coordunits="miles"))
stopifnot(all.equal(z1, z2))


## again, projection onto a plane first, but now using the
## gnomonic projection
## here the (standard) units are [km]
(z1 <- RFcov(RMS(model, s= 6350 / 180 * pi), x, y, grid=FALSE,
             coord_sys="earth", new_coord_sys="gnomonic"))

## as above, but in miles
(z2 <- RFcov(RMS(model, s= 6350 / 1.609344 / 180 * pi), x, y, grid=FALSE,
             coord_sys="earth", new_coord_sys="gnomonic",
             new_coordunits="miles"))
stopifnot(all.equal(z1, z2, tol=1e-5))




\dontshow{FinalizeExample()}}


\keyword{spatial}

