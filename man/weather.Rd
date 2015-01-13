\name{weather}
\docType{data}
\alias{weather}
\alias{GKS11}
\title{Pressure and temperature forecast errors over the Pacific Northwest}
\usage{data(weather)}
\description{
 Meteorological dataset, which consists of difference between forecasts
 and observations (forcasts minus observations) of temperature and
 pressure at 157 locations in the North American Pacific Northwest. 
}
\format{
  The data frame \code{weather} contains the following columns:
  \describe{
    \item{pressure}{in units of Pascal}
    \item{temperature}{in units of degree Celcius}
    \item{lon}{longitudinal coordinates of the locations}
    \item{lat}{latitude coordinates of the locations}
  }
  
  Furthermore, some results obtained from the data analysis in
  \code{\link{jss14}} are delivered that are
  \code{pars.model}, \code{pars}, \code{whole.model}, \code{whole}.
  
  Finally, the variable \code{information} contains packing information
  (the date and the version of \pkg{RandomFields})
}
\details{
 The forecasts are from the GFS member of the University of Washington
regional numerical weather prediction ensemble (UWME; Grimit and Mass
2002; Eckel and Mass 2005); they were valid on December 18, 2003 at 4 pm
local time, at a forecast horizon of 48 hours.
}
\source{
The data were obtained from Cliff Mass and Jeff Baars in the University
of Washington Department of Atmospheric Sciences.
}
\references{
  \itemize{
    \item
    Eckel, A. F. and Mass, C. F. (2005) Aspects of effective mesoscale,
    short-range ensemble forecasting \emph{Wea. Forecasting} \bold{20},
    328-350.

    \item
    Gneiting, T., Kleiber, W. and Schlather, M. (2010) Matern
    cross-covariance functions for multivariate random fields 
    \emph{J. Amer. Statist. Assoc.} \bold{105}, 1167-1177.
    
    \item
    Grimit, E. P. and Mass, C. F. (2002) Initial results of a
    mesoscale short-range forecasting system over the Pacific
    Northwest \emph{Wea. Forecasting} \bold{17}, 192-205.
  }
}

\seealso{
  A reanalysis has been performed in Section 5 of the
  \link{jss14} paper
}

\examples{
## see 'jss14'
\dontshow{FinalizeExample()}
}

\keyword{datasets}




