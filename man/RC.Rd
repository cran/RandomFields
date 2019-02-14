\name{Constants}
\alias{constants}
\alias{RC}
\alias{RC_ISOTROPIC}
\alias{RC_ISO_NAMES}
\alias{RC_TYPE_NAMES}
\alias{RC_DOMAIN_NAMES}
\alias{RC_MONOTONE_NAMES}
\alias{RC_SPACEISOTROPIC}
\alias{RC_CARTESIAN_COORD}
\alias{RC_GNOMONIC_PROJ}
\alias{RC_ORTHOGRAPHIC_PROJ}
\alias{RC_EARTH_COORDS}
\alias{RC_SPHERICAL_COORDS}
\alias{RC_OPTIMISER_NAMES}
\alias{RC_NLOPTR_NAMES}
\alias{RC_LIKELIHOOD_NAMES}
\alias{RC_EARTH_ISOTROPIC}
\alias{RC_DOUBLEISOTROPIC}
\alias{RC_EARTH_ISOTROPIC}
\alias{RC_UNREDUCED}

\title{Constants used in RandomFields (RC constants)}
\description{
  Several constants are provided that might make
  the use of some functions easier, e.g. \command{\link{RFgetModelNames}}.
}
\value{
  \code{RC_TYPE_NAMES = c(
    "tail correlation", "positive definite", "variogram", "negative definite",
    "point-shape function", "shape function", "trend", "distribution or shape",
    "of manifold type", "process", "method for Gauss process", 
    "normed process (non-negative values with maximum value being 0 or 1)",
    "method for Brown-Resnick process",
    "Smith", "Schlather", "Poisson", "PoissonGauss",
    "distribution family", "interface", "mathematical operator",
    "other type")
  }

  \code{RC_DOMAIN_NAMES = c("single variable", "kernel",
    "framework dependent", "submodel dependent",
    "parameter dependent", "<keep copy>", "mismatch")
  }

  \code{RC_ISO_NAMES = c("isotropic", "space-isotropic",
    "vector-isotropic", "symmetric",
    "cartesian system", "gnomonic projection", "orthographic projection",
    "spherical isotropic", "spherical symmetric", "spherical system", 
    "earth isotropic", "earth symmetric",  "earth system",
    "cylinder system", "non-dimension-reducing", "framework dependent", 
    "submodel dependent", "parameter dependent",
    "<internal keep copy>", "<mismatch>")
  }

  \code{RC_MONOTONE_NAMES = c(
    "not set", "mismatch in monotonicity", "submodel dependent monotonicity",
    "previous model dependent monotonicity",
    "parameter dependent monotonicity",
    "not monotone", "monotone", "Gneiting-Schaback class", 
    "normal mixture", 
    "completely monotone",  
    "Bernstein")
  }
 
  \code{RC_ISOTROPIC} gives the numerical code for option
  \code{"isotropic"}

  
  \code{RC_DOUBLEISOTROPIC} gives the numerical code for option
  \code{"space-isotropic"}

  \code{RC_CARTESIAN_COORD} gives the numerical code for option
  \code{"cartesian system"}

  \code{RC_GNOMONIC_PROJ} gives the numerical code for
  the gnomonic projection, see also \code{zenit} in \code{\link{RFoptions}}.

  \code{RC_ORTHOGRAPHIC_PROJ} gives the numerical code for
  the orthographic
  projection, see also \code{zenit} in \code{\link{RFoptions}}.

  \code{RC_EARTH_COORDS} gives the numerical code for option
  \code{"earth coordinates"}
  
  \code{RC_EARTH_ISOTROPIC} gives the numerical code for option
  \code{"earth isotropic"}
  

  \code{RC_SPHERICAL_COORDS} gives the numerical code for option
  \code{"earth coordinates"}

  \code{RC_OPTIMISER_NAMES} and \code{RC_NLOPTR_NAMES} give the
  names for the options \code{optimiser} and \code{algorithm},
  respectively, \code{RFfitoptimiser}.

  \code{RC_LIKELIHOOD_NAMES = c("auto", "full", "composite",
    "tesselation")}
  gives the names of the ML variants: (i) internal
  choice according to the
  number of data, (ii) full likelihood, (iii) (pairwise) composite likelihood,
  and (iv) composite likelihood based on a tessellation of the space.
    
}


\seealso{
  \link{RF},   \link{RM},  \link{RP}, \link{RR}, \link{R.},
  \command{\link{RFgetModelNames}},
  \link{RMmodelgenerator-class},
  \command{\link{RMtrafo}}.
}
\me

\keyword{spatial}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RC_ISO_NAMES
RC_ISO_NAMES[RC_ISOTROPIC:RC_CARTESIAN_COORD + 1]
\dontrun{
RFgetModelNames(isotropy=RC_ISO_NAMES[RC_ISOTROPIC:RC_CARTESIAN_COORD +
1])
}
\dontshow{FinalizeExample()}}

