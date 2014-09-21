\name{Constants}
\alias{constants}
\alias{RC}
\alias{RC_ISOTROPIC}
\alias{RC_ISOTROPY}
\alias{RC_TYPE}
\alias{RC_DOMAIN}
\alias{RC_MONOTONE}
\alias{RC_SPACEISOTROPIC}
\alias{RC_CARTESIAN_COORD}

\title{Constants used in RandomFields (RC constants)}
\description{
  Several constants are provided that might make
  the use of some functions easier, e.g. \command{\link{RFgetModelNames}}
}
\value{
  \code{RC_TYPE = c("tail correlation function", "positive definite",
    "negative definite", "process", 
    "method for Gauss processes", "method for Brown-Resnick processes",
    "shifted shape function",
    "distribution family", "shape function", "trend", "interface",
    "undefined", "other type")}

  \code{RC_DOMAIN = c("single variable", "kernel", "framework dependent", "mismatch")
  }

  \code{RC_ISOTROPY = c("isotropic", "space-isotropic", "zero-space-isotropic",
    "vector-isotropic", "symmetric", "cartesian system",
    "earth system", "spherical system", "cylinder system",
    "non-dimension-reducing", "parameter dependent",
    "<mismatch>")
  }

  \code{RC_MONOTONE = c("mismatch in monotonicity",
                 "submodel dependent monotonicity",
                "previous model dependent monotonicity",
                "parameter dependent monotonicity",
                "not monotone", "monotone", "Gneiting-Schaback class", 
                "normal mixture", "completely monotone", "Bernstein")
  }
 
  \code{RC_ISOTROPIC} gives the numerical code for option
  \code{"isotropic"}

  
  \code{RC_SPACEISOTROPIC} gives the numerical code for option
  \code{"space-isotropic"}

  \code{RC_CARTESIAN_COORD} gives the numerical code for option
  \code{"cartesian system"}



}

\seealso{
  \link{RF},   \link{RM},  \link{RP}, \link{RR}, 
  \command{\link{RFgetModelNames}},
  \link{RMmodelgenerator-class},
  \command{\link{RMtrafo}}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\keyword{spatial}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RC_ISOTROPY[1:5]
RFgetModelNames(isotropy=RC_ISOTROPY[1:5])
\dontshow{FinalizeExample()}
}

