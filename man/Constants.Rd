\name{Constants}
\alias{constants}
\alias{ZF_ISOTROPY}
\alias{ZF_TYPE}
\alias{ZF_DOMAIN}
\alias{ISOTROPIC}
\alias{SPACEISOTROPIC}
\alias{CARTESIAN_COORD}

\title{Constants used in RandomFields}
\description{
  Several constants are provided that might make
  the use of some functions easier, e.g.\ \command{\link{RFgetModelNames}}
}
\value{
  \code{ZF_TYPE = c("tail correlation function", "positive definite",
    "negative definite", "process", 
    "method for Gauss processes", "method for Brown-Resnick processes",
    "shifted shape function",
    "distribution family", "shape function", "trend", "interface",
    "undefined", "other type")}

  \code{ZF_DOMAIN = c("single variable", "kernel", "calling model", "mismatch")
  }

  \code{ZF_ISOTROPY = c("isotropic", "space-isotropic", "zero-space-isotropic",
    "vector-isotropic", "symmetric", "cartesian system",
    "earth system", "spherical system", "cylinder system",
    "non-dimension-reducing", "parameter dependent",
    "<mismatch>")
  }
  
 
  \code{ISOTROPIC} gives the numerical code for option
  \code{"isotropic"}

  
  \code{SPACEISOTROPIC} gives the numerical code for option
  \code{"space-isotropic"}

  \code{CARTESIAN_COORD} gives the numerical code for option
  \code{"cartesian system"}

}

\seealso{
  \command{\link{RFgetModelNames}},
  \link{RMmodelgenerator-class},
  \command{\link{RMtrafo}}
}
\examples{
RFgetModelNames(isotropy=ZF_ISOTROPY[1:5])
\dontshow{RFoptions(seed=NA)}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}

