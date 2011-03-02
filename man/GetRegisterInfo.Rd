\name{GetRegisterInfo}
\alias{GetRegisterInfo}
\alias{GetModelInfo}
\alias{GetModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal information}
\description{
  The function returns internal information about the simulation
  of a random field and the calculation of covariance functions
}
\usage{
GetRegisterInfo(register, ignore.active = FALSE, max.elements=10^6)
GetModelInfo(register, modelreg, level=3, gatter=FALSE)
GetModel(register, modelreg, modus=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{register}{-2, -1, 0:9; (-1 and -2 only work with \command{GetModelInfo})
    place where intermediate calculations are stored;
    the numbers \code{0:9} are aliases for 10 internal registers, see also
    \command{\link{GaussRF}}}
  \item{modelreg}{value in  0..4 with the following meaning:
    \itemize{
      \item 0 : register for the functions \command{Covariance},
      \command{CovarianceFct}, ect
      \item 1 : register for simulation
      \item 2 : register for covariance calculation within Kriging
      \item 3 : register for MLE (likelihood)
      \item 4 : register for MLE (bounds)
    }
    Positive values refer to the registers, see \command{\link{GaussRF}}.
  }
  \item{ignore.active}{logical. If \code{FALSE} and the register has
    non-active flag (because of an error or it is deleted) then a list
    is returned that contains only the element \code{active=FALSE}.
    Otherwise, a more complete list is returned -- only for internal
    debugging purposes.
    \bold{Warning:} \code{ignore.active=TRUE} may cause failure of R.    
  }
  \item{max.elements}{integer; since \command{GetRegisterInfo}
    might be a dump of simulation that needs a huge amount of
    memory and since all entries are copied, the maximal amount
    of available memory might be easily exceeded.  Therefore,
    only the size of the critical parts are returned and not the
    vector or matrix itself, if the number of elements exceeds
    \code{max.elements}. 
  }
  \item{level}{integer; level of details}
  \item{modus}{integer; see details}
  \item{gatter}{boolean; only for very advanced interests.}
}
\details{
  \code{GetRegisterInfo(register, ignore.active=TRUE)} is
  useful for debugging and specialists' need to control
  the algorithm, see the examples in \code{\link{RFparameters}}
  and \code{\link{GaussRF}}.
  
  If \command{\link{RFparameters}}\code{()$Storage=FALSE} then values of
  the internal registers are not kept if \command{\link{GaussRF}} or
  \command{\link{DoSimulateRF}} has been called.
  Hence \code{GetRegisterInfo} cannot
  provide any information.

  \code{GetModelInfo} returns the complete information on the internal
  model structure. It allows for \code{register=-1} returning the
  model structure for the last use of \command{\link{CovarianceFct}},
 \link{sophisticated} models
  or \command{\link{Variogram}} or similar commands.
  \code{register=-2} is for internal use only.

  \code{GetModel} returns only the parts of the internal model structure
  that have been defined by the user.
  The \code{modus} is ignored if \code{PracticalRange=FALSE}
  in \command{\link{RFparameters}}. In case \code{PracticalRange=TRUE}
  the \code{modus} has the following effects
  \enumerate{
    \item{0}{the model is returned without the explicite scale transformation}
    \item{1}{the model is returned the way it is stored, including the
      scale transformation}
    \item{2}{the scale
      transformation of the the practical range is included,
      but a simplified model will be returned
    }
  }
}
\note{
  Put \code{Storing=TRUE}, see \command{\link{RFparameters}}
  if you like to have more
  internal information in case of an expected failure of an initialisation of
  a random field simulation.
}
\value{
  List of internal information is returned.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}

\seealso{\command{\link{GaussRF}}}
\examples{
GaussRF(1:4, grid=TRUE, model="exp", param=c(1,2,3,4), Storing=TRUE)
Print(GetRegisterInfo(0))

# finally, a register that could not be successfully initialised
try(GaussRF(runif(4), grid=FALSE, model="exp", param=c(1,2,3,4),
            me="ci", Storing=TRUE))
if (interactive()) Print(GetRegisterInfo(0, TRUE))
}
\keyword{spatial}
