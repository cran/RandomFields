\name{DeleteRegister}
\alias{DeleteRegister}
\alias{DeleteAllRegisters}
\title{Deleting Intermediate Results}
\description{
  \code{DeleteRegister} deletes an internal register of the random
        field simulation

  \code{DeleteAllRegisters} deletes all internal registers
}
\usage{DeleteRegister(nr=0)
DeleteAllRegisters()
}
\arguments{
  \item{nr}{number of the register to be deleted; in 0:9}
}
\details{
  \code{DeleteRegister} and \code{DeleteAllRegister} should be used if
  \enumerate{
    \item the simulations are finished and there
    is not enough memory available to continue the R session,
    and \command{\link{RFparameters}}\code{()$Storing} has been
    \code{TRUE}.
    
    \item one of the parameters
    \command{\link{RFparameters}}\code{()$Storing} or
    \command{\link{RFparameters}}\code{()$PracticalRange}
    ((or \command{\link{RFparameters}}\code{()$TBM.method},
    once this parameter is available)) has been changed by
    the user, and if these options should become
    global.
  }
  See \link{GaussRF} and \link{RFparameters} for further details.
}
\value{
  \code{NULL}.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de},
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}

\seealso{\command{\link{GaussRF}}, \code{\link{RandomFields}} and \command{\link{RFparameters}}.}

\examples{
 DeleteAllRegisters() ## no effect visible for the user
}
\keyword{spatial}





