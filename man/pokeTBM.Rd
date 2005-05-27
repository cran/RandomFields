\name{pokeTBM}
\alias{pokeTBM}
\title{Transfer of initialisation details}
\description{
  The function allows to transfer the initialisation details of the Out
  register into the In register in case both registers are initialised
  by the same TBM method.

  The function is intended to be replaced by a more general, user
  friendly method.
}
\usage{
pokeTBM(Out, In)
}
\arguments{
  \item{Out}{register from which the initialisation is taken}
  \item{In}{register to which the initialisation parameters are transfered}
}
\value{
  integer value of the error code. If no error has occurred, 0 is
  returned.

  Side effect is that the Out register is deleted.
}
\details{
  This function is intended to specialists use only.
  The function is useful if a random field is to be simulated
  for various sets of points (usually grids)
  such that the value of the random field at identical points of the
  various sets are the same.

  The use of pokeTBM avoids the tedious definition of the ensemble of
  points obtained by a superposition of
  grids as an ensemble of arbitrary points.
  See the example below.
}
\note{
  Both grids should be centred at the origin and have the same grid
  length. These conditions are not checked.
  
  
  Further, it is not checked whether the simulated
  segment on the turning bands of the Out register is large enough
  to cover all the points given in the In register.
  This has to be ensured by the user. If this condition is not
  satisfied
  R may crash.
}
\seealso{
  \command{\link{GaussRF}},  \command{\link{RFparameters}} \command{\link{RFMethods}}
}
\examples{
###### initialisation
runif(1)
rs <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
DeleteAllRegisters()
col <- "red"
model <- "exponential"
param <- c(0, 1, 0, 10)
meth <- "TBM2"

x2 <- seq(-50, 50, 1)
x1 <- seq(-150,150,1)
y1 <- seq(-15, 15, 1)
r1x <- range(x2)
r1y <- range(y1)
r2x <- range(x2)
r2y <- range(y1)

###### simulation of a random field on long thing stripe
z1 <- GaussRF(x1, y1, model=model, param=param, grid=TRUE, register=1,
              method=meth)
get(getOption("device"))(height=1.55, width=12)
par(mar=c(2.2, 2.2, 0.1, 0.1))
image(x1, y1, z1, col=rainbow(100))
polygon(r1x[c(1,2,2,1)], r1y[c(1,1,2,2)], border=col)


###### definition of a random field on a square
InitGaussRF(x2, x2, model=model, param=param, grid=TRUE, register=2,
            method=meth)

###### carry over the TBM details from the stripe
if (err <- pokeTBM(1, 2)) stop(paste(err))

##### simulate with the same random seed and the identical TBM details
assign(".Random.seed", rs, envir=.GlobalEnv)
z2 <- DoSimulateRF(register=2)
get(getOption("device"))(height=4.3, width=4.3)
par(mar=c(2.2, 2.2, 0.1, 0.1))
image(x2, x2, z2, zlim=range(z1), col=rainbow(100))
polygon(r2x[c(1,2,2,1)], r2y[c(1,1,2,2)], border=col)
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\keyword{spatial}
