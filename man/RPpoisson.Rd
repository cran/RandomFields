\name{RPpoisson}
\alias{RPpoisson}

\title{Simulation of Poisson Random Fields}
\description{
Shot noise model, which is also called moving
 average model, trigger process, dilution random field, and by several
 other names.
}

\usage{
RPpoisson(phi, intensity)
}

\arguments{
 \item{phi}{the model, \command{\link{RMmodel}}, gives the shape
 function to be used}
 \item{intensity}{the intensity of the underlying stationary Poisson
 point process
 }
 }
 
\me

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

# example 1: Posson field based on disks with radius 1
x <- seq(0,25, 0.02)
model <- RMball()
z <- RFsimulate(RPpoisson(model), x, intensity = 2)
plot(z)
par(mfcol=c(2,1))
plot(z@data[,1:min(length(z@data), 1000)], type="l")
hist(z@data[,1], breaks=0.5 + (-1 : max(z@data)))


# example 2: Poisson field based on the normal density function
# note that
# (i) the normal density as unbounded support that has to be truncated
# (ii) the intensity is high so that the CLT holds
x <- seq(0, 10, 0.01)
model <- RMtruncsupport(radius=5, RMgauss())
z <- RFsimulate(RPpoisson(model), x, intensity = 100)
plot(z)

\dontshow{FinalizeExample()}}

\seealso{
  \command{\link{RMmodel}},
  \link{RP},
 \command{\link{RPcoins}}.
}




\keyword{spatial}


