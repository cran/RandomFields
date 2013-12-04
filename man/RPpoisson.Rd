\name{RPpoisson}
\alias{RPpoisson}

\title{Simulation of Random Fields}
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

%\details{ }

 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}


\examples{
set.seed(0)

# example 1
x <- seq(0,25, if (interactive()) 0.02 else 5)
model <- RPpoisson(RMball(), intensity = 2)
z <- RFsimulate(x=x, model=model)
plot(z)
par(mfcol=c(2,1))
plot(z@data[,1:min(length(z@data), 1000)], type="l")
hist(z@data[,1], breaks=0.5 + (-1 : max(z@data)))



# example 2
x <- seq(0,10, if (interactive()) 0.01 else 1)
model <- RMgauss()
z <- RFsimulate(model=RPpoisson(RMtruncsupport(radius=5, model)),
                 x=x, mpp.intensity = 100)
plot(z)

}

\seealso{
  \command{\link{RMmodel}}
  \link{RP},
 \command{\link{RPcoins}}
}




\keyword{spatial}


