\name{RMangle}
\alias{RMangle}
\title{Anisotropy matrix given by angle}
\description{
  \command{RMangle} delivers an anisotropy matrix for the
  argument \code{Aniso} in \command{\link{RMmodel}} in two dimensions.
  \command{RMangle} requires one or two stretching values, passed by
  \code{ratio} or \code{diag}, and an \command{angle}.

  In two dimensions and with
  \code{angle} equal to \eqn{a} and \code{diag} equal to \eqn{(d1, d2)} the
  anisotropy matrix \eqn{A} is
  
  \code{A = diag(d1, d2) \%*\% matrix(ncol=2, c(cos(a), sin(a), -sin(a), cos(a)))}

  In three dimensions and with
  \code{angle} equal to \eqn{a}, second angle \eqn{L}
  and \code{diag} equal to \eqn{(d1, d2, d3)} the
  anisotropy matrix \eqn{A} is
 
  \code{A = diag(d1, d2, d3) \%*\% matrix(ncol=3, c(cos(a) * cos(L), sin(a) * cos(L), sin(L), -sin(a), cos(a), 0, -cos(a) * sin(L), -sin(a) * sin(L), cos(L)
    ))}
  i.e. \eqn{Ax} turns a vector x first in the \eqn{x-z} plane, then
  in the \eqn{x-y} plane.
}
 

\usage{
RMangle(angle, lat.angle, ratio, diag) 
}

\arguments{
  \item{angle}{angle \code{a}}
  \item{lat.angle}{second angle; in 3 dimensions only}
 \item{ratio}{equivalent to \code{diag=c(1, 1/ratio)}; in 2 dimensions only}
 \item{diag}{the diagonal components of the matrix}
}

\value{
 \command{\link{RMangle}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}

\me
\seealso{
 \command{\link{RMtrafo}}, \command{\link{RMmodel}}
  }

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMexp(Aniso=RMangle(angle=pi/4, ratio=3))
plot(model, dim=2)

x <- seq(0, 2, 0.05)
z <- RFsimulate(x, x, model=model)
plot(z)


model <- RMexp(Aniso=RMangle(angle=pi/4, lat.angle=pi/8, diag=c(1,2,3)))
x <- seq(0, 2, 0.2)
z <- RFsimulate(x, x, x, model=model)
plot(z, MARGIN.slices=3)


## next model gives an example how to estimate the parameters back
n <- 20
x <- runif(n, 0, 10)
y <- runif(n, 0, 10)
coords <- expand.grid(x, y)
model <- RMexp(Aniso=RMangle(angle=pi/4, diag=c(1/4, 1/12)))
d <- RFsimulate(model, x=coords[, 1], y=coords[, 2], n=10)
estmodel <- RMexp(Aniso=RMangle(angle=NA, diag=c(NA, NA)))
system.time(RFfit(estmodel, data=d, modus_operandi='sloppy'))

\dontshow{FinalizeExample()}}


\keyword{spatial}
\keyword{models}




