\name{RMidmodel}
\alias{RMidmodel}
\title{Identical Model}
\description{
  \command{RMidmodel} is the identical operator on models, i.e.
  for objects of class \code{\link[=RMmodel-class]{RMmodel}}.
 }
 \usage{
RMidmodel(phi, vdim, var, scale, Aniso, proj)
}

\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}}
 \item{vdim}{for internal purposes}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}

 
\value{
 \command{\link{RMidmodel}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}


\me

\seealso{
  \command{\link{RMmodel}}.
  \command{RMid},
  \command{RMtrafo},
  \command{RMprod}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMexp()
x <- 0:10
z <- RFsimulate(model, x)

model2 <- RMidmodel(model)
z2 <- RFsimulate(model, x)
sum(abs(as.vector(z)- as.vector(z2))) == 0 # TRUE
\dontshow{FinalizeExample()}}


\keyword{spatial}
\keyword{models}

