\name{RMtruncsupport}
\alias{RMtruncsupport}
\title{Truncating the Support of a Shape Function}
\description{
  \command{\link{RMtruncsupport}} may be used to truncate the support of a shape function when
 Poisson fields or M3 processes are created.
 }
\usage{
RMtruncsupport(phi, radius)
}

\arguments{
 \item{phi}{function of class \code{\link[=RMmodel-class]{RMmodel}}}
 \item{radius}{truncation at \code{radius}}
}

 
\value{
 \command{\link{RMtruncsupport}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}


\me
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RMmatrix}},
 \command{\link{RPpoisson}}.
 }

\keyword{spatial}
\keyword{models}

\references{
  \itemize{
    \item Schlather, M. (2002) Models for stationary max-stable
    random fields. \emph{Extremes} \bold{5}, 33-44.
    }
}
 
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again 

model <- RMgauss()
model1 <- RMtruncsupport(model, radius=1)
plot(model)
lines(model1, col="red")

## For a real application of 'RMtruncsupport' see example 2 of 'RPpoisson'.

\dontshow{\dontrun{
model <- RPcoins(shape = RMtruncsupport(RMgauss(), radius=1.5), intensity=1)
x <- seq(0, 10, 0.02)
plot(RFsimulate(model, x=x))
}}
\dontshow{FinalizeExample()}}