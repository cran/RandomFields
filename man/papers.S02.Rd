\name{S02}
\alias{S02}
\title{Models for stationary max-stable random fields }
\description{
 Here the code of the paper on \sQuote{Models for stationary max-stable
 random fields} is given.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de};
}
\references{
  \itemize{ 
    \item
    Schlather, M. (2002) Models for stationary max-stable
    random fields. \emph{Extremes} \bold{5}, 33-44.
  } 
}

\examples{
set.seed(0)

\dontshow{\dontrun{
% library(RandomFields, lib="~/TMP"); runif(1)
col=grey((100:0)/100)
pts <- if (interactive()) 512 else 32
x <- (1:pts) / pts * 10
scalegauss <- 1.5
runif(1)
root <- 1/2
RFoptions(mpp.plus=2*scalegauss, CE.force=TRUE)

save.seed <- .Random.seed
ms1 <- MaxStableRF(x, x, model="gauss", param=c(0, 1, 0, scalegauss), 
 maxstable="Bool", grid=TRUE, Print=3)
image(x, x, ms1^root, col=col)

.Random.seed <- save.seed 
scalecirc<- 3.3 # 0.16 0.00 0.18 0.00 0.00 
ms2 <- MaxStableRF(x, x, model="sph", param=c(0, 1, 0, scalecirc), 
 maxstable="Bool", grid=TRUE, Print=3)
image(x, x, ms2^root, col=col)
 
.Random.seed <- save.seed 
ms3 <- MaxStableRF(x, x, model="exponen", param=c(0, 1, 0, 1),
 maxstable="extr", grid=TRUE, Print=3, method="ci")
image(x, x, ms3^root, col=col)

.Random.seed <- save.seed 
ms4 <- MaxStableRF(x, x, model="gauss", param=c(0, 1, 0, 1),
 maxstable="extr", grid=TRUE, Print=3, method="ci")
image(x, x, ms4^root, col=col)
}}

}


\keyword{spatial}
