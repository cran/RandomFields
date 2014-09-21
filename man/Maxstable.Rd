\name{Max-stable random fields}
\alias{maxstable}
\alias{Maxstable}
\alias{RPmaxstable}

\title{Simulation of Max-Stable Random Fields}
\description{
  Here, a list of models and methods for simulating max-stable random
  fields is given.

 See also \link{maxstableAdvanced} for more advanced examples.
}


\section{Implemeted models and methods}{
  Models
\tabular{ll}{
 \command{\link{RPbrownresnick}} \tab Brown-Resnick process
 using an automatic choice of the below 3 \code{RPbr*} methods\cr
 \command{\link{RPopitz}} \tab  extremal t
 process \cr 
 \command{\link{RPschlather}} \tab extremal Gaussian
 process \cr 
 \command{\link{RPsmith}} \tab M3 processes \cr
}

Methods
\tabular{ll}{
 \command{\link{RPbrmixed}} \tab simulation of Brown-Resnick processes
 using M3 representation\cr
 \command{\link{RPbrorig}} \tab simulation of Brown-Resnick processes using
 the original definition\cr
 \command{\link{RPbrshifted}} \tab simulation of Brown-Resnick processes using
 a random shift\cr
}
}


\references{
  \itemize{
  \item Kabluchko, Z., Schlather, M. & de Haan, L (2009)
  Stationary max-stable random fields associated to negative definite functions
  \emph{Ann. Probab.} \bold{37}, 2042-2065.

  \item Schlather, M. (2002) Models for stationary max-stable
  random fields. \emph{Extremes} \bold{5}, 33-44.

  \item Smith, R.L. (1990) Max-stable processes and spatial extremes
  Unpublished Manuscript.
  }
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMfbm(alpha=1.5)
x <-  seq(0, 8, if (interactive()) 0.02 else 4) 
z <- RFsimulate(RPbrownresnick(model), x=x, n=if (interactive()) 4 else 1)
plot(z)
\dontrun{\dontshow{
## to do : seq(0, 10, 0.02) oben ist furchtbar langsam. Warum?
}}

\dontrun{\dontshow{
model <- RMball()
x <- seq(0, 10, 5) # nice for   x <- seq(0, 10, 0.02)
n <- if (interactive()) 1000 else 1
z <- RFsimulate(RPsmith(model, xi=0), x, n=n, every=1000)
plot(z)
hist(unlist(z@data), 150, freq=FALSE) #not correct; to do; sqrt(2) wrong
curve(exp(-x) * exp(-exp(-x)), from=-3, to=8, add=TRUE, col=3)
}}

model <- RMgauss()
n <- if (interactive()) 100 else 1
x <- seq(0, 10, if (interactive()) 0.05 else 5)
z <- RFsimulate(RPschlather(model, xi=0), x, n=n)
plot(z)
hist(unlist(z@data), 50, freq=FALSE)
curve(exp(-x) * exp(-exp(-x)), from=-3, to=8, add=TRUE)


## for some more sophisticated models see maxstableAdvanced
\dontshow{FinalizeExample()}

}

\seealso{
 \link{RP},
 \command{\link{RMmodel}},
 \command{\link{RPgauss}},
 \command{\link{RPbernoulli}}
 \command{\link{maxstableAdvanced}}
}




\keyword{spatial}


