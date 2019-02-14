\name{RMpolygon}
\alias{RMpolygon}
\title{RMpolygon}
\description{
  \command{RMpolygon} refers to the indicator function of
  a typical Poisson polygon, used for instance in the
  (mixed) Storm process.
}
\usage{
RMpolygon(lambda)
}
\arguments{
  \item{lambda}{
    intensity of the hyperplane process creating the random shape function.

    The default value is \code{1}.
  }
}



\seealso{
  \command{\link{RMball}},
  \command{\link{RMspheric}},
  \command{\link{RFsimulate}},
  \command{\link{RMmodel}}.
}

\author{
  Felix Ballani,
  \url{http://www.mathe.tu-freiberg.de/sto/mitarbeiter/felix-ballani}
  
  \martin
}

\references{
  Poisson polygons / Poisson hyperplane tessellation
  \itemize{
    \item Lantuejoul, C. (2002)
    \emph{Geostatistical Simulation: Models and Algorithms.}
    Springer.
  }

  Poisson storm process
  \itemize{
    \item Lantuejoul, C., Bacro, J.N., Bel L. (2011)
    Storm processes and stochastic geometry. \emph{Extremes},
    \bold{14}(4), 413-428.
  }

  Mixed Poisson storm process
  \itemize{
    \item  Strokorb, K., Ballani, F., and  Schlather, M. (2014)
    Tail correlation functions of max-stable processes: Construction
    principles, recovery and diversity of some mixing max-stable
    processes with identical TCF. 
    \emph{Extremes}, \bold{} Submitted.
  }
}

\keyword{spatial}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

\dontshow{\dontrun{
model <- RPcoins(shape=RMpolygon())
x <- seq(0, 10, 0.02)
plot(RFsimulate(model, x=x))
}}
\dontshow{FinalizeExample()}}