\name{S02}
\alias{S02}
\title{Models for stationary max-stable random fields }
\description{
 Here the code of the paper on \sQuote{Models for stationary max-stable
 random fields} is given.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software};
}
\references{
  \itemize{ 
    \item
    Schlather, M. (2002) Models for stationary max-stable
    random fields. \emph{Extremes} \bold{5}, 33-44.
  } 
}

\examples{
RFoptions(seed=0, xi=1)
## seed = 0 : *ANY* simulation will have the random seed 0; set
##            RFoptions(seed=NA) to make them all random again
## xi = 0.5: Frechet margins with alpha=2
\dontshow{StartExample()}
## Due to change in the handling the seeds here are different from the
## seed in the paper.

x <- seq(0, 10, length=128)
\dontshow{if(RFoptions()$internal$examples_red){warning("reduced 'x'");x<-seq(0,10,length=4)}}
# Fig. 1-4
\dontrun{\dontshow{plot(RFsimulate(RPsmith(RMgauss(s=1.5)), x, x))   # < 1 sec
plot(RFsimulate(RPsmith(RMball(s=RRspheric(2, 3,
R=3.3))), x, x)) # 30  sec
plot(RFsimulate(RPschlather(RMexp()), x, x))      #   1 sec
plot(RFsimulate(RPschlather(RMgauss()), x, x))    #  17 sec
}}

\dontshow{FinalizeExample()}
}


\keyword{spatial}
