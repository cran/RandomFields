\name{RRdistr}
\alias{RRdistr}
\title{RRdistr}
\description{
 \command{RRdistr} defines of distribution family given by \code{fct}.
 It is used to introduce \link[=RR]{Random parameters} based on
 distributions defined on R.
}
\usage{
RRdistr(fct, nrow, ncol, envir)
}
\arguments{
  \item{fct}{an arbitrary family of distribution. E.g.
    \code{norm()} for the family \command{dnorm}, \command{pnorm},
    \command{qnorm}, \command{rnorm}.
  }
  \item{nrow, ncol}{The matrix size (or vector if \code{ncol=1})
    the family returns. Except for very advanced modelling we always have
    \code{nrow=ncol=1}, which is the default.}
  \item{envir}{an environment; defaults to \code{\link[base]{new.env}()}.
  }
}

\note{
   \command{\link{RRdistr}} is the generic model introduced
   automatically when distribution families in R are used in the model
   definition. See the examples below.
}

\details{
 \command{\link{RRdistr}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\note{
 The use of \command{RRdistr} is completely on the risk of the user. There is no
 way to check whether the expressions of the user are statistically
 correct.

 Further, \command{\link{RRdistr}} may not be used in connection of obsolete
 commands of RandomFields.

 
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de}
}
\seealso{
  \command{\link{RMmodel}},
  \link{RR}, \command{\link{RRdistr}},
  \command{\link{RFsimulate}},
  \command{\link{RRdistr}}
}

\examples{
set.seed(0)
## here a model with random scale parameter
model <- RMgauss(scale=exp(rate=1))
x <- seq(0,10,0.02)
n <- if (interactive()) 10 else 1
for (i in 1:n) {
  readline(paste("Simulation no.", i, ": press return", sep=""))
  plot(RFsimulate(model, x=x))
}

## another possibility to define exactly the same model above is
## model <- RMgauss(scale=exp())

## note that however, the following two definitions lead
## to covariance models with fixed scale parameter:
## model <- RMgauss(scale=exp(1))   # fixed to 2.7181
## model <- RMgauss(scale=exp(x=1)) # fixed to 2.7181


## here, just two other examples:
## fst
model <- RMmatern(nu=unif(min=0.1, max=2)) # random
for (i in 1:n) {
  readline(paste("Simulation no.", i, ": press return", sep=""))
  plot(RFsimulate(model, x=x))
}

## snd
## note that the fist 'exp' refers to the exponential function,
## the second to the exponential distribution.
model <- RMgauss(var=exp(3), scale=exp(rate=1))
Print(model)
plot(z <- RFsimulate(model=model, x=1:100/10))

%#all.equal(parseModel(model), parseModel(list2RMmodel(parseModel(model))))


}

\keyword{spatial}
\keyword{models}

