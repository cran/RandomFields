\name{RRdistr}
\alias{RRdistr}
\alias{RM_DISTR}
\title{Definition of Distribution Families}
\description{
 \command{RRdistr} defines a distribution family given by \code{fct}.
 It is used to introduce \link[=RR]{random parameters} based on
 distributions defined on R.
}
\usage{
RRdistr(name, nrow, ncol, %ddistr, pdistr, qdistr, rdistr,
        envir, ...)
}
\arguments{
  \item{name}{an arbitrary family of distributions. E.g.
    \code{norm()} for the family \command{dnorm}, \command{pnorm},
    \command{qnorm}, \command{rnorm}. See examples below.
  }
  \item{nrow, ncol}{The matrix size (or vector if \code{ncol=1})
    the family returns. Except for very advanced modelling we always have
    \code{nrow=ncol=1}, which is the default.}
%  \item{ddistr, pdistr, qdistr, rdistr}{
%    only for very advanced users. Instead of passing the family name of
%    the distribution, an own family name can be given and all four
%    functions \code{ddistr}, \code{pdistr}, \code{qdistr}, \code{rdistr}
%    are given explicitely.
%    If any of these arguments is passed,
%    all arguments of the function \command{RRdistr} must be passed. 
%  }
  \item{envir}{an environment; defaults to
    \code{\link[base]{new.env}()}.
  }
  \item{...}{Second possibility to pass the distribution family is to
    pass a character string as \code{name} and to give the argument
    within \code{...}. See examples below.
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
  See \link{Bayesian Modelling} for a less technical introduction to
  hierarchical modelling. 

 The use of \command{RRdistr} is completely on the risk of the user. There is no
 way to check whether the expressions of the user are mathematically
 correct.

 Further, \command{\link{RRdistr}} may not be used in connection with obsolete
 commands of RandomFields.
}

\me
\seealso{
  \command{\link{RMmodel}},
  \link{RR},
  \command{\link{RFsimulate}},
  \command{\link{RFdistr}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## here a model with random scale parameter
model <- RMgauss(scale=exp(rate=1))
x <- seq(0,10,0.02)
n <- 10
\dontshow{if(RFoptions()$internal$examples_reduced){warning("reduced 'n'"); n<-1}} 
for (i in 1:n) {  
  readline(paste("Simulation no.", i, ": press return", sep=""))
  plot(RFsimulate(model, x=x, seed=i))
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
  plot(RFsimulate(model, x=x, seed=i))
}

## snd, part 1
## note that the fist 'exp' refers to the exponential function,
## the second to the exponential distribution.
(model1 <- RMgauss(var=exp(3), scale=exp(rate=1)))
x <- 1:100/10
plot(z1 <- RFsimulate(model=model, x=x))

## snd, part 2
## exactly the same result as in the previous example
(model2 <- RMgauss(var=exp(3), scale=RRdistr("exp", rate=1)))
plot(z2 <- RFsimulate(model=model, x=x))
all.equal(model1, model2)

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}

