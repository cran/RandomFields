\name{RMuser}
\alias{RMuser}
\title{User Defined Function}
\description{
  \command{RMuser} allows for a user defined covariance function, variogram
  model, or arbitrary function.

  \bold{RMuser is very slow -- users should avoid this model whenever possible.}
  
}
\usage{
RMuser(type, domain, isotropy, vdim, beta,
       variab.names = c("x", "y", "z", "T"), fctn, fst, snd, envir,  
       var, scale, Aniso, proj)
}
\arguments{
  \item{type}{ See \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}} for
      the range of values of the arguments.

    Default: \code{"shape function"}.
  }
  \item{domain}{ See \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}} for
    the range of values of the arguments.

    Default:  \code{XONLY}.
  }
  \item{isotropy}{See
    \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}} for 
    the range of values of the arguments.
    
    Default:
    \itemize{
      \item  \code{'isotropic'} if \code{type}
      equals \code{'tail correlation function'}, \code{'positive definite'}
      or \code{'negative definite'};
    
      \item
      \code{'cartesian system'} if
      \code{type} indicates a process or simulation method or a shape
      function. 
    }
  }
  \item{vdim}{multivariability.
    
      Default:  \code{vdim} is identified from \code{beta} if given;
      otherwise the default value is \code{1}.
    }
  \item{beta}{a fixed matrix that is multiplied to the return value of the
    given function; the dimension must match.
    Defining a vector valued function and \code{beta} as a vector,
    an arbitrary linear model can be defined. Estimation of \code{beta}
    is, however,  not established yet.
  }
  \item{variab.names}{
    Just the names of the variables. More variable names might
    be given here than used in the function.
    See Detail for the interpretation of variables.
  }
  \item{fctn, fst, snd}{a user defined function and its first,
    second and third derivative, given as
    \code{quote(myfunction(x))} or as 
    \code{quote(myfunction(x, y))}, see Details an Examples below.
  }
  \item{envir}{the environment where the given function shall be evaluated}
  \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
    \command{\link{RMmodel}}. If not passed, the above
    covariance function remains unmodified.}
}

\details{
  Primarily, a function is expected that depends on
  a vector whose components, \eqn{x, y, z, T}, are given separately as
  scalar quantities.
 
  Alternatively, the function might depend only on the first argument
  given by \code{variab.names}.

  A kernel should depend on the first two arguments
  given by \code{variab.names}.
}  

\value{
 \command{\link{RMuser}} returns an object of class
 \code{\link[=RMmodel-class]{RMmodel}}.

 
}
\note{
  \itemize{
    \item
    The use of \command{RMuser} is completely on the risk of the
    user. There is no 
    way to check whether the expressions of the user are
    correct in any sense.

    \item
    Note that \code{x}, \code{y}, \code{z} and \code{T} are reserved
    argument names that define solely the coordinates.
    Hence, none of these names might be used for
    other arguments within these functions.
    \item
    In the user-defined functions the models of
    \pkg{RandomFields} are not recognised, so they cannot be included
    in the function definitions.
    
    \item
    \command{\link{RMuser}} may not be used in connection with obsolete
    commands of RandomFields.
  }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{
  \command{\link{RFfit}},
  \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
  \command{\link{RMmodel}},
  \command{\link{RFsimulate}},
  \link{RC_ISOTROPY}, \link{RC_DOMAIN}
}

\keyword{spatial}
\keyword{models}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


## Alternatively to 'model <- RMexp()' one may define the following
## (which is, however, much slower and cannot use all features of
## RandomFields)
\dontshow{if (RFoptions()$general$storing)
cat("This examples requires 'RFoptions(storing=FALSE)', what is set now.\n");
RFoptions(storing=FALSE)
}
## user-defined exponential covariance model
model <- RMuser(type="positive definite", domain="single variable",
                iso="isotropic", fctn=exp(-x))
x <- y <- seq(1, 10, len=10^(1+interactive()))
plot(model)
z <- RFsimulate(RPcirculant(model), x=x, y=y)
plot(z)

## the kernel, which is the scalar product (not programmed (yet) in
## RandomFields) %to do
model <- RMnugget(var=1e-5) +
         RMuser(type="positive definite", domain="kernel",
                iso="symmetric", fctn=sum(x * y))
x <- y <- seq(1, 10, len=7^(1+interactive()))
z <- RFsimulate(model, x=x, y=y, n=6)
plot(z)
\dontshow{FinalizeExample()}
}
