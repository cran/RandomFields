\name{RMuser}
\alias{RMuser}
\title{User Defined Function}
\description{
  \command{RMuser} allows for a user defined covariance function, variogram
 model, or arbitrary function
}
\usage{
RMuser(type, domain, isotropy, vdim, beta, fctn, fst, snd, envir,
 var, scale, Aniso, proj)
}
\arguments{
  \item{type, domain, isotropy, vdim}{see \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}.}
  \item{beta}{a fixed matrix that is multiplied to the return value of the
    given function; the dimension must match.
    Defining a vector valued function and \code{beta} as a vector,
    an arbitrary linear model can be defined. Estimation of \code{beta}
    is, however,  not established yet.
  }
  \item{envir}{the environment where the given function shall be evaluated}
  \item{fctn, fst, snd}{a user defined function and its first,
    second and third derivative, given as
    \code{quote(myfunction(x))} or as 
    \code{quote(myfunction(x, y))}, see Details an Examples below.
  }
  \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
    \command{\link{RMmodel}}. If not passed, the above
    covariance function remains unmodified.}
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
    parameter names that define solely the coordinates.
    Hence, none of these names might be used for passed functions
    or as other parameters within these functions.
    If they were used, \pkg{RandomFields} would overwrite them and hence
    cause something unpredictable.

    \item
    In the user-defined functions the models of
    \pkg{RandomFields} are not recognised, so they cannot be included.
    
    \item
    \command{\link{RMuser}} may not be used in connection with obsolete
    commands of RandomFields.
  }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de}
}
\seealso{
  \command{\link{RMmodel}},
  \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
  \command{\link{RFsimulate}},
  \command{\link{RFfit}}.
}

\examples{
RFoptions(seed=0)
## Alternatively to 'model <- RMexp()' one may define the following
## (which is, however, much slower and cannot use all features of
## RandomFields)

\dontshow{RFoptions(storing = FALSE)}

model <- RMuser(type="positive definit", domain="single variable",
                iso="isotropic", fctn=exp(-x))
x <- y <- seq(1, 10, len=10^(1+interactive()))
plot(model)
z <- RFsimulate(RPcirculant(model), x=x, y=y, grid=TRUE)
plot(z)

% plot(model)
model <- RMnugget(var=1e-5) +
         RMuser(type="positive definit", domain="kernel",
                iso="symmetric", fctn=sum(x * y))
x <- y <- seq(1, 10, len=7^(1+interactive()))
z <- RFsimulate(model, x=x, y=y, grid=TRUE, n=6)
plot(z)
\dontshow{RFoptions(seed=NA)}
}

\keyword{spatial}
\keyword{models}



%x <- y <- seq(1, 10, len=10)
%model <- RMexp()
%
%z <- RFsimulate(model = model, x=x, y=y, grid=TRUE, seed=0)
%plot(z)
%
%model <- RMuser(type="positive definit", domain="single variable",
%                iso="isotropic",  fctn=quote(exp(-x)))
%
%z1 <- RFsimulate(model = RPcirculant(model), x=x, y=y, grid=TRUE, seed=0)
%plot(z1)
%stopifnot(sum(abs(z@data-z1@data)) == 0) ## 0

