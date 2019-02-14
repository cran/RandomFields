% This file has been created automatically by 'rfGenerateMaths'
\name{Mathematial C functions}
\alias{R.models}
\alias{RFcalc}
\alias{math.c}
\alias{R.}
\alias{R.plus}
%\alias{+}
\alias{R.minus}
\alias{-}
\alias{R.mult}
%\alias{*}
\alias{R.div}
\alias{/}
\alias{R.const}
\alias{R.c}
\alias{c}
\alias{R.is}
\alias{R.p}
\alias{R.lon}
\alias{R.lat}
\alias{R.gamma}
\alias{gamma}
\alias{R.acos}
\alias{R.asin}
\alias{asin}
\alias{R.atan}
\alias{atan}
\alias{R.atan2}
\alias{atan2}
\alias{R.cos}
\alias{cos}
\alias{R.sin}
\alias{sin}
\alias{R.tan}
\alias{tan}
\alias{R.acosh}
\alias{acosh}
\alias{R.asinh}
\alias{asinh}
\alias{R.atanh}
\alias{atanh}
\alias{R.cosh}
\alias{cosh}
\alias{R.sinh}
\alias{sinh}
\alias{R.tanh}
\alias{tanh}
\alias{R.exp}
\alias{exp}
\alias{R.log}
\alias{log}
\alias{R.expm1}
\alias{expm1}
\alias{R.log1p}
\alias{log1p}
\alias{R.exp2}
\alias{exp2}
\alias{R.log2}
\alias{log2}
\alias{R.pow}
\alias{^}
\alias{R.sqrt}
\alias{sqrt}
\alias{R.hypot}
\alias{hypot}
\alias{R.cbrt}
\alias{cbrt}
\alias{R.ceil}
\alias{ceiling}
\alias{R.fabs}
\alias{abs}
\alias{R.floor}
\alias{floor}
\alias{R.fmod}
\alias{\%\%}
\alias{R.round}
\alias{round}
\alias{R.trunc}
\alias{trunc}
\alias{R.erf}
\alias{erf}
\alias{R.erfc}
\alias{erfc}
\alias{R.lgamma}
\alias{lgamma}
\alias{R.nextafter}
\alias{R.nexttoward}
\alias{R.remainder}
\alias{R.ilogb}
\alias{R.fdim}
\alias{R.fmax}
\alias{max}
\alias{R.fmin}
\alias{min}
%\alias{R.nearbyint}
%\alias{R.logb}
%\alias{R.lrint}
%\alias{R.llrint}
%\alias{R.lround}
%\alias{R.llround}
%\alias{R.copysign}
%\alias{R.tgamma}
%\alias{R.rint}
\alias{\%\%,ANY,RMmodel-method}
\alias{\%\%,RMmodel,ANY-method}
\alias{*,RMmodel,character-method}
\alias{*,character,RMmodel-method}
\alias{+,RMmodel,character-method}
\alias{+,RMmodel,factor-method}
\alias{+,RMmodel,list-method}
\alias{+,character,RMmodel-method}
\alias{+,data.frame,RMmodel-method}
\alias{+,factor,RMmodel-method}
\alias{-,RMmodel,character-method}
\alias{-,character,RMmodel-method}
\alias{/,RMmodel,character-method}
\alias{/,character,RMmodel-method}
\alias{^,ANY,RMmodel-method}
\alias{^,RMmodel,ANY-method}
\alias{^,RMmodel,character-method}
\alias{^,character,RMmodel-method}
\alias{abs,RMmodel-method}
\alias{acosh,RMmodel-method}
\alias{asin,RMmodel-method}
\alias{asinh,RMmodel-method}
\alias{atan2,ANY,RMmodel-method}
\alias{atan2,RMmodel,ANY-method}
\alias{atan,RMmodel-method}
\alias{atanh,RMmodel-method}
\alias{ceiling,RMmodel-method}
\alias{cos,RMmodel-method}
\alias{cosh,RMmodel-method}
\alias{exp,RMmodel-method}
\alias{expm1,RMmodel-method}
\alias{floor,RMmodel-method}
\alias{lgamma,RMmodel-method}
\alias{log1p,RMmodel-method}
\alias{log2,RMmodel-method}
\alias{log,RMmodel-method}
\alias{round,RMmodel,missing-method}
\alias{sin,RMmodel-method}
\alias{sinh,RMmodel-method}
\alias{sqrt,RMmodel-method}
\alias{tan,RMmodel-method}
\alias{tanh,RMmodel-method}
\alias{trunc,RMmodel-method}
\title{Transformation of coordinate systems}
\description{
  The functions provide mathematical c functions as \link{RMmodels}
}
 
\usage{
RFcalc(model, params, ...)
R.minus(x, y, factor)
R.plus(x, y, factor)
R.div(x, y, factor)
R.mult(x, y, factor)
R.const(x)
R.c(a, b, c, d, e, f, g, h, i, j, l, m, n, o, p, q, ncol, factor)
R.p(proj, new, factor)
R.is(x, is, y)
R.lon()
R.lat()
%c(...)
R.gamma(x)
R.acos(x)
R.asin(x)
R.atan(x)
R.atan2(y, x)
R.cos(x)
R.sin(x)
R.tan(x)
R.acosh(x)
R.asinh(x)
R.atanh(x)
R.cosh(x)
R.sinh(x)
R.tanh(x)
R.exp(x)
R.log(x)
R.expm1(x)
R.log1p(x)
R.exp2(x)
R.log2(x)
R.pow(x, y)
R.sqrt(x)
R.hypot(x, y)
R.cbrt(x)
R.ceil(x)
R.fabs(x)
R.floor(x)
R.fmod(x, y)
R.round(x)
R.trunc(x)
R.erf(x)
R.erfc(x)
R.lgamma(x)
R.remainder(x, y)
R.fdim(x, y)
R.fmax(x, y)
R.fmin(x, y)

\S4method{\%\%}{ANY,RMmodel}(e1,e2)
\S4method{\%\%}{RMmodel,ANY}(e1,e2)
\S4method{*}{RMmodel,character}(e1,e2)
\S4method{*}{character,RMmodel}(e1,e2)
\S4method{+}{RMmodel,character}(e1,e2)
\S4method{+}{RMmodel,factor}(e1,e2)
\S4method{+}{RMmodel,list}(e1,e2)
\S4method{+}{character,RMmodel}(e1,e2)
\S4method{+}{data.frame,RMmodel}(e1,e2)
\S4method{+}{factor,RMmodel}(e1,e2)
\S4method{-}{RMmodel,character}(e1,e2)
\S4method{-}{character,RMmodel}(e1,e2)
\S4method{/}{RMmodel,character}(e1,e2)
\S4method{/}{character,RMmodel}(e1,e2)
\S4method{^}{ANY,RMmodel}(e1,e2)
\S4method{^}{RMmodel,ANY}(e1, e2)
\S4method{^}{RMmodel,character}(e1,e2)
\S4method{^}{character,RMmodel}(e1,e2)
\S4method{abs}{RMmodel}(x)
\S4method{acosh}{RMmodel}(x)
\S4method{asin}{RMmodel}(x)
\S4method{asinh}{RMmodel}(x)
\S4method{atan2}{ANY,RMmodel}(y,x)
\S4method{atan2}{RMmodel,ANY}(y,x)
\S4method{atan}{RMmodel}(x)
\S4method{atanh}{RMmodel}(x)
\S4method{ceiling}{RMmodel}(x)
\S4method{cos}{RMmodel}(x)
\S4method{cosh}{RMmodel}(x)
\S4method{exp}{RMmodel}(x)
\S4method{expm1}{RMmodel}(x)
\S4method{floor}{RMmodel}(x)
\S4method{lgamma}{RMmodel}(x)
\S4method{log1p}{RMmodel}(x)
\S4method{log2}{RMmodel}(x)
\S4method{log}{RMmodel}(x)
\S4method{round}{RMmodel,missing}(x,digits)
\S4method{sin}{RMmodel}(x)
\S4method{sinh}{RMmodel}(x)
\S4method{sqrt}{RMmodel}(x)
\S4method{tan}{RMmodel}(x)
\S4method{tanh}{RMmodel}(x)
\S4method{trunc}{RMmodel}(x)
}


\arguments{
  \item{model,params}{\argModel
    \code{model} is usually a \command{R.model} given here.}
  \item{e1,e2,x,y,a, b, c, d, e, f, g, h, i, j, l, m, n, o, p, q}{
    constant or object of class
    \code{\link[=RMmodel-class]{RMmodel}}, in particular
    \command{R.model}}
  \item{ncol}{in contrast to \link[base]{c}, \command{R.c} also allows for
   defining matrices; \code{ncol} gives the number of columns}
  \item{factor}{constant factor multiplied with the function. This is
    useful when linear models are built}
  \item{is}{one of \code{"=="}, \code{"!="}, \code{"<="}, \code{"<"},
    \code{">="}, \code{">"}} 
   \item{proj}{selection of a component of the vector giving the
    location. Default value is 1.}
  \item{new}{\link{coordinate system} or other 
    \code{\link[=RC_ISO_NAMES]{kind of isotropy}}
    which is supposed to be present at this model.
    It shold always be given if the coordinates are not cartesian.
  }
  \item{digits}{number of digits. Does not work with a \link{RMmodel}}
  \item{...}{\argDots}
%  \item{\code{signature(x = "RMmodel")}}{an \link{RMmodel}, especially a
%  function from this man page. }
%  \item{\code{signature(x = "ANY")}}{numeric or  RMmodel }
% \item{\code{signature(x = "character")}}{no characters allowed. }
}
    
\details{
  \describe{ 
  \item{R.plus}{adds two values}
  \item{R.minus}{substracts two values}
  \item{R.mult}{multiplies two values}
  \item{R.div}{devides two values}
  \item{R.const}{defines a constant }
  \item{R.c}{builds a vector}
  \item{R.is}{evaluates equalities and inequalities;
	note that \code{TRUE} is returned if the equality or
	 inequality holds up to a tolerance given by 
	\code{\link{RFoptions}()$nugget$tol}}
  \item{R.p}{takes a component out of the vector giving the location}
  \item{R.lon, R.lat}{longitudinal and latitudinal coordinate, given
        in the \emph{spherical system}, i.e. in radians.
        (earth system is in degrees).}
  \item{R.round}{Note that \command{R.round} rounds away from zero.
  }
  }
  For the remaining models see the corresponding C functions for
  their return value. 
  (For any \sQuote{R.model} type \sQuote{man model} under Linux.)
}

\value{
Formally, the functions return an object of class
\code{\link[=RMmodel-class]{RMmodel}}, except for
\command{RFcalc} that returns a scalar.
Neither vectors nor parentheses are allowed.
}
%\references{}
  
\note{
  Instead of \code{R.model} the standard function can be used in case
  there is no ambiguity, i.e.,
\code{c(...)},\code{asin(x)}, \code{atan(x)}, \code{atan2(y, x)}, \code{cos(x)},
  \code{sin(x)}, \code{tan(x)}, \code{acosh(x)}, \code{asinh(x)},
  \code{atanh(x)}, \code{cosh(x)}, \code{sinh(x)}, \code{tanh(x)},
  \code{exp(x)}, \code{log(x)}, \code{expm1(x)}, \code{log2(x)},
  \code{log1p(x)}, \code{exp2(x)}, \code{^}, \code{sqrt(x)},
  \code{hypot(a,b)}, \code{cbrt(x)}, \code{ceiling(x)}, \code{abs(x)},
  \code{floor(x)}, \code{round(x)}, \code{trunc(x)},
  \code{erf(x)}, \code{ erfc(x)}, \code{lgamma(x)}. See the examples
  below. 

  
The function \command{RFcalc} is intended for simple calculations only
and it is not excessively tested. Especially, binary
operators should be used with caution.

Note that all the functions here are NOT recognized as being
positive definite (or negative definite), e.g. \command{cos} in
\eqn{R^1}:
\enumerate{
  \item please use the functions given in \command{\link{RMmodels}}
       for definite functions (for \code{cos} see \command{\link{RMbessel}})
  \item Using uncapsulated substraction to build up a covariance
	function is ambiguous, see the example in \command{\link{RMtrend}}
      }      
}


\me
\seealso{
 \command{\link{RMmodel}}, \command{\link{RFfctn}},
	 \command{\link{RMtrend}}	
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


## simple calculation
RFcalc(3 + R.sin(pi/4))

## calculation performed on a field
RFfctn(R.p(1) + R.p(2), 1:3, 1:3) 
RFfctn(10 + R.p(2), 1:3, 1:3) 

## calculate the distances between two vectors
print(RFfctn(R.p(new="iso"), 1:10, 1:10))

## simulation of a non-stationary field where
## anisotropy by a transform the coordinates (x_1^2, x_2^1.5)
x <- seq(0.1, 6, 0.12)
Aniso <- R.c(R.p(1)^2, R.p(2)^1.5)
z <- RFsimulate(RMexp(Aniso=Aniso), x, x)


## calculating norms can be abbreviated:
x <- seq(-5, 5, 5) #0.1)
z2 <- RFsimulate(RMexp() + -40 + exp(0.5 * R.p(new="isotropic")), x, x)
z1 <- RFsimulate(RMexp() + -40 + exp(0.5 * sqrt(R.p(1)^2 + R.p(2)^2)), x, x)
stopifnot(all.equal(z1, z2))
plot(z1)



\dontshow{FinalizeExample()}}


\keyword{spatial}
\keyword{models}



