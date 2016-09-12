% This file has been created automatically by 'rfGenerateMaths'
\name{Mathematial C functions}
\alias{RFcalc}
\alias{R.models}
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
\alias{R.logb}
\alias{logb}
\alias{R.exp2}
\alias{R.log2}
\alias{log2}
\alias{R.pow}
\alias{R.sqrt}
\alias{sqrt}
\alias{R.hypot}
\alias{R.cbrt}
\alias{R.ceil}
\alias{R.fabs}
\alias{abs}
\alias{R.floor}
\alias{floor}
\alias{R.fmod}
\alias{R.nearbyint}
\alias{R.round}
\alias{round}
\alias{R.trunc}
\alias{trunc}
\alias{R.lrint}
\alias{R.llrint}
\alias{R.lround}
\alias{R.llround}
\alias{R.copysign}
\alias{R.erf}
\alias{R.erfc}
\alias{R.tgamma}
\alias{gamma}
\alias{R.lgamma}
\alias{lgamma}
\alias{R.rint}
\alias{R.nextafter}
\alias{R.nexttoward}
\alias{R.remainder}
\alias{R.fdim}
\alias{R.fmax}
\alias{max}
\alias{R.fmin}
\alias{min}
\title{Transformation of coordinate systems}
\description{
  The functions provide mathematical c functions as \link{RMmodels}
}
 
\usage{
RFcalc(model)
R.minus(a, b, factor)
R.plus(a, b, factor)
R.div(a, b, factor)
R.mult(a, b, factor)
R.const(a)
R.c(a, b, c, d, e, f, g, h, i, j, factor)
R.p(proj, new, factor)
R.is(a, is, b)
R.lon()
R.lat()
%c(...)
R.acos(a)
asin(x)
R.asin(a)
atan(x)
R.atan(a)
atan2(y, x)
R.atan2(a, b)
cos(x)
R.cos(a)
sin(x)
R.sin(a)
tan(x)
R.tan(a)
acosh(x)
R.acosh(a)
asinh(x)
R.asinh(a)
atanh(x)
R.atanh(a)
cosh(x)
R.cosh(a)
sinh(x)
R.sinh(a)
tanh(x)
R.tanh(a)
exp(x)
R.exp(a)
log(x)
R.log(a)
expm1(x)
R.expm1(a)
log1p(x)
R.log1p(a)
logb(x)
R.logb(a)
R.exp2(a)
log2(x)
R.log2(a)
R.pow(a, b)
sqrt(x)
R.sqrt(a)
R.hypot(a, b)
R.cbrt(a)
R.ceil(a)
abs(x)
R.fabs(a)
floor(x)
R.floor(a)
R.fmod(a, b)
R.nearbyint(a)
round(x, ...)
R.round(a)
trunc(x)
R.trunc(a)
R.lrint(a)
R.llrint(a)
R.lround(a)
R.llround(a)
R.copysign(a, b)
R.erf(a)
R.erfc(a)
gamma(x)
R.tgamma(a)
lgamma(x)
R.lgamma(a)
R.rint(a)
R.nextafter(a, b)
R.nexttoward(a, b)
R.remainder(a, b)
R.fdim(a, b)
max(...)
R.fmax(a, b)
min(...)
R.fmin(a, b)
}

\arguments{
  \item{model}{object of class \code{\link[=RMmodel-class]{RMmodel}}, 
	in particular \command{R.model}}
 \item{x,y,a, b, c, d, e, f, g, h, i, j, ...}{constant or object of class \code{\link[=RMmodel-class]{RMmodel}}, 
	in particular \command{R.model}}
  \item{is}{one of \code{"=="}, \code{"!="}, \code{"<="}, \code{"<"}, \code{">="}, \code{">"}}
  \item{factor}{constant factor multiplied with the function. This is
  useful when linear models are built}
  \item{proj}{selection of a component of the vector giving the
  location. Default value is 1.}
  \item{new}{\link{coordinate system} or other 
	\code{\link[=RC_ISONAMES]{kind of isotropy}}
	 which is supposed to be present at this model.
	It shold always be given if the coordinates are not cartesian.
	}
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


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}
\seealso{
 \command{\link{RMmodel}}, \command{\link{RFfctn}},
	 \command{\link{RMtrend}}	
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
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



\dontshow{FinalizeExample()}
}


\keyword{spatial}
\keyword{models}



