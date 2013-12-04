\name{RMtrend}
\alias{RMtrend}
\title{Trend Model}
\description{
 \command{\link{RMtrend}} is a pure trend model with covariance 0.
}
\usage{
RMtrend(mean, plane, polydeg, polycoeff, arbitraryfct, fctcoeff)
}
\arguments{
 \item{mean}{optional; should be a vector of length \eqn{p}{p}, where
 \eqn{p}{p} is the number of variables taken into account by the
 corresponding multivariate random field
 \eqn{(Z_1(\cdot),\ldots,Z_p(\cdot))}{(Z_1(.),\ldots,Z_p(.))}; 
 the \eqn{i}{i}-th component of \code{mean} is interpreted as constant
 mean of \eqn{Z_i(\cdot)}{Z_i(.)}.
 }
 \item{plane}{optional; should be a \eqn{d \times p}{d x p}-matrix where
 \eqn{d} is the dimension of the random field and \eqn{p} the number
 of variables considered; corresponds to a trend described by \eqn{p}{p}
 hyperplanes. The mean of \eqn{Z_i(x)} with \eqn{x= (x_1,\ldots,x_d)}
 is given by
 \deqn{plane_{1i} x_1 + plane_{2i} x_2 + \ldots + plane_{di} x_d.}
 }
 \item{polydeg}{optional; should be an integer vector of length \eqn{p};
 indicates that the mean of \eqn{Z_i(\cdot)}{Z_i(.)} is given by a
 multivariate polynomial of degree less than or equal to the
 \eqn{i}-th component of \code{polydeg}; the coefficients are either assumed
 to be unknown or to be given by \code{polycoeff}.
 }
 \item{polycoeff}{optional; should be a vector of length
 \deqn{{polydeg_1+d\choose d} + \ldots + {polydeg_p+d\choose d}}
 which is the number of monomial basis functions up to degree
 \code{polydeg_1}, \ldots, \code{polydeg_p} in a
 \eqn{d}{d}-dimensional space. Each component of \code{polycoeff}
 gives the coefficient to one basis function; these products are added
 providing \eqn{p}{p} trend polynomials.\cr
 For the order of the monomial basis functions see details.
 CAUTION: This parameter should be used by advanced users only as the
 order might be sophisticated. In many cases it is recommended to use
 \code{arbitraryfct} and \code{fctcoeff} instead.
 } 
 \item{arbitraryfct}{optional; should be a \eqn{p}-variate
 function; the \eqn{i}-th component of \code{arbitraryfct}
 describes the trend surface of \eqn{Z_i(\cdot)}{Z_i(.)};
 the arguments of this function should be location (and
 time) corresponding to the random field to be modelled. If
 \command{\link{RMtrend}} is used in the model definition of the functions
 \command{\link{RFsimulate}}, \command{\link{RFfit}} or
 \command{\link{RFinterpolate}}, the names of the arguments should be
 given by \code{x, y, z, T}.
 }
 \item{fctcoeff}{optional; should be numerical; determines the
 coefficient belonging to \code{arbitraryfct}, i.e. the trend
 is given by \code{fctcoeff * arbitraryfct}. Note that the
 coefficient is the same for each component of
 \code{arbitraryfct};
 \code{fctcoeff} is ignored if \code{arbitraryfct=NULL}.
 } 
}
\details{
 If different trend parameters are given, the corresponding trend
 components are added. Equivalently, \code{+} (see
 \command{\link{RFformula}}) or \command{\link{RMplus}} can be used
 to add trend terms.

 Note that this function refers to trend surfaces in the geostatistical
 framework. Fixed effects in the mixed models framework are
 implemented, as well (see \command{\link{RFformula}}).

 The order of the monomial basis functions and the corresponding
 coefficients given by \code{polycoeff} is the following:\cr
 The first \eqn{{polydeg_1+d\choose d}} components belong to the
 trend polynomial of the first variable, the following
 \eqn{{polydeg_2+d\choose d}} ones to the second one, and so on.
 Within one trend polynomial the monomial basis functions are ordered
 by the powers in an ascending way such that the power of the first
 component varies fastest; e.g. the monomial basis functions up to
 degree \eqn{k}{k} in a two-dimensional space are given by
 \deqn{1, x, \ldots, x^k, y, xy, \ldots, x^{k-1}y, y^2, \ldots, x
 y^{k-1}, y^k.} 
 
}
\value{
 \command{\link{RMtrend}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{Chiles, J. P., Delfiner, P. (1999) \emph{Geostatistics:
 Modelling Spatial Uncertainty.} New York: John Wiley & Sons.
}
\author{
Marco Oesting, \email{oesting@math.uni-mannheim.de}

Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFformula}},
 \command{\link{RFsimulate}},
 \command{\link{RMplus}}
}

\examples{
set.seed(0)
##################################################
# Example 1: # 
# Simulate from model with a plane trend surface #
##################################################

#trend: 1 + x - y, cov: exponential with variance 0.01
model <- ~ RMtrend(mean=1, plane = c(1,-1)) + RMexp(var=0.01)
#equivalent model:
model <- ~ RMtrend(polydeg=1,polycoeff=c(1,1,-1)) + RMexp(var=0.01)
#Simulation
x <- 0:10
simulated <- RFsimulate(model=model, x=x, y=x, grid=TRUE)
plot(simulated)


####################################################################
#
# Example 2: Simulate and fit a multivariate geostatistical model
#
####################################################################
 
# Simulate a bivariate Gaussian random field with trend
# m_1((x,y)) = x + 2*y and m_2((x,y)) = 3*x + 4*y
# where m_1 is a hyperplane describing the trend for the first response
# variable and m_2 is the trend for the second one;
# the covariance function is the sum of a bivariate Whittle-Matern model
# and a multivariate nugget effect
x <- y <- 0:10
model <- RMtrend(plane=matrix(c(1,2,3,4), ncol=2)) + 
         RMparswm(nu=c(1,1)) + RMnugget(var=0.5)
multi.simulated <- RFsimulate(model=model, x=x, y=y, grid=TRUE)
plot(multi.simulated)

\dontshow{\dontrun{
# Fit the Gaussian random field with unknown trend for the second
# response variable and unknown variances
model.na <- RMtrend(plane=matrix(c(1, 2, NA, NA), ncol=2)) + 
            RMparswm(nu=c(1,1), var=NA) + RMnugget(var=NA)
fit <- RFfit(model=model.na, data=multi.simulated)
}}

\dontshow{\dontrun{
##################################################
#
# Example 3:  Simulation and estimation for model with 
#             arbitrary trend functions
#
##################################################

#Simulation
# trend: 2*sin(x) + 0.5*cos(y), cov: spherical with scale 3
model <- ~ RMtrend(arbitraryfct=function(x) sin(x),
 fctcoeff=2) +
 RMtrend(arbitraryfct=function(y) cos(y),
 fctcoeff=0.5) +
 RMspheric(scale=3)
x <- seq(-4*pi, 4*pi, pi/10)
simulated <- RFsimulate(model=model, x=x, y=x, grid=TRUE)
plot(simulated)

################# ?? !!
#Estimation, part 1
# estimate coefficients and scale
model.est <- ~ RMtrend(arbitraryfct=function(x) sin(x), fctcoeff=1) +
 RMtrend(arbitraryfct=function(y) cos(y), fctcoeff=1) +
 RMspheric(scale=NA)
estimated <- RFfit(model=model.est, x=x, y=x, grid=TRUE,
 data=simulated@data, mle.methods="ml")


#Estimation
# estimate coefficients and scale
model.est <- ~ RMtrend(arbitraryfct=function(x) sin(x)) +
 RMtrend(arbitraryfct=function(y) cos(y)) +
 RMspheric(scale=NA)
estimated <- RFfit(model=model.est, x=x, y=x, grid=TRUE,
 data=simulated@data, mle.methods="ml")



##################################################
#
# Example 4: Simulation and estimation for model with 
#            polynomial trend 
#
##################################################

#Simulation
# trend: 2*x^2 - 3 y^2, cov: whittle-matern with nu=1,scale=0.5
model <- ~ RMtrend(arbitraryfct=function(x) 2*x^2 - 3*y^2,
 fctcoeff=1) + RMwhittle(nu=1, scale=0.5)
# equivalent model:
model <- ~ RMtrend(polydeg=2, polycoeff=c(0,0,2,0,0,-3))
x <- 0:20		 
simulated <- RFsimulate(model=model, x=x, y=x, grid=TRUE)
plot(simulated)

#Estimation
# estimate nu and the trend term assuming that it is a polynomial
# of degree 2
model.est <- ~ RMtrend(polydeg=2) + RMwhittle(nu=NA, scale=0.5)
estimated <- RFfit(model=model.est, x=x, y=x, grid=TRUE,
 data=simulated@data, mle.methods="ml")
}}

}
		 
\keyword{spatial}
\keyword{models}



