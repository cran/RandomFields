\name{RFformula}
\alias{RFformula}

\title{RFformula - syntax to design random field models with trend
 or linear mixed models}

\description{It is described how to create a formula, which can
 e.g. be used as an argument of \command{\link{RFsimulate}} and
 \command{\link{RFfit}} to simulate and to fit data accordingly to the
 model described by the formula.

 In general, the created formula serves two purposes:
 \itemize{
 \item to describe models in the \dQuote{Linear Mixed
 Models}-framework including fixed and random effects
 \item to define models for random fields including trend surfaces from
 a geostatistical point of view.
 }

 Thereby, fixed effects and trend surfaces are adressed via
 the expression \command{\link{RMfixed}} and the function
 \command{\link{RMtrend}}; the covariance structures of the zero-mean
 multivariate normally distributed random effects and random field
 components are adressed by objects of class \code{\link[=RMmodel-class]{RMmodel}}, which
 allow for a very flexible covariance specification.
}

\details{ 
 The formula should be of the type
 \deqn{response ~ fixed effects + random effects + error term}
 or 
 \deqn{response ~ trend + zero-mean random field + nugget effect,}
 respectively.
 
 Thereby:
 
 \itemize{
 \item response\cr
 optional; name of response variable

 \item fixed effects/trend:\cr
 optional, should be a sum (using \command{\link[=RMplus]{+}})%%check link
 of components either of the form \code{X@RMfixed(beta)} or
 \code{\link{RMtrend}(...)} with \eqn{X} being a design matrix 
 and \eqn{\beta} being a vector of coefficients (see
 \command{\link{RMfixed}} and \command{\link{RMtrend}}).\cr
 Note that a fixed effect of the form \eqn{X} is interpreted as
 \code{X@RMfixed(beta=NA)} by default (and \eqn{\beta} is estimated
 provided that the formula is used in \command{\link{RFfit}}).
 
 \item random effects/zero-mean random field:\cr
 optional, should be a sum (using \command{\link[=RMplus]{+}})%%check link
 of components of the form \code{Z@model}
 where \eqn{Z} is a design matrix and \code{model} is an object of
 class \code{\link[=RMmodel-class]{RMmodel}}.\cr
 \code{Z@model} describes a vector of random effects which is
 normally distributed with zero mean and covariance matrix \eqn{Z
 \Sigma Z^T} where \eqn{Z^T} is the transpose of \eqn{Z} and
 \eqn{\Sigma} is the covariance matrix according to \code{model}.\cr
 Note that a random effect/random fluctuation of the form
 \code{model} is viewed as \code{I@model} where \eqn{I} is the
 identity matrix of corresponding dimension.
 
 \item error term/nugget effect\cr
 optional, should be of the form \code{\link{RMnugget}(...)}.
 \command{\link{RMnugget}} describes a vector
 of iid Gaussian random variables.

 Please note that the character \dQuote{@} in the RFformula-context can only
 be used to multiply design-matrices with corresponding vectors of
 fixed or random effects, whereas in the context of S4-classes \dQuote{@} is
 used to access slots of corresponding objects.
 }
 
}

\note{
  (additional) argument names should always start with a capital
  letter. Small initial letters are reserved  for \command{\link{RFoptions}}.
}

\references{
 \itemize{
 \item Chiles, J.-P. and P. Delfiner (1999) \emph{Geostatistics. Modeling
 Spatial Uncertainty.} New York, Chichester: John Wiley & Sons.
 
 \item McCulloch, C. E., Searle, S. R. and Neuhaus, J. M. (2008)
 \emph{Generalized, linear, and mixed models.} Hoboken, NJ: John
 Wiley & Sons.
 
 \item Ruppert, D. and Wand, M. P. and Carroll, R. J. (2003)
 \emph{Semiparametric regression.} Cambridge: Cambridge University
 Press.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \code{\link[=RandomFields-package]{RandomFields}}.}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFoptions(modus_operandi="sloppy")

\dontshow{\dontrun{
##############################################################
#
# Example 1: Simulation and fitting of a two-dimensional
# Gaussian random field with exponential covariance function
#
###############################################################

V <- 10
S <- 0.3
model <- RMexp(var=V, scale=S)
# Simulate from model
x <- y <- seq(1,3, by=0.1) # grid coordinates
simulated <- RFsimulate(model = model, x=x, y=y, grid=TRUE)
#returns an object of class RFspatialPointsDataFrame including
#coordinates and the simulated response vector
plot(simulated)

# add overall mean/trend of 3 to the response/simulated data 
simulated@data <- simulated@data + 3

# alternative:
simulated <- RFsimulate(model = ~ 1@RMfixed(beta=3) +
 RMexp(var=V, scale=S),
 x=x, y=y, grid=TRUE, V=V, S=S)

# Estimate parameters of underlying covariance function via
# maximum likelihood
model.na <- ~ 1@RMfixed(beta=NA) + RMexp(var=NA, scale=NA)
fitted <- RFfit(model=model.na, data=simulated)

# compare sample mean of data with ML estimate:
mean(simulated@data[,1]) 
fitted@ml@trend
}}

\dontshow{\dontrun{
##############################################################
#
# Example 2: Fitting a standard linear mixed model using maximum
# likelihood to estimate the variance components
#
###############################################################

# Y = W*beta + Z*u + e
# where u ~ N(0, (sigma_u)^2 A) and e ~ N(0, (sigma_e)^2)
W <- rep(1,times=10)
Z <- matrix(rnorm(150) ,nrow=10, ncol=15)
A <- RFcovmatrix(0:14, model=RMexp()) 
response <- W*5 + Z%*%matrix(1:225, nrow=15)%*%rnorm(15, sd=10) + rnorm(10, sd=3)

# Estimate beta, (sigma_u)^2 and (sigma_e)^2:
model <- response ~ W@RMfixed(beta=NA) +
 Z@RMconstant(A, var=NA) +
 RMnugget(var=NA)
fitted <- RFfit(model=model, data=response, W=W, Z=Z, A=A)
}}


\dontshow{\dontrun{
#### THIS EXAMPLE IS NOT PROGRAMMED YET


###########################################################
#
# Example 3: Simulate and fit a geostatistical model
#
###########################################################

# Simulate a Gaussian random field with trend m((x,y)) = 2 + 1.5 x - 3 y
# with vector of coordinates (x,y)
# and covariance function C(s,t) = 3*exp(-||(2*(s_1-t_1),s_2-t_2)||) +
# 1.5*exp(-||(2*(s_1-t_1),s_2-t_2)||^2)
# for s=(s_1,s_2) and t=(t_1,t_2)

model <- ~ RMtrend(mean=2) +
 RMtrend(arbitraryfct=function(x) x, fctcoeff=1.5) +
 RMtrend(arbitraryfct=function(y) y, fctcoeff=-3) + 
 RMplus(RMexp(var=3), RMgauss(var=1.5),
 Aniso=matrix(c(2,0,0,1),nrow=2))

simulated <- RFsimulate(model=model,x=seq(0,2,0.1),y=seq(-1,3,0.2),
 grid=TRUE)

# equivalent model
model <- RMtrend(polydeg=1, polycoeff=c(2, 1.5, .3)) +
 RMplus(RMexp(var=3), RMgauss(var=1.5),
 Aniso=matrix(c(2,0,0,1), nrow=2))

simulated <- RFsimulate(model=model, x=seq(0,2,0.1), y=seq(-1,3,0.2),
 grid=TRUE)

# Estimate trend (polynomial of degree 1) and variance components such
# that var_exp = 2*var_gauss as in the true model used for simulation
model.na <- ~ RMtrend(polydeg=1) +
 RMplus(RMexp(var=2),RMgauss(var=2),var=NA,
        Aniso=matrix(c(NA,0,0,NA),nrow=2))
fit <- RFfit(model=model.na, data=simulated)
}}



\dontshow{\dontrun{
#### THIS EXAMPLE IS NOT PROGRAMMED YET

####################################################################
#
# Example 4: Simulate and fit a multivariate geostatistical model
#
####################################################################

# Simulate a bivariate Gaussian random field with trend
# m_1((x,y)) = x + 2*y and m_2((x,y)) = 3*x + 4*y
# where m_1 is a hyperplane describing the trend for the first response
# variable and m_2 is the trend for the second one;
# the covariance function is a multivariate nugget effect
x <- y <- 0:10
model <- ~ RMtrend(plane=matrix(c(1,2,3,4), ncol=2)) +
           RMnugget(var=0.5, vdim=2)
multi.simulated <- RFsimulate(model=model, x=x, y=y, grid=TRUE)

# Fit the Gaussian random field with unknown trend for the second
# response variable and unknown variances
model.na <- ~ RMtrend(plane=matrix(c(NA,NA,NA,NA), ncol=2)) + 
              RMnugget(var=NA, vdim=2)
fit <- RFfit(model=model.na, data=multi.simulated)
}}



\dontshow{RFoptions(modus_operandi="normal")}

\dontshow{FinalizeExample()}
}

\keyword{spatial}
