\name{RMmodelsMultivariate}
\alias{RMmodelsMultivariate}
\alias{Multivariate RMmodels}
\title{Multivariate models}
\description{
 Here, multivariate and vector-valued covariance models are presented.
}

\details{
\bold{Bivariate covariance models}
\tabular{ll}{
 \command{\link{RMbicauchy}} \tab a bivariate Cauchy model\cr
 \command{\link{RMbiwm}} \tab full bivariate Whittle-Matern model (stationary and isotropic)\cr
 \command{\link{RMbigneiting}} \tab bivariate Gneiting model (stationary and isotropic)\cr
 \command{\link{RMbistable}} \tab a bivariate stable model\cr
 }
 
 \bold{Physically motivated, vector valued covariance and variogram  models}
\tabular{ll}{
\command{\link{RMcurlfree}} \tab curlfree (spatial) vector-valued field (stationary and anisotropic)\cr
\command{\link{RMdivfree}} \tab divergence free (spatial) vector-valued
field (stationary and anisotropic)\cr
\command{\link{RMkolmogorov}} \tab Kolmogorov's model of turbulence\cr
\command{\link{RMvector}} \tab vector-valued field (combining \command{\link{RMcurlfree}} and \command{\link{RMdivfree}})
}

\bold{Multivariate covariance models}
\tabular{ll}{
\command{\link{RMdelay}} \tab  delay effect model
(stationary)\cr
\command{\link{RMderiv}} \tab field and its gradient\cr
\command{\link{RMmatrix}} \tab linear model of coregionalization\cr
\command{\link{RMparswm}} \tab multivariate Whittle-Matern model (stationary and isotropic)\cr
}

\bold{Operators}
\tabular{ll}{
\command{\link{RMcov}} \tab covariance structure given by a multivariate variogram\cr 
\command{\link{RMexponential}} \tab functional returning \eqn{e^C}{exp(C)}\cr
\command{\link{RMmatrix}} \tab linear model of coregionalization\cr
\command{\link{RMmqam}} \tab multivariate quasi-arithmetic mean (stationary)\cr
\command{\link{RMschur}} \tab element-wise product with a positive definite
matrix\cr 
\command{\link{RMtbm}} \tab turning bands operator\cr
}

\bold{Trend models}
\tabular{ll}{
  \command{\link{RMtrend}} \tab for explicit trend modelling\cr
  \command{\link{R.models}} \tab for implicit trend modelling\cr
  \command{\link{R.c}} \tab binding univariate trend models into a vector\cr
}
 

}

\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 % \item Gneiting, T. and Schlather, M. (2004)
 % Statistical modeling with covariance functions.
 % \emph{In preparation.}
 \item Schlather, M. (2011) Construction of covariance functions and
 unconditional simulation of random fields. In Porcu, E., Montero, J.M.
 and Schlather, M., \emph{Space-Time Processes and Challenges Related
 to Environmental Problems.} New York: Springer.
 % \item Schlather, M. (2002) Models for stationary max-stable
 % random fields. \emph{Extremes} \bold{5}, 33-44.
 \item
   Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and
    Strokorb, K. (2015) 
    Analysis, simulation and prediction of multivariate
    random fields with package \pkg{RandomFields}. \emph{
      Journal of Statistical Software}, \bold{63} (8), 1-25,
    url          = \sQuote{http://www.jstatsoft.org/v63/i08/}
  \item
  Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
  Springer, 3rd edition. 
}
}

\seealso{\command{\link{RFformula}}, \command{\link{RMmodels}},
   \command{\link{RM}},
 \command{\link{RMmodelsAdvanced}}

 \sQuote{\href{../doc/multivariate_jss.pdf}{multivariate}}, a vignette
 for multivariate geostatistics
}

\me
\keyword{spatial}

\examples{\dontshow{StartExample()} %   library(RandomFields)
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

n <- 100
x <- runif(n=n, min=1, max=50)
y <- runif(n=n, min=1, max=50)

rho <- matrix(nc=2, c(1, -0.8, -0.8, 1))
model <- RMparswmX(nudiag=c(0.5, 0.5), rho=rho)

## generation of artifical data
dta <- RFsimulate(model = model, x=x, y=y, grid=FALSE)
\dontshow{if (!interactive()) .dataorig <- dta}
## introducing some NAs ...
dta@data$variable1[1:10] <- NA
if (interactive()) dta@data$variable2[90:100] <- NA
\dontshow{if (!interactive()) {print("no NAs introduced"); dta <- .dataorig}}
plot(dta)

% StartExample NICHT weiter vorne!
## co-kriging
x <- y <- seq(0, 50, 1)
\dontshow{if (!interactive()) x <- y <- seq(0, 5, 1)}
k <- RFinterpolate(model, x=x, y=y, data= dta)
plot(k, dta)

## conditional simulation
z <- RFsimulate(model, x=x, y=y, data= dta) ## takes a while
plot(z, dta)

\dontshow{FinalizeExample()}}
