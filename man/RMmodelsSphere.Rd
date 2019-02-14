\name{Spherical models}
\alias{RMmodelsSphere}
\alias{sphere}
\alias{spherical models}
\alias{Spherical models}
\alias{earth models}
\alias{Earth models}
\title{Covariance models valid on a sphere}
\description{
  This page summarizes the covariance models that can be used
  for spherical coordinates (and earth coordinates).
}

\details{
  The following models are available:
  
  \bold{Completely monotone functions allowing for arbitrary scale}
  \tabular{ll}{
     \command{\link{RMbcw}} \tab Model bridging stationary and
    intrinsically stationary processes for \eqn{\alpha \le 1}
    and \eqn{\beta < 0}\cr
       \command{\link{RMcubic}} \tab cubic model\cr
   \command{\link{RMdagum}} \tab Dagum model with \eqn{\beta < \gamma}
    and \eqn{\gamma \le 1}\cr
    \command{\link{RMexp}} \tab exponential model \cr
    \command{\link{RMgencauchy}} \tab generalized Cauchy family with
    \eqn{\alpha \le 1} (and arbitrary \eqn{\beta> 0})\cr
    \command{\link{RMmatern}} \tab Whittle-Matern model with
    \eqn{\nu \le 1/2}\cr
    %multiquadric todo
    %sine power todo
    \command{\link{RMstable}} \tab symmetric stable family or powered
    exponential model with \eqn{\alpha \le 1}\cr
    \command{\link{RMwhittle}} \tab Whittle-Matern model, alternative
    parametrization with \eqn{\nu \le 1/2}\cr 
  }
  
  \bold{Other isotropic models with arbitrary scale}
  \tabular{ll}{
   \command{\link{RMconstant}} \tab spatially constant model \cr
   \command{\link{RMnugget}} \tab nugget effect model \cr
  }
  
  \bold{Compactly supported covariance functions allowing for scales up to
    \eqn{\pi} (or \eqn{180} degrees)}
    \tabular{ll}{
      \command{\link{RMaskey}} \tab Askey's model\cr
      \command{\link{RMcircular}} \tab circular model\cr
     \command{\link{RMgengneiting}} \tab Wendland-Gneiting model;
      differentiable models with compact support \cr
   \command{\link{RMgneiting}} \tab differentiable model with compact
    support \cr 
    \command{\link{RMspheric}} \tab spherical model \cr
  }
  
  \bold{Anisotropic models}
  \tabular{ll}{
    None up to now.
  }
  
  \bold{Basic Operators}
  \tabular{ll}{
    \command{\link{RMmult}}, \code{*} \tab product of covariance models \cr
    \command{\link{RMplus}}, \code{+} \tab sum of covariance models
    or variograms\cr
  }

 \bold{See \link{RMmodels} for cartesian models.}
}

\me

\seealso{
  \link{coordinate systems},
  \command{\link{RMmodels}},
  \command{\link{RMtrafo}}.
}

\keyword{spatial}
\keyword{models}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFgetModelNames(isotropy=c("spherical isotropic"))

## an example of a simple model valid on a sphere
model <- RMexp(var=1.6, scale=0.5) + RMnugget(var=0) #exponential + nugget
plot(model)


## a simple simulation
l <- seq(0, 85, 1.2)
coord <- cbind(lon=l, lat=l)


z <- RFsimulate(RMwhittle(s=30, nu=0.45), coord, grid=TRUE) # takes 1 min
plot(z)


z <- RFsimulate(RMwhittle(s=500, nu=0.5), coord, grid=TRUE,
                new_coord_sys="orthographic", zenit=c(25, 25)) 
plot(z)


z <- RFsimulate(RMwhittle(s=500, nu=0.5), coord, grid=TRUE,
                new_coord_sys="gnomonic", zenit=c(25, 25)) 
plot(z)


## space-time modelling on the sphere
sigma <- 5 * sqrt((R.lat()-30)^2 + (R.lon()-20)^2)
model <- RMprod(sigma) * RMtrafo(RMexp(s=500, proj="space"), "cartesian") *
  RMspheric(proj="time") 
z <- RFsimulate(model, 0:10, 10:20, T=seq(0, 1, 0.1),
                coord_system="earth", new_coordunits="km")
plot(z, MARGIN.slices=3)


\dontshow{FinalizeExample(); }}