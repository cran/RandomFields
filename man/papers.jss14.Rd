\name{jss14}
\alias{jss14}
\alias{information}
\alias{pars}
\alias{pars.model}
\alias{whole}
\alias{whole.model}
\title{Covariance models for multivariate and vector valued fields}
\description{
 Here the code of the paper on \sQuote{Analysis, simulation and prediction of
 multivariate random fields with package RandomFields}
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software};
}
\references{
 \itemize{ 
 \item
  Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and Strokorb, K. (2014)
 Analysis, simulation and prediction of
 multivariate random fields with package \pkg{RandomFields} \emph{JSS}
 \emph{Submitted}

 }
}

 

\examples{

\dontshow{
if (!interactive()) RFoptions(modus_operandi="sloppy")
}

\donttest{

# This skript is based on version 3.0.28 available at
# http://ms.math.uni-mannheim.de/de/publications/software/



               ###########################################
               ##  SECTION 4: UNCONDITIONAL SIMULATION  ##
               ###########################################


RFoptions(seed = 0, height = 4, always_close_screen = FALSE)
## seed=0:  *ANY* simulation will have the random seed 0; set
##          RFoptions(seed=NA) to make them all random again
## height : height of X11
## always_close_screen=FALSE: the pictures are kept on the screen


## Fig. 1: linear model of coregionalization
M1 <- c(0.9, 0.6)
M2 <- c(sqrt(0.19), 0.8)
model <- RMmatrix(M = M1, RMwhittle(nu = 0.3)) + 
         RMmatrix(M = M2, RMwhittle(nu = 2))
x <- y <- seq(-10, 10, 0.2)
simu <- RFsimulate(model, x, y)
plot(simu)


## Fig. 2: Wackernagel's delay model
model <- RMdelay(RMstable(alpha = 1.9, scale = 2), s = c(4, 4))
simu <- RFsimulate(model, x, y)
plot(simu, zlim = 'joint')
%dev.copy2pdf(file="delay0.pdf")

## Fig. 3: extended Wackernagel's delay model
model <- RMdelay(RMstable(alpha = 1.9, scale = 2), s = c(0, 4)) + 
         RMdelay(RMstable(alpha = 1.9, scale = 2), s = c(4, 0))
simu <- RFsimulate(model, x, y)
plot(simu, zlim = 'joint')
plot(model, dim = 2, xlim = c(-5, 5), main = "Covariance function", 
     cex = 1.5, col = "brown")
%dev.copy2pdf(file="delay_cov.pdf")

## Fig. 4: latent dimension model
## MARGIN.slices has the effect of choosing the third dimension
##               as latend dimension
## n.slices has the effect of choosing a bivariate model
model <- RMgencauchy(alpha = 1.5, beta = 3)
simu <- RFsimulate(model, x, y, z = c(0,1))
plot(simu, MARGIN.slices = 3, n.slices = 2)


## Fig. 5: Gneiting's bivariate Whittle-Matern model
model <- RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 1, cdiag = c(1, 5), 
                s = c(1, 1, 2))
simu <- RFsimulate(model, x, y)
plot(simu)


## Fig. 6: anisotropic linear model of coregionalisaton
M1 <- c(0.9, 0.6)
M2 <- c(sqrt(0.19), 0.8)
A1 <- RMangle(angle = pi/4, diag = c(0.1, 0.5))
A2 <- RMangle(angle = 0, diag = c(0.1, 0.5))
model <- RMmatrix(M = M1, RMgengneiting(kappa = 0, mu = 2, Aniso = A1)) +
         RMmatrix(M = M2, RMgengneiting(kappa = 3, mu = 2, Aniso = A2))
x <- y <- seq(-10, 10, 0.2)
simu <- RFsimulate(model, x, y)
plot(simu)


## Fig. 7: random vector field whose path are curl free
## A 4-variate field is simulated, where the first component
## refers to the potential field, the second and third component
## to the curl free vector field and the forth component to the
## field of sinks and sources
model <- RMcurlfree(RMmatern(nu = 5), scale = 4)
simu <- RFsimulate(model, x, y)
plot(simu, select.variables = list(1, 2 : 3, 4))
plot(model, dim = 2, xlim = c(-3, 3), main = "", cex = 2.3, col="brown") 
%dev.copy2pdf(file="curlfree_cov.pdf")


## Fig. 8: Kolmogorov's model of turbulence
## See the physical literature for its derivation and details
x <- y <- seq(-2, 2, len = 20)
model <- RMkolmogorov()
simu <- RFsimulate(model, x, y, z = 1)
plot(simu, select.variables = list(1 : 2), col = c("red"))
plot(model, dim = 3, xlim = c(-3, 3), MARGIN = 1 : 2, cex = 2.3,
     fixed.MARGIN = 1.0, main = "", col = "brown")
%dev.copy2pdf(file="kolmogorov_cov.pdf")





               ###########################################
               ## SECTION 5: DATA ANALYSIS              ##
               ###########################################

## get the data     
data("weather")
PT <- weather[ , 1 : 2]
## transformation of earth coordinates to Euclidean coordinates
Dist.mat <- as.vector(RFearth2dist(weather[ , 3 : 4]))


#################################
## model definition            ##
#################################
## bivariate pure nugget effect:
nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())
## parsimonious bivariate Matern model
pars.model <- nug + RMbiwm(nudiag = c(NA, NA), scale = NA, cdiag = c(NA, NA),
                           rhored = NA)
## whole bivariate Matern model
whole.model <- nug + RMbiwm(nudiag = c(NA, NA), nured = NA, s = rep(NA, 3),
                            cdiag = c(NA, NA), rhored = NA)



#################################
## model fitting and testing   ## 
#################################
## 'parsimoneous model'
## fitting takes a while; 'pars' is already available within 'data(weather)'
pars <- RFfit(pars.model, distances = Dist.mat, dim = 3, data = PT)
print(pars)
print(pars, full=TRUE)
RFratiotest(pars)
RFcrossvalidate(pars, x = weather[ , 3 : 4], data = PT, full = TRUE)

## 'whole model'
## fitting takes a while; 'whole' is already available within 'data(weather)'
whole <- RFfit(whole.model, distances = Dist.mat, dim = 3, data = PT)
print(whole, full=TRUE)
RFratiotest(whole)
RFcrossvalidate(whole, x = weather[ , 3 : 4], data = PT, full = TRUE)

## compare parsimonous and whole
RFratiotest(nullmodel = pars, alternative = whole)


#################################
## kriging                     ##
#################################
## First the coordindates are projected on a plane
a <- colMeans(weather[ , 3 : 4]) * pi / 180
P <- cbind(c(-sin(a[1]), cos(a[1]), 0),
           c(-cos(a[1]) * sin(a[2]), -sin(a[1]) * sin(a[2]), cos(a[2])),
           c(cos(a[1]) * cos(a[2]), sin(a[1]) * cos(a[2]), sin(a[2])))
x <- RFearth2cartesian(weather[ , 3 : 4])
plane <- (x \%*\%P)[ , 1 : 2]

## here, kriging is performed on a rectangular that covers the 
## the projected points above. The rectangular is approximated
## by a grid of length 200 in each direction.
n <- 200
r <- apply(plane, 2, range)
data <- cbind(plane, weather[ , 1 : 2])
z <- RFinterpolate(pars, data = data, dim = 2,
                   x = seq(r[1, 1], r[2, 1], length = n),
                   y = seq(r[1, 2], r[2, 2], length = n),
                   variab_units = c("Pa", "K"), spConform = TRUE)
plot(z)

}

\dontshow{\dontrun{ # OK
dip <- dimnames(installed.packages())
version <- installed.packages()[which(dip[[1]] == "RandomFields")[1],
                                dip[[2]] == "Version"]
information <- list(date=date(), version=version) # version auf /usr/local !

save(file="/home/schlather/svn/RandomFields/RandomFields/data/weather.rda",# OK
     pars.model, pars, whole.model, whole,
     weather, information=information)
}}


\dontshow{RFoptions(modus_operandi="normal")}
\dontshow{FinalizeExample()}
}

\seealso{
 \link{weather}, \link{SS12}, \link{S10}
}

\keyword{spatial}

