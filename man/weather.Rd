\name{weather}
\alias{weather}
\title{Pressure and temperature forecast errors over the Pacific Northwest}
\usage{data(weather)}
\description{
  Meteorological dataset, which consists of difference between forecasts
  and observations (forcasts minus observations) of temperature and
  pressure at 157 locations in the North American Pacific Northwest. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{pressure}{in units of Pascal}
    \item{temperature}{in units of degree Celcius}
    \item{lon}{longitudinal coordinates of the locations}
    \item{lat}{latitude coordinates of the locations}
  }
}
\details{
  The forecasts are from the GFS member of the University of Washington
regional numerical weather prediction ensemble (UWME; Grimit and Mass
2002; Eckel and Mass 2005); they were valid on December 18, 2003 at 4 pm
local time, at a forecast horizon of 48 hours.
}
\source{
The data were obtained from Cliff Mass and Jeff Baars in the University
of Washington Department of Atmospheric Sciences.
}
\references{
  \itemize{
    \item
    Eckel, A. F. and Mass, C. F. (2005) Aspects of effective mesoscale,
   short-range ensemble forecasting \emph{Wea. Forecasting} \bold{20},
   328-350.

   \item
    Gneiting, T., Kleiber, W.  and  Schlather, M. (2010) Matern
    cross-covariance functions for multivariate random fields 
    \emph{J. Amer. Statist. Assoc.} \bold{105}, 1167-1177.

    \item
    Grimit, E. P. and Mass, C. F. (2002) Initial results of a
    mesoscale short-range forecasting system over the Pacific
    Northwest \emph{Wea. Forecasting} \bold{17}, 192-205.
  }
  }
\examples{

\dontrun{

# % library(RandomFields, lib="~/TMP")

############################################################
##                                                        ##
##    T. Gneiting, W. Kleiber, M. Schlather, JASA 2010    ##
##                                                        ##
############################################################

## Here the analysis in the above mentioned paper is replicated.
## The results differ slightly from those in the paper, as the algorithm
## has been improved, and the estimation has been nearly fully automated.
## In particular, user supplied lower and upper bound for estimating
## the independent model are no longer required.
## As a result, the corresponding MLE fit deteriorates, whereas
## the other fits improve slightly.




#################################
## get the data                ##
#################################
library(fields)
miles2km <- 1.608

data(weather)
## P and T are so different in scale so that they have
## to be normalized first:
sdP <- sd(weather[, 1])
sdT <- sd(weather[, 2])
PT <- cbind( weather[, 1] / sdP, weather[, 2] / sdT )

Dist.mat <- rdist.earth(weather[,3:4])
Dist.mat <- Dist.mat[lower.tri(Dist.mat)] ## in miles


#################################
## first: marginal estimation  ##
#################################
uni.est <- 
  list("+",
       list("$", var=NA, list("nugget")),
       list("$", var=NA, scale=NA, list("whittle", nu=NA))
       )
fvP <- fitvario(Distances=Dist.mat, truedim=2, data=PT[, 1],
                model=uni.est, grid=FALSE, ml="ml", Print=1)$ml # -153.9

    
fvT <- fitvario(Distances=Dist.mat, truedim=2, data=PT[, 2],
                model=uni.est, grid=FALSE, ml="ml", Print=1)$ml # -138.45

# % file <- "fvPT.rda"; load(file); fvT <- fvT$ml; fvP  <- fvP$ml


#################################
## second: parsimoninous model ##
#################################
puni2bi <- function(mP, mT, lower) {
  nugg1 <- mP[[2]]$var
  nugg2 <- mT[[2]]$var
  nu1 <- mP[[3]][[4]]$nu
  nu2 <- mT[[3]][[4]]$nu
  s <- mean(c(mP[[3]]$scale, mT[[3]]$scale))
  c1 <- mP[[3]]$var
  c2 <- mT[[3]]$var
  if (missing(lower)) {
    rhored <- 0
    f <- 1
  } else if (lower) {
    rhored <- -1
    f <- 0.2
    nugg1 <- nugg2 <- 0
  } else {
    rhored <- 1
    f <- 4
    nugg1 <- c1 / 2
    nugg2 <- c2 / 2
  }
  return(list("+",
              list("M", M=matrix(nc=2, c(sqrt(nugg1), 0, 0, sqrt(nugg2))),
                   list("nugget")),
              list("parsbiWM",
                   nu = c(nu1 * f, nu2 * f), 
                   s = s * f,
                   c = c(c1 * f, c2 * f), rhored=rhored 
                   )
              ) )
}

p.est.model <- 
  list("+",
       list("M", M=matrix(nc=2, c(NA, 0, 0, NA)),
            list("nugget")),
       list("parsbiWM",
            nu = c(NA, NA), 
            s = NA,
            c = c(NA, NA), rhored=NA
            ))

## takes a rather long time (15 min in 2011)
pars <- fitvario(Distances=Dist.mat, truedim = 2, data=PT,
                 model=p.est.model, grid=FALSE, trend=list(0,0),
                 lower= puni2bi(fvP$model, fvT$model, lower=TRUE),
                 upper= puni2bi(fvP$model, fvT$model, lower=FALSE),
                 users.guess=puni2bi(fvP$model, fvT$model),
                 Print=2, ml="ml")$ml ## -280.9791

# % save(file="pars", pars)
# % Print(puni2bi(fvP$model, fvT$model, lower=TRUE),puni2bi(fvP$model, fvT$model, lower=FALSE),puni2bi(fvP$model, fvT$model))


#################################
## third: full biwm model      ##
#################################
pars2full <- function(model, lower){
  s <- model[[3]]$s
  nuP <- model[[3]]$nu[1]
  nuT <- model[[3]]$nu[2]
  tauP <- model[[2]]$M[1]
  tauT <- model[[2]]$M[4]
  cP <- model[[3]]$c[1]
  cT <- model[[3]]$c[2]
  Rhored <- model[[3]]$rhored
  nured <- 1

  if (missing(lower)) {
    f <- 1
  } else if (lower) {
    nured <- 0.1
    f <- 0.5
    Rhored <- -1
    tauP <- tauT <- 0
  } else {
    f <- 2
    tauP <- max(f * tauP, cP / 10)
    tauT <- max(f * tauT, cT / 10)
    Rhored <- 0 ## as estimated Rhored is negativ
  }

  list("+",
       list("M", M=matrix(nc=2, c(tauP, 0, 0, tauT)),
            list("nugget")),
       list("biWM",
            nu = c(nuP * f, nuT * f), nured=nured,
            s = c(s * f, s * f), s12=s * f,
            c = c(cP * f, cT * f), rhored=Rhored 
            ) )
}

est.model <- 
  list("+",
       list("M", M=matrix(nc=2, c(NA, 0, 0, NA)),
            list("nugget")),
       list("biWM",
            nu = c(NA, NA), nured=NA,
            s = c(NA, NA), s12=NA,
            c = c(NA, NA), rhored=NA
            ))

fv <- fitvario(Distances=Dist.mat, truedim = 2, data=PT, Print=2,
               model=est.model, grid=FALSE, trend=list(0,0),
               lower=pars2full(pars$model, lower=TRUE),
               upper=pars2full(pars$model, lower=FALSE),
               users.guess=pars2full(pars$model), ml="ml")$ml # -280.6910
## %                      now even  -280.5419
# %  Print(pars2full(pars$model, lower=TRUE),pars2full(pars$model,lower=FALSE),pars2full(pars$model))


#################################
## output                      ##
#################################


Factor <- function(nu1, nu2, nu12, a1, a2, a12, d) {
  gamma(nu1 + d/2) * gamma(nu2 + d/2) / gamma(nu1) / gamma(nu2) *
  (gamma(nu12) / gamma(nu12+d/2))^2 * (a1^nu1 * a2^nu2 / a12^(2*nu12) )^2
}


InfQ <- function(nu1, nu2, nu12, a1, a2, a12, d) {
  gamma <- (2*nu12+d) * a1^2 * a2^2 - (nu2 +d/2)*a1^2*a12^2 -
    (nu1 +d/2) *a2^2*a12^2
  beta <- (2 * nu12 - nu2 + d/2) * a1^2 + (2 * nu12 - nu1 + d/2) * a2^2  -
    (nu1 + nu2 + d) *a12^2
  alpha <- 2 * nu12 - nu1 -nu2
  rednu12 <- 0.5 * (nu1 + nu2) / nu12
 
  if (rednu12 == 1.0) {
    t1sq <- t2sq <- max(0, if (beta==0) 0 else -gamma / beta)
    inf <- 1
  } else {
    inf <- Inf
    discr <- beta^2 - 4 * alpha * gamma
    t1sq <- t2sq <- 0
    if (discr >= 0) {
      discr <- sqrt(discr)
      t1sq = max(0, (-beta + discr) / (2.0 * alpha))
      t2sq = max(0, (-beta - discr) / (2.0 * alpha))
    }  
  }
  tsq <- c(0, t1sq, t2sq)
  return(min(inf, (a12^2 + tsq)^(2.0 * nu12 + d) /
            (a1^2 + tsq)^(nu1 + d/2) / (a2^2 + tsq)^(nu2 + d/2) ))
}


uni2full <- function(mP, mT) {
  nuggP <- mP[[2]]$var
  nuggT <- mT[[2]]$var
  nuP <- mP[[3]][[4]]$nu
  nuT <- mT[[3]][[4]]$nu
  sP <- mP[[3]]$scale
  sT <- mT[[3]]$scale
  c1 <- mP[[3]]$var
  c2 <- mT[[3]]$var
  return(list("+",
             list("M", M=matrix(nc=2, c(sqrt(nuggP), 0, 0, nuggT)),
                 list("nugget")),
             list("biWM",
                  nu = c(nuP, nuT), nured=1,
                  s = c(sP, sT), s12=1,
                  c = c(c1, c2), rhored=0)
            ) )
}




PaperOutput <- function(m, sdP, sdT) {
  m[[2]]$M <- m[[2]]$M * c(sdP, 0, 0, sdT)
  m[[3]]$c <- m[[3]]$c * c(sdP, sdT)^2
  sP <- m[[3]]$s[1] * miles2km
  sT <- m[[3]]$s[2] * miles2km
  sPT <- m[[3]]$s12 * miles2km
  m[[3]]$s <- c(sP, sT)
  m[[3]]$s12 <- sPT
  d <- 2

  ml <- fitvario(Distances=Dist.mat * miles2km,
                 truedim = d, data=t(t(PT) * c(sdP, sdT)),
                 m=m, grid=FALSE, trend=list(0,0), ml="ml")$ml$ml

  nuP <- m[[3]]$nu[1]
  nuT <- m[[3]]$nu[2]
  nuPT <- 0.5 * (nuP + nuT) / m[[3]]$nured
 
  f <- Factor(nuP, nuT, nuPT, 1/sP, 1/sT, 1/sPT, d)
  infQ <- InfQ(nuP, nuT, nuPT, 1/sP, 1/sT, 1/sPT, d)
  return(list(
         model = m,
         sigmaP = sqrt(m[[3]]$c[1]),
         sigmaT = sqrt(m[[3]]$c[2]),
         nuP  = nuP,
         nuT  = nuT,
         nuPT = nuPT,
         inv.aP = sP,
         inv.aT = sT,
         inv.aPT= sPT,
         rhoPT  = m[[3]]$rhored * sqrt(f * infQ),
         tauP = m[[2]]$M[1],
         tauT = m[[2]]$M[4],
         ml = ml
  ))
}


Print(PaperOutput(uni2full(fvP$model, fvT$model), sdP, sdT)) # -1277.11
Print(PaperOutput(pars2full(pars$model), sdP, sdT)) # -1265.73
Print(PaperOutput(fv$model, sdP, sdT)) # -1265.30


}
}
\keyword{datasets}




