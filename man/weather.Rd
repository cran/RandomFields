\name{weather}
\docType{data}
\alias{weather}
\alias{GKS11}
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
 Gneiting, T., Kleiber, W. and Schlather, M. (2010) Matern
 cross-covariance functions for multivariate random fields 
 \emph{J. Amer. Statist. Assoc.} \bold{105}, 1167-1177.

 \item
 Grimit, E. P. and Mass, C. F. (2002) Initial results of a
 mesoscale short-range forecasting system over the Pacific
 Northwest \emph{Wea. Forecasting} \bold{17}, 192-205.
 }
 }
\examples{

\dontshow{if (!.C("isAuthor", a=integer(1))$a) {}}

\dontshow{\dontrun{

# % library(RandomFields, lib="~/TMP")

############################################################
## ##
## T. Gneiting, W. Kleiber, M. Schlather, JASA 2010       ##
## ##
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
data(weather)
sdP <- sd(weather[, 1])
sdT <- sd(weather[, 2])
PT <- cbind( weather[, 1] / sdP, weather[, 2] / sdT )
Dist.mat <- as.vector(RFearth2dist(weather[, 3:4]))

RFoptions(spConform=FALSE)
  
nug <- RMmatrix(M=matrix(nc=2, c(NA, 0, 0, NA)), RMnugget())

#################################
## first: independent model ##
#################################
indep.model <- nug + RMbiwm(nudiag=c(NA,NA), s=c(NA, 1, NA), c=c(NA, 0, NA))
indep <- RFfit(indep.model, distances=Dist.mat, dim = 3, data=PT)

#################################
## second: parsimoninous model ## (takes a rather long time: 15 min in 2012)
#################################
pars.model <- nug + RMbiwm(nudiag=c(NA,NA), scale=NA, c=rep(NA, 3))
pars <- RFfit(pars.model, distances=Dist.mat, dim = 3, data=PT)

#################################
## third: full biwm model ##
#################################
full.model <- nug + RMbiwm(nudiag=c(NA, NA), nured=NA,
                           s=rep(NA, 3), c=rep(NA, 3))
full <- RFfit(full.model, distances=Dist.mat, dim = 3, data=PT)


#################################
## output                      ##
#################################

PaperOutput <- function(model, fit, sdP, sdT) {
  ## globally used: Dist.mat, PT
  m <- fit$model
  if (pars <- !is.null(m$C1$scale)) { ## ! parsimonious
    m$C1$phi$cdiag <- m$C1$phi$cdiag * c(sdP, sdT)^2
    m$C1$phi$c <- m$C1$phi$c * c(sdP^2, sdP * sdT, sdT^2)
    biwm <- m$C1$phi
  } else {
    m$C1$cdiag <- m$C1$cdiag * c(sdP, sdT)^2
    m$C1$c <- m$C1$c * c(sdP^2, sdP * sdT, sdT^2)
    biwm <- m$C1
  }
  m$C0$M <- m$C0$M * c(sdP, 0, 0, sdT)
  sigmaP <- sqrt(biwm$c[1])
  sigmaT <- sqrt(biwm$c[3])

  options(warn=0)
  ml <- RFfit(model, distances=Dist.mat, dim = 3, 
              modus_operandi = "careless", fit.split=FALSE,sub.methods=NULL,
              users.guess = m, only_users=TRUE,
              data=t(t(PT) * c(sdP, sdT)),
              grid=FALSE, meth="ml",
              print=1, cPr=1, spC=TRUE)["ml"]

  return(list(list(sigmaP = sigmaP,
                   sigmaT = sigmaT,
                   nuP = biwm$nu[1],
                   nuT = biwm$nu[3],
                   nuPT = if(biwm$c[2]==0) NA else biwm$nu[2],
                   inv.aP = if (pars) m$C1$s else biwm$s[1],
                   inv.aT = if (pars) m$C1$s else biwm$s[3] ,
                   inv.aPT= if (pars) m$C1$s else
                            if(biwm$c[2]==0) NA else biwm$s[2],
                   rhoPT = biwm$c[2] / (sigmaP * sigmaT), 
                   tauP = m$C0$M[1],
                   tauT = m$C0$M[4],
                   ml = ml@likelihood,
                   AICc = ml@AICc
                   ),
              model = ml
           ))
}

i <- PaperOutput(indep.model, indep$ml, sdP, sdT)#% ever best: -1276.75
p <- PaperOutput(pars.model, pars$ml, sdP, sdT) #% ever best: -1265.73
f <- PaperOutput(full.model, full$ml, sdP, sdT) #% ever best: -1265.20
table <- rbind(indep=unlist(i[[1]]), pars=unlist(p[[1]]), full=unlist(f[[1]]))
print(table, digits=3)
print(table[, ncol(table)-(1:0)], digits=8)



\dontshow{RFoptions(spConform=TRUE)}


}


}

\dontshow{RFoptions(seed=NA)}
}
\keyword{datasets}




