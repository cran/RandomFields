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
## output ##
#################################

PaperOutput <- function(m, sdP, sdT) {
 ## correct for miles2l=km and for sdP and sdT
 if (pars <- !is.null(m$C1$scale)) { ## ! parsimonious
 m$C1$scale <- m$C1$scale * miles2km 
 m$C1$phi$cdiag <- m$C1$phi$cdiag * c(sdP, sdT)^2
 m$C1$phi$c <- m$C1$phi$c * c(sdP^2, sdP * sdT, sdT^2)
 biwm <- m$C1$phi
 } else {
 m$C1$cdiag <- m$C1$cdiag * c(sdP, sdT)^2
 m$C1$c <- m$C1$c * c(sdP^2, sdP * sdT, sdT^2)
 m$C1$s <- m$C1$s * miles2km
 biwm <- m$C1
 }
 m$C0$M <- m$C0$M * c(sdP, 0, 0, sdT)
 
 ml <- RFfit(distances=Dist.mat * miles2km,
 dim = 2, data=t(t(PT) * c(sdP, sdT)),
 model=m, grid=FALSE, meth="ml",
 spC=FALSE)$ml$ml
 
 sigmaP <- sqrt(biwm$c[1])
 sigmaT <- sqrt(biwm$c[3])
 
 return(list(#model = m,
 sigmaP = sigmaP,
 sigmaT = sigmaT,
 nuP = biwm$nu[1],
 nuT = biwm$nu[3],
 nuPT = if(biwm$c[2]==0) NA else biwm$nu[2],
 inv.aP = if (pars) m$C1$s else biwm$s[1],
 inv.aT = if (pars) m$C1$s else biwm$s[3] ,
 inv.aPT= if (pars) m$C1$s else if(biwm$c[2]==0) NA else biwm$s[2],
 rhoPT = biwm$c[2] / (sigmaP * sigmaT), 
 tauP = m$C0$M[1],
 tauT = m$C0$M[4],
 ml = ml
 ))
}

 

#################################
## get the data ##
#################################
library(fields)
miles2km <- 1.608

data(weather)
## P and T are so different in scale so that they have
## to be normalized first:
sdP <- sd(weather[, 1])
sdT <- sd(weather[, 2])
PT <- cbind( weather[, 1] / sdP, weather[, 2] / sdT )

Dist.mat <- rdist.earth(weather[, 3:4])
Dist.mat <- Dist.mat[lower.tri(Dist.mat)] ## in miles

RFoptions(use_naturalscaling=FALSE)
 
nug <- RMmatrix(M=matrix(nc=2, c(NA, 0, 0, NA)), RMnugget())

modes <- c("normal", "easy", "sloppy")

for (mode in modes) {
 cat("mode\n")
 file <- paste(mode, ".rda", sep="")
 #% if (FileExists(file)) next;

 #################################
 ## first: independent model ##
 #################################
 indep.model <- nug + RMbiwm(nudiag=c(NA,NA), s=c(NA, 1, NA), c=c(NA, 0, NA))
 indep <- RFfit(indep.model, distances=Dist.mat, dim=2, data=PT,
                meth="ml", modus_op=mode, spC=FALSE)
 
 #################################
 ## second: parsimoninous model ## (takes a rather long time: 15 min in 2012)
 #################################
 pars.model <- nug + RMbiwm(nudiag=c(NA,NA), scale=NA, c=rep(NA, 3))
 pars <- RFfit(pars.model, distances=Dist.mat, dim = 2, data=PT,
               meth="ml", modus_op=mode, spC=FALSE)

 #################################
 ## third: full biwm model ##
 #################################
 full.model <- nug + RMbiwm(nu=rep(NA, 3), s=rep(NA, 3), c=rep(NA, 3))
 full <- RFfit(full.model, distances=Dist.mat, dim = 2, data=PT,
               meth="ml", modus_op=mode, spC=FALSE)

 #################################
 ## output ##
 #################################

 i <- PaperOutput(indep$ml$model, sdP, sdT)#% ever best: -1276.75
 p <- PaperOutput(pars$ml$model, sdP, sdT) #% ever best: -1265.73
 f <- PaperOutput(full$ml$model, sdP, sdT) #% ever best: -1265.30
 table <- rbind(indep=unlist(i), pars=unlist(p), full=unlist(f))
 print(table, digits=3)
 print(table[, ncol(table)], digits=8)

 save(file=file, indep, pars,full, full.model, pars.model, indep.model,
 sdP, sdT, table, modes)


 file <- paste(mode, ".tex", sep="")
 model.names <- c("indep.{}", "pars.{}", "full")
 write(file=file,
 paste("\\begin{table}[p]\\small\n\\centering\n\\caption{Comparison ",
 "between the bivariate Whittle-Mat\'ern models 'independent',",
 " 'parsimonious' and 'full'. Here $\\sigma_P=\\sqrt{c_{pp}}$ and",
 " $\\nu_P$ instead of $\\nu_{PP}$ etc., for short.", 
 " Optimization mode is '", mode, "'. Here, \\label{tab:", mode,
 "}}\n\\begin{tabular}{l",
 paste(rep("c", ncol(table)), collapse=""), "}\n",
 " & $\\sigma_P$ & $\\sigma_T$ & $\\nu_P$ & $\\nu_T$ & ",
 "$\\nu_{PT}$ & $s_P$ & $s_T$ & $s_{PT}$ & $\\rho_{PT}$ &",
 " $\\tau_P$ & $\\tau_T$ & ML \\\\\\hline",
 sep="")
 )
 for (nr in 1:length(modes)) {
 write(file=file, append=TRUE,
 paste(model.names[nr], " & ",
 paste(signif(table[nr, (0:1) - ncol(table)], 3),
 collapse=" & "),
 " & ", signif(table[nr, ncol(table)-1], 1),
 " & ", signif(table[nr, ncol(table)], 7), "\\\\"))
 }
 write(file=file, append=TRUE, "\\end{tabular}\n\\end{table}\n", sep="")
}


}


}
}
\keyword{datasets}




