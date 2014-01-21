\name{jss13}
\alias{jss13}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{\R package \pkg{RandomFields}: Analysis and simulation of
 multivariate random fields and more}
\description{
 Here the code of the paper on \sQuote{\R package \pkg{RandomFields}:
 Analysis and simulation of 
 multivariate random fields and more} is given. 
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de};
}
\references{
 \itemize{ 
 \item
 Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and Strokorb, K. (2013)
 \R package \pkg{RandomFields}: Analysis and simulation of
 multivariate random fields and more. \emph{Submitted}
 }
}
 

\examples{

RFoptions(seed=0)

\dontrun{

RFoptions(seed = 0, height=4, always_close_screen=FALSE)


## Fig. 1
model <- RMmatrix(M = c(0.9, 0.43), RMwhittle(nu = 0.3)) + 
         RMmatrix(M = c(0.6, 0.8), RMwhittle(nu = 2))
x <- y <- seq(-10, 10, 0.2)
simu <- RFsimulate(model, x, y, grid=TRUE)
plot(simu)


## Fig. 2
model <- RMdelay(RMstable(alpha=1.9, scale=2), s=c(4, 4))
simu <- RFsimulate(model, x, y, grid=TRUE)
plot(simu, zlim='joint')


## Fig. 3
model <- RMdelay(RMstable(alpha=1.9, scale=2), s=c(0, 4)) + 
         RMdelay(RMstable(alpha=1.9, scale=2), s=c(4, 0))
simu <- RFsimulate(model, x, y, grid=TRUE)
plot(simu, zlim='joint')
plot(model, dim=2, xlim=c(-5, 5), main="Covariance function", 
     cex=1.25, cex.lab=1, xlab=" ")


## Fig. 4
model <- RMgencauchy(alpha=1.5, beta=3)
simu <- RFsimulate(model, x, y, z=c(0,1), grid=TRUE)
plot(simu, MARGIN.slices=3, n.slices=2)


## Fig. 5
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
simu <- RFsimulate(model, x, y, grid=TRUE)
plot(simu)


## Fig. 6
model <- RMcurlfree(RMmatern(nu=5), scale=4)
simu <- RFsimulate(model, x, y, grid=TRUE)
plot(simu, select.variables=list(1, 2:3, 4))
plot(model, dim=2, xlim=c(-3, 3), main="")


## Fig. 7
x <- y <- seq(-2, 2, len=20)
model <- RMkolmogorov()
simu <- RFsimulate(model, x, y, z=1, grid=TRUE)
plot(simu, select.variables=list(1:2), col=c("red"))
plot(model, dim=3, xlim=c(-3, 3), MARGIN=1:2, fixed.MARGIN=1, main="")



## Section 5 : Inference
## kriging takes about 14 min on one
## Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz

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



## Section 5 : Spatial prediction (kriging)
## kriging takes about 17sec on one
## Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz
RFoptions(spConform = TRUE)
## project coordinates
a <- colMeans(weather[, 3:4]) * pi / 180
P <- cbind(c(-sin(a[1]), cos(a[1]), 0),
           c(-cos(a[1]) * sin(a[2]), -sin(a[1]) * sin(a[2]), cos(a[2])),
           c(cos(a[1]) * cos(a[2]), sin(a[1]) * cos(a[2]), sin(a[2])))
x <- RFearth2cartesian(weather[, 3:4])
plane <- (x \%*\% P)[, 1:2]
r <- apply(plane, 2, range)

n <- 200
data <- cbind(plane, weather[, 1:2])
z <- RFinterpolate(p$model,  data=data, dim = 2,
                   x=seq(r[1, 1], r[2, 1], length=n),
                   y=seq(r[1, 2], r[2, 2], length=n),                  
                   variab_units=c("Pa", "K"))
plot(z)
dev.copy2pdf(file=paste(sep="/", path, "kriging.pdf"))

save(file=paste(sep="/", path, "jss.kriging.rda"), z)



## Section 5: Likelihood ratio tests
## the likelihood ratio test takes about 18 hours on one
## Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz

pars.full <- RFratiotest(pars.model, full.model, alpha=0.05,
                         distances=Dist.mat, dim = 3, data=PT)

indep.full <- RFratiotest(indep.model, full.model, alpha=0.05,
                          distances=Dist.mat, dim = 3, data=PT)

lapply(list(pars.full, indep.full), function(x) x$msg)



## Section 5: cross validation
## the cross-validation takes 35 hours on a one
## Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz
if (!exists("pred.indep")) {
  pred.indep <- pred.pars <- matrix(nc=2, nr=nrow(weather))
  res <- list()
}

fromto <- 1:nrow(weather)
for (out in fromto) {
  Print(out)
  if (all(!is.na(pred.pars[out, ]))) next
  sdP <- sd(weather[-out, 1])
  sdT <- sd(weather[-out, 2])
  PT <- cbind( weather[-out, 1] / sdP, weather[-out, 2] / sdT )
  Dist.mat <- as.vector(RFearth2dist(weather[-out, 3:4]))
  data <- cbind(x[-out, ], weather[-out, 1:2])
  
  indep <- RFfit(indep.model, distances=Dist.mat, dim = 3, data=PT, spC=FALSE)
  i <- PaperOutput(indep.model, indep$ml, sdP, sdT)#% ever best: -1276.75
  zi <- RFinterpolate(i$model, x=x[out, 1], y=x[out, 2], z=x[out, 3],
                     data = data, spC=FALSE)
  pred.indep[out, ] <- zi
  
  pars <- RFfit(pars.model, distances=Dist.mat, dim = 3, data=PT, spC=FALSE)
  p <- PaperOutput(pars.model, pars$ml, sdP, sdT) #% ever best: -1265.73
  zp <- RFinterpolate(p$model, x=x[out, 1], y=x[out, 2], z=x[out, 3],
                     data = data, spC=FALSE)
  pred.pars[out, ] <- zp

  res[[out]] <- list(indep, i, zi, pars, p, zp)
  Print(out, x[out, ], as.double(zi), as.double(zp),  weather[out, 1:2])

  save(file=paste(sep=".", "jss/pdf/cross","rda"), pred.indep, pred.pars, res)
}
apply(pred.indep - weather[, 1:2], 2, mean)     # -7.28249190  0.04199676 
apply(pred.indep - weather[, 1:2], 2, sd)       # 132.626224    1.641583 
sqrt(apply((pred.indep - weather[, 1:2])^2, 2, mean)) # 132.403601    1.636885
apply(pred.pars - weather[, 1:2], 2, mean)      # -8.40425406  0.02442434
apply(pred.pars - weather[, 1:2], 2, sd)        # 127.111010    1.588013 
sqrt(apply((pred.pars - weather[, 1:2])^2, 2, mean))  # 126.983968    1.583136 
}

\dontshow{RFoptions(seed=NA)}
}

\seealso{
 \link{weather}, \link{SS12}, \link{S10}
}

\keyword{spatial}

