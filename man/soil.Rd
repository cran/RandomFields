\name{soil}
\alias{soil}
\title{Soil data of North Bavaria, Germany}
\usage{data(soil)}
\description{
  Soil physical and chemical data collected on a field in the
  Weissenstaedter Becken, Germany
}
\format{
This data frame contains the following columns:
\describe{
\item{x.coord}{x coordinates given in cm}
\item{y.coord}{y coordinates given in cm}
\item{nr}{number of the samples, which were taken in this order}
\item{moisture}{moisture content [Kg/Kg * 100\%]}
\item{NO3.N}{nitrate nitrogen [mg/Kg]}
\item{Total.N}{total nitrogen [mg/Kg]}
\item{NH4.N}{ammonium nitrogen [mg/Kg]}
\item{DOC}{dissolved organic carbon [mg/Kg]}
\item{N20N}{nitrous oxide [mg/Kg dried substance]}
}
}
\details{
  For technical reasons some of the data were obtained as differences
  of two measurements (which are not available anymore). Therefore,
  some of the data have negative values.
}
\source{
  The data were collected by Wolfgang Falk,
  Soil Physics Group,
  \url{http://www.geo.uni-bayreuth.de/bodenphysik/Welcome.html},
  University of Bayreuth, Germany.
}
\references{
  Falk, W. (2000)
  \emph{Kleinskalige raeumliche Variabilitaet von Lachgas und bodenchemischen
  Parameters [Small Scale Spatial Variability of Nitrous Oxide and
  Pedo-Chemical Parameters]},
  Master thesis, University of Bayreuth, Germany.
}
\examples{
% library(RandomFields);
% source("~/article/R/NEW.RF/RandomFields/R/rf.R")
% source("~/article/R/NEW.RF/RandomFields/R/modelling.R")

################################################################
##                                                            ##
##         a geostatistical analysis that demonstrates        ##
##         features of the package `RandomFields'             ##
##                                                            ##
##        *** THIS EXAMPLE IS NOT RUN AUTOMATICALLY ***       ##
##        ***      SINCE IT IS TIME CONSUMING       ***       ##
##                                                            ##
################################################################

\dontrun{

data(soil)
names(soil)
pts <- soil[,c(1,2)]
d <- soil$moisture

## define some graphical parameters first
close.screen(close.screen())
maxbin <- max(dist(pts)) / 2
(zlim <- range(d))
cn <- 100
colour <- rainbow(cn)
par(bg="white",cex=1, cex.lab=1.4, cex.axis=1.4, mar=c(4.3,4.3,0.8,0.8))
lu.x <- min(soil$x)
lu.y <- max(soil$y)
y <- x <- seq(min(soil$x), max(soil$x), l=121) 

## ... and a certain appearance of the legend
my.legend <- function(lu.x, lu.y, zlim, col, cex=1) {
  ## uses already the legend code of R-1.3.0
  cn <- length(col)
  filler <- vector("character", length=(cn-3)/2)
  legend(lu.x, lu.y, y.i=0.03, x.i=0.1, 
         legend=c(format(zlim[2], dig=2), filler,
         format(mean(zlim), dig=2), filler,
         format(zlim[1], dig=2)),
         lty=1, col=col[length(col):1],cex=cex)
}

## plot the data first
plot(pts, col=colour[1+99*((d-min(d))/diff(zlim))], pch=16,
     xlab="x [cm]", ylab="y [cm]")
my.legend(lu.x, lu.y, zlim=zlim, col=colour, cex=1.4)
% dev.copy2eps(file="data.eps", width=5, height=5)

% source("/home/schlather/article/R/NEW.RF/RandomFields/R/rf.R")
## empirical variogram
ev <- EmpiricalVariogram(pts, data=d, grid=FALSE,
                         bin=c(-1,seq(0,maxbin,l=30)))

## show all models,
%source("/home/schlather/article/R/NEW.RF/RandomFields/R/ShowModels.R")
by.eye <- ShowModels(x=x, y=y, emp=ev, col=colour, zlim=zlim,
                     Mean=mean(d))
%by.eye<-list(param=c(13.23566, 6.1301215, 2.1624282, 19.9774468, 0.5268569))
%dev.copy2eps(file="ShowModels.eps", width=7, height=7)

## fit parameters of the whittlematern model by MLE
%source("/home/schlather/article/R/NEW.RF/RandomFields/R/MLE.R")
mle <- fitvario(x=pts, data=d, model="whittle", par=rep(NA,5))

## plot the fitted model and the empirical variogram
plot(ev$c, ev$emp.var, ylim=c(0,11), ylab="variogram", xlab="lag")
gx <- seq(0.001, max(ev$c), l=100)
if(!is.null(by.eye)) lines(gx, Variogram(gx, model=by.eye)) 
lines(gx, Variogram(gx, model=mle$mle), col=2)
lines(gx, Variogram(gx, model=mle$lsq), col=3)
lines(gx, Variogram(gx, model=mle$nlsq), col=4)
lines(gx, Variogram(gx, model=mle$slsq), col=5)
legend(120, 4, c("empirical", "by eye", "MLE", "lsq", "sqrt(n) lsq",
               "sd^-0.5 lsq"),
       lty=c(-1, rep(1, 5)), pch=c(1, rep(-1, 5)),
       col=c(1, 1, 2, 3, 4, 5), cex=1.4)
%dev.copy2eps(file="mle.eps", width=5, height=5)

## map of expected values
%source("/home/schlather/article/R/NEW.RF/RandomFields/R/modelling.R")
k <- Kriging("O", x=x, y=y, grid=TRUE, model=mle$mle, given=pts, data=d)
image(x, y, k, col=colour, zlim=zlim, xlab="x [cm]", ylab="y [cm]")
par(bg="white")
my.legend(lu.x, lu.y, zlim=zlim, col=colour, cex=1.4)
%dev.copy2eps(file="kriging.eps", width=5, height=5)

%source("/home/schlather/article/R/NEW.RF/RandomFields/R/modelling.R")
## what is the probability that at no point of the
## grid given by x and y the moisture is greater than 24 percent?
cs <- CondSimu("O", x=x, y=y, grid=TRUE, model=mle$mle, given=pts,
               data=d, n=10) # better n=100 or n=1000
%!!!!!!!!! cs ist keine MATRIX !!!!!! ?? 1.10.03
image(x, y, cs[,,1], col=colour, zlim=zlim, xlab="x [cm]",
      ylab="y [cm]")
my.legend(lu.x, lu.y, zlim=zlim, col=colour, cex=1.4)
%dev.copy2eps(file="condsimu.eps", width=5, height=5)
mean(apply(cs<=24, 3, all)) ## about 40 percent ...

}% dontrun

} 

################################################################
##                                                            ##
##        *** THIS EXAMPLE IS NOT RUN AUTOMATICALLY ***       ##
##        ***      SINCE IT IS TIME CONSUMING       ***       ##
##                                                            ##
################################################################
\keyword{datasets}






