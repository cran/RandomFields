\name{SS12}
\alias{SS12}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{Covariance Models for Random Vector Fields}
\description{
 Here the code of the paper on \sQuote{Covariance Models for Random
 Vector Fields} is given. 
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software};
}
\references{
 \itemize{ 
 \item
 Scheuerer, M. and Schlather, M. (2012) Covariance Models for Random
 Vector Fields. \emph{Stochastic Models}, \bold{82}, 433-451.
 
 }
 
 % Schlather, M. (2001) Simulation of stationary and isotropic random
% fields. \emph{R-News} \bold{1} (2), 18-20.


}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}

\dontshow{\dontrun{
my.legend <- function(lu.x, lu.y, zlim, col, cex=1) {
 ## uses already the legend code of R-1.3.0
 cn <- length(col)
 filler <- vector("character", length=(cn-3)/2)
 legend(lu.x, lu.y, y.i=0.03, x.i=0.1, 
 legend=c(format(zlim[2], dig=2), filler,
 format(mean(zlim), dig=2), filler,
 format(zlim[1], dig=2)),
 lty=1, col=rev(col),cex=cex)
}

my.arrows <- function(xy, z, r, thinning) {
 startx <- as.vector(xy[,1] - r/2*z[as.integer(dim(z)[1]/3) + 1,,])
 starty <- as.vector(xy[,2] - r/2*z[as.integer(dim(z)[1]/3) + 2,,])
 endx <- as.vector(xy[,1] + r/2*z[as.integer(dim(z)[1]/3) + 1,,])
 endy <- as.vector(xy[,2] + r/2*z[as.integer(dim(z)[1]/3) + 2,,])
 startx[c(rep(TRUE, thinning), FALSE)] <- NA
 starty[c(rep(TRUE, thinning), FALSE)] <- NA
 endx[c(rep(TRUE, thinning), FALSE)] <- NA
 endy[c(rep(TRUE, thinning), FALSE)] <- NA
 arrows(x0=startx, y0=starty, x1=endx, y1=endy, length=0.03)
}

x <- c(-3, 3, 0.049)
nu <- 3
col <- grey(seq(0, 1, 0.01)) 
thinning <- 21
length.arrow <- 1.5 / thinning
runif(1)
seed <- .Random.seed
eps <- FALSE # true falls eps/pdf drucken
\dontshow{if(RFoptions()$internal$examples_red){warning("modified 'x'"); x <- c(-3, 3, 0.5)}}
for (modelname in c("divfree", "curlfree")) {
 cat(modelname, "\n")
 model <- list(modelname, list("matern", nu=nu))
 xx <- seq(x[1], x[2], x[3])
 RFoptions(print=2)

 if (!eps) {
 cf <- RFcov(model=model, x=cbind(xx, 0))
 do.call(getOption("device"), list(height=5, width=5))
  j <- 3
 plot(xx, cf[(j-1) * length(xx) + (1 : length(xx)), j])
 }

 .Random.seed <- seed
 z <- RFsimulate(x, x, model=model, n=1, CE.trial=2,
 Stor=TRUE, me="ci", print=3, CE.force=!TRUE)
 ## z[1,,] : Potentialfeld
 ## z[2:3,,] : vectorfeld
 ## z[4,,] : div bzw. rot
 
 if (eps) {
 ScreenDevice(height=5, width=5)
  } else {
 ScreenDevice(height=5, width=10)
 par(mfcol=c(1,2))
 }
 par(mar=c(2.2,2.5,0.5,0.5), cex.axis=2, bg="white")
 for (no.vectors in c(TRUE,FALSE)) for (i in c(1,2,4)) {
 ## image i=1: Potential feld + Vektorfeld
 ## image i=4: div/rot feld + Vektorfeld
 if (i==1 || i==4) {
 image(xx, xx, col=col, z[i,,] )
 if (no.vectors)
 my.legend(max(xx) - 0.3 * diff(range(xx)),
 min(xx) + 0.3 * diff(range(xx)),
 zlim=range(z[i,,]), col=col, cex=1.5)
 } else plot(Inf, Inf, xlim=range(xx), ylim=range(xx))

 if (!no.vectors || i!=1 && i!=4) {
 xy <- as.matrix(expand.grid(xx, xx))
 my.arrows(xy, z, length.arrow, thinning)
 }

 if (all(par()$mfcol==1)) {
 name <- paste(modelname, "_", nu, "_", !no.vectors, "_", i, sep="")
 cat(name,"\n")
 dev.copy2eps(file=paste(name, ".eps", sep=""))
 dev.copy2pdf(file=paste(name, ".pdf", sep=""))
 }
 }
}

par(mfcol=c(1,1))

}}

\dontshow{FinalizeExample()}
}


\keyword{spatial}

