\name{Papers}
\alias{papers}
\alias{Papers}
\alias{schlather}
\alias{Schlather}
\title{Code used in the papers co-authored by M. Schlather}
\description{
  Code is provided that provides the results in the
  references below.
 }
 
 \author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather};
  }
  \references{
    \itemize{
      \item
    Gneiting, T., Kleiber, W.  and  Schlather, M. (2011) Matern
    cross-covariance functions for multivariate random fields 
    \emph{J. Amer. Statist. Assoc.} \bold{105}, 1167-1177.
    
     \item
    Gneiting, T., Sevcikova, H., Percival, D.B., Schlather, M., Jiang,
    Y. (2006) Fast and Exact Simulation of Large Gaussian Lattice Systems
    in R2: Exploring the Limits. \emph{J. Comput. Graph. Stat.},
    \bold{15}, 483-501.

      \item
   Schlather, M. (2002) Models for stationary max-stable
   random fields. \emph{Extremes} \bold{5}, 33-44.

      \item
  Scheuerer, M. and Schlather, M. (2011) Covariance Models for Random
   Vector Fields.

      \item
  Schlather, M. (2010)
   On some covariance models based on normal scale mixtures.
   \emph{Bernoulli}, \bold{16}, 780-797.
 }
   
     %  Schlather, M. (2001) Simulation of stationary and isotropic random
%  fields. \emph{R-News} \bold{1} (2), 18-20.


}

\examples{

% library(RandomFields, lib="~/TMP"); runif(1)
############################################################
##                                                        ##
##     M. Schlather, Extremes 5:1, 33-44, 2002            ##
##                                                        ##
############################################################

col=grey((100:0)/100)
pts <- if (interactive()) 512 else 32
x <- (1:pts) / pts * 10
scalegauss <- 1.5
runif(1)
root <- 1/2
RFparameters(mpp.plus=2*scalegauss, CE.force=TRUE)

save.seed <- .Random.seed
ms1 <- MaxStableRF(x, x, model="gauss", param=c(0, 1, 0, scalegauss), 
                  maxstable="Bool", grid=TRUE, Print=3)
image(x, x, ms1^root, col=col)

.Random.seed <- save.seed  
scalecirc<- 3.3 # 0.16 0.00 0.18 0.00 0.00 
ms2 <- MaxStableRF(x, x, model="sph", param=c(0, 1, 0, scalecirc), 
                  maxstable="Bool", grid=TRUE, Print=3)
image(x, x, ms2^root, col=col)
  
.Random.seed <- save.seed  
ms3 <- MaxStableRF(x, x, model="exponen", param=c(0, 1, 0, 1),
           maxstable="extr", grid=TRUE, Print=3, method="ci")
image(x, x, ms3^root, col=col)

.Random.seed <- save.seed 
ms4 <- MaxStableRF(x, x, model="gauss", param=c(0, 1, 0, 1),
           maxstable="extr", grid=TRUE, Print=3, method="ci")
image(x, x, ms4^root, col=col)






############################################################
##                                                        ##
##     M. Schlather, Bernoulli, 16, 780-797, 2010         ##
##                                                        ##
############################################################


# %  library(RandomFields, lib="~/TMP")

jpg <- jpeg
Plot <- function(basicname, x, y, z, pixels=200, col=col, delay=1,
                 basicn=10000, T, redo=TRUE, height=5, zi=1, speed=0.3,
                 zlim = range(z)) {
  Print(try(dir.create(basicname)))
  if (length(T)==3) T <- seq(T[1], T[2], T[3])
  cex = 2
  par(mar=c(4,4.3, 0.8, 0.8), cex=0.7)
  for (file in list(TRUE, "jpg")) {
    h <- if (is.logical(file)) height else pixels
    args <- list(TRUE, file, ps=paste(basicname, "/X2", basicname, sep=""),
                        quiet=FALSE, height=h, width=h)
    do.call("Dev", args)
    image(x, T, z[,zi,], zlim=zlim, col=col, cex.axis=cex, cex.lab=cex,
          ylab="time")
    Dev(FALSE)
    
    args <- list(TRUE, file, ps=paste(basicname, "/X3", basicname, sep=""),
                        quiet=FALSE, height=h, width=h)
    do.call("Dev", args)
    image(x, y, z[,,1], zlim=zlim, col=col, cex.axis=cex, cex.lab=cex)
    Dev(FALSE)
  }
  k <- 0
  time <- length(T)
  if (redo) {
    system(paste("rm ", basicname, "/", basicname, "*.jpg", sep=""))

    par(mar=rep(0, 4))
    Dev(TRUE, "jpg", ps=paste(basicname, "/", basicname, sep=""),
        quiet=FALSE, height=pixels * height, width=pixels * height)
    plot(Inf, Inf, axes=FALSE, xlim=c(0,1), ylim=c(0,1))
    Dev(FALSE)
    
    for (i in 1:(time)) {
      for (j in 1:delay) {
        k <- k + 1
        name <- paste(basicname, "/", basicname, k + basicn, sep="")
        Dev(TRUE, "jpg", ps=name, quiet=FALSE, height=pixels, width=pixels)
        image(x, y, z[,,i], zlim=zlim, col=col, axes=FALSE)
        Dev(FALSE)
      }
    }
  }
  system(paste("cd ", basicname, ";mencoder -mf fps=30 -ffourcc DX50 -ovc lavc ",
               " -speed ", speed,
               " -lavcopts vcodec=mpeg4:vbitrate=9800:aspect=4/3:vhq:keyint=15",
               " -vf scale=720:576 -o ", basicname, ".avi mf://",
               basicname, "*.jpg &", sep=""))
}




### Example 10 in Schlather (2010)
basicname <- "coxisham"
y <- x <- seq(0, 10, len=if (interactive()) 256 else 5)
T <- c(0, 25.5, if (interactive()) 0.1 else 5)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)])
model <- list("coxisham", mu=c(1, 1), D=matrix(nr=2, c(1, 0.5, 0.5, 1)),
              list("whittle", nu=1)
              )
z <- GaussRF(x, y, T=T, grid =TRUE, spectral.lines=1500, every=10,
             model = model, Print=2)

zlim <- range(z)
time <- dim(z)[3]
for (i in 1:time) {
  cat(i,"\n")
  image(x, y, z[, , i], add=i>1, col=col, zlim=zlim)
} 
# load(paste(basicname, ".dat", sep=""))
Plot(basicname, x, y, z, col=col, T=T)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z)


### Example 13 in Schlather (2010)
basicname <- "moving"
y <- x <- seq(0, 10,len = if (interactive()) 256 else 10) 
T <- c(0, 25.5, if (interactive()) 0.1 else 10)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)])        
model <- list("ave1", 
              A = matrix(nc=2,c(0.5, 0, 0, 1)), z= c(2,0),
              list("whittle", nu=1)
              )
z <- GaussRF(x, x, T=T, model=model, grid=TRUE,  Print=2, me="av",
             every = 10, CE.trial=2, CE.force=TRUE, CE.maxmem=16777216*8)

zlim <- range(z)
time <- dim(z)[3]
for (i in 1:time) {
  Print(i)
  image(x, y, z[, , i], add=i>1, col=col, zlim=zlim)
}
# load(paste(basicname, ".dat", sep=""))
Plot(basicname, x, y, z, col=col, T=T)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z)



### Example 16 in Schlather (2010)
intens <- if (interactive()) 100 else 3 ## 1000 takes a huge amount
##                                       ## of time; take smaller values
basicname <- "cyclone"
len <- if (interactive()) 81 else 10
y <- x <- seq(-3, 3, len=len)
T <- seq(0, 6, len=len)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)])
model <- list("stp",  
              M=matrix(nc=3, rep(0, 9)),
              Sx=list("EtAxxA",
                E=c(1, 1, 1),
                alpha = -2 * pi,
                A=t(matrix(nc=3, c(2, 0, 0, 1, 1 , 0, 0, 0, 0)))),
              H = list("Rotat", phi= -2 * pi),
              phi = list("nonstWM", nu = 1)
              )
z <- GaussRF(x, x, z=T, model=model, grid=TRUE,  me="coin", every = 10, 
             mpp.intens=intens, mpp.p = 0.1, mpp.beta=3.5, mpp.plus = 8, Print=3) 
zlim <- c(-3.5, 3.5)
time <- dim(z)[3]
for (i in 1:time) {Print(i);image(x, y, z[,,i], add=i>1, col=col, zlim=zlim)}
# load(paste(basicname, ".dat", sep=""))
Plot(basicname, x, y, z, col=col, T=T, pixels=256, zi=0.5 + dim(z)[2]/2,
     speed=0.2, zlim=zlim)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z)




############################################################
##                                                        ##
##    T. Gneiting, W. Kleiber, M. Schlather, JASA 2011    ##
##                                                        ##
############################################################

##    see   help("weather") 



############################################################
##                                                        ##
##                     Figure 1 in                        ##
##           Gneiting, T., Sevcikova, H.,                 ##
##   Percival, D.B., Schlather, M. and Jiang, Y. (2005)   ##
##                                                        ##
############################################################

stabletest <- function(alpha, theta, size=if (interactive()) 512 else 20) {
  RFparameters(CE.trials=1, CE.tolIm = 1e-8, CE.tolRe=0, CE.force = FALSE,
               CE.useprimes=TRUE, CE.strategy=0, skipchecks=!FALSE,
               Print=2, Storing=TRUE)
  model <- list("cutoff", diameter=theta, a=1, 
                list(model="stable", alpha=alpha))
  CovarianceFct(0, model=model, dim=2)
  r <- GetModelInfo(-1)$sub[[1]]$internalq[5] # theor  R
  x <- c(0, r, r / (size - 1)) * theta
  err <- try(InitGaussRF(x, x, grid=TRUE, model=model, meth="circulant"))
  return(if (!is.numeric(err) || err != 0) NA else r)
}

alphas <- seq(1.52, 2.0, if (interactive()) 0.02 else 0.1) 
thetas <-
   if (interactive()) seq(0.05, 3.5, 0.05) else seq(0.1, 3.5, 0.2) 

m <- matrix(NA, nrow=length(thetas), ncol=length(alphas))
for (it in 1:length(thetas)) {
  theta <- thetas[it]
  for (ia in 1:length(alphas)) {
    alpha <- alphas[ia]
    cat("alpha=", alpha, "theta=", theta,"\n")
    print(unix.time(m[it, ia] <- stabletest(alpha=alpha, theta=theta)))
    if (is.na(m[it, ia])) break
  }
  if (any(is.finite(m))) image(thetas, alphas, m, col=rainbow(100))
}



############################################################
##                                                        ##
##          M. Scheuerer, M. Schlather, 2011              ##
##                                                        ##
############################################################
# %  library(RandomFields, lib="~/TMP")

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

x <- c(-3, 3, if (interactive()) 0.049 else 0.5) 
nu <- 3
col <- grey(seq(0, 1, 0.01)) 
thinning <- 21
length.arrow <- 1.5 / thinning
runif(1)
seed <- .Random.seed
eps <- interactive()  # true falls eps/pdf drucken

for (modelname in c("divfree", "curlfree")) {
  cat(modelname, "\n")
  model <- list(modelname, list("matern", nu=nu))
  xx <- seq(x[1], x[2], x[3])
  RFparameters(Print=2)

  if (!eps) {
    cf <- Covariance(x=cbind(xx,0), model=model)
    do.call(getOption("device"), list(height=5, width=5))
    par(mfcol=c(1,1))
    j <- 3
    plot(xx, cf[(j-1) * length(xx) + (1 : length(xx)), j])
  }

  .Random.seed <- seed
  z <- GaussRF(x, x,
               grid=TRUE, gridtr=TRUE, model=model, n=1, CE.trial=2,
               Stor=TRUE, me="ci", Print=3, CE.force=!TRUE)
  ## z[1,,] : Potentialfeld
  ## z[2:3,,] : vectorfeld
  ## z[4,,] : div bzw. rot
       
  if (eps) {
    do.call(getOption("device"), list(height=5, width=5))
    par(mfcol=c(1,1)) 
  } else {
    do.call(getOption("device"), list(height=5, width=10))
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




}



\keyword{spatial}


