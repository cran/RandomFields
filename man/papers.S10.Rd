\name{S10}
\alias{S10}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{On some covariance models based on normal scale mixtures}
\description{
 Here the code of the paper on \sQuote{On some covariance models based on normal scale mixtures} is given.
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software};
}
\references{
 \itemize{ 
 \item
 Schlather, M. (2010)
 On some covariance models based on normal scale mixtures.
 \emph{Bernoulli}, \bold{16}, 780-797.
 }

}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

### Example 10 in Schlather (2010).
## The field below has more than 80 million points. So the simulation
## takes a while
y <- x <- seq(0, 10, len=if (interactive()) 256 else 5)
T <- c(0, 0.02, if (interactive()) 1275 else 4)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)])
model <- RMcoxisham(mu=c(1, 1), D=matrix(nr=2, c(1, 0.5, 0.5, 1)),
                    RMwhittle(nu=1))
z <- RFsimulate(model, x, y, T=T, sp_lines=1500, every=10)
plot(z, MARGIN.slices=3, col=col)




\dontshow{\dontrun{

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
  system(paste("cd ", basicname,
               ";mencoder -mf fps=30 -ffourcc DX50 -ovc lavc ",
               " -speed ", speed,
               " -lavcopts vcodec=mpeg4:vbitrate=9800:aspect=4/3:vhq:keyint=15",
               " -vf scale=720:576 -o ", basicname, ".avi mf://",
               basicname, "*.jpg &", sep=""))

}

basicname <- "coxisham"
# load(paste(basicname, ".dat", sep=""))
Plot(basicname, x, y, z, col=col, T=T)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z)

}}


\dontshow{\dontrun{
### Example 13 in Schlather (2010)
y <- x <- seq(0, 10,len = if (interactive()) 256 else 10) 
T <- c(0, 0.02, if (interactive()) 1275 else 4)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)]) 
model <- RMave(A = matrix(nc=2,c(0.5, 0, 0, 1)), z= c(2,0),
               RMwhittle(nu=1))
z <- RFsimulate(model, x, y, T=T, every = 10,
                trials=2, force=TRUE, maxmem=16777216*8)
plot(z)
}}

\dontshow{\dontrun{
basicname <- "moving"
# load(paste(basicname, ".dat", sep=""))
Plot(basicname, x, y, z, col=col, T=T)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z)
}}


\dontshow{\dontrun{

### Example 16 in Schlather (2010)
intens <- if (interactive()) 100 else 3 ## 1000 takes a huge amount
## ## of time; take smaller values
len <- if (interactive()) 81 else 4
y <- x <- seq(-3, 3, len=len)
T <- seq(0, 0.075, len=len)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)])
model <- RMstp(M=matrix(nc=3, rep(0, 9)),
               S=RMetaxxa(E=c(1, 1, 1), alpha = -2 * pi,
                          A=t(matrix(nc=3, c(2, 0, 0, 1, 1 , 0, 0, 0, 0)))),
               Aniso = RMrotation(phi= -2 * pi),
               phi = RMwhittle(nu = 1) )
z <- RFsimulate(model, x, x, z=T, me="coin", every = 10, 
 mpp.intens=intens, mpp.p = 0.1, mpp.beta=3.5, mpp.plus = 8, Print=3) 
zlim <- c(-3.5, 3.5)
time <- dim(z)[3]
for (i in 1:time) {Print(i);image(x, y, z[,,i], add=i>1, col=col, zlim=zlim)}

}}

\dontshow{\dontrun{
basicname <- "cyclone"
# load(paste(basicname, ".dat", sep=""))
Plot(basicname, x, y, z, col=col, T=T, pixels=256, zi=0.5 + dim(z)[2]/2,
 speed=0.2, zlim=zlim)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z)
}}


\dontshow{FinalizeExample()}


}


\keyword{spatial}

