\name{S10}
\alias{S10}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{On some covariance models based on normal scale mixtures}
\description{
 Here, the code of the paper on \sQuote{On some covariance models based on normal scale mixtures} is given.
}
 
\me
\references{
 \itemize{ 
 \item
 Schlather, M. (2010)
 On some covariance models based on normal scale mixtures.
 \emph{Bernoulli}, \bold{16}, 780-797.
 }

}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

### Example 10 in Schlather (2010).
## The field below has more than 80 million points. So the simulation
## takes a while
y <- x <- seq(0, 10, len=256) ## currently does not work
T <- c(0, 0.02, 1275)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2,
         101:150)])
%\dontshow{if(RFoptions()$internal$examples_red){warning("reduced x,y,T")
y <- x <- seq(0, 10, len=5)
T <- c(0, 0.02, 4)
%}}
model <- RMcoxisham(mu=c(1, 1), D=matrix(nr=2, c(1, 0.5, 0.5, 1)),
                    RMwhittle(nu=1))
z <- RFsimulate(model, x, y, T=T, sp_lines=1500, every=10)
plot(z, MARGIN.slices=3, col=col)
plot(z, MARGIN.movie=3) # add 'file="ci.avi"' to get it stored


\dontshow{\dontrun{
### Example 13 in Schlather (2010)
y <- x <- seq(0, 10,len = 256)
T <- c(0, 0.02, 1275)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)]) 
\dontshow{if(RFoptions()$internal$examples_red){warning("reduced x,y,T"); 
y <- x <- seq(0, 10,len = 10) 
T <- c(0, 0.02, 4)
}}
model <- RMave(A = matrix(nc=2,c(0.5, 0, 0, 1)), z= c(2,0),
               RMwhittle(nu=1))
z <- RFsimulate(model, x, y, T=T, every = 10,
                trials=2, force=TRUE, maxmem=16777216*8)
plot(z)
}}

\dontshow{\dontrun{
basicname <- "moving"
Plot(basicname, x, y, z, col=col, T=T)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z) # OK
}}


\dontshow{\dontrun{

### Example 16 in Schlather (2010)
intens <- 1000 ## takes a huge amount of time; take smaller values
len <- 81
\dontshow{if(RFoptions()$internal$examples_red){warning("modified intens & len");intens<-3;len<-4}}
y <- x <- seq(-3, 3, len=len)
T <- seq(0, 0.075, len=len)
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2,
101:150)])
\dontshow{if(RFoptions()$internal$examples_red){warning("reduced 'intens'"); intens<-3}}
model <- RMstp(M=matrix(nc=3, rep(0, 9)),
               S=RMetaxxa(E=c(1, 1, 1), alpha = -2 * pi,
                          A=t(matrix(nc=3, c(2, 0, 0, 1, 1 , 0, 0, 0, 0)))),
               Aniso = RMrotation(phi= -2 * pi),
               phi = RMwhittle(nu = 1) )
z <- RFsimulate(model, x, x, z=T, me="coin", every = 10, 
 mpp.intens=intens, mpp.p = 0.1, mpp.beta=3.5, mpp.plus = 8, print=3) 
zlim <- c(-3.5, 3.5)
time <- dim(z)[3]
for (i in 1:time) {Print(i);image(x, y, z[,,i], add=i>1, col=col, zlim=zlim)}

}}

\dontshow{\dontrun{
basicname <- "cyclone"
Plot(basicname, x, y, z, col=col, T=T, pixels=256, zi=0.5 + dim(z)[2]/2,
 speed=0.2, zlim=zlim)
save(file=paste(basicname, "/", basicname, ".dat", sep=""), z) # OK
}}

\dontshow{FinalizeExample()}}


\keyword{spatial}

