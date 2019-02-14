\name{RFoptionsAdvanced}
\alias{RFoptionsAdvanced}
\title{Setting control arguments of \pkg{RandomFields} -- advanced examples}
\description{
  Some more complex examples for the use of \command{\link{RFoptions}}
  are given.
}

\me

\examples{\dontshow{StartExample()}
%% NOTE: THE BELOW EXAMPLE NEEDS seed = NA
 #############################################################
 ##                      EXAMPLE 1                          ##
 ## The following gives an example on the advantage of      ##
 ## dependent=TRUE for simulating with RPcirculant if, in a ##
 ## study, most of the time is spent with simulating the    ##
 ## Gaussian random fields. Here, the covariance at a pair  ##
 ## of points is estimated for n independent repetitions    ##
 ## and 2*n locally dependent repetitions.                  ##
 ## To get the precision, the procedure is repeated m times.##
 #############################################################

# In the example below, local.dependent speeds up the simulation
# by about factor 16 at the price of an increased variance of
# factor 1.5

RFoptions(seed=NA)
len <- 10
\dontshow{if(!interactive()){warning("reduced 'len'");len<-2;}}%ok
x <- seq(0, 1, len=len)
y <- seq(0, 1, len=len)
grid.size <- c(length(x), length(y))
meth <- RPcirculant
model <- RMexp(var=1.1, Aniso=matrix(nc=2, c(2,0.1,1.5,1)))

m <- 5 
n <- 100
\dontshow{if(!interactive()){warning("reduced 'm'");m<-2;}}%ok

# using local.dependent=FALSE (which is the default)
c1 <- numeric(m)
time <- system.time(
  for (i in 1:m) {
    cat("", i, "out of", m, "\n")
    z <- RFsimulate(meth(model), x, y, n=n, pch="", 
                    dependent=FALSE, spConform=FALSE, trials=5, force=TRUE)
    c1[i] <- cov(z[1, dim(z)[2], ], z[dim(z)[1], 1, ])
}) # many times slower than with local.dependent=TRUE below
\dontshow{cat("\n")}
true.cov <- RFcov(model, t(y[c(1, length(y))]), t(x[c(length(x), 1)]))
print(time)
Print(true.cov, mean(c1), sd(c1), empty.lines=1)## true mean is zero

# using local.dependent=TRUE ...
c2 <- numeric(m)
time <- system.time(
  for (i in 1:m) {
    cat("", i)
    z <- RFsimulate(meth(model), x, y, n=2 * n, pch="", 
                    dependent=TRUE, spConform=FALSE, trials=5, force=TRUE)
    c2[i] <- cov(z[1, dim(z)[2], ], z[dim(z)[1], 1, ])
})
\dontshow{cat("\n")}
print(time)                                      ## 20 times faster
Print(true.cov, mean(c2), sd(c2), empty.lines=1) ## much better results

## the sd is smaller (using more locally dependent realisations)
## but it is (much) faster! Note that for n=n2 instead of n=2 * n, 
## the value of sd(c2) would be larger due to the local dependencies 
## in the realisations.




 #############################################################
 ##                      EXAMPLE 2                          ##
 ## This example shows that the same realisation can be     ##
 ## obtained on different grid geometries (or point         ##
 ## configurations, i.e. grid, non-grid) using TBM          ##
 #############################################################

RFoptions(seed=0)
step <- 1
x1 <- seq(-150,150,step)
y1 <- seq(-15, 15, step)
x2 <- seq(-50, 50, step)
model <- RPtbm(RMexp(scale=10))
\dontshow{if(!interactive()){warning("values for x1, y1, x2 changed'");x1<-y1<-x2 <- -1:1}}%ok

RFoptions(storing=TRUE)
mar <- c(2.2, 2.2, 0.1, 0.1)
points <- 700

###### simulation of a random field on long thin stripe
z1 <- RFsimulate(model, x1, y1, center=0, seed=0,
                 points=points, storing=TRUE, spConform=FALSE)
ScreenDevice(height=1.55, width=12)
par(mar=mar)
image(x1, y1, z1, col=rainbow(100))
polygon(range(x2)[c(1,2,2,1)], range(y1)[c(1,1,2,2)],
        border="red", lwd=3)


###### definition of a random field on a square of shorter diagonal
z2 <- RFsimulate(model, x2, x2, register=1, seed=0,
                 center=0, points=points, spConform=FALSE)
ScreenDevice(height=4.3, width=4.3)
par(mar=mar)
image(x2, x2, z2, zlim=range(z1), col=rainbow(100))
polygon(range(x2)[c(1,2,2,1)], range(y1)[c(1,1,2,2)],
        border="red", lwd=3)
tbm.points <- RFgetModelInfo(level=3)$loc$totpts
Print(tbm.points, empty.lines=0) # number of points on the line


\dontshow{FinalizeExample()}}
