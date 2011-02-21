


library(RandomFields, lib="~/TMP")
#
source("RandomFields/tests/source.R")



my.arrows <- function(xy, z, r, thinning, col=1) {
  startx <- as.vector(xy[,1] - r/2*z[, 1])
  starty <- as.vector(xy[,2] - r/2*z[ ,2])
  endx <- as.vector(xy[,1] + r/2*z[, 1])
  endy <- as.vector(xy[,2] + r/2*z[, 2])

#  startx[c(rep(TRUE, thinning), FALSE)] <- NA
#  starty[c(rep(TRUE, thinning), FALSE)] <- NA
#  endx[c(rep(TRUE, thinning), FALSE)] <- NA
#  endy[c(rep(TRUE, thinning), FALSE)] <- NA
#  Print(startx, starty, endx, endy)
  arrows(x0=startx, y0=starty, x1=endx, y1=endy, length=0.03,col=col)
}



x <- c(-3, 3, 0.2) / 10
T = c(0, 15, 0.1)

#
x <- c(-3, 3, 3) 
T = c(0, 15, 10)  

xx <- seq(x[1], x[2], x[3])
tt <- seq(T[1], T[2], T[3])
nu <- 3
col <- grey(seq(0, 1, 0.01)) 
thinning <- 0
length.arrow <- 1.5 / max(thinning, 10)
runif(1)
seed <- .Random.seed
eps <- interactive()  # true falls eps/pdf drucken

a <- 1
M <- rbind(diag(2), 0)
model <- list("vector", a=a,
              list("coxisham", mu=c(1,1), rho=0.5, beta=1,
                   list("whittle", nu=nu)))

r <- c(1,1,1)
model <- list("delay", r=r, list("exp")  )

xy <- as.matrix(expand.grid(xx, xx))

# .Random.seed <- seed
Z <- GaussRF(x, x, T=T,
             grid=TRUE, gridtr=TRUE, model=model, n=1, CE.trial=1,
             Stor=TRUE, me="ci", Print=2, CE.force=TRUE)

tt <- 0
z <- Z[,,,1]
zz <- t(matrix(z, nrow=2, ncol=prod(dim(z))/2))
xyT <- as.matrix(expand.grid(xx, xx, tt))
idx <- sample(nrow(xyT), 5)
estimx <- xyT[idx, ,drop=FALSE]
datax <- zz[idx, ,drop=FALSE]
givenx <- xyT[-idx, ,drop=FALSE]
givendata <- zz[-idx, ,drop=FALSE]

givenx <- xyT
givendata <- zz

CM <- function(xy, model, r) {
  nx <- nrow(xy)
  res <- matrix(nr=nx, nc=nx)
  for (i in 1:nx) {
    res[, i] <- Covariance(t(t(xy)-xy[i,]+r), model=model)
  }
  return(res)
}

options(digits=3)
(C <- rbind(cbind(CM(xyT, model=list("exp"), 0), CM(xyT, model=list("exp"), -r)),
           cbind(CM(xyT, model=list("exp"), r), CM(xyT, model=list("exp"), 0))))
invcov <- solve(C, as.double(zz))

Kx <-  function(xy, x, model, r) {
   cbind( c(Covariance(t(t(xy)-x), model=model),
            Covariance(t(t(xy)-x+r), model=model)),
         c(Covariance(t(t(xy)-x-r), model=model),
           Covariance(t(t(xy)-x), model=model))
         )
}

for (j in 1:1) { # {length(idx)) {
  i <- idx[j]
  Print(invcov, Kx(xyT, x=estimx[1,], model=list("exp"), r),
        invcov %*% Kx(xyT, x=estimx[j,], model=list("exp"), r),
        zz[i,], xyT[i,], estimx[j,])
}

#source("RandomFields/tests/source.R")
SK <- Kriging("S", x=estimx, grid=FALSE, model=model, given=givenx,
              data=givendata)
Print( SK, datax, SK-datax) ## muss NULL rauskommen

xxxxxxxx




my.arrows <- function(xy, z, r, thinning, time, col=1) {
  startx <- as.vector(xy[,1] - r/2*z[1, ,, time])
  starty <- as.vector(xy[,2] - r/2*z[2,,, time])
  endx <- as.vector(xy[,1] + r/2*z[1, ,, time])
  endy <- as.vector(xy[,2] + r/2*z[2, ,, time])
  startx[c(rep(TRUE, thinning), FALSE)] <- NA
  starty[c(rep(TRUE, thinning), FALSE)] <- NA
  endx[c(rep(TRUE, thinning), FALSE)] <- NA
  endy[c(rep(TRUE, thinning), FALSE)] <- NA
  arrows(x0=startx, y0=starty, x1=endx, y1=endy, length=0.03,col=col)
}


x <- c(-3, 3, 0.2)
xx <- seq(x[1], x[2], x[3])
T = c(0, 15, 0.1)

tt <- seq(T[1], T[2], T[3])
nu <- 3
col <- grey(seq(0, 1, 0.01))
thinning <- 0
length.arrow <- 1.5 / max(thinning, 10)
runif(1)
seed <- .Random.seed
eps <- interactive()  # true falls eps/pdf drucken

modelname <- "vector"
a <- 1

model <-
  list(modelname, a=a,
       list("coxisham", mu=c(1,1), D=matrix(c(1,0.5,0.5,1), nr=2),
            beta=2,# in(0,2]
            list("whittle", nu=nu)))

RFparameters(Print=5)
xy <- as.matrix(expand.grid(xx, xx))

.Random.seed <- seed
z <- GaussRF(x, x, T=T,
             grid=TRUE, gridtr=TRUE, model=model, n=1, CE.trial=1,
             Stor=TRUE, me="ci", Print=2, CE.force=TRUE)

par(mar=c(2.2,2.5,0.5,0.5), cex.axis=2, bg="white")



jpg <- jpeg
Plot <- function(basicname, x, y, z, pixels=length(x), col=col, delay=1,
                 basicn=10000, T, redo=TRUE, height=5, zi=1, speed=0.3,
                 zlim = range(z)) {
  Print(try(dir.create(basicname)))
  if (length(T)==3) T <- seq(T[1], T[2], T[3])
  cex = 2
  par(mar=c(4,4.3, 0.8, 0.8), cex=0.7)

  k <- 0
  time <- length(T)
  if (redo) {
    system(paste("rm ", basicname, "/", basicname, "*.jpg", sep=""))

    for (i in 1:(time)) {
      for (j in 1:delay) {
        k <- k + 1
        name <- paste(basicname, "/", basicname, k + basicn, sep="")
        Dev(TRUE, "jpg", ps=name, quiet=FALSE,
            height=pixels * height, width=pixels * height,
            type="Xlib")
        par(mar=c(0,0,0,0))

        plot(Inf, Inf, xlim=range(xx), ylim=range(xx))
        my.arrows(xy, z, length.arrow, thinning, i,
                  col=rainbow(length(tt))[time])
        Dev(FALSE)
      }
    }
  }
  system(paste("cd ", basicname, ";mencoder -mf fps=30 -ffourcc DX50
-ovc lavc ",
               " -speed ", speed,
               " -lavcopts
vcodec=mpeg4:vbitrate=9800:aspect=4/3:vhq:keyint=15",
               " -vf scale=720:576 -o ", basicname, ".avi mf://",
               basicname, "*.jpg &", sep=""))
}


# Plot("vectors", x, x, z, col=col, T=T, height=500)


for (i in 1:length(tt)) {
  plot(Inf, Inf, xlim=range(xx), ylim=range(xx))
  my.arrows(xy, z, length.arrow , thinning, i)
  readline()
}



xxxxxxxxxxxxxxxxxxxxxxxxxxxx




x <- c(-3, 3, 0.1)
xx <- seq(x[1], x[2], x[3])
Print(length(xx))

r <- c(1, 1)

model <- list("delay", r
              =r, list("exp"))
xy <- as.matrix(expand.grid(xx, xx, 0))

M <- CovMatrix(xy, model=model)
Print(dim(M))

Kx <- Covariance(t(t(xy) - xy[1,]), model=model)
M[,1] - Kx[,1]  ## should vanish 
M[,10] - Kx[,2] ## should vanish ?

solve(M, Kx) ## shoud give unit vector


d <- as.matrix(dist(xy))

dd <- NULL
rr <- r
for (i in 1:nrow(xy)) {
  dd <- cbind(dd, colSums((t(xy) - xy[i, ] + rr)^2) )
}
dd <- sqrt(dd)

n <- length(xx)^2 
print(M[1:n, 1:n] - exp(-d))  ## should vanish
print(M[n+(1:n), 1:n  ] - exp(-dd))  ## should vanish




xxxxxxxxxxxxxxxxxxxxxxxxxxxx



my.arrows <- function(xy, z, r, thinning, time, col=1) {
  startx <- as.vector(xy[,1] - r/2*z[, 1 ])
  starty <- as.vector(xy[,2] - r/2*z[, 2])
  endx <- as.vector(xy[,1] + r/2*z[, 1 ])
  endy <- as.vector(xy[,2] + r/2*z[, 2 ])
  startx[c(rep(TRUE, thinning), FALSE)] <- NA
  starty[c(rep(TRUE, thinning), FALSE)] <- NA
  endx[c(rep(TRUE, thinning), FALSE)] <- NA
  endy[c(rep(TRUE, thinning), FALSE)] <- NA
  arrows(x0=startx, y0=starty, x1=endx, y1=endy, length=0.03,col=col)
}


x <- c(-3, 3, 1) / 5
xx <- seq(x[1], x[2], x[3])
Print(length(xx))

T = c(0, 15, 15)

r <- c(1, 1, 1)

model <- list("delay", r=r, list("exp"))
xy <- as.matrix(expand.grid(xx, xx, 0))


Z <- GaussRF(x, x, T=T, grid=TRUE, gridtr=TRUE, model=model, n=1,
CE.trial=1, Stor=TRUE, me="ci", Print=2, CE.force=TRUE)

z <- Z[,,,1, drop=FALSE]
tt <- 0

zz <- t(matrix(z, nrow=2, ncol=prod(dim(z))/2))
xyT <- as.matrix(expand.grid(xx, xx, tt))

idx <- sample(nrow(xyT), 5)
estimx <- xyT[idx ,,drop=FALSE]

datax <- zz[idx, ,drop=FALSE]
givenx <- xyT[-idx ,,drop=FALSE]
givendata <- zz[-idx, ,drop=FALSE]

alen <- diff(range(xx)) / 20
plot(Inf, Inf, xlim=range(xx), ylim=range(xx))

estimx <- estimx[5:1,]
datax <- datax[5:1,]
my.arrows(givenx, givendata, alen, 0, col="black")
my.arrows(estimx, datax, alen, 0, col="green")

source("RandomFields/tests/source.R")
OK <- Kriging("S", x=estimx, grid=FALSE, model=model, given=givenx,
data=givendata)

my.arrows(estimx, OK, alen, 0, col="red")




xxxxxxxxxxxxxxxxxxxx
