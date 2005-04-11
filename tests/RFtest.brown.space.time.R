# source("RFtest.brown.space.time.R")

###############################################################################
##   Brownian motion, test
###############################################################################

if (file.exists("source.R")) source("source.R")

#if (file.exists("random")) load("random")

model <- "fractalB"
## 3d
RFparameters(pch="")
kappa <- 2 * runif(1)
scale <- 5 * runif(1)
var   <- 5 * runif(1)

n <- 1000; TBM2.lines<-60
n <- 2; TBM2.lines<-2
#kappa <- scale <- var <-  1

x <- c(1,5,1)#x <- c(1,4,1)
z <- GaussRF(x=x, y=x, z=x, grid=TRUE, gridtriple=TRUE, model=model,
		me="intr", param=c(0,var,0,scale,kappa), n=5 * n)

cat(Variogram(x=cbind(1,1,1), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1,1, ]- z[2,2,2, ])/2,"\n")

cat(Variogram(x=cbind(1,1,0), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1,1, ]- z[2,2,1,])/2,"\n")

cat(Variogram(x=cbind(0,0,1), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1,1, ]- z[1,1,2,])/2,"\n")

cat(Variogram(x=cbind(1,3,4), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1,1, ]- z[2,4,5,])/2,"\n")

# stop("")
print("A")


## 2d
RFparameters(pch="")
kappa <- 2 * runif(1)
scale <- 5 * runif(1)
var   <- 5 * runif(1)
#scale <- kappa <- var <- 1

x <- c(1,5,1)
z <- GaussRF(x=x, y=x, grid=TRUE, gridtriple=TRUE, model=model,
		me="intr", param=c(0,var,0,scale,kappa) ,n=10 * n)

cat(Variogram(x=cbind(1,1), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1, ]- z[2,2,])/2,"\n")

cat(Variogram(x=cbind(1,0), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1, ]- z[2,1,])/2,"\n")

cat(Variogram(x=cbind(0,1), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1, ]- z[1,2,])/2,"\n")

cat(Variogram(x=cbind(3,4), model=model, param=c(0,var,0,scale,kappa)),
   "",var(z[1,1, ]- z[4,5,])/2,"\n")

# stop("")


print("B")


DeleteRegister(0)
RFparameters(Print=7, TBMCE.trials=4, TBMCE.force=TRUE, pch="")
RFparameters(Print=4, TBMCE.force=TRUE, pch="#",
             CE.mmin=8, TBMCE.strategy=1, TBMCE.trials=8,
             TBMCE.mmin=c(-1,-3),
             TBM2.linesimufactor=50,
             TBM2.lines=TBM2.lines
             )

######################################################################
## space-time
######################################################################


x <- c(1,2,1)
y <- c(4,5,1)
T <- c(1,2,1)

runif(1)
rs <- .Random.seed
#if (file.exists("random")) load("random") else
#{ runif(1); save(file="random",.Random.seed) }
var <- runif(1)
aniso <- cbind(matrix(runif(6),ncol=2), c(0,0,runif(1)))
#aniso <- diag(3) * 0.5
#aniso <- diag(2) * 10
kappa <- c(0.5, 2, runif(1,0,2), runif(1,0,1),
           1 + 00000 * floor(runif(1,0,3)), 2 + rpois(1,2))
var <- 1 ###
model <- list(list(model="nsst", var=var, aniso=aniso, k=kappa))
#model <- list(list(model="sph", var=var, aniso=diag(c(runif(2), 0))),"*",
#            list(model="exp", var=var, aniso=diag(c(0, 0, runif(1))))
#              )
.Random.seed <- rs
print(model)
#x <- c(1,4,1)

print("C")

## stop("")
print(CovarianceFct(x=cbind(1,1,1), model=model))
print(Variogram(x=cbind(1,1,1), model=model))

z <- GaussRF(x=x, y=y, T=T, grid=TRUE, gridtriple=TRUE, model=model, me="TBM2",
             n=n)

print(CovarianceFct(x=cbind(1,1,1), model=model))
print(cov(z[1,1,1, ], z[2,2,2,]))
print("")

print(CovarianceFct(x=cbind(1,0,0), model=model))
print(cov(z[1,1,1, ], z[2,1,1,]))
print("")

print(CovarianceFct(x=cbind(0,0,1), model=model))
print(cov(z[1,1,1, ], z[1,1,2,]))


print("D")

# stop("")

## does spherical work fine?
DeleteRegister(0)
RFparameters(Print=2, TBMCE.force=TRUE, pch="#",
             TBMCE.strategy=1, TBMCE.trials=8,
             TBMCE.mmin=c(-1,-3),
             TBM2.linesimufactor=50,
             CE.mmin=c(-1,-3),
             CE.userfft=TRUE
             )
aniso <- matrix(c(0.01,0,0,0.5),ncol=2)
model <- list(list(model="sph", var=var, aniso=aniso))
z <- GaussRF(x=c(1,75,1), y=y, grid=TRUE, gridtriple=TRUE, model=model, me="ci"
             ,n=3 * n)

print(CovarianceFct(x=cbind(1,1), model=model))
print(cov(z[1,1, ], z[2,2,]))
print("")

print(CovarianceFct(x=cbind(1,0), model=model))
print(cov(z[1,1, ], z[2,1,]))
print("")

print(CovarianceFct(x=cbind(0,1), model=model))
print(cov(z[1,1, ], z[1,2,]))

print(var(z))

print("D")

# stop("")

## this is used in combination with variable xx in RFcircembed.cc
## are all the elements of the FFT-matrix reached?
repeat {
  nn <- 4 ## 1..4
  x <- rbind(1,as.integer(runif(nn,2,2+rpois(1,100)^(1/nn))),1)
  #x <- rbind(1,as.integer(runif(nn,2,2+rpois(1,10000)^(1/nn))),1)
  print(x)
  z <- GaussRF(x=x, grid=TRUE, gridtriple=TRUE, model="exp",
               me="ci", param=c(0,1,0,0.0001) ,n=1)
  cat("OK\n\n")
  break; ##
}


# stop("")
