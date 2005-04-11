# source("test.mle.R")

q()

load("rs")

repeat{
#runif(11)
##
#
  save(file="rs", .Random.seed)

options(warn=4)

library(RandomFields, lib="~/TMP");
source("~/article/R/NEW.RF/RandomFields/tests/mle.old.R");
source("~/article/R/NEW.RF/RandomFields/R/MLE.R");
source("~/article/R/NEW.RF/RandomFields/R/convert.R");
source("~/article/R/NEW.RF/RandomFields/R/modelling.R");

model <-"gencauchy"
param <- c(0, 1, 0, 1, 1, 2)
estparam <- c(0, NA, 0, NA, NA, 2) ## NA means: "to be estimated"

estparam <- param ## NA means: "to be estimated"
estparam[runif(length(estparam)) < 0.5] <- NA
print(estparam)

##estparam <- param; estparam[2] <- NA

## sequence in `estparam' is
## mean, variance, nugget, scale, (+ further model parameters)
## So, mean, variance, and scale will be estimated here.
## Nugget is fixed and equals zero.
n <- 10
#n <- 1
points <- 100
points <- 15
x <- runif(points,0,3)
y <- runif(points,0,3) ## 100 random points in square [0, 3]^2
##
d <- GaussRF(x=x, y=y, grid=FALSE, model=model, param=param, n=n)
str(fitvario(Print=3,
             x=cbind(x,y), data=d, model=model, param=estparam,
             lower=c(0.1, 0.1), upper=c(1.9, 5),
             lsq.methods=c("self", "plain", "sqrt.nr", "sd.inv"),
             #cross="cross.ign"
             )[c("variogram", "values")], vec.len=6, digits.d=5)
str(try(Ofitvario(x=cbind(x,y), data=d, model=model, param=estparam,
              lower=c(0.1, 0.1), upper=c(1.9, 5))[c("mle", "mle.value")]),
    vec.len=6, digits.d=5)
print(estparam)
}


## The next two estimations give about the same result.
## For the first the sill is fixed to 1.5. For the second the sill
## is reached if the estimated variance is smaller than 1.5
estparam <-  c(0, NA, NA, NA, NA, NA) 
str(fitvario(x=cbind(x,y), data=d, model=model, param=estparam, sill=1.5))
str(Ofitvario(x=cbind(x,y), data=d, model=model, param=estparam, sill=1.5))

estparam <-  c(0, NA, NaN, NA, NA, NA) 
parampositions(model=model, param=estparam)
f <- function(param) {
   param[5] <- max(0, 1.5 - param[1])
   return(param)
}
str(fitvario(x=cbind(x,y), data=d, model=model, param=estparam, 
    sill=1, transform=f))
str(Ofitvario(x=cbind(x,y), data=d, model=model, param=estparam, 
    sill=1, transform=f))

## the next call gives a warning, since the user may programme
## strange things in this setup, and the program cannot check it.
estparam <- c(0, NA, NA, NA, NA, NaN) 
parampositions(model=model, param=estparam)
f <- function(param) {param[3] <- param[2]; param}
unix.time(str(fitvario(x=cbind(x,y), data=d, model=model,
      param=estparam, transform=f, standard.style=TRUE), vec.len=6))
unix.time(str(Ofitvario(x=cbind(x,y), data=d, model=model,
      param=estparam, transform=f, standard.style=TRUE), vec.len=6))

## much better programmed, but also much slower:
estmodel <- list(list(model="gencauchy", var=NA, scale=NA,
                      kappa=list(NA, function(m) m[[1]]$kappa[1])))
unix.time(str(fitvario(x=cbind(x,y), data=d, model=estmodel),
              vec.len=6))
unix.time(str(Ofitvario(x=cbind(x,y), data=d, model=estmodel),
              vec.len=6))

