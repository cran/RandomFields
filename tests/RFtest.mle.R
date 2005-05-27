# source("RFtest.MLE.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

note <- readline
note <- print
max.le <- Inf
max.le <- 30
max.le <- 5

n <- 30
model <- list(list(model="exp", var=5, aniso=c(1,0,0,2)), "+",
              list(model="whittle", var=3, kappa=1.5, aniso=c(1/4,0,0,1/2)),
              mean=2)
x <- 1:3
RFparameters(Print=2)

# runif(1); save(file="rs", .Random.seed)
# if (file.exists("rs")) load("rs")

z <- GaussRF(x, x, grid=TRUE, model=model, n=n)

EST <- PrepareModel(model, timespace=2)
le <- min(max.le, 2^(length(EST$param)) - 1)
mle <- NULL
for (i in 1:le) {
  for (me in FALSE:TRUE) {
    est <- PrepareModel(model, timespace=2)
    if (me) est$mean <- NA
    ii <- i
    for (j in 1:length(est$param)) {
      if (ii %% 2 ==1) est$param[j] <- NA
      ii <- as.integer(ii / 2)
    }
  if (i < le || !me) next ###################################
    est <- convert.to.readable(est)
    cat("\n\n i=", i, "m=", me, "\n")
    str(est)
    str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, nphi=2,
                       Print=1, cross.m=NULL)$vario$ml)
    FV <- PrepareModel(fv, timespace=2)
    mle <- rbind(mle, FV$param)
    if (any(abs(FV$param - EST$param) > 0.05)) {
      print(rbind(EST$param, FV$param))
      note("WARNING -- deviations detected! Press return")
    }
  }
}


######################################################################
######################################################################

x <- 1:3
coord <- cbind(c(1, 3, 6, 10),
           c(0.2, 0.5, 0.9, 0))
coord <- expand.grid(x,x)
lc <- nrow(coord)
distances <- matrix(nrow=ncol(coord), ncol=lc * (lc-1) / 2)
dd <- 0
for (i in 1:ncol(coord)) {
  dx <- outer(coord[,i], coord[, i], "-")
  print(dx)
  distances[i, ] <- as.double(dx[lower.tri(dx)])
  ## x <- c(1,6,0.1)
  ## outer(x, x, "-
  ##     [,1] [,2] [,3]
  ##[1,]  0.0 -5.0  0.9
  ##[2,]  5.0  0.0  5.9
  ##[3,] -0.9 -5.9  0.0
  dd <- dd + distances[i ,]^2
}
DD <- as.double(distances)

y <- t(matrix(DD, nrow=2))
len <- (1 + sqrt(1 + 8 * nrow(y))) / 2 
m1 <- model
m1[[1]]$var <- 5.19
cov <- CovarianceFct(y, model=m1)
cm <- diag( CovarianceFct(cbind(0,0), model=m1)/2, ncol=len, nrow=len)
cm[lower.tri(cm)] <- cov
cm <- cm + t(cm)
eigen(cm)$val

CovarianceFct(y, model=m1, fct="CovarianceMatrix")

##################################################################
## dim = 2, old definitions
##################################################################

n <- 3000
n <- 30
model <- list(model="exp", param=c(4, 3 , 1, 5))
x <- 1:10
if (!exists("z") || any(dim(z) != c(rep(length(x), 2), n)))
  z <- GaussRF(x, x, grid=TRUE, model=model, n=n)
est <- PrepareModel(model, timespace=2)
if (TRUE) {
  est$param[1 + runif(2, max=length(est$param))] <- NA
  if (runif(1) < 0.5) est$mean <- NA
  ## est$param[1] <- NA
} else {
  est <- PrepareModel(model, timespace=2)
  est$param[c(1,2,4)] <- NA
  est$mean <- NA
}
est <- convert.to.readable(est, allowed="standard")

str(est)
str(xxx <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1,
                                       cross.m=NULL))
print(xxx$value$ml)
fitvario(x=expand.grid(x,x), data=z, model=xxx$vario$ml, Print=1,
                                       cross.m=NULL)

le <- min(max.le, 2^(length(model$param)) - 1)
mle <- NULL
for (i in 1:le) {
  est <- model
  ii <- i
  for (j in 1:length(est$param)) {
    if (ii %% 2 ==1) est$param[j] <- NA
    ii <- as.integer(ii / 2)
  }
  cat("\n\n i=", i, "\n")
  str(est)
  str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1,
                     cross.m=NULL))
  mle <- rbind(mle, fv$vario$ml$param)
  if (any(abs(fv$vario$ml$param-model$param) > 0.05))
    note("WARNING -- deviations detected! Press return")
  fv2 <- fitvario(x=expand.grid(x,x), data=z, model=fv$vario$ml,
                                       cross.m=NULL)$value$ml
  if ( abs(fv$value$ml - fv2) > abs(fv2) * 10^(-14)){
    print(fv$value$ml)
    print(fv2)
    print(fv$value$ml - fv2)
    stop("mle.values unequal")
  }
}


######################################################################
if (EXTENDED.TESTING) {
######################################################################

n <- 3000
model <- list(model="whittle", param=c(4, 3 , 1, 5, 1.5))
x <- 1:10
#if (!exists("z") || any(dim(z) != c(rep(length(x), 2), n)))
z <- GaussRF(x, x, grid=TRUE, model=model, n=n)

le <- min(max.le, 2^(length(model$param)) - 1) ## -1 wegen scale von nugget
mle <- NULL
for (i in 1:le) {
  est <- model
  ii <- i
  for (j in 1:length(est$param)) {
    if (ii %% 2 ==1) est$param[j] <- NA
    ii <- as.integer(ii / 2)
  }
  cat("\n\n i=", i,  "\n")
  str(est)
  str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1,
                                       cross.m=NULL)$vario)
  mle <- rbind(mle, fv$ml$param)
  if (any(abs(fv$ml$param-model$param) > 0.05))
    note("WARNING -- deviations detected! Press return")
}





##################################################################
## isotropic, new definition
##################################################################

n <- 3000
model <- list(list(model="exp", var=5, sc=1), "+",
              list(model="whittle", var=3, sc=4, kappa=2),
              mean=2)
x <- 1:10
z <- GaussRF(x, x, grid=TRUE, model=model, n=n)

EST <- PrepareModel(model, timespace=2)
le <- min(max.le, 2^(length(EST$param)) - 1)
mle <- NULL
for (i in 1:le) {
  for (me in FALSE:TRUE) {
    est <- PrepareModel(model, timespace=2)
    if (me) est$mean <- NA
    ii <- i
    for (j in 1:length(est$param)) {
      if (ii %% 2 ==1) est$param[j] <- NA
      ii <- as.integer(ii / 2)
    }
    est <- convert.to.readable(est)
    cat("\n\n i=", i, "m=", me, "\n")
    str(est)
    str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1,
                                       cross.m=NULL)$vario$ml)
    FV <- PrepareModel(fv, timespace=2)
    mle <- rbind(mle, FV$param)
    if (any(abs(FV$param - EST$param) > 0.05)) {
      print(rbind(EST$param, FV$param))
      note("WARNING -- deviations detected! Press return")
    }
  }
}

## funktioniert nicht : # 23, var=NA, scale=NA; var=NA, scale=NA
##                        27, var=NA, scale=NA; kappa=NA, scale=NA
##                        28                  ; NA, NA, NA
##                        31  NA, NA; NA, NA, NA
### das scheint aber eher daran zu liegen, dass das Problem schlecht gestellt
### ist


##################################################################
## anisotropic
##################################################################
n <- 3000
model <- list(list(model="exp", var=5, aniso=c(1,0,0,2)), "+",
              list(model="whittle", var=3, kappa=1.5, aniso=c(1/4,0,0,1/2)),
              mean=2)
x <- 1:10
z <- GaussRF(x, x, grid=TRUE, model=model, n=n)
EST <- PrepareModel(model, timespace=2)
le <- min(max.le, 2^(length(EST$param)) - 1)
mle <- NULL
for (i in 1:le) {
  for (me in FALSE:TRUE) {
    est <- PrepareModel(model, timespace=2)
    if (me) est$mean <- NA
    ii <- i
    for (j in 1:length(est$param)) {
      if (ii %% 2 ==1) est$param[j] <- NA
      ii <- as.integer(ii / 2)
    }
    est <- convert.to.readable(est)
    cat("\n\n i=", i, "m=", me, "\n")
    str(est)
    fv <- fitvario(x=expand.grid(x,x), data=z, model=est, nphi=2,
                       Print=1, cross.m=NULL)
    FV <- PrepareModel(fv$vario$ml, timespace=2)
    mle <- rbind(mle, FV$param)
    if (any(abs(FV$param - EST$param) > 0.05)) {
      print(rbind(EST$param, FV$param))
      note("WARNING -- deviations detected! Press return")
    }
  }
}

## i=8
######################################################################
} ## end: if (EXTENDED.TESTING)
######################################################################
