# source("RFtest.MLE.R")

if (file.exists("source.R")) source("source.R")

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
    str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, nphi=2,
                       Print=1)$mle)
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
str(xxx <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1))
print(xxx$mle.value)
fitvario(x=expand.grid(x,x), data=z, model=xxx$mle, Print=1)$mle.value

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
  str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1))
  mle <- rbind(mle, fv$param)
  if (any(abs(fv$mle$param-model$param) > 0.05))
    note("WARNING -- deviations detected! Press return")
  fv2 <- fitvario(x=expand.grid(x,x), data=z, model=fv$mle)$mle.value
  if ( abs(fv$mle.value - fv2) > abs(fv2) * 10^(-14)){
    print(fv$mle.value)
    print(fv2)
    print(fv$mle.value - fv2)
    stop("mle.values unequal")
  }
}


######################################################################
######################################################################

if (file.exists("source.R")) source("source.R")
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
  str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1)$mle)
  mle <- rbind(mle, fv$mle$param)
  if (any(abs(fv$mle$param-model$param) > 0.05))
    note("WARNING -- deviations detected! Press return")
}





##################################################################
## isotropic, new definition
##################################################################

if (file.exists("source.R")) source("source.R")
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
    str(fv <- fitvario(x=expand.grid(x,x), data=z, model=est, Print=1)$mle)
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
if (file.exists("source.R")) source("source.R")
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
    Xfv <- fitvario(x=expand.grid(x,x), data=z, model=est, nphi=2,
                       Print=3)
    str(fv <- Xfv[c("mle", "mle.lower", "mle.upper")])
    FV <- PrepareModel(fv$mle, timespace=2)
    mle <- rbind(mle, FV$param)
    if (any(abs(FV$param - EST$param) > 0.05)) {
      print(rbind(EST$param, FV$param))
      note("WARNING -- deviations detected! Press return")
    }
  }
}

## i=8
