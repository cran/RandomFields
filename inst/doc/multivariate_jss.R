### R code from vignette source 'multivariate_jss.Rnw'

###################################################
### code chunk number 1: multivariate_jss.Rnw:93-98
###################################################
library(RandomFields)
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)

options(device="jpeg")
RFoptions(resolution=60)


###################################################
### code chunk number 2: multivariate_jss.Rnw:377-378
###################################################
RFoptions(seed = 0, height = 4, always_close_screen = FALSE)


###################################################
### code chunk number 3: multivariate_jss.Rnw:381-385
###################################################
if (!interactive()) {
  #options(device="jpeg")
  RFoptions(always_close_screen = TRUE, pch="")
}


###################################################
### code chunk number 4: multivariate_jss.Rnw:483-484
###################################################
RFoptions(file="lmc",onefile=TRUE)


###################################################
### code chunk number 5: multivariate_jss.Rnw:486-493
###################################################
M1 <- c(0.9, 0.6)
M2 <- c(sqrt(0.19), 0.8)
model <- RMmatrix(M = M1, RMwhittle(nu = 0.3)) +
            RMmatrix(M = M2, RMwhittle(nu = 2))
x <- y <- seq(-10, 10, 0.3)
simu <- RFsimulate(model, x, y)
plot(simu)


###################################################
### code chunk number 6: multivariate_jss.Rnw:523-524
###################################################
RFoptions(file="delay0",onefile=TRUE)


###################################################
### code chunk number 7: multivariate_jss.Rnw:526-529
###################################################
model <- RMdelay(RMstable(alpha = 1.9, scale = 2), s = c(4, 4))
simu <- RFsimulate(model, x, y)
plot(simu, zlim = 'joint')


###################################################
### code chunk number 8: multivariate_jss.Rnw:565-567
###################################################
model <- RMdelay(RMstable(alpha = 1.9, scale = 2), s = c(0, 4)) +
            RMdelay(RMstable(alpha = 1.9, scale = 2), s = c(4, 0))


###################################################
### code chunk number 9: multivariate_jss.Rnw:569-573
###################################################
RFoptions(file="delay", onefile=!TRUE)
simu <- RFsimulate(model, x, y)
plot(simu, zlim = 'joint')
#options(device="pdf")


###################################################
### code chunk number 10: multivariate_jss.Rnw:576-578
###################################################
plot(model, dim = 2, xlim = c(-5, 5), main = "Covariance function", 
        cex = 1.5, col = "brown")


###################################################
### code chunk number 11: multivariate_jss.Rnw:611-613
###################################################
RFoptions(file="latent",onefile=TRUE)
#options(device="jpeg")


###################################################
### code chunk number 12: multivariate_jss.Rnw:615-618
###################################################
model <- RMgencauchy(alpha = 1.5, beta = 3)
simu <- RFsimulate(model, x, y, z = c(0, 1))
plot(simu, MARGIN.slices = 3, n.slices = 2)


###################################################
### code chunk number 13: multivariate_jss.Rnw:688-689
###################################################
RFoptions(file="biwm",onefile=TRUE)


###################################################
### code chunk number 14: multivariate_jss.Rnw:691-695
###################################################
model <- RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 1, cdiag = c(1, 5),
                   s = c(1, 1, 2))
simu <- RFsimulate(model, x, y)
plot(simu)


###################################################
### code chunk number 15: multivariate_jss.Rnw:727-730
###################################################
RMparswmX <- function(nudiag, rho, var, scale, Aniso, proj) {
   RMschur(M = rho, RMparswm(nudiag, var, scale, Aniso, proj))
}  


###################################################
### code chunk number 16: multivariate_jss.Rnw:737-738
###################################################
RFoptions(file="aniso",onefile=TRUE)


###################################################
### code chunk number 17: multivariate_jss.Rnw:740-746
###################################################
A1 <- RMangle(angle = pi / 4, diag = c(0.1, 0.5))
A2 <- RMangle(angle = 0, diag = c(0.1, 0.5))
model <- RMmatrix(M = M1, RMgengneiting(kappa = 0, mu = 2, Aniso = A1)) +
            RMmatrix(M = M2, RMgengneiting(kappa = 3, mu = 2, Aniso = A2))
simu <- RFsimulate(model, x, y)
plot(simu)


###################################################
### code chunk number 18: multivariate_jss.Rnw:895-897
###################################################
RFoptions(file="curlfree",onefile=!TRUE)
#options(device="pdf")


###################################################
### code chunk number 19: multivariate_jss.Rnw:899-903
###################################################
model <- RMcurlfree(RMmatern(nu = 5), scale = 4)
simu <- RFsimulate(model, x, y)
plot(simu, select.variables = list(1, 2 : 3, 4))
plot(model, dim = 2, xlim = c(-3, 3), main = "", cex = 2.3, col = "brown")


###################################################
### code chunk number 20: multivariate_jss.Rnw:953-954
###################################################
RFoptions(file="kolmogorov",onefile=!TRUE)


###################################################
### code chunk number 21: multivariate_jss.Rnw:956-962
###################################################
x <- y <- seq(-2, 2, len = 20)
model <- RMkolmogorov()
simu <- RFsimulate(model, x, y, z = 1) 
plot(simu, select.variables = list(1 : 2), col = c("red"))
plot(model, dim = 3, xlim = c(-3, 3), MARGIN = 1 : 2, cex = 2.3,
            fixed.MARGIN = 1.0, main = "", col = "brown")


###################################################
### code chunk number 22: multivariate_jss.Rnw:1008-1011
###################################################
nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())
pars.model <- nug + RMbiwm(nudiag = c(NA, NA), scale = NA, 
                              cdiag = c(NA, NA), rhored = NA)


###################################################
### code chunk number 23: multivariate_jss.Rnw:1016-1018
###################################################
data(weather)
Dist.mat <- as.vector(RFearth2dist(weather[ , 3 : 4]))


###################################################
### code chunk number 24: multivariate_jss.Rnw:1023-1024
###################################################
PT <- weather[ , 1 : 2]


###################################################
### code chunk number 25: multivariate_jss.Rnw:1026-1027 (eval = FALSE)
###################################################
## pars <- RFfit(pars.model, distances = Dist.mat, dim = 3, data = PT)


###################################################
### code chunk number 26: multivariate_jss.Rnw:1029-1030
###################################################
print(pars)


###################################################
### code chunk number 27: multivariate_jss.Rnw:1068-1069
###################################################
print(pars, full=TRUE)


###################################################
### code chunk number 28: multivariate_jss.Rnw:1090-1091
###################################################
RFratiotest(pars)


###################################################
### code chunk number 29: multivariate_jss.Rnw:1104-1105 (eval = FALSE)
###################################################
## RFcrossvalidate(pars, x = weather[ , 3 : 4], data = PT, full = TRUE)


###################################################
### code chunk number 30: multivariate_jss.Rnw:1121-1123
###################################################
whole.model <- nug + RMbiwm(nudiag = c(NA, NA), nured = NA, 
                           s = rep(NA, 3), cdiag = c(NA, NA), rhored = NA)


###################################################
### code chunk number 31: multivariate_jss.Rnw:1125-1126 (eval = FALSE)
###################################################
## whole <- RFfit(whole.model, distances = Dist.mat, dim = 3, data = PT)


###################################################
### code chunk number 32: multivariate_jss.Rnw:1130-1131
###################################################
RFratiotest(nullmodel = pars, alternative = whole)


###################################################
### code chunk number 33: multivariate_jss.Rnw:1141-1147
###################################################
a <- colMeans(weather[ , 3 : 4]) * pi / 180
P <- cbind(c(-sin(a[1]), cos(a[1]), 0),
              c(-cos(a[1]) * sin(a[2]), -sin(a[1]) * sin(a[2]), cos(a[2])),
              c(cos(a[1]) * cos(a[2]), sin(a[1]) * cos(a[2]), sin(a[2])))
x <- RFearth2cartesian(weather[ , 3 : 4])
plane <- (x %*% P)[ , 1 : 2]


###################################################
### code chunk number 34: multivariate_jss.Rnw:1152-1154
###################################################
RFoptions(file="kriging", onefile=TRUE)
#options(device="jpeg")


###################################################
### code chunk number 35: multivariate_jss.Rnw:1156-1164
###################################################
n <- 200
r <- apply(plane, 2, range)
data <- cbind(plane, weather[ , 1 : 2])
z <- RFinterpolate(pars, data = data, dim = 2,
                      x = seq(r[1, 1], r[2, 1], length = n),
                      y = seq(r[1, 2], r[2, 2], length = n),
                      varunits = c("Pa", "K"))
plot(z)


###################################################
### code chunk number 36: multivariate_jss.Rnw:1240-1248
###################################################
# Sweave does not work in BATCH mode without the following lines.
# Sweave memorizes somehow that the above plot's have first closed
# Sweave's graphical output before opening one's own output. As a
# revanche, Sweave kills the last figure.
if (!interactive()) {
  RFoptions(file="schnulli_schnulli",onefile=TRUE)
  plot(simu)
}


