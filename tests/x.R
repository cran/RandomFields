# source("x.R")


if (file.exists("source.R")) source("source.R") 


RFparameters(Print=6)

x <- seq(0, 1, 1/24)
model <- list(
              list(model="stable", var=1, scale=1, kappa=2),
              "+",
              list(model="cauchy", var=1, scale=1, kappa=2)
              )
try(z <- GaussRF(x, x, grid=TRUE, gridtriple=FALSE,
        model=model, n=1, exactness=NA, stationary.only=NA))
