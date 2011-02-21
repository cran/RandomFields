## source("geoR.R")


if (EXTENDED.TESTING <- file.exists("source.R")) {
  source("source.R")
} else if (file.exists(f <- "~/R/RF/RandomFields/tests/RFtest.R")) source(f) 

if (EXTENDED.TESTING) {
  data(soil250, package="geoR")
  soil <- soil250
  x <- soil$Linha
  y <- soil$Coluna
  z <- soil$pHAgua
  
  ##x <- 1:5
  ##y <- 1:8
  ##z <- runif(length(x) * length(y))
  ##xy <- as.matrix(expand.grid(x,y))

  ## xset b off
  ## gauss
  
  bin <- c(-1, seq(0, 20, len=40))
  (emp0 <- EmpiricalVariogram(x, y, data=z, grid=FALSE, bin=bin))
  (emp1 <- EmpiricalVariogram(unique(y), unique(x), data=z, grid=TRUE, bin=bin))

  by <- 0.15
  RFparameters(Print=1, CE.force=TRUE, CE.trials=1, CE.useprimes=TRUE)
  ShowModels(seq(min(y), max(y), by=by), seq(min(x), max(x), by=by), emp=emp1,
             me="ci", fixed.rs=TRUE)
}
