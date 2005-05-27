# source("RFtest.info.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

str(GetRegisterInfo(0))

cat("-------------------------------------------\n")
try(GaussRF(1:4, grid=TRUE, model="exp", param=c(1,2,3,4)))
str(GetRegisterInfo(0))

try(GaussRF(runif(4), grid=FALSE, model="exp", param=c(1,2,3,4), me="ci"))
str(GetRegisterInfo(0, TRUE))

