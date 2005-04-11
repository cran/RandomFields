## check, hana

## library(RandomFields, lib="~/TMP")


if (file.exists("source.R")) source("source.R")

x <- 1:1000
x <-  c(0, 1, 0.001)
model <- "fractalB"
param <- c(0,1,0,1, 1.70)
RFparameters(Storing=FALSE, Print=5)
#plot(seq(x[1], x[2], x[3]),
#     z<-GaussRF(x, model=model, param=param, gridtriple=TRUE))
e <- EmpiricalVariogram(x, data= GaussRF(x, model=model, param=param, n=1000,
                             gridtriple=TRUE),
                        bin=(0.5 + (0 : 20)) / 40, gridtriple=TRUE)
plot(e$c, e$e - param[2] * (e$c / param[4])^param[5])



plot(e$c, e$e / ( param[2] * (e$c / param[4])^param[5]))
