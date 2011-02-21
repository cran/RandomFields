
##  source("RFtest.Precision.R")

if (EXTENDED.TESTING <- file.exists("source.R")) {
  source("source.R")
  source("RFtest.R")
} else if (file.exists(f <- "~/R/RF/RandomFields/tests/RFtest.R")) source(f) 


#RFparameters(Print=2, direct.method=1)
#GaussRF(x=1:3, model="exp", param=c(0,1,0,1), me="di")
#xxxxx

if (EXTENDED.TESTING) {

load("xx.dat") 
repetitions <- 10
rep.factor <- 1000

model[[1]]$model <- "spherical"

#xx <- c(0.9853106,  0.1743256)
#tt <- c(2.414391, 7.243173, 4.828782)

len <- 3
xx <- seq(xx[1], xx[2], len=len)
tt <- c(tt[1:2], diff(tt[1:2])/(len-1))

seq(tt[1], tt[2], tt[3])

setparameters <- function(n.sf) {
  sf <- c(5,10) # linesimufactor
  prnt <- c(1, 6) # >6 gives too much output
  trials <- c(3, 6)
  RFparameters(TBM2.linesimufactor=sf[n.sf],
               TBM3.linesimufactor=sf[n.sf],
               TBMCE.force=TRUE,
               TBMCE.trials=trials[n.sf],
               TBM2.lines=120,#120,
               CE.trials=trials[n.sf],
               TBM3.linesimufactor=sf[n.sf], pch="#", 
               Print=prnt[n.sf], Storing=TRUE, TBM.method="di",
               direct.method=0)
}

ENVIR = environment()
refined.simulation <- function(tol, xx, tt, grid, model, repetitions, method,
                               v, rep.factor) {              
  DeleteAllRegisters()
  assign("z", NULL, envir=ENVIR)
  zaehler <- 0
  dev <- Inf
  setparameters(2)  
  while (abs(dev) > tol) {
    zaehler <- zaehler + 1
    cat("\nlarge deviance ")
    assign("z", cbind(matrix(GaussRF(xx, T=tt, grid=grid, gridtriple=TRUE,
                              model=model, n=repetitions, method=method)
                      , ncol=repetitions), z), envir=ENVIR)
    RFparameters(Print=1)
    
    e <- 0.5 * mean((z[1,] - z[nrow(z), ])^2)
    dev <- if (abs(v)>1e-10) (e-v)/v else e * 100000
    cat("; simu=", e, " (", formatC(dev * 100, dig=2),
        "%, ",  formatC((0.5*(diff(z[c(1,4), 1]))^2-v)/v * 100, dig=2),
        "; ", ncol(z), ") ",
        paste(formatC(z[c(1,4), 1], dig=3), collapse=", "), sep="")
    if (zaehler > rep.factor && abs(dev) > tol) {
      str(model)
      cat("too large deviance:", abs(e-v), ">", tol * v,"\n")
      DeleteAllRegisters()
      return(FALSE)
    }
  } # while dev > tol
  DeleteAllRegisters()
  return(TRUE)
}

e <- NULL
tol <- 0
vv <- v
ee <- e
v <- Variogram(cbind((0:(len-1)) * (xx[2] - xx[1]),
                     (0:(len-1)) * tt[3] / (len-1)),
               model=model)
cat("\nRECHECK: \n", model$model, "grid=",grid, "; d=", d,
    "; T=", !is.null(tt),
    "; ani=", anisotropy, "; k=", paste(nk, collapse=","),
     " ", method, sep="")
print(v)
str(model); str(xx); str(tt)
print(ncol(xx))

if (refined.simulation(tol, xx, tt, grid, model, repetitions, method,
                       v, rep.factor)) {
  ## file.remove("xx.dat")
  stop("OK")
}



} ## if (FALSE)
