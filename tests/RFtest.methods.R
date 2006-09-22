## R --no-save < RFtest.methods.R
# source("RFtest.methods.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")
RFparameters(PrintLevel=6, TBM2.num=EXTENDED.TESTING) ## 0
       
ERRORCOVNOTALLOWED <- 4 ## check with RFsimu.h!!!
ERRORCOVFAILED <- 14

firstmodel <- if (interactive()) 1 else 18
#firstmodel <- 23
firstmethod <- 1
#firstmethod <- 3

models<- list(list("bessel",1),
              list("cauchy",1),
              list("cauchytbm",1:3),
              list("circular",NULL),
              list("cubic",NULL),
              list("cone",1:3),
              list("dampedcosine",1),
              list("exponential",NULL),
              ##  "expP", 
              list("fractgauss",1),
              list("fractalB",1),
              list("gauss",NULL),
              list("gencauchy",1:2),
              list("gengneiting",1:2),
               list("gneiting",NULL),
              ##list("gneitingdiff",1:2),
              list("hyper",1:3),
              list("nugget",NULL),
              list("penta",NULL),
              list("power",1),
              list("qexponen",1),
              list("spherical",NULL),
              list("stable",1),
              list("wave",NULL),
              list("whittle",1)
              )
              
methods <- c("cir", "cutoff", "intrinsic",
             "TBM2", "TBM3", "sp", "dir", "add", "hyper")

working <- matrix(0, nrow=length(models), ncol=length(methods))
for (scale in c(0.3,1,3)) for (kappa1 in c(0.5,1,2,10)) {
  for (mo in firstmodel:length(models)) {
    param <- c(0,1,1,scale, c(kappa1,2,2)[models[[mo]][[2]]])
    param <- c(0,1,1,scale, c(kappa1,2,3)[models[[mo]][[2]]])
    ## cauchytbm : kappa3 >= 3  for tbm2
    ##cat("\n", scale, kappa1, models[[mo]][[1]])
    for (me in firstmethod:length(methods)) {
      cat(">>>>", scale, kappa1, models[[mo]][[1]], methods[me],"\n")
      ##cat("\n\nSTART")
      error <- InitGaussRF(x=1:10,y=1:10,grid=TRUE,model=models[[mo]][[1]],
                             param=param,method=methods[me])
      ##cat(" E=",error)
      if (error==ERRORCOVNOTALLOWED || error==ERRORCOVFAILED)
        error <- InitGaussRF(x=1:10,grid=TRUE,model=models[[mo]][[1]],
                             param=param,method=methods[me])
      ##cat(" E1=",error)
      if (error==0)
        working[mo,me]  <- working[mo,me] +1      
    } 
  }
}

{
  sign <- c(paste(1:9),"X")
  print(paste(formatC("",width=14),paste(formatC(methods,width=4),collapse="")))
  for (mo in 1:length(models)) {
    txt <- ""
    for (me in 1:length(methods)) {
      if (working[mo,me]==0) txt <- paste(txt," . ")
      else txt <- paste(txt,"  ",sign[min(10,working[mo,me])]," ",sep="")
    }
    print(paste(formatC(models[[mo]][[1]],width=15),txt),quote=FALSE)
  }
  PrintModelList()
}





