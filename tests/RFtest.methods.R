## R --no-save < RFtest.methods.R
# source("RFtest.methods.R")

#library(RandomFields,lib="~/TMP")
library(RandomFields)
        
ERRORCOVNOTALLOWED <- 4 ## check with RFsimu.h!!!

models<- list(list("bessel",1),
              list("cauchy",1),
              list("cauchytbm",1:3),
              list("circular",NULL),
              list("cubic",NULL),
              list("cone",1:3),
              list("exponential",NULL),
              ##  "expP", 
              list("gauss",NULL),
              list("gencauchy",1:2),
              list("gneiting",NULL),
              list("gneitingdiff",1:2),
              list("gengneiting",1:2),
              list("holeeffect",1),
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
              
methods <- c("cir","loc","TBM2","TBM3","sp","dir","add")

working <- matrix(0, nrow=length(models), ncol=length(methods))
RFparameters(PrintLevel=0) ## 0
for (scale in c(0.3,1,3)) for (kappa1 in c(0.5,1,2,10)) {
  for (mo in 1:length(models)) {
    param <- c(0,1,1,scale, c(kappa1,2,2)[models[[mo]][[2]]])
    param <- c(0,1,1,scale, c(kappa1,2,3)[models[[mo]][[2]]])
    ## cauchytbm : kappa3 >= 3  for tbm2
    print(models[[mo]][[1]])
    for (me in 1:length(methods)) {

     error <- InitGaussRF(x=1:10,y=1:10,grid=TRUE,model=models[[mo]][[1]],
                             param=param,method=methods[me])
     if (error==ERRORCOVNOTALLOWED)
       error <- InitGaussRF(x=1:10,grid=TRUE,model=models[[mo]][[1]],
                               param=param,method=methods[me])
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





