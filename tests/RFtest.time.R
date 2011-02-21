## R --no-save < RFtest.time.R
# source("RFtest.time.R")

#library(RandomFields)


## unten sind 2 Beispiele ausgeblendet, da zu lange brauchen --
## was ist los??


if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

wait <- TRUE; ##wait <- FALSE;  

RFsimu <- function(cov="expo", param =c(1000,2,0,2/3,1),
                   x=seq(0, 10, if (interactive()) 0.2 else 1),
                   y=seq(0, 10, if (interactive()) 0.1 else 1),
                   key=0)
{
  #str(RFparameters())
  ut <- unix.time(res<-try(GaussRF(x=x, y=y, z=NULL, grid=TRUE, model=cov,
                               param=param, method=NULL, n=1, reg=key)))
  cat("TIME=",  paste(format(ut, dig=3)), "\n")
  if (is.null(res) || !is.numeric(res)) {
    print("Simulation has failed",quote=FALSE)
  } else {
    image(x,y,matrix(res,ncol=length(y)),col=grey(0:100 / 100));
    if (!interactive()) title(cov);
  }
}

RFparameters(PrintLevel=2,
             CE.force=FALSE,CE.mmin=0,CE.tolRe=-1e-7,CE.tolIm=1e-3,
             CE.trials=3)

#####################
# Time checks
####################

print("",quote=FALSE);print("Press RETURN to start",quote=FALSE);
if (wait) readline()

print("",quote=FALSE);
print("STORING TRUE VS FALSE -- NOTE THE DIFFERENCE IN TIME FOR THE RESPECTIVE SECOND SIMULATION",quote=FALSE)
print("Storing=false:",quote=FALSE)
RFparameters(Storing=FALSE)
RFsimu("cubic",param=c(0,1,0,2));
RFsimu("cubic",param=c(0,1,0,2));

print("",quote=FALSE);print("Storing=true:",quote=FALSE)
RFparameters(Storing=TRUE)
RFsimu("cubic",param=c(0,1,0,2));
RFsimu("cubic",param=c(0,1,0,2));

print("Press RETURN",quote=FALSE);if (wait) readline()
print("",quote=FALSE);
print("STORING TWO SIMULATIONS-- NOTE THE DIFFERENCE IN TIME FOR 3RD & 4TH SIMULATION",quote=FALSE);
print("two simulations at the same time are SLOW, if `key' is not used not:",quote=FALSE)
RFsimu("spherical",param=c(0,1,0,1));
RFsimu("cubic",param=c(0,1,0,2));
RFsimu("spherical",param=c(0,1,0,1));
RFsimu("cubic",param=c(0,1,0,2));

print("",quote=FALSE);print("",quote=FALSE);
print("two simulations at the same time, FAST:",quote=FALSE)
RFsimu("spherical",param=c(0,1,0,1));
RFsimu("cubic",param=c(0,1,0,2),key=1);
RFsimu("spherical",param=c(0,1,0,1));
RFsimu("cubic",param=c(0,1,0,2),key=1);

print("Press RETURN",quote=FALSE);if (wait) readline()
print("",quote=FALSE);
print("EXCEPTION HANDLING, 2 EXAMPLES; LOOK AT THE (ENABLED) TRACING",quote=FALSE);
print("Exapmle 1; Waves, once by circular embedding, once by TBM ",quote=FALSE)
RFparameters(Storing=TRUE,PrintLevel=2)
RFsimu("wave",param=c(0,1,0,0.01));
RFsimu("wave",param=c(0,1,0,2));
print("Press RETURN for second example",quote=FALSE);if (wait) readline()

print("",quote=FALSE);print("",quote=FALSE);
print("Example 2: Gneiting, 1) working directly 2)enlarging necessary 3)failed",quote=FALSE)
RFparameters(Storing=TRUE,PrintLevel=2)
RFsimu("gneiting",param=c(0,1,0,10));
RFsimu("gneiting",param=c(0,1,0,15));

RFparameters(Storing=TRUE,PrintLevel=5)
RFsimu("gneiting",param=c(0,1,0,20));







