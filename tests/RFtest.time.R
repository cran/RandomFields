## R --no-save < RFtest.time.R
# source("RFtest.time.R")

wait <- TRUE; ##wait <- FALSE;  

RFsimu <- function(cov="expo", param =c(1000,2,0,2/3,1),x=seq(0,10,0.2),
                   y=seq(0,10,0.1),key=0)
{
  print(paste("TIME=",paste(format(unix.time(res<-GaussRF(x=x,y=y,z=NULL,
                                                          grid=TRUE,
                                                          model=cov,
                                                          param=param,
                                                          method=NULL,
                                                          n=1,reg=key)
                                             ),dig=3))),quote=FALSE)
  if (is.null(res)) {
    print("Simulation has failed",quote=FALSE)
  } else {
    image(x,y,matrix(res,ncol=length(y)),col=grey(0:100 / 100));
    title(cov);
  }
}

library(RandomFields)
RFparameters(PrintLevel=4,
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
print("Example 2: Gneiting, 1) working directly 2)enlarging necessay 3)failed",quote=FALSE)
RFparameters(Storing=TRUE,PrintLevel=2)
RFsimu("gneiting",param=c(0,1,0,10));
RFsimu("gneiting",param=c(0,1,0,20));
RFsimu("gneiting",param=c(0,1,0,60));







