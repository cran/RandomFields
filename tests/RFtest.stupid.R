## RX -o "testlevel<-2" -i RFtest.stupid.R  | R > stupid.out2
## R --no-save < RFtest.stupid.R
# source("RFtest.stupid.R")

source("./RFtest.R")

if (!exists("testlevel")) testlevel <- 1
if (!exists("ps")) ps <- NULL

maxtestlevel <- 5
switch (testlevel,   
        {NP <- 1; repeatscript <-1; },
        {NP <- 2; repeatscript <-500; }
        )


pointrepet <- 1;
valuerepet <- 3;
histo <- FALSE;
save <- FALSE;
RFparameters(PrintLevel=3)


getxyz <- function(grid,dim,pointnumber,fieldsize,quadraticgrid) {
  locations <- list()
  if (grid) {      
    locations$lx <- 3;
    pp <- trunc(pointnumber^(1/dim))
    if (runif(1)<0.1) step<-0 else step<-runif(1,-1,4);
    ## step<-bin[2]-bin[1];  ### step <- 0.3 ###
    locations$x<-c(0,step*(pp-0.5),step);
    locations$y<-NULL;  locations$z<-NULL;
    if (dim>1) {
      if (quadraticgrid) {
        locations$y<-locations$x
        if (dim>2) { locations$z <- locations$x }
      } else {
        if (runif(1)<0.1) step<-0 else step<-runif(1,-1,4);
        ## step<-bin[2]-bin[1]; ### step <- 0.3
        locations$y<-c(0,step*(pp-0.5),step);
        if (dim>2) {
          if (runif(1)<0.1) step<-0 else step<-runif(1,-1,4);
          ## step<-bin[2]-bin[1]; ### step <- 0.3
          locations$z<-c(0,step*(pp-0.5),step);
        }
      }
    }  
  } else {
    locations$x<-runif(pointnumber,max=fieldsize);
    ##locations$x<-(1:pointnumber)*bin[2]-bin[1]; ###
    locations$y<-NULL;  locations$z<-NULL;
    if (dim>1) {
      locations$y <-runif(pointnumber,max=fieldsize)
      if (dim>2) { locations$z <- runif(pointnumber,max=fieldsize) }
    }
  }
  if (runif(1)<0.15) {locations$x <- NULL; print("X IS NULL") }
  if (dim>1) {
    if (runif(1)<0.15) {locations$y <- NULL; print("Y IS NULL") }
    if ((dim>2) && (runif(1)<0.15)) {locations$z<- NULL; print("Z IS NULL") }
  } else {
    if (runif(1)<0.1) { locations$z <- locations$x; print(" Z=X, Y=NULL")}
    else if (runif(1)<0.1) {locations$y <- runif(runif(1,0,5),-1,5); print(" NONSENSE Y")}                 
  }
  print(locations)
  return(locations)
}

ENVIR <- environment()
randomize <- function(){
  assign("quadraticgrid", runif(1)<0.5, envir=ENVIR) ##quadraticgrid <- FALSE
  assign("Mean", runif(1,-10,10), envir=ENVIR)
  assign("scaling", if (runif(1)<0.5) runif(1,1,10) else runif(1,-1,1),
         envir=ENVIR)
  assign("variance", if (runif(1)<0.5) runif(1,1,10) else runif(1,-1,1),
         envir=ENVIR)
  assign("nugget", if (runif(1)<0.5) runif(1,1,10) else runif(1,-1,1),
         envir=ENVIR)
  assign("fieldsize", if (runif(1)<0.5) runif(1,1,10) else runif(1,-1,1),
         envir=ENVIR)
  assign("endofbins", if (runif(1)<0.5) runif(1,1,10) else runif(1,-1,1),
         envir=ENVIR)
  assign("numberbins", as.integer(runif(1,0,5)), envir=ENVIR)
  print(paste("qg=",format(quadraticgrid,dig=3),
              " m=",format(Mean,dig=3),
              " s=",format(scaling,dig=3),
              " v=",format(variance,dig=3),
              " n=",format(nugget,dig=3),
              " f=",format(fieldsize,dig=3),
              " eb=",format(endofbins,dig=3),
              " nb=",format(numberbins,dig=3)))
}

simplemodels <- c("gauss","penta");

models <- list(
 list(model="whittle",     kappa1=runif(NP,-0.5,3),   kappa2=NULL),
 list(model="cauchy",      kappa1=runif(NP,-0.5,3),   kappa2=NULL),
)

largemodels <-
  list(
       list(model="cauchytbm", kappa1=runif(NP,-0.5,2.5),
            kappa2=runif(NP,-0.5,3),  kappa3=round(runif(NP,0,50)/10)),
       list(model="hyper",    kappa1=runif(NP,-0.5,3),
            kappa2=runif(NP,-0.5,3), kappa3=runif(NP,-0.5,3))
       )

for (i in 1:repeatscript) {
  for (naturalscaling in FALSE:TRUE) for (pointnumber in 0:5) {
    print(c("nsc,pn",naturalscaling,pointnumber))
    RFparameters(PracticalRange=naturalscaling)
    
      if (FALSE) for (model in models) {
      cat("\n ********* models "); str(model)
      randomize()
      x <- try(RFcontrol(model$model,pointnumber=pointnumber,valuerepet=valuerepet,
                pointrepet=pointrepet,kappa1=model$kappa1,kappa2=model$kappa2,
                scaling=scaling,var=variance,nug=nugget,mean=Mean,
                field=fieldsize,endof=endofbins,numberb=numberbins,
                histo=histo,ps=ps,quadraticgrid=quadraticgrid,save=save)
          )
      if (!is.null(x)) { readline();}
    }
    
   for (model in simplemodels) {
      cat("\n ********* simplemodels "); str(model)
      randomize()
      x <- try(RFcontrol(model,pointnumber=pointnumber,valuerepet=valuerepet,
                pointrepet=pointrepet,scaling=scaling,var=variance,nug=nugget,
                mean=Mean,field=fieldsize,endof=endofbins,numberb=numberbins,
                histo=histo,ps=ps,quadraticgrid=quadraticgrid,save=save)
          )
      if (!is.null(x)) { readline();}
    }
    
    for (model in largemodels) {
      cat("\n ********* largemodels "); str(model)
      randomize()
      x<-try(RFcontrol(model$model,pointnumber=pointnumber,valuerepet=valuerepet,
                pointrepet=pointrepet,kappa1=model$kappa1,kappa2=model$kappa2,
                kappa3=model$kappa3,scaling=scaling,var=variance,nug=nugget,
                mean=Mean,field=fieldsize,endof=endofbins,numberb=numberbins,
                histo=histo,ps=ps,quadraticgrid=quadraticgrid,save=save)
          )
        if (!is.null(x)) { readline();}
    }    
  }
}

q()







