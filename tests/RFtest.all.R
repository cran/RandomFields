### RX -o 'ps<-"x"' -i RFtest.all.R  | R  > xa.out &
### RX -o 'testlevel<-3;ps<-"x"' -i RFtest.all.R | R 
## R --no-save < RFtest.all.R
# source("RFtest.all.R")

if (file.exists(f <- "RFtest.R")) source(f) else
if (file.exists(f <- "~/R/RF/RandomFields/tests/RFtest.R")) source(f) 

if (!exists("testlevel"))
  testlevel <- 1

maxtestlevel <- 5
switch (testlevel,   
        {pointnumber<-20; valuerepet<-3;  pointrepet<-1; NP<-1;
         duration<-"70 sec"; random.parameters <- TRUE; PrintLevel <- 1;
         lines<-100;repeatscript <- 1 }, 
        {pointnumber<-50; valuerepet<-30;  pointrepet<-4; NP<-1;
         duration<-"2.5 h"; random.parameters <- TRUE; PrintLevel <- 2;
         lines<-300;repeatscript <- 1 }, 
        {pointnumber<-50; valuerepet<-200;  pointrepet<-8; NP<-3;
         duration<-"?? d"; random.parameters <- TRUE; PrintLevel <- 3;
         lines<-500;repeatscript <- 1},
        {pointnumber<-200; valuerepet<-300;  pointrepet<-12; NP<-3;
         duration<-"?? w"; random.parameters <- TRUE; PrintLevel <- 3;
         lines<-1000;repeatscript <- 3 },
        {pointnumber<-200; valuerepet<-300;  pointrepet<-12; NP<-3;
         duration<-"?? w"; random.parameters <- FALSE; PrintLevel <- 3;
         lines<-1000;repeatscript <- 3 }
        # further test: random point number range [0,10]
        )

#******************************************************************
#*                                                                *
#* This test tries to simulate 1, 2, and 3 dimensional random     *
#* fields of any defined covariance function with any defined     *
#* method. The random fields may be defined on a grid or not.     *
#*                                                                *
#* For each combination, several random fields are generated      *
#* and the empirical variogram is calculated and (graphically)    *
#* compared with the true variogram.                              *
#* In some dimensions or for some (randomly chosen) parameters    *
#* the variogram (thus the random field) does not exist. In this  *
#* case, a single red cross in the middle of the picture is       *
#* shown.                                                         *
#* For some parameter values it may happen that none of the       *
#* simulation methods are suitable. In this case the true         *
#* variogram is shown in form of red crosses.                     *
#* In case the random field can be simulated by a certain method  *
#* the estimated variogram is shown. The latter may nonetheless   *
#* deviate largely from the true variogram (depending on the      *
#* parameters and the precision chosen)                           *
#*                                                                *
#* Please ignore the strange titles and labels, which code al l   *
#* the relevant parameter values, and some summary statistic      *
#* of the simulation.                                             *
#*                                                                *
#* As there are about 20 covariance functions, some of them       *
#* allowing for additional parameters, and 5 methods implemented  *
#* the number of combinations is quite large and the running      *
#* time can be pretty high!                                       *
#* Therefore different test levels exist. Low level means quick   *
#* testing, but the estimated variogram may be very eratic.       *
#*                                                                *
#* The estimated time for running the test on a 500Mhz PC is      *
#* given by the variable `duration' above                         *
#*                                                                *
#******************************************************************"


if (exists("ps") && !is.null(ps)) { ## for author's use only
  ps <- paste(ps,pid(),sep=".")
  save <- TRUE;
} else {  ## for users
  ps <- NULL;
  histo<-FALSE;
  save <- FALSE;
}


### minor adjustments ...
fieldsize<-3 
endofbins <- 1.3
numberbins <-31
if (!is.null(ps))  {histo <- FALSE}  ## Da R CMD check die Variable `ps' nicht
                                    ## definiert hat,
if (histo) {X11(); X11();}         ## wird X11() nicht mehr aufgerufen
RFparameters(Storing=TRUE,PrintLevel=PrintLevel,
             TBM2.lines=6/5*lines,TBM2.linesimufactor=0,TBM2.linesimustep=0.01,
             TBM3D2.lines=lines,TBM3D2.linesimufactor=0,TBM3D2.linesimustep=0.01,
             TBM3D3.lines=lines,TBM3D3.linesimufactor=0,TBM3D3.linesimustep=0.01,
             spectral.lines=lines)

ENVIR <- environment()
randomize <- function(){
  assign("quadraticgrid", runif(1)<0.5, envir=ENVIR) ##quadraticgrid <- FALSE
  assign("Mean", runif(1,0,10), envir=ENVIR)
  assign("scaling", if (runif(1)<0.5) runif(1,1,10) else runif(1,0.2,1),
         envir=ENVIR)
  assign("variance", if (runif(1)<0.5) runif(1,1,10) else runif(1,0.2,1),
         envir=ENVIR)
  assign("nugget", if (runif(1)<0.5) runif(1,1,10) else runif(1,0.2,1),
         envir=ENVIR)
}

simplemodels <- c("circular","cubic",
                  "exponential","gauss",
                  "gneiting","penta",
                  "spherical","wave");

if (random.parameters){
models <-
  list(
       list(model="bessel",  kappa1=runif(NP,-1,3),      kappa2=NULL),
       list(model="cauchy",  kappa1=runif(NP,-0.5,3),    kappa2=NULL),
       list(model="damped",  kappa1=runif(NP,-0.5,2.5),  kappa2=NULL),
       list(model="FD",      kappa1=runif(NP,-1,1),      kappa2=NULL),
       list(model="fractg",  kappa1=runif(NP,-0.5,2.5),  kappa2=NULL),
       list(model="gencau",  kappa1=runif(NP,-0.5,2.5),
            kappa2=runif(NP,-0.5,3)),
       list(model="gengn",   kappa1=round(runif(NP,-0.5,4)),
            kappa2=runif(NP,-0.5,3)),
       list(model="lgd1",    kappa1=runif(NP,-0.5,3.5),
            kappa2=runif(NP,-0.5,1.5)),
       list(model="power",   kappa1=runif(NP,-0.5,3),   kappa2=NULL),
       list(model="qexpone", kappa1=runif(NP,-0.2,1.2), kappa2=NULL),
       list(model="stable",  kappa1=runif(NP,-0.5,2.5), kappa2=NULL),
       list(model="whittle", kappa1=runif(NP,-0.5,3),   kappa2=NULL),
       list(model="fractalB",kappa1=runif(NP,-0.5,2.5), kappa2=NULL),
)
largemodels <-
  list(
       list(model="cone", kappa1=runif(NP,0.5,1),
            kappa2=runif(NP,0.5,3),  kappa3=runif(NP,0.5,3)),
       list(model="cone", kappa1=runif(NP,-0.5,2.5),
            kappa2=runif(NP,-0.5,3),  kappa3=runif(NP,-0.5,3)),
       list(model="cauchytbm", kappa1=runif(NP,-0.5,2.5),
            kappa2=runif(NP,-0.5,3),  kappa3=round(runif(NP,0,50)/10)),
        list(model="hyper",    kappa1=runif(NP,-0.5,3),
            kappa2=runif(NP,-0.5,3), kappa3=runif(NP,-0.5,3))
       )

## very large models are not tested yet:
##      nsst, nsst2

} else { ## fixed random parameters
models <-
  list(
       list(model="bessel",  kappa1=c(-1,-0.5,0,0.5,4),kappa2=NULL),
       list(model="cauchy",  kappa1=c(-1,0,0.1,5),     kappa2=NULL),
       list(model="damped",  kappa1=c(-1,0,1,2),       kappa2=NULL),
       list(model="FD",      kappa1=c(-1,-0.5,0,0.5),  kappa2=NULL),
       list(model="fractg",  kappa1=c(0,1,2,3),        kappa2=NULL),
       list(model="gencau",  kappa1=c(-1,0,1,2),       kappa2=c(-1,0,0.1,5)),
       list(model="gengn",   kappa1=c(-1,0,0.1,1,8),   kappa2=c(-1,0,0.1,4)),
       list(model="lgd1",    kappa1=runif(0,1,2,3,4),  kappa2=c(-1,0,1)),
       list(model="power",   kappa1=c(-1,0.1,1,2),     kappa2=NULL),
       list(model="qexpon",  kappa1=c(-1,0,0.5,1),     kappa2=NULL),
       list(model="stable",  kappa1=c(-1,0,0.1,1,2),   kappa2=NULL),
       list(model="whittle", kappa1=c(-1,0,0.1,1,8),   kappa2=NULL),
       list(model="fractalB",kappa1=c(-1,0,1,2,3),     kappa2=NULL),
      )
largemodels <-
  list(
       list(model="cauchytbm", kappa1=c(-1,0,0.1,1,8), kappa2=c(-1,0,0.1,5),
            kappa3=c(0.5,1,1.5,2,3)),
       list(model="cone", kappa1=c(-1,0,8), kappa2=c(-0.5,0,5),
            kappa3=c(-0.5,0,3)),
       list(model="hyper",  kappa1=c(-0.1,0,0.1,5), kappa2=c(-5,-0.1,0,0.1,5),
            kappa3=c(-0.1,0,0.1,5))
       )

## very large models are not tested yet:
##      nsst, nsst2

} 

for (i in 1:repeatscript) {
  randomize()
  print("nugget only simulates the nugget effect. Therefore the value of variance is not taken into account");
  
  RFcontrol("nugget",pointnumber=pointnumber,valuerepet=valuerepet,
            pointrepet=pointrepet,scal=1,var=0,nug=nugget,mean=Mean,
            field=fieldsize,endof=endofbins,numberb=numberbins,histo=histo,
            ps=ps,quadraticgrid=quadraticgrid,save=save)
    
  for (naturalscaling in FALSE:TRUE) {
    RFparameters(PracticalRange=as.logical(naturalscaling))

    for (model in models) {##
      randomize()
      RFcontrol(model$model,pointnumber=pointnumber,valuerepet=valuerepet,
                pointrepet=pointrepet,kappa1=model$kappa1,kappa2=model$kappa2,
                scaling=scaling,var=variance,nug=nugget,mean=Mean,
                field=fieldsize,endof=endofbins,numberb=numberbins,
                histo=histo,ps=ps,quadraticgrid=quadraticgrid,save=save)
    }

    for (model in largemodels) {
      randomize()
      RFcontrol(model$model,pointnumber=pointnumber,valuerepet=valuerepet,
                pointrepet=pointrepet,kappa1=model$kappa1,kappa2=model$kappa2,
                kappa3=model$kappa3,scaling=scaling,var=variance,nug=nugget,
                mean=Mean,field=fieldsize,endof=endofbins,numberb=numberbins,
                histo=histo,ps=ps,quadraticgrid=quadraticgrid,save=save)
      print("OK")
    }

    for (model in simplemodels) {
      randomize()
      RFcontrol(model,pointnumber=pointnumber,valuerepet=valuerepet,
                pointrepet=pointrepet,scaling=scaling,var=variance,nug=nugget,
                mean=Mean,field=fieldsize,endof=endofbins,numberb=numberbins,
                histo=histo,ps=ps,quadraticgrid=quadraticgrid,save=save)
    }
        
  }
}







