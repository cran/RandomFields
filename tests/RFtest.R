
# R <  RFtest.all.R && R <  RFtest.all.R &&R <  RFtest.all.R &&R <  RFtest.all.R &&R <  RFtest.all.R


if (EXTENDED.TESTING <- file.exists("source.R")) {  source("source.R")
} else if (file.exists(f <- "~/R/RF/RandomFields/tests/source.R")) source(f) 
       
## how many simulation methods are available?
MAXMETHODNUMBER <- length(GetMethodNames())

ENVIR <- environment()
assign("zaehler", 0, envir=ENVIR)

pid <- function() 1

getxyz <- function(grid,dim,pointnumber,fieldsize,quadraticgrid) {
  locations <- list()
  if (grid) {      
    locations$lx <- 3;
    pp <- trunc(pointnumber^(1/dim))
    step<-runif(1,max=fieldsize/pp);
    ## step<-bin[2]-bin[1];  ### step <- 0.3 ###
    locations$x<-c(0,step*(pp-0.5),step);
    locations$y<-NULL;  locations$z<-NULL;
    if (dim>1) {
      if (quadraticgrid) {
        locations$y<-locations$x
        if (dim>2) { locations$z <- locations$x }
      } else {
        step<-runif(1,max=fieldsize/pp);
        ## step<-bin[2]-bin[1]; ### step <- 0.3
        locations$y<-c(0,step*(pp-0.5),step);
        if (dim>2) {
          step<-runif(1,max=fieldsize/pp);
          ## step<-bin[2]-bin[1]; ### step <- 0.3
          locations$z<-c(0,step*(pp-0.5),step);
        }
      }
    }
  } else {
    locations$lx <- pointnumber
    locations$x<-runif(pointnumber,max=fieldsize);
    ##locations$x<-(1:pointnumber)*bin[2]-bin[1]; ###
    locations$y<-NULL;  locations$z<-NULL;
    if (dim>1) {
      locations$y <-runif(pointnumber,max=fieldsize)
      if (dim>2) { locations$z <- runif(pointnumber,max=fieldsize) }
    }
  }
  return(locations)
}


RFcontrol <- function (model,kappa1=NULL,kappa2=NULL,kappa3=NULL,
                       nugget=0,mean=0,
                       variance=1,
                       endofbins=1.3, numberbins=30, scaling=1,
                       pointrepet=5,valuerepet=5,pointnumber=100,fieldsize=3,
                       histo=FALSE, wait=FALSE,clean=FALSE, waitfinally=FALSE,
                       dimensions=1:3,grids=c(FALSE,TRUE),
                       #methods=0:MAXMETHODNUMBER,
                       methods=c(0,2,3,4,5,7),
                       ps=NULL,save=TRUE,quadraticgrid=TRUE){
  ## model : covariance model
  ## kappaX: parameters of the model
  ## nugget: amount of nugget effect (variance of the field without nugget
  ##         is always 1)   
  ## endofbins: when calculating the empirical variogram, what is the largest
  ##            distance considered? This is given as endofbins
  ## numberbins: number of bins;
  ## scaling   : if null then automatic detection of scale, so that
  ##             cov(1/scale)~=0.05
  ## pointrepet: number of times, a set locations is simulated
  ## valuerepet: number of times, a set of values is simulated for a given set
  ##             of locations thus the total number of simulated random fields
  ##             is pointrepet * valuerepet
  ## pointnumber: number of locations in the set of locations (if grid, then
  ##              this is an approximate number)
  ## fieldsize : in units of scaling, in each direction
  cat("\n", model)
  MethodNames <- GetMethodNames()

  
  
  model<-model;kappa1<-kappa1;kappa2<-kappa2;kappa3<-kappa3;
  nugget<-nugget; mean<-mean; variance<-variance;
  endofbins<-endofbins; numberbins<-numberbins; scaling<-scaling;
  pointrepet<-pointrepet;valuerepet<-valuerepet;
  pointnumber<-pointnumber;fieldsize<-fieldsize;
  histo<-histo; wait<-wait;clean<-clean; waitfinally<-waitfinally;
  dimensions<-dimensions;grids<-grids;methods<-methods;
  ps<-ps; quadraticgrid<-quadraticgrid;

  if (isnullpar <- is.null(kappa3)) {
    kappa3 <-0
    if (is.null(kappa2)) {
      isnullpar <- 2
      kappa2 <- 0
      if (is.null(kappa1)) {
        isnullpar <- 3
        kappa1 <- 0
      }
    }
  }
  if (isnullpar==3) par.index <- NULL
  else par.index <- 1:(3-isnullpar)
  
  if (!is.numeric(methods)) {
    methnr <- pmatch(methods, MethodNames, dup=TRUE) - 1
    if (any(is.na(methnr))) {
      stop(paste("Method `", methods[is.na(methnr)],
                 "' not identifiable...",sep=""))
      ERR
    }
    methods <- methnr
  }
  
  Pid <- pid()


  if (!is.null(ps) && save) {
    runif(1) # to make sure the random number generator is initialized
    RFpar <- RFparameters()
    save(file=paste(ps,".seed",sep=""),.Random.seed,
         model,kappa1,kappa2,kappa3,nugget, mean,variance,
         endofbins, numberbins, scaling,pointrepet,valuerepet,pointnumber,
         fieldsize,histo, wait,clean, waitfinally, dimensions,grids,methods,
         ps,RFpar,Pid
         )
  }
  
  totalrepet <- pointrepet * valuerepet;
  meanvar<-"" ## to avoid an error in the final output in case everything fails
  for (k1 in kappa1) for (k2 in kappa2) for (k3 in kappa3) for (nug in nugget){
    print(modeletc<-paste(model," m=",format(mean,dig=2),
                          ",v=",format(variance,dig=2),
                          ",ng=", format(nugget,dig=2),
                          ",scl=",format(scaling,dig=2),
                          ",k1=",format(k1,dig=2),
                          ",k2=",format(k2,dig=2),
                          ",k3=",format(k3,dig=2),sep=""))
    param<-as.double(c(mean=mean,variance=variance,nugget=nug,scaling,
                       c(k1,k2,k3)[par.index] ))
    
    bin <- c(-1, seq(0, endofbins, length=numberbins))
    
    midbin <- as.double(0.5 * (bin[-length(bin)] + bin[-1]))
    midbinP0 <- as.double(c(0,midbin))    

    for (dim in dimensions) {
      print("dim")
      print(dim)
      for (grid in grids) {
        print("grids")
        print(grid)
        print(dimgrid<-paste("dim=",dim,",grid=",grid,", pid=",Pid,
                            ", PR=",RFparameters()$PracticalRange,
                            ", qg=",quadraticgrid,sep=""))
        v <- matrix(0, nrow=numberbins, ncol=MAXMETHODNUMBER+1)
        repet <- matrix(0, nrow=numberbins, ncol=MAXMETHODNUMBER+1)
        MethodIgnoreList <- -99
#
#        MethodIgnoreList <- c(0,2,3,4,5)
        for (i in 1:pointrepet) {
          print(c("pointrepet", pointrepet, i))
          locations <- getxyz(grid,dim,pointnumber,fieldsize,quadraticgrid)
          for (method in methods) {
            #print(method); print(methods); print(MethodIgnoreList); readline()
            if (sum(MethodIgnoreList==method)==0) {
              methodname <- MethodNames[method + 1]      
              cat("\n !! ", methodname)
              ##print(date())
##seed <-  get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)             
##xx <- list(x=locations$x,y=locations$y,z=locations$z,
##           grid=grid,model=model,
##            param=param,meth=methodname,
##           reg=0,gridtriple=TRUE,
##           PracticalRange=RFparameters()$PracticalRange,
##           dim = dim, midbinP0=midbinP0,
##           seed=seed)
##save(file="xx", xx)                      
              error <- try(InitGaussRF(locations$x,y=locations$y,z=locations$z,
                                   grid=grid,model=model,
                                   param=param,meth=methodname,
                                   reg=0,gridtriple=TRUE))
              if (!is.numeric(error)) error <- 10000
              if (save) {
                savename<-paste(ps,".save",sep="")              
                save(file=savename,locations,grid,model,param,methodname,
                     valuerepet,.Random.seed)
              }
              if (error) {
                print("Initialisation failed");
                if (error<=4) {
                  print(paste("error=",error,
                              "; Ignoring further point generations"))
                  MethodIgnoreList <- c(MethodIgnoreList,method)
                } else print(error)
              } else {            
                ##print(date())                  
                print("do simulate") ##################
                allres <- DoSimulateRF(n=valuerepet,reg=0)
                if (!is.null(allres)) {
                  if (any(is.na(allres))) {
                    print(allres)
                    stop("some NA in allres"); ERR
                  }
                  meanvar<-paste("m=",format(mean(allres),dig=2),
                                 ", sill=",format(var(as.double(allres)),dig=2))
                  print(paste(i,", Method=",method,":",meanvar))
                  if (histo) {dev.set(dev.next()); hist(as.vector(allres))}
                  binresult<-
                    EmpiricalVariogram(locations$x, locations$y, locations$z,
                                       data=allres,
                                       grid=grid, bin=bin, gridtriple=TRUE)
                  print("done")
                  index <- is.finite(binresult$e);
                  print("Done")
                  v[index,method+1] <- v[index, method+1] + binresult$e[index];
                  print("done.")
                  repet[index,method+1] <- repet[index,method+1] + 1;
                  print("Done.")
                } else { print(" Simulation failed "); }
              }
            } # method
          } # if method should not be ignored ...
        } # pointrep
        v <- v / repet;
    print(model)
    print(param)
    print(dim)
        truevariogram  <-
          try(Variogram(midbinP0, model, param, dim))
        if (!is.numeric(truevariogram)) next
        delta <- colSums(abs(v-truevariogram[-1]),na.rm=TRUE)
        delta[apply(is.na(v), 2, all)]<-NA;
        assign("zaehler", zaehler + 1, envir=ENVIR)
        filename <- paste(ps,".",zaehler,".ps",sep="")
        print(filename)
        print(modeletc)
        if (!is.null(ps)) {
          Dev(TRUE,TRUE,ps=filename,height=6,width=6,paper="special")
        }
        xlab <- ""
        for (i in 1:length(delta)) {          
          if (!is.na(delta[i])) {
            xlab <-
              paste(xlab,i-1,":",
                    paste(repet[1:(min(numberbins,10)),i],collapse=","),
                          " ",sep="")
          }
        }
        if (any(!is.na(delta))) {
          if (any(!is.na(truevariogram)))
            ymax <- max(truevariogram,na.rm=TRUE)*1.05
          else
            ymax <- max(v,na.rm=TRUE)*1.05
          matplot(midbin,v,col=c(1,2,3,4,5,6,7,1,2,3), type="b",
                  pch=c("0","1","2","3" ,"4","5","6","7","8","9","A","B","C"),
                  lty=c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3),
                  ylim=c(0,ymax),
                  ylab=paste("(",dimgrid,"    ",meanvar,")",sep=""),
                  xlab=xlab)
          points(midbinP0,truevariogram,pch="*")
          legend(x=midbin[length(midbin)],y=0,xjust=1,yjust=0,
                 legend=c("true","CE","local","TBM2","TBM3","spectral","direct",
                   "nugget","add.MPP","HypPl","other"),
                 col=c(1,1,2,3,4,5,6,7,1,2,3),
                 pch=c("*","0","1","2","3" ,"4","5","6","7","8","9","A","B","C")
                 )
        } else {
          print("no graphic")
          if (any(is.finite(truevariogram))) {
            print(truevariogram)
            plot(midbinP0,truevariogram,pch="X",col=2,
                 ylab=paste("(",dimgrid,"    ",meanvar,")",sep=""),
                 xlab=xlab)
            text(midbinP0[length(midbinP0)]/2,max(truevariogram)/2,
                 pos=1,"all simulation methods failed");
          } else {
            plot(0,0,pch="X",col=2,
                 ylab=paste("(",dimgrid,"    ",meanvar,")",sep=""),xlab=xlab)
            text(0,0,pos=1,"variogramm does not exist!");
          }
        }
        if (!interactive())
          title(sub=paste("delta=",paste(format(delta,dig=2),collapse=","),
                  sep=""), main=paste(modeletc))
        if (!is.null(ps)) { Dev(FALSE,FALSE) } else if (wait) readline()
      }
    }
  }
  if (waitfinally && !wait)  { print("press key");readline(); }
    if (clean) while (!is.null(dev.list())) dev.off();
  return(NULL);
}





