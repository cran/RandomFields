
## PrintLevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information

mleRF <-
function(coord,data,model,param,
         lower.kappa=NULL,upper.kappa=NULL,sill=NA,
         ...) {
  mleRF.default(coord=coord,data=data,model=model,param=param,
                lower.kappa=lower.kappa,upper.kappa=upper.kappa,sill=sill,...)
}


## to.do: extension to data measured several times, i.e. data equals matrix!

mleRF.default <-
function(coord,data,model,param,
         lower.kappa=NULL,upper.kappa=NULL,sill=NA,
         use.naturalscaling=TRUE,
         PrintLevel=0, trace.optim=0,
         bins=20,distance.factor=0.5,
         upperbound.scale.factor=10,
         lowerbound.scale.factor=20,
         lowerbound.scale.LS.factor=5,
         upperbound.var.factor=10,
         lowerbound.var.factor=100,
         lowerbound.sill=1E-10,
         scale.max.relative.factor=1000,
         minbounddistance=0.001, minboundreldist=0.02,
         approximate.functioncalls=50,
         pch="*"
) {
  index <- is.na(param)
  if (!any(index)) return(param)
  coord <- as.matrix(coord)
  
  ENVIR <- environment()
  ip<-.C("GetParameterIndices", MEAN= integer(1), VARIANCE= integer(1),
          NUGGET= integer(1), SCALE= integer(1), 
          KAPPA= integer(1), LASTKAPPA= integer(1),integer(1),
          SILL= integer(1), DUP = FALSE)
  MEAN <- ip$MEAN+1
  VARIANCE <- ip$VARIANCE+1
  NUGGET <- ip$NUGGET+1
  SCALE <- ip$SCALE+1
  KAPPA <- ip$KAPPA+1
  LASTKAPPA <- ip$LASTKAPPA+1
  SILL <- ip$SILL+1

  stopifnot(is.character(model),
            all(is.na(param)) || is.numeric(param[!is.na(param)]),
            all(is.finite(coord)),
            all(is.finite(as.matrix(data))))
  covnr <- .C("GetModelNr", as.character(model), nr=integer(1))$nr
  if ((covnr==.C("GetModelNr", as.character("nugget"), nr=integer(1))$nr) &&
      (any(index[-c(MEAN,NUGGET)]) || (index[VARIANCE]!=0.0)))
    stop("if model=nugget effect then the variance must be zero, and only mean and/or nugget can be NA")   
  if (covnr<0) stop("model not ok")
  storage.mode(covnr) <- "integer"
   if (length(param)>LASTKAPPA) stop("parameter vector too long!")
  lower <- -Inf * seq(1, 1, l=length(param))
  lower[NUGGET] <- 0 
  if (!is.null(lower.kappa)) {
    lower[KAPPA:(KAPPA+length(lower.kappa)-1)]<-lower.kappa
  }
  upper <- Inf * seq(1,1,l=length(param)) 
  if (!is.null(upper.kappa)) {
    upper[KAPPA:(KAPPA+length(lower.kappa)-1)]<-upper.kappa
  }
  if (length(param)>=KAPPA) {
    l.kappa <- lower[KAPPA:length(param)]
    u.kappa <- upper[KAPPA:length(param)]
    nice.kappa <- (l.kappa+u.kappa)/2
    if (!all(is.finite(nice.kappa)))
      warning("limits for kappa are not all finite")
    nice.kappa[is.na(nice.kappa)] <- 1.0
    nice.kappa[nice.kappa==Inf] <- l.kappa + 1.0
    nice.kappa[nice.kappa==-Inf] <- u.kappa - 1.0
  } else nice.kappa <- NULL
  lc<- nrow(coord)
  storage.mode(lc) <- "integer"
  lpar <- length(param)
  storage.mode(lpar) <- "integer"
  distances <- as.double(dist(coord))
  storage.mode(distances) <- "double"
  mindistances <- min(distances[distances!=0])
  maxdistances <- max(distances)
  vardata <- var(data)
  variables <- c(0,vardata,0,maxdistances,nice.kappa)
  PARAM <- param ## just to make clear in MLEtarget and LStarget what is global
                 ## and not overwrite param

  if (is.na.mean <- is.na(PARAM[MEAN])) { PARAM[MEAN] <- mean(data) }
  MLEtargetV <- data-PARAM[MEAN] ## needed in the  next line
  ##                                and in MLE if !is.na.mean

  ev <- EmpiricalVariogram(coord, data=MLEtargetV, grid=FALSE, 
                     bin=c(-1,seq(0, distance.factor * maxdistances,
                       len=bins+1)))
  bin.centers <-  ev$c
  binned.variogram <- ev$e
  binned.n <- ev$n
  
  bins <- length(ev$n)
  max.bin.vario <- max(binned.variogram)
  if ( (dim <- as.integer(ncol(coord))) > 3) {
    warning("only the first three columns are considered as coordinates")
    dim <- as.integer(3)
  }

  
  varnugNA <- FALSE 
  zeronugget <- FALSE
  sillbounded <- !is.na(sill)
  if (sum(index[-c(MEAN,VARIANCE)])==0) {
    if (index[VARIANCE]) PARAM[VARIANCE] <- 1
    cov.matrix <- chol(matrix(.C("CovarianceMatrixNatSc",
                                 distances,
                                 lc,
                                 covnr,
                                 as.double(PARAM),
                                 lpar,
                                 dim,
                                 cov.matrix=double(lc * lc),
                                 as.integer(0),
                                 DUP=FALSE)$cov.matrix     
                              ,ncol=lc))
    stopifnot(all(diag(cov.matrix)>=0))
    cov.matrix<- chol2inv(cov.matrix)
    cc <- crossprod(rep(1,lc), cov.matrix)
    PARAM[MEAN] <- (cc %*% data) / sum(cc)
    if (!index[VARIANCE]) return(PARAM)
    else {
      MLEtargetV <- data-PARAM[MEAN] 
      PARAM[VARIANCE] <-
        sqrt((crossprod(MLEtargetV, cov.matrix) %*% MLEtargetV)/lc)
      return(PARAM)
    }
  }
  if (sillbounded) {
    ## only VARIANCE need to be optimised
    ## NUGGET = SILL - VARIANCE
    if (xor(is.na(PARAM[VARIANCE]),is.na(PARAM[NUGGET])))
      stop("Sill fixed. Then variance and nugget should be given or unknown simultaneously.")
    if (!is.na(PARAM[VARIANCE]) && (PARAM[VARIANCE]+PARAM[NUGGET]!=sill))
      stop("sill !=  variance + nugget") 
    index[NUGGET] <- FALSE;
    PARAM[NUGGET] <- 0; ## just to avoid trouble with NAs later on (*)
    variables[VARIANCE]<-sill/2
    upper[NUGGET] <- upper[VARIANCE]<-sill
    lower[VARIANCE] <- lower[NUGGET] <- 0
  } else {
    if (is.na(PARAM[VARIANCE])) {
      if (is.na(PARAM[NUGGET])) {
         ## both NUGGET and VARIANCE have to be optimised
         ## instead optimise the SILL (which can be done by
         ## a formula) and estimate the part of NUGGET
         ## (which will be stored in NUGGET) !
         ##  consequences:
         ## * estimtation space increased only by 1, not 2
         ## * real nugget: NUGGET * SILL
         ## *  variance  : (1-NUGGET) * SILL
         varnugNA <- TRUE
         lower[NUGGET] <- 0
         upper[NUGGET] <- 1
	 variables[NUGGET] <- 0.5
         index[VARIANCE] <- FALSE
         PARAM[VARIANCE] <- 0 ## does not have any meaning, just to get rid
         ##                       of NA since later (*) checked whether VARIANCE
         ##                        in [0,oo] 
      } else {       
	if (PARAM[NUGGET]==0) {
          ## here SILL==VARIANCE and therefore, variance
          ## can be estimated by a formula without increasing the dimension
          
          ## more than only the variance has to be estimated
          zeronugget <- TRUE
          index[VARIANCE] <- FALSE
          PARAM[VARIANCE] <- 0 ## dito
        } else {
          upper[VARIANCE] <- upperbound.var.factor * max(binned.variogram)
          lower[VARIANCE] <- (vardata-PARAM[NUGGET])/lowerbound.var.factor
          if (lower[VARIANCE]<lowerbound.sill) {
            if (PrintLevel>0)
              cat("low.var=",lower[VARIANCE]," low.sill",lowerbound.sill,"\n")
            warning("param[NUGGET] might not be given correctly.")
            lower[VARIANCE] <- lowerbound.sill
          }
        }
      }
    }
    else {
      if (PARAM[VARIANCE]==0.0) {
        
      }
      if (is.na(PARAM[NUGGET])) { ## and !is.na(param[VARIANCE])
        lower[NUGGET]<-(var(data)-PARAM[VARIANCE])/lowerbound.var.factor
        if (lower[NUGGET]<lowerbound.sill) {
          if (PrintLevel>0)
            cat("\nlow.nug=",lower[NUGGET]," low.sill",lowerbound.sill,"\n")
          warning("param[VARIANCE] might not be given correctly!")
          lower[NUGGET]<-lowerbound.sill
        }
        upper[NUGGET] <- upperbound.var.factor * max.bin.vario
      }
    }
  }  

  if (index[SCALE]) {
    lower[SCALE]<- mindistances / lowerbound.scale.LS.factor
    upper[SCALE]<- maxdistances * upperbound.scale.factor
    if (use.naturalscaling) {
      ##  11: exact or numeric
      ##  13: MLE or numeric
      scalingmethod <- as.integer(11)  ## check with RFsimu.h, NATSCALEMLE!
      ##          or   3 ?
      if (.C("GetNaturalScaling",covnr,as.double(variables),
             as.integer(length(variables)),
             scalingmethod,double(1),error=integer(1),DUP=FALSE)$error)
        {
          scalingmethod <- as.integer(0)
          if (PrintLevel>0) cat("No natural scaling.")
        }
    } else {
      scalingmethod <- as.integer(0)## check with RFsimu.h, NATSCALEMLE!
    }
  } else {
    scalingmethod <- as.integer(0)## check with RFsimu.h, NATSCALEMLE!
  }

  stopifnot(length(index)==length(variables))
  if (any((PARAM[!index]>upper[!index]) | (PARAM[!index]<lower[!index]))){
    ## (*) : this part makes trouble if index has been reset (TRUE->FALSE)
    ##       and PARAM is not changed accordingly
    stop("fixed parameters out of range")
  }  
  index[MEAN] <- FALSE ## directly estimated
  variab <- variables[index]
  
  LB  <- lower[index]
  UB  <- upper[index]
  if (PrintLevel>2) {print("LB,UB");print(cbind(LB,UB))}

  weights <- cbind(sqrt(bins:1 * binned.n), bins:1, sqrt(bins:1), sqrt(binned.n))

  show <- function(text,col,param) {
    if (PrintLevel>3) {
      cat(text,format(param,dig=3),"\n")
      if (PrintLevel>4) {
        plot(bin.centers,binned.variogram)
        model.values <- .C("CovarianceNatSc",bin.centers,bins,
                          covnr,as.double(param),lpar,
                          dim,model.values=double(bins),
                          scalingmethod,DUP=FALSE)$model.values
        lines(bin.centers,param[VARIANCE]-model.values,col=col)
      }
    }
  }
  
 
  CoVariate <- matrix(1,ncol=1,nrow=lc)
  ## note "global" variables MLEMIN, MLEPARAM
  MLEtarget<-function(variab) {
    if (PrintLevel>2) print(variab,dig=10)
    if (any((variab<LB) | (variab>UB))) {
      if (PrintLevel>1) cat("WARNING! forbidden variab values !\n")  
      penalty <- variab 
      variab<-pmax(LB,pmin(UB,variab)) 
      penalty <- sum(variab-penalty)^2
      res <- MLEtarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      if (PrintLevel>3) {
        cat("penalty",format(c(variab,penalty + res * (1+penalty))),"\n")
      }
      return(penalty + res * (1+penalty))
    }
    param <- PARAM
    param[index]<-variab  
    if (sillbounded) param[NUGGET]<-sill - param[VARIANCE]
    if (varnugNA) {
      param[VARIANCE] <- 1.0 - param[NUGGET]
    }
    if (zeronugget) param[VARIANCE] <- 1.0
    
    options(show.error.messages = FALSE)
    cov.matrix <- try(chol(matrix(.C("CovarianceMatrixNatSc",
                                     distances,
                                     lc,
                                     covnr,
                                     as.double(param),
                                     lpar,
                                     dim,
                                     cov.matrix=double(lc * lc),
                                     scalingmethod,
                                     DUP=FALSE)$cov.matrix     
                                 ,ncol=lc)))
    options(show.error.messages = TRUE)
    if (is.numeric(cov.matrix)) {
      if (any(diag(cov.matrix)<0)) {stop("chol det<0!")}
      logdet <- 2 * sum(log(diag(cov.matrix))) 
      cov.matrix<- chol2inv(cov.matrix)
    } else {
      return(1E300)
    }
    if (!is.finite(logdet)) logdet<-1E300 ## the best ?! 
    
    if (is.na.mean)  {
      dummy <- crossprod(CoVariate,cov.matrix)
      m <- solve(dummy %*% CoVariate, dummy %*% data)
      MLEtargetV  <- data - CoVariate %*% m
    }
    quadratic  <- crossprod(MLEtargetV, cov.matrix) %*% MLEtargetV
    
    if (varnugNA || zeronugget) {
      res <- logdet + lc * log(quadratic) 
    } else {
      res <- logdet + quadratic
    }   
    if (res<MLEMIN) {
      if (is.na.mean) param[MEAN] <- m
      if (sillbounded) param[NUGGET]<-sill-param[VARIANCE]
      if (varnugNA) {
        sill <- quadratic / lc
        param[VARIANCE] <- sill * (1.0 - param[NUGGET])
        param[NUGGET] <- sill * param[NUGGET]
      }
      if (zeronugget) {
        param[VARIANCE] <- quadratic / lc
      }
      assign("MLEMIN",res,envir=ENVIR)
      assign("MLEPARAM",param,envir=ENVIR)
    }
    if (PrintLevel>2) cat("result MLE",res,"\n")
    return(res)
  }
    
  ## Note
  ## creating "global" variables LSMIN, LSPARAM
  ## using WEIGHTS, BINNEDSQUARE
  LStarget <- function(variab) {
    if (PrintLevel>2) print(variab,dig=20)
    if (any((variab<LB) | (variab>UB))) {
      if (PrintLevel>1) cat("WARNING! forbidden variab values !\n")
      penalty <- variab 
      variab<-pmax(LB,pmin(UB,variab)) 
      penalty <- sum(variab-penalty)^2
      res <- LStarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      return(penalty + res * (1+penalty))
    }
    param <- PARAM
    param[index]<-variab
    if (sillbounded) param[NUGGET]<-sill-param[VARIANCE]
    if (varnugNA) param[VARIANCE] <- 1.0 - param[NUGGET]
    if (zeronugget) param[VARIANCE] <- 1.0
    
    ## if CoVariates
    ##binned.variogram <- .C("binnedvariogram",Xcoord,Ycoord,
    ##     as.double(data-f(coord,parameter)),   ###
    ##     lc,step,binned.variogram=double(bins),integer(bins),bins,
    ##     DUP=FALSE)$binned.variogram
    model.values <-
      .C("VariogramNatSc", bin.centers, bins, covnr, as.double(param), lpar,
         dim,model.values=double(bins), scalingmethod, DUP=FALSE)$model.values
    if (varnugNA || zeronugget) {
      bgw <- sum(binned.variogram * model.values * WEIGHTS)
      g2w <- sum(model.values^2 * WEIGHTS)
      ergb <- BINNEDSQUARE - bgw^2/g2w
    } else {
      ergb<-sum((binned.variogram - model.values)^2 * WEIGHTS) 
    }
    if (ergb<LSMIN) {
      if (sillbounded) param[NUGGET]<-sill-param[VARIANCE]
      if (varnugNA) {
        sill <- bgw/g2w
        param[VARIANCE] <- sill * (1.0 - param[NUGGET])
        param[NUGGET] <- sill * param[NUGGET]
      }
      if (zeronugget) {
        param[VARIANCE] <- bgw/g2w
      }
      assign("LSMIN",ergb,envir=ENVIR)
      assign("LSPARAM",param,envir=ENVIR)
      show("LSX ",1,LSPARAM)
    }
    return(ergb)
  }

  ## Note
  ## creating "global" variable WEIGHTS, LSMIN, BINNEDSQUARE, LSpMIN,LSpMAX
  ## using LSPARAM, LSMIN
  LSoptimize <- function(index,variab) {
    minimum <- Inf
    LSpmin <- rep(Inf,length(variab))
    LSpmax <- rep(-Inf,length(variab))
    for (W in 1:ncol(weights)) {
      assign("WEIGHTS",weights[,W],envir=ENVIR)
      if (varnugNA || zeronugget) {
        assign("BINNEDSQUARE",sum(binned.variogram^2*WEIGHTS),envir=ENVIR)
      }
      assign("LSMIN", Inf,envir=ENVIR)
      
      options(show.error.messages = FALSE) ##
      neuvariab <- try(optim(variab, LStarget, method ="L-BFGS-B", lower = LB,
                            upper = UB, control=list(trace=trace.optim))$par)
      ## side effect: minimum so far is in LSMIN and LSPARAM
      ## even if the algorithm finally fails
      options(show.error.messages = TRUE) ##
      
      if (LSMIN>1E299) next ### ==1E300 or Inf,,see LStarget,
      ##                         and assign(LSMIN) above

      neuvariab <- LSPARAM[index]
      LSpmin <- pmin(LSpmin,neuvariab)
      LSpmax <- pmax(LSpmax,neuvariab)

      mlet <- MLEtarget(neuvariab) ## necessary since LStarget is optimised!

      if (mlet<minimum) {
        minimum.param <- LSPARAM
        minimum.variab <- LSPARAM[index]
        minimum <- mlet
        minimum.index <- W## only for debugging reasons.
        ##                   Which LS weights have been the best?
      }
    }    
    if (!is.finite(minimum)) {stop("Minimum not found")}
    if (PrintLevel>3)
      cat("LS. . . ",format(minimum.param,dig=3),"(",format(minimum.index),")")
    return(param=minimum.param,variab=minimum.variab,minimum=minimum,
           LSpmin=LSpmin,LSpmax=LSpmax)
  }

  ## find a good initial value for MLE using  weighted least squares
  ## and binned variogram
  ##
  ## background: if the number of observations (and the observation
  ## field) tends to infinity then any least square algorithm should
  ## yield the same result as MLE
  ## so the hope is that for a finite number of points the least squares
  ## find an acceptable initial values  
  assign("MLEMIN",Inf,envir=ENVIR) ## necessary here, as LSoptim calls MLEtarget!
  cat(pch)
  LSopt <- LSoptimize(index,variab)
  PARAM <- LSopt$param    ## could use also LSPARAM directly...
  variab <- LSopt$variab
  LSpmin <- LSopt$LSpmin
  LSpmax <- LSopt$LSpmax
   
  ## if rather close to nugget and nugget==0 (fixed), do exception handling.
  ## This part is active only if
  ## scale.max.relative.factor < lowerbound.scale.LS.factor
  ## By default it is not, i.e., the following if-condition
  ## will "always" be FALSE.
  if ((PARAM[SCALE] < mindistances/scale.max.relative.factor) &&
      (!index[NUGGET]))
    {
      if (PrintLevel>1) cat("Looks like if a nugget effect is included...\n")
      if (PARAM[NUGGET]!=0) stop("Chosen model does not seem to be appropriate!")
      param[NUGGET] <- NA;
      return(mleRF(coord=coord, data=data, model=model, param=param,
                   lower.kappa=lower.kappa, upper.kappa=upperkappa, sill=sill,
                   use.naturalscaling=use.naturalscaling,
                   PrintLevel=PrintLevel, trace.optim=trace.optim,
                   bins=bin, distance.factor=distance.factor,
                   upperbound.scale.factor=upperbound.scale.factor,
                   lowerbound.scale.factor=lowerbound.scale.factor,
                   lowerbound.scale.LS.factor=lowerbound.scale.LS.factor,
                   upperbound.var.factor=upperbound.var.factor,
                   lowerbound.var.factor=lowerbound.var.factor,
                   lowerbound.sill=lowerbound.sill,
                   scale.max.relative.factor=scale.max.relative.factor,
                   minbounddistance=minbounddistance,
                   minboundreldist=minboundreldist,
                   approximate.functioncalls=approximate.functioncalls,
                   pch=pch
                   )
             )
    }

  show("Final",3,PARAM)
  
  ## now MLE
  assign("MLEMIN",LSopt$minimum,envir=ENVIR)
  assign("MLEPARAM",PARAM,envir=ENVIR)

  ## lowerbound.scale.LS.factor <  lowerbound.scale.factor, usually
  ## LS optimisation should not run to a boundary (what often happens
  ## for the scale) since a boundary value is usually a bad initial
  ## value for MLE (heuristic statement). Therefore a small
  ## lowerbound.scale.LS.factor is used for LS optimisation.
  ## For MLE estimation we should include the true value of the scale;
  ## so the bounds must be larger. Here lower[SCALE] is corrected
  ## to be suitable for MLE estimation
  lower[SCALE] <-
    lower[SCALE] * lowerbound.scale.LS.factor / lowerbound.scale.factor
  LB  <- lower[index]

  ##save(file="debug.save",MLEMIN,MLEPARAM,variab,LB,UB,trace.optim,covnr,index,
  ##     PrintLevel,param,sillbounded,varnugNA,zeronugget,MEAN,NUGGET,VARIANCE,
  ##     is.na.mean,SCALE,lc,distances,lpar,dim,scalingmethod,CoVariate,data)

  distances <- as.double(distances)
  options(show.error.messages = FALSE) ##

  cat(pch)
  variab <- try(optim(variab, MLEtarget, method="L-BFGS-B", lower = LB,
                      upper = UB, control=list(trace=trace.optim))$par)
  
  options(show.error.messages = TRUE) ##
  
  if (!is.numeric(variab) && (PrintLevel>0)) cat("MLEtarget I failed.\n")

  PARAM <- MLEPARAM
  variab <- PARAM[index]
 
  mindistance <- pmax(minbounddistance, minboundreldist * abs(variab))
  onborderline <- any(abs(variab-LB) <
                      pmax(mindistance,              ## absolute difference
                           minboundreldist * abs(LB) ## relative difference
                           )) ||
                  any(abs(variab-UB) <
                      pmax(mindistance, minboundreldist * abs(UB))) 
  ## if the MLE result is close to the border, it usually means that
  ## the algorithm has failed, especially because of a bad starting
  ## value (least squares do not always give a good starting point, helas)
  ## so the brutal method:
  ## calculate the MLE values on a grid and start the optimization with
  ## the best grid point. Again, there is the believe that the
  ## least square give at least a hint what a good grid is
  if ( onborderline ) {
    if (PrintLevel>3) cat("MLE a ", format(c(MLEPARAM,MLEMIN),dig=4),"\n")
    gridlength <- round(approximate.functioncalls^(1/sum(index)))
    ## grid is given by the extremes of the LS results
    ## so, therefore we should examine above at least 4 different sets
    ## of weights 
    step <- (LSpmax-LSpmin) / (gridlength-2) # grid starts bit outside

    LSpmax <- pmin(LSpmax + step/2,UB)       # the extremes of LS
    LSpmin <- pmax(LSpmin - step/2,LB)
    step <- (LSpmax-LSpmin) / (gridlength-1)
    i <- 1
    zk <-  paste("LSpmin[",i,"] + step[",i,"] * (0:",gridlength-1,")")
    if (length(step)>1)
      for (i in 2:length(step))
        zk <- paste(zk,",LSpmin[",i,"] + step[",i,"] * (0:",gridlength-1,")")
    zk <- paste("expand.grid(",zk,")")
    startingvalues <- eval(parse(text=zk))
    MLEMIN.old <- MLEMIN
    MLEPARAM.old <- MLEPARAM
    MLEMIN <- Inf
    cat(pch)
    apply(startingvalues,1,MLEtarget) ## side effect: Minimum is in
    ##                                   MLEMIN !
    PARAM <- MLEPARAM
    variab <- PARAM[index]
    options(show.error.messages = FALSE) ##

    cat(pch)
    variab <- try(optim(variab, MLEtarget,
                        method ="L-BFGS-B",
                        lower = LB, upper = UB,
                        control=list(trace=trace.optim))$par)
    options(show.error.messages = TRUE) ##
    if (!is.numeric(variab) && (PrintLevel>0)) cat("MLEtarget II failed.\n")
    ## do not check anymore whether there had been convergence or not.
    ## just take the best of the two strategies (initial value given by
    ## LS, initial value given by a grid), and be happy.
    if (MLEMIN<MLEMIN.old)  PARAM <- MLEPARAM
    else PARAM <- MLEPARAM.old
    variab <- PARAM[index]        
  }
  show("MLE",col=2,MLEPARAM)

  ## if the covariance functions use natural scaling, just
  ## correct the final output by GNS$natscale
  ## (GNS$natscale==1 if no rescaling was used)
  GNS <- .C("GetNaturalScaling",
            covnr,
            as.double(PARAM),
            as.integer(length(PARAM)),
            scalingmethod,
            natscale=double(1),
            error=integer(1),DUP=FALSE)
  if (GNS$error)
    stop(paste("Error",error,"occured whilst rescaling"))
  PARAM[SCALE] <- PARAM[SCALE] * GNS$natscale
  if (PrintLevel>2)
    cat("MLE(end)",format(c(PARAM,,MLEMIN),dig=4),"\n")
  if (pch!="") cat("\n")
  return(PARAM)
}


