

## PrintLevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information

fitvario <-
function(x, y=NULL, z=NULL, T=NULL, data, model,param,
         lower=NULL, upper=NULL, sill=NA,
         ...) {
  
  fitvario.default(x=x, y=y, z=z, T=T, data=data, model=model, param=param,
                lower=lower, upper=upper, sill=sill,
                ...)
}


mleRF <-
function(x, y=NULL, z=NULL, T=NULL, data, model,param,
         lower=NULL, upper=NULL, sill=NA,
         ...) {
  warning("mleRF is obsolete -- use fitvario instead")
  fitvario.default(x=x, y=y, z=z, T=T, data=data, model=model, param=param,
                lower=lower, upper=upper, sill=sill,
                ...)
}



fitvario.default <-
function(x, y=NULL, z=NULL, T=NULL, data, model, param,
         lower=NULL, upper=NULL, sill=NA,
         trend,
         use.naturalscaling=TRUE,
         PrintLevel=RFparameters()$Print, trace.optim=0,
         bins=20, nphi=1, ntheta=1, ntime=20,
         distance.factor=0.5,
         upperbound.scale.factor=10,
         lowerbound.scale.factor=20,
         lowerbound.scale.LS.factor=5,
         upperbound.var.factor=10,
         lowerbound.var.factor=100,
         lowerbound.sill=1E-10,
         scale.max.relative.factor=1000,
         minbounddistance=0.001, minboundreldist=0.02,
         approximate.functioncalls=50,
         pch="*", # no.mle=FALSE,
         var.name="X", time.name="T",
         transform=NULL, standard.style=NULL
) {
  
  ENVIR <- environment()
  save.transform <- transform
  save.bins <- bins;
  save.upper <- upper
  save.lower <- lower
  OLD.NUGGET.POSITION <- 3
                
######################################################################
###                function definitions                            ###
######################################################################

  ## function not used yet -- for trend definitions
  formulaToMatrix <- function(formula, coord, var.name="X", time.name="T") {
    coord <- as.matrix(coord)
    if (is.null(time.name)) co <- ""
    else co <- paste(time.name,"=coord[,ncol(coord)]")
    for (i in 1:ncol(coord))
      co <- paste(co, ",", var.name, i, "=coord[,", i, "]", sep="")
    eval(parse(text=paste("l <- list(",co,")")))
    fo <- as.list(formula)[[2]]
    variables <- NULL
    while (length(fo)==3) {
      dummy <- as.list(fo)
      stopifnot(dummy[[1]]=="+")
      fo <- dummy[[2]]
      variables <- cbind(eval(dummy[[3]],l), variables)
    }
    variables <- cbind(eval(fo, l), variables) 
    return(variables)
  }


  ## Note
  ## creating "global" variables LSMIN, LSPARAM
  ## using WEIGHTS
  LStarget <- function(variab) {
    if (PrintLevel>3) {cat(LEVEL); print(variab,dig=20)}
    if (any((variab<LB) | (variab>UB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 13.12.03 still happens ...
      if (PrintLevel>1) 
        cat("WARNING! forbidden variab values !",variab,"[",LB,",",UB,"]\n")
      penalty <- variab 
      variab<-pmax(LB,pmin(UB,variab)) 
      penalty <- sum(variab-penalty)^2
      res <- LStarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      return(penalty + res * (1+penalty)) ## not the best ....
    }
    param <- PARAM
    param[index] <- variab
    param <- transform(param)
    
    model.values <-
     .C("VariogramNatSc", bin.centers, bins, covnr, as.double(param), lpar,
        logicaldim, xdim,
        as.integer(length(covnr)), as.integer(pm$anisotropy),         
        as.integer(pm$op),model.values=double(bins), scalingmethod,
        PACKAGE="RandomFields", DUP=FALSE)$model.values

    if (any(is.na(model.values))) {
      if (PrintLevel>1) {
        cat("\nmissing values!")
        print(variab, dig=20)
        print(model.values)
        print("end model.values")
      }
      return(1E300)
    }
    
    if (SELF.WEIGHING) {
      ## weights = 1/ model.values^2
        gx <- binned.variogram / model.values
        gx <- gx[!is.na(gx)]
        if (varnugNA || zeronugget) {
          bgw <- sum(gx^2)
          g2w <- sum(gx)
          ergb <- length(gx) - g2w^2 / bgw
        } else {
          ergb <- sum((gx - 1)^2)
        }
    } else {
      if (varnugNA || zeronugget) {
        ## in this case the calculated variogram model is not the one
        ## to which we should compare, but:
        bgw <- sum(binned.variogram * model.values * WEIGHTS)
        g2w <- sum(model.values^2 * WEIGHTS)
        ergb <- BINNEDSQUARE - bgw^2/g2w
      } else ergb<-sum((binned.variogram - model.values)^2 * WEIGHTS)
    }
    
    if (ergb<LSMIN) {
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
    for (W in 0:ncol(weights)) {
      assign("LEVEL", paste("LSQ",W,":"), envir=ENVIR)
      assign("SELF.WEIGHING", W==0, envir=ENVIR)
      if (!SELF.WEIGHING) {
        assign("WEIGHTS",weights[, W],envir=ENVIR)
        if (varnugNA || zeronugget) {
          assign("BINNEDSQUARE",sum(binned.variogram^2*WEIGHTS),envir=ENVIR)
        }
      }
      assign("LSMIN", Inf, envir=ENVIR)
      
      options(show.error.messages = FALSE) ##
      if (length(variab)==1) {
        neuvariab <- try(optimize(LStarget, lower = LB, upper = UB)$minimum)
      } else {
        neuvariab <- try(optim(variab, LStarget, method ="L-BFGS-B", lower = LB,
                               upper = UB, control=optim.control)$par)
      }      
      ## side effect: minimum so far is in LSMIN and LSPARAM
      ## even if the algorithm finally fails
      options(show.error.messages = TRUE) ##

      if (W<=3) assign("lsqPARAM", cbind(lsqPARAM, LSPARAM), envir=ENVIR)## only
      ## used in return(...) -- first set of weights is the set of sqrt(nbin)
      
      if (LSMIN>1E299) next ### ==1E300 or Inf, see LStarget,
      ##                         and assign(LSMIN) above

      neuvariab <- LSPARAM[index]
      LSpmin <- pmin(LSpmin, neuvariab) ## used in case it is not clear
      ##      that later on the MLE was successful. Then ML is 
      ##      calculated on a grid, and best value for the grid is
      ##      used as starting value for another ML optimisation
      ##      LSpmin and LSpmax will give the approximate range of the grid
      LSpmax <- pmax(LSpmax, neuvariab)

      assign("LEVEL", paste("LSQ - MLE :"), envir=ENVIR)
      mlet <- MLEtarget(neuvariab) ## necessary since LStarget not MLE
      ##                              is optimised!

      if (mlet<minimum) {
        minimum.param <- LSPARAM
        minimum.variab <- LSPARAM[index]
        minimum <- mlet
        minimum.index <- W## only for debugging reasons.
        ##                   Which set of LS weights has been the best?
     }
    }    
    if (!is.finite(minimum)) {stop("Minimum not found")}
    if (PrintLevel>3)
      cat("LS. . . ",format(minimum.param,dig=3),"(",format(minimum.index),")")
    return(list(param=minimum.param, variab=minimum.variab, minimum=minimum,
                LSpmin=LSpmin, LSpmax=LSpmax))
  }

      
  ## note "global" variables MLEMIN, MLEPARAM
  MLEtarget<-function(variab) {
    if (PrintLevel>3) {cat(LEVEL); print(variab,dig=10)}
    if (any((variab<LB) | (variab>UB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 23.12.03 : still happens
      if (PrintLevel>0)
        cat("WARNING! forbidden variab values ! -- if there are too many warnings try narrower lower and upper bounds for the variables.",variab,"[",LB,",",UB,"]\n")
      penalty <- variab 
      variab <- pmax(LB,pmin(UB,variab)) 
      penalty <- sum(variab-penalty)^2 ## not the best ....
      res <- MLEtarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      if (PrintLevel>3) {
        cat("penalty",format(c(variab,penalty + res * (1+penalty))),"\n")
      }
      return(penalty + res * (1+penalty))
    }
    param <- PARAM
    param[index] <- variab
    param <- transform(param)
    options(show.error.messages = FALSE)
    cov.matrix <-
      try(chol(matrix(.C("CovarianceMatrixNatSc",
                         distances, lc,
                         covnr, as.double(param), lpar,
                         logicaldim, ##  space time dim of random field
                         xdim,       ##  dimension of the distance (vector)
                         ##          so for isotropic fields it might be 1
                         ##          although that the random field is 3
                         as.integer(length(covnr)),
                         as.integer(pm$anisotropy),
                         as.integer(pm$op),
                         cov.matrix=double(lc * lc),
                         scalingmethod,
                         PACKAGE="RandomFields", DUP=FALSE)$cov.matrix     
                      ,ncol=lc)))
    options(show.error.messages = TRUE)
    if (!is.numeric(cov.matrix) || any(is.na(cov.matrix))) {
      if (PrintLevel>1) {
        cat("\nerror in cholesky decomposition -- matrix pos def?")
        print(variab, dig=20)
        str(cov.matrix)
      }
      return(1E300)
    }    
    if (any(diag(cov.matrix)<0)) {stop("chol det<0!")}
    
    ## der faktor zwei erklaert sich wie folgt:
    ## logarithmus der wahrscheinlichkeitsdichte ergibt
    ## - lc/2 log(2pi) - 1/2 log(det C) - 1/2 (D-X m)^T C^{-1} (D-X m)
    ## nur wird statt obiges maximimiert :  (-2) * () minimiert ohne die 
    ##     konstante lc/2 log(2pi):
    ##     log(det C) + (D-X m)^T C^{-1} (D-X m),
    ## genauer, fuer wiederholte daten:
    ## repet * log(det C) + sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
    ##
    ## somit hat log(det(C)) vorfaktor 1
    ## nun erhaelt man aus chol(matrix)=wurzel der matrix,
    ## dass sum(log(diag(cov.matrix))) = 1/2 log(det C)
    ##
    ## repet is number of (independent) repeted measurement of data
    logdet <- 2 * sum(log(diag(cov.matrix))) * repet #  repet * log(det C)
    if (!is.finite(logdet)) logdet <- 1E300 ## the best ?! 
     cov.matrix<- chol2inv(cov.matrix) # La.chol2inv, LIN=TRUE
    if (givenCoVariates)  {
    ## m minimimiert, falls
      ##      sum_i C^{-1} (D_i - X m) = 0
      ##  <=> \sum_i (X^T C^{-1} D_i) = repet X^T C^-1 X m
      ##  <=> m = repet^{-1} (X^T C^-1 X)^{-1} X^T C^{-1} \sum_i D_i
      ##  <=> m =  (X^T C^-1 X)^{-1} X^T C^{-1} meandata
      XC <- crossprod(CoVariates, cov.matrix) # X^T C^-1
      m <- solve(XC %*% CoVariates, XC %*% meandata)
      MLEtargetV  <- data - as.double(CoVariates %*% m) # (sumD-X m)
    } ## otherwise MLEtargetV had been set by the calling function

    quadratic <- sum(MLEtargetV * (cov.matrix %*% MLEtargetV))
    ##               sum_i (D_i - Xm)^T C^{-1} (D_i - X m)
    if (varnugNA || zeronugget) {
      ## varnugNA : sill = s^2, var=p * s^2, nugget= (1-p) * s^2
      ##            and s^2 is treated as it was a variance parameter
      ##            and p is an internal parameter stored in NUGGET
      ##
      ## so define C_1 =  s^2 C
      ## Then
      ##    repet * log(det C) + sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##  = repet * log(det C_1) + 2 * repet * lc * log(s) + s^{-2} *
      ##                                      sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##
      ##  Hence the optimal s is given by
      ## s^2 = (2 * repet * lc)^{-1} * 2 * sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##     =  (lc *  repet)^{-1} * sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##
      ## plugging in gives
      ##     repet * log(det C_1) -  repet * lc * log(repet * lc) +
      ##     repet * lc * log(sum_i (D_i-Xm)^T C^{-1} (D_i-X m)) + repet * lc
      ## = logdet + loglcrepet + lcrepet * log(sum_i (D_i-Xm)^T C^{-1} (D_i-X m))
      res <- sum(logdet + loglcrepet + lcrepet * log(quadratic))
    } else {
      res <- sum(logdet + quadratic)
    }   
    if (res<MLEMIN) {
      if (givenCoVariates) assign("MLETREND", m, envir=ENVIR)
      if (varnugNA) {
        sill <- quadratic / lcrepet
        param[VARIANCE] <- sill * (1.0 - param[NUGGET])
        param[NUGGET] <- sill * param[NUGGET]
      } else
      if (zeronugget) param[VARIANCE] <- quadratic / lcrepet
      assign("MLEMIN",res,envir=ENVIR)
      assign("MLEPARAM",param,envir=ENVIR)
    }
    if (PrintLevel>3) cat("result MLE",res,"\n")
    return(res)
  }

  
######################################################################
###                Initial settings                                ###
######################################################################
  if (PrintLevel>2) cat("\ninitial settings...")
  
  new <- CheckXT(x, y, z, T, grid=FALSE)
  new$x <- as.matrix(new$x)

  
###########################      transform     #######################
  ## in case functions are given within model, they are extracted and
  ## replaced NaN; `transform' is set instead
  ## must be the first one, since model cannot be Prepare-d
  if (is.list(model) && is.list(unlist(model))) {
    if (!any(sapply(unlist(model), is.function)))
      stop("model does not seem to be defined correctly")
    if (!is.null(transform)) 
      stop("if functions are given in `model', `transform' must be NULL") 
    transform <- function(param) {
      p <- pm
      p$param <- param
      p <- convert.to.readable(p, allow="list")$model
      for (i in 1:length(p.transform))
        param[p.transform[i]] <- f.transform[[i]](p)
      return(param)
    }
    functionToNaN <- function(p) {
      ## replaces functions by NaN
      if (is.list(p)) for (i in 1:length(p)) p[[i]] <- functionToNaN(p[[i]])
      else if (is.function(p)) {
        assign("f.transform", c(f.transform, p), envir=ENVIR)
        return(NaN)
      }
      else if (is.nan(p))
        stop("NaN are not allowed if functions are given in `model'")
      return(p)
    }
    functionNum <- function(p) {
      ## replaces functions by an enumeration
      if (is.list(p)) for (i in 1:length(p)) p[[i]] <- functionNum(p[[i]])
      else if (is.function(p)) {
        assign("zaehler", zaehler + 1, envir=ENVIR)
        return(zaehler)
      }
      return(if (is.numeric(p)) rep(100000, length(p)) else p) # 100000 > anzahl
      ##                                                             parameter
    }
    ## aniso could have been a list, this is corrected by the following command
    f.transform <- NULL
    model.nan <- functionToNaN(model)
    model.nan <-
      convert.to.readable(PrepareModel(model.nan,NULL, # new$spacedim+new$Time
                                       ), allow="list")
    if (any(sapply(model.nan[-1], function(x) !is.null(x) && is.nan(x))))
      stop("transform functions may not given (yet) for mean and trend")
    zaehler <- 0
    fct.num <- PrepareModel(functionNum(model))$param
    p.transform <- PrepareModel(parampositions(model.nan, print=FALSE))$param
    stopifnot(length(fct.num)==length(p.transform)) ## if so, programming error!
    p.transform <-
      p.transform[order(fct.num)][1:sum(fct.num <= length(p.transform),
                                        na.rm=TRUE)]
    model <- model.nan
    rm("model.nan")
  }


############################    model       ##########################
  stopifnot(all(is.finite(as.matrix(data))))
  data <- as.vector(data)
  pm <- PrepareModel(model, param, new$spacedim+new$Time, trend)
  ctr <- convert.to.readable(pm)
  if (is.null(standard.style))
    standard.style <- is.vector(ctr$param) && is.null(transform)
  else if (!is.vector(ctr$param) && standard.style) {
    standard.style <- FALSE
    warning("standard.style must be FALSE for the given model specification.")
  }
  if (!is.null(ctr$param))
    ## to make sure that the model is  given in a standardised way
    ## if not truely a composed covariance function
    pm <- PrepareModel(ctr$model, ctr$param, new$spacedim+new$Time, trend)


##############              Coordinates               #################
  if (PrintLevel>2) cat("\ncoordinates...")
  spacedim <- ncol(new$x) 
  logicaldim <- as.integer(spacedim + !is.null(T))
  if (pm$anisotropy) {
    xdim <- as.integer(logicaldim)
    dd <- 0
    
    coord <- new$x
    if (!is.null(T)) {
      y <- expand.grid(1:nrow(new$x), seq(new$T[1], new$T[2], new$T[3]))
      coord <- cbind(coord[y[,1], ], y[,2])
    }
    
    lc <- nrow(coord)
    ## note: the direct C call needs matrix where points are given column-wise
    ##       whereas the R function CovarianceFct need them row-wise,
    ##                   except for fctcall==CovarianceMatrix
    distances <- matrix(nrow=logicaldim, ncol=lc * (lc-1) / 2)
    ## distances are used to calculate bounds for the scale parameter(s)
    ## to be estimated -- further it is used in MLEtarget

    for (i in 1:logicaldim) {
      dx <- outer(coord[,i], coord[,i], "-")
      distances[i, ] <- as.double(dx[lower.tri(dx)])
      ## x <- c(1,6,0.1)
      ## outer(x, x, "-")
      ##     [,1] [,2] [,3]
      ##[1,]  0.0 -5.0  0.9
      ##[2,]  5.0  0.0  5.9
      ##[3,] -0.9 -5.9  0.0
      if (i<=spacedim) dd <- dd + distances[i ,]^2
    }
    dd <- sqrt(dd)
    coord <- NULL
   
    maxdistances <- max(dd)
    mindistances <- min(dd[dd!=0])

    rm("dd")
  } else {
    lc <- nrow(new$x)
    xdim <- as.integer(1)
    distances <- as.double(dist(new$x))
    mindistances <- min(distances[distances!=0])
    maxdistances <- max(distances)
  }
  gc() ## necessairy, since otherwise too much memory is being taken because
  ##      of dd[dd!=0]
  
  storage.mode(lc) <- "integer"


##########   Check consistency of NA, NaN and transform  ##########
  if (PrintLevel>2) cat("\nconsistency of NA, NaN, transform...")
  PARAM <- pm$param
  ##PARAM <- param: just to make clear in MLEtarget and LStarget what is global,
  ##and do not overwrite param
  index <- is.na(PARAM) & !is.nan(PARAM)
  if (!is.null(transform)) {
    ## transform <- c("para[1] <- para[2]", "para[3] <- 4 - para[1]")
    pp <- para <- runif(length(PARAM))
    para[is.nan(PARAM)] <- NaN
    para <- transform(para)
    if (length(para)!=length(pp))
      stop("transform does not preserve length of the parameter vector")
    txt <- ""
    if (any(is.nan(para)))
     txt <-
       paste("position(s)", paste(which(is.nan(para)), collapse=","),
             "of param still NaN after application of 'transform'.\n'transform'",
             "is too complex or does not fill the positions with value NaN.\n\n")
    if (any(pp <- para[!is.nan(PARAM)]!=pp[!is.nan(PARAM)])) {
      pos <- which(!is.nan(PARAM))
      if (options()$warn<2) {
        if (any(pq <- is.na(PARAM)[!is.nan(PARAM)] & pp))
          warning(paste("position(s)", paste(pos[pq], collapse=","),
                        "of param has been changed by 'transform'",
                        "that were NA.\n"))
        pp <- pp & !pq
      }
      if (any(pp))
        txt <- paste(txt,"position(s)", paste(pos[pp], collapse=","),
                     "of param has been changed although not being NaN.\n\n")
    }
    if (txt!="") stop(txt)
  } 

  
##############      find upper and lower bounds      #################  
  if (PrintLevel>2) cat("\nbounds...")
  UnList <- function(l) {
    l <- unlist(l)
    paste(names(l), l, sep="=", collapse="|")
  }
  txt <- "lower and upper are both lists or vectors of the same length or NULL"
  lengthmismatch <- "lengths of bound vectors do not match model"
  structuremismatch <- "structures of bounds do not match the model"
  parscale <- start <- convert.to.readable(pm, allowed="list")$model
  ## parscale gives the magnitude of the parameters to be estimated
  ##     passed to optim/optimise so that the optimiser estimates
  ##     values around 1
  ## start will give the starting values for LSQ

  ## press lower and upper into a standard format...
  if (is.list(lower)) {
    if (!is.list(upper)) stop(txt)
    test <- pm
    test$param <- rep(NA, length(test$param))
    test <- Unlist(convert.to.readable(test, allowed="list"))
    if (!is.numeric(lower[[i]]$v)) lower[[i]]$var <- NA ##i.e. not given
    if (pm$anisotropy)
      if (is.null(lower[[i]]$a)) lower[[i]]$aniso <- start$a * NA
      else if (!is.numeric(lower[[i]]$s)) lower[[i]]$scale <- NA
    if (!is.numeric(lower[[i]]$k)) lower[[i]]$kappas <- start$param * NA
    bndtst <- PrepareModel(lower, time=new$spacedim+new$Time, trend=trend)
    bndtst$param <- rep(NA, length(bndtst$param))
    if (test!=Unlist(convert.to.readable(bndtst, allowed="list")))
      stop(structuremismatch)
    lower <- convert.to.readable(bndtst, allowed="list")

    if (!is.numeric(upper[[i]]$v)) upper[[i]]$var <- NA
    if (pm$anisotropy)
      if (is.null(upper[[i]]$a)) upper[[i]]$aniso <- start$a * NA
      else if (!is.numeric(upper[[i]]$s)) upper[[i]]$scale <- NA
    if (!is.numeric(upper[[i]]$k)) upper[[i]]$kappas <- start$param * NA
    bndtst <- PrepareModel(upper, time=new$spacedim+new$Time, trend=trend)
    bndtst$param <- rep(NA, length(bndtst$param))
    if (test!=Unlist(convert.to.readable(bndtst, allowed="list")))
      stop(structuremismatch)
    upper <- convert.to.readable(bndtst, allowed="list")
    rm("test")
  } else if (is.vector(lower)) {
    if (!is.vector(upper) || length(upper)!=length(lower)) stop(txt)
    if (is.vector(param)) {
      nulltrend <- as.integer(missing(trend) || is.null(trend))
      if (length(lower) < length(param)) {
        if (length(param) - 3 - nulltrend != length(lower))
          stop("length of lower does not match param")
        nugg <- if (is.finite(param[OLD.NUGGET.POSITION]) &&
                    param[OLD.NUGGET.POSITION]==0) 0 else NA
        lower <- c(rep(NA, 1 + nulltrend), nugg, NA, lower)
        upper <- c(rep(NA, 1 + nulltrend), nugg, NA, upper)
      } else if (length(param) < length(lower)) stop(lengthmismatch)      
      lower <- convert.to.readable(PrepareModel(model, lower,
                                                new$spacedim+new$Time, trend),
                                   allowed="list")
      upper <- convert.to.readable(PrepareModel(model, upper,
                                                new$spacedim+new$Time, trend),
                                   allowed="list")
    } else if (is.matrix(param)) {
      if (nrow(lower)<length(param)) {
        if (length(param) != 2 + nrow(lower)) stop(lengthmismatch)
        lower <- c(NA, NA, lower)
        upper <- c(NA, NA, upper)
      } else if (nrow(lower)>length(param)) stop(lengthmismatch)     
      lower <- convert.to.readable(PrepareModel(model,
                matrix(lower, ncol=ncol(param), nrow=nrow(param)),
                                                new$spacedim+new$Time, trend),
                                                allowed="list")
      upper <- convert.to.readable(PrepareModel(model,
               matrix(upper, ncol=ncol(param), nrow=nrow(param)),
                                                new$spacedim+new$Time, trend),
                                   allowed="list")
    } else stop("vector valued lower bound only allowed if param is given")
  } else { # !(is.vector(lower)) i.e. not given
    if (!is.null(upper) || !is.null(lower)) stop(txt)
    lower <- pm
    lower$param <- rep(NA, length(lower$param))
    lower <- upper <- convert.to.readable(lower, allowed="list")
  }
  lower <- lower$model
  upper <- upper$model

### now the bounds and starting values are set...
  scale.pos <- pm
  scale.pos$param <- rep(0, length(scale.pos$param))
  scale.pos <- convert.to.readable(scale.pos, allowed="list")$model
  ## scale.pos will give the positions where scale parameters can be found
  ##    this is necessary to change the limits for the scale parameters
  ##    before MLE is performed
  vardata <- var(data)
  truedim <- new$spacedim+new$Time ## could made better by considering
  ##   the anisotropy matrices if there are any
  ##   but too much work with the low profit too allow some more
  ##   models for the estimation (those that are allowed for
  ##   lower dimensions only, and the anisotropy matrix has a rank
  ##   small enough

  for (i in seq(1, length(lower), 2)) {   
    if (is.na(lower[[i]]$var))
      ## lower bound of first model is treated differently!
      ## so the "main" model should be given first!
      ## usually "main" model + (estimated) nugget
      lower[[i]]$var <- if (i==1) vardata / lowerbound.var.factor else 0
    if (is.na(upper[[i]]$var)) upper[[i]]$var <- upperbound.var.factor *vardata
    if (is.na(start[[i]]$var))
      parscale[[i]]$var <- start[[i]]$var <- vardata / (length(lower) + 1) *2
    if (parscale[[i]]$var==0) parscale[[i]]$var <- 1
    if (pm$anisotropy) {
      scale.pos[[i]]$aniso <- rep((i+1)/2, length(scale.pos[[i]]$aniso))      
      diagonal <- as.vector(diag(logicaldim))
      if (upper[[i]]$model=="nugget") {
        ## user may have given an anisotropy matrix
        ## matrix entries are supposed to be within [-1,1], so
        ## -10, 10 as bounds is sufficient
        upper[[i]]$aniso[is.na(upper[[i]]$aniso)] <- 10
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & diagonal] <- 0
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & !diagonal] <- -10
      } else {
        upper[[i]]$aniso[is.na(upper[[i]]$aniso)] <-
          lowerbound.scale.LS.factor / mindistances
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & diagonal] <-
          1 / (upperbound.scale.factor * maxdistances)
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & !diagonal] <-
          -lowerbound.scale.LS.factor / mindistances
      }
      idx <- is.na(start[[i]]$aniso)
      parscale[[i]]$aniso[idx] <- 8 / (maxdistances + 7 * mindistances)
      start[[i]]$aniso[idx & diagonal] <- 8 / (maxdistances + 7 * mindistances)
      start[[i]]$aniso[idx & !diagonal] <- 0      
    } else { ## isotropy
      scale.pos[[i]]$scale <- (i+1)/2
      if (upper[[i]]$model=="nugget") {
        ## user should never have chosen the value for the scale 
        lower[[i]]$scale <- 0
        upper[[i]]$scale <- 1 
      } else {
        if (is.na(lower[[i]]$scale))
          lower[[i]]$scale <- mindistances / lowerbound.scale.LS.factor
        if (is.na(upper[[i]]$scale))
          upper[[i]]$scale <- maxdistances * upperbound.scale.factor
      }
      if (is.na(start[[i]]$scale))
        parscale[[i]]$scale <- start[[i]]$scale <-
          (maxdistances + 7 * mindistances) / 8
    }
    kappas <- parameter.range(start[[i]]$model, truedim)
    if (!is.null(kappas)) {
      if (length(kappas)==1 && is.nan(kappas)) stop(paste("model #", (i+1)/2,
                 "is not allowed for the considered dimension"))
      fixed.kappas <- kappas$th[[1]][2,] == kappas$th[[1]][1,]
      if (any(fixed.kappas)) {
        if (!all(is.finite(start[[i]]$kappa[fixed.kappas])))
          stop(paste("in model #", (i+1)/2,
                     "the parameters that select subclasses are not all given"))
        for (nk in 1:length(kappas$th)) 
          if (all(kappas$th[[nk]][fixed.kappas] ==
                  start[[nk]]$kappa[fixed.kappas])) break
        if (any(kappas$th[[nk]][fixed.kappas]!=start[[nk]]$kappa[fixed.kappas]))
          stop(paste("Could not find the indicated subclass for model #",
                     (i+1)/2))
      } else nk <- 1
      lower[[i]]$kappas[is.na(lower[[i]]$kappas)] <-
        kappas$pr[[nk]][1, is.na(lower[[i]]$kappas)]
      upper[[i]]$kappas[is.na(upper[[i]]$kappas)] <-
        kappas$pr[[nk]][2, is.na(upper[[i]]$kappas)] 
      start[[i]]$kappas <- (upper[[i]]$kappas + lower[[i]]$kappas) / 2
      parscale[[i]]$kappas <- rep(1, length(start[[i]]$kappas))
    } # !is.null(kappas)
  }

  lower <- PrepareModel(model=lower, time=new$spacedim+new$Time, trend=trend,
                        named=PrintLevel>1)$par
  upper <- PrepareModel(model=upper, time=new$spacedim+new$Time, trend=trend,
                        named=PrintLevel>1)$par
  lower[is.nan(PARAM)] <- -Inf ## nan should not be checked -- this is
  ##                              completely up to the user ...
  upper[is.nan(PARAM)] <- Inf
  
  start <- PrepareModel(model=start, time=new$spacedim+new$Time, trend=trend)$par
  parscale <- 
    PrepareModel(model=parscale, time=new$spacedim+new$Time, trend=trend)$par
  scale.pos <-
    PrepareModel(model=scale.pos, time=new$spacedim+new$Time, trend=trend)$par


##########################       data         ########################
  ## coord instead of new$x because of possible time component

  if (length(data)==lc) { # data is as.double(data)
    repet <- 1
    meandata <- data
  } else {
    if (as.integer(length(data) / lc) * lc != length(data))
      stop("length of data not a multiple of number of locations")
    data <- matrix(data, nrow=lc)
    repet <- ncol(data)
    meandata <- (data %*% rep(1, repet)) / repet
  }
  lcrepet <- lc * repet
  loglcrepet <- lcrepet * (1 - log(lcrepet))
  
  
##############           Covariates                   #################
  if (PrintLevel>2) cat("\ncovariates...")
  givenCoVariates <- TRUE
  if (is.null(pm$trend)) {
    if (is.na(pm$mean)) {
      CoVariates <- matrix(1, nrow=nrow(new$x), ncol=1)
      stopifnot(!is.null(CoVariates)) 
    }
    else {
      givenCoVariates <- FALSE 
      MLEtargetV <- data - pm$mean
    }
  } else {
    stop("sorry. not available yet") ## programmed but not tested yet
    if (is.matrix(trend)) {
      stopifnot(nrow(trend)==nrow(new$x))
      CoVariates <- trend
    } else { ## formula
      CoVariates <- formulaToMatrix(trend, cbind(new$x,new$T),
                                 var.name=var.name, time.name=time.name)
    }
  }
  
  if (!any(index)) {
    ## nothing to be estimated (except maybe for the mean)
    C <-  CovarianceFct(fct="CovarianceMatrix", model=ctr,
                         if (is.vector(distances)) distances else t(distances))
    logeigen <- sum(log(eigen(C, symm=TRUE)$values))
    InvC <- solve(C)
    if (givenCoVariates) {
      if (!is.null(pm$trend)) stop("sorry, not implemented yet")
      ## m =  (X^T C^-1 X)^{-1} X^T C^{-1} meandata
      XC <- crossprod(CoVariates, InvC)  
      m <- solve(XC %*% CoVariates, XC %*% meandata)
      pm$mean <- m
      ctr <- convert.to.readable(pm)
      MLEtargetV <- data - as.double(CoVariates %*% m)
    }
    loglikelihood <- -0.5 * (repet * (lc * log(2 * pi) + logeigen) +
                             sum(MLEtargetV * (InvC %*% MLEtargetV))
                             )
    return(list(mle.value=loglikelihood, trend.coff=NULL,
                mle=ctr, ev  = ctr, lsq = ctr, nlsq= ctr, slsq= ctr, flsq=ctr))
  } 

  
##############          Empirical Variogram           #################
  if (PrintLevel>2) cat("\nempirical variogram...")
  if (givenCoVariates) { 
    regr <- lsfit(CoVariates, data, intercept=FALSE)
    TREND <- regr$coeff
    EVtargetV <- regr$residuals
  } else EVtargetV <- data

  if (length(nphi)==1) nphi <- c(0, nphi) # starting angle; lines per half circle
  if (length(ntheta)==1) ntheta <- c(0, ntheta) # see above
  if (length(ntime)==1) ntime <- c(ntime, 1) 
  ntime <- ntime * T[3] ## time endpoint; step
  
  ev <- EmpiricalVariogram(new$x, T=new$T, data=EVtargetV, grid=FALSE, 
                           bin=if (length(bins)>1) bins else 
                           c(-1, seq(0, distance.factor * maxdistances,
                                    len=bins+1)),
                           phi=if ((new$spacedim>=2) && pm$anisotropy) nphi,
                           theta=if ((new$spacedim>=3) && pm$anisotropy) ntheta,
                           deltaT=if (!is.null(T)) ntime
                           )
  
  index.bv <- as.vector(!is.na(ev$e)) ## exclude bins without entry
  if (sum(index.bv)==0)
    stop("only NAs in empirical variogram; check values of bins and distance.factor")
  binned.variogram <- as.double(ev$e[index.bv])
 
  bin.centers <- as.matrix(ev$c)
  if (pm$anisotropy) {
    ## complete the coordinates of bin.centers to a vector of the dimension
    ## considered
    if (!is.null(ev$phi)) {
      if (spacedim<2) stop("x dimension is less than two, but phi is given") 
      bin.centers <- cbind(as.vector(outer(bin.centers, cos(ev$phi))),
                           as.vector(outer(bin.centers, sin(ev$phi))))
    }
    if (!is.null(ev$theta)) {
      if (spacedim<3)
        stop("x dimension is less than three, but theta is given") 
      if (ncol(bin.centers)==1) bin.centers <- cbind(bin.centers, 0)
      bin.centers <- cbind(as.vector(outer(bin.centers[, 1], cos(ev$theta))),
                           as.vector(outer(bin.centers[, 2], cos(ev$theta))),
                           rep(sin(ev$theta), each=nrow(bin.centers)))
    } else {
      if (nrow(bin.centers) < spacedim) # dimension of bincenter vector
        ##                       smaller than dimension of location space
        bin.centers <- 
          cbind(bin.centers, matrix(0, nrow=nrow(bin.centers),
                                    ncol=spacedim - ncol(bin.centers)
                                    ))
    }
    if (!is.null(ev$T)) {
      bin.centers <-
        cbind(matrix(rep(t(bin.centers), length(ev$T)), byrow=TRUE,
                     ncol = ncol(bin.centers)),
              rep(ev$T, each=nrow(bin.centers)))      
    }
  }
 
  bin.centers  <- as.double(t(bin.centers[index.bv, ])) #
  ##  es muessen beim direkten C-aufruf die componenten der Punkte
  ##  hintereinander kommen (siehe auch variable distance). Deshalb t()

  evsd <- as.double(ev$sd)
  evsd[is.na(evsd)] <- 0
  evsd[evsd==0] <- 10 * sum(evsd, na.rm=TRUE) ## == "infinity"
  
  bins             <- length(ev$n)
  binned.n         <- as.integer(ev$n)
  weights <- cbind(rep(1, bins), ## must be first, since lsqPARAM
                   ##                           is returned
                   sqrt(binned.n), ## must be 2nd, since lsqPARAM is returned
                   1 / evsd,## must be 2nd, since lsqPARAM is returned
                   ## see return and W<=3 in LSoptim if further lsq are added
                   sqrt(bins:1 * binned.n), bins:1, sqrt(bins:1)
                   )[index.bv, ]  
  bins <- as.integer(sum(index.bv))
  EVtargetV <- NULL
 
######################################################################
###   check which parameters are NA -- only old style is allowed   ###
###     certain combinations of NA allow for a much faster MLE     ###
######################################################################
  covnr <- as.integer(pm$covnr)
  lpar <- as.integer(length(PARAM))  
  varnugNA <- zeronugget <- sillbounded <- FALSE
  scalingmethod <- as.integer(0)
  if (standard.style) {
    if (PrintLevel>2) cat("\nstandard style settings...")
    VARIANCE <- 1
    KAPPA <- 2
    n.kappas <- .C("GetNrParameters", covnr, as.integer(1), k=integer(1),
                   PACKAGE="RandomFields", DUP=FALSE)$k
    SCALE  <- KAPPA + n.kappas
    if (SCALE<length(PARAM)) {
      NUGGET <- SCALE + 1
      nugget <- PARAM[NUGGET]
    } else {
      NUGGET <- 0  ## nugget effect is zero, so not included in PARAM
      nugget <- 0.0
    }

    if ((covnr==.C("GetModelNr", as.character("nugget"), as.integer(1),
           nr=integer(1), PACKAGE="RandomFields")$nr) &&
        (any(index[-NUGGET]) || (index[VARIANCE]!=0.0))) #works also if NUGG==0??
      stop("if model=nugget effect then the variance must be zero, and only mean and/or nugget can be NA")   
              
    if (sillbounded <- !is.na(sill)) {
      ## only VARIANCE need to be optimised
      ## NUGGET = SILL - VARIANCE
      if (xor(is.na(PARAM[VARIANCE]), is.na(nugget)))
        stop("Sill fixed. Then variance and nugget should be given or unknown simultaneously.")
      if ((is.nan(PARAM[VARIANCE]) || is.nan(nugget)) && !is.null(transform))
        stop("not sure what to do now: sill fixed, but some transformation is given -- try standard.style=FALSE")
      if (!is.na(PARAM[VARIANCE]) && (PARAM[VARIANCE] + nugget!=sill))
        stop("sill !=  variance + nugget") 
      index[NUGGET] <- FALSE;
      PARAM[NUGGET] <- 0; ## just to avoid trouble with NAs later on (*)
      start[VARIANCE] <- sill/2
      upper[NUGGET] <- upper[VARIANCE] <- sill
      lower[NUGGET] <- lower[VARIANCE] <-  0
    } else { ## not sill bounded
      if (is.na(PARAM[VARIANCE]) &&
          (!is.nan(PARAM[VARIANCE]) || is.null(transform))) {
        if (is.na(nugget) && (!is.nan(nugget) || is.null(transform))) {
          ## both NUGGET and VARIANCE have to be optimised
          ## instead optimise the SILL (which can be done by
          ## a formula) and estimate the part of NUGGET
          ## (which will be stored in NUGGET) !
          ##  consequences:
          ## * estimtation space increased only by 1, not 2
          ## * real nugget: NUGGET * SILL
          ## *  variance  : (1-NUGGET) * SILL
          varnugNA <- TRUE
          upper[NUGGET] <- 1
          start[NUGGET] <- 0.5
          index[VARIANCE] <- FALSE
          lower[VARIANCE] <- 0
          PARAM[VARIANCE] <- 0 ## does not have any meaning, just to getrid
          ##                       of NA since later (*) checked whether VARIANCE
          ##                        in [0,oo] 
        } else { ## not sillbounded, is.na(variance), !is.na(nugget)
          if (nugget==0) {
            ## here SILL==VARIANCE and therefore, variance
            ## can be estimated by a formula without increasing the dimension
            
            ## more than only the variance has to be estimated
            zeronugget <- TRUE
            if (sum(!index)>1) index[VARIANCE] <- FALSE ## otherwise we get
            ##           again a special case to treat
            lower[VARIANCE] <- PARAM[VARIANCE] <- 0 ## dito
          } else { ## not sillbounded, is.na(variance), nugget!=0
            lower[VARIANCE] <- (vardata-nugget)/lowerbound.var.factor
            if (lower[VARIANCE]<lowerbound.sill) {
              if (PrintLevel>0)
                cat("low.var=",lower[VARIANCE]," low.sill",lowerbound.sill,"\n")
              warning("param[NUGGET] might not be given correctly.")
              lower[VARIANCE] <- lowerbound.sill
            }
          }
        }
      } else { ## not sillbounded, !is.na(variance)
        if (PARAM[VARIANCE]==0.0) {
          if (any(is.na(PARAM[KAPPA:SCALE])))
            stop("If variance=0, estimating scale or model parameters does not make sense")
          lower[VARIANCE] <- 0
        }
        if (is.na(nugget)) { ## and !is.na(param[VARIANCE])
          lower[NUGGET] <-
            pmax(0, (vardata - PARAM[VARIANCE]) / lowerbound.var.factor)
          if (lower[NUGGET] < lowerbound.sill) {
            if (PrintLevel>0)
              cat("\nlower nugget bound=", lower[NUGGET],
                  " < lower sill bound=", lowerbound.sill,
                  " -- is the variance given correctly?\n",sep="")
            warning("Has param[VARIANCE] been given correctly?!")
            lower[NUGGET] <- lowerbound.sill
          }
          start[NUGGET] <- sqrt(lower[NUGGET] * upper[NUGGET])
        }
      } ##  else { ## not sillbounded, !is.na(variance)
    } ## else { ## not sill bounded
    
    if (index[SCALE]  && use.naturalscaling) {
      ##  11: exact or numeric
      ##  13: MLE or numeric
      scalingmethod <- as.integer(11)  ## check with RFsimu.h, NATSCALEMLE!
      ##          or   3 ?
      if (.C("GetNaturalScaling", covnr, as.double(start[-1:-(KAPPA-1)]),
             scalingmethod,double(1), error=integer(1),
             PACKAGE="RandomFields", DUP=FALSE)$error)
        {
          scalingmethod <- as.integer(0)
          if (PrintLevel>0) cat("No natural scaling.")
        }
    }
    
    stopifnot(length(index)==length(start)) ## simple check
    notidx <- !index & !is.nan(PARAM)
    ## (*)
    if (any((PARAM[notidx]>upper[notidx]) | (PARAM[notidx]<lower[notidx]))){
      if (PrintLevel>1) {cat("\n");print(rbind(notidx, upper, PARAM, lower))}
      warning("fixed parameters out of range\n")
    }     
    stopifnot( varnugNA + zeronugget + sillbounded <= 1)
    trafo <- NULL;
    if (varnugNA)
      trafo <- function(param) {param[VARIANCE] <- 1.0 - param[NUGGET]; param}
    if (sillbounded)
      trafo <- function(param) {param[NUGGET] <- sill-param[VARIANCE]; param}
    if (zeronugget)
      trafo <- function(param) {param[VARIANCE] <- 1.0; param}
    if (is.null(transform)) transform <- trafo
    else if (!is.null(trafo)) {
      warning("standard.style and transform!=NULL may cause strange effects -- internal transformation is performed before the user's one")
      trafoUser <- transform
      transform <- function(param) trafoUser(trafo(param))
    }
  } ## standard.style
  if (is.null(transform)) transform <- function(x) x

 
 
######################################################################
###                     Estimation part itself                     ###
######################################################################

  
###################  LSQ  ########################
  if (PrintLevel>2) cat("\nLSQ optimisation...") else cat(pch)
  variab <- start[index]
  parscale <- parscale[index]
  stopifnot(all(is.finite(parscale)))
  optim.control <- list(trace=trace.optim, parscale=parscale)
  
  # do not place after  LSoptimize(index,variab)!!
  # since LSoptimize may overwrite it
  assign("lsqPARAM", NULL, envir=ENVIR) ## only used for returning the values
  ##                              of LSQ fits by fitvario
  ##                              lsqPARAM is a table
  assign("MLETREND", NULL, envir=ENVIR) ## just to avoid errors due to
  ##                         not being defined in case trend is not estimated    

  LB  <- lower[index]
  UB  <- upper[index]
  if (PrintLevel>1) {
    cat("lower:", lower, "\nupper:", upper, "\n")
    if (PrintLevel>2) {print("LB,UB");print(cbind(LB,UB))}
  }
  
    
  ## find a good initial value for MLE using weighted least squares
  ## and binned variogram
  ##
  ## background: if the number of observations (and the observation
  ## field) tends to infinity then any least square algorithm should
  ## yield the same result as MLE
  ## so the hope is that for a finite number of points the least squares
  ## find an acceptable initial values
 
  distances <- as.double(distances)## need for call in MLE target
  assign("MLEMIN",Inf,envir=ENVIR) ##necessary here, as LSoptim calls MLEtarget!
  LSopt <- LSoptimize(index, variab)
  PARAM <- LSopt$param
  variab <- LSopt$variab
  LSpmin <- LSopt$LSpmin
  LSpmax <- LSopt$LSpmax
  
  ## if rather close to nugget and nugget==fixed, do exception handling.
  ## This part is active only if
  ## scale.max.relative.factor < lowerbound.scale.LS.factor
  ## By default it is not, i.e., the following if-condition
  ## will "always" be FALSE.
  if (standard.style && PARAM[SCALE] < mindistances/scale.max.relative.factor
      && !is.na(nugget)) {
    if (PrintLevel>1) cat("Looks like if a nugget effect is included...\n")
    if (nugget!=0.0)
      stop("Chosen model does not seem to be appropriate!")
    ctr$param[OLD.NUGGET.POSITION] <- NA;
    return(fitvario(x=x, T=T, data=data, model=ctr, param=NULL,
                 lower=save.lower, upper=save.upper, sill=sill,
                 trend=new$trend,
                 use.naturalscaling=use.naturalscaling,
                 PrintLevel=PrintLevel, trace.optim=trace.optim,
                 bins=save.bins, nphi=nphi, ntheta=ntheta, ntime=ntime,
                 distance.factor=distance.factor,
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
                 pch=pch, var.name=var.name, time.name=time.name,
                 transform=save.transform, standard.style=standard.style
                 )
             )
  }

 
###################   optimization by MLE    ########################

  LEVEL <- "1st MLE :"
  assign("MLEMIN", LSopt$minimum,envir=ENVIR)
  assign("MLEPARAM", PARAM,envir=ENVIR)
  if (PrintLevel>2) cat("\nMLEMIN(LSQ)=", MLEMIN, "MLE optimisation...")
  else cat(pch)


  ## lowerbound.scale.LS.factor <  lowerbound.scale.factor, usually
  ## LS optimisation should not run to a boundary (what often happens
  ## for the scale) since a boundary value is usually a bad initial
  ## value for MLE (heuristic statement). Therefore a small
  ## lowerbound.scale.LS.factor is used for LS optimisation.
  ## For MLE estimation we should include the true value of the scale;
  ## so the bounds must be larger. Here lower[SCALE] is corrected
  ## to be suitable for MLE estimation
  if (pm$anisotropy)
    upper[scale.pos>0] <-
      upper[scale.pos>0] * lowerbound.scale.factor / lowerbound.scale.LS.factor
  else 
    lower[scale.pos>0] <-
      lower[scale.pos>0] * lowerbound.scale.LS.factor / lowerbound.scale.factor
  LB  <- lower[index]
  
  options(show.error.messages = FALSE) ##  
  if (length(variab)==1) {
    variab <- try(optimize(MLEtarget, lower = LB, upper = UB)$minimum)
  } else {
    variab <- try(optim(variab, MLEtarget, method="L-BFGS-B", lower = LB,
                        upper = UB, control=optim.control)$par)
  }
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
    if (PrintLevel>2) cat("MLE at ", format(c(MLEPARAM,MLEMIN),dig=4),"\n")
    else cat(pch)
    gridlength <- max(3, round(approximate.functioncalls^(1/sum(index))))
    ## grid is given by the extremes of the LS results
    ## so, therefore we should examine above at least 4 different sets
    ## of weights

    step <- (LSpmax-LSpmin) / (gridlength-2) # grid starts bit outside
    LSpmax <- pmin(LSpmax + step/2, UB)       # the extremes of LS
    LSpmin <- pmax(LSpmin - step/2, LB)
    step <- (LSpmax-LSpmin) / (gridlength-1)
    i <- 1
    zk <-  paste("LSpmin[",i,"] + step[",i,"] * (0:",gridlength-1,")")
    if (length(step)>1)
      for (i in 2:length(step))
        zk <- paste(zk,",LSpmin[",i,"] + step[",i,"] * (0:",gridlength-1,")")
    zk <- paste("expand.grid(",zk,")")
    startingvalues <- eval(parse(text=zk))
    limit <- 10 * approximate.functioncalls
    if ((rn <- nrow(startingvalues)) > limit) {
      if (PrintLevel>2)
        cat("using only a random subset of the", rn, "starting values")
      rand <- runif(rn)
      startingvalues <- startingvalues[rand < quantile(rand,  limit / rn), ]
      gc()
    }
    MLEMIN.old <- MLEMIN
    MLEPARAM.old <- MLEPARAM
    MLETREND.old <- MLETREND
    MLEMIN <- Inf

    LEVEL <- "grid ML :"
    apply(startingvalues, 1, MLEtarget) ## side effect: Minimum is in
    ##                                   MLEMIN !
    if (length(variab)>1) {
      PARAM <- MLEPARAM
      variab <- PARAM[index]
      options(show.error.messages = FALSE) ##
      LEVEL <- "2nd MLE :"
      
      if (PrintLevel>2) cat("second ML optimisation...\n") else cat(pch)
      variab <- try(optim(variab, MLEtarget,
                          method ="L-BFGS-B",
                          lower = LB, upper = UB,
                          control=optim.control)$par)      
      options(show.error.messages = TRUE) ##
      if (!is.numeric(variab) && (PrintLevel>0)) cat("MLEtarget II failed.\n")
      ## do not check anymore whether there had been convergence or not.
      ## just take the best of the two strategies (initial value given by
      ## LS, initial value given by a grid), and be happy.
    }
    if (MLEMIN<MLEMIN.old) {      
      PARAM <- MLEPARAM
    }
    else {
      PARAM <- MLEPARAM.old
      MLEMIN <- MLEMIN.old
      MLETREND <- MLETREND.old
    }
    variab <- PARAM[index]        
  }
  
  #if (PrintLevel>2) cat("\nMLEMIN=", MLEMIN, "nearly finished...")

  ## if the covariance functions use natural scaling, just
  ## correct the final output by GNS$natscale
  ## (GNS$natscale==1 if no rescaling was used)
  ##
  ## currently natural scaling only for standard.style...
  if (standard.style && use.naturalscaling && index[SCALE]) {
    GNS <- .C("GetNaturalScaling",
              covnr,
              as.double(PARAM[-1:-(KAPPA-1)]),
              scalingmethod,
              natscale=double(1),
              error=integer(1),
              PACKAGE="RandomFields", DUP=FALSE)
    if (GNS$error)
      stop(paste("Error", error, "occured whilst rescaling"))
    PARAM[SCALE] <- PARAM[SCALE] * GNS$natscale
    
    for (i in 1:ncol(lsqPARAM)) {
      GNS <- .C("GetNaturalScaling",
                covnr,
                as.double(lsqPARAM[-1:-(KAPPA-1), i]),
                scalingmethod,
                natscale=double(1),
                error=integer(1),
                PACKAGE="RandomFields", DUP=FALSE)
      if (GNS$error)
        stop(paste("Error", error, "occured whilst rescaling"))
      lsqPARAM[SCALE, i] <- lsqPARAM[SCALE, i] * GNS$natscale
    }
  }
  if (pch!="") cat("\n")

  if (is.null(pm$trend)) {
    TC <- NULL
    M <- if (is.na(pm$mean)) MLETREND else pm$mean
  } else {
    M <- NULL
    TC <- MLETREND
  }
  
  l <- list(covnr=pm$covnr, anisotropy=pm$anisotropy,
            op=pm$op, mean=M, trend=pm$trend,
            method=pm$method, timespacedim=pm$timespacedim)
  
  return(list(mle.value = -0.5 * (MLEMIN + lcrepet * log(2 * pi)),
              trend.coff=TC,
              mle=convert.to.readable({l$param <- PARAM; l}),
              ev = ev,
              lsq=convert.to.readable({l$param <- lsqPARAM[, 2]; l}),
              nlsq=convert.to.readable({l$param <- lsqPARAM[, 3]; l}),
              slsq=convert.to.readable({l$param <- lsqPARAM[, 4]; l}),
              flsq=convert.to.readable({l$param <- lsqPARAM[, 1]; l}),
              mle.lower=convert.to.readable({l$param <- lower;
                          if (is.null(M)) {l$trend <- -Inf; l$mean <- NULL;}
                          else {l$mean <- -Inf; l$trend <- NULL}; l}),
              mle.upper=convert.to.readable({l$param <- upper;
                          if (is.null(M)) {l$trend <- Inf; l$mean <- NULL;}
                          else {l$mean <- Inf; l$trend <- NULL}; l})
              ))
}

