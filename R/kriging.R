
## source("modelling.R")

### in Init als Parameter uebergeben; inklusive 
loc.input <- function(x, y, z, T) {
  x.given <- !missing(x) && !is.null(x)
  y.given <- !missing(y) && !is.null(y)
  z.given <- !missing(z) && !is.null(z)
  T.given <- !missing(T) && !is.null(T)
  env = new.env()
  return (function(f) {
    n <- names(as.list(args(f)))
    n <- n[-length(n)]
  #  if (any(is.na(pmatch(n, c("x", "y", "z", "T")))))
  #      stop("function arguments may only be 'x', 'y', 'z', 'T'")
    if (xor(x.given, "x" %in% n) || xor(y.given, "y" %in% n) ||
        xor(z.given, "z" %in% n) || xor(T.given, "T" %in% n))
      stop("mismatch of arguments")    
    
    if (!x.given) {
      if (!T.given || y.given || z.given) 
        stop("if 'x' is not argument only 'T' can be an argument")
      fctcall <- "fct(x)"
      res <- f(0)
    } else {
      if (y.given) {
        fctcall <- "fct(x=x[1], y=x[2]"
        xlen <- 2
        if (z.given) {
          fctcall <- paste(fctcall, ", z=x[3]")
          xlen <- 3
        }
        if (T.given) {
          fctcall <- paste(fctcall, ", T=x[length(x)]")
          xlen <- xlen + 1
        }      
      } else {
        if (z.given) stop("mismatch of arguments!")
        xlen <- if (is.vector(x)) 1 else ncol(x)
        fctcall <- "fct(x"
        if (T.given) {
          fctcall <- paste(fctcall, "[-length(x)], T=x[length(x)]")
          xlen <- xlen + 1
        }
      }
    }

    g <- NULL
    eval(parse(text=paste("g <- function(x, fct) {", fctcall, ") }")))
    res <- g(rep(0, xlen), f)
    if (is.vector(res)) res <- as.matrix(res) else
    if (!is.matrix(res)) stop("the function should return a vector or a matrix")
    vdim <- ncol(res)
    cond <- nrow(res)
    
    return(list(g=g, fct=f, vdim=vdim, cond=cond, env=env))
  })
  
}


rfPrepareNames <- function(model, coords, chkXT, locinfo) {
  coord.names <- if (!missing(coords)) colnames(coords) else NULL
  variab.names <- extractVarNames(model)
                       # gets response part of model, if model is a formula

  if (missing(chkXT)) {
    Time <- locinfo$Time
    ts.xdim <- locinfo$spatialdim + locinfo$Time
  } else {
    Time <- chkXT$Time
    ts.xdim <- chkXT$spacedim + chkXT$Time
  }
  if (is.null(coord.names))
    coord.names <- paste("coords.x", 1:ts.xdim, sep="")
  if (Time) coord.names[ts.xdim] <- "coords.T1"

  return(list(coord.names=coord.names, variab.names=variab.names))
}

data.starting.col <- function(data, halt=TRUE) {
  # str(data)
  if (is.null(colnames(data))) {
    if (halt) stop('colnames of data argument must contain "data" or "variable"')
    else return(-1);
  }
  is.data <- (tolower(substr(colnames(data), 1, 4)) == "data" |
              tolower(substr(colnames(data), 1, 4)) == "value" |
              tolower(substr(colnames(data), 1, 8)) == "variable")
  if (!any(is.data)) {
    if (halt) stop('no colname starts with "data" or "variable"')
    else return(-2);
  }
  return(which(is.data)[1])
}


rfPrepare <- function(model, x, y, z, T, distances, grid, data,
                      fillall, names.only=FALSE,
                      na.rm = "any", # "all", "none"
                      unconditional = FALSE,
                      ...) {  
  missing.x <- missing(x)
  if (missing(data)) stop("missing data")
  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")

  if (missing.orig.model <- missing(model)) model <- list("null")
  else model <- PrepareModel2(model, ...)

  
  T.tmp <- NULL

  if (isRFsp <- is(data, "RFsp")) {
    if (hasArg(given))
      stop("the coordinates for the measured data are given ambigiously (by 'given' and within 'data')")
    data <- selectDataAccordingFormula(data, model=model)
    if (!is.null(data@.RFparams)) {
      vdim <- data@.RFparams$vdim
      repet <- data@.RFparams$n
    } else {
      vdim <- repet <- NULL
    }
    variab.names <- colnames(data@data)

  #  Print(data, coordinates(data), RFspDataFrame2conventional(data)); kkkkkk
    
    grid.data <- isGridded(data)
    if (grid.data) {
      data <- RFspDataFrame2conventional(data) ## zeitkritisch ! to do
      x.tmp <- data$x
      given <- apply(x.tmp, 2, function(x) seq(x[1], by=x[2], length.out=x[3]))
      if (!is.list(given)) {
        given <- unlist(apply(given, 2, list), recursive=FALSE)
      }

      #Print(given)      
      given <- as.matrix(do.call("expand.grid", given))
      
      T.tmp <- data$T
      data <- as.matrix(data$data)
    } else {
      given <- coordinates(data) ## stehend; darf nicht gaendert werden
      data <- as.matrix(data@data)
    }
  } else {
    vdim <- repet <- NULL
 #    if (missing(vdim)) vdim <- 1
 #   if (missing(repet)) repet <- 1
  }

  #  Print(given, data); # ppp
 
  if (is.vector(data)) dimdata <- c(length(data), 1)
  else {
    if (!is.array(data)) data <- as.matrix(data)
    dimdata <-  base::dim(data)
  }
  base::dim(data) <- c(dimdata[1], prod(dimdata[-1]))

   if (missing.x) { ## generate locations for kriging
    ## then data cannot have repeated meseasurements
#
    stopifnot(is.null(y), is.null(z), is.null(T))
    if (isRFsp) {
      ## standard: data where at least 1 component or repetition
      ##           is NA is reestimated
      ## This is overwritten by two facts: fillall=TRUE or no NAs at all
      ## In this case all the data are reestimated.

      
      data.na <- rowSums(is.na(data)) > 0
      ## Print(data.na, length(data.na), given, grid)
      complete <- fillall || !any(data.na)
      if (complete) data.na <- rep(TRUE, length(data.na))
      else grid.data <- FALSE
      compareGridBooleans(grid, grid.data)
      grid <- grid.data
      T <- T.tmp
      if (!grid) x.tmp <- if (fillall) given else given[data.na, , drop=FALSE]
    } else {
      grid <- FALSE    ## Beibehaltung von grid, falls grid waere besser !
      data.col <- data.starting.col(data) : ncol(data)
      data.na <-  rowSums(is.na(data[, data.col, drop=FALSE])) > 0
      if (fillall || !any(data.na)) data.na <- data.na | TRUE
      x.tmp <- data[data.na, -data.col, drop=FALSE]
    }
    # if (!any(data.na)) stop("neither 'x' is given nor are there NAs in 'data' to be filled")
  } else {
     x.tmp <- x
  }

#  Print(missing.x, grid.data)

  neu <- CheckXT(x = x.tmp, y = y, z = z, T = T, grid = grid, 
                 distances=distances, 
                 if (!missing(distances) && length(distances) > 0) spdim=1
                 ##, length.data=length(data)
                 )
  
  
  ts.xdim <- as.integer(neu$spacedim + neu$Time)

  if (!isRFsp) {
    if (hasArg(given)) {
      given <- list(...)$given  ## stehend
      variab.names <- try(silent = TRUE, colnames(data))      
    } else {
      if (unconditional) {
        vdimrepet <- as.integer(length(data) / neu$restotal)
        if (vdimrepet * neu$restotal != length(data))
          stop("dimension of data does not match coordinates")
        dim(data) <- c(neu$restotal, vdimrepet)
        given <- NULL
        variab.names <- try(silent = TRUE, colnames(data))
      } else {
        given <- data[, 1:ts.xdim, drop=FALSE]
        idx <- try(silent = TRUE, data.starting.col(data):ncol(data))
        data2 <- try(silent = TRUE, data[, idx, drop=FALSE])
        variab.names <- try(silent = TRUE, colnames(data[ , idx, drop=FALSE]))
        if (class(data2) == "try-error" || class(variab.names) == "try-error") {
          idx <- -(1:ts.xdim)
          if (class(data2) == "try-error") {
            data2 <- try(data[, idx, drop=FALSE])
            if (class(data2) == "try-error")
              stop("columns of 'data' cannot be identified")
          }
          if (class(variab.names) == "try-error") {
            variab.names <- try(silent = TRUE, colnames(data[,idx, drop=FALSE]))
          }         
        }
        data <- data2
        rm("data2")
      }
    }
  }
  if (ncol(given) != ts.xdim && !unconditional) {
    #Print(given, ts.xdim, unconditional)
    stop("dimension of kriging points does not match dimension of given points")
  }
   
  names <- rfPrepareNames(model=model, coords=x.tmp, chkXT=neu)
  if (is.null(names$variab.names))
    names$variab.names <-
      if (class(variab.names) == "try-error") NULL else variab.names

  if (names.only) return(list(coord.names=names$coord.names,
                              variab.names=names$variab.names))

  ## extract on which data should be conditioned
  if (na.rm == "any") Data.NA <- rowSums(is.na(data)) > 1
  else if (na.rm == "all")  Data.NA <- rowSums(!is.na(data)) == 0
  else if (na.rm =="none") Data.NA <- FALSE
  else stop("unknown option for 'na.rm'.")

  fullgiven <- given
  fulldata <- data
  if (any(Data.NA)) {
    given <- given[!Data.NA, , drop=FALSE] # ts.xdim x number
    data <- data[!Data.NA, , drop=FALSE]    
  }
  total <- length(data)
 

 return(c(neu,                # locations x dim
          ts.xdim=ts.xdim,
          vdim =vdim,
          repet = repet,
          list(given=given, # locations x dim; given is n x dim matrix
               data=data,   # repet & multivariate alles in den Spalten
               ##             zuerst multivariat dann 
               xgr = if (neu$grid) cbind(neu$x, neu$T) else NULL,   # 3 x dim
               coord.names=names$coord.names,
               variab.names=names$variab.names,            
               model = model,
                                        #userdefined = userdefined,
               data.col = if (missing.x && !isRFsp) data.col else NULL,
               data.na  = if (missing.x) which(data.na) else NULL,
               fullgiven= fullgiven,
               fulldata = fulldata),
          total = total,        
          missing.x=missing.x)
         )
}



rfPrepareData <- function(model, x, y, z, T,
                          distances,                          
                          grid, data,
                          cathegories, reg, fillall, names.only=FALSE, ...) {


  
  neu <- rfPrepare(model=model, x=x, y=y, z=z, T=T, distances=distances,
                   grid=grid, data=data, reg=reg, fillall=fillall,
                   names.only = names.only, ...)
  ts.xdim <- neu$ts.xdim

#  Print(neu); str(data); #oooo

  trend <- NULL
  if (!is.null(cathegories)) {
    pm <- rfSplitTrendAndCov(model=neu$model, spatialdim=ts.xdim,
                             xdimOZ=ts.xdim, Time = FALSE,
                             cathegories=cathegories)
     if (!is.null(pm$trend)) {
       xgiv <- if (is.null(x)) NULL else
               if (is.null(y) && is.null(T)) neu$given else
               if (is.null(y)) neu$given[, -ts.xdim] else neu$given[, 1]
       ygiv <- if(is.null(y)) NULL else neu$given[, 2]
       zgiv <- if(is.null(z)) NULL else neu$given[, 3]
       Tgiv <- if(is.null(T)) NULL else neu$given[, ts.xdim]
       ## subtract the trend part that is known:

#       Print(neu$data, RFsimulate(model=pm$trend, x=xgiv,y=ygiv, z=zgiv, 
#                         T=Tgiv, grid=FALSE, register = reg, spConform=FALSE))
       #Print(pm$model)

       
       trend <- pm$trend
       neu$data <- neu$data - RFsimulate(model=trend, x=xgiv,y=ygiv, z=zgiv, 
                         T=Tgiv, grid=FALSE, register = reg, spConform=FALSE)

       ## parameter
    }
  } else {
    pm <- list(cov=model)
  }

  
  # userdefined <- GetParameterModelUser(pm$cov)
  xy <- matrix(nrow=if (neu$grid) 3 else ts.xdim - neu$Time,
                 ncol=if (neu$grid) ts.xdim - neu$Time  else 3,
                 as.double(1:3)) ## nur dummies

  vdim <- 
    .Call("Init",
          as.integer(reg),
          list("CovMatrix", pm$cov),  ## geaendert
          xy, ## nur dummies
          xy, ## nur dummies
          as.double(rep(1, 3)),
          as.integer(neu$spacedim),
          as.logical(neu$grid),
          as.logical(neu$distances),
          as.logical(neu$Time),
          NAOK=TRUE,
          PACKAGE="RandomFields")
  stopifnot(diff(vdim) == 0)
  vdim <- vdim[1]
  if (!is.null(neu$vdim) && vdim != neu$vdim)
    stop("multivariate dimensions of data and model do not match")
 
   # Vorsicht: vorige 2 Zeilen muessen vorher kommen

#  Print(neu$grid)
  
  if (neu$grid) {
    xlist <- mapply(function(a,b,c) a + b * (0:(c-1)),
                    neu$xgr[1,], neu$xgr[2,], neu$xgr[3,], SIMPLIFY=FALSE)
     ## 'x' will be used in apply within kriging
    neu$x <- do.call(base::expand.grid, xlist)
  } else if (neu$Time)
    neu$x <- cbind(matrix(rep(t(neu$x), times=neu$T[3]),
                        ncol=if (is.matrix(neu$x)) ncol(neu$x) else 1,
                        byrow=FALSE),
                 rep(neu$T[1] + neu$T[2] * (0:(neu$T[3]-1)),
                     each=if (is.matrix(neu$x)) nrow(neu$x) else length(neu$x)))

  repet <- neu$total / (vdim * nrow(neu$given))

#  Print(neu, vdim)
  
  if (repet != as.integer(repet)) {
    # Print(length(all$data), vdim)
    stop("number of data not a multiple of the number of locations")
  }
#

  neu$vdim <- vdim
  neu$repet <- repet

  return(c(neu, pm))
}


RFinterpolate <- function(model, x, y=NULL, z=NULL, T=NULL, grid, data,
                          distances, dim, err.model, method="ml", ...) {
  #Print("entering interpolate")
                                        #  str(data)

  if (!missing(err.model)) stop("'err.model' not programmed yet.")
  if (!missing(distances) && length(distances) > 0) stop("'distances' not programmed yet.")
  # if (!missing(dim)) warning("'dim' is ignored.")
  krige.methlist <- c("A", "S", "O", "M", "U", "I")
  ## currently only univariate Kriging possible

  opt <- list(...)
  i <- pmatch(names(opt), c("MARGIN", "given"))
  opt <- opt[is.na(i)]
  
  RFoptOld <- do.call("internal.rfoptions",
                      c(opt, RELAX=isFormulaModel(model)))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  ## eingabe wird anstonsten auch als vdim_close erwartet --
  ## dies ist nocht nicht programmiert! Ausgabe ist schon programmiert
  ## CondSimu ist auch noch nicht programmiert
  if (RFopt$general$vdim_close_together)
    stop("'vdim_close_together' must be FALSE")

  maxn <- RFopt$krige$locmaxn
  split <- RFopt$krige$locsplitn[1]
  fillall <- RFopt$krige$fillall
  spConform <- RFopt$general$spConform

 #maxn <- 70
#  split <- 40

  reg <- RFopt$registers$interpolregister
  return.variance <- RFopt$krige$return_variance
  
  ##  MINV <<- NULL
  SOLVE <- if (RFopt$krige$cholesky_R) {
    function(M, v) {
      sqrtM <- chol(M)
      if (RFopt$general$printlevel>=PL.FCTN.SUBDETAILS)
        Print(range(diag(sqrtM))^2) #
      if (missing(v)) chol2inv(sqrtM) else chol2inv(sqrtM) %*% v
    }
  } else solve
  
  silent = TRUE
  nuggetrange <- 1e-15 + 2 * RFopt$nugget$tol
  krige.meth.nr <- pmatch(RFopt$krige$method, krige.methlist) -1

  if (FALSE && krige.methlist[krige.meth.nr + 1] == "O") { ###### FALSE
    new.model <- PrepareModel2(model, ...)
    if (new.model[[1]] %in% ZF_PLUS) {
      new.model[[length(new.model)+1]] <- list(ZF_TREND[2], mean=NA)
    } else {
      new.model <- list(ZF_SYMBOLS_PLUS, new.model,  list(ZF_TREND[2], mean=NA))
    }
    return(RFinterpolate(new.model, x=x, y=y, z=z, T=T, grid=grid, data=data,
                         ..., krige.method="U"))
  } # end FALSE
  
  if (is.na(krige.meth.nr))
    stop(paste("kriging method not identifiable from the list.",
               "Possible values are", paste(krige.methlist, collapse = ",")))
  krige.mean <- which(krige.methlist == "M") - 1
   
  cathegories <- list(trend=c(DetTrendEffect, DeterministicEffect),
                      estimtrend=c(FixedTrendEffect),
                      random=c(FixedEffect:SpVarEffect))


  #Print(fillall)
#  Print(class(model), model[method])
  if (class(model) == "RF_fit") model <- model[[method]]$model
  else if (class(model) == "RFfit") model <- PrepareModel2(model[method])
 #  Print(model)
  all <- rfPrepareData(model=model, x=x, y=y, z=z, T=T, grid=grid,
                       data=data, cathegories=cathegories, reg=reg,
                       fillall=fillall, ...)

#  Print(all); kkk

   #userdefined <- all$userdefined
  
  ts.xdim <- as.integer(all$ts.xdim)
  ngiven <- as.integer(nrow(all$given)) ## number of given points
  repet <- as.integer(all$repet)
  vdim = all$vdim

  if (!RFopt$general$na_rm_lines && repet > 1 && all$missing.x) {
     base::dim(all$fulldata) <- c(length(all$fulldata) / repet, repet)
    for (i in 1:repet) {
      dummy <- if (hasArg(given)) {
        RFinterpolate(model=model, x=x, y=y, z=z, T=T, all$grid,
                      data = all$fulldata[, i], 
                      ..., spC=if (i==1) spConform else FALSE)
      } else {
        RFinterpolate(model=model, x=x, y=y, z=z, T=T, all$grid,
                      data = all$fulldata[, i],
                      given = all$given,
                      ..., spC=if (i==1) spConform else FALSE)
        
      }
      
      if (i==1) {
        if (return.variance) {
          v <- dummy$var
          res <- dummy$estim
        }
        if (spConform) {
          str(dummy) #
          coords <- res@coords
          stop("spConform output for missing values not programmed yet. Please contact author")
        } else {
          d <- if (is.vector(dummy)) length(dummy)  else  base::dim(dummy)
          res <- array(0, dim=c(d, repet)) 
        }
      }

      if (return.variance) {
        res[, i] <- dummy$estim
        var[is.na(all$fulldata[, i]), i] <- if (i==1) v else dummy$var
      } else {
        res[, i] <- res
      }
    }
  
    if (spConform){
      res <- conventional2RFspDataFrame(res,
                                        coords=coords,
                                        gridTopology=gridTopology,
                                        n=repet,
                                        vdim=vdim,
                                        vdim_close_together =
                                        RFopt$general$vdim_close_together)
      if (return.variance){
        var <- conventional2RFspDataFrame(var,
                                          coords=coords,
                                          gridTopology=gridTopology,
                                          n=1,
                                          vdim=vdim,
                                          vdim_close_together =
                                          RFopt$general$vdim_close_together)
        names(var@data) <- paste("var.", names(var@data), sep="")
        res <- cbind(res, var)
        res@.RFparams$has.variance <- TRUE
      }
      if (is.raster(x)) {
        res <- raster::raster(res)
        projection(res) <- projection(x)
      }
      return(res)
    } else {
      return(if (return.variance) list(res, var) else res)
    }
  }  # end !na_rm_lines   


  trendfct <- NULL
  polydeg <- NULL
  
  if (!is.null(all$random))
    stop("random effects cannot be treated within kriging")
  if (length(all$estimtrend)>0) {
    trendparlist <- c("mean", "plane", "polydeg", "arbitraryfct")
    pref.methlist <- list("O", "I", "I", "U")
    add.methlist <- list(NULL, "U", "U", NULL)
    if (all$estimtrend[[1]] %in% ZF_PLUS) {
      all$estimtrend[[1]] <- NULL
      if (any(sapply(all$estimtrend, 
                     function(ll) length(ll$arbitraryfct)) == 0))
        stop("complex trend may consist of arbitrary functions only")
      trendfct <- lapply(all$estimtrend, function(ll) ll$arbitraryfct)
      if (krige.meth.nr==0)
        krige.meth.nr <- pmatch("U", krige.methlist)-1 #arbitraryfct
      else if (krige.methlist[krige.meth.nr+1] != "U")
        stop("use universal kriging for more complicated trends")
    } else {
      stopifnot(length(all$estimtrend)==2)
      trendparnr <- pmatch(names(all$estimtrend)[2],trendparlist)
      if (is.na(trendparnr))
        stop("unknown trend does not match kriging specification")
      if (krige.meth.nr == 0)
        krige.meth.nr <- pmatch(pref.methlist[[trendparnr]],
                                krige.methlist) - 1
      if (!(krige.methlist[krige.meth.nr+1] %in%
            union(pref.methlist[[trendparnr]], 
                  add.methlist[[trendparnr]])) )
        stop("unknown trend does not match kriging specification")
      
      if (krige.meth.nr !=
          (krige.meth.nr <- pmatch(pref.methlist[[trendparnr]],
                                   krige.methlist) - 1         ))
        warning("use intrinsic kriging for unknown trend plane / polynomial")
      
      switch(trendparnr ,{ ##mean
      },{ ## plane
        polydeg <- rep(1, times=ncol(all$estimtrend$plane))
      },{ ##  polydeg
        polydeg <- all$estimtrend$polydeg
      },{ ## arbitraryfct
        trendfct <- all$estimtrend$arbitraryfct
      })
    }
  } else {
    if (krige.meth.nr==0) {
      info <- RFgetModelInfo(reg, level=3)
      abbr <- if (isNegDef(info$type) && !isPosDef(info$type)) "I" else "S"
      krige.meth.nr <- pmatch(abbr, krige.methlist) - 1
    }
  }
  
  if (krige.meth.nr==0) stop("auto kriging detecting cannot be resolved") 
  ##old.data <- data

  
  ## die zu interpolierenden Orte
  xx <- t(all$x)
  nx <- as.integer(ncol(xx))
  
  exact <- RFopt$general$exact
  if (ngiven > maxn || !is.na(exact) && !exact && ngiven > split) {
    ## neighbourhood kriging !
    if (!is.na(exact) && exact)
      stop("number of conditioning locations too large for an exact result.")
    if (ngiven > maxn && is.na(exact) && RFopt$general$printlevel>=PL.IMPORPANT)
      message("performing neighbourhood kriging")

    ## calculate the boxes for the locations where we will interpolate
    idx <- GetNeighbourhoods (MODEL.INTERN, all,
                              RFopt$krige$locsplitfactor,
                              RFopt$krige$locmaxn,
                              RFopt$krige$locsplitn)
    totalparts <- length(idx[[2]])
  } else {
    idx <- list(list(matrix(1:ngiven, ncol=1)),
                list(1),
                list(1:nx),
                list(nDsplitn=1)
                )
    totalparts <- 1
  }
  if (totalparts > 1) RFoptions(general.pch="")
  pr <- totalparts > 1 && RFopt$general$pch != "" &&  RFopt$general$pch != " "



  ################ KRIGING THE MEAN #################
#  old <- as.matrix(-1)
  
  if (krige.meth.nr == krige.mean) {
    if (return.variance) {
      if (vdim>1)
        stop("kriging variance for multivariate mean not programmed yet")
      if (totalparts>1)
        stop("variance cannot be identified for neighbourhood kriging.\nSet return.variance=FALSE")
    }
    if (hasArg(MARGIN)) {
      ## if given, then the mean is calculated for the "MARGIN"
      ## coordinates (e.g. space) , for each "non.const" coordinate (e.g. time)
      MARGIN <- list(...)$MARGIN
      if (is.character(MARGIN)) {
        stop("character for MARGIN not allowed yet")
         MARGIN <- pmatch(MARGIN, c("x", "y", "z", "T"))
      }
    } else {
      MARGIN <- integer(0)
    }   
    not.const <- MARGIN[MARGIN >= 1 & MARGIN <= ts.xdim]
    constant <- if (length(not.const) == 0) 1:ts.xdim else (1:ts.xdim)[-not.const]
    
  
    U <- list()
    if (length(not.const)==0) {
      totalCond <- where <- 1
    } else {
      ## get first the "non.const" coordinates, e.g. time. As we may not
      ## have a grid, we have to be careful getting them
      for (i in 1:length(not.const)) {
        U[[i]] <- sort(unique(xx[not.const[i], ]))
      }
      nU <- sapply(U, length)
      cnU <- c(1, cumprod(nU)) ## could be more than one not.const coordinates.
      ##                   we work on a grid then, for ease
      totalCond <- prod(nU)
    }
    sigma2 <- res <- array(dim=c(totalparts, repet * vdim, totalCond))

    if (pr) cat(totalparts)
    for (p in 1:totalparts) {
      if (pr && p %% 1==0) cat(RFopt$general$pch)
      givenidx <- unlist(idx[[1]][idx[[2]][[p]]])
      dat <- all$data[givenidx, , drop=FALSE]
      
      XX <- xx[, idx[[3]][[p]], drop=FALSE]
       
      given <- t(all$given[givenidx, , drop=FALSE])
      storage.mode(given) <- "double"
      Ngiven <- as.integer(ncol(given)) 
      ngvdim <- Ngiven * vdim

      notna <- as.vector(is.finite(dat))
      dat <- dat[notna]
  
      base::dim(dat) <- c(length(dat) / repet, repet)# Voraussetzung, dass es aufgeht!
      storage.mode(dat) <- "double"
      covmatrix <- double(ngvdim^2)

      .Call("CovMatrixIntern", reg, given, FALSE, FALSE, # no dist, no grid
            Ngiven, covmatrix, # userdefined,
            PACKAGE="RandomFields")
       
       base::dim(covmatrix) <- c(ngvdim, ngvdim)
      covmatrix <- covmatrix[notna, notna]
      
      if (length(not.const)==0) {
        conditions <- vdim
        onevector <- c(rep(c(rep(1, times=Ngiven), 
                             rep(0, times=Ngiven*vdim)), times=vdim-1),
                       rep(1,times=Ngiven))
         base::dim(onevector) <- c(Ngiven * vdim, vdim)
        covmatrix <- rbind(cbind(covmatrix, onevector), 
                           cbind(t(onevector), matrix(0, nrow=vdim, ncol=vdim)))
      } else {
        if (vdim > 1) stop("multivariate version of marginal kriging not programmed yet")
        u <- v <- list()        
        for (i in 1:length(not.const)) {
          ## die XX definieren mir, wo geschaetzt werden soll
          u[[i]] <- sort(unique(XX[not.const[i], ]))
          v[[i]] <- pmatch(u[[i]], U[[i]]) * cnU[i]
        }
        uall <- expand.grid(u)
        where <- rowSums(expand.grid(v))
        conditions <- nrow(uall)
        lagrange <- matrix(0, nrow=conditions, ncol=ngvdim)
        
        for (i in 1:conditions) {
          sameinstance <-  apply(given[not.const, , drop=FALSE] == uall[i, ],
                                 2, all)
          lagrange[i, sameinstance] <- 1
        }
        covmatrix <-
          cbind(rbind(covmatrix, lagrange),
                rbind(t(lagrange), matrix(0, nrow=conditions, ncol=conditions)))
      }
      
      lambda.mu <-
        try(SOLVE(covmatrix,
                  rbind(matrix(0, nrow=ngvdim, ncol=conditions),
                        diag(conditions))),
            silent = silent)
       if (!(is.numeric(lambda.mu))) {
         stop("Covmatrix is singular .") #
       }
    
      res[p, , where] <- as.vector(crossprod(dat, lambda.mu[1:ngvdim, ]))  
      mu <- lambda.mu[-1:-ngvdim, ]
      sigma2[p, , where] <-
        if (length(mu) == 1) -mu else -diag(mu)  ## unklar ob es passt
    } # p in totalparts

    invsd <- 1/sqrt(sigma2)
    res <- apply(res * invsd, 2:3, sum, na.rm=TRUE) /
      apply(invsd, 2:3, sum, na.rm=TRUE)
    sigma2 <- apply(sigma2, 2:3, mean, na.rm=TRUE)
 
     base::dim(sigma2) <-  base::dim(res) <-
      c(totalCond, if (vdim>1) vdim, if (repet>1) repet)
    
    if (spConform) {
      if (all$missing.x) {
        Res <- all$fulldata
      } else {
        Res <- matrix(nrow=nx, ncol=repet * vdim)
      }
      
      if (length(not.const) > 0) {
        Uall <- expand.grid(U)
        for (i in 1:nrow(Uall)) {
          j <- which(xx[not.const, ]==Uall[i, ])
          Res[j, ] <- rep(res[i], each=length(j))
        }
      } else {
        Res[,] <- res[1]
      }
    } else {
      return(list(x=if (length(not.const) > 0) expand.grid(U),
                  estim = res,
                  var = if (return.variance) sigma2))
    }
  ############## END KRIGING THE MEAN ######
  } else {    
     if (all$missing.x) {
      Res <- all$fulldata
    } else {
      Res <- matrix(nrow=nx, ncol=repet * vdim)
    }


    for (p in 1:totalparts) {
      stopifnot((Nx <- as.integer(length(idx[[3]][[p]]))) > 0)
      if (pr && p %% 5==0) cat(RFopt$general$pch)
      res <- double(Nx * repet * vdim)
      storage.mode(all$data) <- "double"
      givenidx <- unlist(idx[[1]][idx[[2]][[p]]])
      Ngiven <- length(givenidx)
      dat <- all$data[givenidx, , drop=FALSE]
      base::dim(dat) <- c(Ngiven, vdim, repet)
      given <- t(all$given[givenidx,  , drop=FALSE])
      storage.mode(given) <- "double"
      Ngiven <- as.integer(ncol(given)) ## may some point had been deleted
      ngvdim <- Ngiven * vdim
      
      covmatrix <- double(ngvdim^2)

      pos <- integer(Ngiven)
      ## lexicographical ordering of vectors --- necessary to check
      ## whether any location is given twice, but with different value of
      ## the data
      .C("Ordering", given, Ngiven, ts.xdim,
         pos, PACKAGE = "RandomFields", DUP = DUPFALSE)
      pos <- pos + 1
      
      ## are locations given twice with different data values?
      check <- given[ , pos, drop = FALSE]
      if (any(dup <-
              c(FALSE, colSums(abs(check[, -1, drop = FALSE] -
                                   check[, -Ngiven, drop = FALSE])) == 0))) {
        if (RFopt$general$allowdistanceZero) {
          given[1, dup] <- check[1, dup] + rnorm(n=sum(dup), 0, nuggetrange)
          dat <- dat[pos, , , drop = FALSE]
        } else {
          dat <- dat[pos, , , drop = FALSE]##pos is recycled if vdim > 1
          if (any(dat[dup, , ] != dat[c(dup[-1], FALSE), , ]))
            stop("duplicated conditioning points with different data values")
          stop("duplicated values not allowed")
          given <- given[, !dup, drop = FALSE]
          Ngiven <- as.integer(ncol(given)) ## may some point had been deleted
          dat <- dat[!dup, , drop = FALSE]
          stop("duplicated values not allowed")
        }
        if (is.list(trendfct))
          return(matrix(sapply(trendfct, function(f) eval(body(f))),
                        nrow=nfct, byrow=TRUE))
        else return(matrix(eval(body(trendfct)), nrow=nfct, byrow=TRUE))
      }

      notna <- as.vector(is.finite(dat))
      base::dim(notna) <- c(length(notna) / repet, repet)
      if (any(xor(notna, notna[,1]))) {
        stop("kriging with repeated measurements but different NAs not programmed yet")
      } else {
        notna <- notna[ ,1]
      }
      
      dat <- dat[notna]
      base::dim(dat) <- c(length(dat) / repet, repet)
      storage.mode(dat) <- "double"
       
      ## works also for variograms and returns -(gamma(x_i-x_j))_{ij}
      

      .Call("CovMatrixIntern", reg, given, FALSE, FALSE,  Ngiven, covmatrix,
            #userdefined,
            PACKAGE="RandomFields")
       base::dim(covmatrix) <- c(ngvdim, ngvdim)

      

      covmatrix <- covmatrix[notna, notna]

      XX <- as.double(xx[, idx[[3]][[p]], drop=FALSE])


      switch(krige.meth.nr, {
         ## simple kriging

         stopifnot(is.null(all$estimtrend))
         if (return.variance) {
          if (!(is.numeric(try(invcov <- SOLVE(covmatrix), silent = silent)))){
            stop("Covmatrix is singular")
          }

          sigma2 <- double(Nx*vdim)

          .Call("simpleKriging2", reg, given, XX, dat,
             invcov, notna, Nx, Ngiven, ts.xdim, repet,
             res, sigma2, #userdefined,
             PACKAGE = "RandomFields")
        } else {          
          if (!(is.numeric(try(invcov <- SOLVE(covmatrix, dat),
                               silent = silent)))) {
#            print(covmatrix)
 #           print(given)
            
            if (RFopt$general$printlevel>2)
              #Print(covmatrix,all, sum(dat),  base::dim(dat),  base::dim(covmatrix), eigen(covmatrix)$values) 
            stop("Covmatrix is singular.")        
          } 

          #     Print(XX, length(XX), given, notna, Nx, Ngiven, ts.xdim, repet)

         .Call("simpleKriging", reg, given, XX, invcov,
             notna, Nx, Ngiven, ts.xdim, repet, res, #userdefined,             
               PACKAGE = "RandomFields")
        }
       }, {
        ## ordinary kriging
        onevector <- c(rep(c(rep(1, times=Ngiven), rep(0, times=Ngiven*vdim)),
                           times=vdim-1),
                       rep(1,times=Ngiven))
        base::dim(onevector) <- c(Ngiven*vdim, vdim)
        onevector <-  onevector[notna, ]
        #Print(onevector, covmatrix, rep(1, times=Ngiven))
        covmatrix <- rbind(cbind(covmatrix, onevector), 
                           cbind(t(onevector), matrix(0, nrow=vdim, ncol=vdim)))
       
        if (return.variance) {
          if (!(is.numeric(try(invcov <- SOLVE(covmatrix), silent = silent)))) {
            stop("Covmatrix is singular!")
          }
          sigma2 <- double(Nx*vdim)
          .Call("ordinaryKriging2", reg, given, XX, as.double(dat),
             invcov,
             notna, Nx, Ngiven, ts.xdim, repet, res, sigma2,
#             userdefined,
             PACKAGE = "RandomFields")
        }  else  {
          invcov <- try(SOLVE(covmatrix,
                              rbind(dat,matrix(0,nrow=vdim,ncol=repet))),
                      silent = silent)
          if (!(is.numeric(invcov))) {
            stop("Covmatrix is singular .") #
          }
          
          .Call("ordinaryKriging", reg, given, XX, invcov,
             notna, 
             Nx, Ngiven, ts.xdim, repet, res, #userdefined,
             PACKAGE = "RandomFields")
        }
      }, {
        ## kriging the mean
        stop("kriging the mean went wrong -- please contact author")
      }, {
        ## universal kriging
        nfct <- length(trendfct)
        spacedim <- ts.xdim - all$Time
                
        trend_expr <- trendeval <- NULL ## !! Marco ?!
        
           if(nfct > 0) {
          trendargs <- rep(0,times=ts.xdim)
          if(!is.numeric(try(testres <- eval(trend_expr), silent=silent)))
          stop("misspecified trend function")
          else if (length(testres)!=vdim*nfct)
            stop("dimension of trend function does not match model")
        }
        
        if (return.variance) {
          if (!(is.numeric(try(invcov <- SOLVE(covmatrix), silent=silent))))
            stop("covmatrix is singular!")
          sigma2 <- double(Nx*vdim)
          .Call("universalKriging2", reg, given, XX, dat, 
             invcov, notna, Nx, Ngiven, ts.xdim, repet,
             res, sigma2, as.integer(nfct), trend_expr, new.env(),#userdefined, 
             PACKAGE = "RandomFields")
        } else {
           base::dim(dat) <- c(Ngiven * vdim, repet)
          data_null <- rbind(dat, matrix(0.0, nrow=nfct, ncol=repet))
          tFmatrix <- matrix(nrow=nfct, ncol=Ngiven*vdim)
          for (i in 1:Ngiven)
            tFmatrix[,(0:(vdim-1))*Ngiven+i] <- trendeval(args=given[,i])
          KMat <- rbind(cbind(covmatrix,t(tFmatrix)), 
                        cbind(tFmatrix, matrix(0.0, nrow=nfct, ncol=nfct)))
          if (!(is.numeric(try(invcov <-
                               SOLVE(KMat, data_null, silent=silent)))))
            stop("covmatrix is singular..")
          .Call("universalKriging", reg, given, XX, invcov,
             notna,
             Nx, Ngiven, ts.xdim, repet, res, as.integer(nfct), trend_expr,
             new.env(), #userdefined,
                PACKAGE = "RandomFields")
        }
      }, {
        ## intrinsic kriging, notation see p. 169, Chiles 
        ## cross variogram 
        ## nullx <- matrix(0, nrow=ts.xdim, ncol=Ngiven)
        ## nullcov <- double(ngvdim^2)
        ## .Call("CovMatrixIntern", reg, nullx, FALSE,
        ##    FALSE, Ngiven, nullcov, userdefined,
        ##    PACKAGE="RandomFields")
        ##  base::dim(nullcov) <- c(ngvdim, ngvdim)
        ## #2*nu_ij(h) = 2*C_ij(0)-C_ij(h)-C_ij(-h) 
        ## covmatrix <- 0.5*covmatrix + 0.5*t(covmatrix) - nullcov
        ## pseudo cross variogram
        nullx <- matrix(0, nrow=ts.xdim, ncol=Ngiven)
        nullcov <- double(ngvdim^2)
        .Call("CovMatrixIntern", reg, nullx, FALSE, FALSE, Ngiven, nullcov,
#              userdefined,
              PACKAGE="RandomFields")

         base::dim(nullcov) <- c(ngvdim, ngvdim)
        
        nullcov <- rep(diag(nullcov), times=ngvdim)
         base::dim(nullcov) <- c(ngvdim, ngvdim)
        ## 2*gamma_ij(h) = C_ii(0) + C_jj(0) - 2*C_ij(h) 
        covmatrix <- covmatrix - 0.5*nullcov - 0.5*t(nullcov)
        ##
        if (is.null(polydeg)) {
          polydeg <- rep(0, times=vdim)
          message("!!!")
          if (FALSE)
          warning(paste("by default the polynomial degree for",
                        " intrinsic kriging is set to zero", collapse=""))
        } else if (!all(polydeg==(polydeg <- max(polydeg))))
          warning(paste("the degree of the polynomial trend for each component",
                        "is set to", polydeg[1], collapse=""))
        
        polydeg <- polydeg[1] ##!!!!
        
        if (return.variance) {
          ## Inverting CovMatrix and initializing sigma 2

#          Print(eigen(covmatrix)$val, covmatrix[length(covmatrix)])
          
          if (!(is.matrix(try(invcov <- SOLVE(covmatrix), silent = silent))))
            stop("covmatrix is singular.")
          sigma2 <- double(Nx*vdim)
          ##now: Kriging
    
          .Call("intrinsicKriging2", reg, given, XX, dat,
             invcov, notna, Nx, Ngiven, ts.xdim, repet,#9
             res, sigma2, as.integer(polydeg), #userdefined,
             PACKAGE = "RandomFields")
        } else {
          ##Initializing Fmatrix
          nfct <- choose(polydeg + ts.xdim, ts.xdim)
          powers <- integer(ts.xdim * nfct)
          .C("poly_basis_extern", ts.xdim, as.integer(polydeg),
             powers, PACKAGE = "RandomFields", DUP = DUPFALSE)
          powers <- matrix(powers, nrow=nfct, ncol=ts.xdim, byrow=TRUE)
          smallFmatrix <- matrix(0, nrow=Ngiven, ncol=nfct)
          for (j in 1:nfct)
            smallFmatrix[,j] <- apply(given^powers[j,], 2, prod)
          Fmatrix <- matrix(0, nrow=Ngiven*vdim, ncol=nfct*vdim)
          for (v in 1:vdim)
            Fmatrix[(v-1)*Ngiven+1:Ngiven,(v-1)*nfct+1:nfct] <- smallFmatrix
          KMat <- rbind(cbind(covmatrix,Fmatrix),
                        cbind(t(Fmatrix), matrix(0.0, ncol=nfct*vdim, 
                                                 nrow=nfct*vdim)))
           base::dim(dat) <- c(Ngiven * vdim, repet)
          data_null <- rbind(dat, matrix(0.0, ncol=repet, nrow=nfct*vdim))
          if (!(is.numeric(try(invcov <- SOLVE(KMat, data_null, silent=silent)))))
            stop("covmatrix is singular!")
        .Call("intrinsicKriging", reg, given, XX,
           invcov, notna,
           Nx, Ngiven, ts.xdim, repet, res, as.integer(polydeg), #userdefined,
           PACKAGE = "RandomFields")
        }
      })

      if (all$missing.x) {
         
        where <- all$data.na[idx[[3]][[p]]]
#      Print(where, length(where), diff(where), Res)
        isNA <- is.na(Res[where, ])
#      Print(isNA)

 #       Print(where, isNA, Res, res)
        
        Res[where, ][isNA] <- res[isNA]        
      } else {
        Res[idx[[3]][[p]], ] <- res
      }
      if (any(base::dim(Res) == 1)) Res <- as.vector(Res)
    
    } ## for p in totalparts
     
    if (pr) cat("\n")

    dimension <- if (all$grid) all$xgr[3, ]  else  nx
     
     newdim <- c(dimension, if (vdim>1) vdim, if (repet>1) repet)
    if (length(newdim)>1)  base::dim(Res) <- newdim
  
    if (return.variance && length(newdim <- c(dimension, if (vdim>1) vdim)) > 1)
       base::dim(sigma2) <- newdim
    
     
     if (length(all$trend) > 0){
       Res <- Res + RFsimulate(x=x, y=y, z=z, T=T, grid=all$grid,
                              model=all$trend, spC=FALSE, reg=reg,
                              n=repet) ## n=repet von Alex 4.12.2011;
 
    }
  } ## END OF  **NOT**  KRIGING THE MEAN

  
  
  if (!is.null(all$variab.names))
    attributes(Res)$variab.names <- all$variab.names

  if (!spConform && !all$missing.x) {
    if (vdim > 1 && RFopt$general$vdim_close_together) {
      Resperm <- c(length(dimension)+1, 1:length(dimension),
                   if(repet>1) length(dimension)+2)      
      Res <- aperm(Res, perm=Resperm)
      
      if (return.variance)
        sigma2 <- aperm(sigma2, perm=Resperm[1:(length(dimension)+1)])
    }
   
    return(if (return.variance) list(estim = Res, var = sigma2) else Res)
  }

  coord.names.incl.T <- c(if (!is.null(all$coord.names)) all$coord.names
  else paste("coords.x", 1:all$spatialdim, sep=""),
                          if (all$Time) paste("coords.T", 1, sep="")
                          else NULL)

  if (all$missing.x) {
    all$fillall <- fillall
    all$simu <- Res    
    # Print(data, all, 1, spConform, Res) 
    Res <- FinishImputing(data, all, 1, spConform)
    if (return.variance){
      var <- FinishImputing(data, all, 1, spConform)
      if (spConform) {
        names(var@data) <- paste("var.", names(var@data), sep="")
        Res@.RFparams$has.variance <- TRUE
        Res <-  cbind(Res, var)
      } else Res <- list(Res, var)
    }
 #   print(Res)
    return(Res)
  } else {
    if (!all$grid) {
      coords <- all$x
      colnames(coords) <- coord.names.incl.T
      gridTopology <- NULL
    } else { ## grid=TRUE:
      if (!FALSE) {
        coords <- NULL
        xgr <- all$xgr
        colnames(xgr) <- coord.names.incl.T
        xgr[is.na(xgr)] <- 0
        gridTopology <- GridTopology(xgr[1, ], xgr[2, ], xgr[3, ])
      } else {
        info <- RFgetModelInfo(reg, level=3)
        prep <- prepare4RFspDataFrame(model, info, x, y, z, T,
                                      all$grid, data, RFopt)
        attributes(Res)$variab.names <- prep$names$variab.names
      }
    }

    Res <- conventional2RFspDataFrame(Res, coords=coords,
                                      gridTopology=gridTopology,
                                      n=repet, vdim=vdim,
                                      vdim_close_together =
                                      RFopt$general$vdim_close_together)
  
    if (return.variance){
      var <- conventional2RFspDataFrame(sigma2, coords=coords,
                                        gridTopology=gridTopology,
                                        n=1, vdim=vdim,
                                        vdim_close_together =
                                        RFopt$general$vdim_close_together)
      names(var@data) <- paste("var.", names(var@data), sep="")
      Res@.RFparams$has.variance <- TRUE
      Res <-  cbind(Res, var)
    }
  }
  
  Res@.RFparams$krige.method <-
    c("Simple Kriging", "Ordinary Kriging", "Kriging the Mean",
      "Universal Kriging", "Intrinsic Kriging")[krige.meth.nr]
  

  ## Res@.RFparams$var <- sigma2 ## sehr unelegant.
  ## * plot(Res) sollte zwei Bilder zeigen
  ## * var(Res) sollte sigma2 zurueckliefern
  ## * summary(Res) auch summary der varianz, falls vorhanden
  ## * summary(Res) auch die Kriging methode

  if (is.raster(x)) {
    Res <- raster::raster(Res)
    projection(Res) <- projection(x)
  }
   
  return(Res)
}




rfCondGauss <- function(model, x, y=NULL, z=NULL, T=NULL, grid, n=1,
                        data,  # first coordinates, then data
                        err.model=NULL, ...) { # ... wegen der Variablen  
  dots <- list(...)
  if ("spConform" %in% names(dots))
    dots$spConform <- NULL

   ## currently only univariate Kriging possible

  RFoptOld <- internal.rfoptions(..., RELAX=isFormulaModel(model)) 
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  cond.reg <- RFopt$registers$condregister
  interpol.reg <- RFopt$registers$interpolregister

  cathegories <- list(trend=c(DetTrendEffect, DeterministicEffect),
                      random=c(FixedTrendEffect:SpVarEffect))


  all <- rfPrepareData(model=model,
                       x=x, y=y, z=z, T=T, grid=grid, data=data,
                       cathegories=cathegories,
                       reg=interpol.reg, fillall=FALSE, ...)
  ##  userdefined <- all$userdefined

#  Print(all)
  
  ##  krige.mean <- 0 ## all$mean + err$mean ?? 
  stopifnot(length(all$random) ==0)
 
  krige <- all$model
  vdim <- all$vdim
  ts.xdim <- all$ts.xdim
  
  if (vdim > 1) stop("multivariate version not programmed yet")
  
  if (!is.null(err.model)) {
    err <- rfSplitTrendAndCov(model=err.model, spatialdim=ts.xdim,
                              xdimOZ=ts.xdim, Time=FALSE)
    if (!is.null(err$trend)) stop("trend for the error term not allowed yet")
    err <- err$cov

    krige <-
      if (err[[1]] %in% ZF_PLUS) {
        if (krige[[1]] %in% ZF_PLUS) c(krige, err[-1]) else c(err, list(krige))
      } else {
        if (krige[[1]] %in% ZF_PLUS) c(krige, list(err))
        else list(ZF_SYMBOLS_PLUS, err, krige)
      }   
  }

  simu.grid <- all$grid

  if (nrow(all$given) * all$vdim != length(all$data)) {
    #Print(all)
    stop("dimension of 'given' does not match 'data'") # message
    # return(NA)
  }
 
  txt <- "kriging in space time dimensions>3 where not all the point ly on a grid is not possible yet"
  ## if 4 dimensional then the last coordinates should ly on a grid

  ## now check whether and if so, which of the given points belong to the
  ## points where conditional simulation takes place
  simu <- NULL
  if (simu.grid) {
    ind <- 1 + (t(all$given) - all$xgr[1, ]) / all$xgr[2, ]
    index <- round(ind)
    outside.grid <-
      apply((abs(ind-index)>RFopt$general$gridtolerance) | (index<1) |
            (index > 1 + all$xgr[3, ]), 2, any)   
    if (any(outside.grid)) {
      ## at least some data points are not on the grid:
      ## simulate as if there is no grid
      simu.grid <- FALSE
      ll <- NULL ##  otherwise check will give a warning
      l <- ncol(all$xgr)
      if (l>3) stop(txt)
      xx <- if (l==1) ## dim x locations
                 matrix(seq(from=all$xgr[1], by=all$xgr[2], len=all$xgr[3]),
                        nrow=1)
            else eval(parse(text=paste("t(expand.grid(",
                            paste("seq(from=all$xgr[1,", 1:l,
                                  "], by=all$xgr[2,", 1:l,
                                  "], len=all$xgr[3,", 1:l, "])", collapse=","),
                         "))")))     
      ll <- eval(parse(text=paste("c(",
                   paste("length(seq(from=all$xgr[1,", 1:l,
	                 "], by=all$xgr[2,", 1:l, 
		         "], len=all$xgr[3,", 1:l, "]))",
                         collapse=","),
                   ")")))

      new.index <- rep(0,ncol(index))
      ## data points that are on the grid, must be registered,
      ## so that they can be used as conditioning points of the grid
      if (!all(outside.grid)) {
        new.index[!outside.grid] <- 1 +
          colSums((index[, !outside.grid, drop=FALSE]-1) *
                  cumprod(c(1, ll[-length(ll)])))
      }
      index <- new.index
      new.index <- NULL
    } else {  
      ## data points are all lying on the grid
      simu <- do.call(RFsimulate, args=c(list(x=x, y=y, z=z, T=T,
                                    grid=all$grid,
                                    model=all$model, n=n, 
                                    register=cond.reg,
                                    seed = NA,
                                    spConform=FALSE),
                                    dots))
      ## for all the other cases of simulation see, below
      index <- t(index)
    }
  } else { ## not a grid
    xx <- t(all$x)  ## dim x locations
   
    ## the next step can be pretty time consuming!!!
    ## to.do: programme it in C
    ##
    ## identification of the points that are given twice, as points to
    ## be simulated and as data points (up to a tolerance distance !)
    ## this is important in case of nugget effect, since otherwise
    ## the user will be surprised not to get the value of the data at
    ## that point
    one2ncol.xx <- 1:ncol(xx)
    index <- apply(all$given, 1, function(u){
      i <- one2ncol.xx[colSums(abs(xx - u)) < RFopt$general$gridtolerance]
      if (length(i)==0) return(0)
      if (length(i)==1) return(i)
      return(NA)
    })
  }

  
  if (!simu.grid) {
    ## otherwise the simulation has already been performed (see above)
    tol <- RFopt$general$gridtolerance * nrow(xx)
    if (any(is.na(index)))
      stop("identification of the given data points is not unique - `tol' too large?")
    if (any(notfound <- (index==0))) {
      index[notfound] <- (ncol(xx) + 1) : (ncol(xx) + sum(notfound))
    }

#    Print(t(xx), all$given[notfound, , drop=FALSE])

    xx <- rbind(t(xx), all$given[notfound, , drop=FALSE])
    simu <- do.call(RFsimulate,
                    args=c(list(x=xx, grid=FALSE, model=all$model, n=n,
                      register = cond.reg, seed = NA,  spConform=FALSE), dots))

    xx <- NULL
  } else stopifnot(!is.null(simu))
  
  if (is.null(simu)) stop("random field simulation failed")

   
  if (n==1) {
     
    ## simu values at the `given' locations
    simu.given <- simu[index]

    simu <- as.vector(simu[1:all$restotal]) # as.vector is necessary !! Otherwise
    ##                          is not recognized as a vector
     
  } else {
    
    ## this is a bit more complicated since index is either a vector or
    ## a matrix of dimension  base::dim(simu)-1
    simu.given <-
      matrix(ncol=n, apply(simu, length(base::dim(simu)), function(m) m[index]))
    simu <- apply(simu, length(base::dim(simu)), function(m) m[1:all$restotal])
    simu <- as.vector(simu)
   }

  if (!is.null(err.model)) {
      error <- do.call(RFsimulate, args=c(list(model=err.model, all$given,
                                     grid=FALSE, n=n,
                                     register = RFopt$registers$errregister,
                                     seed = NA,
                                     spConform=FALSE),
                                     dots))
    if (is.null(error)) stop("error field simulation failed")
     simu.given <- simu.given + as.vector(error)
     error <- NULL
   }


  ## to do: als Naeherung bei UK, OK:
  ## kriging(data, method="A") + simu - kriging(simu, method="O") !

  simu <- simu + RFinterpolate(krige.method="O",
                               x=if (!all$missing.x) x else all$x,
                               y=if (!all$missing.x) y else all$y,
                               z=if (!all$missing.x) z else all$z,
                               T=if (!all$missing.x) T else all$T, 
                               grid=if (!all$missing.x) grid else all$grid, 
                               model=krige,
                               data= cbind(all$given,
                                 data=as.vector(all$data)-simu.given),
                               spConform=FALSE)
 
  if (!is.null(all$trend)) {
    simu <- simu + RFsimulate(model=all$trend, x=x, y=y,
                              z=z, T=T, grid=all$grid, n=n,
                              register=interpol.reg,
                              seed = NA,
                              spConform=FALSE)
  }


 return(list(simu=simu, coord.names=all$coord.names, vdim=all$vdim, x=all$x,
             fullgiven=all$fullgiven, fulldata=all$fulldata,
             data.na=all$data.na, data.col=all$data.col))
  
}


