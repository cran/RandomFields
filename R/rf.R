# source("rf.R")
# see getNset.R for the variable .methods



## sudo apt-get install tcl8.5-dev
## sudo apt-get install tk8.5-dev
## see R Installation and Administration
##./configure --with-tcl-config=/usr/lib/tcl8.5/tclConfig.sh --with-tk-config=/usr/lib/tk8.5/tkConfig.sh
# sudo tar xvf ~/TMP/bwidget-1.9.5.tar 


# readline <- function(...) return("") 

.onAttach <- function (lib, pkg) {
  #packageStartupMessage("\n")
}

#.onUnload <- function(lib, pkg){
#  }
#Implementierung von Cox & Isham's non-separable model

#individuelle Praeferenzliste:
#  1) erste praeferenz
#  2) if # pkte < 10 / 50 -> Gauss??
#  4) if nugget/variance > 0.01 -> CE (and not spectral, e.g.)
#  3) if C(maxh)/C(0) > 0.05  -> TBM else CE


.onUnload <- function(lib, pkg){
  RFoptions(storing=FALSE) ## delete everything
}
 

RFdistr <- function(model, x, q, p, n, dim=1, ...) {
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x
  RFoptOld <- internal.rfoptions(..., RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  model<- list("Distr", PrepareModel2(model), dim=as.integer(dim));
  if (!missing(x)) model$x <- x 
  if (!missing(q)) model$q <- q 
  if (!missing(p)) model$p <- p
  if (!missing(n)) model$n <- n

  
  rfInitSimulation(model=model, x=matrix(0, ncol=dim, nrow=1),
                   y=NULL, z=NULL, T=NULL, grid=TRUE,
                   reg = MODEL.USER,
                   dosimulate=FALSE)
  result <- .Call("EvaluateModel", double(0), as.integer(MODEL.USER),
                  PACKAGE="RandomFields")
}
RFddistr <- function(model, x, dim=1,...) RFdistr(model=model, x=x, dim=dim,...)
RFpdistr <- function(model, q, dim=1,...) RFdistr(model=model, q=q, dim=dim,...)
RFqdistr <- function(model, p, dim=1,...) RFdistr(model=model, p=p, dim=dim,...)
RFrdistr <- function(model, n, dim=1,...) RFdistr(model=model, n=n, dim=dim,...)



RFcov <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
                  distances, dim, 
#                  dim = ifelse(is.matrix(x), ncol(x), 1),
                   fctcall=c("Cov", "CovMatrix", "Variogram",
                           "Pseudovariogram", "Fctn"), ...) {

  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x
  RFoptOld <- internal.rfoptions(xyz_notation=2*(!is.null(y) && !is.matrix(y)),
                                 ..., RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  fctcall <- match.arg(fctcall)
  
  p <- list(fctcall, PrepareModel2(model));
  
  rfInitSimulation(model=p,
                   x=x, y=y, z=z, T=T, grid=grid,
                   distances=distances, spdim=dim, reg = MODEL.USER,
                   dosimulate=FALSE)

  result <- .Call("EvaluateModel", double(0), as.integer(MODEL.USER),
                  PACKAGE="RandomFields")
 
  return(result)
}

RFcovmatrix <- function(...) {
  RFcov(..., fctcall="CovMatrix")
}

RFvariogram <- function (model, x, y=NULL, ...)
  RFcov(model=model, x=x, y, fctcall="Variogram", ...)


RFpseudovariogram <- function(model, x, y=NULL, ...)
  RFcov(model=model, x=x, y, fctcall="Pseudovariogram", ...)

RFfctn <- function(...) {
  RFcov(..., fctcall="EvalFctn")
}


######################################################################
######################################################################


prepare4RFspDataFrame <- function(model.orig=NULL,
                                  info, x, y, z, T, grid, data, RFopt) {
  
  vdim <- info$vdim
  locinfo <- info$loc
  
  ## names 
  if (missing(data))
    names <- rfPrepareNames(model=model.orig, locinfo=locinfo)
  else
    names <- rfPrepareData(model.orig, x, y, z, T, grid, data,
                           cathegories=NULL, names.only=TRUE)
   
  if (!is.null(names$variab.names)) {
    if (vdim != length(names$variab.names))
      stop(paste("you passed a formula for 'model' with left-hand side '",
                 paste(names$variab.names, collapse=","),
                 "', but vdim of the model equals ", vdim, sep=""))
  }
  
  coord.names.incl.T <- names$coord.names
  
  ## coords or GridTopology 
  if (locinfo$grid) { ## grid=TRUE
    coords <- NULL
    xgr <- cbind(locinfo$xgr, locinfo$T)
    colnames(xgr) <- coord.names.incl.T
    xgr[is.na(xgr)] <- 0
    gridTopology <- GridTopology(xgr[1, ], xgr[2, ], xgr[3, ])
  } else { ## grid == FALSE
    gridTopology <- NULL
    
    # cbind of locations from x-matrix and T (if given)
    coords <- apply(t(locinfo$x), 2, rep,
                    times=(locinfo$totpts/locinfo$spatialtotpts))
    if (locinfo$Time) {
      T <- locinfo$T
      coord <- cbind(coord, rep(seq(T[1], by=T[2], len=T[3]),
                                each=locinfo$spatialtotpts))
    }
    colnames(coords) <- coord.names.incl.T
  }

#  Print(RFopt)
  if (RFopt$general$printlevel>0 && RFopt$warn$newstyle) {
    RFoptions(warn.newstyle = FALSE)
    message("New output format of RFsimulate: object of class 'RFsp';\n",
            "for conventional array format use 'RFoptions(spConform=FALSE)'.")
  }

  return(list(coords=coords, gridTopology=gridTopology, vdim=vdim, names=names))
}


rfDoSimulate <- function(n = 1, reg, spConform) {
  stopifnot(length(n) == 1, n>0, is.finite(n))
  RFopt <- RFoptions()
  if (missing(spConform)) spConform <- RFopt$general$spConform

  if (RFopt$gauss$paired && (n %% 2 != 0))
    stop("if paired, then n must be an even number")

  info <- RFgetModelInfo(RFopt$general$register)

  # Print(info); str(info$loc); xxx

  len <- info$loc$len
  vdim <- info$vdim
  total <- info$loc$totpts
  if (is.null(total) || total <= 0)
    stop("register ", RFopt$general$register, " does not look initialized")
 
  error <- integer(1)
#  result <- double(total * n * vdim)

 # Print(info,  total, vdim, len, info$simu$distr, n);  xxx

  result <- .Call("EvaluateModel", as.double(n), as.integer(reg), #userdefined,
                  PACKAGE="RandomFields", DUP=FALSE)
  
  if (!spConform) return(result)
  
  prep <- prepare4RFspDataFrame(model.orig=NULL, info=info, RFopt=RFopt)
  attributes(result)$variab.names <- prep$names$variab.names
  
  res2 <- conventional2RFspDataFrame(result, coords=prep$coords,
                                     gridTopology=prep$gridTopology,
                                     n=n, vdim=prep$vdim)
  return(res2)
}



rfInitSimulation <- function(model, x, y = NULL, z = NULL, T=NULL, grid=FALSE,
                             distances, spdim,
                             reg, dosimulate=TRUE) {

  stopifnot(xor(missing(x), missing(distances)))
  neu <- CheckXT(x, y, z, T, grid=grid, distances=distances, spdim=spdim,
                 y.ok=!dosimulate)
  ts.xdim <- as.integer(ncol(neu$x) + neu$Time)

  if (missing(reg)) stop("rfInitSimulation may not be called by the user")
 # userdefined <- GetParameterModelUser(model)

#  Print(userdefined)

  ## 
  vdim <-
    .Call("Init",
          as.integer(reg),
          model, 
          if (neu$grid) neu$x else t(neu$x),
          if (is.null(neu$y)) {
            double(0)
          } else if (neu$grid) neu$y else t(neu$y),
          as.double(neu$T),
          as.integer(neu$spacedim),
          as.logical(neu$grid),
          as.logical(neu$distances),
          as.logical(neu$Time),
          NAOK=TRUE,
          PACKAGE="RandomFields")  
}



    
RFsimulate <- function (model, x, y = NULL, z = NULL, T = NULL, grid,
                        data, distances, dim, err.model,
                        n = 1,  ...) {
 
  ### preparations #########################################################  
  RFoptOld <- internal.rfoptions(xyz_notation=!is.null(y),
                                 expected_number_simu=n, ..., 
                                 RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]] 
  cond.simu <- !missing(data) && !is.null(data)  
  reg <- if (cond.simu) RFopt$general$condregister else RFopt$general$register
  
  ### simulate from stored register ########################################
  mc <- as.list(match.call(expand.dots=FALSE))
  if (length(mc)==1 ||
      length(mc)==2 && !is.null(mc$n) ||
      length(mc)==3 && !is.null(mc$n) && "..." %in% names(mc)) {
    if (cond.simu) {
      stop("repeated simulation of conditional simulation not programmed yet")
    } else {
      # userdefined <- GetParameterModelUser(PrepareModel2(model, ...))
      return(rfDoSimulate(n=n, reg=reg, spConform=RFopt$general$spConform
                          #userdefined=userdefined
                          ))
    }
  }
  
  ### preparations #########################################################
  stopifnot(!missing(model) && !is.null(model))

  model.orig <- model
  model <- PrepareModel2(model, ...)
  err.model <-  if (missing(err.model)) NULL else PrepareModel2(err.model, ...)
 
  ### conditional simulation ###############################################
  
  if (cond.simu) {
    #Print(RFoptions()$general)
    if (missing.x <- missing(x)) {
      if (is(data, "RFsp")) {
        if (!missing(dim)) 
          stop("'dim' may not be given when 'data' is an RFsp object")
        x.tmp <- coordinates(data)
      } else { ## data not RFsp
        start <- data.starting.col(data, halt=FALSE)
        if (missing(dim)) {
          if (ncol(data) == 2) dim <- 1
          else if (start<0) stop("Please give 'dim' explicitely since x is missing and data is not an RFsp object.")
          else dim <- start - 1
        } else {  
          if (start < 0) colnames(data)[dim+1] <- "data"
          else if (start - 1 != dim) stop(paste("value of 'dim' does not match the column name convention of the data: the name 'data' is found in position", start, "which does not equal 'dim'+1 (dim=", dim, ")"))
        }
        x.tmp <- data[, 1:dim]
      }   
    } else {
      x.tmp <- x
    }

    p <- list("Simulate", model);
    rfInitSimulation(model=p, x=x.tmp, y=y, z=z, T=T, grid=grid,
                                    distances=distances, spdim=dim,
                                    reg=reg, dosimulate=-1)

    info <- RFgetModelInfo(reg)
    grid <- info$loc$grid

#    xxxx
#    Print(info$type, info)
    
    res <- switch(info$role,
                  Gauss = 
                  rfCondGauss(model=model.orig, x=x, y=y, z=z, T=T,
                              grid=grid, n=n, data=data,
                              err.model=err.model, ...),
                  BrownResnick = stop(paste("conditional simulation for ",
                    "BRprocesses not programmed yet.")),
                  Smith =  stop(paste("conditional simulation for Smith ",
                    "processes not programmed yet.")),
                  Schlather =  stop(paste("conditional simulation for ",
                    "Schlather processes not programmed yet.")),
                  Poisson =  stop(paste("conditional simulation for Poisson ",
                    "processes not programmed yet.")),
                  stop(paste(info$role, "not programmed yet"))
                  )

    ## if (missing.x==TRUE), result can already be returned here
    ## but rfPrepareData has to be called again to get the full vector
    ## of given locations and the given data
    if (missing.x) {
      all <- rfPrepareData(model.orig=model.orig,
                           x, y, z, T, grid, data, cathegories=NULL)
      
      coords <- all$x 
      colnames(coords) <- all$coord.names
      
      tmp.res <- conventional2RFspDataFrame(data=res, coords=coords,
                                            gridTopology=NULL,
                                            n=n, vdim=info$vdim)
      res <- tmp.res@data
      tmp.res@coords <- all$fullgiven
      tmp.res@data <- as.data.frame(all$fulldata)
      tmp.res@data[all$data.na, ] <- res
      
      if (is(tmp.res, "RFspatialPointsDataFrame"))
        try(tmp.res <- as(tmp.res, "RFspatialGridDataFrame"), silent=TRUE)
      if (is(tmp.res, "RFpointsDataFrame"))
        try(tmp.res <- as(tmp.res, "RFgridDataFrame"), silent=TRUE)
      
      return(if (RFopt$general$spConform) tmp.res else tmp.res@data)
    }
    ## end of cond simu 

    
  } else { 

    ## unconditional simulation #############################################
    if(!is.null(err.model))
      warning("model for measurement error is ignored in unconditional simulation")
    stat <-  !missing(distances) || !is.matrix(x) || is.null(y)
   
    p <- list("Simulate", model);
 
    #Print(unix.time(
    rfInitSimulation(model=p, x=x, y=y, z=z, T=T, grid=grid,
                     distances=distances, spdim=dim, 
                     reg=reg)#)) #$userdefined
    info <- RFgetModelInfo(reg)
  
    #Print(info$simu)
    #
    #  Print(info);xxxxxx

    ## perform simulation ###################################################
    #Print(unix.time(
    if (n < 1) return(NULL)
    res <- rfDoSimulate(n=n, reg=reg, spConform=FALSE)#))
  } # end of uncond simu


  ## output: either conventional or RFsp   #################################

  if (!RFopt$general$spConform)
    return(res)

  
  prep <- prepare4RFspDataFrame(model.orig, info, x, y, z, T, grid, data, RFopt)
  attributes(res)$variab.names <- prep$names$variab.names

#  Print(res, coords=prep$coords,
#                                     gridTopology=prep$gridTopology,
#                                     n=n, vdim=prep$vdim)

  res2 <- conventional2RFspDataFrame(res, coords=prep$coords,
                                     gridTopology=prep$gridTopology,
                                     n=n, vdim=prep$vdim)
  return(res2)
}




### ist keine Methode im engeren Sinne. Habe ich aus Methods-RFsp.R
### rausgenommen, da bei jeglicher Aenderung in Methods-RFsp.R ich
### komplett neu installieren muss. Bei rf.R muss ich es nicht.
conventional2RFspDataFrame <-
  function(data, coords=NULL, gridTopology=NULL, n=1, vdim=1){

    # stopifnot(!is.na(data))
 
  #   Print(data, coords, gridTopology, n, vdim)
    
    if (!xor(is.null(coords), is.null(gridTopology)))
      stop("one and only one of 'coords' and 'gridTopology' must be NULL")

    variab.names <- attributes(data)$variab.names
    ## may be NULL, if called from RFsimulate, the left hand side of model, if
    ## model is a formula, is passed to 'variab.names'
    attributes(data)$variab.names <- NULL

    #Print(data, coords, gridTopology, n, vdim)
    
    ## grid case
    if (length(coords) == 0) {# war is.null(coords) -- erfasst coords=list() nicht
      grid <- convert2GridTopology(gridTopology) 
      timespacedim <- length(grid@cells.dim)

  
   #   Print(data, dim(data))
      ## naechste Zeile eingefuegt !! und (Martin 30.6.13) wieder
      ## auskommentiert. s. Bsp in RFsimulate
     # if (!is.null(dim(data)) && all(dim(data)[-1]==1)) data <- as.vector(data)
     
      if (is.null(dim(data))) {
      #Print(1, timespacedim, n, vdim, data, coords)
        stopifnot(1 == timespacedim + (n > 1) + (vdim > 1))
      } else {
        if (length(dim(data)) != timespacedim + (n>1) + (vdim > 1)){          
           stop(paste(length(dim(data)),
                     "= length(dim(data)) != timespacedim + (n>1) + (vdim>1) =",
                      timespacedim, '+', (n>1), '+', (vdim > 1)))
        }
      } 
     
 
      if (vdim>1){
        # new order of dimensions: space-time-dims, vdim, n
        perm <- c( 1+(1:timespacedim), 1, if (n>1) timespacedim+2 else NULL)
        data <- aperm(data, perm=perm)
      }
      if (timespacedim==1)
        call <- "RFgridDataFrame"
      else {
        data <- reflection(data, 2, drop=FALSE)
        call <- "RFspatialGridDataFrame"
      }


     }
    
    ## coords case
    if (is.null(gridTopology)){
      #stopifnot(is.matrix(coords))
      if (vdim>1){
        n.dims <- length(dim(data))
        perm <- c(2:(n.dims - (n>1)), 1, if (n>1) n.dims else NULL)
        data <- aperm(data, perm=perm)
      }
      if (is.null(dim(coords)) || ncol(coords)==1)
        call <- "RFpointsDataFrame"
      else call <- "RFspatialPointsDataFrame"
    }

    ## in both cases:
    dim(data) <- NULL
    data <- as.data.frame(matrix(data, ncol=n*vdim))
    
    if (is.null(variab.names))
      variab.names <- paste("variable", 1:vdim, sep="")
    if (length(variab.names) == n*vdim)
      names(data) <- variab.names
    else
      if (length(variab.names) == vdim)
        names(data) <- paste(rep(variab.names, times=n),
                             if (n>1) ".n", if (n>1) rep(1:n, each=vdim),sep="")
      else names(data) <- NULL

    ##  Print(variab.names, names(data), coords)
    
    if (is.null(coords)){
      do.call(call, args=list(data=data, grid=grid,
                      RFparams=list(n=n, vdim=vdim)))
    } else {
      do.call(call, args=list(data=data, coords=coords,
                      RFparams=list(n=n, vdim=vdim)))
    }
  }


list2RMmodelExt <- function(x, isRFsp=FALSE,
                            coords, gridTopology, data.RFparams) {
  stopifnot(is.list(x),
            all(c("model", "trend", "ml.value", "residuals") %in% names(x)))

  if (isRFsp) {
    stopifnot(!missing(coords) &&
              !missing(gridTopology) &&
              !missing(data.RFparams))

    ## convert residuals to RFsp class
    err <-
      try({
        lres <- length(x$residuals)
        if (lres > 0) {
          for (i in 1:lres) {
            gT <- if (length(gridTopology) < i) NULL else gridTopology[[i]]
            co <- if (length(coords)<i) NULL else coords[[i]]
            if (!is.null(x$residuals[[i]]))
                x$residuals[[i]]<-
                  conventional2RFspDataFrame(data=x$residuals[[i]],
                                             coords=co,
                                             gridTopology=gT,
                                             n=data.RFparams[[i]]$n,
                                             vdim=data.RFparams[[i]]$vdim)
          }
        }
      }
          , silent=!TRUE)
    if(class(err)=="try-error")
      warning(paste("residuals could not be coerced to class'RFsp';",
                    err))
  }
    
  return(new(ZF_MODELEXT,
             list2RMmodel(x$model),
             trend = list2RMmodel(x$trend),
             likelihood = x$ml.value,
             residuals = x$residuals)) 
}
  



