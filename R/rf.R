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
 


xRFloglikelihood <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL, data,
                            distances, dim, likelihood, ...) {
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x
  relax <- isFormulaModel(model)
  RFoptOld <-
    if (missing(likelihood)) internal.rfoptions(..., RELAX=relax)
    else internal.rfoptions(likelihood=likelihood, ..., RELAX=relax)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  Reg <- MODEL_USER
  all <- rfPrepare(model=model, x=x, y=y, z=z, T=T,
                   distances=distances, grid=grid, data=data,
                   reg = Reg, fillall = FALSE, names.only=FALSE, ...)
  
  model <- list("loglikelihood",  all$model, data=double(all$data),
                len=length(all$data))

  #Print(all, model)

  rfInit(model=model, x=as.double(t(all$given)),
         y=NULL, z=NULL, T=all$T, grid=all$grid,
         reg = MODEL_USER, dosimulate=FALSE, old.seed=RFoptOld[[1]]$general$seed)
  
  return(.Call("EvaluateModel", double(0),  Reg, PACKAGE="RandomFields"))
}

rfdistr <- function(model, x, q, p, n, dim=1, ...) {
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x

  RFoptOld <- internal.rfoptions(..., RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  model<- list("Distr", PrepareModel2(model), dim=as.integer(dim));
  if (!missing(x)) {
    model$x <- if (is.matrix(x)) t(x) else x
  }
  if (!missing(q)) {
    model$q <- if (is.matrix(q)) t(q) else q
  }
  if (!missing(p)) {
    model$p <- if (is.matrix(p)) t(p) else p
  }
  if (!missing(n)) {
    model$n <- n
  }

#  Print(model)
  
  rfInit(model=model, x=matrix(0, ncol=dim, nrow=1),
         y=NULL, z=NULL, T=NULL, grid=FALSE, reg = MODEL_USER,
         dosimulate=FALSE, old.seed=RFoptOld[[1]]$general$seed)

  res <-  .Call("EvaluateModel", double(0), as.integer(MODEL_USER),
                PACKAGE="RandomFields")

  if (RFoptOld[[2]]$general$returncall) attr(res, "call") <-
    as.character(deparse(match.call(call=sys.call(sys.parent()))))
  attr(res, "coord_system") <- c(orig=RFoptOld[[2]]$coords$coord_system,
                                 model=RFoptOld[[2]]$coords$new_coord_system)
  return(res)
}

RFdistr <- function(model, x, q, p, n, dim=1, ...) {
   rfdistr(model=model, x=x, q=q, p=p, n=n, dim=dim, ...)
}
RFddistr <- function(model, x, dim=1,...) {
  if (hasArg("q") || hasArg("p") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, x=x, dim=dim,...)
}
RFpdistr <- function(model, q, dim=1,...) {
  if (hasArg("x") || hasArg("p") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, q=q, dim=dim,...)
}
RFqdistr <- function(model, p, dim=1,...){
  if (hasArg("x") || hasArg("q") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, p=p, dim=dim,...)
}
RFrdistr <- function(model, n, dim=1,...) {
  if (hasArg("x") || hasArg("q") || hasArg("p")) stop("unknown argument(s)");
  rfdistr(model=model, n=n, dim=dim,...)
}



rfeval <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                  distances, dim, ..., 
                  ##                  dim = ifelse(is.matrix(x), ncol(x), 1),
                  fctcall=c("Cov", "CovMatrix", "Variogram",
                    "Pseudovariogram", "Fctn")) {
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x
  RFoptOld <- internal.rfoptions(xyz_notation=2*(!is.null(y) && !is.matrix(y)),
                                 ..., RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  fctcall <- match.arg(fctcall)
  
  p <- list(fctcall, PrepareModel2(model));
 
  rfInit(model=p, x=x, y=y, z=z, T=T, grid=grid,
         distances=distances, spdim=dim, reg = MODEL_USER, dosimulate=FALSE,
         old.seed=RFoptOld[[1]]$general$seed)

  res <- .Call("EvaluateModel", double(0), as.integer(MODEL_USER),
               PACKAGE="RandomFields")

   if (RFoptOld[[2]]$general$returncall) attr(res, "call") <-
    as.character(deparse(match.call(call=sys.call(sys.parent()))))
  attr(res, "coord_system") <- .Call("GetCoordSystem", as.integer(MODEL_USER),
              RFoptOld[[2]]$coords$coord_system,
              RFoptOld[[2]]$coords$new_coord_system)
   return(res)
}


RFcov <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
                        distances, dim, ...) {
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Cov")
}


RFcovmatrix <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
                        distances, dim, ...) {
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="CovMatrix")
}

RFvariogram <- function (model, x, y=NULL,  z = NULL, T=NULL, grid,
                        distances, dim,...){
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Variogram")
}

RFpseudovariogram <- function(model, x, y=NULL,  z = NULL, T=NULL, grid,
                        distances, dim,...){
   rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Pseudovariogram")
}

RFfctn <- function(model, x, y=NULL,  z = NULL, T=NULL, grid,
                   distances, dim,...) {
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Fctn")
}


######################################################################
######################################################################


prepare4RFspDataFrame <- function(model=NULL,
                                  info, x, y, z, T, grid=NULL, data, RFopt) {
  
  vdim <- info$vdim
  locinfo <- info$loc
  
  ## names 
  if (missing(data))
    names <- rfPrepareNames(model=model, locinfo=locinfo)
  else
    names <- rfPrepare(model=model, x=x, y=y, z=z, T=T,
                       grid=grid, data=data,
                       names.only=TRUE)
   
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
    gridTopology <- sp::GridTopology(xgr[1, ], xgr[2, ], xgr[3, ])
  } else { ## grid == FALSE
    gridTopology <- NULL
    
    # cbind of locations from x-matrix and T (if given)
    coords <- as.matrix(apply(t(locinfo$x), 2, rep,
                              times=(locinfo$totpts/locinfo$spatialtotpts)))
    if (locinfo$Time) {
      T <- locinfo$T
      coords <- cbind(coords, rep(seq(T[1], by=T[2], len=T[3]),
                                each=locinfo$spatialtotpts))
    }
    if (is.matrix(coords)) colnames(coords) <- coord.names.incl.T
  }

  if (RFopt$general$printlevel>0 && RFopt$internal$warn_newstyle) {
    RFoptions(internal.warn_newstyle = FALSE)
    message("New output format of RFsimulate: S4 object of class 'RFsp';\n",
            "for a bare, but faster array format use 'RFoptions(spConform=FALSE)'.")
  }

  return(list(coords=coords, gridTopology=gridTopology, vdim=vdim, names=names))
}


rfDoSimulate <- function(n = 1, reg, spConform) {
  stopifnot(length(n) == 1, n>0, is.finite(n))
  RFopt <- RFoptions()
  if (missing(spConform)) spConform <- RFopt$general$spConform

  if (RFopt$gauss$paired && (n %% 2 != 0))
    stop("if paired, then n must be an even number")

  info <- RFgetModelInfo(RFopt$registers$register, level=3)

 
  len <- info$loc$len
  vdim <- info$vdim
  total <- info$loc$totpts
  if (is.null(total) || total <= 0)
    stop("register ", RFopt$registers$register, " does not look initialized")
 
  error <- integer(1)

  result <- .Call("EvaluateModel", as.double(n), as.integer(reg), #userdefined,
                  PACKAGE="RandomFields")
  
  if (!spConform) return(result)
  
  prep <- prepare4RFspDataFrame(model=NULL, info=info, RFopt=RFopt)
  attributes(result)$variab.names <- prep$names$variab.names
  
  res2 <- conventional2RFspDataFrame(result, coords=prep$coords,
                                     gridTopology=prep$gridTopology,
                                     n=n, vdim=prep$vdim, T=info$loc$T,
                                     vdim_close_together=RFopt$general$vdim_close_together)
  return(res2)
}



rfInit <- function(model, x, y = NULL, z = NULL, T=NULL, grid=FALSE,
                   distances, spdim, reg, dosimulate=TRUE, old.seed=NULL) {
  stopifnot(xor(missing(x), #|| length(x)==0,
                missing(distances) || length(distances)==0))
 
  RFopt <- RFoptions() 
  if (!is.na(RFopt$general$seed)) {
    allequal <- all.equal(old.seed, RFopt$general$seed)
    allequal <- is.logical(allequal) && allequal
    if (dosimulate && RFopt$general$printlevel >= PL_IMPORTANT &&
        (is.null(old.seed) || (!is.na(old.seed) && allequal)
         )
        ) {
       message("NOTE: simulation is performed with fixed random seed ",
               RFopt$general$seed,
               ".\nSet 'RFoptions(seed=NA)' to make the seed arbitrary.")
     }
    set.seed(RFopt$general$seed)
  }
  ##  if (missing(x) || length(x) == 0) stop("'x' not given")
  
  neu <- CheckXT(x, y, z, T, grid=grid, distances=distances, spdim=spdim,
                 y.ok=!dosimulate)

  ts.xdim <- as.integer(ncol(neu$x) + neu$Time)

#  if (missing(reg)) stop("'rfInit' may not be called by the user")
 # userdefined <- GetParameterModelUser(model)


  ##  Print(userdefined)
  vdim <- .Call("Init",
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
                NAOK=TRUE, # ok
                PACKAGE="RandomFields")
}


FinishImputing <- function(data, res, n, spConform) {
  
  #Print(data, res, tail(res$simu), spConform);
  if (is(data, "RFsp")) {
    if (spConform) {
      rows <- nrow(data@data)
      simu <- res$simu
      sum.idx <- start <- 0
      for (i in 1:length(data@data)) {          
        idx <- which(res$data.na > start & res$data.na <= start + rows)
        len.idx <- length(idx)
        data@data[res$data.na[idx] - start, i] <- simu[sum.idx + 1:len.idx]
        start <- start + rows
        sum.idx <- sum.idx + len.idx
      }
      
      # Print(data) 
      return(data)
    } else {
      values <- as.matrix(data@data)
      values[is.na(values)] <- res$simu
      
      #print(cbind(coordinates(data), values))
      return(cbind(coordinates(data), values))
    }
  } else { ## not RFsp
    #Print("for testing")
    coords <- res$fullgiven
    ##  coords <- res$x
    colnames(coords) <- res$coord.names
    
    values <- data[, res$data.col]
    values[is.na(values)] <- res$simu

    #Print(values)
    
    if (!spConform)  return(cbind(coords, values))
    
    tmp.res <- conventional2RFspDataFrame(data=values, coords=coords,
                                          gridTopology=NULL,
                                          n=n, vdim=res$vdim,
                                          vdim_close_together=FALSE)
    res <- tmp.res@data         
    if (is(tmp.res, "RFspatialPointsDataFrame"))
      try(tmp.res <- as(tmp.res, "RFspatialGridDataFrame"), silent=TRUE)
    if (is(tmp.res, "RFpointsDataFrame"))
      try(tmp.res <- as(tmp.res, "RFgridDataFrame"), silent=TRUE)
    
    return(tmp.res)
  }
}



## todo: finish programming
.RFprobab <- function(model, x, y = NULL, z = NULL, T = NULL, grid,
             t, log=FALSE, ...) {
  .RFprobabdens(model=model, x=x, y = y, z = z, T = T, grid=grid,
               t=t, log=log, ..., type="Probab")
}

.RFdensity <- function(model, x, y = NULL, z = NULL, T = NULL, grid,
             t, log=FALSE, ...) {
  .RFprobabdens(model=model, x=x, y = y, z = z, T = T, grid=grid,
               t=t, log=log, ..., type="Density")
}


.RFprobabdens <- function(model, x, y = NULL, z = NULL, T = NULL, grid=NULL,
             t, log=FALSE, ..., type) {
#  prepareOptions(environment())
#  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFoptOld <- internal.rfoptions(xyz_notation=!is.null(y),
                                 RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  
  reg <- MODEL_UNUSED

 ### preparations #########################################################
  stopifnot(!missing(model) && !is.null(model))
  model <- PrepareModel2(model, ...)

  rfInit(model=list(type, log=log,
             setseed=eval(parse(text="quote(set.seed(seed=seed))")),
             env=.GlobalEnv, model), x=x, y=y, z=z, T=T,
         grid=grid,  reg=reg,
         old.seed=RFoptOld[[1]]$general$seed)

  

  result <- .Call("EvaluateModel", as.double(t), reg, #userdefined,
                  PACKAGE="RandomFields")
 
  return(result)
}

    
RFsimulate <- function (model, x, y = NULL, z = NULL, T = NULL, grid=NULL,
                        distances, dim, data, given = NULL, err.model,
                        n = 1,  ...) {

  mc <- as.character(deparse(match.call()))

### preparations #########################################################  
  if (!missing(distances) && length(distances)  > 0)
    RFoptOld <- internal.rfoptions(xyz_notation=!is.null(y),
                                   expected_number_simu=n, ..., 
                                   general.spConform = FALSE,
                                   RELAX=isFormulaModel(model))
  else 
    RFoptOld <- internal.rfoptions(xyz_notation=!is.null(y),
                                   expected_number_simu=n, ...,
                                   RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  
  cond.simu <- !missing(data) && !is.null(data)  
  reg <- if (cond.simu) MODEL_COND else RFopt$registers$register

  ### simulate from stored register ########################################
  mcall <- as.list(match.call(expand.dots=FALSE))
  if (length(mcall)==1 ||
      length(mcall)==2 && !is.null(mcall$n) ||
      length(mcall)==3 && !is.null(mcall$n) && "..." %in% names(mcall)) {
    if (cond.simu) {
      stop("repeated simulation of conditional simulation not programmed yet")
    } else {
      # userdefined <- GetParameterModelUser(PrepareModel2(model, ...))
      res <- rfDoSimulate(n=n, reg=reg, spConform=RFopt$general$spConform
                          #userdefined=userdefined
                          )
      if (RFopt$general$returncall) attr(res, "call") <- mc
      attr(res, "coord_system") <- .Call("GetCoordSystem", reg,
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
      return(res)
    }
  }
  
  ### preparations #########################################################
  stopifnot(!missing(model) && !is.null(model))

  model.orig <- model
  model <- PrepareModel2(model, ...)
  err.model <-  if (missing(err.model)) NULL else PrepareModel2(err.model, ...)
 
  ### conditional simulation ###############################################
  #Print(cond.simu)
  
  if (cond.simu) {
    #Print(RFoptions()$general)
    if (isSpObj(data)) data <- sp2RF(data)
    res <- switch(GetProcessType(model),
                  RPgauss = 
                  rfCondGauss(model=model.orig, x=x, y=y, z=z, T=T,
                              grid=grid, n=n, data=data, given=given,
                              err.model=err.model, ...),
                  stop(GetProcessType(model),
                       ": conditional simulation of the process not programmed yet")
                  )

    ## if (missing.x==TRUE), result can already be returned here
    ## but rfPrepare has to be called again to get the full vector
    ## of given locations and the given data
    if (missing.x <- missing(x)) {

#      Print(data=data, res=res, n=n,
#                            spConform=RFopt$general$spConform)
      res <- FinishImputing(data=data, res=res, n=n,
                            spConform=RFopt$general$spConform)
      if (RFopt$general$returncall) attr(res, "call") <- mc
      attr(res, "coord_system") <- .Call("GetCoordSystem",
                                     as.integer(reg),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
      return(res)      
    } else {
      res <- res$simu
      rfInit(model=list("Simulate",
               setseed=eval(parse(text="quote(set.seed(seed=seed))")),
               env=.GlobalEnv, model), x=x, y=y, z=z, T=T,
           grid=grid, distances=distances, spdim=dim, reg=reg,
           old.seed=RFoptOld[[1]]$general$seed)
    }        ## end of cond simu

  } else { ## unconditional simulation ####
    if(!is.null(err.model))
      warning("model for measurement error is ignored in unconditional simulation")

    rfInit(model=list("Simulate",
               setseed=eval(parse(text="quote(set.seed(seed=seed))")),
               env=.GlobalEnv, model), x=x, y=y, z=z, T=T,
           grid=grid, distances=distances, spdim=dim, reg=reg,
           old.seed=RFoptOld[[1]]$general$seed)

    if (n < 1) return(NULL)
    
    res <- rfDoSimulate(n=n, reg=reg, spConform=FALSE)    
  } # end of uncond simu
  info <- RFgetModelInfo(reg, level=3)


  ## output: either conventional or RFsp   #################################
  
  if (RFopt$general$returncall) attr(res, "call") <- mc
  attr(res, "coord_system") <- .Call("GetCoordSystem",
                                     as.integer(reg),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
  if (!RFopt$general$spConform) return(res)
  if (length(res) > 1e7) {
    message("Too big data set (more than 1e7 entries) to allow for 'spConform=TRUE'. So the data are returned as if 'spConform=FALSE'")
    return(res)
  }

  
  prep <- prepare4RFspDataFrame(model.orig, info, x, y, z, T,
                                grid, data, RFopt)
  attributes(res)$variab.names <- prep$names$variab.names
 
  res2 <- conventional2RFspDataFrame(data=res, coords=prep$coords,
                                      gridTopology=prep$gridTopology,
                                      n=n, vdim=prep$vdim,
                                      T=info$loc$T,
                                      vdim_close_together=
                                      RFopt$general$vdim_close_together)

      
  if (is.raster(x)) {
    res2 <- raster::raster(res2)
    raster::projection(res2) <- raster::projection(x)
  }

  
  if (RFopt$general$returncall) attr(res2, "call") <- mc
  attr(res, "coord_system") <- .Call("GetCoordSystem",
                                     as.integer(reg),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
  return(res2)
}




### ist keine Methode im engeren Sinne. Habe ich aus Methods-RFsp.R
### rausgenommen, da bei jeglicher Aenderung in Methods-RFsp.R ich
### komplett neu installieren muss. Bei rf.R muss ich es nicht.
conventional2RFspDataFrame <-
  function(data, coords=NULL, gridTopology=NULL, n=1, vdim=1, T=NULL,
           vdim_close_together) {
  
  
  if (!xor(is.null(coords), is.null(gridTopology)))
    stop("one and only one of 'coords' and 'gridTopology' must be NULL")
  
  variab.names <- attributes(data)$variab.names
  ## may be NULL, if called from 'RFsimulate', the left hand side of model, if
  ## model is a formula, is passed to 'variab.names'
  attributes(data)$variab.names <- NULL
  
  ## grid case
  if (length(coords) == 0) {# war is.null(coords) -- erfasst coords=list() nicht
    grid <- convert2GridTopology(gridTopology) 
    timespacedim <- length(grid@cells.dim)
    
    
    
    ## naechste Zeile eingefuegt !! und (Martin 30.6.13) wieder
    ## auskommentiert. s. Bsp in 'RFsimulate'
    ## if (!is.null(dim(data)) && all(dim(data)[-1]==1)) data <- as.vector(data)
     
    if (is.null(dim(data))) {
      stopifnot(1 == timespacedim + (n > 1) + (vdim > 1))
    } else {
     
      if (length(dim(data)) != timespacedim + (n>1) + (vdim > 1)){          
        stop(paste(length(dim(data)),
                   "= length(dim(data)) != timespacedim + (n>1) + (vdim>1) =",
                   timespacedim, '+', (n>1), '+', (vdim > 1)))
      }
    } 
        
    if (vdim>1 && vdim_close_together){
      ## new order of dimensions: space-time-dims, vdim, n
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
    if (vdim>1 && vdim_close_together){
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
                    RFparams=list(n=n, vdim=vdim, T=T)))
  } else {
    #Print(call, args=list(data=data, coords=coords,RFparams=list(n=n, vdim=vdim)))
     do.call(call, args=list(data=data, coords=coords,
                    RFparams=list(n=n, vdim=vdim, T=T)))
  }
}




