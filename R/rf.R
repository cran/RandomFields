# source("rf.R")
# see getNset.R for the variable .methods


xRFloglikelihood <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                             data, distances, dim, likelihood, ...) {
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
#  all <- Check*Data(model=model, x=x, y=y, z=z, T=T,
 #                  distances=distances, grid=grid, data=data, ...)
  
  model <- list("loglikelihood",  all$model, data=double(all$fulldata),
                len=length(all$fulldata))

  #Print(all, model)

  rfInit(model=model, x=as.double(t(all$fullgiven)),
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

RFcalc <- function(model) {
  if (is.numeric(model)) return(model)
  RFfctn(model, 0,
         coord_system="cartesian", new_coord_system="keep", spConform = FALSE)
}


######################################################################
######################################################################


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
  new <- CheckXT(x, y, z, T, grid=grid, distances=distances, dim=spdim,
                 y.ok=!dosimulate)

  ts.xdim <- as.integer(ncol(new$x) + new$Zeit)

#  if (missing(reg)) stop("'rfInit' may not be called by the user")
 # userdefined <- GetParameterModelUser(model)


  #Print(model);

  ##  Print(userdefined)
  vdim <- .Call("Init",
                as.integer(reg),
                model, 
                if (new$grid) new$x else t(new$x),
                if (is.null(new$y)) {
                  double(0)
                } else if (new$grid) new$y else t(new$y),
                as.double(new$T),
                as.integer(new$spatialdim),
                as.logical(new$grid),
                as.logical(new$distances),
                as.logical(new$Zeit),
                NAOK=TRUE, # ok
                PACKAGE="RandomFields")
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
    stopifnot(missing(distances) || is.null(distances))
    res <- switch(GetProcessType(model),
                  RPgauss = 
                  rfCondGauss(model=model.orig, x=x, y=y, z=z, T=T,
                              grid=grid, n=n, data=data, given=given,
                              err.model=err.model, ...),
                  stop(GetProcessType(model),
                       ": conditional simulation of the process not programmed yet")
                  )

    ## if (missing.x==TRUE), result can already be returned here
    ## but Check Data has to be called again to get the full vector
    ## of given locations and the given data
    if (missing.x <- missing(x)) {

#      Print(data=data, res=res, n=n,
#                            spConform=RFopt$general$spConform)
      res <- FinishImputing(data=data, all=res, n=n,
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


