#getusermodel <- function(info, range, cov, derivative) {
#  ## hypermodels not allowed!
#  stopifnot(!is.null(info), !is.null(range), is.null(check))
#  .C("Init")
#  if (!is.null(cov))
#    .C("AddCov")
#  if (!is.null(derivative))
#    .C("AddTBM")
#}


GetRegisterInfo <- function(register=0, ignore.active=FALSE, max.elements=10^6) {
  # ignore.active=TRUE only for internal debugging information!
  if (ignore.active)
     warning("ignore.active = TRUE may cause failure of R -- do not use it!!")
  .Call("GetRegisterInfo", as.integer(register), as.logical(ignore.active),
        as.integer(if (max.elements > .Machine$integer.max)
                   .Machine$integer.max else max.elements),
        PACKAGE="RandomFields")
}

GetModelInfo <- function(register, modelreg, level=3, gatter=FALSE) {
## positive values refer the covariance models in the registers
#define MODEL_USER : 0  /* for user call of Covariance etc. */
#define MODEL_SIMU : 1  /* for GaussRF etc */ 
#define MODEL_INTERN  : 2 /* for kriging, etc; internal call of covariance */
#define MODEL_MLE  : 3
#define MODEL_BOUNDS  4 : - /* MLE, lower, upper */
# [ in RF.h we have the modulus minus 1 ]
  if (!missing(modelreg)) {
    stopifnot(missing(register))
    register <- -modelreg -1
  } else if (missing(register)) register <- 0
  .Call("GetExtModelInfo", as.integer(register), as.integer(level),
        as.integer(gatter), PACKAGE="RandomFields")
}
  

GetModel <- function(register, modelreg, modus=0) {
#define MODEL_USER 0  /* for user call of Covariance etc. */
#define MODEL_SIMU 1  M/* for GaussRF etc */ 
#define MODEL_INTERN 2 /* for kriging, etc; internal call of covariance */
#define MODEL_MLE  3 
#define MODEL_BOUNDS 4 /* MLE, lower, upper */
  if (!missing(modelreg)) {
    stopifnot(missing(register))
    register <- -modelreg -1
  } else if (missing(register)) register <- 0
  .Call("GetModel", as.integer(register), as.integer(modus),
        PACKAGE="RandomFields")
}

GetPracticalRange <- function(model, param, dim=1) {
  vdim <- .Call("CheckModelUser", PrepareModel(model, param)$model,
                dim=as.integer(dim), xdim=as.integer(dim), FALSE,
                PACKAGE="RandomFields")

  natscl <- double(1)
  .C("UserGetNatScaling", natscl, PACKAGE="RandomFields", DUP=FALSE)
  return(1.0 / natscl)
}

#GetPracticalRange("whittle", 1)

GetrfParameters <- function(initcov){
  ## if initcov then InitModelList is called before
  ## values are return (necessary to get covnr right, but should
  ## not be done for RFparameters (debugging reasons)
  maxints <- integer(1)
  .C("GetMaxDims", maxints, DUP=FALSE, PACKAGE="RandomFields")
  
  p <- .C(c("GetrfParameters", "GetrfParametersI")[1 + (initcov != 0) ],
      covmaxchar=integer(1), methodmaxchar=integer(1),
      distrmaxchar=integer(1),
      covnr=integer(1), methodnr=integer(1), distrnr=integer(1),
      maxdim=integer(maxints), maxmodels=integer(1), 
      PACKAGE="RandomFields")
  names(p$maxdim) <- c("cov", "mle", "simu", "ce", "tbm", "mpp", "hyper", "nug",
         "vario")
  return(p)
}

GetMethodNames <- function() {
  assign(".p", GetrfParameters(TRUE))
  l <- character(.p$methodnr)
  for (i in 1:.p$methodnr) {
    l[i] <- .C("GetMethodName", as.integer(i-1),
               n=paste(rep(" ", .p$methodmaxchar), collapse=""),
               PACKAGE="RandomFields")$n
  }
  return(l)
}

GetDistributionNames <- function() {
  assign(".p", GetrfParameters(TRUE))
  l <- character(.p$distrnr)
  for (i in 1:.p$distrnr) {
    l[i] <- .C("GetDistrName", as.integer(i-1),
               n=paste(rep(" ",.p$distrmaxchar), collapse=""),
               PACKAGE="RandomFields")$n
  }
  return(l)
}


GetModelNames <- function(stationary = c("any", "stationary", "variogram", "irf",
                            "auxmatrix", "auxvector",  "gencovariance",
                            "nonstationary", "genvariogram", "prev.model"),
                          isotropy = c("any", "isotropic", "spaceisotropic",
                                  "zerospaceiso", "anisotropic", "prev.model"),
                          multivariate = -9999:9999, ## any !
                          operator,
                          normalmix,
                          finiterange,
                          methods,
                          internal=FALSE,
                          dim = 1) {
  # vdim koennte man auch noch in die args aufnehmen.

  assign(".p", GetrfParameters(TRUE))
  
  par <- as.character(formals()$stationary)[-1]  
  stationary <- pmatch(stationary, par, duplicates.ok=TRUE) - 1
  stationary <- if (any(stationary == 0)) 0:(length(par) - 2) else stationary - 1

  par <- as.character(formals()$isotropy)[-1]  
  isotropy <- pmatch(isotropy, par, duplicates.ok=TRUE) - 1
  isotropy <- if (any(isotropy == 0)) 0:(length(par) - 2) else isotropy - 1

  mn <- GetMethodNames()  
  methods <-
    if (missing(methods)) 1:length(mn) else pmatch(methods, mn,
                                                   duplicates.ok=TRUE)

  empty <- paste(rep(" ", .p$covmaxchar), collapse="")
  l <- character(.p$covnr)
  for (i in 1:.p$covnr) {
    l[i] <- .C("GetModelName",as.integer(i-1),
               n=empty, PACKAGE="RandomFields")$n
  }

  loperator <- integer(.p$covnr)
  lnormalmix <- integer(.p$covnr)
  lfiniterange <- integer(.p$covnr)
  linternal    <- integer(.p$covnr)
  stats <- integer(.p$covnr)
  isos <- integer(.p$covnr)
  maxcdim <- integer(.p$covnr)
  vdim <- integer(.p$covnr)
  .C("GetAttr", loperator, lnormalmix, lfiniterange, linternal,
     stats, isos, maxcdim, vdim,
     DUP=FALSE, PACKAGE="RandomFields")
  
  idx <- integer(.p$covnr * length(mn))
  .C("GetModelList", idx, as.integer(TRUE), PACKAGE="RandomFields", DUP=FALSE)
  dim(idx) <- c(.p$covnr, length(mn))
  
  if (!(hyp <- missing(operator))) hyp <- loperator - !operator
  if (!(nmix <- missing(normalmix))) nmix <- lnormalmix - !normalmix
  if (!(nfr <- missing(finiterange))) nfr <- lfiniterange - !finiterange
  if (!(nint <- missing(internal))) nint <- linternal - !internal

  return(l[(stats %in% stationary) &(isos %in% isotropy) &
           (vdim %in% multivariate) &
           hyp & nmix & nfr & nint & (maxcdim >= dim) &
           apply(idx, 1, function(x) any(x[methods]))
           ])
}



GetModelList <- function(abbr=TRUE, internal=FALSE) {
  names <- GetModelNames(internal=internal)
  methods <- GetMethodNames()
  if (abbr) methods <- substr(methods, 1, if (is.logical(abbr)) 5 else abbr)
 
  idx <- integer(length(names) * length(methods))
  .C("GetModelList", idx, as.integer(internal), PACKAGE="RandomFields", DUP=FALSE)
  t(matrix(as.logical(idx), ncol=length(names), dimnames=list(methods, names)))
}


RFparameters <- function (..., no.readonly=FALSE) {
  assign(".p", GetrfParameters(FALSE))
  assign(".methods",  if (.p$covnr == -1) 0:100 else GetMethodNames())
  
  ## do not add any temporary variable til ## **
  ## do not remove leading "." from .maxdim
  PracticalRange <- integer(1)
  ## always logical returned
  ## PracticalRange also allows for being set to
  ##  0 : no practical range
  ##  1,11 : practical range, evaluated exactly (if given in RFCovfct.cc)
  ##  2,12 : approximative value (if given in RFCovfct.cc)
  ##  3,13 : rough guess (good enough for MLE) (if given in RFCovfct.cc)
  ## >10: and if nothing appropriate given in RFCovfct.cc then numerical approx.
  
  pch <- as.character(" ")
  PrintLevel <- integer(1)
  Storing <- integer(1)
  skipchecks <- integer(1)
  stationary.only <- integer(1)
  exactness <- integer(1)
  every <- integer(1)

#  ave.r <- double(1)
#  ave.dist <- double(1)
  
  CE.force <- integer(1)
  CE.mmin <- double(.p$maxdim["ce"])
  CE.strategy <- integer(1)
  CE.maxmem <- double(1)
  CE.tolIm <- double(1)
  CE.tolRe <- double(1)
  CE.trials <- integer(1)
  CE.useprimes <- integer(1)
  CE.dependent <- integer(1)
  CE.method <- integer(1)
  
  local.force <- integer(1)
  local.mmin <- double(.p$maxdim["ce"])
  local.maxmem <- double(1)
  local.tolIm <- double(1)
  local.tolRe <- double(1)
  local.useprimes <- integer(1)
  local.dependent <- integer(1)

  direct.bestvariables <- integer(1)
  direct.maxvariables <- integer(1)
  direct.method <- integer(1)
  direct.svdtolerance <- double(1)

  markov.neighbours <- integer(1)
  markov.precision <- double(1)
  markov.cyclic <- integer(1)
  markov.maxmem <- integer(1)
 
  nugget.tol <- double(1)
  nugget.meth <- integer(1)
  
  sequ.max <- integer(1)
  sequ.back <- integer(1)
  sequ.initial <- integer(1)

  ## parameters for special functions
  ## currently none
  ##
  
  spectral.lines <- integer(.p$maxdim["tbm"])
  spectral.grid <- integer(1)
  spectral.ergodic <- integer(1)
  spectral.metro <- integer(1)
  spectral.nmetro <- integer(1)
  spectral.sigma <- double(1)

  TBM.method <- integer(1)
  TBM.center <- double(.p$maxdim["tbm"])
  TBM.points <- integer(1)


  TBM2.lines <- integer(1)
  TBM2.linesimufactor <- double(1)
  TBM2.linesimustep <- double(1)
  TBM2.layers <- integer(1)
#  TBM2.num <- integer(1)
 
  TBM3.lines <- integer(1)
  TBM3.linesimufactor <- double(1)
  TBM3.linesimustep <- double(1)
  TBM3.layers <- integer(1)

  TBMCE.force <- integer(1)
  TBMCE.mmin <- double(.p$maxdim["ce"])
  TBMCE.strategy <- integer(1)
  TBMCE.maxmem <- double(1)
  TBMCE.tolIm <- double(1)
  TBMCE.tolRe <- double(1)
  TBMCE.trials <- integer(1)
  TBMCE.useprimes <- integer(1)
  TBMCE.dependent <- integer(1)
 
  mpp.locations <- integer(1)
  mpp.intensity <- double(.p$maxdim["mpp"])
  mpp.plus  <- double(.p$maxdim["mpp"])
  mpp.relRadius <- double(.p$maxdim["mpp"])
#  mpp.scale <- double(.p$maxdim["mpp"])
  mpp.approxzero <- double(1)
  mpp.samplingdist <- double(1)
  mpp.samplingr <- double(1)
  mpp.p <- double(1)
  mpp.beta <- double(1)

  hyper.superpos <- integer(1)
  hyper.maxlines <- integer(1)
  hyper.mar.distr <- integer(1)
  hyper.mar.param <- double(1)
 
  maxstable.maxGauss <- double(1)

  arg.list <- ls()
  ## **

 
  ## first element is the function name
  parameters <- list(...)
  for (m in 0:1) {
    # m = 0 reading parameters
    # m = 1 storing parameters
    storage.mode(m) <- "integer"
    ## "SetParam" more complicated since pch is of character type
 
    pch <- .Call("SetParamPch", m, pch, PACKAGE="RandomFields");
    .C("SetParam", m, Storing, PrintLevel, PracticalRange, skipchecks,
       every,  PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamDecision", m, stationary.only, exactness, 
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamCircEmbed", m, CE.force, CE.tolRe, CE.tolIm, CE.trials, 
       CE.mmin, CE.useprimes, CE.strategy, CE.maxmem, CE.dependent,
       CE.method,
       PACKAGE="RandomFields", DUP=FALSE)
    
    .C("SetParamLocal", m, local.force, local.tolRe, local.tolIm,
       local.mmin, local.useprimes, local.maxmem,
       local.dependent,
       PACKAGE="RandomFields", DUP=FALSE)

    .C("SetParamDirect", m, direct.method,
       direct.bestvariables, direct.maxvariables, direct.svdtolerance,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamMarkov", m, markov.neighbours, markov.precision, markov.cyclic,
       markov.maxmem,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamNugget", m, nugget.tol, nugget.meth,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamSequential", m, sequ.max, sequ.back, sequ.initial,
       PACKAGE="RandomFields", DUP=FALSE)
   
   .C("SetParamSpectral", m, spectral.lines, spectral.grid, spectral.ergodic,
       spectral.metro, spectral.nmetro, spectral.sigma,
       PACKAGE="RandomFields", DUP=FALSE)

    .C("SetParamTBMCE", m, TBMCE.force, TBMCE.tolRe, TBMCE.tolIm, TBMCE.trials,
       TBMCE.mmin, TBMCE.useprimes, TBMCE.strategy, TBMCE.maxmem,
       TBMCE.dependent,
       PACKAGE="RandomFields", DUP=FALSE)
    
    .C("SetParamTBM2", m, TBM2.lines, TBM2.linesimufactor,
       TBM2.linesimustep, TBM2.layers,
       PACKAGE="RandomFields", DUP=FALSE)

     .C("SetParamTBM3", m, TBM3.lines, TBM3.linesimufactor,
       TBM3.linesimustep, TBM3.layers,
       PACKAGE="RandomFields", DUP=FALSE)

    .C("SetParamTBM", m, TBM.method, TBM.center, TBM.points,
       PACKAGE="RandomFields", DUP=FALSE, NAOK=TRUE)
    .C("SetParamMPP", m, mpp.locations, mpp.intensity,
       mpp.plus, mpp.relRadius, # mpp.scale,
       mpp.approxzero,
       mpp.samplingdist, mpp.samplingr, mpp.p, mpp.beta,
       PACKAGE="RandomFields", DUP=FALSE)
  ## in der naehcsten Zeile bringt valgrind einen Fehler
  ## im system, nicht klar wieso

##17264## Conditional jump or move depends on uninitialised value(s)
##17264##    at 0x401620C: (within /lib/ld-2.8.so)
##17264##    by 0x42895C3: (within /lib/libc-2.8.so)
##17264##    by 0x4289989: _dl_sym (in /lib/libc-2.8.so)
##17264##    by 0x4181DE7: (within /lib/libdl-2.8.so)
##17264##    by 0x400DE25: (within /lib/ld-2.8.so)
##17264##    by 0x41820DB: (within /lib/libdl-2.8.so)
##17264##    by 0x4181D72: dlsym (in /lib/libdl-2.8.so)
##17264##    by 0x8119155: R_dlsym (Rdynload.c:796)
##17264##    by 0x8119AC8: R_FindSymbol (Rdynload.c:845)
##17264##    by 0x8176E9C: resolveNativeRoutine (dotcode.c:238)
##17264##    by 0x817782A: do_dotCode (dotcode.c:1651)
##17264##    by 0x819FC0A: Rf_eval (eval.c:495)
##17264##    by 0x81A233A: do_begin (eval.c:1266)
##17264##    by 0x819F9D0: Rf_eval (eval.c:469)
##17264##    by 0x81A296A: do_for (eval.c:1154)
##17264##    by 0x819F9D0: Rf_eval (eval.c:469)
##17264##    by 0x81A233A: do_begin (eval.c:1266)
##17264##    by 0x819F9D0: Rf_eval (eval.c:469)
##17264##    by 0x81A36DB: Rf_applyClosure (eval.c:704)
##17264##    by 0x819F923: Rf_eval (eval.c:513)

    .C("SetParamHyperplane", m, hyper.superpos, hyper.maxlines,
       hyper.mar.distr, hyper.mar.param, NAOK=TRUE, DUP=FALSE)
    .C("SetExtremes", m, maxstable.maxGauss, PACKAGE="RandomFields", DUP=FALSE)
    if (length(parameters)==0)
      return(c(list(
                    PracticalRange=if (PracticalRange<=1)
                    as.logical(PracticalRange) else PracticalRange, 
                    PrintLevel=PrintLevel,
                    pch=pch,
                    Storing=as.logical(Storing),
                    skipchecks = if (skipchecks<=1)
                                 as.logical(skipchecks) else skipchecks,
                    stationary.only = (if(stationary.only==-1) NA else
                                       as.logical(stationary.only)),
                    exactness =if(exactness==-1) NA else as.logical(exactness),
                    every = every,

                    CE.force=as.logical(CE.force),
                    CE.mmin=CE.mmin,
                    CE.strategy=CE.strategy,
                    CE.maxmem=CE.maxmem,
                    CE.tolIm=CE.tolIm,
                    CE.tolRe=CE.tolRe,
                    CE.trials=CE.trials,
                    CE.useprimes=as.logical(CE.useprimes),
                    CE.dependent=as.logical(CE.dependent),
                    CE.method=CE.method,
                    
                    direct.bestvariables=direct.bestvariables,
                    direct.maxvariables=direct.maxvariables,
                    direct.method=direct.method,
                    direct.svdtolerance=direct.svdtolerance,
                    
                    hyper.superpos=hyper.superpos,
                    hyper.maxlines=hyper.maxlines,
                    hyper.mar.distr=hyper.mar.distr,
                    hyper.mar.param=hyper.mar.param,

                    
                    local.force=as.logical(local.force),
                    local.mmin=local.mmin,
                    local.maxmem=local.maxmem,
                    local.tolIm=local.tolIm,
                    local.tolRe=local.tolRe,
                    local.useprimes=as.logical(local.useprimes),
                    local.dependent=as.logical(local.dependent),

                    markov.neighbours = markov.neighbours,
                    markov.precision = markov.precision,
                    markov.cyclic = as.logical(markov.cyclic),
                    markov.maxmem = markov.maxmem,
                    
                    maxstable.maxGauss=maxstable.maxGauss,
                    
                    mpp.locations = mpp.locations,
                    mpp.intensity =mpp.intensity,
                    mpp.plus = mpp.plus,
                    mpp.relRadius = mpp.relRadius,
                   # mpp.scale =mpp.scale,
                    mpp.approxzero=mpp.approxzero,
                    mpp.samplingdist = mpp.samplingdist,
                    mpp.samplingr = mpp.samplingr,
                    mpp.p = mpp.p,
                    mpp.beta = mpp.beta,

                    nugget.tol=nugget.tol,
                    nugget.meth= as.logical(nugget.meth),
                    
                    sequ.max = sequ.max,
                    sequ.back = sequ.back,
                    sequ.initial = sequ.initial,
                    
                    spectral.lines=spectral.lines,
                    spectral.grid=as.logical(spectral.grid),
                    spectral.ergodic=as.logical(spectral.ergodic),
                    spectral.metro = as.logical(spectral.metro),
                    spectral.nmetro = spectral.nmetro,
                    spectral.sigma = spectral.sigma,

                    
                    TBM.method=.methods[TBM.method+1],
                    TBM.center=TBM.center,
                    TBM.points=TBM.points,
                    
                    TBM2.lines=TBM2.lines,
                    TBM2.linesimufactor=TBM2.linesimufactor,
                    TBM2.linesimustep=TBM2.linesimustep,
                    TBM2.layers=if (TBM2.layers==0 || TBM2.layers==1)
                                   as.logical(TBM2.layers) else TBM2.layers,
                    #                    TBM2.num=as.logical(TBM2.num),
                    
                    TBM3.lines=TBM3.lines,
                    TBM3.linesimufactor=TBM3.linesimufactor,
                    TBM3.linesimustep=TBM3.linesimustep,
                    TBM3.layers=if (TBM3.layers==0 || TBM3.layers==1)
                                  as.logical(TBM3.layers) else TBM3.layers,
                    
                    TBMCE.force=as.logical(TBMCE.force),
                    TBMCE.mmin=TBMCE.mmin,
                    TBMCE.strategy=TBMCE.strategy,
                    TBMCE.maxmem = TBMCE.maxmem,
                    TBMCE.tolIm=TBMCE.tolIm,
                    TBMCE.tolRe=TBMCE.tolRe,
                    TBMCE.trials=TBMCE.trials,
                    TBMCE.useprimes=as.logical(TBMCE.useprimes),
                    TBMCE.dependent=as.logical(TBMCE.dependent)
                    ),
               if (!no.readonly)
               list(
                    covmaxchar=.p$covmaxchar,
                    covnr=.p$covnr,
                    distrmaxchar=.p$distrmaxchar,
                    distrnr=.p$distrnr,
                    maxdim=.p$maxdim,
                    maxmodels=.p$maxmodels,
                    methodmaxchar=.p$methodmaxchar,
                    methodnr=.p$methodnr
                    )
               )
             )
    if (m==1) return(invisible(parameters))
   
    orig.name <- names(parameters)
    if (is.null(orig.name) || (orig.name[1]=="")) {
      txt <- "either a single unnamed list must be given or the parameters should be referenced by names"
      if ((length(parameters)!=1)) stop(txt)
      parameters <- parameters[[1]]
      orig.name <- names(parameters)
      if ((length(parameters) != sum(orig.name!="") ||
           (length(parameters)==0))) stop(txt)
    }

    name <- arg.list[pmatch(orig.name, arg.list)]
    if (any(is.na(name)))
      stop("the following parameter(s) could not be matched: ",
           paste(orig.name[is.na(name)], collapse=", "))
    names(parameters) <- name

    for (i in 1:length(parameters)) {
      type <- storage.mode(get(name[i]))
      ## parameters[i] is not sufficient since user give expression,
      ## which have type "language"
      v <- parameters[[i]]
      if (name[i]=="TBM.method" && !is.integer(v)) v <- pmatch(v, .methods) - 1
      if (name[i]=="stationary.only" && is.na(v)) v <- -1
      if (name[i]=="exactness" && is.na(v)) v <- -1
      if (switch(type,
                 character = !is.character(v),
                 integer = !is.finite(v) || (v != as.integer(v)),
                 double = !is.numeric(v)))
        stop("`", orig.name[i], "' is not ", type)
      len <- length(get(name[i]))
      if (length(v) > len)
        stop("`", orig.name[i], "' is a too long vector")
      assign(name[i], rep(v, length=len))
      eval(parse(text=paste("storage.mode(",name[i],") <- type")))
    }
    stopifnot(PracticalRange %in% c(FALSE, TRUE, 2, 3, 11, 12, 13, 999))
  }
}


PrintModelList <-function (operators=FALSE, internal=FALSE) {
   stopifnot(internal <=2, internal >=0, operators<=1, operators>=0)
    .C("PrintModelList", as.integer(internal), as.integer(operators),
       PACKAGE="RandomFields")
    invisible()
}


PrintMethodList <-function () {
    .C("PrintMethods", PACKAGE="RandomFields")
    invisible()
}



parampositions <- function(model, param, trend=NULL, dim, print=1) {
  stopifnot(!missing(dim))
  if (!is.null(trend)) stop("trend not programmed yet")
  model <- PrepareModel(model, param, trend=trend, nugget.remove=FALSE)
  .Call("GetNAPositions", model$model, as.integer(dim), as.integer(dim),
        FALSE, as.integer(print), PACKAGE="RandomFields")
}


.distInt <- function(x) {
  ##
  ## only for gene data where each coordinates takes
  ## only three neighboured integer values !
  stopifnot(is.matrix(x), is.integer(x))
  n <- nrow(x)
  genes <- ncol(x)
  res <- double(n * n)
  .C("distInt", t(x), n, genes, res, PACKAGE="RandomFields", DUP=FALSE)
  dim(res) <- rep(n, 2)
  res
}




parameter.range <- function(model, param, dim=1){
  cat("sorry not programmed yet\n")
  return(NULL)
  
  pm <- PrepareModel(model=model, param=param, nugget.remove=FALSE)        
  storage.mode(dim) <- "integer"
  ResGet <- .Call("MLEGetModelInfo", pm$model, dim, dim,
                  PACKAGE="RandomFields")
  minmax <- ResGet$minmax[, 1:2]
  dimnames(minmax) <- list(attr(ResGet$minmax, "row.names"), c("min", "max"))
  return(minmax)
}
