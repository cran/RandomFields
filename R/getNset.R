
paramextract <- function(p, model=c("cutoff")) {
  i <- .C("GetParamterPos", variance=integer(1), kappa=integer(1),
          lastkappa=integer(1), tbm2num=integer(1), hyperinternal=integer(1),
          lasthyperinternal=integer(1), scale=integer(1), aniso=integer(1),
          hypernr=integer(1), localdiameter=integer(1), localr=integer(1),
          cutoffr=integer(1),
          hyperkappa=integer(1), total=integer(1), PACKAGE="RandomFields")
  i <- lapply(i, function(x) x + 1)
  model <- match.arg(model)
  return(switch(model,
                cutoff = list(hypernr=p[i$hypernr], diameter=p[i$localdiameter],
                  lcoalr=p[i$localr], cutoffr=p[i$cutoffr])
                )
         )
}

CheckAndComplete <- function(model, param, dim) {
    p <- PrepareModel(model, param, dim)
    i <- .C("GetParamterPos", variance=integer(1), kappa=integer(1),
            lastkappa=integer(1), tbm2num=integer(1), hyperinternal=integer(1),
            lasthyperinternal=integer(1), scale=integer(1), aniso=integer(1),
            hypernr=integer(1), localdiameter=integer(1), localr=integer(1),
            cutoffr=integer(1),
            hyperkappa=integer(1), total=integer(1), PACKAGE="RandomFields")
    error = integer(1)
    param = double(i$total * length(p$covnr))
    .C("CheckAndCompleteParameters",
       as.integer(p$covnr),
       as.double(p$param), 
       as.integer(length(p$param)), ## control
       as.integer(dim), 
       as.integer(length(p$covnr)),
       as.integer(p$anisotropy),
       as.integer(p$op),
       param,
       error,
       PACKAGE="RandomFields", DUP = FALSE)$res
    dim(param) = c(i$total, length(p$covnr))
    return(list(error=error, param=param))
}

GetRegisterInfo <- function(register=0, ignore.active=FALSE)
  # ignore.active=TRUE only for internal debugging information!
  .Call("GetExtKeyInfo", as.integer(register), as.logical(ignore.active),
        PACKAGE="RandomFields")

GetPracticalRange <- function(model, kappas=NULL) {
  covnr <-
    as.integer(.C("GetModelNr", as.character(model), as.integer(1),
                  nr = integer(1), PACKAGE="RandomFields")$nr)
  if (covnr < 0) {
    .C("PrintModelList", PACKAGE="RandomFields")
    stop("given model cannot be (uniquely) identified from the above list")
  }
  if (length(kappas)!=.C("GetNrParameters", covnr, as.integer(1),
              k = integer(1), PACKAGE="RandomFields", DUP=FALSE)$k)
    stop("incorrect number of parameters!")
  nat.scl <- double(1)
  error <- integer(1)
  .C("GetNaturalScaling",
     as.integer(covnr),
     as.double(kappas),         ## not stable w.r.t. to changings !!
     as.integer(11),
     nat.scl,
     error,
     PACKAGE="RandomFields", DUP=FALSE)
  if (error) stop("natural scaling could not be obtained")
  return(1.0 / nat.scl)
}

GetMethodNames <- function() {
  assign(".p",
         .C("GetrfParameters", covmaxchar=integer(1), methodmaxchar=integer(1),
            distrmaxchar=integer(1),
            covnr=integer(1), methodnr=integer(1), distrnr=integer(1),
            maxdim=integer(1), maxmodels=integer(1),
            PACKAGE="RandomFields"))
  l <- character(.p$methodnr)
  for (i in 1:.p$methodnr) {
    l[i] <- .C("GetMethodName", as.integer(i-1),
               n=paste(rep(" ", .p$methodmaxchar), collapse=""),
               PACKAGE="RandomFields")$n
  }
  return(l)
}

GetDistributionNames <- function() {
  assign(".p",
         .C("GetrfParameters", covmaxchar=integer(1), methodmaxchar=integer(1),
            distrmaxchar=integer(1),
            covnr=integer(1), methodnr=integer(1), distrnr=integer(1),
            maxdim=integer(1), maxmodels=integer(1),
            PACKAGE="RandomFields"))
  l <- character(.p$distrnr)
  for (i in 1:.p$distrnr) {
    l[i] <- .C("GetDistrName", as.integer(i-1),
               n=paste(rep(" ",.p$distrmaxchar), collapse=""),
               PACKAGE="RandomFields")$n
  }
  return(l)
}


GetModelNames <- function() {
  assign(".p",
         .C("GetrfParameters", covmaxchar=integer(1), methodmaxchar=integer(1),
            distrmaxchar=integer(1),
            covnr=integer(1), methodnr=integer(1), distrnr=integer(1),
            maxdim=integer(1), maxmodels=integer(1),
            PACKAGE="RandomFields"))
  l <- character(.p$covnr)
  for (i in 1:.p$covnr) {
    l[i] <- .C("GetModelName",as.integer(i-1),
               n=paste(rep(" ",.p$covmaxchar), collapse=""),
               PACKAGE="RandomFields")$n
  }
  return(l)
}


GetModelList <- function(abbr=TRUE) {
  names <- GetModelNames()
  methods <- GetMethodNames()
  if (abbr) methods <- substr(methods, 1, if (is.logical(abbr)) 5 else abbr)
  idx <- integer(length(names) * length(methods))
  .C("GetModelList", idx, PACKAGE="RandomFields", DUP=FALSE)
  t(matrix(as.logical(idx), ncol=length(names), dimnames=list(methods, names)))
}

parampositions <- function(model, param, print=TRUE) {
  type <- if (!missing(param) && !is.null(param))
    if (is.matrix(param)) "n" else "s" else "l"
  old.model <- PrepareModel(model, param)
  model <- PrepareModel(convert.to.readable(old.model))
  if (length(old.model$param) != length(model$param))
    stop("The model should be simplified beforehand") 
  model$param <- 1:length(model$param)
  model$mean <- NA
  model <- convert.to.readable(model, allowed=type)
  model$method <-  model$trend <- NULL
  if (type=="l") {
    if (print) str(model)
  } else {
    if (print) cat("model:", model$model, "\nparam: ")
    if (type=="s") { # standard
      NUGGET <- 3
      if (is.finite(param[NUGGET]) && param[NUGGET]==0)
        model$param[NUGGET] <- NA
      if (print) cat(model$param, "\n")
    } else { # nested
      model$param[model$param==0] <- NA
      if (length(model$param) !=
          length(convert.to.readable(old.model, allow="n")$param))
        stop("Model is too complex to be identified")
      if (print) {
        cat("\n")
        print(model$param)
      }
    }
  }
  invisible(model)
}

"RFparameters" <- function (..., no.readonly=FALSE) {
  assign(".p",
         .C("GetrfParameters", covmaxchar=integer(1), methodmaxchar=integer(1),
            distrmaxchar=integer(1),
            covnr=integer(1), methodnr=integer(1), distrnr=integer(1),
            maxdim=integer(1), maxmodels=integer(1),
            PACKAGE="RandomFields"))
  assign(".methods", GetMethodNames())
  ## do not add any temporary variable til ## **
  ## do not remove leading "." from .maxdim
  Storing <- integer(1)
  PrintLevel <- integer(1)
  PracticalRange <- integer(1)
  ## always logical returned
  ## PracticalRange also allows for being set to
  ##  0 : no practical range
  ##  1,11 : practical range, evaluated exactly (if given in RFCovfct.cc)
  ##  2,12 : approximative value (if given in RFCovfct.cc)
  ##  3,13 : rough guess (good enough for MLE) (if given in RFCovfct.cc)
  ## >10: and if nothing appropriate given in RFCovfct.cc then numerical approx.
  pch <- as.character("  ")
  stationary.only <- integer(1)
  exactness <- integer(1)
  
  CE.force <- integer(1)
  CE.tolRe <- double(1)
  CE.tolIm <- double(1)
  CE.trials <- integer(1)
  CE.strategy <- integer(1)
  CE.mmin <- double(.p$maxdim)
  CE.enlarge <- integer(1)
  CE.maxmem <- double(1)
  CE.useprimes <- integer(1)
  
  local.force <- integer(1)
  local.tolRe <- double(1)
  local.tolIm <- double(1)
  local.trials <- integer(1)
  local.strategy <- integer(1)
  local.mmin <- double(.p$maxdim)
  local.enlarge <- integer(1)
  local.maxmem <- double(1)
  local.useprimes <- integer(1)
  
  TBMCE.force <- integer(1)
  TBMCE.tolRe <- double(1)
  TBMCE.tolIm <- double(1)
  TBMCE.trials <- integer(1)
  TBMCE.strategy <- integer(1)
  TBMCE.mmin <- double(.p$maxdim)
  TBMCE.enlarge <- integer(1)
  TBMCE.maxmem <- double(1)
  TBMCE.useprimes <- integer(1)
  
  TBM.method <- integer(1)
  TBM.center <- double(.p$maxdim)
  TBM.points <- integer(1)

  TBM2.lines <- integer(1)
  TBM2.linesimufactor <- double(1)
  TBM2.linesimustep <- double(1)
  TBM2.every <- integer(1)
  TBM2.num <- integer(1)
  
  
  TBM3.lines <- integer(1)
  TBM3.linesimufactor <- double(1)
  TBM3.linesimustep <- double(1)
  TBM3.every <- integer(1)

  spectral.lines <- integer(1)
  spectral.grid <- integer(1)

  direct.method <- integer(1)
  direct.checkprecision <- integer(1)
  direct.requiredprecision <- double(1)
  direct.bestvariables <- integer(1)
  direct.maxvariables <- integer(1)

  MPP.approxzero <- double(1)
  add.MPP.realisations <- double(1)
  MPP.radius <- double(1)

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
    x <- .C("SetParam", m, Storing=Storing, PrintLevel=PrintLevel,
            PracticalRange=PracticalRange, pch=pch, PACKAGE="RandomFields")
    Storing <- x$Storing
    PrintLevel <- x$PrintLevel
    PracticalRange <- x$PracticalRange
    pch <- x$pch
    .C("SetParamDecision", m, stationary.only, exactness,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamCircEmbed", m, CE.force, CE.tolRe, CE.tolIm, CE.trials, 
       CE.mmin, CE.useprimes, CE.strategy, CE.maxmem,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamLocal", m, local.force, local.tolRe, local.tolIm, local.trials,
       local.mmin, local.useprimes, local.strategy, local.maxmem,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamTBMCE", m, TBMCE.force, TBMCE.tolRe, TBMCE.tolIm, TBMCE.trials,
       TBMCE.mmin, TBMCE.useprimes, TBMCE.strategy, TBMCE.maxmem,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamTBM2", m, TBM2.lines, TBM2.linesimufactor,
       TBM2.linesimustep, TBM2.every, TBM2.num,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamTBM3", m, TBM3.lines, TBM3.linesimufactor,
       TBM3.linesimustep, TBM3.every, PACKAGE="RandomFields", DUP=FALSE)
   .C("SetParamTBM", m, TBM.method, TBM.center, TBM.points,
       PACKAGE="RandomFields", DUP=FALSE, NAOK=TRUE)
    .C("SetParamSpectral", m, spectral.lines, spectral.grid,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamDirectGauss", m, direct.method, direct.checkprecision,
       direct.requiredprecision, direct.bestvariables, direct.maxvariables,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetMPP", m, MPP.approxzero, add.MPP.realisations, MPP.radius,
       PACKAGE="RandomFields", DUP=FALSE)
    .C("SetParamHyperplane", m, hyper.superpos, hyper.maxlines,
       hyper.mar.distr, hyper.mar.param, NA.OK=TRUE, DUP=FALSE)
    .C("SetExtremes", m, maxstable.maxGauss, PACKAGE="RandomFields", DUP=FALSE)
    
    if (length(parameters)==0)
      return(c(list(
                    PracticalRange=if (PracticalRange<=1)
                    as.logical(PracticalRange) else PracticalRange, 
                    PrintLevel=PrintLevel,
                    pch=pch,
                    Storing=as.logical(Storing),
                    stationary.only = (if(stationary.only==-1) NA else
                                       as.logical(stationary.only)),
                    exactness = if(exactness==-1) NA else as.logical(exactness),
                    CE.force=as.logical(CE.force),
                    CE.mmin=CE.mmin,
                    CE.strategy=CE.strategy,
                    CE.maxmem=CE.maxmem,
                    CE.tolIm=CE.tolIm,
                    CE.tolRe=CE.tolRe,
                    CE.trials=CE.trials,
                    CE.useprimes=as.logical(CE.useprimes),
                    local.force=as.logical(local.force),
                    local.mmin=local.mmin,
                    local.strategy=local.strategy,
                    local.maxmem=local.maxmem,
                    local.tolIm=local.tolIm,
                    local.tolRe=local.tolRe,
                    local.trials=local.trials,
                    local.useprimes=as.logical(local.useprimes),
                    direct.checkprecision=as.logical(direct.checkprecision),
                    direct.bestvariables=direct.bestvariables,
                    direct.maxvariables=direct.maxvariables,
                    direct.method=direct.method,
                    direct.requiredprecision=direct.requiredprecision,
                    spectral.grid=as.logical(spectral.grid),
                    spectral.lines=spectral.lines,
                    TBM.method=.methods[TBM.method+1],
                    TBM.center=TBM.center,
                    TBM.points=TBM.points,
                    TBM2.every=TBM2.every,
                    TBM2.lines=TBM2.lines,
                    TBM2.linesimufactor=TBM2.linesimufactor,
                    TBM2.linesimustep=TBM2.linesimustep,
                    TBM2.num=as.logical(TBM2.num),
                    TBM3.every=TBM3.every,
                    TBM3.lines=TBM3.lines,
                    TBM3.linesimufactor=TBM3.linesimufactor,
                    TBM3.linesimustep=TBM3.linesimustep,
                    TBMCE.force=as.logical(TBMCE.force),
                    TBMCE.mmin=TBMCE.mmin,
                    TBMCE.strategy=TBMCE.strategy,
                    TBMCE.maxmem = TBMCE.maxmem,
                    TBMCE.tolIm=TBMCE.tolIm,
                    TBMCE.tolRe=TBMCE.tolRe,
                    TBMCE.trials=TBMCE.trials,
                    TBMCE.useprimes=as.logical(TBMCE.useprimes),
                    add.MPP.realisations=add.MPP.realisations,
                    MPP.approxzero=MPP.approxzero,
                    MPP.radius=MPP.radius,
                    hyper.superpos=hyper.superpos,
                    hyper.maxlines=hyper.maxlines,
                    hyper.mar.distr=hyper.mar.distr,
                    hyper.mar.param=hyper.mar.param,
                    maxstable.maxGauss=maxstable.maxGauss,
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
                    methodnr=.p$methodnr,
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
      if (name[i]=="TBM.method") v <- pmatch(v, .methods) - 1
      if (name[i]=="stationary.only" && is.na(v)) v <- -1
      if (name[i]=="exactness" && is.na(v)) v <- -1
      if (switch(type,
                 character = !is.character(v),
                 integer = !is.finite(v) || (v != as.integer(v)),
                 double = !is.numeric(v)))
        stop(paste("`", orig.name[i], "' is not ", type, sep=""))
      len <- length(get(name[i]))
      if (length(v) > len)
        stop(paste("`", orig.name[i], "' is a too long vector", sep=""))
      assign(name[i], rep(v, length=len))
      eval(parse(text=paste("storage.mode(",name[i],") <- type")))
    }
    stopifnot(PracticalRange %in% c(FALSE, TRUE, 2, 3, 11, 12, 13))
  }
}


"PrintModelList" <-function () {
    .C("PrintModelList", PACKAGE="RandomFields")
    invisible()
}


"PrintMethodList" <-function () {
    .C("PrintMethods", PACKAGE="RandomFields")
    invisible()
}


parameter.range <- function(model, dim){
  if (length(model)==0) stop("model not given")
  stopifnot(is.character(model))
  nr <- .C("GetModelNr", as.character(model), as.integer(1), nr=integer(1),
           PACKAGE="RandomFields")$nr
  if (nr < 0) {
    .C("PrintModelList", PACKAGE="RandomFields")
    stop("given model cannot be (uniquely) identified from the above list")
  }
  storage.mode(nr) = "integer"
  storage.mode(dim) = "integer"
  l <- as.integer(4 * .C("GetNrParameters", nr, as.integer(1), k=integer(1),
                         PACKAGE="RandomFields")$k)
  r <- list()
  r$theoretical <- list()
  r$practical <- list()
  index <- as.integer(1)
  while (index>0) {
    R <- double(l)
    index.orig <- as.integer(index) ## without index.orig points to index,
    ## what is a bug in R -- do report! -- Check if bug is still there
    .C("GetRange", nr, dim, index, R, l, PACKAGE="RandomFields", DUP=FALSE)
    R <- matrix(R, nrow=4)
    r$theoretical[[index.orig]] <- R[1:2, , drop=FALSE]
    r$practical[[index.orig]] <- R[3:4, ,drop=FALSE]
  }
  if (index <= -2) {
    if (index==-2) r <- NaN ##  stop("dimension not correct")
    else stop(paste("error: inform maintainer (error nr.",index,")"))
  }
  return(if (is.list(r) && ncol(r$theoretical[[1]])==0) NULL else r)
}
