
summary.RFopt <- function(object, ...) {
  object <- lapply(object, function(z) z[order(names(z))])
  object <- object[c(1, 1 + order(names(object[-1])))]
  class(object) <- "summary.RFopt"
  object
}


print.summary.RFopt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFopt <- function(x, ...) {
  print.summary.RFopt(summary.RFopt(x, ...)) #
  invisible(x)
}

summary.RFoptElmnt <- function(object, ...) {
  object <- object[order(names(object))]
  class(object) <- "summary.RFoptElmt"
  object
}

print.summary.RFoptElmnt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFoptElmnt <- function(x, ...) {
  print.summary.RFoptElmnt(summary.RFoptElmnt(x, ...)) #
  invisible(x)
}


RFoptions <- function(..., no.readonly=TRUE) {
##  on.exit(.C("RelaxUnknownRFoption", FALSE))
##  .C("RelaxUnknownRFoption", TRUE)
  opt <- lapply(.External("RFoptions", ...),
                function(x) {
                  class(x) <- "RFoptElmnt"
                  x
                })
  if (length(opt)!=0) {      
    class(opt) <-  "RFopt"
    if (!no.readonly) {
      opt$readonly <- list(covmaxchar=MAXCHAR,
                           covnr=GetCurrentNrOfModels(FALSE),
                           maxdim=c(cov=MAXCOVDIM, mle=MAXMLEDIM,
                               simu=MAXSIMUDIM, ce=MAXCEDIM, tbm=MAXTBMSPDIM,
                               mpp=MAXMPPDIM, hyper=MAXHYPERDIM,
                               nug=MAXNUGGETDIM, vario=MAXVARIODIM),
                           maxmodels=MAXFIELDS,
                           methodmaxchar=METHODMAXCHAR,
                           methodnr=Forbidden
                           )
    }
  }
  if (length(opt)==0) invisible(opt) else opt
}

GetCurrentNrOfModels <- function(init=TRUE) {
  res<- .C("GetCurrentNrOfModels", as.integer(init), nr=as.integer(1))$nr
}

internal.rfoptions <- function(..., REGISTER=FALSE, COVREGISTER=as.integer(NA),
                               RELAX=FALSE){
  RFopt <- list()
  RFopt[[1]] <- .External("RFoptions")
  if (is.logical(REGISTER)) {
    REGISTER <- if (REGISTER) RFopt[[1]]$registers$register else as.integer(NA)
  }
  RFopt[[1]]$general$storing <-
    c(RFopt[[1]]$general$storing, REGISTER, COVREGISTER)
  l <- list(asList=FALSE, ...)
#  l <- list(...)
  if (length(l) > 1) { ## wegen asList; sonst  > 0 !!
    storing <- (substr(names(l), 1, 3) == "sto" |
                substr(names(l), 1, 9) == "general.sto")
    if (any(storing)) last <- rev(which(storing))[1]
    if (any(storing) && !l[[last]]) {
      for (p in which(storing)) l[[p]] <- c(FALSE, REGISTER, COVREGISTER)
    }
    on.exit(.C("RelaxUnknownRFoption", FALSE))
    .C("RelaxUnknownRFoption", RELAX)
    .External("RFoptions", LIST=l)
    RFopt[[2]] <- .External("RFoptions")    
  } else {
    RFopt[[2]] <- RFopt[[1]]
  }
  return(RFopt)
}


print_RFgetModelInfo <- function(x, max.level=10-attr(x, "level"),
                                 give.attr=FALSE, ...) {
  str(object = x,  max.level=max.level, give.attr=give.attr, ...) #
}

print.RFgetModelInfo <- function(x, ...) {
  print_RFgetModelInfo(x, ...) 
}


RFgetModelInfo <-
  function(register, level=1, 
           spConform=RFoptions()$general$spConform,
           which.submodels = c("user", "internal", "both"),
           modelname=NULL) {  
  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register)) register else
                               deparse(substitute(register)))
  which.submodels <- match.arg(which.submodels)
  ## positive values refer the covariance models in the registers
  ##define MODEL_USER : 0  /* for user call of Covariance etc. */
  ##define MODEL_SIMU : 1  /* for GaussRF etc */ 
  ##define MODEL_INTERN  : 2 /* for kriging, etc; internal call of covariance */
  ##define MODEL_MLE  : 3
  ##define MODEL_BOUNDS  4 : - /* MLE, lower, upper */
  ## level + 10: auch die call fctn !
  ## [ in RF.h we have the modulus minus 1 ]
  
  cov <- .Call("GetExtModelInfo", as.integer(register), as.integer(level),
               as.integer(spConform),
               as.integer(if (which.submodels=="user") 0 else
                          1 + (which.submodels=="both")),
               PACKAGE="RandomFields")

  if (!is.null(modelname)) {
    cov <- search.model.name(cov, modelname, 0)
  }
  class(cov) <- "RFgetModelInfo"
  attr(cov, "level") <- level
 
  return(cov)
}



RFgetModelNames <- function(type = RC_TYPENAMES, domain = RC_DOMAIN_NAMES,
                            isotropy = RC_ISONAMES, operator = c(TRUE, FALSE),
                            monotone = RC_MONOTONE_NAMES,
                            implied_monotonicities = length(monotone) == 1,
                            finiterange = c(TRUE, FALSE, NA),
                            valid.in.dim = c(1, Inf), 
                            vdim = c(1, 5),
                            group.by,
                            simpleArguments = FALSE,
                            internal,
                            newnames
                            ){ #, .internal=FALSE) {

  group.names <- c("type", "domain", "isotropy", "operator",
                   "monotone", "finiterange", "valid.in.dim", "vdim")

  if (hasArg(internal)) {
    return(PrintModelList(operators=operator, internal = internal,
                          newstyle=missing(newnames) || newnames))
  }
  if (!missing(newnames) && !newnames) {
    if (hasArg(type) || hasArg(domain) || hasArg(isotropy) || hasArg(operator)
        || hasArg(monotone) || hasArg(finiterange) || hasArg(valid.in.dim)
        || hasArg(vdim) || hasArg(group.by))
      stop("use 'newnames=FALSE' without further parameters or in combination with 'internal'")
    return (.Call("GetAllModelNames", PACKAGE="RandomFields"))
  }
  

  if (!(length(valid.in.dim) %in% 1:2)) stop("'valid.in.dim' has wrong size.")
  if (length(valid.in.dim) == 1) valid.in.dim <- c(valid.in.dim, Inf)
  
  if (!(length(vdim) %in% 1:2)) stop("'vdim' has wrong size.")
  if (length(vdim) == 1) vdim <- rep(vdim, 2)

  debug <- !TRUE
  if (hasArg(type)) type <- TYPENAMES[pmatch(type, TYPENAMES)]
  if (hasArg(domain)) domain <- DOMAIN_NAMES[pmatch(domain, DOMAIN_NAMES)]
  if (hasArg(isotropy)) isotropy <- ISONAMES[pmatch(isotropy, ISONAMES)]
  if (hasArg(monotone)) monotone <- MONOTONE_NAMES[pmatch(monotone,
                                                           MONOTONE_NAMES)]
  if (any(is.na(pmatch(type, TYPENAMES))))
    stop(paste("'", type, "'", " is not a valid category", sep=""))
  if (any(is.na(pmatch(domain, DOMAIN_NAMES))))
    stop(paste("'", domain, "'", " is not a valid category", sep=""))
  if (any(is.na(pmatch(isotropy, ISONAMES))))
    stop(paste("'", isotropy, "'", " is not a valid category", sep=""))
  if (any(is.na(pmatch(monotone, MONOTONE_NAMES))))
    stop(paste("'", monotone, "'", " is not a valid category", sep=""))


  if (missing(group.by)) {
    if (length(type) == 1) {
      if (type == TYPENAMES[TcfType + 1]) type <- c(type, "undefined") # to do
      else if (type == TYPENAMES[PosDefType + 1])
        type <- c(TYPENAMES[TcfType + 1], type, "undefined") # to do
      else if (type == TYPENAMES[NegDefType + 1])
        type <- c(TYPENAMES[c(TcfType, PosDefType) + 1], type,"undefined")#to do
     }
     if (!hasArg("group.by")) group.by <- if (length(type) == 1) NULL else 'type'
  }
 
  if (length(group.by) > 0) {
    group.idx <- pmatch(group.by, group.names)
    if (any(is.na(group.idx)))
      stop("'group.by' can be equal to '",
           paste(group.names, collapse="', '"), "'")
    group.by <- group.names[group.idx]
    group.idx <- group.idx[1]
  }

  
  if (group <- !is.null(group.by)) {  
    FUN <- function(string){
      args <- list(type=type, domain=domain, isotropy=isotropy,
                   operator=operator, monotone=monotone,
                   implied_monotonicities = implied_monotonicities &&
                   group.by[1] != "monotone",
                   finiterange=finiterange, valid.in.dim=valid.in.dim,
                   vdim=vdim,
                   group.by =
                     if (group && length(group.by) > 1) group.by[-1] else NULL
                   )
      args[[group.idx]] <- string
      list(do.call("RFgetModelNames", args))
    }
    li <- sapply(get(group.by[1]), FUN=FUN)
    if (is.null(names(li)))
      names(li) <- paste(group.by[1], get(group.by[1]), sep="=")
    length <- unlist(lapply(li, FUN=length))
    li <- li[length>0]
    return(li)
  } # matches  if (hasArg(group.by)) {
  
  if (implied_monotonicities) {
    mon <- MONOTONE_NAMES[-1 : (MISMATCH - 1)]
    for (i in MONOTONE:NORMAL_MIXTURE)
      if (mon[i] %in% monotone) monotone <- c(monotone, mon[i + 1])
    if (mon[BERNSTEIN] %in% monotone) monotone <- c(monotone, mon[MONOTONE])
    monotone <- unique(monotone)
  }

 
  envir <- as.environment("package:RandomFields")
  ls.RFmg <- ls(envir=envir)
  ls.RFmg <- ls.RFmg[substr(ls.RFmg, 1, 1) == "R"]
  idx <- logical(len <- length(ls.RFmg))
   for (i in 1:len){
    fun <- get(ls.RFmg[i], envir=envir)
    idx[i] <- is.function(fun) && is(fun, class2="RMmodelgenerator")
    if (!idx[i]) next
     idx[i] <- (!all(is.na(pmatch(fun["type"], type, 
                                        duplicates.ok=TRUE) &
                           pmatch(fun["isotropy"], isotropy, 
                                  duplicates.ok=TRUE))) &&
          #      !any(pt[!is.na(pt)] %in% pi[!is.na(pi)]) &&
                !all(is.na(pmatch(domain, fun["domain"]))) &&
                fun["operator"] %in% operator &&
                !all(is.na(pmatch(monotone, fun["monotone"]))) &&
                (!simpleArguments || fun["simpleArguments"]) &&
                fun["finiterange"] %in% finiterange &&
                (fun["maxdim"] < 0 ||
                 (fun["maxdim"] >= valid.in.dim[1] &&
                  fun["maxdim"] <= valid.in.dim[2]))  &&
                (fun["vdim"] < 0 ||
                 (fun["vdim"] >= vdim[1] && fun["vdim"] <= vdim[2]))  &&
                ls.RFmg[i] != ZF_INTERNALMIXED                
                )
   }

  
  return(unique(sort(ls.RFmg[idx])))
}


RFformula <- function(f)
  return(parseModel(f))


RFgetMethodNames <-function () {
  RFgetModelNames(type=TYPENAMES[c(GaussMethodType, BrMethodType) + 1])
}


GetProcessType <- function(model) {
  stopifnot(is.list(model))
  return(.Call("GetProcessType", MODEL_INTERN, model))
}
