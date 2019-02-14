
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017  Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



### get&set functions using the C interface
### and RFgetNset functions

GetCurrentNrOfModels <- function() {
  .C(C_GetCurrentNrOfModels, nr=as.integer(1))$nr
}

xylabs <- function(x, y, T=NULL, units=NULL) {
  if (is.null(units)) units <- RFoptions()$coords$coordunits
  xlab <- if (is.null(x)) NULL
          else if (units[1]=="") x else paste(x, " [", units[1], "]", sep="")
  ylab <- if (length(y)==0) NULL
          else if (units[2]=="") y else paste(y, " [", units[2], "]", sep="")
  Tlab <- if (length(T)==0) NULL
          else if (units[3]=="") T else paste(T, " [", units[3], "]", sep="")
  return (list(xlab=xlab, ylab=ylab, Tlab=Tlab))
}

add.units <- function(x,  units=NULL) {
    if (is.null(x)) return(NULL)
  if (is.null(units)) units <- RFoptions()$coords$varunits
  return(ifelse(units=="", x, paste(x, " [", units, "]", sep="")))
}


internal.rfoptions <- function(..., REGISTER=FALSE, COVREGISTER=as.integer(NA),
                               RELAX=FALSE){
### TO DO: FLAG GLOBAL (ja/nein) an C_* funktion uebergeben
### an RandomFieldsUtils::RFoptions ; am besten eine PseudoOption aehnlich
### zu "SAVEOPTIONS" und "GETOPTIONS"
  
 # Print(list(...))
  RFopt <- list()
  RFopt[[1]] <- RandomFieldsUtils::RFoptions(#LOCAL=TRUE
                                             )
  if (is.logical(REGISTER)) {
    REGISTER <- if (REGISTER) RFopt[[1]]$registers$register else as.integer(NA)
  }
  RFopt[[1]]$general$storing <-
    c(RFopt[[1]]$general$storing, REGISTER, COVREGISTER)
  l <- list(#LOCAL=TRUE,
            asList=FALSE, ...)
#  l <- list(...)
  if (length(l) > 1) { ## wegen asList; sonst  > 0 !!
    storing <- (substr(names(l), 1, 3) == "sto" |
                substr(names(l), 1, 9) == "general.sto")
    if (any(storing)) last <- rev(which(storing))[1]
    if (any(storing) && !l[[last]]) {
      for (p in which(storing)) l[[p]] <- c(FALSE, REGISTER, COVREGISTER)
    }
    on.exit(.C(C_RelaxUnknownRFoptions, as.integer(FALSE)))
    .C(C_RelaxUnknownRFoptions, as.integer(RELAX))
    RandomFieldsUtils::RFoptions(LIST=l)
    RFopt[[2]] <- RandomFieldsUtils::RFoptions(#LOCAL=TRUE
                                               )
  } else {
    RFopt[[2]] <- RFopt[[1]]
  }
  return(RFopt)
}


InitModel <- function(reg, model, dim, NAOK=FALSE){ # ok
  for (y in list(double(0), matrix(nrow=dim, ncol=3, as.double(1:3)))) {
    vdim <- try(.Call(C_Init, as.integer(reg), model,
                      list(x=matrix(nrow=dim, ncol=3, as.double(1:3)), #0 nur dummies
                           y=y, #1 y 
                           as.double(0), #2 T
                           FALSE, #3 grid
                           as.integer(dim), #4 spatdim
                           FALSE, #5 Zeit
                           FALSE #6 distances
                           ),
                      NAOK=NAOK # ok
                      ), silent=TRUE)
    if (is.numeric(vdim)) return(vdim)
    msg <- strsplit(vdim[[1]], "\n")[[1]][2]
    if (RFoptions()$basic$printlevel >= PL_ERRORS) cat(msg, "\n")
   }
  stop(msg)
  ##  stop("model could not be initialized")
}


resolve.register <- function(register){
  if (missing(register) || length(register) == 0) {
    register <- .C(C_GetCurrentRegister, reg=integer(1))$reg
    if (register < 0) stop("model that has been used right now cannot be determined or no model has been used up to now")
  }
  if (!is.numeric(register)) {
 #   register <- deparse(substitute(register))   
    register <-
      switch(register,
             "RFcov" = MODEL_COV,
             "RFcovmatrix" = MODEL_COVMATRIX,
             "RFfctn" = MODEL_FCTN,
             "RFpseudovariogram" =  MODEL_PSEUDO,
             "RFvariogram" = MODEL_VARIOGRAM,
             
             "RFdistr" = MODEL_DISTR,
             "RFddistr" = MODEL_DISTR,
             "RFpdistr" = MODEL_DISTR,
             "RFqdistr" = MODEL_DISTR,
             "RFrdistr" = MODEL_DISTR,
             
             "RFcrossvalidate" =  MODEL_MLE,
             "RFfit" =  MODEL_MLE,
             "RFgui" =  MODEL_GUI,
             "RFinterpolate" =  MODEL_KRIGE,
             "RFlikelihood" = RFoptions()$register$likelihood_register,
             "RFlinearpart" = MODEL_USER,
             "RFratiotest" =  MODEL_MLE,
             "RFsimulate" =  RFoptions()$registers$register,
             stop("register unknown")
             )
  }
  stopifnot(is.numeric(register))
  if (register < 0) stop("'register' does not have a non-negative value.")
  return(register)
}

shorten.info <- function(x) {
  sub <- which(names(x) == "submodels")
  if (length(sub) == 1) {    
    names(x[[sub]]) <- lapply(x[[sub]], function(y) y$name)
    for (i in 1:length(x[[sub]])) {
      x[[sub]][[i]] <- shorten.info(x[[sub]][[i]])
      x[[sub]][[i]]$name <- NULL ## muss immer das erste sein      
      if (length(x[[sub]][[i]]$coordinates) > 0) {    
        if (length(x[[sub]][[i]]) == 1)
          x[[sub]][[i]] <- x[[sub]][[i]]$coordinates
      }
    }

    if (length(x$submodels) == 1) {
      sub <- which(names(x) == "submodels")
      names(x)[sub] <- names(x[[sub]])
      x[[sub]] <- x[[sub]][[1]]
    }
  }
  
  return(x)
}

print_RFgetModelInfo <- function(x, max.level=99, short=TRUE,
                                 give.attr=FALSE, ...) {
  if (is.null(attr(x, "level"))) {
    y <- x
    y$minmax <- NULL
    str(y, give.attr=FALSE) #
    types <- sort(unique(x$minmax$type))
    if (length(types) > 0) {
      cat(" $ minmax: \n")
      print(x$minmax, justify="right") #
      cat(" pmin/pmax : bound usually met in practice\n",
          "type\n",
          paste("   type =", formatC(types, width=2), ":",
                TYPEOF_PARAM_NAMES[types + 1], "\n"),
          "NAN/bayes : internal\n",
          "min/max   : mathematically valid interval for the parameter\n",
          "omin/omax : whether the interval is open to the left/right\n",
          "col/row   : size of parameter\n")
    }
  } else {
    if (short) x <- shorten.info(x)
    str(object = x,  max.level=max.level, give.attr=give.attr, ...) #
  }  
}

print.RFgetModelInfo <- function(x, ...) {
  print_RFgetModelInfo(x, ...) 
}


RFgetModelInfo <- function(...) {
  x <- list(...)
  if (length(x) > 0 &&
      (is(x[[1]], CLASS_CLIST) || is(x[[1]], "list"))) RFgetModelInfo_model(...)
  else RFgetModelInfo_register(...)
}

RFgetModelInfo_model <- function(model, params, dim = 1, Time = FALSE,
                                 kernel = FALSE, exclude_trend = TRUE, ...) {
  Reg <- MODEL_AUX
  RFoptOld <- internal.rfoptions(..., RELAX=is(model, "formula"))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  if (!exclude_trend) {
    stop("'exclude_trend = FALSE' not programmed yet.")
  }

  dim <- as.integer(dim)
  intern <- try(.Call(C_SetAndGetModelInfo, Reg,
                      list("Dummy", PrepareModel2(model, params=params, ...)),
                          dim, FALSE, as.logical(kernel), as.logical(Time), dim,
                          as.integer(10), ## ehemals RFoptions(short=10)
                          TRUE, as.logical(exclude_trend) ))
  if (is(intern, "try-error")) return(0)
  intern$effect <- intern$xdimOZ <- intern$matrix.indep.of.x <- NULL

  intern$minmax <- as.data.frame(intern$minmax)
  storage.mode(intern$minmax$omin) <- "logical"
  storage.mode(intern$minmax$omax) <- "logical"
  storage.mode(intern$minmax$bayes) <- "logical"
  storage.mode(intern$minmax$NAN) <- "integer"
  storage.mode(intern$minmax$col) <- "integer"
  storage.mode(intern$minmax$row) <- "integer"
 
  class(intern) <- "RFgetModelInfo"
  intern
}

RFgetModelInfo_register <- function(register, level=1, 
   spConform=RFoptions()$general$spConform,
   which.submodels = c("user", "internal",
		       "call+user", "call+internal",
		       "user.but.once", "internal.but.once",
		       "user.but.once+jump", "internal.but.once+jump", "all"),
   modelname=NULL, origin = "original") {    
  register <- resolve.register(if (missing(register)) NULL
			       else if (is.numeric(register)) register
			       else deparse(substitute(register)))
  ##  Print(which.submodels, as.list(RFgetModelInfo)$which.submodels[-1])  
  w <- (if (!hasArg("which.submodels")) 0
        else pmatch(match.arg(which.submodels),
                    as.list(RFgetModelInfo_register)$which.submodels[-1]) - 1)
  if (is.na(w) || length(w) ==0)
    stop("the value for 'which.submodels' does not match. Note that the definitoin has been changed in version 3.0.70")
  
  if (length(origin) != 1) stop("length of 'origin' is not 1.")
  if (is.character(origin)) {
    origin <- pmatch(origin, SORT_ORIGIN_NAMES) - 1
    if (is.na(origin)) stop("unallowed value of 'origin'")
  } else {
    if (!is.numeric(origin) || origin < 0 || origin >=length(SORT_ORIGIN_NAMES))
      stop("unallowed value for 'origin'")
  }
      
  cov <- .Call(C_GetModelInfo, as.integer(register), as.integer(level),
               as.integer(spConform), as.integer(w), as.integer(origin))

  if (!is.null(modelname)) {
    cov <- search.model.name(cov, modelname, 0)
  }
  class(cov) <- "RFgetModelInfo"
  attr(cov, "level") <- level
 
  return(cov)
}


RFgetModel <- function(register, explicite.natscale, show.call=FALSE,
		       origin = "original") {
  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register)) register else
                               deparse(substitute(register)))
  modus <- (if (missing(explicite.natscale)) GETMODEL_AS_SAVED else
            if (explicite.natscale)  GETMODEL_DEL_NATSC else
            GETMODEL_SOLVE_NATSC)

  which <- if (is.logical(show.call)) {
    if (show.call) "call+user" else "user"
  } else {
    show.call
  }
  
  if (is.character(origin)) {
    origin <- pmatch(origin, SORT_ORIGIN_NAMES) - 1
    if (is.na(origin)) stop("unallowed value of 'origin'")
  } else {
    if (!is.numeric(origin) || origin < 0 || origin >=length(SORT_ORIGIN_NAMES))
      stop("unallowed value for 'origin'")
  }
  m <- GetModel(register=register, modus=modus, which.submodels = which,
		as.integer(origin))
  class(m) <-  "RM_model"
  m
}
           
           
GetModel <- function(register, modus=GETMODEL_DEL_NATSC,
                     spConform=RFoptions()$general$spConform,      
                     which.submodels = c("user", "internal",
                         "call+user", "call+internal",
                         "user.but.once", "internal.but.once",
                         "user.but.once+jump", "internal.but.once+jump", "all"),
                     return.which.param=STANDARD,
                     ## STANDARD, ALLPARAMETERS, NOPARAMETERS,
                     ## INTERNALPARAMETER, INCLUDENOTRETURN,
                     solve_random = FALSE,
		     origin = original){
  
  ## modus:
  ##  AS_SAVED : Modell wie gespeichert
  ##  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  ##  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  ##               drauf multipliziert; Rest wie gespeichert)
  ##  DEL_MLE : nur natscale_MLE werden geloescht
  ##  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  ##               drauf multipliziert; Rest wie gespeichert)
  ## 

  ## do.notreturnparam : if true, then also parameters with the flag
  ##                      DONOTRETURN are returned

  ## spConform : only the names of the models
  ##              if > 1 than also "+" is forced to be "RMplus"

  
  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register))  register else
                               deparse(substitute(register)))
   w<-pmatch(match.arg(which.submodels),
            as.list(GetModel)$which.submodels[-1]) - 1
  if (missing(register)) register <- 0
  
  return(.Call(C_GetModel, as.integer(register), as.integer(modus),
	       as.integer(spConform), as.integer(w),
	       as.logical(solve_random), as.integer(return.which.param),
	       as.integer(origin)))
}

GetModelRegister <- function(name) { ## obsolete
  stopifnot(is.character(name))
  .C(C_GetModelRegister, name, integer(1))[[2]]
}


RFgetModelNames <- function(type = RC_TYPE_NAMES, domain = RC_DOMAIN_NAMES,
                            isotropy = RC_ISO_NAMES, operator = c(TRUE, FALSE),
                            monotone = RC_MONOTONE_NAMES,
                            implied_monotonicities = length(monotone) == 1,
                            finiterange = c(TRUE, FALSE, NA),
                            valid.in.dim = c(1, Inf), 
                            vdim = c(1, 5),
                            group.by,
			    exact.match = !missing(group.by),
                            simpleArguments = FALSE,
                            internal,
                            newnames
                            ){ #, .internal=FALSE

  group.names <- c("type", "domain", "isotropy", "operator",
		   "monotone", "finiterange", "valid.in.dim", "vdim")
  
  if (hasArg(internal)) {
    return(PrintModelList(operators=operator, internal = internal,
                          newstyle=missing(newnames) || newnames))
  }
  if (!missing(newnames) && !newnames) {
    if (hasArg(type) || hasArg(domain) || hasArg(isotropy) || hasArg(operator)
        || hasArg(monotone) || hasArg(finiterange) || hasArg(valid.in.dim)
        || hasArg(vdim) || hasArg(group.by) || hasArg(simpleArguments))
      stop("use argument 'newnames' without further arguments or in combination with 'internal'")
    return( .Call(C_GetAllModelNames, as.logical(newnames)) )
  }
  

  if (!(length(valid.in.dim) %in% 1:2)) stop("'valid.in.dim' has wrong size.")
  if (length(valid.in.dim) == 1) valid.in.dim <- c(valid.in.dim, Inf)
  
  if (!(length(vdim) %in% 1:2)) stop("'vdim' has wrong size.")
  if (length(vdim) == 1) vdim <- rep(vdim, 2)


##  Print(type, TYPE_NAMES)
  if (hasArg(type)) type <- TYPE_NAMES[pmatch(type, TYPE_NAMES)]
  if (hasArg(domain)) domain <- DOMAIN_NAMES[pmatch(domain, DOMAIN_NAMES)]
  if (hasArg(isotropy)) isotropy <- ISO_NAMES[pmatch(isotropy, ISO_NAMES)]
  if (hasArg(monotone)) monotone <- MONOTONE_NAMES[pmatch(monotone, MONOTONE_NAMES)]

 ## Print(type)
  if (any(is.na(pmatch(type, TYPE_NAMES))))
    stop(paste("'", type, "'", " is not a valid type category", sep=""))
  if (any(is.na(pmatch(domain, DOMAIN_NAMES))))
    stop(paste("'", domain, "'", " is not a valid domain category", sep=""))
  if (any(is.na(pmatch(isotropy, ISO_NAMES))))
    stop(paste("'", isotropy, "'", " is not a valid isotropy category", sep=""))
  if (any(is.na(pmatch(monotone, MONOTONE_NAMES))))
    stop(paste("'", monotone, "'", " is not a valid monotonicity category", sep=""))

 
  if (!exact.match) {
    if (length(domain) == 1 &&
	any(type %in% TYPE_NAMES[c(TcfType, PosDefType, VariogramType) + 1]) &&
	domain %in% DOMAIN_NAMES[1 + c(XONLY, KERNEL)])
      domain <- c(domain, DOMAIN_NAMES[c(PREVMODEL_D, SUBMODEL_D, PARAMDEP_D) + 1])
    
     if (length(type) == 1) {
      idx <- which(type ==  TYPE_NAMES[1 + c(TcfType, PosDefType, VariogramType,
					     NegDefType, ProcessType,
                                             GaussMethodType, 
                                             BrMethodType, SmithType,
                                             SchlatherType
                                             )])
       if (length(idx) > 0)
	type <- TYPE_NAMES[1 + list(c(TcfType, ManifoldType),
				    c(TcfType, PosDefType, ManifoldType),
				    c(TcfType, PosDefType, VariogramType,
                                      ManifoldType),
				    c(TcfType, PosDefType, VariogramType,
                                      NegDefType,
				      ManifoldType),
				    c(GaussMethodType, BrMethodType, SmithType,
				      SchlatherType, ProcessType),
				    c(BrMethodType, SmithType, SchlatherType)
				    )[[idx]] ]
    }
    
    group.by <- if (length(type) == 1) NULL else 'type'
    if (length(domain) > 1) group.by <- c(group.by, 'domain')
  }


  if (!missing(group.by) && length(group.by) > 0) {
    group.idx <- pmatch(group.by, group.names)
    if (any(is.na(group.idx)))
      stop("'group.by' can be equal to '",
           paste(group.names, collapse="', '"), "'")
    group.by <- group.names[group.idx]
    group.idx <- group.idx[1]
  }
   
  if (group <- !missing(group.by) && !is.null(group.by)) {
    FUN <- function(string) {
      args <- list(type=type, domain=domain, isotropy=isotropy,
                   operator=operator, monotone=monotone,
                   implied_monotonicities = implied_monotonicities &&
                   group.by[1] != "monotone",
                   finiterange=finiterange, valid.in.dim=valid.in.dim,
                   vdim=vdim,
                   group.by = if (group && length(group.by) > 1) group.by[-1] else NULL
                   )
      args[[group.idx]] <- string
      list(do.call("RFgetModelNames", args))
    }
    li <- sapply(get(group.by[1]), FUN=FUN)
    if (is.null(names(li))) names(li) <- paste(group.by[1], get(group.by[1]), sep="=")
    li <- li[sapply(li, FUN=length) > 0]
    return(li)
  } # matches  if (hasArg(group.by)) {}
  
  if (implied_monotonicities) {
    mon <- MONOTONE_NAMES[(UNSET-1) : -1]
    for (i in MONOTONE:NORMAL_MIXTURE)
      if (mon[i] %in% monotone) monotone <- c(monotone, mon[i + 1])
    if (mon[BERNSTEIN] %in% monotone) monotone <- c(monotone, mon[MONOTONE])
    monotone <- unique(monotone)
  }

  all <- list2RMmodel_Names
  internal <- substr(all, 1, 2) == "iR"
  idx <- logical(len <- length(all))
  for (i in 1:len) {
    fun <-  do.call(":::", list("RandomFields", all[i])) 
    idx[i] <- is.function(fun) && is(fun, class2=CLASS_RM)    
    if (!idx[i]) next
    idx[i] <- !all(is.na(pmatch(fun["type"], type, duplicates.ok=TRUE))) &&
 	      !all(is.na(pmatch(fun["isotropy"], isotropy, duplicates.ok=TRUE))) &&
	      !all(is.na(pmatch(domain, fun["domain"]))) &&
              fun["operator"] %in% operator &&
              !all(is.na(pmatch(monotone, fun["monotone"]))) &&
	      (!simpleArguments || fun["simpleArguments"]) &&
	      fun["finiterange"] %in% finiterange &&
	      (fun["maxdim"] < 0 ||
	      (fun["maxdim"] >= valid.in.dim[1] &&
	      fun["maxdim"] <= valid.in.dim[2]))  &&
	      (fun["vdim"] < 0 ||
	      (fun["vdim"] >= vdim[1] && fun["vdim"] <= vdim[2]))  
  }
#  Print(internal, substring(all[internal], 2))
  all[internal] <- substring(all[internal], 2)
  
  return(unique(sort(all[idx])))
}


RFgetMethodNames <- function() {
  RFgetModelNames(type=TYPE_NAMES[c(GaussMethodType, BrMethodType) + 1])
}


RFformula <- function(f) return(parseModel(f))

GetProcessType <- function(model) {
  stopifnot(is.list(model))
  .Call(C_GetProcessType, MODEL_INTERN, list("RFdummy", model))
}


parameter.range <- function(model, param, dim=1){
  cat("sorry not programmed yet\n")
  return(NULL)
 }

mergeWithGlobal <- function(dots) {
  if (exists(par.storage, envir=.RandomFields.env)) {
    current <- dots
    dots <- get(par.storage, envir=.RandomFields.env)
    names.current <- names(current)
    if (length(current) > 0) {
      for (i in 1:length(current))
        dots[[names.current[i]]] <- current[[i]]
    }
  }
  dots
}


RFpar <- function(...) {
  #l <- eval(substitute(list(...)))
  
  l <- list(...)
  if (length(l) == 1 && is.null(l[[1]]) && length(names(l)) ==0) {
    assign(par.storage, list(), envir=.RandomFields.env)
    return(NULL)
  }

  if (!exists(par.storage, envir=.RandomFields.env))
    assign(par.storage, list(), envir=.RandomFields.env)

  par <- get(par.storage, envir=.RandomFields.env)
  if (length(l) == 0) return(par)

  n <- names(l)
  for (i in 1:length(l)) {
    par[[n[i]]] <- l[[i]]
  }
  assign(par.storage, par, envir=.RandomFields.env) 
}
