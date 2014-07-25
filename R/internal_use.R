########################################################################
## Hier sind Funktionen, die vermutlich kaum ein Anwender interessiert
########################################################################

reps <- function(n, sign=",") paste(rep(sign, n), collapse="")

search.model.name <- function(cov, name, level) {
  #Print(cov, name, length(cov)); str(cov)
  if (length(name) == 0 || length(cov) ==0) return(cov);
 # Print("A", cov)
  if (!is.na(pmatch(name[1], cov))) return(search.model.name(cov, name[-1], 1))

  for (i in 1:length(cov$submodels)) {
    found <- search.model.name(cov$submodels[[i]], name, 1)
    if (!is.null(found)) return(found)      
  }
  found <- search.model.name(cov$key, name, 1)
  if (!is.null(found)) return(found)
  if (level == 0) stop("model name not found")
  return(NULL)
}



GetNeighbourhoods <- function(model.nr, all,
                              splitfactor, maxn, split_vec,
                              shared=FALSE)  {
  ## model.nr < 0 : standard scales
  ## MODEL_USER 0  /* for user call of Covariance etc. */
  ## MODEL_SIMU 1  /* for GaussRF etc */ 
  ## MODEL_INTERN 2 /* for kriging, etc; internal call of covariance */
  ## MODEL_SPLIT 3  /* split  covariance model */
  ## MODEL_GUI 4   /* RFgui */
  ## MODEL_MLE  5  /* mle covariance model */
  ## MODEL_MLESPLIT 6 /* ="= */
  ## MODEL_MLETREND 7  /* mle trend model !! */
  ## MODEL_BOUNDS 8 /* MLE, lower, upper */
  
  locfactor <- as.integer(splitfactor * 2 + 1) ## total number of
  ##    neighbouring boxes in each direction, including the current box itself
  maxn <- as.integer(maxn)        ## max number of points, including neighbours
  splitn <- as.integer(split_vec[1]) ## number of location when split up
  minimum <- as.integer(split_vec[2])## min. number of points in a neighbourhood
  maximum <- as.integer(split_vec[3])## maximum number of points when still
  ##                                    neighbours of neighbours are included.
  ##                         Note that, mostly, an additional box is included.

  #  Print(RFoptions(), locfactor, maxn, splitn, minimum, maximum)

  d <- dim(all$given)
  ts.xdim <- as.integer(d[1])
  n <- as.integer(d[2])

  if (model.nr >= 0) {
    natsc <- double(ts.xdim)
    .C("MultiDimRange", as.integer(model.nr), natsc,
       DUP=DUPFALSE, PACKAGE="RandomFields")
  } else natsc <- rep(1, ts.xdim)

 
  u <- numeric(ts.xdim)
  for (i in 1:ts.xdim) {
    u[i] <- length(unique(all$given[i,]))
  }
  
  Range <- apply(all$given, 1, range)
  rd <- apply(Range, 2, diff) 
  len <- pmax(1e-10 * max(rd), rd * (1 + 1e-10))
  units <- pmax(1, len * natsc)
  nDsplitn <- n / splitn

  ## * "gerechte" Aufteilung in alle Richtungen waere nDsplitn
  ## * die Richtung in die viele units sind, soll eher aufgespalten werden
  ## * ebenso : wo viele Werte sind eher aufspalten
  idx <- (nDsplitn / prod(units * u))^{1/ts.xdim} * units * u > 0.5 
  reddim <- sum(idx)
  units <- units[idx]
  zaehler <- 1
  parts <- rep(1, ts.xdim)
  OK <- integer(1)
  
  repeat {
    parts[idx] <- (nDsplitn / prod(units))^{1/reddim} *
                  locfactor * zaehler * units * all$vdim
    parts <- as.integer(ceiling(parts))

    ## zuordnung der coordinaten_Werte zu den jeweiligen "parts"
    ## given ist liegend
    coord.idx <- floor((all$given - Range[1,]) / (len / parts))
    
    cumparts <- cumprod(parts)
    totparts <- cumparts[length(cumparts)]
    Ccumparts <- as.integer(c(1, cumparts))
    cumparts <- Ccumparts[-length(Ccumparts)]
    
    elms.in.boxes <- integer(totparts)    
    neighbours <- integer(totparts)
    cumidx <- as.integer(colSums(coord.idx * cumparts))

 
    .C("countelements", cumidx, n, elms.in.boxes, 
       DUP=DUPFALSE, PACKAGE="RandomFields")

     
    .C("countneighbours", ts.xdim, parts, locfactor, Ccumparts,
       elms.in.boxes, neighbours, OK,
       DUP=DUPFALSE, PACKAGE="RandomFields")

    ## if there too many points within all the neighbours, then split
    ## into smaller boxes
    zaehler <- zaehler * 2

    ## image(neighbours, zlim=c(0:(prod(parts)-1)))
    if (OK) break;
  }
    

#  Print(elms.in.boxes, neighbours)

  l <- list()
  l[[1]] <- .Call("getelements", cumidx, ts.xdim, n, Ccumparts, elms.in.boxes,  
                  PACKAGE="RandomFields")
  l1len <- sapply(l[[1]], length)

  if (length(all$x) > 0) {
    ## now calculate the boxes for the locations where we will interpolate
    i <- pmax(0, pmin(parts-1, floor((t(all$x) - Range[1,]) / (len / parts)) ))
    dim(i) <- rev(dim(all$x))
    i <- as.integer(colSums(i * cumparts))
    res.in.boxes <- integer(totparts)
    
    .C("countelements", i, nrow(all$x), res.in.boxes,
       DUP=DUPFALSE, PACKAGE="RandomFields")
    
    notzeros <- res.in.boxes > 0
    l[[3]] <- .Call("getelements", i, ts.xdim, as.integer(nrow(all$x)),
                    Ccumparts, res.in.boxes, PACKAGE="RandomFields")[notzeros]
  } else {
    notzeros <- TRUE
  }

  ll <- .Call("getneighbours", ts.xdim, parts, locfactor, Ccumparts,
              neighbours, PACKAGE="RandomFields")[notzeros]
  less <- sapply(ll, function(x) sum(elms.in.boxes[x]) < minimum) | !shared
  ##                  if !shared then all(less)==TRUE

 # Print(l, less, minimum,sapply(ll, function(x) sum(elms.in.boxes[x])),
  #sapply(ll, length), elms.in.boxes, sum(elms.in.boxes)) ; 
  
  if (any(less)) {
    not.considered.yet <- sapply(l[[1]], length) > 0   
    newll <- ll
    for (i in which(less)) {
      current <- ll[[i]]
      elements <- sum(elms.in.boxes[current] *
                      (shared | not.considered.yet[current]))# number of pts in a neighbourhood
      while (elements < minimum) {
        new <- unique(unlist(ll[current])) # neighbours of neighbours, but not
        new <- new[which(is.na(pmatch(new, current)))]# neighbours themselves
        nn <- elms.in.boxes[new] * (shared | not.considered.yet[new]) # how many pts are in each of these boxes?
        ordr <- order(nn)
        new <- new[ordr]
        nn <- nn[ordr]
        cs <- elements + cumsum(nn)
        smaller <- sum(cs <= maximum) ## now, check which neighbours of
        ## the neigbours can be included in the list of neighbours of i
        ## to increase the number of points in the kriging neighbourhood
        if (smaller == 0) break; ## none
        if (smaller == length(cs) || cs[smaller] >= minimum ||
            cs[smaller+1] > maxn) {
          if ( (elements <- cs[length(cs)]) <= maxn ) {            
            current <- c(current, new)            
          } else {
            current <- c(current, new[1:smaller])
            elements <- cs[smaller]
          }
          if (smaller != length(cs)) break
        } else {
          ## smaller < length(cs) && cs[smaller]<minimum && cs[smaller+1]<=maxn
          ## i.e., include the next one, but there is no chance to include
          ## more within the rules.
          elements <- cs[smaller+1]
          current <- c(current, new[1:(smaller+1)])
          break;
        }
      }
      current <- current[l1len[current] > 0]
      if (!shared) current <- current[not.considered.yet[current]]
      newll[[i]] <- current
      not.considered.yet[current] <- FALSE                            
    }
    newll <- newll[sapply(newll, length) > 0]
    l[[2]] <- newll
  } else l[[2]] <- ll
 #l[[4]] <- list(Range=Range, len=len, parts=parts, cumparts=Ccumparts,
 #               locfactor=locfactor, nDsplitn=nDsplitn)

  # Print(l)
 # Print(sapply(newll, length), sapply(newll, function(x) sum(elms.in.boxes[x])), length(ll), length(newll), l, not.considered.yet, l[not.considered.yet]) ;
  #hhh

# str(l)
 # cccc
  return(if (shared) l else lapply(l[[2]], function(x) unlist(l[[1]][x])))
}

resolve.register <- function(register){
  stopifnot(!missing(register))
  if (length(register) == 0) {
    register <- .C("GetCurrentRegister", reg=integer(1))$reg
    if (register < 0) stop("model that has been used right now cannot be determined. This happens, for instance, after the use of 'RFfit'")
  }
  if (!is.numeric(register)) {
 #   register <- deparse(substitute(register))   
    register <-
      switch(register,
             "RFcov" = MODEL.USER,
             "RFdistr" = MODEL.USER,
             "RFcovmatrix" = MODEL.USER,
             "RFvariogram" = MODEL.USER,
             "RFpseudovariogram" =  MODEL.USER,
             "RFsimulate" =  RFoptions()$registers$register,
             "RFinterpolate" =  RFoptions()$registers$interpolregister,
             "RFgui" =  GetModelRegister("gui"),
             "RFfit" =  MODEL.MLE,
             "RFratiotest" =  MODEL.MLE,
             stop("register unknown")
             )
  }
  stopifnot(is.numeric(register))
  if (register < 0) stop("'register' does not have a non-negative value.")
  return(register)
}

RFgetModelInfo <-
  function(register, level=3, 
           spConform=RFoptions()$general$spConform,
           which.submodels = c("submodels", "keys", "both"),
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
               as.integer(if (which.submodels=="submodels") 0 else
                          1 + (which.submodels=="both")),
               PACKAGE="RandomFields")

  if (!is.null(modelname)) {
    cov <- search.model.name(cov, modelname, 0)
  }

  return(cov)
}


RFgetModel <- function(register, explicite.natscale, show.call=FALSE)  {
  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register)) register else
                               deparse(substitute(register)))
  modus <- (if (missing(explicite.natscale)) GETMODEL_AS_SAVED else
            if (explicite.natscale)  GETMODEL_DEL_NATSC else
            GETMODEL_SOLVE_NATSC)
  if (show.call) modus <- modus + 10
  m <- GetModel(register=register, modus=modus)
  class(m) <-  "RM_model"
  m
}
           
           
GetModel <- function(register, modus=GETMODEL_DEL_NATSC,
                     spConform=RFoptions()$general$spConform,
                     do.notreturnparam=FALSE,
                     replace.select = FALSE) {
  ## modus:
  ##  AS_SAVED : Modell wie gespeichert
  ##  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  ##  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  ##               drauf multipliziert; Rest wie gespeichert)
  ##  DEL_MLE : nur natscale_MLE werden geloescht
  ##  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  ##               drauf multipliziert; Rest wie gespeichert)
  ## 
  ## modus: 10+ : wie oben, jedoch ohne CALL_FCT zu loeschen 

  ## do.notreturnparam : if true, then also parameters with the flag
  ##                      DONOTRETURN are returned

  ## spConform : only the names of the models

  ## replace.select : if TRUE then model "select " is replaced by model
  ##                  model "plus" -- they should be joined anyway
  

  register <- resolve.register(if (missing(register)) NULL else
                               if (is.numeric(register))  register else
                               deparse(substitute(register)))
  if (missing(register)) register <- 0
  
  model <- .Call("GetModel", as.integer(register), as.integer(modus),
                 as.integer(spConform), as.integer(do.notreturnparam),
                 PACKAGE="RandomFields")

  if (replace.select) {
    if (model[[1]] == ZF_SELECT[1]) {
      model[[1]] <- ZF_SYMBOLS_PLUS
      stopifnot(names(model)[2]=="subnr")
      model[[2]] <- NULL
    }    
  }
  return(model)
}


GetPracticalRange <- function(model, param, dim=1) {
  ## dim=spdim=tsdim
  model <- PrepareModel(model, param)
#  userdefined <- GetParameterModelUser(model)
  InitModel(MODEL.USER, list("Dummy", model), dim) 
  natscl <- double(1)
  .C("UserGetNatScaling", natscl, PACKAGE="RandomFields", DUP=DUPFALSE)
  .C("DeleteKey", MODEL.USER)# to do : nicht sauber
  return(1.0 / natscl)
}

#GetPracticalRange("whittle", 1)

GetrfParameters <- function(initcov=TRUE){
  ## if initcov then InitModelList is called before
  ## values are return (necessary to get covnr right, but should
  ## not be done for RFoptions (debugging reasons)

  maxints <- integer(1)
  .C("GetMaxDims", maxints, DUP=DUPFALSE, PACKAGE="RandomFields")
  name <- c("GetrfParameters", "GetrfParametersI")[1 + (initcov != 0) ]
  
  p <- .C(name, covmaxchar=integer(1), methodmaxchar=integer(1),
          distrmaxchar=integer(1),
          covnr=integer(1), methodnr=integer(1), distrnr=integer(1),
          maxdim=integer(maxints), maxmodels=integer(1),
          type = integer(1),
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
  ## verwendet??
  assign(".p", GetrfParameters(TRUE))
  l <- character(.p$distrnr)
  for (i in 1:.p$distrnr) {
    l[i] <- .C("GetDistrName", as.integer(i-1),
               n=paste(rep(" ",.p$distrmaxchar), collapse=""),
               PACKAGE="RandomFields")$n
  }
  return(l)
}



PrintModelList <-function (operators=FALSE, internal=FALSE,
                           newstyle=TRUE) {
   stopifnot(internal <=2, internal >=0, operators<=1, operators>=0)
    .C("PrintModelList", as.integer(internal), as.integer(operators),
       as.integer(newstyle),
       PACKAGE="RandomFields")
    invisible()
}


PrintMethodList <-function () {
    .C("PrintMethods", PACKAGE="RandomFields")
    invisible()
}



distInt <- function(x) {
  ##
  ## only for gene data where each coordinates takes
  ## only three neighboured integer values !
  stopifnot(is.matrix(x), is.integer(x))
  n <- nrow(x)
  genes <- ncol(x)
  res <- double(n * n)
  .C("distInt", t(x), n, genes, res, PACKAGE="RandomFields", DUP=DUPFALSE)
  dim(res) <- rep(n, 2)
  res
}




parameter.range <- function(model, param, dim=1){
  cat("sorry not programmed yet\n")
  return(NULL)
 

#parampositions < - function(model, param, trend=NULL, dim, print=1) {
#  stopifnot(!missing(dim))
#  model <- PrepareModel(model, param, trend=trend, nugget.remove=FALSE)
#  .Call("Get NA Positions", reg, model, as.integer(dim), as.integer(dim),
#        FALSE, FALSE, as.integer(print), PACKAGE="RandomFields")
#}
#  pm <- PrepareModel(model=model, param=param, nugget.remove=FALSE)        
#  storage.mode(dim) <- "integer"
#  ResGet <- .Call("SetAnd  ?? GetModelInfo",
#                  reg,
#                  pm, dim, Time, dim, FALSE, MaxNameCharacter,
#                  TRUE, TRUE,
#                  PACKAGE="RandomFields")
#  minmax <- ResGet$minmax[, 1:2]
#  dimnames(minmax) <-
#    list(attr(ResGet$minmax, "dimnames")[[1]], c("min", "max"))
#  return(minmax)
}

GetModelRegister <- function(name) {
  stopifnot(is.character(name))
  return(as.integer(.C("GetModelRegister", name, integer(1),
                       PACKAGE="RandomFields")[[2]]))
}


cmplists <- function(l, m) {
  stopifnot(is.list(l) || is.list(m))
  if (length(l) == 0) {
    if (length(m)==0) return(c("empty lists =", names(l)))
    else return(list("empty first list, but not second:", m))
  } else if (length(m) == 0)
    return(list("empty second list, but not first:", m))
  n <- list()
  for (i in 1:min(length(m), length(l))) {
     if (xor(is.list(l[[i]]), is.list(m[[i]])) ||
        xor(is.numeric(l[[i]]), is.numeric(m[[i]])) ||
        xor(is.logical(l[[i]]), is.logical(m[[i]])) ||
        xor(is.matrix(l[[i]]), is.matrix(m[[i]])) ||
        xor(is.array(l[[i]]), is.array(m[[i]])) ||
        xor(is.character(l[[i]]), is.character(m[[i]]))
        ) {
      n[[i]] <- list(first = l[[i]], second = m[[i]])
    } else {
      if (is.list(l[[i]])) n[[i]] <- cmplists(l[[i]], m[[i]]) else 
      if (is.numeric(l[[i]]) | is.logical(l[[i]])) {
        if (length(l[[i]]) != length(m[[i]]) ||
            is.array(l[[i]]) && any(dim(m[[i]]) != dim(n[[i]])))
          n[[i]] <- list("differ in size", names(l)[[i]],
                         first=l[[i]], second=l[[i]])
        else {
          idx <- is.na(l[[i]]) != is.na(m[[i]]) |
                 (!is.na(l[[i]]) && l[[i]] != m[[i]])
        #  Print(idx, l[[i]] , m[[i]]  )
          n[[i]] <-
            if (!any(idx))
              c(if (is.numeric(l[[i]])) "num =" else "logi=", names(l)[[i]])
            else list(c(names(l)[i], which(idx)), l[[i]], m[[i]])
        }
      } else
      if (is.character(l[[i]])) {     
        idx <-  l[[i]] != m[[i]]
     #   Print(idx, l[[i]] , m[[i]])
        n[[i]] <- if (!any(idx)) c("char =", names(l)[[i]], l[[i]])
                  else list(c(names(l)[i], which(idx)), l[[i]], m[[i]])
      } else {
        n[[i]] <- c("cannot be compared", names(l)[[i]])
      }
    }
  }
  if (length(l) > length(m))
    for (i in (length(m)+1):length(l))
      n[[i]] <- list("only first list has entry", l[[i]])
   if (length(m) > length(l))
    for (i in (length(l)+1):length(m))
      n[[i]] <- list("only second list has entry", m[[i]])
  return(n)
}

checkExamples <- function(exclude=NULL, include=1:length(.fct.list),
                          ask=FALSE, echo=TRUE, halt=FALSE, ignore.all=FALSE,
                          path="RandomFields", package="RandomFields",
                          read.rd.files=TRUE) {
  .exclude <- exclude
  .ask <- ask
  .echo <- echo
  .halt <- halt
  .ignore.all <- ignore.all
  .path <- path
  .package <- package
  useDynLib <- importClassesFrom <- import <-
    importFrom <- exportClasses <-
      importMethodsFrom <- exportMethods <- S3method <- function(...) NULL
  .env <- new.env()
  if (!require(.package, character.only=TRUE)) return(NA)
  exportPattern <- function(p) {
    all.pattern <- p %in% c("^[^\\.]", "^[^.]", ".") | get("all.pattern", .env)
    if (!.ignore.all) assign("all.pattern", all.pattern, .env)
    if (all.pattern) return(NULL)
 #   Print(p)
    stopifnot(nchar(p)==2, substr(p,1,1)=="^")

 #   Print(substr(p, 2, 1), c(get("p", .env)))
    
     assign("p", c(get("p", .env), substring(p, 2)), .env)
 #   Print(get("p", .env))
     
  }
  export <- function(...) {
    ## code from 'rm'
    dots <- match.call(expand.dots = FALSE)$...
    z <-deparse(substitute(...))
    if (length(dots) && !all(sapply(dots, function(x) is.symbol(x) || 
                                    is.character(x)))) 
      stop("... must contain names or character strings")
    z <- sapply(dots, as.character)
    assign("export", c(get("export", .env), z), .env)
  }
  assign("export", NULL, .env)
  assign("all.pattern", FALSE, .env)
  assign("p", NULL, .env)
  source(paste(.path, "NAMESPACE", sep="/"), local=TRUE)
  if (is.logical(read.rd.files) && !read.rd.files) {
    .package.env <- parent.env(.GlobalEnv)
    while (attr(.package.env, "name") != paste("package:", .package, sep="")) {
      ##Print(attr(.package.env, "name"), paste("package:", .package, sep=""))
      .package.env <- parent.env(.package.env)
    }
    .orig.fct.list <- ls(envir=.package.env)
    .ok <- (get("all.pattern", .env) |
            substr(.orig.fct.list, 1, 1) %in% get("p", .env) | 
            .orig.fct.list %in% get("export", .env))
    .fct.list <- .orig.fct.list[.ok]
  } else {
    .path <- read.rd.files
    if (is.logical(.path))
      .path <- "/home/schlather/svn/RandomFields/RandomFields/man"
    .files <- dir(.path, pattern="d$")
    .fct.list <- character(length(.files))
    for (i in 1:length(.files)) {
      .content <- scan(paste(.path, .files[i], sep="/") , what=character(),
                       quiet=TRUE)
      .content <- strsplit(.content, "alias\\{")
      .content <- .content[which(lapply(.content, length) > 1)][[1]][2]
      .fct.list[i] <-
        strsplit(strsplit(.content,"\\}")[[1]][1], ",")[[1]][1]
    }
  }
  Print(.fct.list, length(.fct.list)) #
  .include <- include
  .RFopt <- RFoptions()
  .not_working_no <- .not_working <- NULL
  .included.fl <- .fct.list[.include]
  .not.isna <- !is.na(.included.fl)
  .include <- .include[.not.isna]
  .included.fl <- .included.fl[.not.isna]
  Print(.included.fl, get("export", .env), get("p", .env)); #return(NA)
  .max.fct.list <- max(.included.fl)
  for (.idx in .include) {
    if (.idx %in% .exclude) next
    cat("\n\n\n\n\n", .idx, " ", .package, ":", .fct.list[.idx],
        " (total=", length(.fct.list), ") \n", sep="")
    RFoptions(LIST=.RFopt)
    .C("ResetWarnings")
    if (.echo) cat(.idx, "")
    .tryok <- TRUE
    if (.halt) {
      do.call(utils::example, list(.fct.list[.idx], ask=.ask, echo=.echo))
    } else {
      if (is(try(do.call(utils::example, list(.fct.list[.idx], ask=.ask,
                                         echo=.echo))), "try-error")) {
        .not_working_no<- c(.not_working_no, .idx)
        .not_working <- c(.not_working, .fct.list[.idx])
        .tryok <- FALSE
      }
    }
    RFoptions(storing = FALSE);
    cat("****** '", .fct.list[.idx], "' (", .idx, ") done. ******\n")
    if (.tryok && !is.na(RFoptions()$general$seed)) {
      Print(.not_working, paste(.not_working_no, collapse=", ")) #
      stop("seed not NA: ", .fct.list[.idx])
    }
  }
  Print(.not_working, paste(.not_working_no, collapse=", ")) #
  .ret <- list(.not_working, .not_working_no)
  names(.ret) <- c(.package, "")
  return(.ret)
}


FinalizeExample <- function() {
  if (!interactive()) {
    close.screen(all.screens = TRUE)
    n <- length(dev.list()) - 5
    if (n > 0) {
      for (i in 1:n) {      
        dev.off() ## OK
      }
    }
  }
  RFoptions(seed = NA)
}


Dependencies <- function(install= all(pkgs == c(depends, imports, suggests)),
                         check=TRUE,
                         pkgs = c(depends, imports, suggests),
                         lib.loc ="/usr/local/lib64/R/library" ) {
  depends <- c("DSpat", "Geneland", "GeoGenetix", "LS2Wstat", "MarkedPointProcess", "ProbForecastGOP", "Sunder")
  imports <- c("constrainedKriging", "geoR", "georob", "lgcp", "spatsurv")
  suggests <- c("CompRandFld", "fractaldim", "geostatsp", "rpanel", "spatstat")
 
  PKGS <- character(length(pkgs))
  ip <- installed.packages(lib.loc = lib.loc)
  
  rn <- rownames(ip)
  for (i in 1:length(pkgs)) {
    idx <- which(pkgs[i] == rn)
    v <- ip[idx, "Version"]
    PKGS[i] <- paste(pkgs[i], "_", v, ".tar.gz", sep="")    
  }
 
  if (install) {
    for (i in 1:length(pkgs)) {
      x <- system(paste("(cd Dependencies; wget http://cran.r-project.org/src/contrib/", PKGS[i], ")", sep=""))
      if (x != 0) {
        install.packages(pkgs[i], repos="http://ftp5.gwdg.de/pub/misc/cran/")
        ip <- installed.packages(lib.loc = lib.loc)
        rn <- rownames(ip)
        idx <- which(pkgs[i] == rn)
        v <- ip[idx, "Version"]
        PKGS[i] <- paste(pkgs[i], "_", v, ".tar.gz", sep="")
        x <- system(paste("(cd Dependencies; wget http://cran.r-project.org/src/contrib/", PKGS[i], ")", sep=""))
        stopifnot(x == 0)
      }
    }
  }
  
  if (check) {
    for (i in 1:length(pkgs)) {
      x <- system(paste("(cd Dependencies; R CMD check --as-cran", PKGS[i],")"))
      if (x != 0) stop(PKGS[i], "failed")
    }
  }
}

showManpages <- function(path="/home/schlather/svn/RandomFields/RandomFields/man") {
  files <- dir(path)
  result <- character(length(files))
  for (i in 1:length(files)) {
    cat("\n\n\n\n\n")
    system(paste("grep \"\\examples\" -A 10000 ", path, "/", files[i], sep=""))
    result[i] <- readline(files[i])
  }
  result <- cbind(files, result)
  result <- result[result[, 2]!= "", ]
  result[result[, 2]=="t", 2] <- "To Do"
  result[result[, 2]=="s", 2] <- "dontshow"
  return(result)
}
# showManpages()
