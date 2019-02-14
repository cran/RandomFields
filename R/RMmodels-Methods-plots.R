
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017 Alexander Malinowski & Martin Schlather
## Copyright (C) 2017 -- 2017 Martin Schlather
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


## accessing CLASS_CLIST and RMmodelgenerator via '['-operator
# e.g. RMwhittle["domain"]

accessSlotsByName <- function(x,i,j,drop=FALSE) {
  names <- slotNames(x)
  if (!(i %in% names))
    stop(paste(i, "is not a valid slot specification"))
  return(slot(x, i))
}

accessReplaceSlotsByName <- function(x,i,j,value) {
  names <- slotNames(x)
  if (!(i %in% names))
    stop(paste(i, "dis not a valid slot specification"))
  else
    slot(x, i) <- value
  validObject(x)
  return(x)
}

setMethod(f="[", signature = CLASS_CLIST, def = accessSlotsByName)
setMethod(f="[", signature = CLASS_FIT, def = accessSlotsByName)
setMethod(f="[", signature = CLASS_RM, def = accessSlotsByName)

setReplaceMethod(f="[", signature = CLASS_CLIST, accessReplaceSlotsByName)
setReplaceMethod(f="[", signature = CLASS_FIT, accessReplaceSlotsByName)
setReplaceMethod(f="[", signature = CLASS_RM,accessReplaceSlotsByName)


## summing up RMmodels
setMethod('c', signature=c(CLASS_CLIST),
          function(x, ..., recursive = FALSE) R.c(x, ...) )


## funktioniert im speziellen Modus nicht, allgemein schon
## warum ist unklar
#setMethod('diag', signature=c(CLASS_CLIST, 'missing', 'missing', 'missing'),
#	  function(x=1, nrow, ncol, names=TRUE) {
#            ## Print("hier");  str(x)
#            RMmatrix(x)
#          })
          


resolve <- function(e1, e2, sign) {
  ## RMmodels only
  d <- list()
  if (e1@name==sign && (len.e1 <- length(e1@submodels)) < MAXSUB) {             
    for (i in 1:len.e1)  d[[i]] <- e1@submodels[[i]]
  } else {
    len.e1 <- 1
    d[[1]] <- if (is.character(e1))
      stop("characters cannot be combined with 'RMmodels'")
                                        #  do.call(SYMBOL_P, list(e1))
    else e1
  }
  d[[len.e1 + 1]] <-
     (if (is.character(e2))
       stop("characters cannot be combined with 'RMmodels'")
                                        #  do.call(SYMBOL_P, list(e2))
     else e2)
  if (sign == RM_MULT[1]) {
#    print(d)
    model <- do.call(sign, d)
  } else model <- do.call(sign, d)
  return(model)
}

warn.resolve.txt <- "A large vector consists fully of NAs -- the model is probably not correct.\nNote that it is always better to define the covariance model in the first\nsummands and then the trend. Also better use explicitely 'R.c', 'RMcovariate'\nand 'R.const' if the model is more complicated"

resolveRight<- function(e1, e2, sign) {
  ## left constant
  d <- list()
  if (e1@name==sign && (len.e1 <- length(e1@submodels)) < MAXSUB) {
    for (i in 1:len.e1)  d[[i]] <- e1@submodels[[i]]
  } else {
    len.e1 <- 1
    d[[1]] <- (if (is.character(e1)) 
       stop("characters cannot be combined with 'RMmodels'")
                                        #do.call(SYMBOL_P, list(e1))
      else e1)
  }
  d[[len.e1 + 1]] <-
    if (is.list(e2)) do.call(R_C, e2) 
    else if (length(e2)==1) do.call(R_CONST, list(e2))
    else if (!all(is.finite(e2))) {
      if (all(is.na(e2)) && length(e2)>5) warning(warn.resolve.txt)
      do.call(R_C, list(e2))
    } else if (sign == RM_PLUS[1]) {
      tmpList <- list(RM_COVARIATE)
      tmpList[[COVARIATE_C_NAME]] <- e2
      tmpList[[COVARIATE_X_NAME]] <- NULL
      tmpList[[COVARIATE_ADDNA_NAME]] <- TRUE
      do.call(RM_COVARIATE, tmpList)
    } else do.call(R_C, list(e2))
  ## Version brasilienc
 ##  if (length(e2)==1) d[[len.e1 + 1]] <- do.call(R_CONST, list(e2))
##  else {
##    e <- list(e2)
##    e[[COVARIATE_ADDNA_NAME]] <- TRUE
##    d[[len.e1 + 1]] <- do.call(RM_COVARIATE, e)
##  }

  model <- do.call(sign, d)
  return(model)
}

resolveLeft<- function(e1, e2, sign) {
  ## right constant
  d <- list()
  len.e1 <- 1
## Version Brasilien  
# if (length(e1)==1)  d[[1]] <- do.call(R_CONST, list(e1))
#  else {
#    e <- list(e1)
#    e[[COVARIATE_ADDNA_NAME]] <- TRUE   
#    d[[1]] <- do.call(RM_COVARIATE, e)
#  }

  d[[1]] <- if (is.list(e1)) do.call(R_C, e1)  
    else if (length(e1)==1) do.call(R_CONST, list(e1))
    else if (!all(is.finite(e1))) {
      if (all(is.na(e1)) && length(e1)>5) warning(warn.resolve.txt)
      do.call(R_C, list(e1))
    } else if (sign == RM_PLUS[1]) {
      tmpList <- list(RM_COVARIATE)
      tmpList[[COVARIATE_C_NAME]] <- e1
      tmpList[[COVARIATE_X_NAME]] <- NULL
      tmpList[[COVARIATE_ADDNA_NAME]] <- TRUE
      do.call(RM_COVARIATE, tmpList)
    } else do.call(R_C, list(e1))
  
  d[[len.e1 + 1]] <-  (if (is.character(e2)) 
			 stop("characters cannot be combined with 'RMmodels'")
		       ##do.call(SYMBOL_P, list(e2))
		       else e2)
  model <- do.call(sign, d)
  return(model)
}

## summing up RMmodels
setMethod('+', signature=c(CLASS_CLIST, CLASS_CLIST),
          function(e1, e2) resolve(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c(CLASS_CLIST, 'numeric'), 
          function(e1, e2) resolveRight(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c(CLASS_CLIST, 'logical'),
          function(e1, e2) resolveRight(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c(CLASS_CLIST, 'factor'),
          function(e1, e2) resolveRight(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c(CLASS_CLIST, 'list'),
          function(e1, e2) resolveRight(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c('numeric', CLASS_CLIST),
          function(e1, e2) resolveLeft(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c('logical', CLASS_CLIST),
          function(e1, e2) resolveLeft(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c('data.frame', CLASS_CLIST),
          function(e1, e2) resolveLeft(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c('factor', CLASS_CLIST),
          function(e1, e2) resolveLeft(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c(CLASS_CLIST, 'character'),
          function(e1, e2) resolve(e1, e2, RM_PLUS[1]))
setMethod('+', signature=c('character', CLASS_CLIST),
          function(e1, e2) resolve(e1, e2, RM_PLUS[1]))


## multiplying RMmodels
setMethod('*', signature=c(CLASS_CLIST, CLASS_CLIST),
         function(e1, e2) resolve(e1, e2, RM_MULT[1]))
setMethod('*', signature=c('numeric', CLASS_CLIST),
         function(e1, e2) resolveLeft(e1, e2, RM_MULT[1]))
setMethod('*', signature=c('logical', CLASS_CLIST),
         function(e1, e2) resolveLeft(e1, e2, RM_MULT[1]))
setMethod('*', signature=c(CLASS_CLIST, 'logical'),
          function(e1, e2) resolveRight(e1, e2, RM_MULT[1]))
setMethod('*', signature=c(CLASS_CLIST, 'numeric'),
          function(e1, e2) resolveRight(e1, e2, RM_MULT[1]))
setMethod('*', signature=c(CLASS_CLIST, 'character'),
         function(e1, e2) resolve(e1, e2, RM_MULT[1]))
setMethod('*', signature=c('character', CLASS_CLIST),
         function(e1, e2) resolve(e1, e2, RM_MULT[1]))


Xresolve <- function(e1, e2, model) {
  d <- list()
  len.e1 <- 1
  d[[1]] <-
    (if (is.character(e1)) 
     stop("characters cannot be combined with 'RMmodels'")
                                        #do.call(SYMBOL_P, list(e1))
    else e1)
  d[[len.e1 + 1]] <-
    ( if (is.character(e2)) 
      stop("characters cannot be combined with 'RMmodels'")
                                        #do.call(SYMBOL_P, list(e2))
    else e2)
  model <- do.call(model, d)
  return(model)
}
XresolveLeft <- function(e1, e2, model) {
  d <- list()
  len.e1 <- 1
  if (length(e1)==1) d[[1]] <- do.call(R_CONST, list(e1))
  else {
    e <- list(e1)
    d[[1]] <-  do.call(R_CONST, e)
## version brasilien    
##    e[[COVARIATE_ADDNA_NAME]] <- TRUE   
##    d[[1]] <-  do.call(RM_COVARIATE, e)
 }
  d[[len.e1 + 1]] <- (if (is.character(e2)) 
			stop("characters cannot be combined with 'RMmodels'")
                                        #do.call(SYMBOL_P, list(e2))
		      else e2)       
  model <- do.call(model, d)
  return(model)
}
XresolveRight <- function(e1, e2, model) {
  d <- list()
  len.e1 <- 1
  d[[1]] <- (if (is.character(e1)) 
     stop("characters cannot be combined with 'RMmodels'")
                                        #do.call(SYMBOL_P, list(e1))
    else e1)
  if (length(e2)==1) d[[len.e1 + 1]] <- do.call(R_CONST, list(e2))
  else {
    e <- list(e2)
    d[[len.e1 + 1]] <- do.call(R_CONST, e)
  ## version brasilien
  ## e[[COVARIATE_ADDNA_NAME]] <- TRUE   
  ## d[[len.e1 + 1]] <- do.call(RM_COVARIATE, e)
  }
  model <- do.call(model, d)
  return(model)
}
  

## substracting RMmodels
setMethod('-', signature=c(CLASS_CLIST, CLASS_CLIST),
          function(e1, e2) Xresolve(e1, e2, "R.minus"))
setMethod('-', signature=c('numeric', CLASS_CLIST),
          function(e1, e2) XresolveLeft(e1, e2, "R.minus"))
setMethod('-', signature=c('logical', CLASS_CLIST),
          function(e1, e2) XresolveLeft(e1, e2, "R.minus"))
setMethod('-', signature=c(CLASS_CLIST, 'numeric'),
          function(e1, e2) XresolveRight(e1, e2, "R.minus"))
setMethod('-', signature=c(CLASS_CLIST, 'logical'),
          function(e1, e2) XresolveRight(e1, e2, "R.minus"))
setMethod('-', signature=c(CLASS_CLIST, 'character'),
          function(e1, e2) Xresolve(e1, e2, "R.minus"))
setMethod('-', signature=c('character', CLASS_CLIST),
          function(e1, e2) Xresolve(e1, e2, "R.minus"))

## dividing up RMmodels
setMethod('/', signature=c(CLASS_CLIST, CLASS_CLIST),
          function(e1, e2) Xresolve(e1, e2, "R.div"))
setMethod('/', signature=c('numeric', CLASS_CLIST),
          function(e1, e2) XresolveLeft(e1, e2, "R.div"))
setMethod('/', signature=c('logical', CLASS_CLIST),
          function(e1, e2) XresolveLeft(e1, e2, "R.div"))
setMethod('/', signature=c(CLASS_CLIST, 'numeric'),
          function(e1, e2) XresolveRight(e1, e2, "R.div"))
setMethod('/', signature=c(CLASS_CLIST, 'logical'),
          function(e1, e2) XresolveRight(e1, e2, "R.div"))
setMethod('/', signature=c('RMmodel', 'character'),
          function(e1, e2) Xresolve(e1, e2, "R.div"))
setMethod('/', signature=c('character', 'RMmodel'),
          function(e1, e2) Xresolve(e1, e2, "R.div"))


## powering up RMmodels
setMethod('^', signature=c(CLASS_CLIST, CLASS_CLIST),
          function(e1, e2) Xresolve(e1, e2, "R.pow")) 
setMethod('^', signature=c('numeric', CLASS_CLIST),
          function(e1, e2) XresolveLeft(e1, e2, "R.pow"))
setMethod('^', signature=c('logical', CLASS_CLIST),
          function(e1, e2) XresolveLeft(e1, e2, "R.pow"))
setMethod('^', signature=c(CLASS_CLIST, 'numeric'),
          function(e1, e2) XresolveRight(e1, e2, "R.pow"))
setMethod('^', signature=c(CLASS_CLIST, 'logical'),
          function(e1, e2) XresolveRight(e1, e2, "R.pow"))
setMethod('^', signature=c('RMmodel', 'character'),
          function(e1, e2) Xresolve(e1, e2, "R.pow")) 
setMethod('^', signature=c('character', 'RMmodel'),
          function(e1, e2) Xresolve(e1, e2, "R.pow")) 



## str, print and show

## str.RMmodel is basically copied from str.default, but where
## those slots of 'par.general' are excluded from being printed
## for which the value is 'RFdefault', i.e. for which there is no
## explicit value given
str.RMmodel <-
  function(object, max.level = NA, vec.len = strO$vec.len,
           digits.d = strO$digits.d,
           nchar.max = 128, give.attr = TRUE, give.head = TRUE,
           give.length = give.head, 
           width = getOption("width"), nest.lev = 0,
           indent.str = paste(rep.int(" ", 
             max(0, nest.lev + 1)), collapse = ".."),
           comp.str = "$ ", no.list = FALSE, envir = baseenv(),
           strict.width = strO$strict.width, 
           formatNum = strO$formatNum, list.len = 99, ...) 
{
  oDefs <- c("vec.len", "digits.d", "strict.width", "formatNum")
  strO <- getOption("str")
  if (!is.list(strO)) {
    warning("invalid options('str') -- using defaults instead")
    strO <- strOptions()
  }
  else {
    if (!all(names(strO) %in% oDefs)) 
      warning("invalid components in options('str'): ", 
              paste(setdiff(names(strO), oDefs), collapse = ", "))
    strO <- modifyList(strOptions(), strO)
  }
  strict.width <- match.arg(strict.width, choices = c("no", "cut", "wrap"))
  if (strict.width != "no") {
    ss <- capture.output(str(object, max.level = max.level, #
                             vec.len = vec.len, digits.d = digits.d,
                             nchar.max = nchar.max, 
                             give.attr = give.attr,
                             give.head = give.head,
                             give.length = give.length, 
                             width = width, nest.lev = nest.lev,
                             indent.str = indent.str, 
                             comp.str = comp.str,
                             no.list = no.list || is.data.frame(object),
                             envir = envir, strict.width = "no", ...))
    if (strict.width == "wrap") {
      nind <- nchar(indent.str) + 2
      ss <- strwrap(ss, width = width, exdent = nind)
    }
    if (any(iLong <- nchar(ss) > width)) 
      ss[iLong] <- sub(sprintf("^(.{1,%d}).*", width - 2), "\\1..", ss[iLong])
    cat(ss, sep = "\n")
    return(invisible())
  }
  oo <- options(digits = digits.d)
  on.exit(options(oo))
  le <- length(object)
  P0 <- function(...) paste(..., sep = "")
  `%w/o%` <- function(x, y) x[is.na(match(x, y))]
  nfS <- names(fStr <- formals())
  strSub <- function(obj, ...) {
    nf <- nfS %w/o% c("object", "give.length", "comp.str", 
                      "no.list", names(match.call())[-(1:2)], "...")
    aList <- as.list(fStr)[nf]
    aList[] <- lapply(nf, function(n) eval(as.name(n)))
    if ("par.general" %in% names(obj)){
      is.RFdefault <-
        unlist(lapply(obj$par.general,
                      FUN=function(x){
                        !is(x, CLASS_CLIST) && !is.na(x) && x==RM_DEFAULT
                      }))
      obj$par.general[is.RFdefault] <- NULL
      if (all(is.RFdefault)) obj$par.general <- list()
    }
    do.call(utils::str, c(list(object = obj), aList, list(...)), 
            quote = TRUE)
  }
  v.len <- vec.len
  std.attr <- "names"
  cl <- if ((S4 <- isS4(object))) 
    class(object)
  else oldClass(object)
  has.class <- S4 || !is.null(cl)
    mod <- ""
  char.like <- FALSE
  if (give.attr) 
    a <- attributes(object)
  if (is.null(object)) 
    cat(" NULL\n")
  else if (S4) {
    a <- sapply(methods::.slotNames(object), methods::slot, 
                object = object, simplify = FALSE)
    cat("Formal class", " '", paste(cl, collapse = "', '"), 
        "' [package \"", attr(cl, "package"), "\"] with ", 
        length(a), " slots\n", sep = "")
    strSub(a, comp.str = "@ ", no.list = TRUE,
           give.length = give.length, 
           indent.str = paste(indent.str, ".."),
           nest.lev = nest.lev + 1)
    return(invisible())
  }
}

summary.RMmodel <- function(object, max.level=5, ...)
  summary(PrepareModel2(object, ...), max.level=max.level)

summary.RM_model <- function(object, ...) {
  class(object) <- "summary.RMmodel"
  object
}

print.summary.RMmodel <- function(x, max.level=5, ...) {
  str(x, no.list=TRUE, max.level = max.level, give.attr=FALSE) #
  invisible(x)
}

print.RM_model <- function(x, max.level=5,...) {
  print.summary.RMmodel(summary.RM_model(x, max.level=max.level,...),
                        max.level=max.level)#
}
print.RMmodel <- function(x, max.level=5,...) {
  print.summary.RMmodel(summary.RMmodel(x, max.level=max.level, ...),
                        max.level=max.level)#
}


setMethod(f="show", signature=CLASS_CLIST,
          definition=function(object) print.RMmodel(object))#

	  
print.RMmodelgenerator <- function(x, ...) {
  cat("*** object of Class '", CLASS_RM, "' ***\n", sep="")
  str(args(x@.Data)) ## ok
  cat("  type : \t", paste(x@type, collapse=", "), "\n")
  cat("  domain : \t", paste(x@domain, collapse=", "), "\n")
  cat("  isotropy : \t", paste(x@isotropy, collapse=", "), "\n")
  cat("  monotoniciy :\t", paste(x@monotone, collapse=", "), "\n")
  cat("  multivariate:\t",
      if (x@vdim >= 0) x@vdim
      else if (x@vdim == PARAM_DEP) "parameter dependent"
      else if (x@vdim == PREVMODEL_DEP) "depends on calling model"
      else if (x@vdim == SUBMODEL_DEP) "submodel dependent"
      else "specification unclear -- please contact maintainer",
      "\n")
  cat("  max. dimen. :\t",
      if (x@maxdim >= 0) x@maxdim
      else if (x@maxdim == PARAM_DEP) "parameter dependent"
      else if (x@maxdim == PREVMODEL_DEP) "depends on calling model"
      else if (x@maxdim == SUBMODEL_DEP) "submodel dependent"
      else "specification unclear -- please contact maintainer",
      "\n")
  cat("  finite range:\t", x@finiterange, "\n")
  cat("  operator : \t", x@operator, "\n")
  cat("  simple fctn :\t", x@simpleArguments, "\n")
}

setMethod("show", signature=CLASS_RM,
	  definition=function(object) print.RMmodelgenerator(object))

rfConvertRMmodel2string <- function(model){
  if (!is(model, class2=CLASS_CLIST))
    stop("model must be of class '", CLASS_CLIST, "'")
  par <- c(model@par.model, model@par.general)
  idx.random <- unlist(lapply(par, FUN=isRMmodel))
  if (is.null(idx.random)){
    param.string <- ""
    param.random.string <- ""
  } else {
    idx.default <- par[!idx.random] == RM_DEFAULT
    param.string <- paste(names(par[!idx.random][!idx.default]),
                          par[!idx.random][!idx.default],
                          sep="=", collapse=", ")

    #if (sum(idx.random) > 0){
    string.vector <- lapply(par[idx.random], FUN=rfConvertRMmodel2string)
    param.random.string <- paste(names(par[idx.random]), string.vector,
                                 sep="=", collapse=", ")
  }
  
  if (length(model@submodels) > 0){
    string.vector <- lapply(model@submodels, FUN=rfConvertRMmodel2string)
    submodel.string <- paste(names(model@submodels), string.vector,
                             sep="=", collapse=", ")
    if (!(nchar(param.string)==0))
      submodel.string <- paste(submodel.string, ", ", sep="")
  }
  else submodel.string <- ""
  
  string <- paste(model@name, "(", submodel.string, param.random.string,
                  param.string, ")", sep="")
  return(string)
}


## Plot
preparePlotRMmodel <- function(x, xlim, ylim, n.points, dim, fct.type,
                               MARGIN, fixed.MARGIN, ...){
  types <- c("Cov", "Variogram", "Fctn")

##  types <- "Fctn"; print("remove")
  
  verballist <- paste("'", types, "'", sep="", collapse="")
  if (!missing(fct.type) && length(fct.type) > 0) {
    if (!(fct.type %in% types))
      stop("fct.type must be NULL or of the types ", verballist)
    types <- fct.type
  }

  all.fct.types <- character(length(x))
  all.vdim <- numeric(length(x))

  for (i in 1:length(x)) {
    fct.type <- types
    m <- list("", PrepareModel2(x[[i]]), ...)
    while (length(fct.type) > 0 &&
           { m[[1]] <- fct.type[1];
           !is.numeric(vdim <- try( InitModel(MODEL_AUX, m, dim),
				   silent=TRUE))
##	   Print(CCC)
	   })
      fct.type <- fct.type[-1]
    if (!is.numeric(vdim)) {
      ##    Print(vdim, model)
      ##   stop("'vdim' is not numeric")
      stop(attr(vdim, "condition")$message)
    }
    if (vdim[1] != vdim[2]) stop("only simple models can be plotted")
    all.vdim[i] <- vdim[1]
    all.fct.types[i] <- fct.type[1]
  }

  if (!all(all.vdim == all.vdim[1]))
    stop("models have different multivariability")
   

#  fctcall.vec <- c("Cov", "Variogram")
  #covinfo <- RFgetModelInfo(register=MODEL_AUX, level=3, spConform=TRUE)
  ### Print((covinfo))
#  Print(covinfo$domown, TRANS_INV , covinfo$type, isPosDef(covinfo$type),
 #       is.null(fct.type), fct.type,"Variogram") ; xxx

  if (is.null(xlim)) {
    xlim <- if (dim > 1 || all.vdim[1] > 1) c(-1, 1) * 1.75 else c(0, 1.75)    
  }
  if (is.null(ylim) && dim > 1) ylim <- xlim
 
  
  distance <- seq(xlim[1], xlim[2], length=n.points)
  
  if (prod(xlim) <= 0)
    distance <- sort(c(if (!any(distance==0)) 0, 1e-5, distance))

  if (dim > 1) {
    distanceY <- seq(ylim[1], ylim[2], length=n.points)
    if (prod(ylim) < 0 & !(any(distanceY==0)))
      distanceY <- sort(c(0, distanceY))
  }

#  Print(covinfo)

#  if(covinfo$domown != TRANS_INV)
#    stop("method only implemented for models that are not kernels")

 
  if (all(all.fct.types[1] == all.fct.types)) {
    switch(all.fct.types[1],
         "Cov" = {
           main <-"Covariance function"
           ylab <- "C(distance)"
         },
         "Variogram" = {
           main <- "Variogram"
           ylab <- expression(gamma(distance))
         },
         "Fctn" = {
           main<- ""
           ylab <- "f(distance)"
         },
         stop("method only implemented for ", verballist)
       )
  } else {
    main <- ""
    ylab <- "f(distance)"
  }

  if (length(x) == 1) {
    for.what <- rfConvertRMmodel2string(x[[1]])
    if (nchar(for.what) > 20) for.what <- strsplit(for.what,"\\(")[[1]][1]
  } else {
    for.what <- "various models"
  }
   
 # main <- paste(main, " plot for ", for.what, "\n", sep="")
  main <- paste("plot for ", for.what, "\n", sep="")
  
  if (dim >= 3)
    main <- paste(sep="", main, ";  component", if (dim>3) "s", " ",
                  paste((1:dim)[-MARGIN], collapse=", "), 
                  " fixed to the value", if (dim>3) "s", " ",
                  paste(format(fixed.MARGIN, digits=4), collapse=", "))
  
  return(list(main=main, fctcall=all.fct.types, ylab=ylab, distance=distance,
              distanceY = if (dim > 1) distanceY, xlim=xlim, ylim=ylim))
}


singleplot <- function(cov, dim, li, plotmethod, dots, dotnames) {    
  if (dim==1) {
    D <- li$distance             
    iszero <- D == 0
    if (plotmethod == "matplot")  {
        plotpoint <- any(iszero) && diff(range(dots$ylim)) * 1e-3 <
          diff(range(cov)) - diff(range(cov[!iszero]))
        if (plotpoint) D[iszero] <- NA
    } else plotpoint <- FALSE
    liXY <-
      if (plotmethod=="plot.xy") list(xy = xy.coords(x=D, y=cov)) ## auch D ?
      else list(x=D, y=cov)
    
    ##       Print(plotmethod, args=c(dots, liXY))
    
    do.call(plotmethod, args=c(dots, liXY))   
    if (plotpoint) {
      for (i in 1:ncol(cov))
        points(0, cov[iszero, i], pch=19 + i, col = dots$col[i])
    }
    
  } else {
    if (!("zlim" %in% dotnames)) dots$zlim <- range(unlist(cov), finite=TRUE)
    addgiven <- "add" %in% dotnames
    local.dots <- dots
    local.dots$col <- NULL
    local.dots$lty <- NULL
    is.contour <- is.character(plotmethod) && plotmethod == "contour"
    if (!is.contour) {
      if (ncol(cov) > 1)
        stop("several models can be plotted at once only with 'contour'")
      col <- default.image.par(dots$zlim, NULL)$data$default.col
    }
    for (i in 1:ncol(cov)) {
      if (!addgiven && is.contour) local.dots$add <- i > 1
      do.call(plotmethod,
              args=c(list(x=li$distance, y=li$distanceY,
                  col=if (is.contour) dots$col[i] else col,
                  lty = dots$lty[i],
                  z=matrix(cov[, i], nrow=length(li$distance))), local.dots))
    }
  }
}


RFplotModel <- function(x, y, dim=1,
                        n.points=
                        if (dim==1 || is.contour) 200 else 100,
                        fct.type=NULL,
                        MARGIN, fixed.MARGIN,
                        maxchar=15, ...,
                        plotmethod=if (dim==1) "matplot" else "contour") {
  
  is.contour <- is.character(plotmethod) && plotmethod == "contour"
  RFopt <- RFoptions()
 
  if (ex.red <- RFopt$internal$examples_reduced)
    n.points <- as.integer(min(n.points, ex.red - 2))
  
  stopifnot(length(dim)==1)
  if (!(dim %in% 1:10))
    stop("only 'dim==1', 'dim==2' and 'dim==3' are allowed")
  if (dim==3)
    if (missing(MARGIN) || missing(fixed.MARGIN))
      stop("'MARGIN' and 'fixed.MARGIN' must be given if dim >=3")
  if ((!missing(MARGIN)) || (!missing(fixed.MARGIN))) {
    stopifnot((!missing(MARGIN)) && (!missing(fixed.MARGIN)))
    if (dim < 3)
      stop("'MARGIN' and 'fixed.MARGIN' should only be given for dim>=3")
    stopifnot(is.numeric(MARGIN) && length(MARGIN)==2)
    stopifnot(is.numeric(fixed.MARGIN) && (length(fixed.MARGIN) == dim - 2))
  }
  dots <- list(...)
  dotnames <- names(dots)
  models <- substr(dotnames, 1, 5) == "model"
  x <- c(list(x), dots[models])
  mnames <- c("", substr(dotnames[models], 6, 100))
  idx <- substr(mnames, 1, 1) == "."
  mnames[idx] <- substr(mnames[idx], 2, 1000)
  for (i in which(!idx)) {
    mnames[i] <- rfConvertRMmodel2string(x[[i]])
    nmn <- nchar(mnames[i]) - 1 
    if (substr(mnames[i], nmn, nmn) == "(") {
      mnames[i] <- substr(mnames[i], 1, nmn - 1)
    }
    msplit <- strsplit(mnames[i], "RM")[[1]]
    if (length(msplit) == 2 && msplit[1]=="")  mnames[i] <- msplit[2]
  }
  mnames <- substr(mnames, 1, maxchar)

  dots <- mergeWithGlobal(dots[!models])
  dotnames <- names(dots)  

  if (!("type" %in% dotnames)) dots$type <- "l"
  li <- preparePlotRMmodel(x=x, xlim=dots$xlim, ylim=dots$ylim,
                           n.points=n.points, dim=dim,
                           fct.type=fct.type,
                           MARGIN=MARGIN, fixed.MARGIN=fixed.MARGIN)

  dots$xlim <- li$xlim
  if (!is.null(li$ylim)) dots$ylim <- li$ylim
  if (!("main" %in% dotnames)) dots$main <- li$main
  if (!("cex" %in% dotnames)) dots$cex <- 1
  if (!("cex.main" %in% dotnames)) dots$cex.main <- 1.3 * dots$cex
  if (!("cex.axis" %in% dotnames)) dots$cex.axis <- 1.0 * dots$cex
  if (!("cex.lab" %in% dotnames)) dots$cex.lab <- 1.0 * dots$cex
  if (!("col" %in% dotnames)) dots$col <- rep(1:7, length.out=length(x))

  cov <- list()
  if (dim==1) {
    for (i in 1:length(x)) 
      cov[[i]] <- rfeval(x=li$distance, model=x[[i]], fctcall=li$fctcall[i])

    #print(cov, li$distance)
    lab <- xylabs("distance", li$ylab)
    ## strokorb monotone ist z.B. nicht ueberall finite:
    if (!("ylim" %in% dotnames)) dots$ylim <- range(0, unlist(cov), finite=TRUE)
    if (!("xlab" %in% dotnames)) dots$xlab <- lab$x
    if (!("ylab" %in% dotnames)) dots$ylab <- lab$y
  } else {  
    lab <- xylabs("", "")    
    if (!("xlab" %in% dotnames)) dots$xlab <- lab$x
    if (!("ylab" %in% dotnames)) dots$ylab <- lab$y
    dots$type <- NULL
 
    if (dim==2) {
      di <- as.matrix(expand.grid(li$distance, li$distanceY))

       for (i in 1:length(x)) {
        cov[[i]] <- rfeval(x=di, model=x[[i]], fctcall=li$fctcall[i])      
     }
    } else if (dim>=3) {
      m1 <- expand.grid(li$distance, li$distanceY)
      m2 <- matrix(NA, ncol=dim, nrow=nrow(m1))
      m2[,MARGIN] <- as.matrix(m1)
      m2[,-MARGIN] <- rep(fixed.MARGIN, each=nrow(m1))
      for (i in 1:length(x))
        cov[[i]] <- rfeval(x=m2, model=x[[i]], fctcall=li$fctcall[i])
    } else stop("this error should never appear")
  }

  if ((is.null(dots$xlab) || dots$xlab=="") &&
      (is.null(dots$ylab) || dots$ylab=="")) {
    margins <- c(3, 3, if (dots$main=="") 0 else 2, 0) + 0.2
  } else {
    margins <- c(5, 5, if (dots$main=="") 0 else 2, 0) + 0.2
  }

  dimcov <-  dim(cov[[1]])
  graphics <- RFopt$graphics  
  if (plotmethod != "plot.xy") { ### points and lines !!
     if (is.null(dimcov)) { ## i.e. vdim == 1
      ArrangeDevice(graphics, c(1,1))
    } else {
      figs <- dimcov[2:3]
      ArrangeDevice(graphics, figs)
    }
  }
    
  
  scr <- NULL
  if (is.null(dimcov)) { ## i.e. vdim == 1
    cov <- sapply(cov, function(x) x)
    if (!("lty" %in% dotnames)) dots$lty <- 1:5
    singleplot(cov=cov, dim=dim, li=li, plotmethod, dots=dots,
               dotnames=dotnames)
 #   Print(mnames, dots$col, dots$lty)
    
    if (length(x) > 1) legend(x="topright", legend=mnames,
                             col=dots$col, lty=dots$lty)
  } else { ## multivariate 
    scr <- matrix(split.screen(figs=figs), ncol=dimcov[2], byrow=TRUE)
    par(oma=margins)
    title(main=dots$main, xlab=dots$xlab, ylab=dots$ylab,
          outer=TRUE, cex.main=dots$cex.main)
    dots.axis = dots[names(dots) != "type"]
    dots.axis$col = dots.axis$col.axis
    if (!("axes" %in% dotnames)) dots$axes <- FALSE
    
    for (i in 1:ncol(scr)) {
      for (j in 1:nrow(scr)) {
        dots$main <-
          eval(parse(text=paste("expression(C[", i, j,"])", sep="")))
        screen(scr[i,j])
        par(mar=c(0,0,2,1))

        singleplot(cov = sapply(cov, function(x) x[, i,j]),
                   dim = dim, li=li, plotmethod=plotmethod, dots=dots,
                   dotnames = dotnames)

        box()
        if (i==1) do.call(graphics::axis,
              args=c(dots.axis, list(side=1, outer = TRUE, line=1)))
        if (j==1) do.call(graphics::axis,
              args=c(dots.axis, list(side=2, outer = TRUE, line=1)))
      }
    }
  }

  if (graphics$split_screen && graphics$close_screen) {
    close.screen(scr)
    scr <- NULL
  }

  return(scr)
}


points.RMmodel <- function(x, ..., type="p")  
  RFplotModel(x, ...,type=type,  plotmethod="plot.xy")
lines.RMmodel <- function(x, ..., type="l")
  RFplotModel(x, ...,type=type,  plotmethod="plot.xy")



# @FUNCTION-STARP**********************************************************
# @NAME		list2RMmodel
# @REQUIRE	a list
# @AUTHOR	A Malinowski <malinows@math.uni-goettingen.de>
# @DATE		26.09.2011
# @FUNCTION-END************************************************************


list2RMmodel <- function(x) { 
  if (!is.list(x)) return(x)

  name <- x[[1]]
  if (!is.character(name)) return(x)
  if (name == RM_DECLARE) return(NULL)
  
  len <- length(x)
  
  if (name %in% DOLLAR) return(list2RMmodel(c(x[[len]], x[-c(1, len)])))
#  if (name == RM_DISTR[1]){
#    x <- x[-1]
#    m <- do.call(RM_DISTR[1], args=list())
#    m@par.model <- x
#    return(m)
#  }
  
  if (name==SYMBOL_PLUS) name <- RM_PLUS[1] else
  if (name==SYMBOL_MULT) name <- RM_MULT[1] else
  if (!(name %in% list2RMmodel_Names)) {
    if (!(name %in% list2RMmodel_oldNames))
      stop(paste("'", name, "' is not the name of a valid model", sep=""))
  }
  ## hier noch mehr Ausnahmen?


  if (len==1) return(eval(parse(text=paste(name, "()", sep=""))))
  else {
    x <- lapply(x, FUN=list2RMmodel)
    if (length(idx <- which("anisoT" == names(x))) == 1){
      names(x)[idx] <- "Aniso"
      x[[idx]] <- t(x[[idx]])
    }
    return(do.call(name, args=x[-1]))
  }
}

setMethod(f="plot", signature(x=CLASS_CLIST, y="missing"),
          function(x, y, ...) RFplotModel(x, ...))
setMethod(f="lines", signature(x=CLASS_CLIST),
          function(x, ..., type="l")
          RFplotModel(x, ..., type=type, plotmethod="plot.xy"))
setMethod(f="points", signature(x=CLASS_CLIST),
          function(x, ..., type="p")
          RFplotModel(x, ..., type=type, plotmethod="plot.xy"))
setMethod(f="persp", signature(x=CLASS_CLIST),
          function(x, ..., dim=2, zlab="")
          RFplotModel(x,...,dim=dim,zlab=zlab,plotmethod="persp"))
setMethod(f="image", signature(x=CLASS_CLIST),
          function(x, ..., dim=2) RFplotModel(x,...,dim=dim,plotmethod="image"))
