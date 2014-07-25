# accessing 'RMmodel' and RMmodelgenerator via '['-operator
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
    stop(paste(i, "is not a valid slot specification"))
  else
    slot(x, i) <- value
  validObject(x)
  return(x)
}

setMethod(f="[", signature = 'RMmodel', def = accessSlotsByName)
setMethod(f="[", signature = 'RMmodelFit', def = accessSlotsByName)
setMethod(f="[", signature = 'RMmodelgenerator', def = accessSlotsByName)

setReplaceMethod(f="[", signature = 'RMmodel',
                 def = accessReplaceSlotsByName)
setReplaceMethod(f="[", signature = 'RMmodelFit',
                 def = accessReplaceSlotsByName)
setReplaceMethod(f="[", signature = 'RMmodelgenerator',
                 def = accessReplaceSlotsByName)





## summing up RMmodels
setMethod('+', signature=c('RMmodel', 'RMmodel'),
          function(e1, e2) {
            d <- list()
            if (e1@name==ZF_PLUS[1]){
              len.e1 <- length(e1@submodels)              
              for (i in 1:len.e1)  d[[i]] <- e1@submodels[[i]]
            } else {
              len.e1 <- 1
              d[[1]] <- e1
            }
            d[[len.e1 + 1]] <- e2            
            model <- do.call(ZF_PLUS[1], d)
            return(model)
          })


## summing up RMmodels
setMethod('*', signature=c('RMmodel', 'RMmodel'),
          function(e1, e2) {
            d <- list()
            if (e1@name==ZF_MULT[1]){
              len.e1 <- length(e1@submodels)              
              for (i in 1:len.e1)  d[[i]] <-  e1@submodels[[i]]
            } else {
              len.e1 <- 1
              d[[1]] <- e1
            }
            d[[len.e1 + 1]] <- e2            
            model <- do.call(ZF_MULT[1], d)
            return(model)
          })



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
                        !is(x, ZF_MODEL) && !is.na(x) && x==ZF_DEFAULT_STRING
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

summary.RMmodel <- function(object, ...) summary(PrepareModel2(object, ...))

setMethod(f="show", signature='RMmodel',
          definition=function(object) print.RMmodel(object, max.level=8))#

summary.RM_model <- function(object, ...) {
  class(object) <- "summary.RMmodel"
  object
}

print.summary.RMmodel <- function(x, max.level=3, ...) {
  str(x, no.list=TRUE, max.level = max.level) #
  invisible(x)
}

print.RM_model <- function(x, ...) {
  print.summary.RMmodel(summary.RM_model(x, ...))#
}
print.RMmodel <- function(x, ...) {
  print.summary.RMmodel(summary.RMmodel(x, ...))#
}

setMethod(f="show", signature='RMmodel',
          definition=function(object) print.RMmodel(object))#







print.RMmodelgenerator <- function(x, ...) {
  cat("*** object of Class '", ZF_MODEL_FACTORY,
      "' ***\n  ** str(object):\n  ", sep="") #
  str(x)#
  cat("*** end of '", ZF_MODEL_FACTORY, "' ***\n", sep="")#
}
  
setMethod(f="show", signature='RMmodelgenerator',
          definition=function(object) print.RMmodelgenerator(object))#



## ## buildCovList
## setGeneric(name = "rfBuildCovList", 
##            def = function(model) standardGeneric("rfBuildCovList"))
## setMethod("rfBuildCovList", signature=c('RMmodel'),
##           function(model) buildCovList(model))
##           ## see convert_new.R for fct buildCovList



rfConvertRMmodel2string <- function(model){
  if (!is(model, class2=ZF_MODEL))
    stop("model must be of class 'RMmodel'")
  par <- c(model@par.model, model@par.general)
  idx.random <- unlist(lapply(par, FUN=isModel))
  if (is.null(idx.random)){
    param.string <- ""
    param.random.string <- ""
  } else {
    idx.default <- par[!idx.random] == ZF_DEFAULT_STRING
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

setMethod(f="plot", signature(x='RMmodel', y="missing"),
	  definition=function(x, y, dim=1, n.points=200, fct.type=NULL,
            MARGIN, fixed.MARGIN, ...)
          plotRMmodel(x, dim=dim, n.points=n.points, fct.type=fct.type,
                      MARGIN=MARGIN, fixed.MARGIN=fixed.MARGIN, ...))
#setMethod(f="points", signature(x='RMmodel', y="missing"),
#	  definition=function(x, y, ...) pointsRMmodel(x, ...))
#setMethod(f="lines", signature(x='RMmodel', y="missing"),
#	  definition=function(x, y, ...) linesRMmodel(x, ...))

preparePlotRMmodel <- function(x, xlim, ylim, n.points, dim, fct.type,
                               MARGIN, fixed.MARGIN){
  types <- c("Cov", "Variogram", "Fctn")
  verballist <- paste("'", types, "'", sep="", collapse="")
  if (missing(fct.type) || length(fct.type) == 0) {
    fct.type <- types
  } else {
    if (!(fct.type %in% types))
      stop("fct.type must be NULL or of the types ", verballist)
  }

  model <- list("", PrepareModel2(x))
  while (length(fct.type) > 0 &&
         { model[[1]] <- fct.type[1];
           !is.numeric(vdim <- try( InitModel(MODEL.USER, model, dim),
                           silent=TRUE))})
    fct.type <- fct.type[-1]

  fct.type <- fct.type[1]
  
  if (!is.numeric(vdim)) {
#    Print(vdim, model)
 #   stop("'vdim' is not numerical")
    stop(attr(vdim, "condition")$message)
  }

#  fctcall.vec <- c("Cov", "Variogram")
  covinfo <- RFgetModelInfo(register=MODEL.USER, level=2, spConform=TRUE)
  ### Print((covinfo))
#  Print(covinfo$domown, TRANS_INV , covinfo$type, isPosDef(covinfo$type),
 #       is.null(fct.type), fct.type,"Variogram") ; xxx
  
  distance <- seq(xlim[1], xlim[2], length=n.points)
  
  if (xlim[1] * xlim[2] < 0 & !(any(distance==0)))
    distance <- sort(c(0, distance))

  if (dim > 1) {
    distanceY <- seq(ylim[1], ylim[2], length=n.points)
    if (ylim[1]*ylim[2] < 0 & !(any(distanceY==0)))
      distanceY <- sort(c(0, distanceY))
  }

#  Print(covinfo)

#  if(covinfo$domown != TRANS_INV)
#    stop("method only implemented for models that are not kernels")

  switch(fct.type,
         "Cov" = {
           main<- paste("Covariance fct plot for\n", rfConvertRMmodel2string(x))
           ylab <- "C(distance)"
         },
         "Variogram" = {
           main <- paste("Variogram plot for\n", rfConvertRMmodel2string(x))
           ylab <- expression(gamma(distance))
         },
         "Fctn" = {
           main<- paste("plot for shape function\n", rfConvertRMmodel2string(x))
           ylab <- "f(distance)"
         },
         stop("method only implemented for ", verballist)
       )
  if (dim >= 3)
    main <- paste(sep="", main, ";  component", if (dim>3) "s", " ",
                  paste((1:dim)[-MARGIN], collapse=", "), 
                  " fixed to the value", if (dim>3) "s", " ",
                  paste(format(fixed.MARGIN, digits=4), collapse=", "))

  .C("DeleteKey", MODEL.USER)

  
  return(list(main=main, fctcall=fct.type, ylab=ylab, distance=distance,
              distanceY = if (dim > 1) distanceY))
}

plotRMmodel <- function(x, y, dim, n.points, fct.type, MARGIN, fixed.MARGIN,
                        ...) {
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
  
  if (any(par()$mfcol != c(1,1))) par(mfcol=c(1,1))
  plotfctcall <- if ("plotfctcall" %in% dotnames)
    dots$plotfctcall else "plot.default"
  dots$plotfctcall <- NULL

  if (!("xlim" %in% dotnames)) {
    dots$xlim <- c(-1, 1)
  }
  if (!("ylim" %in% dotnames) && dim > 1) dots$ylim <- dots$xlim
  if (!("type" %in% dotnames)) dots$type <- "l"

  
  li <- preparePlotRMmodel(x=x, xlim=dots$xlim, ylim=dots$ylim,
                           n.points=n.points, dim=dim,
                           fct.type=fct.type,
                           MARGIN=MARGIN, fixed.MARGIN=fixed.MARGIN)
  
  if (!("main" %in% dotnames)) dots$main <- li$main
  if (!("cex" %in% dotnames)) dots$cex <- 1
  if (!("cex.main" %in% dotnames)) dots$cex.main <- 0.8 * dots$cex
  if (!("cex.axis" %in% dotnames)) dots$cex.axis <- 0.8 * dots$cex
  if (!("col" %in% dotnames)) 
    dots$col <- default.image.par(dots$ylim, NULL, FALSE)$data$default.col[1]
 
  
  if (dim==1) {
    cov <- RFcov(x=li$distance, model=x, fctcall=li$fctcall)

    #print(cov, li$distance)
    lab <- xylabs("distance", li$ylab)
    ## strokorb monotone ist z.B. nicht ueberall finite:
    if (!("ylim" %in% dotnames)) dots$ylim <- range(0, cov[is.finite(cov)])
    if (!("xlab" %in% dotnames)) dots$xlab <- lab$x
    if (!("ylab" %in% dotnames)) dots$ylab <- lab$y
  } else {  
    lab <- xylabs("", "")    
    if (!("xlab" %in% dotnames)) dots$xlab <- lab$x
    if (!("ylab" %in% dotnames)) dots$ylab <- lab$y
    dots$type <- NULL
 
    if (dim==2)
      cov <- RFcov(x=as.matrix(expand.grid(li$distance, li$distanceY)),
                   model=x, fctcall=li$fctcall)

    if (dim>=3) {
      m1 <- expand.grid(li$distance, li$distanceY)
      m2 <- matrix(NA, ncol=dim, nrow=nrow(m1))
      m2[,MARGIN] <- as.matrix(m1)
      m2[,-MARGIN] <- rep(fixed.MARGIN, each=nrow(m1))
      cov <- RFcov(x=m2, model=x, fctcall=li$fctcall)      
    }
  }

  if ((is.null(dots$xlab) || dots$xlab=="") &&
      (is.null(dots$ylab) || dots$ylab=="")) {
    margins <- c(3, 3, if (dots$main=="") 0 else 2, 0) + 0.2
  } else {
    margins <- c(5, 5, if (dots$main=="") 0 else 2, 0) + 0.2
  }

  if (is.null(dim(cov))) { ## i.e. vdim == 1
    if (dim==1){
       D <- li$distance             
       if (plotfctcall == "plot.default")  {
         iszero <- D == 0
         D[iszero] <- NA
         plotpoint <- any(iszero)
       }  
       liXY <-
        if (plotfctcall=="plot.xy") list(xy = xy.coords(x=D, y=cov)) ## auch D ?
        else list(x=D, y=cov)

  #       Print(plotfctcall, args=c(dots, liXY))
      
       do.call(plotfctcall, args=c(dots, liXY))
      if (plotpoint) points(0, cov[iszero], pch=20, col = dots$col)
      
    } else {
      do.call(graphics::contour, args=c(list(x=li$distance, y=li$distanceY,
                         z=matrix(cov, nrow=length(li$distance))), dots))
    }
  } else { ## multivariate 
    graphics <- RFoptions()$graphics
    ArrangeDevice(graphics, dim(cov)[2:3])
    scr <- matrix(split.screen(figs=dim(cov)[2:3]),
                  ncol=dim(cov)[2], byrow=TRUE)
    par(oma=margins)
    title(main=dots$main, xlab=dots$xlab, ylab=dots$ylab,
          outer=TRUE, cex.main=dots$cex.main)
    dots.axis = dots[names(dots) != "type"]
    dots.axis$col = dots.axis$col.axis
    for (i in 1:ncol(scr)) {
      for (j in 1:nrow(scr)) {
        dots$main <-
          eval(parse(text=paste("expression(C[", i, j,"])", sep="")))
        screen(scr[i,j])
        par(mar=c(0,0,2,1))
        covij <- cov[,i,j]
        if (dim==1) {
          liXY <-
            if (plotfctcall=="plot.xy")
              list(xy = xy.coords(x=li$distance, y=covij))
            else
              list(x=li$distance, y=covij)
          do.call(plotfctcall, args=c(dots, liXY, axes=FALSE))
        } else {
          do.call(graphics::contour, args=c(dots, list(axes=FALSE,
                             x=li$distance, y=li$distanceY,
                             z=matrix(covij, nrow=length(li$distance)))))
        }
        box()
        if (i==1) do.call(graphics::axis,
              args=c(dots.axis, list(side=1, outer = TRUE, line=1)))
        if (j==1) do.call(graphics::axis,
              args=c(dots.axis, list(side=2, outer = TRUE, line=1)))
      }
    }
    if (graphics$always_close_screen)
      close.screen(all.screens=TRUE)
  } 
}


points.RMmodel <- function(x, y, n.points = 200, fct.type = NULL, ...) {
  dots <- list(...)
  if (!("type" %in% names(dots))) dots$type <- "p"
  do.call(plotRMmodel, args=c(dots, list(
                         x=x, dim=1, n.points=n.points,
                         fct.type=fct.type, plotfctcall="plot.xy"))
          )
  ##  li <- preparePlotRMmodel(x=x, xlim=xlim, n.points=n.points, dim=dim,
  ##                           fct.type=fct.type)
  ##  cov <- RFcov(x=li$distance, model=x, fctcall=li$fctcall)
  ##  plot.xy(xy.coords(x=li$distance, y=cov), xlim=xlim, type="p",
  ##          pch=pch, lty=lty, col=col, bg=bg, cex=cex, lwd=lwd, ...)
}

lines.RMmodel <- function(x, y, n.points = 200, fct.type = NULL, ...) {
  dots <- list(...)
  if (!("type" %in% names(dots))) dots$type <- "l"
  do.call(plotRMmodel, args=c(dots, list(
                         x=x, dim=1, n.points=n.points,
                         fct.type=fct.type, plotfctcall="plot.xy"))
          )
}


# @FUNCTION-STARP**********************************************************
# @NAME		list2RMmodel
# @REQUIRE	a list
# @AUTHOR	A Malinowski <malinows@math.uni-goettingen.de>
# @DATE		26.09.2011
# @FUNCTION-END************************************************************


list2RMmodel <- function(x, skip.prefix.check=FALSE) { 
  #Print(x)

  if (!is.list(x))
    return(x)

  name <- x[[1]]
  if (!is.character(name))
    return(x)
  
  len <- length(x)
  
  if (name %in% DOLLAR)
    return(list2RMmodel(c(x[[len]], x[-c(1, len)])))

  if (name == ZF_DISTR[1]){
    x <- x[-1]
    m <- do.call(ZF_DISTR[1], args=list())
    m@par.model <- x
    return(m)
  }
  
  if (name==ZF_SYMBOLS_PLUS) name <- ZF_PLUS[1] else
  if (name==ZF_SYMBOLS_MULT) name <- ZF_MULT[1] else
  if (name %in% ZF_MIXED) name <- ZF_INTERNALMIXED


  ## add prefix 'RM' if necessary
  if (!skip.prefix.check){
    if (!(name %in% c(RFgetModelNames(), ZF_INTERNALMIXED, ZF_TREND))) {
      name2 <- paste(ZF_MODEL_PREFIX, name, sep="")
      if (!(name2 %in% c(RFgetModelNames(), ZF_INTERNALMIXED, ZF_TREND)))
        stop(paste("'", name, "' is not the name of a valid cov model", sep=""))
      else
        name <- name2
    }
  }

  
  if (len==1)
    return(eval(parse(text=paste(name, "()", sep=""))))
  else {
    x <- lapply(x, FUN=list2RMmodel, skip.prefix.check=skip.prefix.check)
    #argnames <- names(x[-1])
    if (length(idx <- which("anisoT" == names(x))) == 1){
      #argnames[idx] <- "Aniso"
      names(x)[idx] <- "Aniso"
      x[[idx]] <- t(x[[idx]])
    }
    #if (!is.null(argnames)) {
    #  isEmpty <- argnames=="" 
    #  argnames[!isEmpty] <- paste(argnames[!isEmpty], "=")
    #}
    return(do.call(name, args=x[-1]))
    
    #return(eval(parse(text =
    #                  paste(name, "(",
    #                        paste(argnames, " x[[", 2:len, "]]", 
    #                              collapse=","),
    #                        ")", sep=""))))
  }
}
