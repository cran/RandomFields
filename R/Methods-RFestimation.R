## Methods for classes 'RFempVario' and 'RFfit'  #######################

## accessing slots
accessByNameOrNumber <- function(x, i, j, drop=FALSE) {
  stopifnot(length(i)==1)
  if (is.numeric(i))  i <- slotNames(x)[i]
  return(accessSlotsByName(x=x, i=i, j=j, drop=drop))
}

setMethod("[", signature = "RFfit", def=accessByNameOrNumber)



anova.RFfit <- function(object, ...)  RFratiotest(nullmodel=object, ...)
anova.RF_fit <- function(object, ...) RFratiotest(nullmodel=object, ...)
anova.RMmodelFit <- function(object, ...) RFratiotest(nullmodel=object, ...)
anova.RM_modelFit <- function(object, ...) RFratiotest(nullmodel=object, ...)

setMethod(f="anova", signature='RMmodelFit',
          definition=function(object, ...) anova.RFfit(object, ...))#

boundary_values <- function(variab) {
  upper.bound <- variab[4, , drop=FALSE]
  lower.bound <- variab[3, , drop=FALSE]
  # sd <- variab[2, ]
  variab <- variab[1, , drop=FALSE]
  lidx <- variab < lower.bound + 1e-8
  uidx <- variab > upper.bound - 1e-8
  nl <- sum(lidx, na.rm=TRUE)
  nu <- sum(uidx, na.rm=TRUE)
  if (nl + nu > 0) {
    lidx[is.na(lidx)] <- FALSE
    uidx[is.na(uidx)] <- FALSE
    txt <-
      paste(sep="", "Note that the (internal) fitted variable",
            if (nl > 0)
            paste(if (nl > 1) "s", " ",
                  paste("'", names(variab[lidx]), "'",
                        sep="", collapse=", "),
                  " ", if (nl == 1)  "is" else "are",                  
                  " close to or on the effective lower boundary", sep=""),
            if (nl > 0 && nu > 0) " and the variable",
            if (nu > 0)
            paste(if (nu > 1) "s", " ",
                  paste("'", names(variab[uidx]), "'",
                        sep="", collapse=", "),
                  " ", if (nu == 1) "is" else "are",
                  " close to or on the effective upper boundary"),
            ".\nHence the gradient of the likelihood function might not be zero and none of the\nreported 'sd' values might be reliable.")
  } else txt <- NULL
  return(txt)
}


summary.RMmodelFit <- function(object, ..., isna.param) {
  model <- PrepareModel2(object, ...)
  residuals <- object@residuals
  l <- list(model=model,
            loglikelihood=object@likelihood,
            AIC = object@AIC,
            AICc= object@AICc,
            BIC = object@BIC,
            residuals=if (length(residuals) == 1) residuals[[1]] else residuals)
  if (missing(isna.param)) isna.param <- any(is.na(object@param)) 
  l$boundary <- boundary_values(object@variab)
  if (!any(is.na(object@param[1, ])) )
    l$param <- cbind(object@param, object@covariab)
  if (isna.param || !is.null(l$boundary)) {
#   print(object@variab)
#    print(object@covariab)
#   str(object)
    l$variab <- cbind(object@variab,
                      if (length(object@covariab) > 0)
                        rbind(object@covariab, NA, NA))
  }
  class(l) <- "summary.RMmodelFit"
  l
}

setMethod(f="summary", signature='RMmodelFit',
          definition=function(object) summary.RMmodelFit(object))#


summary.RM_modelFit <- function(object, ..., isna.param) {
  model <- object$model
  residuals <- object$residuals
  l <- list(model=model,
            loglikelihood=object$ml.value,
            AIC = object$AIC,
            AICc= object$AICc,
            BIC = object$BIC,
            residuals=if (length(residuals) == 1) residuals[[1]] else residuals)
  if (missing(isna.param)) isna.param <- any(is.na(object$param))
  l$boundary <- boundary_values(object$variab)
  if (!any(is.na(object$param[1,])))
      l$param <- cbind(object$param, object$covariab)
  if (isna.param || !is.null(l$boundary))
    l$variab <-
      cbind(object$variab,
            if (length(object$covariab) > 0) rbind(object$covariab, NA, NA))
  class(l) <- "summary.RMmodelFit"
  l
}



print.summary.RMmodelFit <- function(x, ...) {

  printVariab <- function(x) {
    cat("Internal variables:\n")
    if (is.null(x$boundary)) print(x$variab[1:2, , drop=FALSE], ..., na.print="-")#
    else print(x$variab, ..., na.print="-")#
    cat("\n")
    return(ncol(x$variab))
  }


  printParam <- function(param) {
    cat("User's variables:\n")
    print(param, ..., na.print="-")#
    return(ncol(param))
  }
  
  printRest <- function(...) {
    x <- unlist(list(...))
    stopifnot(length(x) == 3)
    names(x) <- c("#variab", "loglikelihood", "AIC")
    cat("\n")
    print(x) #
    cat("\n")
  }

 if (RFoptions()$general$detailed_output) str(x$model, no.list=TRUE) #
  cat("\n")
  np <- AIC <- ll <- nm <- NA  
  if (!is.null(x$submodels)) {
    cur_name <- ""
    len <- length(x$submodels)
    
    for (i in 1:len) {
      sm <- x$submodels[[i]]
      n <- sm$report
      nnxt <- if (i==len) "" else x$submodels[[i+1]]      
      if (n != cur_name) {
        if (i > 1) {         
          if (!is.null(sm$param)) printParam(cparam)
          printRest(np, ll, AIC) #
          if (!is.null(sm$boundary)) cat(sm$boundary, "\n\n")
        }

        if (nnxt != n && length(sm$fixed) > 0) {
          
          nX <- paste(sep="", n, " (",
                      paste(c(if (length(sm$fixed$zero) > 0)
                              paste(colnames(x$param)[sm$fixed$zero], "= 0"),
                              if (length(sm$fixed$one) > 0)
                              paste(colnames(x$param)[sm$fixed$one], "= 1")),
                            sep=", "),
                      ")")
        } else nX <- n
        cat(if (!is.na(nm)) cat("\n"), nX, "\n",
            paste(rep("=", min(80, nchar(nX))), collapse=""),
            "\n", sep="")
        np <- 0
        AIC <- 0
        ll <- 0
        cparam <- NULL
        nm <- 1
      }
      if (!is.null(sm$variab)) {
        if (nm > 1 || (i<len && n==nnxt)) cat("model", nm, ", ")
        printVariab(sm)
      }
      if (!is.null(sm$param)) {
        param <- x$param * NA
#        Print(i, sm$p.proj, param, sm, sm$param)
        
        param[, sm$p.proj] <- sm$param
        fixed <- sm$fixed
        if (length(fixed) > 0) {
          param[1, fixed$zero] <- 0
          param[1, fixed$one] <- 1
        }
      
      #  if (!is.null(cparam)) cparam <- rbind(cparam, NA)
        cparam <- rbind(cparam, param)
  
 #       Print(param, cparam)
        
      }
      np <- np + length(sm$p.proj)
      ll <- ll + sm$loglikelihood
      AIC <- AIC + sm$AIC
      nm <- nm + 1;
      cur_name <- n
   }

    if (!is.null(sm$param)) printParam(param)
    printRest(np, ll, AIC) #
    if (!is.null(sm$boundary)) cat(sm$boundary, "\n\n")
     
    cat("\nuser's model\n", paste(rep("=", 12), collapse=""), "\n", sep="")
  }
 
  if (!is.null(x$variab)) np <- printVariab(x)
  if (!is.null(x$param)) np <- printParam(x$param)
  
  printRest(np, x[c("loglikelihood", "AIC")])#
  if (!is.null(x$boundary)) cat(x$boundary, "\n\n")
  
  invisible(x)
}

print.RMmodelFit <- function(x, ...)
  print.summary.RMmodelFit(summary.RMmodelFit(x, ...))#
print.RM_modelFit <- function(x, ...)
  print.summary.RMmodelFit(summary.RM_modelFit(x, ...))#

setMethod(f="show", signature='RMmodelFit',
          definition=function(object) print.RMmodelFit(object))#

  
summary.RFfit <- function(object, ...,  method="ml", full=FALSE) {
  s <- summary.RMmodelFit(object[method])
  len <-  length(object@submodels)
#  str(object@submodels, 2)
  if (full && length(object@submodels) > 0) {
    submodels <- list()
    for (i in 1:len) {
      ## war summary.RM_modelFit
#      Print(object@submodels[[i]][[method]])
      
      submodels[[i]] <- summary(object@submodels[[i]][[method]],# 'summary'
                                isna.param=is.null(s$param))    # nicht
      submodels[[i]]$report <- object@submodels[[i]]$report     # spezifizieren!
      submodels[[i]]$p.proj <- object@submodels[[i]]$p.proj
      submodels[[i]]$fixed  <- object@submodels[[i]]$fixed
    }     
    s$submodels <- submodels
  }
  s
}

summary.RF_fit <- function(object, ..., method="ml", full=FALSE) {
  s <- summary.RM_modelFit(object[[method]])
  len <-  length(object$submodels)
  if (full && len > 0) {
    submodels <- list()
    for (i in 1:len) {      
      submodels[[i]] <- summary.RM_modelFit(object$submodels[[i]][[method]],
                                            isna.param=is.null(s$param))
      submodels[[i]]$report <- object$submodels[[i]]$report
      submodels[[i]]$p.proj <- object$submodels[[i]]$p.proj
      submodels[[i]]$fixed  <- object$submodels[[i]]$fixed
    }
    s$submodels <- submodels
   }
  s
}



print.RFfit <- function(x, ...,  method="ml", full=FALSE) {
  print.summary.RMmodelFit(summary.RFfit(x, ..., method=method, full=full))
}

setMethod(f="show", signature='RFfit',
          definition=function(object) print.RFfit(object))#

print.RF_fit <- function(x, ...,  method="ml", full=FALSE) {
   print.summary.RMmodelFit(summary.RF_fit(x, ..., method=method, full=full))
}



logLik.RF_fit <- function(object, REML = FALSE, ..., method="ml") {
  if (hasArg("REML")) stop("parameter 'REML' is not used. Use 'method' instead")
  ## according to geoR
  val <- object[[method]]$ml.value  
  attr(val, "df") <- object$number.of.parameters
  attr(val, "method") <- method
  class(val) <- "logLik"
  return(val)
}

logLik.RFfit <- function(object, REML = FALSE, ..., method="ml") {
  if (hasArg("REML")) stop("parameter 'REML' is not used. Use 'method' instead")
 ## according to geoR
  val <- object[method]@likelihood
  attr(val, "df") <- object@number.of.parameters
  attr(val, "method") <- method
  class(val) <- "logLik"
  return(val)
}



print.AICRFfit<- function(x, ..., digits=3) {
  ## nur deshalb 
  fstcol <- 3
  sndcol <- 55
  trdcol <- 4
  forthcol<-9
  leer <- formatC("", width=fstcol)
  size <- max(abs(x[[2]]))
  size <- if (size>0) ceiling(log(size) / log(10)) else 1
  cat(leer, formatC("model", flag="-", width=sndcol), " ",
      formatC(names(x)[1], width=trdcol),
      formatC(names(x)[2], width=forthcol), "\n", sep="")
  names <- attr(x, "row.names")
  for (i in 1:length(names)) {
    cat(formatC(i, width=fstcol, flag="-"))
    if (nchar(xx <- names[i]) <= sndcol)
      cat(formatC(xx, width=sndcol, flag="-"))
    else {
      yy <- strsplit(xx, " \\* ")[[1]]
      for (j in 1:length(yy)) {
        ncyy <- nchar(yy[j])
        if (ncyy <= sndcol && j==length(yy))
          cat(format(yy[j], width=sndcol, flag="-"))
        else {
          if (ncyy <= sndcol - 2) {
            cat(yy[j])
          } else {
            zz <- strsplit(yy[j], ", ")[[1]]
            ncyy <-  0
            lenzz <- length(zz)
            for (k in 1:lenzz) {
              len <- nchar(zz[k])
              if (k > 1 && len > sndcol - 1) {
                cat("\n", leer, zz[k], sep="")
                if (k < lenzz)
                  cat(formatC(",", flag="-",  width=pmax(1, sndcol-len)))
              } else {
                if (ncyy + len > sndcol - 1) {
                  cat("\n", leer, sep="")
                  ncyy <- len
                } else {
                  ncyy <- ncyy + len
                }
                cat(zz[k])
                if (k < lenzz) {
                  cat(", ")
                  ncyy <- ncyy + 2
                }
              }
            } # for k 1:lenzz
          } # split according to commata
          if (j < length(yy)) cat(" *\n", leer, sep="")
          else if (ncyy < sndcol) cat(formatC("", width=sndcol-ncyy))
        }
      } # for 1:products
    } ## not be written in a single line
    cat("",
        formatC(x[[1]][i], width=trdcol),
        formatC(x[[2]][i], format="f", width=size + digits + 1,
                    digits=digits),"\n")
  }
}



fullAIC <- function(x, method="ml", AIC="AIC") {
  ats <- approx_test_single(x, method=method)$result
  values <- c("name", "df", AIC)
  model2 <- paste("model2.", values, sep="")
  ats2 <- ats[ !is.na(ats[, model2[2]]), model2]
  colnames(ats2) <- values
  ats <- ats[, paste("model1.", values, sep="")]
  colnames(ats) <- values
  ats <- unique(rbind(ats, ats2))
  dimnames(ats) <- list(1:nrow(ats), dimnames(ats)[[2]])

  names  <- as.character(ats$name)
  ats <- ats[-1]
  attr(ats, "row.names") <- names  
  class(ats) <- "AICRFfit"
  ats
}

AIC.RFfit <- function(object, ..., k=2, method="ml", full=TRUE) {
  if (full) {
    fullAIC(object, method=method)
  } else {
    AIC <- object[method]@AIC
    names(AIC) <- "AIC"
    AIC
  }
}
AIC.RF_fit <- function(object, ..., k=2, method="ml", full=TRUE) {
  if (full) {
    fullAIC(object, method=method)
  } else {
    AIC <- object[[method]]$AIC
    names(AIC) <- "AIC"
    AIC
  }
}

AICc.RFfit <- function(object, ...,  method="ml", full=FALSE) {
  if (full) {
    stop("for 'AICc' the option 'full=TRUE' has not been programmed yet.")
    fullAIC(object, method=method)
  } else {
    AIC <- object[method]@AIC
    names(AIC) <- "AICc"
    AIC
  }
}
AICc.RF_fit <- function(object, ..., method="ml", full=TRUE) {
  if (full) {
    stop("for 'AICc' the option 'full=TRUE' has not been programmed yet.")
    fullAIC(object, method=method)
  } else {
    AIC <- object[[method]]$AIC
    names(AIC) <- "AICc"
    AIC
  }
}

BIC.RFfit <- function(object, ..., method="ml", full=TRUE) {
  if (full) {
    fullAIC(object, method=method, AIC="BIC")
  } else {
    BIC <- object[method]@BIC
    names(BIC) <- "BIC"
    BIC
  }
}

BIC.RF_fit <- function(object, ..., method="ml", full=TRUE) {
  if (full) {
    fullAIC(object, method=method, AIC="BIC")
  } else {
    BIC <- object[[method]]$BIC
    names(BIC) <- "BIC"
    BIC
  }
}


resid.RFfit <- function(object, ..., method="ml") {
  resid <- object[method]@residuals
  names(resid) <- "residuals"
  resid
}
resid.RF_fit <- function(object, ..., method="ml") {
  resid <- object[[method]]$residuals
  names(resid) <- "residuals"
  resid
}
residuals.RFfit <- function(object, ..., method="ml")
  resid.RFfit(object=object, method=method)
residuals.RF_fit <- function(object, ..., method="ml")
  resid.RF_fit(object=object, method=method)


summary.RFempVariog <- function(object, ...) {
  variogram <- object@emp.vario
  if (is.array(variogram)) {
    dims <- dim(variogram)
    dims <- dims[dims > 1]
    dim(variogram) <- dims
  }
  l <- list(centers=object@centers,
            variogram=variogram)
  if (length(object@sd) > 0) {
    l$sd <- object@sd
    if (is.array(variogram)) dim(l$sd) <- dims
  }
  
  if (length(object@phi.centers) > 1) l$phi <- object@phi.centers
  if (length(object@theta.centers) > 1) l$theta <- object@theta.centers
  if (length(object@T)>0) l$T <- object@T
  if (object@vdim > 1) l$vdim <- object@vdim  
  class(l) <- "summary.RFempVariog"
  l
}

summary.RF_empVariog <- function(object, ...) {
  variogram <- object$emp.vario
  if (is.array(variogram)) {
    dims <- dim(variogram)
    dims <- dims[dims > 1]
    dim(variogram) <- dims
  }
  l <- list(centers=object$centers,
            variogram=variogram)
  if (length(object$sd) > 0) {
    l$sd <- object$sd
    if (is.array(variogram)) dim(l$sd) <- dims
  }
  
  if (length(object$phi.centers) > 1) l$phi <- object$phi.centers
  if (length(object$theta.centers) > 1) l$theta <- object$theta.centers
  if (length(object$T)>0) l$T <- object$T
  if (object$vdim > 1) l$vdim <- object$vdim
  class(l) <- "summary.RFempVariog"
  l
}


print.summary.RFempVariog <- function(x, ...) {
  cat("Object of class 'RFempVariog'\n")
  str(x, no.list=TRUE) #
  invisible(x)
}

print.RFempVariog <- function(x, ...) {
  print.summary.RFempVariog(summary.RFempVariog(x, ...))
}
print.RF_empVariog <- function(x, ...) {
  print.summary.RFempVariog(summary.RF_empVariog(x, ...))
}


setMethod(f="show", signature="RFempVariog",
          definition=function(object) print.RFempVariog(object)) #

## coercion methods
setAs("RFempVariog", "list",
      def=function(from) {
        li <- list()
        for (name in slotNames(from)) {
          li[[name]] <- eval(parse(text=paste("from@", name, sep="")))
        }
        return(li)
      })

setAs("RFfit", "RFempVariog", def=function(from) from@ev)

list2RFempVariog <- function(li) {
  RFopt.coords <- RFoptions()$coords
  return(new("RFempVariog",
             centers=li$centers,
             emp.vario=li$emp.vario,
             var=var,
             sd=li$sd,
             n.bin=li$n.bin,
             phi.centers=li$phi.centers,
             theta.centers=li$theta.centers,
             T=li$Tbins,
             vdim = li$vdim,
             coord.units = RFopt.coords$coordunits,
             variab.units = RFopt.coords$varunits,
             call=li$call ))
}


## plot method

setMethod(f="plot", signature(x="RFempVariog", y="missing"),
	  definition=function(
            x, model = NULL, nmax.phi=NA, nmax.theta=NA, nmax.T=NA,
            plot.nbin=TRUE, plot.sd=FALSE, variogram=TRUE,
            boundaries = TRUE,...){ 
          plotRFempVariog(x, model = model,
                          nmax.phi=nmax.phi, nmax.theta=nmax.theta,
                          nmax.T=nmax.T, plot.nbin=plot.nbin, plot.sd=plot.sd,
                          variogram=variogram, boundaries = boundaries,
                          ...)
                   })

setMethod(f="plot", signature(x="RFfit", y="missing"),
	  definition=function(
            x, model = NULL,
            fit.method="ml", nmax.phi=NA, nmax.theta=NA, nmax.T=NA,
            plot.nbin=TRUE, plot.sd=FALSE, variogram=TRUE,
            boundaries = TRUE,...)
          
          plotRFempVariog(x, fit.method=fit.method,
                          nmax.phi=nmax.phi, nmax.theta=nmax.theta,
                          nmax.T=nmax.T, plot.nbin=plot.nbin, plot.sd=plot.sd,
                          model = model, variogram=variogram,
                          boundaries = boundaries,
                          ...))

setMethod(f="persp", signature(x="RFempVariog"),
	  definition=function(
            x, model = NULL, nmax.phi=NA, nmax.theta=NA, nmax.T=NA,
            plot.nbin=TRUE, plot.sd=FALSE, variogram=TRUE,
            boundaries = TRUE,...){ 
            plotRFempVariog(x, model = model,
                          nmax.phi=nmax.phi, nmax.theta=nmax.theta,
                          nmax.T=nmax.T, plot.nbin=plot.nbin, plot.sd=plot.sd,
                          variogram=variogram, boundaries = boundaries,
                          ..., plotmethod=graphics::persp)
                   })

setMethod(f="persp", signature(x="RFfit"),
	  definition=function(
            x, model = NULL,
            fit.method="ml", nmax.phi=NA, nmax.theta=NA, nmax.T=NA,
            plot.nbin=TRUE, plot.sd=FALSE, variogram=TRUE,
            boundaries = TRUE,...)          
          plotRFempVariog(x, fit.method=fit.method,
                          nmax.phi=nmax.phi, nmax.theta=nmax.theta,
                          nmax.T=nmax.T, plot.nbin=plot.nbin, plot.sd=plot.sd,
                          model = model, variogram=variogram,
                          boundaries = boundaries,
                          ..., plotmethod="persp"))




contour.RFfit <- contour.RFempVariog <- 
  function(x, model = NULL,
           fit.method="ml", nmax.phi=NA, nmax.theta=NA, nmax.T=NA,
           plot.nbin=TRUE, plot.sd=FALSE, variogram=TRUE,
           boundaries = TRUE,...) {
    stopifnot(!( (is(x, "RFfit") && is.list(x@ev[[1]]@centers))
                || (is(x, "RFempVariog") && is.list(x@centers))
                ))
    plotRFempVariog(x, fit.method=fit.method,
                    nmax.phi=nmax.phi, nmax.theta=nmax.theta,
                    nmax.T=nmax.T, plot.nbin=plot.nbin, plot.sd=plot.sd,
                    model = model, variogram=variogram,
                    boundaries = boundaries,
                    ..., plotmethod="contour")
  }


plotRFempVariogUnbinned <- function(x, coord.units, variab.units, variab.names,
                                    ..., plotmethod="image") {
  dots = list(...)
  dotnames <- names(dots)
  coords <- GridTopology2gridVectors(cbind(x@centers$x, x@centers$T))
  if (length(coords)>1) {
    coords[[1]] <- sort(unique((coords[[1]] - min(coords[[1]])) *
                               rep(c(-1, 1), each=length(coords[[1]]))))
    lab.names <- dimnames(x@centers$x)[[2]]
    
    if (!(x@centers$spacedim==1 && x@centers$Time) &&
        !(x@centers$x[3,2]==dim(x@emp.vario)[2]))
      coords[[2]] <- sort(unique((coords[[2]] - min(coords[[2]])) *
                                 rep(c(-1, 1), each=length(coords[[2]]))))
    else
      lab.names[2] <- "T"
  } else {
    x@centers <- coords[[1]]-min(coords[[1]])
    do.call(graphics::plot, args=c(dots, list(x=x, plot.nbin=FALSE)))
    return()
  }
  
  if (!("main" %in% dotnames)) {
    main <- "Variogram image plot"
    if (length(variab.names)>0) main <- paste(main, "for", variab.names)
    dots$main <- main
  }
  lab.names <- paste(lab.names, "-distance", sep="")
  
  idx <- lab.names != "T-distance"
  if (any(idx) && all(coord.units[idx] != ""))
      lab.names[idx] <-
        paste(lab.names[idx], " [", coord.units[idx], "]", sep="")
  if (!all(idx) && all(coord.units[!idx] != ""))
    lab.names[!idx] <-
      paste(lab.names[!idx], " [", coord.units[!idx], "]", sep="")
  
  if (!("xlab" %in% dotnames)) dots$xlab <- lab.names[1]
  if (!("ylab" %in% dotnames)) dots$ylab <- lab.names[2]
  if (!("xlim" %in% dotnames)) dots$xlim <- range(coords[[1]])
    if (!("ylim" %in% dotnames)) dots$ylim <- range(coords[[2]])
  
  idx1 <- coords[[1]] >= dots$xlim[1] & coords[[1]] <= dots$xlim[2]
  idx2 <- coords[[2]] >= dots$ylim[1] & coords[[2]] <= dots$ylim[2]
  coords[[1]] <- coords[[1]][idx1]
  coords[[2]] <- coords[[2]][idx2]
  
  dims <- dim(x@emp.vario)
  ev.plot <- matrix(x@emp.vario, nrow=dims[1], ncol=dims[2])
  ev.plot <- ev.plot[idx1, idx2]
  
  range.ev <- range(ev.plot)#x@emp.vario)
  col <- if ("col" %in% dotnames) dots$col else
  default.image.par(NULL, NULL)$default.col
  
  split.screen( rbind(c(0,.85,0,1), c(.85,1,0,1)))
  screen(2)
  par(mar=c(4,0,3,3))
  z.legend <- seq(range.ev[1], range.ev[2], length=length(col))
    image(y=z.legend, x=c(0,1), z=rbind(z.legend, z.legend), axes=F, col=col,
          xlab="")
  axis(4)
  box()
  screen(1)
  par(mar=c(4,4,3,1))
  dots$type <- NULL

  do.call(plotmethod,
          args=c(dots, list(x=coords[[1]], y=coords[[2]], z=ev.plot)))
                                        #      xlab=lab.names[1],
                                        #      ylab=lab.names[2],
                                        #      main="Variogram image plot")
  close.screen(all.screens=TRUE)
  return(invisible())
}

plotRFempVariog <- function(x, model, nmax.phi, nmax.theta, nmax.T,
                            plot.nbin, plot.sd, fit.method = NULL,
                            variogram=variogram,
                            coord.units=c(""),
                            variab.units=c(""),
                            boundaries = TRUE,
                            ...) {
  
  if (!variogram)
    stop("plot of estimated covariance functions not programmed yet.")

  newx <- list()
  methodnames <- double()
  if (is(x, "RFfit")) {    
    if(length(x@ev)==0) stop("The fit does not contain an empirical variogram.")
    newx$autostart <- x@autostart
    newx$self <- x@self
    newx$plain <- x@plain
    newx$sqrt.nr <- x@sqrt.nr
    newx$sd.inv <- x@sd.inv
    newx$internal1 <- x@internal1
    newx$internal2 <- x@internal2
    newx$internal3 <- x@internal3
    newx$ml <- x@ml
    x <- x@ev[[1]]
    if( is.null(fit.method) )
      fit.method <- if (length(newx$ml@name) > 0) "ml" else "plain"

    methodidx <- fit.method %in% names(newx[unlist(lapply(newx, FUN = function(x) is(x, ZF_MODELEXT)))])
    methodnames <- fit.method[methodidx]
    nomethodnames <- fit.method[!methodidx]

    if( !all(methodidx) )
      warning( paste("The following method does not exist: ", nomethodnames) )
  } else if (!is(x, "RFempVariog"))
    stop("method only for objects of class 'RFempVariog' or 'RFfit'")

  variab.names <-
    if (is.matrix(x@emp.vario)) dimnames(x@emp.vario)[[2]][1]
    else names(x@emp.vario)[1]
 
  ## case without binning
  if (is.list(x@centers))
    return(plotRFempVariogUnbinned(x=x,                         
                                   coord.units=coord.units,
                                   variab.units=variab.units,
                                   variab.names=variab.names, ...))

             
  dots = list(...)
  dotnames <- names(dots)
   if (!("type" %in% dotnames)) dots$type <- "b" 
  cex <- if ("cex" %in% dotnames) dots$cex else .8
  par(cex=cex, xaxs="i")
  dots$cex <- NULL
  if (!("pch" %in% dotnames)) dots$pch <- 19
             
  if(!("xlim" %in% dotnames)) dots$xlim <- range(x@centers)
  ylim.not.in.dotnames <- !("ylim" %in% dotnames)
  xlab <- if ("xlab" %in% dotnames) dots$xlab else "distance"
  dots$xlab <- NULL
  ylab.ev <- if ("ylab" %in% dotnames) dots$ylab else "semivariance"
  dots$ylab <- NULL
  
  main0 <- if ("main" %in% dotnames) dots$main else "Variogram plot"  
  dots$main <- NULL
  if (!is.null(main0)) oma.top <- 2 else oma.top <- 0
  

  has.sd <- !is.null(x@sd)
  model <- c(model)
  if(!is.null(model))
    if (!all(unlist(lapply(model, FUN=function(x) is(x, ZF_MODEL)))))
      stop("model must be (a list of elements) of class 'ZF_MODEL'")
  if(!is.null(model)){
    modelnames <-
      if(length(names(model)) > 0) names(model)
      else paste("model", 1:length(model))
    methodnames <- c(methodnames, modelnames)
    names(model) <- modelnames
    newx <- c(newx, model)
  }
  n.methods <- length(methodnames)
  
    
  n.phi <- min(nmax.phi, l.phi <- max(1,length(x@phi.centers)), na.rm=TRUE)
  n.theta <- min(nmax.theta, l.theta <- max(1,length(x@theta.centers)),
                 na.rm=TRUE)
  n.T <- min(nmax.T, l.T <- max(1,length(x@T)), na.rm=TRUE)
  vdim <- dim(x@emp.vario)
  vdim <- vdim[length(vdim)]

  if (n.phi > 6 || n.theta > 3 || n.T > 3)
    message("'If you feel the picture is overloaded, set the parameters 'nmax.phi', 'nmax.theta' and 'nmax.T'")
  
  halfphi.sector <- pi/(2*l.phi)
  halftheta.sector <- pi/(2*l.theta)
  phi.angles <- c(halfphi.sector, 0, -halfphi.sector) 
  theta.angles <- seq(-halftheta.sector, halftheta.sector, len=5) # len=odd!
  if (n.phi>1 && boundaries) {
    phi.angles <- phi.angles * 0.96
    theta.angles <- theta.angles * 0.96
  } 
  
  TandV <- n.T > 1 && vdim > 1
  if (vdim>1 && length(variab.names)==0)
    variab.names <- paste("v", 1:vdim, sep="")

  range.nbin <- range(c(0, x@n.bin), na.rm=TRUE)
  ylim.nbin <- range.nbin * c(1,1.2)
  
  col.v <- col <- 
    if ("col" %in% dotnames) rep(dots$col, len=n.phi) else 1:max(n.phi)
  dots$col <- NULL

  if (n.methods > 0){
    dotsRFfit <- dots
    dotsRFfit$type <- "l"
    dotsRFfit$lwd <- 2
    ltyRFfit.v <- 1:n.methods
    dotsRFfit$lty <- NULL
  }


  dir2vario <- function(dir.vec, x.eval, x.time, method.model, v1, v2){
    x.space <- as.matrix(x.eval) %*% t(dir.vec)            
    for (j in (1:4)) { ## orig 20
      vario.vals <- try(RFvariogram(x = cbind(x.space, x.time),
                                    model = method.model,
                                    grid = FALSE),
                        silent = FALSE)
      
      if(!is(vario.vals, "try-error")) {
        if (is.vector(vario.vals)) {
          return(vario.vals)
        } else {
          return(vario.vals[, v1, v2])
        }
      }      
      x.space <- cbind(x.space, 0)    
    }
    return(NULL)
  }

  
  oma.left <- 6

  Screens <- if (!TandV) c(n.T * vdim * vdim, n.theta) else  c(n.T, n.theta)
  n.all <- prod(Screens)
  if (any(par()$mfcol != c(1,1))) par(mfcol=c(1,1))

  for (v1 in 1:vdim) {
    for (v2 in 1:vdim) {
      if (TandV || (v1==1 && v2==1)) scr <- split.screen(Screens)
      if (vdim == 1) {
        main <-
          if (!is.null(main0) && length(variab.names)>0)
            paste(main0, "for", variab.names) else main0
      } else {
        main <-
          if (!TandV) main0
          else paste(main0, "for", variab.names[v1], "vs.",  variab.names[v2])
      }      
      if (ylim.not.in.dotnames)
        dots$ylim <- range(x@emp.vario[,,,, v1, v2], na.rm=TRUE)
      for (iT in 1:n.T) {
        for (ith in 1:n.theta) {
          ## plot n.bin
          if (plot.nbin) {
            screen(scr[1])
            par(oma=c(4,oma.left,oma.top,0))
            scr2 <- split.screen(rbind(c(0,1,.2,1), c(0,1,0,.2)), screen=scr[1])
            screen(scr2[2])
            par(mar=c(0,.5,0,.5))
            for (iph in 1:n.phi) {
              if (n.phi > 1) col <- col.v[iph]
              if (iph==1) {
                lab <- xylabs("bin centers", NULL, units=x@coord.units)
                plot(x@centers, x@n.bin[ ,iph, ith, iT, v1, v2],
                     xlim=dots$xlim, ylim=ylim.nbin,
                     type=if (n.phi>1) "p" else "h",
                     col =if (n.phi>1) col else "darkgray", lwd=4,
                     pch=16, axes=FALSE, ylab="", xlab = lab$x
                     )
                box()
                at <- seq(range.nbin[1], range.nbin[2], len=3)
                if (ith==1)
                  axis(2, las=1, at=at, labels=format(at, scient=-1, dig=2),
                       outer=TRUE)
                else  axis(2, las=1, at=at, labels=FALSE)
                
                if (iT==n.T && (n.T > 1 || (v1==vdim && v2==vdim))) axis(1)
                if (ith==1) title(ylab="n.bin", line=5, outer=TRUE, adj=0)
                
              } else {
                points(x@centers, x@n.bin[ ,iph, ith, iT, v1, v2],
                       type="p", col=col, pch=16)
              }
            }
            screen(scr2[1])
          } else {
            screen(scr[1])
            par(oma=c(4,oma.left,oma.top,0))
          }
          
          ## plot emp.vario
          ##if (ith==1) par(mar=c(0,6,1,1)) else
          par(mar=c(0,.5,1,.5))
          
          plotted.meth <- NULL  # needed for legend
          
          for (iph in 1:n.phi) {
            if (n.phi>1) col <- col.v[iph]
            if (iph==1) {
              do.call(graphics::plot, args=c(dots, list(
                              x=x@centers, y=x@emp.vario[ ,iph, ith, iT, v1,v2],
                                        #ylim=ylim.ev, type=type, pch=19,
                              col=col, axes=FALSE, ylab=""))
                      )
              box()
              axis(2, las=1, labels=(ith==1), outer=(ith==1))
              if (!plot.nbin) axis(1)          
              if (l.theta > 1 || l.T > 1 || vdim > 1) {
                L <- character(3)
                if (!TandV && vdim > 1) L[1] <- paste(variab.names[c(v1,v2)],
                                                      collapse=":")
                if (l.T>1) L[2] <- paste("T=",signif(x@T[iT],3)," ",sep="")
                if (l.theta>1) L[3] <- paste(sep="", "theta=",
                                             signif(x@theta.centers[ith],3))
                legend("topleft", legend=paste(L[L!=""], collapse=","))
              }
              if (ith == 1) title(ylab=ylab.ev, line=5, outer=TRUE)
              if (has.sd && plot.sd)
                legend("topright", bty="n", #col=NA, lwd=NA,
                       legend="arrows represent +-1 sd intervals")
              
            } else {  ## iph > 1
              do.call(graphics::points, args=c(dots, list(
                                x=x@centers,
                                y=x@emp.vario[ ,iph, ith, iT, v1, v2], col=col))
                      )   #type=type,  pch=19)
            } ## if iph 
            
            k <- 1  
            if (n.methods > 0) {          
              for(i in 1:n.methods) {        
                method.model <- newx[[methodnames[i]]]
                if(length(method.model@name) == 0){
                  warning("The method '", methodnames[i], "' was not fitted.")
                  next
                }
                if (!is.null(x@phi.centers)){
                  x.radial <- cbind(cos(x@phi.centers[iph]+phi.angles),
                                    sin(x@phi.centers[iph]+phi.angles))
                  if(!is.null(x@theta.centers))
                    x.radial <- cbind(x.radial,
                                      cos(x@theta.centers[ith] +
                                          rep(theta.angles,
                                              each= length(phi.angles))))
                  
                } else x.radial <- matrix(1, nrow=1, ncol=1)
                
                x.time <- NULL
                if(!is.null(x@T)) x.time <- x@T[iT]
                x.eval <- seq(from = max(dotsRFfit$xlim[1],1e-3),
                              to = dotsRFfit$xlim[2], len = 150) 
                
                ## sehr genau Abschaetzung, indem mehrere (3)
                ## Winkel angeschaut werden und dann der mittlere
                ## Wert der Variogramme angeschaut wird, da ja auch
                ## das emp. Variogramm ueber ein Winkelintervall gemittelt wird
                dummy.vals <- apply(x.radial, 1, FUN = dir2vario,
                                    x.eval=x.eval, x.time=x.time,
                                    method.model=method.model,
                                    v1=v1, v2=v2)                
                if(!is.null(dummy.vals)){
                  if (n.phi>1 && boundaries)
                    do.call(graphics::matplot,
                            args=c(dotsRFfit, list(x=x.eval,
                              y=t(apply(dummy.vals, 1, range)),
                              add = TRUE, col=col, lty = 3)))
                  else {
                    do.call(graphics::points,
                            args=c(dotsRFfit, list( x=x.eval,
                              y = rowMeans(as.matrix(dummy.vals))  , 
                              col=col, lty = ltyRFfit.v[k])))
                  }                
                  k <- k+1
                  if(iph == 1) plotted.meth <- c(plotted.meth, methodnames[i])
                }
                
              }
            }  
            
            
            if (has.sd && plot.sd) {
              sdnot0 <-  x@sd[ ,iph, ith, iT] != 0
              arrows(x@centers[sdnot0],
                     x@emp.vario[sdnot0 ,iph, ith, iT] - x@sd[sdnot0 ,iph, ith, iT],
                     x@centers[sdnot0],
                     x@emp.vario[sdnot0 ,iph, ith, iT] + x@sd[sdnot0 ,iph, ith, iT],
                     code=2, angle=90, length=0.05, col=col)
            }
            
          } # nphi
          
          
          pifactor <- signif((x@phi.centers[1:n.phi]) / pi, 2)
          
          len.mnames <- length(plotted.meth)
          string.emp <- "empirical"
           if(len.mnames > 0) {
            labels <-
              if (n.phi>1) paste('"phi=",', rep(pifactor, each = len.mnames+1),
                            ', pi, ', '", ')
              else '" '
            labels <- paste("c(",
                            paste('expression(paste(',labels,
                                  rep(c(string.emp, plotted.meth), l.phi),
                                  '"', '))', collapse=","),
                            ")")
            labels <- eval(parse(text=labels))
          } else
            labels <- if (n.phi > 1) 
              eval(parse(text=paste(
                           "c(",
                           paste('expression(paste("phi=",', pifactor, ', pi))',
                                 collapse=","),
                           ")")))
          
          
          if (l.phi > 1 || len.mnames > 0)# && iT==1 && ith==1) # auskommentiert auf Sebs wunsch        
            legend("bottomright", col=rep(col.v, each = len.mnames+1),
                   lwd=1, pch=rep(c(19,rep(NA, len.mnames)), l.phi), bty="n",
                   legend=labels, lty =
                   rep(c(1, if(len.mnames==0) NULL else ltyRFfit.v[1:len.mnames]),
                   l.phi))
          
          scr <- scr[-1]
        } # n.theta
      } # T
    } # vdim 1
  } # vdim 2
  dots$type <- NULL
  
  if (!is.null(main))
    do.call(graphics::title, args=c(dots, main=main, outer=TRUE))
  if (!is.null(xlab))
    do.call(graphics::title, args=c(dots, xlab=xlab, outer=TRUE))

  close.screen(all.screens=TRUE)
}
