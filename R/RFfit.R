
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017 Martin Schlather
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



### !!!!!!!!!!!! ACHTUNG !!!!!!!!!!!! TREND als cov-fct muss
### noch programmiert werden !!!

##   source("~/R/RF/RandomFields/R/MLES.R")

## PrintLevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information

## jetzt nur noch global naturalscaling (ja / nein)
## spaeter eine Funktion schreibbar, die den naturscaling umwandelt;
##   im prinzipt CMbuild, aber ruechwaers mit 1/newscale und eingefuegt
##   in eventuell schon vorhandene $ operatoren


#Beim paper lesen im Zug nach Muenchen heute morgen ist mir eine Referenz zu einem R Paket "mlegp: Maximum likelihood estimates of Gaussian processes" aufgefallen. Ist Dir aber sicher schon bekannt! 

#  stop("")
  # problem: natscale; im moment 2x implementiert, 1x mal ueber
  # scale/aniso (user) und einmal gedoppelt -- irgendwas muss raus

## LSQ variogram fuer trend = const.
## kann verbessert werden, insb. fuer fixed effects, aber auch eingeschraenkt
## fuer random effects -> BA/MA


## REML fehlt

## users.guess muss in eine List von meheren Vorschlaegen umgewandelt werden !!! Und dann muss RFfit recursiver call mit allen bisherigen Werden laufen !!


## zentrale C -Schnittstellen
##    .C(C_PutValuesAtNA, RegNr, param)

## bins bei Distances automatisch


## bei repet sind die Trends/fixed effects gleich, es muessen aber die
## random effects unterschiedlich sein.
## bei list(data) werden auch trend/fixed effects unterschiedlich geschaetzt.


## Erweiterungen: Emilio's Bi-MLE, Covarianz-Matrix-INversion per fft oder
## per INLA, grosse Datensaetze spalten in kleinere "unabhaengige".


###################################
## !!! Mixed Model Equations !!! ##
###################################

## accessing slots
accessByNameOrNumber <- function(x, i, j, drop=FALSE) {
  stopifnot(length(i)==1)
  if (is.numeric(i))  i <- slotNames(x)[i]
  return(accessSlotsByName(x=x, i=i, j=j, drop=drop))
}

setMethod("[", signature = "RFfit", def=accessByNameOrNumber)


### to do : ask Paulo
#effects_RFfit <- function(OP, object, method) {
#  eff <- RFrandef(object=object, method=method, OP=OP)
#  linpart <- fitted_RFfit(OP=OP, object=object, method=method)
#  stop("unclear how these two results should be combined in the output")
#}
#effects_RMmodelFit <- function(...) stop("'effects' can only be used with the original and sp_conform output of 'RFfit'.")
#setMethod(f="effects", signature='RFfit',
#          definition=function(object, method="ml")
#          effects_RFfit("@", object=object, method=method))#
#setMethod(f="effects", signature='RMmodelFit',
#          definition=function(object, newdata=NULL) effects_RMmodelFit())#
#effects.RM_modelFit <- function(object, ...) effects_RMmodelFit()
#effects.RF_fit <- function(object, method="ml") effects_RMmodelFit()


simulate_RFfit <- function(OP, object, newdata, conditional, method) {
  Z <- do.call(OP, list(object, "Z"))
  L <- length(Z$data)
  ans <- rep(list(NULL), L)
  m <- ModelParts(object[method], effects=Z$effect, complete=FALSE) ## no params
  if (conditional) {
    for (i in 1:L) {
      ans[[i]] <- RFsimulate(model=m$model,
                             data = Z$data[[i]],
                             given = Z$coord[[i]],
                             x = if (!is.null(newdata)) newdata,
                             err.model = m$err.model)
        if (is.list(ans[[i]])) stop("the case with NAs not completely programmed yet. Please let the maintainer now that it is needed ")
    }
    return (if (L == 1) ans[[1]] else ans)
  } else {
    for (i in 1:L) {
       ans[[i]] <- RFsimulate(model=,
                              x = if (length(newdata)==0) Z$coord[[i]]
                              else newdata)
      if (is.list(ans[[i]])) stop("the case with NAs not completely programmed yet. Please let the maintainer now that it is needed ")
    }
  }
}
simulate_RMmodelFit <- function(...) stop("'simulate' can only be used with the original and sp_conform output of 'RFfit'.")
setMethod(f="simulate", signature='RFfit',
          definition=function(object, newdata=NULL,
              conditional=!is.null(newdata),  method="ml")
          simulate_RFfit("@", object=object, conditional=conditional,
                         method=method))#
setMethod(f="simulate", signature='RMmodelFit',
          definition=function(object, newdata=NULL) simulate_RMmodelFit())#
simulate.RM_modelFit <- function(object, ...) simulate_RMmodelFit()
simulate.RF_fit <- function(object, method="ml") simulate_RMmodelFit()


predict_RFfit <- function(OP, object, newdata, impute, method) {
  Z <- do.call(OP, list(object, "Z"))
  L <- length(Z$data)
  ans <- rep(list(NULL), L)
  if (impute) {
    if (length(newdata) > 0) stop("for imputing, 'newdata' may not be given")
    for (i in 1:L) {
      ans[[i]] <- RFinterpolate(model=object[method],
                                data = Z$data[[i]],
                                given = Z$coord[[i]],
                                err.model = NA)
      if (is.list(ans[[i]])) stop("the case with NAs not completely programmed yet. Please let the maintainer now that it is needed ")
    }
  } else {
    for (i in 1:L) {
      ans[[i]] <- RFinterpolate(model=object[method],
                                data = Z$data[[i]],
                                given = Z$coord[[i]],
                                x = if (!is.null(newdata)) newdata else
                                     Z$coord,
                                err.model = NA)
      if (is.list(ans[[i]])) stop("the case with NAs not completely programmed yet. Please let the maintainer now that it is needed ")
    }
  }
  return (if (L == 1) ans[[1]] else ans)
}
predict_RMmodelFit <- function(...) stop("'predict' can only be used with the original and sp_conform output of 'RFfit.")
setMethod(f="predict", signature='RFfit',
          definition=function(object, newdata=NULL, impute=FALSE, method="ml")
          predict_RFfit("@", object=object, newdata=newdata,
                        impute=impute, method=method))#
setMethod(f="predict", signature='RMmodelFit',
          definition=function(object, newdata=NULL) predict_RMmodelFit())#
predict.RM_modelFit <- function(object, ...)
  predict_RMmodelFit(object=object, ...)
predict.RF_fit <- function(object, method="ml")
  predict_RMmodelFit(object, method=method)



coef_RMmodelFit <- function(OP, object) {
  covariat <- do.call(OP, list(object, "covariat"))
  glbl.var <- do.call(OP, list(object, "globalvariance"))
  p <- do.call(OP, list(object, "param"))
  if (length(covariat) > 0)  covariat <- as.matrix(covariat)
  nr_p <- nrow(p)
  if (length(glbl.var) > 0)
    glbl.var <- c(glbl.var, rep(NA, nr_p - length(glbl.var)))
  p <- cbind(p, glbl.var,
             if (length(covariat) > 0)
             rbind(covariat, matrix(NA, ncol=ncol(covariat),
                                    nrow= nr_p - nrow(covariat))))
  #class(p) <- "coef.RMmodelFit"
  p[1, ]
}
setMethod(f="coef", signature='RMmodelFit',
          definition=function(object) coef_RMmodelFit("@", object))#
setMethod(f="coef", signature='RFfit',
          definition=function(object, method="ml")
          coef_RMmodelFit("@", object[method]))#
coef.RM_modelFit <- function(object, ...)
  coef_RMmodelFit("$", object)
coef.RF_fit <- function(object, method="ml")
  coef_RMmodelFit("$", object[method])



residuals_RMmodelFit <- function(OP, object) {
  resid <- do.call(OP, list(object, "residuals"))
  message("Note that 'residuals' equals the difference between the data and the linear part (fixed effects).")
  if (length(resid) == 1) resid[[1]] else resid
}
setMethod(f="residuals", signature='RMmodelFit',
          definition=function(object) residuals_RMmodelFit("@", object))#
setMethod(f="residuals", signature='RFfit',
          definition=function(object, method="ml")
          residuals_RMmodelFit("@", object[method]))#
residuals.RM_modelFit <- function(object, ...)
  residuals_RMmodelFit("$", object)
residuals.RF_fit <- function(object, method="ml")
  residuals_RMmodelFit("$", object[method])



fitted_RFfit <- function(OP, object, method) {
  data <- do.call(OP, list(object, "Z"))$data
  resid <- do.call(OP, list(object[method], "residuals"))
  for (i in 1:length(data)) data[[i]] <- data[[i]] - resid[[i]]
  message("Note that 'fitted' equals the linear part (fixed effects).")
  if (length(data) > 1) data else
  if (ncol(data[[1]]) > 1)  data[[1]] else as.vector(data[[1]])
}
fitted_RMmodelFit <- function(...) stop("'fitted' can only be used with the original output of 'RFfit', not with some of its extraction.")
setMethod(f="fitted", signature='RMmodelFit',
          definition=function(object) fitted_RMmodelFit())#
setMethod(f="fitted", signature='RFfit',
          definition=function(object, method="ml")
          fitted_RFfit("@", object=object, method=method))#
fitted.RM_modelFit <- function(object, ...) fitted_RMmodelFit()
fitted.RF_fit <- function(object, method="ml")
  fitted_RFfit("$", object=object, method=method)



RFhessian <- function(model) {
  method <- "ml"
  if (is(model, "RF_fit")) return(model[[method]]@hessian)
  else if (is(model, "RFfit")) return(model[method]$hessian)
  else stop("'model' is not an output of 'RFfit'")
}

anova.RFfit <- function(object, ...)  RFratiotest(nullmodel=object, ...)
anova.RF_fit <- function(object, ...) RFratiotest(nullmodel=object, ...)
anova.RMmodelFit <- function(object, ...) RFratiotest(nullmodel=object, ...)
anova.RM_modelFit <- function(object, ...) RFratiotest(nullmodel=object, ...)

setMethod(f="anova", signature=CLASS_FIT, anova.RFfit)#
setMethod(f="anova", signature='RFfit', anova.RFfit)#

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
      paste(sep="", "Note that the (possibly internal) fitted variable",
            if (nl > 0)
            paste(if (nl > 1) "s " else " ",
                  paste("'", colnames(variab)[lidx], "'", sep="", collapse=", "),
                  if (nl == 1)  " is " else " are ",                  
                  "close to or on the effective lower boundary", sep=""),
            if (nl > 0 && nu > 0) " and the variable",
            if (nu > 0)
            paste(if (nu > 1) "s " else " ", 
                  paste("'", colnames(variab)[uidx], "'",
                        sep="", collapse=", "),
                  if (nu == 1) "is" else "are",
                  "close to or on the effective upper boundary"),
            ".\nHence the gradient of the likelihood function might not be zero and none of the\nreported 'sd' values might be reliable.")
  } else txt <- NULL
  return(txt)
}


  
summary_RMmodelFit <- function(OP, object, ..., isna.param) {
  model <- if (OP == "@") PrepareModel2(object, ...) else object$model
  covariat <- do.call(OP, list(object, "covariat"))
  glbl.var <- do.call(OP, list(object, "globalvariance"))
  p <- do.call(OP, list(object, "param"))
  r <- do.call(OP, list(object, "residuals"))
  v <- do.call(OP, list(object, "variab"))
  l <- list(model=model,
            loglikelihood=do.call(OP, list(object, "likelihood")),
            AIC = do.call(OP, list(object, "AIC")),
            AICc= do.call(OP, list(object, "AICc")),
            BIC = do.call(OP, list(object, "BIC")),
            residuals=if (length(r) == 1) r[[1]] else r)
  if (missing(isna.param)) isna.param <- any(is.na(p)) 
  l$boundary <- boundary_values(v)
  if (length(covariat) > 0)  covariat <- as.matrix(covariat)
  if (!any(is.na(p[1, ]))) {
    nr_p <- nrow(p)
    if (length(glbl.var) > 0)
      glbl.var <- c(glbl.var, rep(NA, nr_p - length(glbl.var)))
    l$param <- cbind(p, glbl.var,
                     if (length(covariat) > 0)
                     rbind(covariat, matrix(NA, ncol=ncol(covariat),
                                            nrow= nr_p - nrow(covariat))))
  }
  if (isna.param || !is.null(l$boundary)) {
    nr_v <- nrow(v)
    if (length(glbl.var) > 0)
      glbl.var <- c(glbl.var, rep(NA, nr_v - length(glbl.var)))
    l$variab <- cbind(v, glbl.var,
                      if (length(covariat) > 0)
                      rbind(covariat, matrix(NA, ncol=ncol(covariat),
                                             nrow=nr_v - nrow(covariat)))
                      )
   }
  class(l) <- "summary.RMmodelFit"
  l
}

summary.RMmodelFit <- function(object, ..., isna.param) {
  summary_RMmodelFit("@", object, ..., isna.param=isna.param)
}

setMethod(f="summary", signature=CLASS_FIT, summary.RMmodelFit)#

summary.RM_modelFit <- function(object, ..., isna.param) {
  summary_RMmodelFit("$", object, ..., isna.param=isna.param)
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
  if (length(x$submodels) > 0) {
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
        param[, sm$p.proj] <- sm$param
        fixed <- sm$fixed
        if (length(fixed) > 0) {
          param[1, fixed$zero] <- 0
          param[1, fixed$one] <- 1
        }
      
      #  if (!is.null(cparam)) cparam <- rbind(cparam, NA)
        cparam <- rbind(cparam, param)
  
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

  np <- NA
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

setMethod(f="show", signature=CLASS_FIT,
          definition=function(object) print.RMmodelFit(object))#

  
summary.RFfit <- function(object, ...,  method="ml", full=FALSE) {
  s <- summary.RMmodelFit(object[method])
  len <-  length(object@submodels)
  if (full && length(object@submodels) > 0) {
    submodels <- list()
    for (i in 1:len) {
      ## war summary.RM_modelFit
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
  val <- object[[method]]$likelihood  
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
  if (ats2$df < 0) ats2 <- NULL
  ats <- ats[, paste("model1.", values, sep="")]
  colnames(ats) <- values
  if (ats$df < 0) ats <- NULL
  ats <- unique(rbind(ats, ats2))
  dimnames(ats) <- list(1:nrow(ats), colnames(ats))

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



coef_RMmodelFit <- function(OP, object, ...) {
  covariat <- do.call(OP, list(object, "covariat"))
  glbl.var <- do.call(OP, list(object, "globalvariance"))
  p <- do.call(OP, list(object, "param"))
  if (length(covariat) > 0)  covariat <- as.matrix(covariat)
  nr_p <- nrow(p)
  if (length(glbl.var) > 0)
    glbl.var <- c(glbl.var, rep(NA, nr_p - length(glbl.var)))
  p <- cbind(p, glbl.var,
             if (length(covariat) > 0)
             rbind(covariat, matrix(NA, ncol=ncol(covariat),
                                    nrow= nr_p - nrow(covariat))))
  #class(p) <- "coef.RMmodelFit"
  p[1, ]
}

setMethod(f="coef", signature='RMmodelFit',
          definition=function(object) coef_RMmodelFit("@", object))#

setMethod(f="coef", signature='RFfit',
          definition=function(object) coef_RMmodelFit("@", object["ml"]))#

coef.RM_modelFit <- function(object, ...) coef_RMmodelFit("$", object, ...)
coef.RF_fit <- function(object, ...) coef_RMmodelFit("$", object["ml"], ...)



setMethod(f="plot", signature(x="RFfit", y="missing"),
          function(x, y, ...) RFplotEmpVariogram(x, ...))
setMethod(f="persp", signature(x="RFfit"),
	  function(x, ...) RFplotEmpVariogram(x, ..., plotmethod="persp"))


contour.RFfit <- contour.RFempVariog <- 
  function(x,...) {
    stopifnot(!( (is(x, "RFfit") && is.list(x@ev@centers))
                || (is(x, "RFempVariog") && is.list(x@centers))
                ))
    RFplotEmpVariogram(x, ..., plotmethod="contour")
  }


ExpliciteGauss <- function(model) {
  if (model[[1]] != "RPgauss" && model[[1]] != "gauss.process") {
    boxcox <- RFoptions()$gauss$boxcox
    if (any(is.na(boxcox)) || any(boxcox[c(TRUE, FALSE)] != Inf))
      return(list("RPgauss", boxcox=boxcox, model))
  }
  return(model)
}

RFfit <-
  function(model, x, y=NULL, z=NULL, T=NULL,  grid=NULL, data, 
           lower=NULL, upper=NULL, 
           methods, # "reml", "rml1"),
           sub.methods,
           ## "internal" : name should not be changed; should always be last
           ##              method!
           optim.control=NULL,
           users.guess=NULL,  
           distances=NULL, dim,
           transform=NULL,
           params=NULL,
            ##type = c("Gauss", "BrownResnick", "Smith", "Schlather",
           ##             "Poisson"),
           ...)
{

  .C(C_NoCurrentRegister)
  
  RFoptOld <- internal.rfoptions(xyz=length(y)!=0,...,
                                 internal.examples_reduced = FALSE,
                                 RELAX=is(model, "formula"))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  if (length(params) > 0) {
    if ((!is.na(RFopt$fit$estimate_variance_globally) &&
         RFopt$fit$estimate_variance_globally) &&
        RFopt$basic$printlevel > 0)
      message("Value of option 'hestimate_variance_globally' is ignored.")
    RFopt$fit$estimate_variance_globally <- FALSE
    RFoptions(fit.estimate_variance_globally = FALSE)
  }
  fit <- RFopt$fit
  
  if (RFopt$general$vdim_close_together)
    stop("'vdim_close_together' must be FALSE")

  if (is.data.frame(data)) {
    name <- "RFfit.user.dataset"
    do.call("attach", list(what=data, name=name))
    on.exit(detach(name, character.only = TRUE), add=TRUE)
  }

  ## in UnifyData the further.models that contain only the parameter data
  ## are turned into genuine models 
  further.models <- list()
  models <- c("lower", "upper", "users.guess", "parscale")
  if (paramlist <- length(params) > 0) {
    parscale <- optim.control$parscale
    for (m in models) {
      fm <- get(m)
      if (!is.null(fm) && !is.numeric(fm))
        further.models[[m]] <- PrepareModel2(fm, ...)
    }
  }  


##  Print(further.models, model)
  
  Z <- UnifyData(model=model, x=x, y=y, z=z, T=T, grid=grid,
                 data=data, distances=distances, dim=dim,
                 RFopt=RFopt,
                 mindist_pts = RFopt$fit$smalldataset / 2,
		 further.models = further.models, params=params, ...)

  ## Print(Z); kkk

  Z <- BigDataSplit(Z, RFopt)
  
 if (!hasArg("transform")) transform <- NULL
  if (paramlist) {
    for (m in models) 
      if (!is.null(get(m)) && !is.numeric(get(m)))
        assign(m, Z$further.models[[m]])
    optim.control$parscale <- parscale
    if (!is.null(Z$transform)) {
      if (!is.null(transform))
        stop("argument 'transform' may not be given if 'params' is given")
      transform <- Z$transform
    }
  } else {    
    parscale <- optim.control$parscale
    for (m in models)
      if (!is.null(get(m)) && !is.numeric(get(m)))
       assign(m, ReplaceC(PrepareModel2(get(m), ...)))
    optim.control$parscale <- parscale
  }


  new.model <- Z$model
  if (new.model[[1]] %in% c("RPpoisson", "poisson")) {
    res <- fit.poisson()
  } else if (new.model[[1]] %in% c("BRmixed", "BRshifted", "BRmixedIntern",
                               "RFbrownresnick")) {
    res <- fit.br()
  } else if (new.model[[1]] %in% c("RPschlather", "extremalgauss")) {
    res <- fit.extremal.gauss()
  } else if (new.model[[1]] %in% c("RPsmith", "smith")) {
    res <- fit.smith()
  } else if (new.model[[1]] %in% c("RPbernoulli", "binaryprocess")) {
    res <- fit.bernoulli()    
  } else {   
    Z$model <- ExpliciteGauss(ReplaceC(Z$model))
    res <- do.call("rffit.gauss",
            c(list(Z, lower=lower, upper=upper, users.guess=users.guess,  
                   optim.control=optim.control,
                   transform=transform,
                   recall = FALSE),
              if (!missing(methods))  list(mle.methods = methods),
              if (!missing(sub.methods)) list(lsq.methods=sub.methods)
              ## "internal" : name should not be changed; should always
                ## be last method!
              ))
  }
  if (RFopt$general$returncall)
    attr(res, "call") <- as.character(deparse(match.call()))
  attr(res, "coord_system") <- .Call(C_GetCoordSystem,
                                     as.integer(MODEL_MLE),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
   return(res)
}
