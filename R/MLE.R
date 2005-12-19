## PrintLevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information

fitvario <-
function(x, y=NULL, z=NULL, T=NULL, data, model, param,
         lower=NULL, upper=NULL, sill=NA,
         ...) {
  
  fitvario.default(x=x, y=y, z=z, T=T, data=data, model=model, param=param,
                lower=lower, upper=upper, sill=sill,
                ...)
}

fitvario.default <-
function(x, y=NULL, z=NULL, T=NULL, data, model, param,
         lower=NULL, upper=NULL, sill=NA,
         trend,
         use.naturalscaling=TRUE,
         PrintLevel=RFparameters()$Print, optim.control=NULL,
         bins=20, nphi=1, ntheta=1, ntime=20,
         distance.factor=0.5,
         upperbound.scale.factor=10,
         lowerbound.scale.factor=20,
         lowerbound.scale.LS.factor=5,
         upperbound.var.factor=10,
         lowerbound.var.factor=100,
         lowerbound.sill=1E-10,
         scale.max.relative.factor=1000,
         minbounddistance=0.001, minboundreldist=0.02,
         approximate.functioncalls=50,
         refine.onborder=TRUE,
         pch=RFparameters()$pch, # no.mle=FALSE,
         var.name="X", time.name="T",
         transform=NULL, standard.style=NULL,
         lsq.methods=c("self", "plain", "sqrt.nr", "sd.inv", "internal"),
         ## "internal" : name should not be changed; should always be last
         ##              method!
         mle.methods=c("ml", "reml"),
         cross.methods=NULL,
 #       cross.methods=c("cross.sq", "cross.abs", "cross.ign", "cross.crps"),
         users.guess=NULL, users.beta=NULL, table.format=FALSE
) {
  ####################################################################
  ##all the following save.* are used for debugging only
  debug <- FALSE
  ##
  ## debug <- TRUE
  save.bins <- bins;
  save.upper <- upper
  save.lower <- lower
  save.optim.control <- optim.control
  save.standard.style <- standard.style
  missing.param <- missing(param)
  missing.trend <- missing(trend)
  trend <- NULL
  OLD.NUGGET.POSITION <- 3

  detailpch <- if (pch=="") "" else "+" 
  show.error.message <- FALSE
 ##
  show.error.message <- TRUE
  save.RFparameters <- RFparameters(no.readonly=TRUE)
  save.options <- options()
  on.exit({RFparameters(save.RFparameters);
           options(save.options)
         })
  if (length(cross.methods) > 0)
    if (is.null(standard.style)) standard.style <- FALSE
    else if (!standard.style)
      stop("if cross.methods then standard.style may not be used")
  ENVIR <- environment()
  save.model <- model
               
######################################################################
###                function definitions                            ###
######################################################################

  ## function not used yet -- for trend definitions
  formulaToMatrix <- function(formula, coord, var.name="X", time.name="T") {
    coord <- as.matrix(coord)
    if (is.null(time.name)) co <- ""
    else co <- paste(time.name,"=coord[,ncol(coord)]")
    for (i in 1:ncol(coord))
      co <- paste(co, ",", var.name, i, "=coord[,", i, "]", sep="")
    eval(parse(text=paste("l <- list(",co,")")))
    fo <- as.list(formula)[[2]]
    variables <- NULL
    while (length(fo)==3) {
      dummy <- as.list(fo)
      stopifnot(dummy[[1]]=="+")
      fo <- dummy[[2]]
      variables <- cbind(eval(dummy[[3]],l), variables)
    }
    variables <- cbind(eval(fo, l), variables) 
    return(variables)
  }


  ## Note
  ## creating "global" variables LSMIN, LSPARAM
  ## using LSQ.WEIGHTS
  LSQsettings <- function(M) {
    assign("LSQ.SELF.WEIGHING", M=="self", envir=ENVIR)
    if (!LSQ.SELF.WEIGHING) {
      assign("LSQ.WEIGHTS", weights[[M]], envir=ENVIR)
      if (varnugNA || zeronugget)
        assign("LSQ.BINNEDSQUARE", sum(binned.variogram^2 * LSQ.WEIGHTS),
               envir=ENVIR)
    }
  }

  LStarget <- function(variab) {
    if (PrintLevel>4) {cat("LSMIN=", LSMIN); print(variab,dig=20)}
    if (any((variab<LSQLB) | (variab>LSQUB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 13.12.03 still happens ...
      if (PrintLevel>1) 
        cat("LSQ WARNING! Forbidden values!", variab, "[", LSQLB, ",",
            LSQUB, "]\n")
      penalty <- variab 
      variab<-pmax(LSQLB, pmin(LSQUB, variab)) 
      penalty <- sum(variab-penalty)^2
      res <- LStarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      return(penalty + res * (1+penalty)) ## not the best ....
    }
    param <- PARAM
    param[LSQINDEX] <- variab
    param <- LSQTRANSFORM(param)

#    str(list("VariogramNatSc", bin.centers, bins, covnr, as.double(param), lpar,
#         logicaldim, xdim,
#         as.integer(length(covnr)), as.integer(pm$anisotropy),         
#         as.integer(pm$op),model.values=double(bins), scalingmethod,
#         PACKAGE="RandomFields", DUP=FALSE))
    
    model.values <-
      .C("VariogramNatSc", bin.centers, bins, covnr, as.double(param), lpar,
         logicaldim, xdim,
         as.integer(length(covnr)), as.integer(pm$anisotropy),         
         as.integer(pm$op),model.values=double(bins), scalingmethod,
         PACKAGE="RandomFields", DUP=FALSE)$model.values
#    print("OK")
    
    if (any(is.na(model.values))) {
      if (PrintLevel>1) {
        cat("\nmissing values!")
        print(variab, dig=20)
        print(model.values)

        print(PARAM)
        print(LSQINDEX)

        param <- PARAM
        param[LSQINDEX] <- variab
        print(param)
        print(LSQTRANSFORM)

        param <- LSQTRANSFORM(param)
        print(param)

        str(.C("VariogramNatSc", bin.centers, bins, covnr, as.double(param), lpar,
               logicaldim, xdim,
               as.integer(length(covnr)), as.integer(pm$anisotropy),         
               as.integer(pm$op),model.values=double(bins), scalingmethod,
               PACKAGE="RandomFields", DUP=FALSE))
       
        print("end model.values")
      }
      return(1E300)
    }
    
    if (LSQ.SELF.WEIGHING) {
      ## weights = 1/ model.values^2
        gx <- binned.variogram / model.values
        gx <- gx[is.finite(gx)]
        if (varnugNA || zeronugget) {
          bgw <- sum(gx^2)
          g2w <- sum(gx)
          ergb <- length(gx) - g2w^2 / bgw
        } else {
          ergb <- sum((gx - 1)^2)
        }
    } else {
      if (varnugNA || zeronugget) {
        ## in this case the calculated variogram model is not the one
        ## to which we should compare, but:
        bgw <- sum(binned.variogram * model.values * LSQ.WEIGHTS)
        g2w <- sum(model.values^2 * LSQ.WEIGHTS)
        ergb <- LSQ.BINNEDSQUARE - bgw^2/g2w
      } else ergb <- sum((binned.variogram - model.values)^2 * LSQ.WEIGHTS)
    }
    
    if (ergb<LSMIN) {
      if (varnugNA) {
        sill <- bgw/g2w
        param[VARIANCE] <- sill * (1.0 - param[NUGGET])
        param[NUGGET] <- sill * param[NUGGET]
      } else if (zeronugget) {
        param[VARIANCE] <- bgw/g2w
      }
      assign("LSMIN",ergb,envir=ENVIR)
      assign("LSPARAM",param,envir=ENVIR)
    }
    return(ergb)
  }


  MLEsettings <- function(M) {
    switch(M,
       "ml" = {
         assign("ML.lcrepet", lc * repet, envir=ENVIR)
         assign("ML.loglcrepet", ML.lcrepet * (1 - log(ML.lcrepet)),
                envir=ENVIR)
         if (!givenCoVariates)
           assign("MLtargetV", data - pm$mean, envir=ENVIR) #else calculated
         ##                   in MLtarget itself
         assign("ML.twopilcrepet", ML.lcrepet * log(2 * pi), envir=ENVIR)
         assign("MLEtarget", MLtarget, envir=ENVIR)
       },
       "reml" = {
         if (givenCoVariates) { ## same as MLE
           assign("REML.lc", nrow(CoVariates) - ncol(CoVariates), envir=ENVIR)
           assign("REML.lcrepet", REML.lc * repet, envir=ENVIR)
               assign("REML.loglcrepet", REML.lcrepet * (1 - log(REML.lcrepet)),
                      envir=ENVIR)
           assign("REML.twopilcrepet", REML.lcrepet * log(2 * pi), envir=ENVIR)
           assign("REML.A", eigen(CoVariates %*% t(CoVariates)), envir=ENVIR)
           assign("REML.A", t(REML.A$vectors[order(abs(REML.A$value))
                                             [1:REML.lc], ]), envir=ENVIR)
           assign("REML.data", crossprod(REML.A, data), envir=ENVIR)
           assign("MLEtarget", REMLtarget, envir=ENVIR)
         } else {
           assign("MLEtarget", REMLtarget, envir=ENVIR) #!
         } 
       },
       stop(paste(M, "unknown"))
      )
  }
    
     
  ## note "global" variables REMLMIN, REMLPARAM
  REMLtarget <- function(variab) {
    if (PrintLevel>4) {print(variab, dig=10)}
    if (any((variab < MLELB) | (variab > MLEUB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 23.12.03 : still happens
      if (PrintLevel>0)
        cat("REML WARNING! Forbidden values! -- if there are too many",
            "warnings try narrower lower and upper bounds for the variables.",
            variab, "[", MLELB, ",", MLEUB, "]\n")
      penalty <- variab 
      variab <- pmax(MLELB, pmin(MLEUB, variab)) 
      penalty <- sum(variab - penalty)^2 ## not the best ....
      res <- REMLtarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      if (PrintLevel>4) {
        cat("penalty", format(c(variab, penalty + res * (1 + penalty))),"\n")
      }
      return(penalty + res * (1+penalty))
    }
    param <- PARAM
    param[MLEINDEX] <- variab
    param <- MLETRANSFORM(param)
    options(show.error.messages = show.error.message)
    cov.matrix <- try(chol(crossprod(REML.A,
                         matrix(.C("CovarianceMatrixNatSc",
                                   distances, lc,
                                   covnr, as.double(param), lpar,
                                   logicaldim,#space time dim of random field
                                   xdim,      #dimension of the distance (vector)
                                   ##      so for isotropic fields it might be 1
                                   ##      although that the random field is 3
                                   as.integer(length(covnr)),
                                   as.integer(pm$anisotropy),
                                   as.integer(pm$op),
                                   cov.matrix=double(lc * lc),
                                   scalingmethod,
                                   PACKAGE="RandomFields", DUP=FALSE)$cov.matrix
                      ,ncol=lc) %*% REML.A)
               ))
    options(show.error.messages = TRUE)
    if (!is.numeric(cov.matrix) || any(is.na(cov.matrix))) {
      if (PrintLevel>1) {
        cat("\nerror in cholesky decomposition -- matrix pos def?")
        print(variab, dig=20)
        str(cov.matrix)
      }
      return(1E300)
    }    
    if (any(diag(cov.matrix)<0)) {stop("chol det<0!")}
    
    ## der faktor zwei erklaert sich wie folgt:
    ## logarithmus der wahrscheinlichkeitsdichte ergibt
    ## - lc/2 log(2pi) - 1/2 log(det C) - 1/2 (D-X m)^T C^{-1} (D-X m)
    ## nun wird zunaechst  (-2) * () betrachtet
    ##     konstante lc/2 log(2pi):
    ##     log(det C) + (D-X m)^T C^{-1} (D-X m),
    ## genauer, fuer wiederholte daten:
    ## repet * log(det C) + sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
    ##
    ## somit hat log(det(C)) vorfaktor 1
    ## nun erhaelt man aus chol(matrix)=wurzel der matrix,
    ## dass sum(log(diag(cov.matrix))) = 1/2 log(det C)
    ##
    ## repet is number of (independent) repeted measurement of data
    ##
    ## bei der Rueckgabe der MLE-werte wird der Faktor und die Konstante
    ## korrigiert, s.u. !!
    logdet <- 2 * sum(log(diag(cov.matrix))) * repet #  repet * log(det C)
    if (!is.finite(logdet)) logdet <- 1E300 ## the best ?! 
    cov.matrix <- chol2inv(cov.matrix, LIN = TRUE) # La.chol2inv, LIN=TRUE
    quadratic <- sum(REML.data * (cov.matrix %*% REML.data))
    ##               sum_i (D_i - Xm)^T C^{-1} (D_i - X m)
    if (varnugNA || zeronugget) {
      ## varnugNA : sill = s^2, var=p * s^2, nugget= (1-p) * s^2
      ##            and s^2 is treated as it was a variance parameter
      ##            and p is an internal parameter stored in NUGGET
      ##
      ## so define C_1 =  s^2 C
      ## Then
      ##    repet * log(det C) + sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##  = repet * log(det C_1) + 2 * repet * lc * log(s) + s^{-2} *
      ##                                      sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##
      ##  Hence the optimal s is given by
      ## s^2 = (2 * repet * lc)^{-1} * 2 * sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##     =  (lc *  repet)^{-1} * sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##
      ## plugging in gives
      ##     repet * log(det C_1) -  repet * lc * log(repet * lc) +
      ##     repet * lc * log(sum_i (D_i-Xm)^T C^{-1} (D_i-X m)) + repet * lc
      ## = logdet + loglcrepet + lcrepet * log(sum_i (D_i-Xm)^T C^{-1} (D_i-X m))
      res <- sum(logdet + REML.loglcrepet + REML.lcrepet * log(quadratic))
    } else {
      res <- sum(logdet + quadratic)
    }
    res <- -0.5 * (res + REML.twopilcrepet)
    if (is.na(res) || is.na(MLEMAX) || debug) {
      filename <- "RandomFields.fitvario.reml.bug.rda"
      txt <- paste("algorithm has failed for unknown reasons -- please contact the author: schlather@cu.lu; please send the data that have been written on the file '", filename, "'.", sep="")
      fitvario.bug <-
        list(x=x, T=T, data=data, model=ctr, param=NULL,
             lower=lower, upper=upper, sill=sill,
             trend=new$trend,
             use.naturalscaling=use.naturalscaling,
             PrintLevel=PrintLevel, optim.control=save.optim.control,
             bins=save.bins, nphi=nphi, ntheta=ntheta, ntime=ntime,
             distance.factor=distance.factor,
             upperbound.scale.factor=upperbound.scale.factor,
             lowerbound.scale.factor=lowerbound.scale.factor,
             lowerbound.scale.LS.factor=lowerbound.scale.LS.factor,
             upperbound.var.factor=upperbound.var.factor,
              lowerbound.var.factor=lowerbound.var.factor,
             lowerbound.sill=lowerbound.sill,
              scale.max.relative.factor=scale.max.relative.factor,
             minbounddistance=minbounddistance,
             minboundreldist=minboundreldist,
             approximate.functioncalls=approximate.functioncalls,
             refine.onborder=refine.onborder,
             pch=pch, var.name=var.name, time.name=time.name,
             transform=MLETRANSFORM, standard.style=standard.style,
             lsq.methods=lsq.methods, cross.methods=cross.methods,
             users.guess=users.guess, users.beta=users.beta,
             table.format=table.format,
             old=list(y=y, z=z, model=model,
               param=if (!missing.param) get("param", envir=ENVIR),
               lower=save.lower,
               upper=save.upper, trend=if (!missing.trend) trend,
               transform=transform, standard.style=save.standard.style),
             )
      if (file.exists(filename) && !debug) {
        cat(txt, "\n")
        stop("file already exists")
      }
      save(file=filename, fitvario.bug, txt, filename, variab, MLELB, MLEUB,
           PARAM, param, MLEINDEX, distances, lc, covnr, lpar, logicaldim, xdim,
           pm, cov.matrix, scalingmethod, REML.A, logdet, quadratic, REML.data,
           varnugNA, zeronugget, REML.twopilcrepet, res, MLEMAX)
      if (is.na(res) || is.na(MLEMAX)) stop(txt)
    }
    if (res > MLEMAX) {
      if (varnugNA) {
        sill <- quadratic / REML.lcrepet
        param[VARIANCE] <- sill * (1.0 - param[NUGGET])
        param[NUGGET] <- sill * param[NUGGET]
      } else {
        if (zeronugget) param[VARIANCE] <- quadratic / REML.lcrepet
      }
      assign("MLEMAX", res, envir=ENVIR)
      assign("MLEPARAM", param, envir=ENVIR)
    }
    if (PrintLevel>3) cat("result REML",res,"\n")
    return(res)
  }

  ## note "global" variables MLEMAX, MLEPARAM
  MLtarget<-function(variab) {
    if (PrintLevel>4) {print(variab, dig=10)}
    if (any((variab < MLELB) | (variab > MLEUB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 23.12.03 : still happens
      if (PrintLevel>0)
        cat("ML WARNING! Forbidden values! -- if there are too many",
            "warnings try narrower lower and upper bounds for the variables.",
            variab, "[", MLELB, ",", MLEUB, "]\n")
      penalty <- variab 
      variab <- pmax(MLELB, pmin(MLEUB, variab)) 
      penalty <- sum(variab - penalty)^2 ## not the best ....
      res <- MLtarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      if (PrintLevel>4) {
        cat("penalty", format(c(variab, penalty + res * (1 + penalty))),"\n")
      }
      return(penalty + res * (1+penalty))
    }
    param <- PARAM
    param[MLEINDEX] <- variab
    param <- MLETRANSFORM(param)
    options(show.error.messages = show.error.message)
    cov.matrix <- try(chol(CC <- matrix(.C("CovarianceMatrixNatSc",
                         distances, lc,
                         covnr, as.double(param), lpar,
                         logicaldim, ##  space time dim of random field
                         xdim,       ##  dimension of the distance (vector)
                         ##          so for isotropic fields it might be 1
                         ##          although that the random field is 3
                         as.integer(length(covnr)),
                         as.integer(pm$anisotropy),
                         as.integer(pm$op),
                         cov.matrix=double(lc * lc),
                         scalingmethod,
                         PACKAGE="RandomFields", DUP=FALSE)$cov.matrix     
                      ,ncol=lc)), silent=!debug)
    options(show.error.messages = TRUE)
    if (!is.numeric(cov.matrix) || any(is.na(cov.matrix))) {
      if (PrintLevel>1) {
        cat("\nerror in cholesky decomposition -- matrix pos def?")
        print(variab, dig=20)
        str(cov.matrix)
      }
      return(1E300)
    }
    if (any(diag(cov.matrix)<0)) {stop("chol det<0!")}
    
    ## der faktor zwei erklaert sich wie folgt:
    ## logarithmus der wahrscheinlichkeitsdichte ergibt
    ## - lc/2 log(2pi) - 1/2 log(det C) - 1/2 (D-X m)^T C^{-1} (D-X m)
    ## nun wird zunaechst  (-2) * () betrachtet:
    ##     konstante lc/2 log(2pi):
    ##     log(det C) + (D-X m)^T C^{-1} (D-X m),
    ## genauer, fuer wiederholte daten:
    ## repet * log(det C) + sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
    ##
    ## somit hat log(det(C)) vorfaktor 1
    ## nun erhaelt man aus chol(matrix)=wurzel der matrix,
    ## dass sum(log(diag(cov.matrix))) = 1/2 log(det C)
    ##
    ## repet is number of (independent) repeted measurement of data
    ##
    ## bei der Rueckgabe der MLE-werte wird der Faktor und die Konstante
    ## korrigiert, s.u. !!
    logdet <- 2 * sum(log(diag(cov.matrix))) * repet #  repet * log(det C)
    if (!is.finite(logdet)) logdet <- 1E300 ## the best ?! 
    cov.matrix <- chol2inv(cov.matrix, LIN = TRUE) # La.chol2inv, LIN=TRUE
    m <- XC <- NULL
    if (givenCoVariates)  {
    ## m minimimiert, falls
      ##      sum_i C^{-1} (D_i - X m) = 0
      ##  <=> \sum_i (X^T C^{-1} D_i) = repet X^T C^-1 X m
      ##  <=> m = repet^{-1} (X^T C^-1 X)^{-1} X^T C^{-1} \sum_i D_i
      ##  <=> m =  (X^T C^-1 X)^{-1} X^T C^{-1} meandata
      XC <- crossprod(CoVariates, cov.matrix) # X^T C^-1
      m <- solve(XC %*% CoVariates, XC %*% meandata)    
      MLtargetV  <- data - as.double(CoVariates %*% m) # (sumD-X m)
    } ## otherwise MLtargetV had been set by the calling function

    quadratic <- sum(MLtargetV * (cov.matrix %*% MLtargetV))
    ##               sum_i (D_i - Xm)^T C^{-1} (D_i - X m)
    if (varnugNA || zeronugget) {
      ## varnugNA : sill = s^2, var=p * s^2, nugget= (1-p) * s^2
      ##            and s^2 is treated as it was a variance parameter
      ##            and p is an internal parameter stored in NUGGET
      ##
      ## so define C_1 =  s^2 C
      ## Then
      ##    repet * log(det C) + sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##  = repet * log(det C_1) + 2 * repet * lc * log(s) + s^{-2} *
      ##                                      sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##
      ##  Hence the optimal s is given by
      ## s^2 = (2 * repet * lc)^{-1} * 2 * sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##     =  (lc *  repet)^{-1} * sum_i (D_i-Xm)^T C^{-1} (D_i-X m)
      ##
      ## plugging in gives
      ##     repet * log(det C_1) -  repet * lc * log(repet * lc) +
      ##     repet * lc * log(sum_i (D_i-Xm)^T C^{-1} (D_i-X m)) + repet * lc
      ## = logdet + loglcrepet + lcrepet * log(sum_i (D_i-Xm)^T C^{-1} (D_i-X m))
      res <- sum(logdet + ML.loglcrepet + ML.lcrepet * log(quadratic))
    } else {
      res <- sum(logdet + quadratic)
    }
    res <- -0.5 * (res + ML.twopilcrepet)
    if (is.na(res) || is.na(MLEMAX) || debug) {
      filename <- "RandomFields.fitvario.bug.rda"
      txt <- paste("algorithm has failed for unknown reasons -- please contact the author: schlather@cu.lu; please send the data that have been written on the file '", filename, "'.", sep="")
      fitvario.bug <-
        list(x=x, T=T, data=data, model=ctr, param=NULL,
             lower=lower, upper=upper, sill=sill,
             trend=new$trend,
             use.naturalscaling=use.naturalscaling,
             PrintLevel=PrintLevel, optim.control=save.optim.control,
             bins=save.bins, nphi=nphi, ntheta=ntheta, ntime=ntime,
             distance.factor=distance.factor,
             upperbound.scale.factor=upperbound.scale.factor,
             lowerbound.scale.factor=lowerbound.scale.factor,
             lowerbound.scale.LS.factor=lowerbound.scale.LS.factor,
             upperbound.var.factor=upperbound.var.factor,
              lowerbound.var.factor=lowerbound.var.factor,
             lowerbound.sill=lowerbound.sill,
              scale.max.relative.factor=scale.max.relative.factor,
             minbounddistance=minbounddistance,
             minboundreldist=minboundreldist,
             approximate.functioncalls=approximate.functioncalls,
             refine.onborder=refine.onborder,
             pch=pch, var.name=var.name, time.name=time.name,
             transform=MLETRANSFORM, standard.style=standard.style,
             lsq.methods=lsq.methods, cross.methods=cross.methods,
             users.guess=users.guess, users.beta=users.beta,
             table.format=table.format,
             old=list(y=y, z=z, model=model,
               param=if (!missing.param) get("param", envir=ENVIR),
               lower=save.lower,
               upper=save.upper, trend=if (!missing.trend) trend,
               transform=transform, standard.style=save.standard.style),
             )
      if (file.exists(filename) && !debug) {
        cat(txt, "\n")
        stop("file already exists")
      }
      save(file=filename, fitvario.bug, txt, filename, variab, MLELB, MLEUB,
           PARAM, param, MLEINDEX, MLETRANSFORM, distances, lc, covnr,
           lpar, logicaldim, xdim,
           pm, cov.matrix, scalingmethod, repet, logdet, quadratic,
           givenCoVariates, XC, m, meandata, MLtargetV, ML.loglcrepet,
           ML.lcrepet, ML.twopilcrepet, CC,
           varnugNA, zeronugget, res, MLEMAX)
      if (is.na(res) || is.na(MLEMAX)) stop(txt)
    }
    if (res > MLEMAX) {
      if (varnugNA) {
        sill <- quadratic / ML.lcrepet
        param[VARIANCE] <- sill * (1.0 - param[NUGGET])
        param[NUGGET] <- sill * param[NUGGET]
      } else{
        if (zeronugget) param[VARIANCE] <- quadratic / ML.lcrepet
      }
      assign("MLEMAX", res, envir=ENVIR)
      assign("MLEPARAM", param, envir=ENVIR)
    }
    if (PrintLevel>5) cat("result MLE", res, "\n")
    return(res)
  }

  
     
  crosssettings <- function(M) {
    stopifnot(is.null(trend))
    assign("CROSS.KRIGE", if (is.na(pm$mean)) "O" else "S",# ordinary/simple
           envir=ENVIR)
    switch(M,
           "cross.sq" = {
              assign("CROSS.VAR", FALSE, envir=ENVIR)
              ## d can have more than 1 element
              CROSS.DIST <- function(x, d) sum((x-d)^2)
              environment(CROSS.DIST) <- NULL
            },
           "cross.abs" = {
              assign("CROSS.VAR", FALSE, envir=ENVIR)
              ## d can have more than 1 element
              CROSS.DIST <- function(x, d) sum(abs(x-d))
              environment(CROSS.DIST) <- NULL
           },
           "cross.ign" = {
             assign("CROSS.VAR", TRUE, envir=ENVIR)
             CROSS.DIST <- function(x, d) {
               sum(0.5 * log(2 * pi * x$var) + (d - x$est)^2/ (2 * x$var))
             }
             environment(CROSS.DIST) <- NULL
           },
           "cross.crps" = {
             assign("CROSS.VAR", TRUE, envir=ENVIR)
             CROSS.DIST <- 
               function(x, d) {
                 stopifnot(x$var>0)
                 sigma <- sqrt(x$var)
                 n <- (d - x$est) / sigma
                 sigma * sum(n * (2 * pnorm(n) - 1) + 2 * dnorm(n) - 1/sqrt(pi))
               }
             environment(CROSS.DIST) <- .GlobalEnv
           },
           stop(paste(M, "unknown"))
           )
    assign("CROSS.DIST", CROSS.DIST, envir=ENVIR)
  }
 
  crosstarget <- function(variab) {
    if (PrintLevel>4) {print(variab, dig=10)}
    if (any((variab<CROSSLB) | (variab>CROSSUB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 23.12.03 : still happens
      if (PrintLevel>0)
        cat("CROSS WARNING! forbidden variab values ! -- if there are too many warnings try narrower lower and upper bounds for the variables.",variab,"[",CROSSLB,",",CROSSUB,"]\n")
      penalty <- variab 
      variab <- pmax(CROSSLB, pmin(CROSSUB, variab)) 
      penalty <- sum(variab-penalty)^2 ## not the best ....
      res <- crosstarget(variab)
      if (res<=0) return(penalty + res * (1-penalty))
      if (PrintLevel>4) {
        cat("penalty",format(c(variab,penalty + res * (1+penalty))),"\n")
      }
      return(penalty + res * (1+penalty))
    }
    param <- PARAM
    model <- pm
    
    if (givenCoVariates) {
      stopifnot(is.null(trend))
      model$mean <- variab[nCROSSINDEX + 1:nCoVariates]
      model$trend <- NULL
      stopifnot(is.null(trend)) ## not programmed yet,
      ##                           currently only nCoVariates==1
      param[CROSSINDEX] <- variab[1:nCROSSINDEX]
    } else param[CROSSINDEX] <- variab
    model$param <- CROSSTRANSFORM(param)

    res <- 0.0
    ## data is matrix!
    for (d in 1:lc) {
 #if (PrintLevel>3) cat(d, "")   
      res <- res +
        CROSS.DIST(Kriging(CROSS.KRIGE, x=new$x[d, , drop=FALSE], grid=FALSE,
                           model=model, given=new$x[-d, , drop=FALSE],
                           data=data[-d, , drop=FALSE], trend = trend, pch="",
                           return.variance=CROSS.VAR, internal = TRUE),
                   data[d, ])
    }
    res <- res / CROSS.lcrepet
    if (res<CROSSMIN) {
      assign("CROSSMIN", res, envir=ENVIR)
      assign("CROSSMODEL", model, envir=ENVIR)
    }
    if (PrintLevel>6) cat("result cross val", res, "\n")
    return(res)
  }
  
  lsq.covariates <- function(modelparam) {
     if (givenCoVariates) {
      if (!is.null(pm$trend)) stop("sorry, not implemented yet")
      XC <- solve(matrix(.C("CovarianceMatrixNatSc",
                            distances, lc,
                            covnr, as.double(modelparam), lpar,
                            logicaldim,#space time dim of random field
                            xdim,      #dimension of the distance (vector)
                            ##      so for isotropic fields it might be 1
                            ##      although that the random field is 3
                            as.integer(length(covnr)),
                            as.integer(pm$anisotropy),
                            as.integer(pm$op),
                            cov.matrix=double(lc * lc),
                            scalingmethod,
                            PACKAGE="RandomFields", DUP=FALSE)$cov.matrix
                         ,ncol=lc), CoVariates)
      solve(crossprod(XC, CoVariates), crossprod(XC, meandata))
    } else {
      stopifnot(!is.na(pm$mean))
      pm$mean
    }
  }

show <- function(nr, M, OPT, PARAM)
  cat("\n ", M, ", ", switch(nr, "start", "grid ", "re-do"), ": value=",
      format(OPT, dig=6), ", param=", format(PARAM, dig=2), sep="")

######################################################################
###                Initial settings                                ###
######################################################################
  if (PrintLevel>4) cat("\ninitial settings...")
  new <- CheckXT(x, y, z, T, grid=FALSE)
  new$x <- as.matrix(new$x)
  
###########################      transform     #######################
  ## in case functions are given within model, they are extracted and
  ## replaced NaN; `transform' is set instead
  ## must be the first one, since model cannot be Prepare-d
  
  users.transform <- transform 
  if (is.list(model) && is.list(unlist(model))) {
    if (!any(sapply(unlist(model), is.function)))
      stop("model does not seem to be defined correctly")
    if (!is.null(users.transform)) 
      stop("if functions are given in `model', `transform' must be NULL") 
    users.transform <- function(param) {
      p <- pm
      p$param <- param
      p <- convert.to.readable(p, allow="list")$model
      for (i in 1:length(p.transform))
        param[p.transform[i]] <- f.transform[[i]](p)
      return(param)
    }
    functionToNaN <- function(p) {
      ## replaces functions by NaN
      if (is.list(p)) for (i in 1:length(p)) p[[i]] <- functionToNaN(p[[i]])
      else if (is.function(p)) {
        assign("f.transform", c(f.transform, p), envir=ENVIR)
        return(NaN)
      }
      else if (is.nan(p))
        stop("NaN are not allowed if functions are given in `model'")
      return(p)
    }
    functionNum <- function(p) {
      ## replaces functions by an enumeration
      if (is.list(p)) for (i in 1:length(p)) p[[i]] <- functionNum(p[[i]])
      else if (is.function(p)) {
        assign("zaehler", zaehler + 1, envir=ENVIR)
        return(zaehler)
      }
      return(if (is.numeric(p)) rep(0, length(p)) else p) # 100000 > anzahl
      ##                                                             parameter
    }
    ## aniso could have been a list, this is corrected by the following command
    f.transform <- NULL
    model.nan <- functionToNaN(model)
    model.nan <-
      convert.to.readable(PrepareModel(model.nan,NULL, # new$spacedim+new$Time
                                       ), allow="list")
    if (any(sapply(model.nan[-1], function(x) !is.null(x) && is.nan(x))))
      stop("transform functions may not given (yet) for mean and trend")
    zaehler <- -10000
    ## PrepareModel kann nicht direkt aufgerufen werden, da im Modell
    ## Funktionen vorhanden sein koennen!
    ## fct.num hat negative werte, wo funktionen waren;
    ##         0 (oder 1.0 fuer fehlender scale bei nugget modell
    ##         an Stellen mit numerischen Werten
    fct.num <- PrepareModel(functionNum(model))$param
    ## p.transform enthaelt zunaechst alle param-Stellen
    p.transform <- PrepareModel(parampositions(model.nan, print=FALSE))$param
    stopifnot(length(fct.num)==length(p.transform)) ## if so, programming error!
    ## order schiebt die Positionen aller negativen Werte nach vorne
    ## somit enthaelt p.transform geordnet die Positionen, wo bisher
    ## Funktionen standen.
    p.transform <-
      p.transform[order(fct.num)][1:sum(fct.num < 0, na.rm=TRUE)]
    model <- model.nan
    rm("model.nan")
  }


############################    model       ##########################
  stopifnot(all(is.finite(as.matrix(data))))
  data <- as.vector(data)
  spacedim <- new$spacedim
  logicaldim <- as.integer(spacedim + new$Time)
  pm <- PrepareModel(model, param, logicaldim, trend)
  ctr <- convert.to.readable(pm)
  if (!is.null(ctr$param))
    ## to make sure that the model is  given in a standardised way
    ## if not truely a composed covariance function
    pm <- PrepareModel(ctr$model, ctr$param, new$spacedim+new$Time, trend)
  covnr <- as.integer(pm$covnr)
  lpar <- as.integer(length(pm$param))  

  
##############              Coordinates               #################
  if (PrintLevel>4) cat("\ncoordinates...")
  if (pm$anisotropy) {
    xdim <- as.integer(logicaldim)
    dd <- 0
    
    coord <- new$x
    if (!is.null(T)) {
      y <- expand.grid(1:nrow(new$x), seq(new$T[1], new$T[2], new$T[3]))
      coord <- cbind(coord[y[,1], ], y[,2])
    }
    
    lc <- nrow(coord)
    ## note: the direct C call needs matrix where points are given column-wise
    ##       whereas the R function CovarianceFct need them row-wise,
    ##                   except for fctcall==CovarianceMatrix
    distances <- matrix(nrow=logicaldim, ncol=lc * (lc-1) / 2)
    ## distances are used to calculate bounds for the scale parameter(s)
    ## to be estimated -- further it is used in MLtarget

    for (i in 1:logicaldim) {
      dx <- outer(coord[,i], coord[,i], "-")
      distances[i, ] <- as.double(dx[lower.tri(dx)])
      ## x <- c(1,6,0.1)
      ## outer(x, x, "-")
      ##     [,1] [,2] [,3]
      ##[1,]  0.0 -5.0  0.9
      ##[2,]  5.0  0.0  5.9
      ##[3,] -0.9 -5.9  0.0
      if (i<=spacedim) dd <- dd + distances[i, ]^2
    }
    dd <- sqrt(dd)
    coord <- NULL
   
    maxdistances <- max(dd)
   ## mindistances <- min(dd[dd!=0])
    mindistances <- min(dd)
    rm("dd")
  } else {
    lc <- nrow(new$x)
    xdim <- as.integer(1)
    distance.matrix <- dist(new$x)
    distances <- as.double(distance.matrix)
    distance.matrix <- as.matrix(distance.matrix)
    mindistances <- min(distances[lower.tri(distances)])
    maxdistances <- max(distances)
  }
  if (mindistances==0) stop("two points given twice")
  gc() ## necessairy, since otherwise too much memory is being taken because
  ##      of dd[dd!=0]
  
  storage.mode(lc) <- "integer"


##########   Check consistency of NA, NaN and users.transform  ##########
  if (PrintLevel>4) cat("\nconsistency of NA, NaN, transform...")
  PARAM <- pm$param
  ##PARAM <- param: just to make clear in MLtarget and LStarget what is global,
  ##and do not overwrite param
  index <- is.na(PARAM) & !is.nan(PARAM)
  nindex <- sum(index)
  if (!is.null(users.transform)) {
    pp <- para <- runif(length(PARAM))
    para[is.nan(PARAM)] <- NaN
    para <- users.transform(para)
    if (length(para)!=length(pp))
      stop("transform does not preserve length of the parameter vector")
    txt <- ""
    if (any(is.nan(para)))
     txt <-
       paste("position(s)", paste(which(is.nan(para)), collapse=","),
             "of param still NaN after application of 'transform'.\n'transform'",
             "is too complex or does not fill the positions with value NaN.\n\n")
    if (any(pp <- para[!is.nan(PARAM)]!=pp[!is.nan(PARAM)])) {
      pos <- which(!is.nan(PARAM))
      if (options()$warn<2) {
        if (any(pq <- is.na(PARAM)[!is.nan(PARAM)] & pp))
          warning(paste("position(s)", paste(pos[pq], collapse=","),
                        "of param has been changed by 'transform'",
                        "that were NA.\n"))
        pp <- pp & !pq
      }
      if (any(pp))
        txt <- paste(txt,"position(s)", paste(pos[pp], collapse=","),
                     "of param has been changed although not being NaN.\n\n")
    }
    if (txt!="") stop(txt)
  }

 
##############      find upper and lower bounds      #################  
  if (PrintLevel>4) cat("\nbounds...")
  UnList <- function(l) {
    l <- unlist(l)
    paste(names(l), l, sep="=", collapse="|")
  }
  txt <- "lower and upper are both lists or vectors of the same length or NULL"
  lengthmismatch <- "lengths of bound vectors do not match model"
  structuremismatch <- "structures of bounds do not match the model"
  parscale <- autostart <- convert.to.readable(pm, allowed="list")$model
  ## parscale gives the magnitude of the parameters to be estimated
  ##     passed to optim/optimise so that the optimiser estimates
  ##     values around 1
  ## autostart will give the starting values for LSQ

  ## press lower and upper into a standard format...
  if (is.list(lower)) {
    if (!is.list(upper)) stop(txt)
    dummy <- pm
    dummy$param <- rep(NA, length(dummy$param))
    dummy <- Unlist(convert.to.readable(dummy, allowed="list"))
    if (!is.numeric(lower[[i]]$v)) lower[[i]]$var <- NA ##i.e. not given
    if (pm$anisotropy)
      if (is.null(lower[[i]]$a)) lower[[i]]$aniso <- autostart$a * NA
      else if (!is.numeric(lower[[i]]$s)) lower[[i]]$scale <- NA
    if (!is.numeric(lower[[i]]$k)) lower[[i]]$kappas <- autostart$param * NA
    bndtst <- PrepareModel(lower, time=new$spacedim+new$Time, trend=trend)
    bndtst$param <- rep(NA, length(bndtst$param))
    if (dummy!=Unlist(convert.to.readable(bndtst, allowed="list")))
      stop(structuremismatch)
    lower <- convert.to.readable(bndtst, allowed="list")

    if (!is.numeric(upper[[i]]$v)) upper[[i]]$var <- NA
    if (pm$anisotropy)
      if (is.null(upper[[i]]$a)) upper[[i]]$aniso <- autostart$a * NA
      else if (!is.numeric(upper[[i]]$s)) upper[[i]]$scale <- NA
    if (!is.numeric(upper[[i]]$k)) upper[[i]]$kappas <- autostart$param * NA
    bndtst <- PrepareModel(upper, time=new$spacedim+new$Time, trend=trend)
    bndtst$param <- rep(NA, length(bndtst$param))
    if (dummy!=Unlist(convert.to.readable(bndtst, allowed="list")))
      stop(structuremismatch)
    upper <- convert.to.readable(bndtst, allowed="list")
    rm("dummy")
  } else if (is.vector(lower)) {
    if (!is.vector(upper) || length(upper)!=length(lower)) stop(txt)
    if (is.vector(ctr$param)) {
      nulltrend <- as.integer(missing(trend) || is.null(trend))
      if (length(lower) < length(ctr$param)) {
        if (length(ctr$param) - 3 - nulltrend != length(lower))
          stop("length of lower does not match param")
        nugg <- if (is.finite(ctr$param[OLD.NUGGET.POSITION]) &&
                    ctr$param[OLD.NUGGET.POSITION]==0) 0 else NA
        lower <- c(rep(NA, 1 + nulltrend), nugg, NA, lower)
        upper <- c(rep(NA, 1 + nulltrend), nugg, NA, upper)
      } else if (length(ctr$param) < length(lower)) stop(lengthmismatch)
      lower <- convert.to.readable(PrepareModel(ctr$model, lower,
                                                new$spacedim+new$Time, trend),
                                   allowed="list")
      upper <- convert.to.readable(PrepareModel(ctr$model, upper,
                                                new$spacedim+new$Time, trend),
                                   allowed="list")
    } else if (is.matrix(ctr$param)) {
      if (nrow(lower)<length(ctr$param)) {
        if (length(ctr$param) != 2 + nrow(lower)) stop(lengthmismatch)
        lower <- c(NA, NA, lower)
        upper <- c(NA, NA, upper)
      } else if (nrow(lower)>length(ctr$param)) stop(lengthmismatch)     
      lower <- convert.to.readable(PrepareModel(ctr$model,
                matrix(lower, ncol=ncol(ctr$param), nrow=nrow(ctr$param)),
                                                new$spacedim+new$Time, trend),
                                                allowed="list")
      upper <- convert.to.readable(PrepareModel(ctr$model,
               matrix(upper, ncol=ncol(ctr$param), nrow=nrow(ctr$param)),
                                                new$spacedim+new$Time, trend),
                                   allowed="list")
    } else stop("vector valued lower bound only allowed if param is given")
  } else { # !(is.vector(lower)) i.e. not given
    if (!is.null(upper) || !is.null(lower)) stop(txt)
    lower <- pm
    lower$param <- rep(NA, length(lower$param))
    lower <- upper <- convert.to.readable(lower, allowed="list")
  }
  lower <- lower$model
  upper <- upper$model

### now the bounds and starting values are set...
  scale.pos <- pm
  scale.pos$param <- rep(0, length(scale.pos$param))
  scale.pos <- convert.to.readable(scale.pos, allowed="list")$model
  ## scale.pos will give the positions where scale parameters can be found
  ##    this is necessary to change the limits for the scale parameters
  ##    before MLE is performed
  vardata <- var(data)
  truedim <- new$spacedim+new$Time ## could made better by considering
  ##   the anisotropy matrices if there are any
  ##   but too much work with the low profit too allow some more
  ##   models for the estimation (those that are allowed for
  ##   lower dimensions only, and the anisotropy matrix has a rank
  ##   small enough

  for (i in seq(1, length(lower), 2)) {   
    if (is.na(lower[[i]]$var))
      ## lower bound of first model is treated differently!
      ## so the "main" model should be given first!
      ## usually "main" model + (estimated) nugget
      lower[[i]]$var <- if (i==1) vardata / lowerbound.var.factor else 0
    if (is.na(upper[[i]]$var)) upper[[i]]$var <- upperbound.var.factor *vardata
    if (is.na(autostart[[i]]$var))
      parscale[[i]]$var <- autostart[[i]]$var <- vardata / (length(lower) + 1) *2
    if (parscale[[i]]$var==0) parscale[[i]]$var <- 1
    if (pm$anisotropy) {
      scale.pos[[i]]$aniso <- rep((i+1)/2, length(scale.pos[[i]]$aniso))      
      diagonal <- as.vector(diag(logicaldim))
      if (upper[[i]]$model=="nugget") {
        ## user may have given an anisotropy matrix
        ## matrix entries are supposed to be within [-1,1], so
        ## -10, 10 as bounds is sufficient
        upper[[i]]$aniso[is.na(upper[[i]]$aniso)] <- 10
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & diagonal] <- 0
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & !diagonal] <- -10
      } else {
        upper[[i]]$aniso[is.na(upper[[i]]$aniso)] <-
          lowerbound.scale.LS.factor / mindistances
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & diagonal] <-
          1 / (upperbound.scale.factor * maxdistances)
        lower[[i]]$aniso[is.na(lower[[i]]$aniso) & !diagonal] <-
          -lowerbound.scale.LS.factor / mindistances
      }
      idx <- is.na(autostart[[i]]$aniso)
      parscale[[i]]$aniso[idx] <- 8 / (maxdistances + 7 * mindistances)
      autostart[[i]]$aniso[idx & diagonal] <- 8 / (maxdistances + 7 * mindistances)
      autostart[[i]]$aniso[idx & !diagonal] <- 0      
    } else { ## isotropy
      scale.pos[[i]]$scale <- (i+1)/2
      if (upper[[i]]$model=="nugget") {
        ## user should never have chosen the value for the scale 
        lower[[i]]$scale <- 0
        upper[[i]]$scale <- 1 
      } else {
        if (is.na(lower[[i]]$scale))
          lower[[i]]$scale <- mindistances / lowerbound.scale.LS.factor
        if (is.na(upper[[i]]$scale))
          upper[[i]]$scale <- maxdistances * upperbound.scale.factor
      }
      if (is.na(autostart[[i]]$scale))
        parscale[[i]]$scale <- autostart[[i]]$scale <-
          (maxdistances + 7 * mindistances) / 8
    }
    kappas <- parameter.range(autostart[[i]]$model, truedim)
    if (!is.null(kappas)) {
      if (length(kappas)==1 && is.nan(kappas)) stop(paste("model #", (i+1)/2,
                 "is not allowed for the considered dimension"))
      fixed.kappas <- kappas$th[[1]][2,] == kappas$th[[1]][1,]
      if (any(fixed.kappas)) {
        if (!all(is.finite(autostart[[i]]$kappa[fixed.kappas])))
          stop(paste("in model #", (i+1)/2,
                     "the parameters that select subclasses are not all given"))
        for (nk in 1:length(kappas$th)) 
          if (all(kappas$th[[nk]][fixed.kappas] ==
                  autostart[[nk]]$kappa[fixed.kappas])) break
        if (any(kappas$th[[nk]][fixed.kappas]!=autostart[[nk]]$kappa[fixed.kappas]))
          stop(paste("Could not find the indicated subclass for model #",
                     (i+1)/2))
      } else nk <- 1
      lower[[i]]$kappas[is.na(lower[[i]]$kappas)] <-
        kappas$pr[[nk]][1, is.na(lower[[i]]$kappas)]
      upper[[i]]$kappas[is.na(upper[[i]]$kappas)] <-
        kappas$pr[[nk]][2, is.na(upper[[i]]$kappas)] 
      autostart[[i]]$kappas <- (upper[[i]]$kappas + lower[[i]]$kappas) / 2
      parscale[[i]]$kappas <- rep(1, length(autostart[[i]]$kappas))
    } # !is.null(kappas)
  }

  lower <- PrepareModel(model=lower, time=new$spacedim+new$Time, trend=trend,
                        named=PrintLevel>1)$par
  upper <- PrepareModel(model=upper, time=new$spacedim+new$Time, trend=trend,
                        named=PrintLevel>1)$par
  lower[is.nan(PARAM)] <- -Inf ## nan should not be checked -- this is
  ##                              completely up to the user ...
  upper[is.nan(PARAM)] <- Inf
  
  autostart <- PrepareModel(model=autostart, time=new$spacedim+new$Time, trend=trend)$par
  parscale <- 
    PrepareModel(model=parscale, time=new$spacedim+new$Time, trend=trend)$par
  scale.pos <-
    PrepareModel(model=scale.pos, time=new$spacedim+new$Time, trend=trend)$par
  
  
###################      check optim.control      ####################
  stopifnot(all(is.finite(parscale[index])))
  if (length(optim.control)>0) {
    stopifnot(is.list(optim.control))
    forbidden.param <- c("parscale", "fnscale")
    if (any(!is.na(pmatch(names(optim.control), forbidden.param))))
      stop(paste(paste(forbidden.param, collapse=" or "),
                 "may not be given in the optim.contol list"))
    ## fnscale=-1 turns the problem into a maximisation problem
  } else optim.control <- list()
  

##########################       data         ########################
  ## coord instead of new$x because of possible time component
  
  if (length(data)==lc) { # data is as.double(data)
    repet <- 1
    meandata <- data
  } else {
    if (as.integer(length(data) / lc) * lc != length(data))
      stop("length of data not a multiple of number of locations")
    data <- matrix(data, nrow=lc)
    repet <- ncol(data)
    meandata <- (data %*% rep(1, repet)) / repet
  }
  
  
##############           Covariates                   #################
  if (PrintLevel>4) cat("\nCoVariates...")
  givenCoVariates <- TRUE
  nCoVariates <- 1
 if (is.null(pm$trend)) {
    if (is.na(pm$mean)) {
      CoVariates <- matrix(1, nrow=nrow(new$x), ncol=1)
      stopifnot(!is.null(CoVariates)) 
    }
    else {
      givenCoVariates <- FALSE 
    }
  } else {
    stop("sorry. not available yet") ## programmed but not tested yet
    if (is.matrix(trend)) {
      stopifnot(nrow(trend)==nrow(new$x))
      CoVariates <- trend
    } else { ## formula
      CoVariates <- formulaToMatrix(trend, cbind(new$x,new$T),
                                 var.name=var.name, time.name=time.name)
    }
    nCoVariates <- ncol(CoVariates)
  }


######################################################################
###                                                                ###
###   check which parameters are NA -- only old style is allowed   ###
###                                                                ###
###     certain combinations of NA allow for faster algorithms     ###
###                                                                ###
###     !is.na(sill) needs special treatment, hence must be        ###
###     identified                                                 ###
###                                                                ###
###     scaling method must be identified                          ###
###                                                                ###
###     some autostart values are calculated                       ###
###                                                                ###
###                                                                ###
######################################################################
  varnugNA <- ## is.na(nugget) and is.na(var)
    zeronugget <- ## nugget==0
      sillbounded <- ## !is.na(sill)
        FALSE
  scalingmethod <- as.integer(0)
  if (is.null(standard.style)) ## useful for lsq and mle methods
    standard.style <- is.vector(ctr$param) && is.null(users.transform)
  else if (!is.vector(ctr$param) && standard.style) {
    standard.style <- FALSE
    warning("standard.style must be FALSE for the given model specification.")
  }
  if (standard.style) {
    if (PrintLevel>4) cat("\nstandard style settings...")
    VARIANCE <- 1
    KAPPA <- 2
    n.kappas <- .C("GetNrParameters", covnr, as.integer(1),
                   as.integer(logicaldim), k=integer(1),
                   PACKAGE="RandomFields", DUP=FALSE)$k
    SCALE  <- KAPPA + n.kappas
    if (SCALE<length(PARAM)) {
      NUGGET <- SCALE + 1
      nugget <- PARAM[NUGGET]
    } else {
      NUGGET <- 0  ## nugget effect is zero, so not included in PARAM
      nugget <- 0.0
    }
    if ((covnr==.C("GetModelNr", as.character("nugget"), as.integer(1),
           nr=integer(1), PACKAGE="RandomFields")$nr) &&
        (any(index[-NUGGET]) || (index[VARIANCE]!=0.0))) #works also if NUGG==0??
      stop("if model=nugget effect then the variance must be zero, and only mean and/or nugget can be NA")   
    
    if (sillbounded <- !is.na(sill)) {
      ## only VARIANCE need to be optimised
      ## NUGGET = SILL - VARIANCE
      if (xor(is.na(PARAM[VARIANCE]), is.na(nugget)))
        stop("Sill fixed. Then variance and nugget should be given or unknown simultaneously.")
      if ((is.nan(PARAM[VARIANCE]) || is.nan(nugget)) &&
          !is.null(users.transform))
        stop("not sure what to do now: sill fixed, but some transformation is given -- try standard.style=FALSE")
      if (!is.na(PARAM[VARIANCE]) && (PARAM[VARIANCE] + nugget!=sill)) {
        cat("var=", PARAM[VARIANCE] , "nug=", nugget, "sill=", sill,
            "sill-var-nug=", sill - PARAM[VARIANCE] - nugget,"\n")
        stop("sill !=  variance + nugget")
      }
      autostart[VARIANCE] <- parscale[VARIANCE] <- sill/2
      upper[NUGGET] <- upper[VARIANCE] <- sill
      lower[NUGGET] <- lower[VARIANCE] <-  0
    } else { ## not sill bounded
      if (is.na(PARAM[VARIANCE]) &&
          (!is.nan(PARAM[VARIANCE]) || is.null(users.transform))) {
        if (is.na(nugget) && (!is.nan(nugget) || is.null(users.transform))) {
          ## of interest currently 
          ## both NUGGET and VARIANCE have to be optimised
          ## instead optimise the SILL (which can be done by
          ## a formula) and estimate the part of NUGGET
          ## (which will be stored in NUGGET) !
          ##  consequences:
          ## * estimtation space increased only by 1, not 2
          ## * real nugget: NUGGET * SILL
          ## *     variance: (1-NUGGET) * SILL
          varnugNA <- TRUE
          autostart[NUGGET] <- autostart[VARIANCE] <- vardata / 2
          parscale[NUGGET] <- 0.5
          lower[VARIANCE] <- 0
        } else { ## not sillbounded, is.na(variance), !is.na(nugget)
          if (nugget==0) {
            ## here SILL==VARIANCE and therefore, variance
            ## can be estimated by a formula without increasing the dimension
            
            ## more than only the variance has to be estimated
            zeronugget <- TRUE
            lower[VARIANCE] <- PARAM[VARIANCE] <- 0 ## dito
            autostart[VARIANCE] <- vardata
          } else { ## not sillbounded, is.na(variance), nugget!=0
            lower[VARIANCE] <- pmax(0, (vardata-nugget)/lowerbound.var.factor)
            if (lower[VARIANCE] < lowerbound.sill) {
              if (PrintLevel>1)
                cat("low.var=",lower[VARIANCE]," low.sill",lowerbound.sill,
                    "\ estimated variance from data=", vardata,
                    "nugget=", nugget, "\n")
              warning("param[NUGGET] might not be given correctly.")
              lower[VARIANCE] <- lowerbound.sill
            }
            parscale[VARIANCE] <- sqrt(lower[VARIANCE] * vardata)
            autostart[VARIANCE]<- lower[VARIANCE]
          }
        }
      } else { ## not sillbounded, !is.na(variance)
        if (PARAM[VARIANCE]==0.0) {
          if (any(is.na(PARAM[KAPPA:SCALE])))
            stop("If variance=0, estimating scale or model parameters does not make sense")
          lower[VARIANCE] <- 0
        }
        if (is.na(nugget)) { ## and !is.na(param[VARIANCE])
          lower[NUGGET] <-
            pmax(0, (vardata - PARAM[VARIANCE]) / lowerbound.var.factor)
          if (lower[NUGGET] < lowerbound.sill) {
            if (PrintLevel>1)
              cat("\nlower nugget bound=", lower[NUGGET],
                  " < lower sill bound=", lowerbound.sill,
                  " -- is the variance given correctly?\n",sep="")
            ## warning("Has param[VARIANCE] been given correctly?!")
            lower[NUGGET] <- lowerbound.sill
          }
          autostart[NUGGET] <- lower[NUGGET]
          parscale[NUGGET] <- sqrt(lower[VARIANCE] * vardata)
        }
      } ##  else { ## not sillbounded, !is.na(variance)
    } ## else { ## not sill bounded
    
    if (index[SCALE] && use.naturalscaling) {
      ##  11: exact or numeric
      ##  13: MLE or numeric
      scalingmethod <- as.integer(11)  ## check with RFsimu.h, NATSCALEMLE!
      ##          or   3 ?
      if (.C("GetNaturalScaling", covnr, as.double(autostart[-1:-(KAPPA-1)]),
             scalingmethod, double(1), error=integer(1),
             PACKAGE="RandomFields", DUP=FALSE)$error)
        {
          scalingmethod <- as.integer(0)
          if (PrintLevel>1) cat("No natural scaling.")
        }
      ## used in cross validation method only via kriging
    }
    
    stopifnot(length(index)==length(autostart)) ## simple check
    notidx <- !index & !is.nan(PARAM)
    incons.idx <- PARAM[notidx] > upper[notidx] | PARAM[notidx] < lower[notidx]
    if (any(incons.idx)){
      inconsist <- rep(NA, length(notidx))
      inconsist[notidx] <- incons.idx
      if (PrintLevel>1) {
        cat("\ninconsistent boundaries\n")
        print(rbind("fixed"=notidx, "inconsist?"=inconsist, lower, PARAM, upper))
      }
      warning("fixed parameters out of range\n")
    }
    stopifnot( varnugNA + zeronugget + sillbounded <= 1)
  }
  RFparameters(PracticalRange=scalingmethod)


  
######################################################################
###                     Estimation part itself                     ###
######################################################################

###################  preparation  ################
  ## methods
  formals <- formals()
  allprimmeth <- c("autostart", "users.guess")
  nlsqinternal <- 3 ## cross checked after definition of weights below
  lsq.orig.methods <- eval(formals$lsq.methods)
  alllsqmeth <- c(lsq.orig.methods[-length(lsq.orig.methods)],
                  paste("internal", 1:nlsqinternal, sep=""))
  allmlemeth <- eval(formals$mle.methods)
  allcrossmeth <- eval(formals$cross.methods)
  allmethods <- c(allprimmeth, alllsqmeth, allmlemeth, allcrossmeth)

  ## how preceding methods have been condidered ?
  ## note cm is used again at the very end when error checking
  cm <- cumsum(c(0, length(allprimmeth), length(alllsqmeth),
                     length(allmlemeth), length(allcrossmeth)))
  cm <- cbind(cm[-length(cm)] + 1, cm[-1])
  cm <- apply(cm, 1, function(x) x[1] : x[2])
  names(cm) <- c("prim", "lsq", "mle", "cross")
  methodprevto <- list(lsq=c(cm$prim),
                       mle=c(cm$prim, cm$lsq),
                       cross=c(cm$prim, cm$lsq, cm$cross)
                       )

  ## index (start, end) to the various cathegories of
  ## information to be stored
  tblidx <- cumsum(c(0, length(PARAM), nCoVariates,
                     rep(1, length(allmethods) - length(allprimmeth)),
                     nindex, givenCoVariates * nCoVariates,
                     nindex, givenCoVariates * nCoVariates))
  tblidx <- rbind(tblidx[-length(tblidx)] + 1, tblidx[-1])
  dimnames(tblidx) <- list(c("start", "end"),
                           c("variab", "covariab",
                             allmethods[-1:-length(allprimmeth)],
                             "lower", "lowbeta", "upper", "upbeta"))
  maxtblidx <- tblidx[length(tblidx)]
  tblidx <- data.frame(tblidx)

  ## table of all information; col:various method; row:information to method
  varnames <- names(PrepareModel(ctr, named =TRUE)$param)
  var.idx <- which("var" == varnames)  # "var"=variance
  stopifnot(length(var.idx) > 0, var.idx[1]==1) ## assumed that "var" is always the first 
  var.idx <- c(var.idx, length(varnames) +1)
  for (i in 1:(length(var.idx)-1))
    varnames[var.idx[i] : (var.idx[i+1] - 1)] <-
      paste(varnames[var.idx[i] : (var.idx[i+1] - 1)],i, sep=".")
  betanames <- paste("beta", 1:nCoVariates, sep=".")

  tablenames <- c(varnames,
                  betanames,
                  allmethods[-1:-length(allprimmeth)],
                  ## do not try to join the next two lines, since both
                  ## varnames and betanames may contain nonsense if
                  ## nindex==0 and !givenCoVariates, respectively
                  if (nindex>0) paste("lower", varnames[index], sep=":"),
                  if (givenCoVariates) paste("lower", betanames, sep=":"),
                  if (nindex>0) paste("upper", varnames[index], sep=":"),
                  if (givenCoVariates) paste("upper", betanames, sep=":")
                  )
  
  param.table <- matrix(1 * NA, nrow=maxtblidx, ncol=length(allmethods),
                        dimnames=list(tablenames, allmethods))
  param.table <- data.frame(param.table)

  
##########  Trafo def + bounds for LSQ, also used for autostart  ###############

  LSQTRANSFORM <- users.transform  ## note: LSQTRANSFORM is used also in MLE
  LSQinvTrafo <- function(param) param
  LSQINDEX <- index
  nLSQINDEX <- sum(LSQINDEX)
  lsqtrafo <- NULL
  lsqlower <- lower
  lsqupper <- upper
  if (varnugNA) {
    lsqupper[NUGGET] <- 1
    LSQINDEX[VARIANCE] <- FALSE
    lsqtrafo <- function(param) {param[VARIANCE] <- 1.0 - param[NUGGET]; param}
    LSQinvTrafo <- function(param) {
      param[NUGGET] <- param[NUGGET] / (param[NUGGET] + param[VARIANCE])
      param[VARIANCE] <- NA
      param
    }
  } else if (sillbounded) {
    LSQINDEX[NUGGET] <- FALSE;
    lsqtrafo <- function(param) {param[NUGGET] <- sill - param[VARIANCE]; param}
  } else if (zeronugget) {
    if (sum(!LSQINDEX)>1) LSQINDEX[VARIANCE] <- FALSE ## otherwise we get
    ##           again a special case to treat
    lsqtrafo <- function(param) {param[VARIANCE] <- 1.0; param}
  }
  if (is.null(LSQTRANSFORM)) LSQTRANSFORM <- lsqtrafo
  else if (!is.null(lsqtrafo)) {
    warning("standard.style and transform!=NULL may cause strange effects -- internal transformation is performed before the user's one")
    lsqtrafoUser <- LSQTRANSFORM
    LSQTRANSFORM <- function(param) lsqtrafoUser(lsqtrafo(param))
  }
  if (is.null(LSQTRANSFORM)) LSQTRANSFORM <- function(x) x
  LSQLB  <- lsqlower[LSQINDEX]
  LSQUB  <- lsqupper[LSQINDEX]
  ixdLSQINDEX <- LSQINDEX[index]

  
##################################################
###############    PRIMITIVE METHODS   ###########
  ##****************    autostart    *****************
  ## for historical reasons, autostart uses the coding of LSQtarget;
  ## hence, LSQTRANSFORM must be called to get the true parameter values
  M <- "autostart"
  idx <- tblidx[["variab"]]
#print(autostart)
  if (!is.null(users.transform)) autostart <- users.transform(autostart)
#print(autostart)

  
#  print(LSQTRANSFORM)
#  print(lsqtrafoUser)
 # print(lsqtrafo)
#  xxx
  if (nLSQINDEX > 0) param.table[[M]][idx[1]:idx[2]] <- autostart 
  idx <- tblidx[["covariab"]]
  param.table[[M]][idx[1]:idx[2]] <- lsq.covariates(autostart)
  default.param <- param.table[["autostart"]]


  ##****************    user's guess    *****************
  ## user's guess is already clear text
  if (!is.null(users.guess)) {
    M <- "users.guess"
    ug <- if (is.list(users.guess)) {
      PrepareModel(users.guess, NULL, new$spacedim+new$Time, trend=users.beta)
    } else {
        if (is.null(param))
          stop("cannot interpret users.guess -- better use the list definition for model to define user's guess")
        PrepareModel(save.model, users.guess, new$spacedim+new$Time,
                     trend=users.beta)
      }
    ug <- convert.to.readable(ug)
    ug <- PrepareModel(ug$model, ug$param, new$spacedim+new$Time,
                       trend=ug$trend)$param
    if (!is.null(users.transform) &&  up != users.transform(ug))
      stop("user's guess does not satisfy the user's transformation function") 
    if (length(ug)!=length(PARAM))
      stop("model given by users.guess does not match 'model'")
    idx <- tblidx[["variab"]]
    param.table[[M]][idx[1]:idx[2]] <- ug
    idx <- tblidx[["covariab"]]
    if (is.null(trend)) {
      param.table[[M]][idx[1]:idx[2]] <-
        if (is.na(pm$mean)) {
          stopifnot(!is.null(users.beta))
          users.beta
        } else {
          stopifnot(is.null(users.beta))
          pm$mean
        }
    } else {
      stop("not programmed yet")
      stopifnot(!is.null(users.beta))
      users.beta <- as.double(users.beta)
      stopifnot(length(users.beta) == ncol(CoVariates))
      param.table[[M]][idx[1]:idx[2]] <- users.beta
    }
  }

  if (standard.style && use.naturalscaling && index[SCALE]) {
      idx <- tblidx[["variab"]]
      for (i in 1:length(allprimmeth)) if (!is.na(param.table[1, i])) {
      GNS <- .C("GetNaturalScaling",
                covnr,
                as.double(param.table[idx[1]:idx[2], i][-1:-(KAPPA-1)]),
                scalingmethod,
                natscale=double(1),
                error=integer(1),
                PACKAGE="RandomFields", DUP=FALSE)
      if (GNS$error)
        stop(paste("Error", error, "occured whilst rescaling"))
      param.table[idx[1]:idx[2], i][SCALE] <-
        param.table[idx[1]:idx[2], i][SCALE] / GNS$natscale
    }
  }

                                  
##################################################
###################  LSQ  ########################
  ## see above for the trafo definitions
  ##
  ## iterativer fit des trends: zuerst regressions fit,
  ## dann schaetzung der Parameter, dann lsq.covariates fit
  ##************   Empirical Variogram    ***********
  if (PrintLevel>4) cat("\nempirical variogram...")
  if (givenCoVariates) { 
    regr <- lsfit(CoVariates, data, intercept=FALSE)
    TREND <- regr$coeff
    EVtargetV <- regr$residuals
  } else EVtargetV <- data

  if (length(nphi)==1) nphi <- c(0, nphi) # starting angle; lines per half circle
  if (length(ntheta)==1) ntheta <- c(0, ntheta) # see above
  if (length(ntime)==1) ntime <- c(ntime, 1) 
  ntime <- ntime * T[3] ## time endpoint; step
  
  ev <- EmpiricalVariogram(new$x, T=new$T, data=EVtargetV, grid=FALSE, 
                           bin=if (length(bins)>1) bins else 
                           c(-1, seq(0, distance.factor * maxdistances,
                                    len=bins+1)),
                           phi=if ((new$spacedim>=2) && pm$anisotropy) nphi,
                           theta=if ((new$spacedim>=3) && pm$anisotropy) ntheta,
                           deltaT=if (!is.null(T)) ntime
                           )
  
  index.bv <- as.vector(!is.na(ev$e)) ## exclude bins without entry
  if (sum(index.bv)==0)
    stop("only NAs in empirical variogram; check values of bins and distance.factor")
  binned.variogram <- as.double(ev$e[index.bv])
 
  bin.centers <- as.matrix(ev$c)
  if (pm$anisotropy) {
    ## complete the coordinates of bin.centers to a vector of the dimension
    ## considered
    if (!is.null(ev$phi)) {
      if (spacedim<2) stop("x dimension is less than two, but phi is given") 
      bin.centers <- cbind(as.vector(outer(bin.centers, cos(ev$phi))),
                           as.vector(outer(bin.centers, sin(ev$phi))))
    }
    if (!is.null(ev$theta)) {
      if (spacedim<3)
        stop("x dimension is less than three, but theta is given") 
      if (ncol(bin.centers)==1) bin.centers <- cbind(bin.centers, 0)
      bin.centers <- cbind(as.vector(outer(bin.centers[, 1], cos(ev$theta))),
                           as.vector(outer(bin.centers[, 2], cos(ev$theta))),
                           rep(sin(ev$theta), each=nrow(bin.centers)))
    } else {
      if (nrow(bin.centers) < spacedim) # dimension of bincenter vector
        ##                       smaller than dimension of location space
        bin.centers <- 
          cbind(bin.centers, matrix(0, nrow=nrow(bin.centers),
                                    ncol=spacedim - ncol(bin.centers)
                                    ))
    }
    if (!is.null(ev$T)) {
      bin.centers <-
        cbind(matrix(rep(t(bin.centers), length(ev$T)), byrow=TRUE,
                     ncol = ncol(bin.centers)),
              rep(ev$T, each=nrow(bin.centers)))      
    }
  }
 
  bin.centers  <- as.double(t(bin.centers[index.bv, ])) #
  ##  es muessen beim direkten C-aufruf die componenten der Punkte
  ##  hintereinander kommen (siehe auch variable distance). Deshalb t()

  evsd <- as.double(ev$sd)
  evsd[is.na(evsd)] <- 0
  evsd[evsd==0] <- 10 * sum(evsd, na.rm=TRUE) ## == "infinity"
  
  bins             <- length(ev$n)
  binned.n         <- as.integer(ev$n)
  weights <- cbind(NA,                      # self
                   rep(1, bins),            # plain 
                   sqrt(binned.n),          # sqrt(#)
                   1 / evsd,                # sd^-1
                   sqrt(bins:1 * binned.n), # internal1
                   bins:1,                  # internal2
                   sqrt(bins:1)             # internal3
                   )[index.bv, ]
  stopifnot(ncol(weights)==length(alllsqmeth))
  dimnames(weights) <- list(NULL, alllsqmeth)
  weights <- data.frame(weights)
  bins <- as.integer(sum(index.bv))
  EVtargetV <- NULL

 
  ##***********   estimation part itself   **********     
  ## find a good initial value for MLE using weighted least squares
  ## and binned variogram
  ##
  ## background: if the number of observations (and the observation
  ## field) tends to infinity then any least square algorithm should
  ## yield the same result as MLE
  ## so the hope is that for a finite number of points the least squares
  ## find an acceptable initial values
 

  ## advantage of the following way is that the for-loop is run through
  ## in an ordered sense -- this might be useful in case partial results
  ## are reused
  methods <- (if (is.null(lsq.methods)) NULL else
              lsq.orig.methods[pmatch(lsq.methods, lsq.orig.methods)])
  if (any(is.na(methods))) stop("not all lsq.methods could be matched")
  if ("internal" %in% methods)
    methods <- c(methods, paste("internal", 1:nlsqinternal, sep=""))
  
  for (M in c(alllsqmeth)) {
    if (!(M %in% methods)) next;
    if (PrintLevel>2) cat("\n", M) else cat(pch)
    param.table[[M]] <- default.param
    LSQsettings(M)
    LSMIN <- Inf ## must be before next "if (nLSQINDEX==0)"
    LSPARAM <- NA 
    if (nLSQINDEX == 0) {
      LStarget(numeric(0))
    } else {
      idx <- tblidx[["lower"]]
      param.table[[M]][idx[1]:idx[2]][ixdLSQINDEX] <- LSQLB
      idx <- tblidx[["upper"]]
      param.table[[M]][idx[1]:idx[2]][ixdLSQINDEX] <- LSQUB
      options(show.error.messages = show.error.message) ##
      if (nLSQINDEX == 1) {
        variab <- try(optimize(LStarget, lower = LSQLB, upper = LSQUB)$minimum,
                      silent=!debug)
      } else {
        min <- Inf
        for (i in methodprevto$lsq) { ## ! -- the parts that change if
          ##                               this part is copied for other methods
          idx <- tblidx[["variab"]]
          if (!is.na(param.table[1, i])) {
            variab <- LSQinvTrafo(param.table[idx[1]:idx[2], i])[LSQINDEX]
            value <- LStarget(variab) ## !
            if (is.finite(value)) {
              param.table[tblidx[[M]][1], i] <- value
              if (value < min) {
                min.variab <- variab
                min <- value
              }
            } else param.table[tblidx[[M]][1], i] <- NaN
          }
        }
        stopifnot(length(min.variab) == length(LSQLB))
        
        lsq.optim.control <-
          c(optim.control, parscale=list(parscale[LSQINDEX]), fnscale=min)

#        print(parscale)
#        print(parscale[LSQINDEX])
        
        variab <- ## fnscale=1: minimisation
          try(optim(min.variab, LStarget, method ="L-BFGS-B", lower = LSQLB,
                    upper = LSQUB, control= lsq.optim.control)$par, silent=!debug)
      }
    }
    options(show.error.messages = show.error.message)  
    ## side effect: minimum so far is in LSMIN and LSPARAM
    ## even if the algorithm finally fails
    if (is.finite(LSMIN)) {
      idx <- tblidx[["variab"]]
      param.table[[M]][tblidx[[M]][1]] <- LSMIN
      param.table[[M]][idx[1]:idx[2]] <- LSPARAM
      idx <- tblidx[["covariab"]]
      param.table[[M]][idx[1]:idx[2]] <- lsq.covariates(LSPARAM)
    } else {
      param.table[[M]] <- NaN
    }
  } # for M


##################################################
### optional parameter grid for MLE and CROSS  ###

  idx <- tblidx$variab
  gridmax <- as.matrix(param.table[idx[1]:idx[2], cm$lsq])
  gridmin <- apply(gridmax, 1, min, na.rm=TRUE)
  gridmax <- apply(gridmax, 1, max, na.rm=TRUE)
  gridbound <- lower
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  abase <- 0.25
  a <- is.na(gridmin[idx]) * (1-abase) + abase
  ## maybe there have not been any lsq estimate; then a=1
  gridmin[idx] <- (1-a) * gridmin[idx] + a * gridbound[idx]
  gridbound <- upper
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  a <- is.na(gridmax[idx]) * (1-abase) + abase
  gridmax[idx] <- (1-a) * gridmax[idx] + a * gridbound[idx]


  
##################################################
###################   MLE    #####################
  
  methods <- (if (is.null(mle.methods)) NULL else
              allmlemeth[pmatch(mle.methods, allmlemeth)])
  if ("reml" %in% methods && !givenCoVariates) methods <- c(methods, "ml")
  mlelower <- lower
  mleupper <- lsqupper
  MLEINDEX <- LSQINDEX
  nMLEINDEX <- sum(MLEINDEX)
  MLETRANSFORM <- LSQTRANSFORM 
  MLEinvTrafo <- LSQinvTrafo
  ## lowerbound.scale.LS.factor <  lowerbound.scale.factor, usually
  ## LS optimisation should not run to a boundary (what often happens
  ## for the scale) since a boundary value is usually a bad initial
  ## value for MLE (heuristic statement). Therefore a small
  ## lowerbound.scale.LS.factor is used for LS optimisation.
  ## For MLE estimation we should include the true value of the scale;
  ## so the bounds must be larger. Here lower[SCALE] is corrected
    ## to be suitable for MLE estimation
  if (pm$anisotropy)
    mleupper[scale.pos>0] <- mleupper[scale.pos>0] *
      lowerbound.scale.factor / lowerbound.scale.LS.factor
  else 
    mlelower[scale.pos>0] <- mlelower[scale.pos>0] *
      lowerbound.scale.LS.factor / lowerbound.scale.factor
  MLELB  <- mlelower[MLEINDEX]
  MLEUB  <- mleupper[MLEINDEX]
  ixdMLEINDEX <- MLEINDEX[index]

  ## fnscale <- -1 : maximisation
  for (M in c(allmlemeth)) {
    if (!(M %in% methods)) next;
    if (PrintLevel>2) cat("\n", M) else cat(pch)
    param.table[[M]] <- default.param
    if (M=="reml" && !givenCoVariates) { ## same as MLE
      param.table[[M]] <- param.table[["ml"]]
      param.table[[M]][tblidx[[M]][1]] <- param.table[[M]][tblidx[["ml"]][1]]
      next
    }
    MLEsettings(M)
    MLEMAX <- -Inf ## must be before next "if (nMLEINDEX==0)"
    if (nMLEINDEX == 0) {
      MLEtarget(numeric(0))
    } else {
      idx <- tblidx[["lower"]]
      param.table[[M]][idx[1]:idx[2]][ixdMLEINDEX] <- MLELB
      idx <- tblidx[["upper"]]
      param.table[[M]][idx[1]:idx[2]][ixdMLEINDEX] <- MLEUB
      options(show.error.messages = show.error.message) ##
      if (nMLEINDEX == 1) {
        variab <- try(optimize(MLEtarget, lower = MLELB,
                               upper = MLEUB, maximum=TRUE)$maximum,
                      silent=!debug)
      } else {
        min <- Inf
        for (i in methodprevto$mle) { ## ! -- the parts that change if
          ##                               this part is copied for other methods
          idx <- tblidx[["variab"]]
          if (!is.na(param.table[1, i])) {
            variab <- MLEinvTrafo(param.table[idx[1]:idx[2], i])[MLEINDEX]
            value <- MLEtarget(variab) ## !
            if (is.finite(value)) {
              param.table[tblidx[[M]][1], i] <- value
              if (value < min) {
                min.variab <- variab
                min <- value
              }
            } else param.table[tblidx[[M]][1], i] <- NaN
          }
        }
        stopifnot(length(min.variab) == length(MLELB))        
        mle.optim.control <-
          c(optim.control, parscale=list(parscale[MLEINDEX]), fnscale=-min)
        variab <-
          try(optim(min.variab, MLEtarget, method="L-BFGS-B", lower = MLELB,
                    upper = MLEUB, control=mle.optim.control)$par,
              silent=!debug)
      }
      options(show.error.messages = TRUE) ##
      variab <- MLEPARAM[MLEINDEX] ## to check onborderline
      mindistance <- pmax(minbounddistance, minboundreldist * abs(variab))
      onborderline <- (any(abs(variab - MLELB) <
                           pmax(mindistance,              ## absolute difference
                                minboundreldist * abs(MLELB)##relative difference
                                )) ||
                       any(abs(variab - MLEUB) <
                           pmax(mindistance, minboundreldist * abs(MLEUB))))
    }
    idx <- tblidx[["variab"]]
    if (is.finite(MLEMAX)) {
      param.table[[M]][tblidx[[M]][1]] <- MLEMAX
      param.table[[M]][idx[1]:idx[2]] <- MLEPARAM
    } else {
      if (PrintLevel>0) cat(M, "MLEtarget I failed.\n")
      param.table[[M]] <- MLEPARAM <- NaN
      variab <- MLELB ## to call for onborderline
    }

    if (nMLEINDEX > 0 && onborderline && refine.onborder) {
      ## if the MLE result is close to the border, it usually means that
      ## the algorithm has failed, especially because of a bad starting
      ## value (least squares do not always give a good starting point, helas)
      ## so the brutal method:
      ## calculate the MLE values on a grid and start the optimization with
      ## the best grid point. Again, there is the believe that the
      ## least square give at least a hint what a good grid is
      MLEgridmin <- MLEinvTrafo(gridmin)[MLEINDEX]
      MLEgridmax <- MLEinvTrafo(gridmax)[MLEINDEX]
      if (any(is.na(MLEgridmin)) || any(is.na(MLEgridmax))) {
        warning(paste(M, "converged to a boundary value -- better performance might be obtained when allowing for more lsq.methods"))
      } else {
        if (PrintLevel>5) show(1, M, MLEMAX, MLEPARAM) else cat(detailpch)
        MLEgridlength <- max(3, round(approximate.functioncalls ^ (1/nMLEINDEX)))
        ## grid is given by the extremes of the LS results
        ## so, therefore we should examine above at least 4 different sets
        ## of weights
        ## wichtig: gridmin/max basiert auf den reduzierten variablen
        step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-2) # grid starts
        ##                                                      bit outside
        MLEgridmin <- pmax(MLEgridmin - step/2, MLELB)     # the extremes of LS
        MLEgridmax <- pmin(MLEgridmax + step/2, MLEUB)       
        step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-1)
        i <- 1
        zk <-  paste("MLEgridmin[",i,"] + step[",i,"] * (0:",MLEgridlength-1,")")
        if (length(step)>1)
          for (i in 2:length(step))
            zk <- paste(zk,",MLEgridmin[",i,"] + step[",i,"] * (0:",
                        MLEgridlength-1,")")
        zk <- paste("expand.grid(",zk,")")
        startingvalues <- eval(parse(text=zk))
        limit <- 10 * approximate.functioncalls
        if ((rn <- nrow(startingvalues)) > limit) {
          if (PrintLevel>4)
            cat("using only a random subset of the", rn, "grid points")
          rand <- runif(rn)
          startingvalues <- startingvalues[rand < quantile(rand, limit / rn), ]
          gc()
        }
        
        MLEMAX <- -Inf
        apply(startingvalues, 1, MLEtarget) ## side effect: Maximum is in MLEMAX!
        ##  optimal parameter is in MLEPARAM
        if (PrintLevel>5) show(2, M, MLEMAX, MLEPARAM)
        if (nMLEINDEX > 1) {
          cat(detailpch)
          variab <- MLEinvTrafo(MLEPARAM)[MLEINDEX]
          options(show.error.messages = show.error.message) ##
          variab <-
            try(optim(variab, MLEtarget, method ="L-BFGS-B",lower = MLELB,
                      upper = MLEUB, control=mle.optim.control)$par,
                silent=!debug)
          options(show.error.messages = TRUE) ##
          if (!is.numeric(variab) && (PrintLevel>0)) cat("MLtarget II failed.\n")
          ## do not check anymore whether there had been convergence or not.
          ## just take the best of the two strategies (initial value given by
          ## LS, initial value given by a grid), and be happy.
          if (PrintLevel>5) show(3, M, MLEMAX, MLEPARAM)
        }

        idx <- tblidx[["variab"]]
        if (is.finite(MLEMAX) && MLEMAX > param.table[[M]][tblidx[[M]][1]]) {
          param.table[[M]][tblidx[[M]][1]] <- MLEMAX
          param.table[[M]][idx[1]:idx[2]] <- MLEPARAM
        }
      } # (is.na(MLEgridmin[1]))
    } # onborderline
    
    if (is.finite(param.table[[M]][tblidx[[M]][1]])) {
      varidx <- tblidx[["variab"]]
      idx <- tblidx[["covariab"]]
      param.table[[M]][idx[1]:idx[2]] <-
        lsq.covariates(param.table[[M]][varidx[1]:varidx[2]])
    }
  } ## M

    
########  estimation by cross validation  ########                      
  data <- as.matrix(data)
  CROSSINDEX <- index
  nCROSStotINDEX <- nCROSSINDEX <- sum(CROSSINDEX)
###############################################################
##  besser: lsqlower oder mlelower, dann andere Trafo, etc!  ## 
###############################################################
  crosslower <- lower
  crossupper <- upper
  CROSSTRANSFORM <- users.transform
  CROSSinvTrafo <- function(param) param
  crosstrafo <- NULL;
  if (sillbounded) {
    CROSSINDEX[NUGGET] <- FALSE;
    crosstrafo <- function(param) {
      param[NUGGET] <- sill - param[VARIANCE]
      param
    }
  }
  if (is.null(CROSSTRANSFORM)) CROSSTRANSFORM <- crosstrafo
  else if (!is.null(crosstrafo)) {
    warning("standard.style and transform!=NULL may cause strange effects -- internal transformation is performed before the user's one")
    crosstrafoUser <- CROSSTRANSFORM
    CROSSTRANSFORM <- function(param) crosstrafoUser(crosstrafo(param))
  }
  if (is.null(CROSSTRANSFORM)) CROSSTRANSFORM <- function(x) x
  CROSSLB <- crosslower[CROSSINDEX] 
  CROSSUB <- crossupper[CROSSINDEX]
  ixdCROSSINDEX <- CROSSINDEX[index]

  methods <- (if (is.null(cross.methods)) NULL else
              allcrossmeth[pmatch(cross.methods, allcrossmeth)])
  CROSS.lcrepet <- lc * repet
  cross.optim.control <-
    c(optim.control, parscale=list(parscale[CROSSINDEX]), fnscale=1)
  crossLBcovariates <-  crossUBcovariates <- NULL
  if (givenCoVariates) {
    stopifnot(is.null(trend))
    crossLBcovariates <- min(data)
    crossUBcovariates <- max(data)
    ## die neuen Grenzen sind 1. katastrophal schlecht; 2. werden
    ## sie nicht in die Tabelle uebernommen!!; 3. gibt es chaos mit
    ## den schnellen loesungen
    CROSSLB <- c(CROSSLB, crossLBcovariates)
    CROSSUB <- c(CROSSUB, crossUBcovariates)
    nCROSStotINDEX <- nCROSSINDEX + nCoVariates
    cross.optim.control$parscale <- ### auch nicht gut !!!
      c(cross.optim.control$parscale, rep(1, nCoVariates))
  }
  for (M in c(allcrossmeth)) {
    if (!(M %in% methods)) next;
    if (PrintLevel>2) cat("\n", M) else cat(pch)
    stopifnot(is.null(trend)) ## vuniversal kriging not programmed yet
    crosssettings(M)    
    CROSSMIN <- Inf
    param.table[[M]] <- default.param ## in case the covariates are not estimated

    ## optimisation part
    RFparameters(Print=0)
    if (nCROSStotINDEX == 0) {
      crosstarget(numeric(0))
    } else {
      if (length(ixdCROSSINDEX) > 0) {
        idx <- tblidx[["lower"]]
        param.table[[M]][idx[1]:idx[2]][ixdCROSSINDEX] <- CROSSLB[1:nCROSSINDEX]
        idx <- tblidx[["upper"]]
        param.table[[M]][idx[1]:idx[2]][ixdCROSSINDEX] <- CROSSUB[1:nCROSSINDEX]
      }
      if (givenCoVariates) {
        idx <- tblidx[["lowbeta"]]
        param.table[[M]][idx[1]:idx[2]] <- CROSSLB[nCROSSINDEX + 1:nCoVariates]
        idx <- tblidx[["upbeta"]]
        param.table[[M]][idx[1]:idx[2]] <- CROSSUB[nCROSSINDEX + 1:nCoVariates]
      }
      options(show.error.messages = show.error.message) ##  
      if (nCROSStotINDEX==1) {
        variab <-
          try(optimize(crosstarget, lower = CROSSLB, upper = CROSSUB)$minimum,
              silent=!debug)
      } else {
        min <- Inf
        for (i in methodprevto$cross) { ## ! -- the parts that change if
          ##                               this part is copied for other methods
          idx <- tblidx[["variab"]]
          idxCovar <-  tblidx[["covariab"]]
          if (!is.na(param.table[1, i])) {
            variab <- CROSSinvTrafo(param.table[idx[1]:idx[2], i])[CROSSINDEX]
            if (givenCoVariates)
              variab <- c(variab, param.table[idxCovar[1]:idxCovar[2], i])
            value <- crosstarget(variab) ## !
            if (is.finite(value)) {
              param.table[tblidx[[M]][1], i] <- value
              if (value < min) {
                min.variab <- variab
                min <- value
              }
            } else param.table[tblidx[[M]][1], i] <- NaN
          }
        }
        stopifnot(length(min.variab) == length(CROSSLB))
        variab <-
          try(optim(min.variab, crosstarget, method ="L-BFGS-B", lower = CROSSLB,
                    upper = CROSSUB, control = cross.optim.control)$par,
              silent=!debug)
      } # nCROSStotINDEX > 1

      ## check onborderline
      variab <- CROSSMODEL$param[CROSSINDEX]
      if (givenCoVariates) variab <- c(variab, CROSSMODEL$mean)

      mindistance <- pmax(minbounddistance, minboundreldist * abs(variab))
      onborderline <-
        (any(abs(variab - CROSSLB) <
             pmax(mindistance,              ## absolute difference
                  minboundreldist * abs(CROSSLB)##relative difference
                  )) ||
         any(abs(variab - CROSSUB) <
             pmax(mindistance, minboundreldist * abs(CROSSUB))))
    } # nCROSStotINDEX > 0
    options(show.error.messages = TRUE) ##
    if (is.finite(CROSSMIN)) {
      idx <- tblidx[["variab"]]
      param.table[[M]][tblidx[[M]][1]] <- CROSSMIN
      param.table[[M]][idx[1]:idx[2]] <- CROSSMODEL$param
      if (givenCoVariates) {
        stopifnot(is.null(trend))
        idx <- tblidx[["covariab"]]
        param.table[[M]][idx[1]:idx[2]] <- CROSSMODEL$mean
      }
    } else {
      if (PrintLevel>0) cat(M, "target I failed.\n")
      param.table[[M]] <- NaN
    }
    if (nCROSStotINDEX > 0 && onborderline && refine.onborder) {
      CROSSgridmin <- c(CROSSinvTrafo(gridmin)[CROSSINDEX], crossLBcovariates)
      CROSSgridmax <- c(CROSSinvTrafo(gridmax)[CROSSINDEX], crossUBcovariates)
      if (any(is.na(CROSSgridmin)) || any(is.na(CROSSgridmax))) {
        warning(paste(M, "converged to a boundary value -- better performance might be obtained when allowing for more lsq.methods"))
      } else {
        if (PrintLevel>5) show(1, M, CROSSMIN, CROSSMODEL$param)
        else cat(detailpch)
        CROSSgridlength <-
          max(3, round(approximate.functioncalls ^
                       (1/(nCROSStotINDEX + givenCoVariates * nCoVariates))))
        step <- (CROSSgridmax - CROSSgridmin) / (CROSSgridlength - 2) 
        CROSSgridmin <- pmax(CROSSgridmin - step/2, CROSSLB)    
        CROSSgridmax <- pmin(CROSSgridmax + step/2, CROSSUB)       
        step <- (CROSSgridmax - CROSSgridmin) / (CROSSgridlength - 1)
        i <- 1
        zk <-  paste("CROSSgridmin[", i, "] + step[", i, "] * (0:",
                     CROSSgridlength - 1,")")
        if (length(step)>1)
          for (i in 2:length(step))
            zk <- paste(zk,",CROSSgridmin[", i, "] + step[", i, "] * (0:",
                        CROSSgridlength - 1,")")
        zk <- paste("expand.grid(", zk, ")")
        startingvalues <- eval(parse(text=zk))
        limit <- 10 * approximate.functioncalls
        if ((rn <- nrow(startingvalues)) > limit) {
          if (PrintLevel>4)
            cat("using only a random subset of the", rn, "starting values")
          rand <- runif(rn)
          startingvalues <- startingvalues[rand < quantile(rand, limit / rn), ]
          gc()
        }
        CROSSMIN <- Inf
        apply(startingvalues, 1, crosstarget) ## side effect: Min. in C.-MODEL!
        if (PrintLevel>2) show(2, M, CROSSMIN, CROSSMODEL$param)
        if (nCROSStotINDEX>1) {
          cat(detailpch)
          variab <- CROSSinvTrafo(CROSSMODEL$param)[CROSSINDEX]
          if (givenCoVariates) {
            variab <- c(variab, param.table[idxCovar[1]:idxCovar[2], i])
          }
          options(show.error.messages = show.error.message) ##
          variab <-
            try(optim(variab, crosstarget, method ="L-BFGS-B",lower = CROSSLB,
                      upper = CROSSUB, control=cross.optim.control)$par,
                silent=!debug)
          options(show.error.messages = TRUE) ##
          if (!is.numeric(variab) && (PrintLevel>0))
            cat("cross target II failed.\n")
        }
        if (PrintLevel>2) show(3, M, CROSSMIN, CROSSMODEL$par) 

        if (is.finite(CROSSMIN) && CROSSMIN < param.table[[M]][tblidx[[M]][1]]) {
          idx <- tblidx[["variab"]]
          param.table[[M]][tblidx[[M]][1]] <- CROSSMIN
          param.table[[M]][idx[1]:idx[2]] <- CROSSMODEL$param
          if (givenCoVariates) {
            stopifnot(is.null(trend))
            idx <- tblidx[["covariab"]]
              param.table[[M]][idx[1]:idx[2]] <- CROSSMODEL$mean
          }
        }
      } # (is.na(CROSSgridmin[1]))
    } # onborderline
    RFparameters(Print=save.RFparameters$PrintLevel)
    
    if (is.finite(param.table[[M]][tblidx[[M]][1]])) {
      varidx <- tblidx[["variab"]]
      idx <- tblidx[["covariab"]]
      param.table[[M]][idx[1]:idx[2]] <-
        lsq.covariates(param.table[[M]][varidx[1]:varidx[2]])
    }
  } ## for M

  
######################################################################
###     calculate all target values for all optimal parameters     ###
######################################################################
  if (table.format) {
    for (i in 1:length(allmethods)) if (!is.na(param.table[1, i])) {
      for (M in alllsqmeth) {
        cur <- param.table[tblidx[[M]][1], i]
        if (is.na(cur) && !is.nan(cur)) {
          LSQsettings(M)
          param.table[tblidx[[M]][1], i] <-
            LStarget(LSQinvTrafo(param.table[, i])[LSQINDEX])
        }
      }

      for (M in allmlemeth) {
        cur <- param.table[tblidx[[M]][1], i]
        if (is.na(cur) && !is.nan(cur)) {
          MLEsettings(M)
          param.table[tblidx[[M]][1], i] <-
            MLEtarget(MLEinvTrafo(param.table[, i])[MLEINDEX])
        }
      }
      
      idx <- tblidx[["variab"]]
      idxCovar <-  tblidx[["covariab"]]
      for (M in allcrossmeth) {
        cur <- param.table[tblidx[[M]][1], i]
        if (is.na(cur) && !is.nan(cur)) {
          crosssettings(M)
          variab <- CROSSinvTrafo(param.table[idx[1]:idx[2], i])[CROSSINDEX]
          if (givenCoVariates) {
            variab <- c(variab, param.table[idxCovar[1]:idxCovar[2], i])
          }
          param.table[tblidx[[M]][1], i] <- crosstarget(variab)
        }
      }
    }
  }
  if (pch!="") cat("\n")

  
######################################################################
###             rescaling in case of old.style                     ###
######################################################################
  ## if the covariance functions use natural scaling, just
  ## correct the final output by GNS$natscale
  ## (GNS$natscale==1 if no rescaling was used)
  ##
  ## currently natural scaling only for standard.style...
  if (standard.style && use.naturalscaling && index[SCALE]) {
    idx <- tblidx[["variab"]]
    for (i in 1:length(allmethods)) if (!is.na(param.table[1, i])) {
      GNS <- .C("GetNaturalScaling",
                covnr,
                as.double(param.table[idx[1]:idx[2], i][-1:-(KAPPA-1)]),
                scalingmethod,
                natscale=double(1),
                error=integer(1),
                PACKAGE="RandomFields", DUP=FALSE)
      if (GNS$error)
        stop(paste("Error", error, "occured whilst rescaling"))
      param.table[idx[1]:idx[2], i][SCALE] <-
        param.table[idx[1]:idx[2], i][SCALE] * GNS$natscale
    }
  }

  
######################################################################
###                     error checks                               ###
######################################################################
  ## if rather close to nugget and nugget==fixed, do exception handling.
  ## This part is active only if
  ## scale.max.relative.factor < lowerbound.scale.LS.factor
  ## By default it is not, i.e., the following if-condition
  ## will "always" be FALSE.
  if (standard.style && !is.na(nugget)) {
    idx <- tblidx[["variab"]]
    alllsqscales <- param.table[idx[1]:idx[2], cm$lsq][SCALE, ]
    if (any(alllsqscales < mindistances/scale.max.relative.factor, na.rm=TRUE))
      warning(paste(sep="",
                    "Chosen model does not seem to be appropriate!\n Probably a",
                  if (nugget!=0.0) "larger",
                    "nugget effect should be considered.\n")
              )
  }
  
######################################################################
###                   format and return values                     ###
######################################################################
  # print(param.table, digits=3)
  if (table.format) return(param.table)
  
  ## else table.format=FALSE :

  r <- list(covnr=pm$covnr, anisotropy=pm$anisotropy, op=pm$op, mean=NA,
            trend=pm$trend, method=pm$method, timespacedim=pm$timespacedim)
  idx <- tblidx[["variab"]]
  idxCovar <- tblidx[["covariab"]]
  idx.meth <- rep(FALSE, length(allmethods))
  res <- values.res <- list()
  for (i in 1:length(allmethods)) {
    M <- allmethods[i]
    if (idx.meth[i] <- !is.na(param.table[1, i]) || is.nan(param.table[1,i])) {
      r$param <- param.table[idx[1]:idx[2], i]
      if (is.null(trend)) {
        r$mean <- param.table[idxCovar[1]:idxCovar[2], i]
      } else {
        stopifnot(is.null(trend))
      }
      res <- c(res, list(convert.to.readable(r)))
      values.res[[M]] <- param.table[[M]][tblidx[[M]][1]]
    }
  }
  names(res) <- names(param.table)[idx.meth]
  
  return(c(list(ev = ev),
           variogram=list(res),
           values=list(values.res)))
}

