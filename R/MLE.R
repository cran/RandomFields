
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



#  stop("")
  # problem: natscale; im moment 2x implementiert, 1x mal ueber
  # scale/aniso (user) und einmal gedoppelt -- irgendwas muss raus

## LSQ variogram fuer trend = const.
## kann verbessert werden, insb. fuer fixed effects, aber auch eingeschraenkt
## fuer random effects -> BA/MA


## REML fehlt


## NAs in data mit mixed model grundsaetzlich nicht kombinierbar !
## NAs in data mit trend (derzeit) nicht kombinierbar

## zentrale C -Schnittstellen
##    .C("PutValuesAtNA", param, PACKAGE="RandomFields", DUP=FALSE)

## bins bei Distances automatisch




fitvario <-
  function(x, y=NULL, z=NULL, T=NULL, data, model, param,
           lower=NULL, upper=NULL, sill=NA, grid=!missing(gridtriple),
           gridtriple,
           ...) {
    fitvario.default(x=x, y=y, z=z, T=T, data=data, model=model, param=param,
                     lower=lower, upper=upper, sill=sill,
                     grid=grid, gridtriple,     
                     ...)
  }

fitvario.default <-
  function(x, y=NULL, z=NULL, T=NULL, data, model, param,
           grid=!missing(gridtriple), gridtriple=FALSE,     
           trend = NULL,
##         BoxCox ((x+c)^l - 1) / l; log(x+c); with c==1
##         BC.lambda : NA / value
##         BC.c: nur value
           BC.lambda, ## if missing then no BoxCox-Trafo
           BC.lambdaLB=-10, BC.lambdaUB=10,
           lower=NULL, upper=NULL, sill=NA,
           use.naturalscaling=FALSE,
           ## speed=FALSE, i.e. coordniates will be saved in GATTER
           PrintLevel, optim.control=NULL,
           bins=20, nphi=1, ntheta=1, ntime=20,
           distance.factor=0.5,
           upperbound.scale.factor=3, lowerbound.scale.factor=3, 
           lowerbound.scale.LS.factor=5,
           upperbound.var.factor=10, lowerbound.var.factor=100,     
           lowerbound.sill=1E-10, scale.max.relative.factor=1000,
           minbounddistance=0.001, minboundreldist=0.02,
           approximate.functioncalls=50, refine.onborder=TRUE,
           minmixedvar=1/1000, maxmixedvar=1000,
           pch=RFparameters()$pch, 
           transform=NULL,  standard.style=NULL,
           var.name="X", time.name="T",
           lsq.methods=c("self", "plain", "sqrt.nr", "sd.inv", "internal"),
           ## "internal" : name should not be changed; should always be last
           ##              method!
           mle.methods=c("ml", "reml"), #, "reml"),
           cross.methods=NULL,
  #       cross.methods=c("cross.sq", "cross.abs", "cross.ign", "cross.crps"),
           users.guess=NULL,  only.users = FALSE,
           Distances=NULL, truedim,
           solvesigma = NA, # if NA then use algorithm -- ToDo
           allowdistanceZero = FALSE,
           na.rm = TRUE) { ## do not use FALSE for mixed models !
    
####################################################################
###                      Preludium                               ###
####################################################################
  mle.methods <- "ml"
    
 
  if (length(cross.methods) > 0) {
    stop("not programmed yet")
    if (is.null(standard.style)) standard.style <- FALSE
    else if (!standard.style)
      stop("if cross.methods then standard.style may not be used")
  }


  save.options <- options()
  
##  library(RandomFields, lib="~/TMP")
  old.param <- RFparameters(no.readonly=TRUE) ## valgrind bringt system fehler
  if (old.param$Practical && use.naturalscaling)
    stop("PracticalRange must be FALSE if use.naturalscaling=TRUE")
  if (use.naturalscaling) RFparameters(PracticalRange = 3)
  if (missing(PrintLevel)) PrintLevel <- old.param$PrintLevel
  else RFparameters(PrintLevel = PrintLevel)

  if (PrintLevel>1) cat("\nfunction defintions...\n")

 
  if (!FALSE) { ### on.exit currently causes errors !!
    on.exit(RFparameters(old.param))
    on.exit(options(save.options), add=TRUE)
  }
  
  ##all the following save.* are used for debugging only
  silent <- PrintLevel < 3   # optim/ize
  show.error.message <- TRUE # options

  detailpch <- if (pch=="") "" else "+"
  nuggetrange <- 1e-15 + 2 * RFparameters()$nugget.tol


### definitions from 
  VARPARAM <- 0  ## to be consistent with the C definitions in RF.h
  SCALEPARAM <- 1
  DIAGPARAM <- 2
  ANISOPARAM <- 3
  INTERNALPARAM <- 4
  ANYPARAM <- 5
  TRENDPARAM <- 6
  NUGGETVAR <- 7
  MIXEDVAR <- 8
  REGRESSION <- 9

  VARIANCE <- 1
  NUGGET <- 2
  SCALE <- 3
  KAPPA <- 4

  DeterministicEffect <- 0
  FixedEffect <- 1 ## trend is also converted to FixedEffect
  RandomEffect <- 2## b is random, no variance is estimated; covariance matrix
  RVarEffect <- 3 ## b is random, variance is estimated; covariance matrix
  LargeEffect <- 4## as RandomEffect, but huge covariance matrix
  LVarEffect <- 5## as RandomEffect, but huge covariance matrix
  SpaceEffect <- 6 ## complex covariance model
  RemainingError <- 7
  EffectName <- c("Determ", "FixEff", "RanEff", "RanEff",
                  "XRnEff", "XRnEff", "SpcEff")

  ModelNr <- 3  ## == MODEL_MLE in RF.h

  
  ## Note: upper case variables are global variables !!

  
######################################################################
###                function definitions                            ###
######################################################################
 
  BoxCox <- function(x, lambda)
     if (abs(lambda) < 10^-20) log(x) else (x^lambda - 1) / lambda


  show <- function(nr, M, OPT, PARAM)
    cat("\n ", M, ", ", switch(nr, "start", "grid ", "re-do"), ": value=",
        format(OPT, dig=6), ", param=", format(PARAM, dig=2), sep="")

##used for trend functions
  formulaToMatrix <- function(formula, coord, var.name="X", time.name="T") {
    l <- NULL ## otherwise check will give a warning
    coord <- as.matrix(coord)
    if (is.null(time.name)) co <- "time.name<-NULL"
    else co <- paste(time.name,"=coord[nrow(coord),]")
    for (i in 1:(nrow(coord)-!is.null(time.name)))
      co <- paste(co, ",", var.name, i, "=coord[", i, ",]", sep="")
    eval(parse(text=paste("l <- list(",co,")")))
    fo <- as.list(formula)[[length(as.list(formula))]]
    variables <- NULL
    while (length(fo)==3) {
      dummy <- as.list(fo)
      if(dummy[[1]]!="+") break
      fo <- dummy[[2]]
      if (class(dummy[[3]])=="numeric") {
	 text <- paste("~", dummy[[3]], "*", var.name, "1^0", sep="")
	 dummy[[3]] <- as.list(eval(parse(text=text)))[[2]]
      }
      variables <- cbind(eval(dummy[[3]],l), variables)
    }
    if (class(fo)=="numeric") {
	text <- paste("~", fo, "*", var.name, "1^0", sep="")
	fo <- as.list(eval(parse(text=text)))[[2]]
    }
    variables <- cbind(eval(fo,l), variables)
    return(variables)
  }

  EmptyEnv <- function() baseenv()

  LSQsettings <- function(M) {
    assign("LSQ.SELF.WEIGHING", M=="self", envir=ENVIR)
    if (!LSQ.SELF.WEIGHING) {
      assign("LSQ.WEIGHTS", weights[[M]], envir=ENVIR)
      if (globalsigma)
        assign("LSQ.BINNEDSQUARE", sum(binned.variogram^2 * LSQ.WEIGHTS),
               envir=ENVIR)
    }
  }

  WarningMessage <- function (variab, LB, UB, txt) {
    cat("Note:", txt, ": forbidden values -- if there are too many warnings",
        "try narrower lower and upper bounds for the variables: (",
        paste(variab, collapse=","), ") not in [(",
        paste(LB, collapse=", "),  ") ; (",
        paste(UB, collapse=", "), ")]\n")
  }

  
  LStarget <- function(variab) {
    variab <- variab + 0## unbedingt einfuegen, da bei R Fehler
    ##                     der Referenzierung !! 16.2.10
    if (PrintLevel>4) Print(LSMIN, format(variab, dig=20))

    if (n.variab==0) return(NA)		##trivial case

    if (any((variab<LSQLB) | (variab>LSQUB))) {            
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 13.12.03 still happens ...
      if (PrintLevel>2) WarningMessage(variab, LSQLB, LSQUB, "LSQ")
      penalty <- variab
      variab <- pmax(LSQLB, pmin(LSQUB, variab))
      penalty <- + sum(variab-penalty)^2
      res <- LStarget(variab)
      return(res + penalty * (1 + abs(res)))
    }

    param <- as.double(transform(variab[1:nvarWithoutBC]))
    .C("PutValuesAtNA", param, PACKAGE="RandomFields", DUP=FALSE)

    model.values <- double(bins * vdim^2)
    .C("VariogramMLE", bin.centers, NULL, bins, model.values,
       PACKAGE="RandomFields", DUP=FALSE)

    if (any(!is.finite(model.values))) {
      if (PrintLevel>1) {
        cat("\nLSQ missing values!")
        Print(format(variab, dig=20), format(param, dig=20))
        Print(cbind(bin.centers, as.matrix(model.values, ncol=vdim))) #
      }
      return(1E300)
    } 

    if (LSQ.SELF.WEIGHING) {
      ## weights = 1/ model.values^2
        gx <- binned.variogram / model.values
        gx <- gx[is.finite(gx)]
        if (globalsigma) {
          bgw <- sum(gx^2)
          g2w <- sum(gx)
          ergb <- length(gx) - g2w^2 / bgw

         } else {
          ergb <- sum((gx - 1)^2)
       }
    } else {
      if (globalsigma) {
        ## in this case the calculated variogram model is not the one
        ## to which we should compare, but:
        bgw <- sum(binned.variogram * model.values * LSQ.WEIGHTS)
        g2w <- sum(model.values^2 * LSQ.WEIGHTS)
        ergb <- LSQ.BINNEDSQUARE - bgw^2/g2w
      } else {
        ergb <- sum((binned.variogram - model.values)^2 * LSQ.WEIGHTS)
       }
    }
 
    if (ergb<LSMIN) { 
      if (varnugNA) {
        sill <- bgw/g2w
        param[var.global] <- sill * (1.0 - param[nugget.idx])
        param[nugget.idx] <- sill * param[nugget.idx]
      } else if (globalsigma) {
         param[var.global] <- bgw/g2w ## variance
      }
      assign("LSMIN", ergb, envir=ENVIR)
      assign("LSPARAM", param, envir=ENVIR)
      assign("LSVARIAB", variab, envir=ENVIR)
    }
    return(ergb)
  }

  

  MLEsettings <- function(M) {
    switch(M,
      "ml" = {
        base <- 0
        ML.totalRV <- ntotdata * sum(effect == RemainingError)

        if (globalsigma)
          base <- base + ML.totalRV * (1 - log(ML.totalRV))

        dummy <-
          rowSums(nCoVar[, effect >= RandomEffect && effect <= SpaceEffect,
                         drop=FALSE])
        if (solvesigma) {
          base <- base + sum(dummy * (1 - log(dummy)))
        }
 
        base <-  base + (sum(dummy) +  ML.totalRV) * log(2 * pi)
    
        assign("ML.totalRV", ML.totalRV, envir=ENVIR)
        assign("ML.base", base, envir=ENVIR)
        assign("MLEtarget", MLtarget, envir=ENVIR)
       },
      "reml" = {
        stop("not programmed yet")
        if (any(effect>2)) stop("mixed model not yet programmed in REML.")
        reml.lc <- sapply(CoVariates, function(x) nrow(x) - ncol(x))
        assign("REML.totalRV", reml.lc * repet, envir=ENVIR)
        sumREML.totalRV <- sum(REML.totalRV)
        assign("REML.loglcrepet", sumREML.totalRV * (1 -log(sumREML.totalRV)),
               envir=ENVIR)
        assign("REML.twopilcrepet", sumREML.totalRV * log(2 * pi),
               envir=ENVIR)
        reml.a0 <- lapply(CoVariates, function(x) eigen(tcrossprod(x)))
        reml.a <- reml.data <- vector("list", length(reml.a0))
        for (i in 1:length(reml.a0)) {
          reml.a[[i]] <-
            reml.a0[[i]]$vectors[ ,
                                 order(abs(reml.a0[[i]]$value))[1:reml.lc[i]],
                                 drop=FALSE]
          reml.data[[i]] <- crossprod(reml.a[[i]], werte[[i]]) 
        }
        assign("REML.A", reml.a, envir=ENVIR)
        assign("REML.data", reml.data, envir=ENVIR)
        assign("MLEtarget", REMLtarget, envir=ENVIR)        
       },
       stop(M, " unknown")
           ) # switch
  }


  linearpart <- function(sigma2, return.z=FALSE) {
    z <- list()
       
    for (i in 1:sets) {
      D <- matrix(0, ncol=nCoVarSets[i], nrow=nCoVarSets[i])
      j <- 1
      for (k in mmcomponents) {
        if (effect[k] > DeterministicEffect && effect[k] < RemainingError) {
          idx <-  startCoVar[i, k] : endCoVar[i, k]

          if (effect[k] == RVarEffect) {
            D[idx, idx] <- Sinv[[i]][[k]] / sigma2[j]
            j = j + 1
          } else if (effect[k] > RVarEffect) {
            D[idx, idx] <- Sinv[[i]][[k]]
          }
        }
      }
       
      rS <- 0
      for (m in 1:length(idx.repet[[i]])) {
        Spalte <- crossprod(Xges[[i]][[m]], Sinv[[i]][[idx.error]][[m]]) # all k
        for (k in 1:nrow(nCoVar)) if (nCoVar[i, k] > 0) {
          idx <- startCoVar[i, k] : endCoVar[i, k]
                    
          D[, idx] <- D[, idx] +
            Spalte %*% Xges[[i]][[m]][, idx, drop=FALSE] * idx.repet[[i]][m]
        }
        rS <- rS + Spalte %*% sumdata[[i]][[m]]
      }      
      z[[i]] <- try(solve(D, rS))
      if (!is.numeric(z[[i]]))
        stop("Design matrix shows linear dependence. Check the design matrix and \n  make sure that you do not use `trend' and the `mixed' model at the same time.")
    }
 
    if (return.z) {
      return(z)
    }

    cumsq <- j <- 0
    for (k in mmcomponents) {
      s2est <- 0
      if (effect[k] == RVarEffect) {
        idx <- startCoVar[i, k] : endCoVar[i, k]
        for (i in 1:sets) {
          zi <- z[idx, i]
          s2est <- s2est + crossprod(zi[idx], Sinv[[i]][[k]] %*% zi[idx])
        }
        cumsq <-  cumsq +  (sigma2[j] - s2est / nTotalComp[k])^2
        j <- j + 1
      }
    }

    if (cumsq < LINEARPART) {
      assign("LINEARPART", cumsq, envir=ENVIR)
      assign("LINEARPARTZ", z, envir=ENVIR)
      assign("LINEARPARTS2", sigma2, envir=ENVIR)
    }    
    return(cumsq)
  }

  
          ## note "global" variables MLEMAX, MLEPARAM
  MLtarget <- function(variab) {
    if (length(variab)>0) variab <- variab + 0  ## unbedingt einfuegen, da bei R Fehler der Referenzierung !! 16.2.10

    if (n.variab > 0) {
      if (PrintLevel>4)  Print(format(variab, dig=20))
      if (any((variab < MLELB) | (variab > MLEUB))) {
        ## for safety -- should not happen, older versions of the optimiser
        ## did not stick precisely to the given bounds
        ## 23.12.03 : still happens
        if (PrintLevel>2) WarningMessage(variab, MLELB, MLEUB, "MLE")   
        penalty <- variab
        variab <- pmax(MLELB, pmin(MLEUB, variab)) 
        penalty <-  - sum(variab - penalty)^2 ## not the best ....
        res <- MLtarget(variab)

        return(res + penalty * (1+ abs(res)))
      }

      param <- as.double(transform(variab[1:nvarWithoutBC]))
      .C("PutValuesAtNA", param, PACKAGE="RandomFields", DUP=FALSE)
      options(show.error.messages = show.error.message)
    }
    if (PrintLevel > 5) {
      cat("\n\nAufruf von MLtarget\n===================\n")
      Print(GetModelInfo(modelreg=ModelNr, level=0))
    }
        
    logdet <- ML.base # Konstante; wird mit -0.5 multipliziert
    Werte <- list() ## nur bei BoxCox Schaetzung verwendet

    stopifnot(is.finite(logdet))
    
    for (i in 1:sets) {
      S <- double((lc[i] * vdim)^2)
      for (k in allcomponents) {

        
        if (effect[k] >= SpaceEffect) { ## ansonsten muss schon vorher gesetzt
          ##                              werden

          .C("CovMatrixMLE", Xdistances[[i]], as.integer(is.dist),
             as.integer(lc[i]),
             as.integer(i),
             as.integer(if (anyeffect) k else 0),  
             S,
             PACKAGE="RandomFields", DUP=FALSE)
          dim(S) <- rep(lc[i] * vdim, 2)         

          if (effect[k] == SpaceEffect) {
            Si <- try(chol(S), silent=silent)

            
            logdet <- logdet + 2 * sum(log(diag(Si)))
            Sinv[[i]][[k]] <<- chol2inv(Si, LINPACK=TRUE)
          } else {
            for (m in 1:length(idx.repet[[i]])) {
              Si <- S[idx.na[[i]][[1]][, m], idx.na[[i]][[1]][, m]] 
              Si <- try(chol(Si), silent=silent)
              
              if (!is.numeric(Si) || any(is.na(Si))) {
                
                assign("MLEINF", TRUE, envir=ENVIR)

                if (FALSE) {
                  Si <- S[idx.na[[i]][[1]][, m], idx.na[[i]][[1]][, m]]
                  for (u in 1:lc[i]) {
                    minor <- det(Si[1:u,1:u, drop=FALSE])
                    if (PrintLevel>5) Print(minor) #
                    if (minor < -1e-10) stop("Covmatrix not positive definite")
                  }
                }
                if (PrintLevel>1 || is.na(MLEVARIAB)) {
                  cat("\nMLE: error in cholesky decomp. -- matrix pos def?")
                }
                return(1E300)
              }
              
              logdet <- logdet + 2 * sum(log(diag(Si))) * idx.repet[[i]][m]
              
              ## Si ist hier Si^{1/2} !! Zum Schluss wird mit -0.5 multipliziert

              Sinv[[i]][[k]][[m]] <<- chol2inv(Si, LINPACK=TRUE)
            }
          }
        }
      }
      S <- NULL
      
      Werte[[i]] <- list()
      for (m in 1:length(idx.repet[[i]])) {
        if (lambdaest) {
          
          Werte[[i]][[m]] <-
            BoxCox(werte[[i]] [ idx.na[[i]][[1]][, m], , idx.na[[i]][[2]][m, ],
                               drop = FALSE ], variab[n.variab])
          for (k in det.effect) {
            if (length(modelinfo$sub[[k]]$param$X)==1)
              Werte[[i]][[m]] <- Werte[[i]][[m]] -
                modelinfo$sub[[k]]$param$X[[i]]
            else 
             Werte[[i]][[m]] <- Werte[[i]][[m]] -
               modelinfo$sub[[k]]$param$X[[i]][idx.na[[i]][[1]][, m]]
            }
        } else {
          Werte[[i]][[m]] <-
            werte[[i]] [ idx.na[[i]][[1]][, m], , idx.na[[i]][[2]][m, ],
                        drop = FALSE ]
        }

        sumdata[[i]][[m]] <<- as.vector(t(apply(Werte[[i]][[m]], 1:2, sum)))
        d <- dim(Werte[[i]][[m]])
        dim(Werte[[i]][[m]]) <- c(d[1] * d[2], d[3])
      } # for m
    } # for i

    if (nCoVarAll > 0) {
      assign("LINEARPART", Inf, envir=ENVIR)
      if (!solvesigma) {
        ## !solvesigma: alle Parameter der Kovarianzen auf
        ## globaler Ebene geschaetzt       
        linearpart(LINEARPARTS2)     
      } else if (length(mixed.idx) == 1) {
        try(optimize(linearpart,
                     lower = LINEARPARTS2 / 10,
                     upper = LINEARPARTS2 * 10),
            silent=silent)
      } else {
        try(optim(LINEARPARTS2, linearpart, method ="L-BFGS-B",
                  lower = LINEARPARTS2 / 10,
                  upper = LINEARPARTS2 * 10,
                  control= lsq.optim.control),
            silent=silent)
      }
    }
    
    quadratic <- 0
    quad.effect <- rep(0, max(1, length(mmcomponents)))
    MLtargetV <- list()
    for (i in 1:sets) {
      MLtargetV[[i]] <- list()
      for (m in 1:length(idx.repet[[i]])) {
        
        MLtargetV[[i]][[m]] <- Werte[[i]][[m]] -
          if (nCoVarAll>0)
            as.vector(Xges[[i]][[m]] %*% LINEARPARTZ[[i]])
          else 0
                                        # y-\sum A_i z_i
        quadratic <- quadratic +
          sum(MLtargetV[[i]][[m]] *
              (Sinv[[i]][[idx.error]][[m]] %*% MLtargetV[[i]][[m]]))       
      }
      
      for (k in mmcomponents) if (effect[k] >= RandomEffect) {
        idx <- startCoVar[i, k] : endCoVar[i, k]
        quad.effect[k] <- quad.effect[k] +
          sum(LINEARPARTZ[[i]][idx] *
              (Sinv[[i]][[k]] %*% LINEARPARTZ[[i]][idx]))
      }
    }

    if (solvesigma) {
      idx <- effect == RVarEffect
      nc <- rowSums(nCoVar[, idx])
      quad.effect[idx] <- log(quad.effect[idx] / nc) * nc
    }

     res <- -0.5 * (logdet + sum(quad.effect) +
                    if (globalsigma) ntotdata * log(quadratic) else quadratic)
    
    if (is.na(res)) {
      warning("NA results appeared")
    }
    
    if (PrintLevel > 6) {
      Print(globalsigma, variab, param, var.global, ML.totalRV,
            quadratic / ML.totalRV)
      readline(paste(res, "Bitte return druecken."))
    }

    if (res > MLEMAX) {
      if (varnugNA) { ## then globalsigma is true as well. So never change oder of
        ##               if's
        sill <- quadratic / ML.totalRV ### check!!!
        param[var.global] <- sill * (1.0 - param[nugget.idx])        
        param[nugget.idx] <- sill * param[nugget.idx]
      } else{
        if (globalsigma)
          param[var.global] <- quadratic / ML.totalRV  ### check!!!
      }

      assign("MLEMAX", res, envir=ENVIR)
      assign("ML.RESIDUALS", MLtargetV, envir=ENVIR)
      if (n.variab > 0)  {
        assign("MLEPARAM", param, envir=ENVIR)
        assign("MLEVARIAB", variab, envir=ENVIR)
      }     
    }

    if (PrintLevel > 3) Print(logdet, quadratic, res)
    
    return(res)
  } # mltarget


  

  
  REMLtarget <- function(variab) {
    variab <- variab + 0  ## unbedingt einfuegen, da bei R Fehler der Referenzierung !! 16.2.10
    if (n.variab > 0) {
      if (PrintLevel>4) {Print(format(variab, dig=10))} #
      if (n.variab == 0) return(NA)
      if (any((variab < MLELB) | (variab > MLEUB))) {
        ## for safety -- should not happen, older versions of the optimiser
        ## did not stick precisely to the given bounds
        ## 23.12.03 : still happens
        if (PrintLevel>2) WarningMessage(variab, MLELB, MLEUB, "REML")
         penalty <- variab 
        variab <- pmax(MLELB, pmin(MLEUB, variab)) 
        penalty <-  - sum(variab - penalty)^2 ## not the best ....
        res <- REMLtarget(variab)
        return(res + penalty * (1+ abs(res)))
      }

      param <- as.double(transform(variab[1:nvarWithoutBC]))

      .C("PutValuesAtNA", param, PACKAGE="RandomFields", DUP=FALSE)
    }

    logdet <- numeric(length(coord))
    quadratic <- numeric(length(coord))

    options(show.error.messages = show.error.message)
    if(length(coord) >1) warning("REML not programmed yet")
    for (i in 1:length(coord)) {
      cov.matrix <- double((vdim*lc[i])^2)

      .C("CovMatrixMLE", Xdistances[[i]], as.integer(is.dist),
         as.integer(lc[i]), as.integer(i-1),
         as.integer(length(mmcomponents)), cov.matrix, 
         PACKAGE="RandomFields", DUP=FALSE)

       dim(cov.matrix) <- rep(vdim*lc[i], 2)

       if (any(eigen(cov.matrix)$val <0) && PrintLevel>1) {
         Print(param, eigen(cov.matrix)$val) #
      }

     cov.matrix <-
        try(chol(crossprod(REML.A[[i]], cov.matrix %*% REML.A[[i]])),
            silent=silent)
      options(show.error.messages = TRUE)
      if (!is.numeric(cov.matrix) || any(is.na(cov.matrix))) {
        if (PrintLevel>1 || is.na(MLEVARIAB)) {
          cat("\nerror in cholesky decomposition -- matrix pos def?")
          Print( format(variab, dig=20), format(param, dig=20), cov.matrix)
        }
        return(1E300)
      }
      stopifnot(all(diag(cov.matrix) > 0))

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
      logdet[i] <- 2 * sum(log(diag(cov.matrix)))# log(det C)
      if (!is.finite(logdet[i])) logdet[i] <- 1E300 ## the best ?! 
      cov.matrix <- chol2inv(cov.matrix, LINPACK=TRUE) # La.chol2inv, LIN=TRUE
      if (lambdaest)  {
          REML.data.transform <- crossprod(REML.A[[i]], BoxCox(as.matrix(werte[[i]]), variab[n.variab]))
	  quadratic[i] <- sum(REML.data.transform * (cov.matrix %*% REML.data.transform))
      }
      else
      	  quadratic[i] <- sum(REML.data[[i]] * (cov.matrix %*% REML.data[[i]]))
    }

    ##               sum_i (D_i - Xm)^T C^{-1} (D_i - X m)
    if (globalsigma) {
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
      res <- REML.loglcrepet + sum(repet *logdet + REML.totalRV * log(quadratic))
    } else {
      res <- sum(repet * logdet + quadratic)  ## 23.9.08 Marco
      ##res <- sum(repet * (logdet + quadratic)) 
    }
    res <- -0.5 * (res + REML.twopilcrepet)
    if (is.na(res) || is.na(MLEMAX) #|| debug
        ) {
      filename <- "RandomFields.fitvario.reml.bug.rda"
      txt <- paste("algorithm has failed for unknown reasons -- please contact the author; please send the data that have been written on the file '",
             filename, "'.", sep="")
      fitvario.bug <- NULL
      if (is.na(res) || is.na(MLEMAX)) stop(txt)
    }

    
    if (res > MLEMAX) {
      if (varnugNA) {
        sill <- sum(quadratic) / sum(REML.totalRV)
        param[var.global] <- sill * (1.0 - param[nugget.idx])
        param[nugget.idx] <- sill * param[nugget.idx]
      } else {
        if (globalsigma) param[var.global] <- quadratic / REML.totalRV
      }
      assign("MLEMAX", res, envir=ENVIR)
      if (n.variab > 0) {
        assign("MLEPARAM", param, envir=ENVIR)
        assign("MLEVARIAB", variab, envir=ENVIR)
      }
    }
    if (PrintLevel>4) cat("result REML",res,"\n")
    return(res)
  }



  
  crosssettings <- function(M) {
    stopifnot(is.numeric(trend) && length(trend)==1)

    stop("die unterscheidung funktioniert nicht mehr!!")
    assign("CROSS.KRIGE", if (is.na(pm$mean)) "O" else "S",# ordinary/simple
           envir=ENVIR)
    switch(M,
           "cross.sq" = {
              assign("CROSS.VAR", FALSE, envir=ENVIR)
              ## d can have more than 1 element
              CROSS.DIST <- function(x, d) sum((x-d)^2)
              environment(CROSS.DIST) <-EmptyEnv()
            },
           "cross.abs" = {
              assign("CROSS.VAR", FALSE, envir=ENVIR)
              ## d can have more than 1 element
              CROSS.DIST <- function(x, d) sum(abs(x-d))
              environment(CROSS.DIST) <- EmptyEnv()
           },
           "cross.ign" = {
             assign("CROSS.VAR", TRUE, envir=ENVIR)
             CROSS.DIST <- function(x, d) {
               sum(0.5 * log(2 * pi * x$var) + (d - x$est)^2/ (2 * x$var))
             }
             environment(CROSS.DIST) <- EmptyEnv()
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
             environment(CROSS.DIST) <- EmptyEnv()
           },
           stop(M, " unknown")
           )
    assign("CROSS.DIST", CROSS.DIST, envir=ENVIR)
  }


  crosstarget <- function(variab) {
    variab <- variab + 0 ## unbedingt einfuegen, da bei R Fehler der Referenzierung !! 16.2.10
    stop("cross not programmed")

    if (PrintLevel>4) {Print(format(variab, dig=10))} #
    if (any((variab<CROSSLB) | (variab>CROSSUB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 23.12.03 : still happens
      if (PrintLevel>2) WarningMessage(variab, CROSSLB, CROSSUB, "Cross")
      penalty <- variab 
      variab <- pmax(CROSSLB, pmin(CROSSUB, variab)) 
      penalty <- + sum(variab-penalty)^2 ## not the best ....
      res <- crosstarget(variab)
      return(res + penalty * (1+ abs(res)))
    }


    param <- as.double(transform(variab[1:nvarWithoutBC]))
    model <- pm

    if (nCoVarAll > 0) {
      stopifnot(is.numeric(trend) && length(trend)==1)
      model$mean <- variab[nCROSSINDEX + 1:nCoVarAll]
      model$trend <- NULL
      stopifnot(is.numeric(trend) && length(trend)==1) ## not programmed yet,
      ##                           currently only nCoVarAll==1
      param[CROSSINDEX] <- variab[1:nCROSSINDEX] # variabX ?!
    } else param[CROSSINDEX] <- variab # variabX ?!
    model$param <- transform(param)

    res <- 0.0
    ## data is matrix!

    stop("internal Kriging not programmed yet!")

    i <- NA ### neu zu schreiben

    for (d in (1:lc[i])) {
      res <- res +
        CROSS.DIST(Kriging(CROSS.KRIGE, x=coord[d, , drop=FALSE], grid=FALSE,
                           model=model, given=coord[-d, , drop=FALSE],
                           data=werte[-d, , drop=FALSE], trend = trend, pch="",
                           return.variance=CROSS.VAR),
                   werte[d, ])
    }
    res <- res / CROSS.lcrepet
    if (res<CROSSMIN) {
      assign("CROSSMIN", res, envir=ENVIR)
      assign("CROSSMODEL", model, envir=ENVIR)
    }
    if (PrintLevel>6) cat("result cross val", res, "\n")
    return(res)
  }

  get.covariates <- function(Variab) {
    if (nCoVarAll == 0) return(NULL)
    max <- MLEMAX
    param <- MLEPARAM
    variab <- MLEVARIAB
    inf <- MLEINF
    resid <- ML.RESIDUALS
    linpart <- LINEARPART
    linz <- LINEARPARTZ
    lins2 <- LINEARPARTS2
    
    assign("MLEMAX", -Inf, envir=ENVIR)
    assign("LINEARPART", Inf, envir=ENVIR)
    MLEsettings("ml")
    MLtarget(Variab)
    result <- unlist(LINEARPARTZ)
    
    assign("MLEMAX", max, envir=ENVIR)
    assign("ML.RESIDUALS", resid, envir=ENVIR)
    assign("MLEPARAM", param, envir=ENVIR)
    assign("MLEVARIAB", variab, envir=ENVIR)
    assign("MLEINF", inf, envir=ENVIR)
    assign("LINEARPART", linpart, envir=ENVIR)
    assign("LINEARPARTZ", linz, envir=ENVIR)
    assign("LINEARPARTS2", lins2, envir=ENVIR)
    return(result)
  }
     
  IDX <- function(name) {idx <- tblidx[[name]]; idx[1]:idx[2]}


  ## to avoid warning on "no visible binding" we define the following
  ## variable that are used in the local functions:
  ENVIR <- environment()
  LSQ.SELF.WEIGHING <- LSQ.WEIGHTS <- LSQ.BINNEDSQUARE <- 
    REML.A <- REML.data <- REML.loglcrepet <- REML.totalRV <-REML.twopilcrepet <-
      ML.base <- ML.totalRV <- MLEtarget <-
        ML.RESIDUALS <- MLEMAX <- MLEVARIAB <- MLEINF <- MLEPARAM <- 
      CROSS.DIST <- CROSS.KRIGE <- CROSS.VAR <- CROSSMODEL <-
        LINEARPART <- LINEARPARTS2 <- LINEARPARTZ  <- 
          NULL       

######################################################################
###              End of definitions of local functions             ###
######################################################################
  
######################################################################
###    Initial settings, coord part I (without model info)         ###
######################################################################

  if (PrintLevel>1) cat("\ninitial settings...\n")

  
  if (!missing(x) && is.data.frame(x)) x <- as.matrix(x)
  if (is.data.frame(data)) data <- as.matrix(data)
  
    
  if (xgiven <- is.null(Distances)) {
    stopifnot(missing(truedim))
    if (is.list(x)) {
      if (!is.null(y) & !is.list(y) |
          !is.null(z) & !is.list(z) |
          !is.null(T) & !is.list(T) |
          !is.list(data))
        stop("either all coordinates and data must be given by a list or none.")
    } else {
       x <- list(x)
       if (!is.null(y)) {
         stopifnot(!is.list(y))
         y <- list(y)
       }
       if (!is.null(z)) {
         stopifnot(!is.list(z))
         z <- list(z)
       }
       if (!is.null(T)) {
         stopifnot(!is.list(z))
         T <- list(T)
       }
       stopifnot(!is.list(data))
       data <- list(data)
    }

    j <- 1
    coord <- list()
    for (i in 1:length(x)) {
      neu <- CheckXT(x[[i]], y[[i]], z[[i]], T[[i]], grid, gridtriple)
      if (PrintLevel>2) Print(neu)
      if (i==1) {
        spacedim <- neu$spacedim
        xdim <- truedim <- as.integer(spacedim + !is.null(T))
        storage.mode(truedim) <- "integer"
        
      } else if (spacedim != neu$spacedim |
                 truedim != as.integer(spacedim + !is.null(T)))
        stop("dimensions of all sets must be the same")
      coord[[i]] <- neu$x
      neu <- NULL
      
      if (grid) {
        ## note: checkxt returns always gridtriple format !
        ## ToDo: expand.grid might be avoided as in GaussRF
        coord[[i]] <- apply(coord[[i]], 2, function(z) seq(z[1], z[2], z[3]))
        args <- lapply(apply(coord[[i]], 2, list),
                       function(z) c(z, recursive=TRUE))
        coord[[i]] <- as.matrix(do.call("expand.grid", args=args))
      }
      if (!is.null(T)) {
        tt <- seq(T[[i]][1], T[[i]][2], T[[i]][3])
        coord[[i]] <- cbind(coord[[i]][rep(1:nrow(coord[[i]]), length(tt)), ],
                        rep(tt, each=nrow(coord[[i]])))
      }
      coord[[i]] <- t(coord[[i]]) ## !!!!!
    }  
    neu <- x <- y <- z <- T <- NULL    
    lc <- sapply(coord, ncol)
  } else { ## x coordinates, not distances
    if (!is.list(Distances)) {
      Distances <- list(Distances)
      if (is.list(data))
        stop("if list of data is given then also for Distances ")
      data <- list(data)
    } else if (!is.list(data)) {
      stop("if list of Distances is given then also for data ")
      if (length(data) != length(Distances))
          stop("length of Distances does not match length of data")
    }
   
    
    for (i in 1:length(Distances)) {
      if (any(is.na(data)))
          stop("missing data are not allowed if Distances are used.")
    }
    coord <- Distances
    stopifnot(missing(x) || is.null(x), is.null(y), is.null(z))
    xdim <- as.integer(1)
    spacedim <- truedim <- as.integer(truedim)
    lcc <- sapply(Distances, function(x) 0.5 * (1 + sqrt(1 + 8 * length(x))) )
    lc <- as.integer(lcc)
    
    if (any(lc != lcc))
      stop("number of distances not of form k(k-1)/2")
  
  }
  sets <- length(lc)

  ## check box-cox transformation 
  if (lambdaest <- !missing(BC.lambda) && is.na(BC.lambda)) {
    if (any(unlist(data) <= 0))   warning("negative data may cause errors")
    if (BC.lambdaLB > 1 || BC.lambdaUB<1)
      stop("bounds for Box-Cox lambda should satisfy lambdaLB <= 1 <= lambdaUB")
  }



################    analyses of orginal model        ###############
##### variables needed for analysis of trend, upper and lower input --
##### user cannot know what the internal represenatation is
  
  pm <- PrepareModel(model, param, nugget.remove=FALSE)

  ResMLEGet <- .Call("MLEGetModelInfo", pm$model, truedim, xdim,
                     PACKAGE="RandomFields")    
  orig.anyeffect <- length(ResMLEGet$effect) > 0
  stationary <- ResMLEGet$stationary ## note: only with respect to the
  ##                                    coordinates !!
  isotropy <- ResMLEGet$isotropy
  stopifnot(xgiven || stationary && isotropy)

  vdim <- GetModelInfo(modelreg=ModelNr, level=1)$vdim
  
  minmax <- ResMLEGet$minmax # 4 Spalten: 1:min, 2:max, 3:type, 4:is.nan, not na
  ##VARPARAM, SCALEPARAM, DIAGPARAM, ANISOPARAM, // 0..3
  ##INTERNALPARAM, ANYPARAM, TRENDPARAM, NUGGETVAR, MIXEDVAR, // 4..8
  ## REGRESSION // 9
 # mm.attr <- attr(minmax, "row.names")
#  idx <- minmax[, 3] != REGRESSION
#  minmax <- minmax[idx, , drop=FALSE]
#  attr(minmax, "row.names") <- mm.attr[idx]
  

  ncovparam <- nrow(minmax)
  ptype <- minmax[, 3]

 ################    analyse trend and include it in model    ###############
  
  ## note: only vdim is needed here. MLEGetModelInfo and GetModelInfo will be
  ## called again with a slightly different model
  
  nCoVar <- 0
  if (PrintLevel>1) cat("\npreparing trend ...")
  
  CoVarblockmatrix <- function(n, A, B) {
     if (is.logical(A) && is.na(A) && is.logical(B) && is.na(B))
       NA
     else if (is.logical(A) && is.na(A))
       rbind(matrix(0, nrow=nrow(B)*(n-1), ncol=ncol(B)), B)
     else if (is.logical(B) && is.na(B))
       rbind(A, matrix(0, nrow=nrow(A), ncol=ncol(A)))
     else
       rbind(cbind(A, matrix(0, nrow=nrow(A), ncol=ncol(B))),
             cbind(matrix(0, nrow=nrow(B), ncol=ncol(A)), B)
             )
   }
  
  if (!missing(param)) {
    if (!is.null(trend)) stop("param and trend cannot be given simultaneously.")
    trend <- pm$trend 
    if (length(trend) > 0) param <- param[-1] # delete mean
  }
  
  if (is.list(trend)) {  ### list
    PreCoVariates <-
      eval(parse(text=paste("list(",
                   paste(rep(NA, times=vdim * sets), collapse=","), ")", sep="")))
    if (is.matrix(trend[[1]]) && length(trend) == sets * vdim) {
      ## list of matrices: 
      ## 1st realisation of 1st component, 1st real. of 2nd comp, ...,
      ##                                             2nd real. of 1st comp,...
      stopifnot(all(sapply(trend, is.matrix)))
      PreCoVariates <- trend
    } else {  ### list of univariate trends
      stopifnot(length(trend) == vdim)
      for (k in 1:vdim) {
        if (is.list(trend[[k]])) {  ### list
          if(length(trend[[k]]) == sets)  { ## list of matrices
            stopifnot(all(sapply(trend[[k]], is.matrix)))
            for (i in 1:sets)
              PreCoVariates[[vdim*(i-1)+k]] <- trend[[k]][[i]]
          } else if(length(trend[[k]]) == sets + 1){# list of matrices & formula
            if (!xgiven) stop("distance matrix cannot be passed along trends")
            trend_form <- trend[[k]][[sets+1]]
            trend_mat <- vector("list", sets)
            for(i in 1:length(trend_mat))
              trend_mat[[i]] <- trend[[k]][[i]]
            if (class(trend_form) !=  "formula")
              stop(paste("Last entry of trend should be a formula",
                         "if length(trend) > number of repetitions."))
            if (is.null(T)) time.name <- NULL
            CoVariates_form <-
              lapply(coord, function(x) formulaToMatrix(trend_form, x,
                                                        var.name=var.name,
                                                        time.name=time.name))
            if (!all(sapply(trend_mat, ncol) == ncol(trend_mat[[1]]) ))
              stop("CoVariates do not have all the same dimension.")
            stopifnot(all(sapply(trend_mat, is.matrix)))
            for (i in 1:length(trend_mat))
              PreCoVariates[[vdim*(i-1)+k]] <- cbind(trend_mat[[i]],
                                                     CoVariates_form[[i]])
          } else stop("length of trend should equal the # of repetitions (+1)")
        } else {  # no list of lists
          if (any( (is.logical(trend[[k]]) || is.numeric(trend[[k]]))
                  && is.na(trend[[k]]))  )     {	### NA
            if (length(trend[[k]]) == 1) {
              if (PrintLevel > 0)
                cat("fitvario assumes that the mean (of component no.",
                    k, ") has to be estimated.\n")
              for (i in 1:sets)
                PreCoVariates[[vdim*(i-1)+k]] <- matrix(1, nrow=lc[i], ncol=1)
            } else stop("NAs are only allowed in list mode")
          } else if (length(trend[[k]]) == 1 && is.numeric(trend[[k]])) {
### number
            if (PrintLevel > 0)
              cat("fitvario assumes that the mean is constant in",
                  "component no.", k,"and equals", trend[[k]],"\n")
            ##lapply(sumdata, function(x) x[,k] <- x[,k] - trend[[k]])
            ##lapply(werte, function(x) x[,seq(k, ncol(x), vdim)]
            ##       <- x[,seq(k, ncol(x), vdim)] - trend[[k]])
            trend[[k]] <- 0
          } else if (is.null(trend[[k]])) { ### NULL
            trend[[k]] <- 0			
          } else if (is.matrix(trend[[k]])) {  ### matrix
            if (sets > 1 && PrintLevel > 0)
              cat("it is assumed that the trend matrix (of component ",
                  "no. ", k,") is the same for each measurement;",
                  " otherwise use a list of matrices as trend\n")
            for (i in 1:sets)
              PreCoVariates[[vdim*(i-1)+k]] <- trend[[k]]
          } else if (class(trend[[k]])=="formula") { ### formula
            if (is.null(T)) time.name <- NULL
            for(i in 1:sets)
              PreCoVariates[[vdim*(i-1)+k]] <-
                formulaToMatrix(trend[[k]], coord[[i]],
                                var.name=var.name, time.name=time.name)
          } else stop("incorrect trend")
        }
      }
    }
    
    CoVariates <- eval(parse(text=paste("list(", paste(rep(NA, times=sets),
                               collapse=","), ")",sep="")))
    for (i in 1:sets) {
      for (k in 1:vdim) {
        CoVariates[[i]] <-
          CoVarblockmatrix(k, CoVariates[[i]], PreCoVariates[[vdim*(i-1)+k]])
      }
      if (!(is.logical(CoVariates[[i]]) && is.na(CoVariates[[i]]))) {
        stopifnot(as.integer(nrow(CoVariates[[i]])/ncol(coord[[i]])) ==
                  nrow(CoVariates[[i]])/ncol(coord[[i]]))
        idx <- is.na(rowSums(matrix(rowSums(CoVariates[[i]]), 
                                    nrow=ncol(coord[[i]]))))
        if (any(idx)) {
          stop("there are NA Covariates") ## 31.12.10: warning 
          CoVariates[[i]] <- CoVariates[[i]][rep(!idx, times=vdim), ,
                                             drop=FALSE]
          coord[[i]] <- coord[[i]][, !idx, drop=FALSE]
          werte[[i]] <- werte[[i]][!idx, , drop=FALSE]
          sumdata[[i]] <- sumdata[[i]][!idx, ,drop=FALSE]
        }
      }
    }
   ## if (!(is.logical(CoVariates[[1]]) && is.na(CoVariates[[1]]))) {
   ##   if (sum(abs(sapply(CoVariates, ncol) - ncol(CoVariates[[1]]))) == 0)
   ##     nCoVar = ncol(CoVariates[[1]])
   ## else nCoVar <- sapply(CoVariates, ncol)
   ## }
      ## 5.4.10:
      ##nCoVar <- sapply(CoVariates, ncol)
   
  } else { ## no list => univariate model
    if (vdim!=1)
      stop("Trend definition does not match multivariate covariance model.")
    
    if (any( (is.logical(trend) || is.numeric(trend)) && is.na(trend))) {# NA
      if (length(trend) == 1) {
        if (PrintLevel > 0)
          cat("fitvario assumes that the mean must be estimated.\n")
        CoVariates <- lapply(coord, function(x) matrix(1, nrow=lc[1], ncol=1))
        nCoVar <- 1
      }
      else stop("NAs are only allowed in list mode.")
    } else if (length(trend) == 1 && is.numeric(trend)) {   ### number
      if (PrintLevel > 0)
        cat("fitvario assumes that the mean is constant and equals", trend,".\n")
      ## lapply(sumdata, function(x) x - trend)
      ## lapply(werte, function(x) x - trend)       
      ## trend <- 0
    } else if (!is.null(trend)) { ### NULL
      if (is.matrix(trend)) { ### matrix
        if (sets > 1 && PrintLevel > 0)
          cat("It is assumed that the trend matrix is the same for each",
              "measurement; otherwise use a list of matrices as trend.\n")
        CoVariates <- vector("list", sets)
        for (i in 1:sets)
          CoVariates[[i]] <- trend[, , drop=FALSE]
        ##if (sum(abs(sapply(CoVariates, ncol) - ncol(CoVariates[[1]]))) == 0)
        ##  nCoVar = ncol(CoVariates[[1]])
        ##else nCoVar <- sapply(CoVariates, ncol)
        ## 5.4.10; rueck 7.1.11
        nCoVar <- sapply(CoVariates, ncol)
      } else if (class(trend)=="formula") { ### formula
        if (is.null(T)) time.name <- NULL
        CoVariates <- lapply(coord, function(x)
                             formulaToMatrix(trend, x, var.name=var.name,
                                             time.name=time.name))
        nCoVar <- sapply(CoVariates, ncol)
      } else stop("incorrect trend")
      for(i in 1:sets) {
        if (any(nCoVar>0) && any(idx <- is.na(CoVariates[[i]]))) {
          stop("there are NA Covariates") # 31.12.10: warning
          CoVariates[[i]] <- CoVariates[[i]][!idx, , drop=FALSE]
          coord[[i]] <- coord[[i]][, !idx, drop=FALSE]
          werte[[i]] <- werte[[i]][!idx, , drop=FALSE]
          sumdata[[i]] <- sumdata[[i]][!idx, ,drop=FALSE]
        }
      }
    }
  }



######################################################################
## model specific settings
######################################################################

  
#######################   upper,lower,user     ########################
  if (PrintLevel>1) cat("\n lower and upper ...\n")

  users.lower <- users.upper <- NULL
  if (missing(param))  {
    if (is.list(lower))  {
      pm.l <- PrepareModel(lower, nugget.remove=FALSE)
      
      users.lower <- lower <-
        .Call("Take2ndAtNaOf1st", pm$model, pm.l$model,
               truedim, xdim, stationary, ncovparam, PACKAGE="RandomFields")
    } else if (is.null(transform)) stopifnot(is.null(lower))
    if (is.list(upper))  {
      pm.u <- PrepareModel(upper, nugget.remove=FALSE)      
      users.upper <- upper <-
        .Call("Take2ndAtNaOf1st", pm$model, pm.u$model,
               truedim, xdim, stationary, ncovparam, PACKAGE="RandomFields")
    } else if (is.null(transform)) stopifnot(is.null(upper))
  } else {
    if (!is.null(lower)) {
      if (length(lower) != ncovparam)
        stop("length of lower does not match the number of parameters to be estimated")
      users.lower <- lower
    }
    if (!is.null(upper)) {
      if (length(upper) != ncovparam)
        stop("length of lower does not match the number of parameters to be estimated")
      users.upper <- upper
    }
  }
  if(!is.null(users.guess))  {
    if (!missing(param)) 
      pm.u <- PrepareModel(model=model, users.guess, nugget.remove=FALSE)
    else 
      pm.u <- PrepareModel(model=users.guess, nugget.remove=FALSE)

#    Print(ncovparam)
    
    users.guess <-
      .Call("Take2ndAtNaOf1st", pm$model, pm.u$model,
            truedim, xdim, stationary, ncovparam, PACKAGE="RandomFields")
  }


#######################   include trend        ########################

  if (trend.input <- (nCoVar > 0)) {
    if (nCoVar == 1 && is.numeric(trend) && !is.na(trend)) {
      mm <- list("mixed", X=CoVariates, b=1)
    } else {
      mm <- list("mixed", X=CoVariates, b=rep(NA, ncol(CoVariates[[1]])))
      ## b is always of the same length
    }
    if (orig.anyeffect) {
      if (PrintLevel>0) cat("The model is considered as a mixed effect model\n")
      trend.input <- FALSE
      pm$model[[length(pm$model) + 1]] <- mm
    } else {
      pm$model <- list("+", mm, pm$model)
      if (PrintLevel>0) {
        if (.C("MLEanymixed", anymixed=integer(1))$anymixed)
          cat("The model is not considered as a mixed effect model although `mixed` appears as a submodel\n")
        else
          cat("The model is not a mixed effect model.\n")
      }
     }    
    rm("CoVariates")
    if (PrintLevel >2) Print(pm)
  }
  
  
#######################         model          ########################

  ResMLEGet <- .Call("MLEGetModelInfo", pm$model, truedim, xdim,
                     PACKAGE="RandomFields")

   
  NAs <-  ResMLEGet$NAs
  effect <- ResMLEGet$effect

  anyeffect <- length(effect) > 0
  if (anyeffect) {
    idx <- effect == RemainingError    
    if (sum(idx) != 1)
      stop("a model with fixed or random effects must have exactly one error component")
    allcomponents <- 1:length(effect)
    mmcomponents <- which(!idx)
  } else {
    mmcomponents <- integer(0)
    allcomponents <- 1
    effect <- RemainingError
  }
  idx.error<- which(effect == RemainingError)
  if (!any(effect == RVarEffect)) solvesigma <- FALSE    


  modelinfo <- GetModelInfo(modelreg=ModelNr, level=1)
  vdim <- modelinfo$vdim
  
  xdim <- as.integer(if (isotropy) 1 else truedim + 0)


#  Print(sum(NAs), nrow(minmax), NAs, minmax)
  
  if (anyeffect) stopifnot(sum(NAs) == nrow(minmax))
  if (anyeffect) {
    eff <- rep(FALSE, nrow(minmax))
    csNAs <- cumsum(c(0, NAs))
    for (k in 1:length(effect)) {
      if (effect[k] > FixedEffect && effect[k] <= SpaceEffect && NAs[k] > 0)
        eff[(csNAs[k]+1) : csNAs[k]] <- TRUE
    }
    minmax[ptype == VARPARAM & eff, 1] <- minmixedvar
    minmax[ptype == VARPARAM & eff, 2] <- maxmixedvar
  }


 
###########################      transform     #######################
## either given bu users.transform + users.min, users.max
  delete.idx <- NULL
   if (is.null(transform)) {
    if (any(minmax[,4]==1)) ## nan, not na
      stop("NaN only allowed if transform is given.")
    transform <- function(x) x;
    if (is.null(lower)) lower <- minmax[,1] 
    if (is.null(upper)) upper <- minmax[,2]
    
  } else {
    standard.style <- solvesigma <- FALSE
    stopifnot(!is.null(lower), !is.null(upper) || is.logical(lower),
              missing(param))
    if (is.logical(lower)) {
      stopifnot(length(lower) == nrow(minmax))
      delete.idx <- which(!lower)
      lower <- minmax[,1]
      upper <- minmax[,2]
      if (any(minmax[,4]==1))
        stop("logical lower and NaN may not given at the same time")
    } else {
      stopifnot(is.function(transform))
      stopifnot(!is.null(users.guess))
      ptype <- NULL  ## if NULL, then better bounds are not searched for
    }
  }
  if (!is.null(lower)) {
    stopifnot(!is.null(upper),
              typeof(lower) == typeof(upper),
              length(lower) == length(upper))
  }
  stopifnot(length(lower) == length(upper))



#################################################################
##############     prepare constants in S, X,etc      ###########
#################################################################

##############         distances              #################

  ## note: the direct C call needs matrix where points are given column-wise
  ##       whereas the R function CovarianceFct need them row-wise,
  ##                   except for fctcall==CovarianceMatrix
  ## distances are used to calculate bounds for the scale parameter(s)
  ## to be estimated -- further it is used in MLtarget

###  30.3.06 folgenden code durch nachfolgenden ersetzt.
###     dd <- 0   
###     for (i in 1:truedim) {
###      dx <- outer(coord[,i], coord[,i], "-")
###      distances[i, ] <- as.double(dx[lower.tri(dx)])
###      ## x <- c(1,6,0.1)
###      ## outer(x, x, "-")
###      ##     [,1] [,2] [,3]
###      ##[1,]  0.0 -5.0  0.9
###      ##[2,]  5.0  0.0  5.9
###      ##[3,] -0.9 -5.9  0.0
###      if (i<=spacedim) dd <- dd + distances[i, ]^2 
###    }
###    dd <- sqrt(dd)
###    coord <- NULL
###    maxdistances <- max(dd)
###    mindistances <- min(dd[dd!=0])  ## notwendig falls !is.null(T)
###    rm("dd")
###
        
  Xdistances <- vector("list", sets)
  mindistances <- matrix(ncol=sets, nrow=2)
  for (i in 1:sets) {       
    if (stationary && isotropy) {
      Xdistances[[i]] <-
        if (xgiven) as.vector(dist(t(coord[[i]]))) else Distances[[i]]
      idx <- Xdistances[[i]] < nuggetrange
      if (any(idx)) {
        if (!allowdistanceZero)
          stop("distance with value 0 identified -- use allowdistanceZero=T?")
        Xdistances[[i]][idx] <- nuggetrange
      }

      mindistances[, i] <- range(Xdistances[[i]][!idx])
      
    } else {
      d <- double(truedim * lc[i] * (lc[i] - 1) / 2)
      anyidx <- TRUE
      while (anyidx) {
        .C("vectordist",
           as.double(coord[[i]]),
           as.integer(dim(coord[[i]])),
           d,
           as.integer(FALSE),
           PACKAGE="RandomFields", DUP=FALSE)
        
        dim(d) <- c(truedim, lc[i] * (lc[i] - 1) / 2) 
        absdist <-
          sqrt(colSums(if (is.null(T)) d^2 else 
                       d[1:spacedim, d[truedim, ]==0, drop=FALSE]^2))
        idx <- absdist < nuggetrange
        if (anyidx <- any(idx)) {
          if (!allowdistanceZero) 
            stop("distance with value 0 identified - use allowdistanceZero=TRUE?")
          coord[[i]] <- coord[[i]] *
            (1 + rnorm(length(coord[[i]]), 0, nuggetrange))
        }
      }
      mindistances[, i] <- range(absdist[!idx])
      
      if (PrintLevel ==3) Print(stationary, isotropy)
      
      if (stationary) {
        if (isotropy) {
          Xdistances[[i]] <- absdist
          if (PrintLevel > 8) hist(absdist)
        } else
        Xdistances[[i]] <- d
      } else {
        Xdistances[[i]] <- coord[[i]]
      }
      storage.mode(Xdistances[[i]]) <- "double"
    }
  }
  maxdistances <- max(mindistances[2, ])
  mindistances <- min(mindistances[1, ])
  is.dist <- TRUE ## might be changed in future when fitvario is speeded up
  ##                 by intermediate results in the knots
  
##############         Coordinates & data    #################
  if (PrintLevel>1) cat("\ndistances and data...")
    ## note: the direct C call needs matrix where points are given column-wise
    ##       whereas the R function CovarianceFct need them row-wise,
    ##                   except for fctcall==CovarianceMatrix

  if (lambdaest && vdim > 1) {
      ## ToDo: should be possible in future
    stop("lambda can be estimated only in the univariate case")
  }    

  ntotdata <- sum(sapply(data, length)) 
  idx.na <- werte <- sumdata <- vector("list", length(coord))
  repet <- numeric(length(coord))
  idx.repet <-  vector("list", length(sets))

  det.effect <- which(effect == DeterministicEffect)
  for (i in 1:sets) {
    repet <- length(data[[i]]) /  (lc * vdim)
    
    if (repet != as.integer(repet))
      stop("number of data does not match number of coordinates" )
    werte[[i]] <- array(data[[i]], c(lc, vdim, repet))
   
    for (k in det.effect) {
      size <- nrow(modelinfo$sub[[k]]$param$X[[i]])
      if (size != 1) {
        
        if (size != prod(dim(werte[[i]])[1:2]))
          stop("matrix of data does not match deterministic effect")
      #  dim(modelinfo$sub[[k]]$param$X[[i]]) <- dim(werte[[i]])[2:1]
     #   modelinfo$sub[[k]]$param$X[[i]] <-
      #    t(modelinfo$sub[[k]]$param$X[[i]])
      }
    }

    if (!lambdaest) {
      if (!missing(BC.lambda)) werte[[i]] <- BoxCox(werte[[i]], BC.lambda)
      for (k in det.effect) {
        werte[[i]] <- werte[[i]] - as.vector(modelinfo$sub[[k]]$param$X[[i]])
      }
    }


    isnawerte <- is.na(werte[[i]])
    idx <- apply(isnawerte, 1, any)
    idx.na[[i]] <- sumdata[[i]] <- list()
    if (na.rm || !any(idx)) {
      idx.na[[i]][[1]] <- matrix(ncol=1, nrow=nrow(werte[[i]]), TRUE)
      idx.na[[i]][[2]] <- matrix(rep(TRUE, repet), nrow=1)
    } else { 
      ## ToDo: durch Gruppierung deutlich verbesserbar
      ## Gruppierung wird bereits in den targets verwendet
     idx.na[[i]][[1]] <- isnawerte
     idx.na[[i]][[2]] <- diag(ncol(isnawerte))
       storage.mode(idx.na[[i]][[1]]) <- "logical"
    }
    
    idx.repet[[i]] <- rowSums(idx.na[[i]][[2]])
    
    for (m in 1:length(idx.repet)) {
       sumdata[[i]][[m]] <-
        as.vector(t(apply(werte[[i]][idx.na[[i]][[1]][, m], ,
                                     idx.na[[i]][[2]][m, ],
                                     drop =FALSE],  1:2, sum)))
    }
  }

  if (vdim>1 && PrintLevel>0)
    cat("Due to the covariance model a ", vdim,
        "\b-variate random field is expected. Therefore the data matrix",
        "is assumed to consist of ", repet,
        "independent measurements for each point.",
        "each realisation is given as the entries of ", vdim,
        "consecutive columns.\n")


##############      find upper and lower bounds      #################
  if (PrintLevel>1) cat("\nbounds...")

  txt <- "lower and upper are both lists or vectors of the same length or NULL"
  lengthmismatch <- "lengths of bound vectors do not match model"
  structuremismatch <- "structures of bounds do not match the model"

  var.idx <- which(ptype == VARPARAM)
  MIXED.IDX <- (ptype == MIXEDVAR) & solvesigma
  mixed.idx <-  which(MIXED.IDX) 
  solvesigma <- length(mixed.idx) > 0
  nugget.idx <- which(ptype == NUGGETVAR)
  SCALE.IDX <- ptype == SCALEPARAM  ## large capitals 
  varnames <- minmax.names <- attr(minmax, "row.names")

  vardata <- var(unlist(werte))
  if (vardata==0) stop("data values are identically constant")
  ## autostart will give the starting values for LSQ
  if (is.null(ptype)) {
    ## appears if transform is given. Then better do not search for
    ## automatic bounds
    autostart <- users.guess
    SCALE.IDX <- NUGGET.IDX <- VAR.IDX <- FALSE
  } else {
    autostart <- numeric(nrow(minmax))

    ## now the bounds and starting values are set...
    idx <- c(var.idx, nugget.idx)
    if (length(idx) > 0) {

      ## lower bound of first model is treated differently!
      ## so the "main" model should be given first!             !!!!!
      if (!anyeffect) {
        lower[idx[1]] <- vardata / lowerbound.var.factor
        if (length(idx) > 1) lower[idx[-1]] <- 0
      }

      upper[idx] <- vardata * upperbound.var.factor
      autostart[idx] <- vardata / length(idx)
    }
     
    if (any(idx <- ptype == DIAGPARAM)) {
      lower[idx] <- 1 / (upperbound.scale.factor * maxdistances)
      upper[idx] <- lowerbound.scale.LS.factor / mindistances
      autostart[idx] <- 8 / (maxdistances + 7 * mindistances)
    }

    if (any(idx <- ptype == ANISOPARAM)) {
      lower[idx] <- -lowerbound.scale.LS.factor / mindistances
      upper[idx] <- lowerbound.scale.LS.factor / mindistances
      autostart[idx] <- 0
    }
    
    if (any(SCALE.IDX)) {
      idx <- which(SCALE.IDX)
      lower[idx] <- mindistances / lowerbound.scale.LS.factor
      upper[idx] <- maxdistances * upperbound.scale.factor
      autostart[idx] <- (maxdistances + 7 * mindistances) / 8      
    }

    if (any(idx <- ptype == ANYPARAM)) {
      neg <- idx & (lower <= 0)
      autostart[neg] <-  0.5 * (lower[neg] + upper[neg])
      pos <- idx & (lower > 0)
      autostart[pos] <-  sqrt(lower[pos]*upper[pos])
    }
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
  if (PrintLevel>1) cat("\nauto...")
  autopar <- autostart
  autopar[autopar == 0] <- 1
  varnugNA <- ## is.na(nugget) and is.na(var)
    globalsigma <- ## nugget==0
      sillbounded <- ## !is.na(sill)
        FALSE
  if (is.null(standard.style)) { ## useful for lsq and mle methods
    standard.style <- !missing(param)
  } else if (missing(param) && standard.style) {
    standard.style <- FALSE
    warning("standard.style must be FALSE for the given model specification.")
  }

  var.global <- var.idx
  if (standard.style) {
    if (PrintLevel>1) cat("\nstandard style settings...")
    variance <- param[VARIANCE]
    nugget <- param[NUGGET]
    covnr <- .C("GetModelNr", as.character(model), nr=integer(1),
                PACKAGE="RandomFields")$nr
    n.kappas <- .C("GetNrParameters", covnr, k=integer(1),
                   PACKAGE="RandomFields", DUP=FALSE)$k

    stopifnot(n.kappas >= 0)

    if ((covnr==.C("GetModelNr", as.character("nugget"), 
           nr=integer(1), PACKAGE="RandomFields")$nr) &&
        any(!is.finite(param[-VARIANCE]) | param[-VARIANCE] != 0.0))
        stop("if model=nugget effect then the nugget variance must be zero, and only mean and/or variance can be NA")   
    if (sillbounded <- !is.na(sill)) {
      ## only VARIANCE need to be optimised
      ## NUGGET = SILL - VARIANCE
      stopifnot(sill > 0)

      if (xor(is.na(variance), is.na(nugget)))
        stop("Sill fixed. Then variance and nugget should be given or unknown simultaneously.")
      if (!is.na(variance) && (variance + nugget!=sill)) {
        cat("var=", variance , "nug=", nugget, "sill=", sill,
            "sill-var-nug=", sill - variance - nugget,"\n")
        stop("sill !=  variance + nugget")
      }

      if (is.na(variance)) {
        autopar[c(nugget.idx, var.idx)] <-
          autostart[var.idx] <- sill/ length(c(nugget.idx, var.idx))
        upper[var.idx] <- sill
        delete.idx <- c(delete.idx, nugget.idx)
        lower[var.idx] <- 0
        transform <- function(x) {
          y <- numeric(length(x) + 1 + length(mixed.idx) )
          y[-c(nugget.idx, mixed.idx)] <- x
          y[nugget.idx] <- sill - y[var.idx]
          y[mixed.idx] <- 1
          y
        }
      }
    } else { ## not sill bounded
      if (is.na(variance)) {
        if (varnugNA <- is.na(nugget)) {
          ## of interest currently 
          ## both NUGGET and VARIANCE have to be optimised
          ## instead optimise the SILL (which can be done by
          ## a formula) and estimate the part of NUGGET
          ## (which will be stored in NUGGET) !
          ##  consequences:
          ## * estimtation space increased only by 1, not 2
          ## * real nugget: NUGGET * SILL
          ## *     variance: (1-NUGGET) * SILL
          autostart[nugget.idx] <- 1 / 2 ## lassen, da nicht standard Alg.
 
          upper[nugget.idx] <- 1
          delete.idx <- c(delete.idx, var.idx)
          transform <- function(x) {
            y <- numeric(length(x) + 1 + length(mixed.idx))
            y[-c(var.idx, mixed.idx)] <- x
            y[var.idx] <- 1.0 - y[nugget.idx]
            y[mixed.idx] <- 1
            y
          }
        } else { ## not sillbounded, is.na(variance), !is.na(nugget)
          if (globalsigma <- nugget==0.0) {
            ## here SILL==VARIANCE and therefore, variance
            ## can be estimated by a formula without increasing the dimension

            ## more than only the variance has to be estimated
            delete.idx <- c(delete.idx, var.idx)
            # if (length(lower) == 0) stop("trivial case not programmed yet")
            transform <- function(x) {
              y <- numeric(length(x) + 1 + length(mixed.idx))
              y[-c(var.idx, mixed.idx)] <- x
              y[var.idx] <-  1.0 ## 30.12.09 ex: vardata 
              y[mixed.idx] <- 1
              y
            }
          } else { ## not sillbounded, is.na(variance), nugget!=0
            lower[var.idx] <-
              pmax(lower[var.idx], 0, (vardata-nugget)/lowerbound.var.factor)
            if (lower[var.idx] < lowerbound.sill) {
              if (PrintLevel>1)
                cat("low.var=",lower[VARIANCE]," low.sill",lowerbound.sill,
                    "\ estimated variance from data=", vardata,
                    "nugget=", nugget, "\n")
              warning("param[NUGGET] might not be given correctly.")
              lower[var.idx] <- lowerbound.sill
            }
            autostart[var.idx] <-  autopar[var.idx] <- lower[var.idx]
          }
        }
      } else { ## not sillbounded, !is.na(variance)
        if (variance==0.0) {
          ## TODO: muss auf globalsigma und var.global gesetzt werden
          if (any(is.na(param[-c(VARIANCE, SCALE)])))
            stop("If variance=0, estimating scale or model parameters does not make sense")
          lower[var.idx] <- 0
        }
        if (is.na(nugget)) { ## and !is.na(param[VARIANCE])
          lower[nugget.idx] <-
            pmax(0, (vardata - variance) / lowerbound.var.factor)
          if (lower[nugget.idx] < lowerbound.sill) {
            if (PrintLevel>1)
              cat("\nlower nugget bound=", lower[nugget.idx],
                  " < lower sill bound=", lowerbound.sill,
                  " -- is the variance given correctly?\n",sep="")
            ## warning("Has param[VARIANCE] been given correctly?!")
            lower[nugget.idx] <- lowerbound.sill
            autostart[nugget.idx] <- autopar[nugget.idx] <- lower[nugget.idx]
          }
        }
      } ##  else { ## not sillbounded, !is.na(variance)
    } ## else { ## not sill bounded

    stopifnot( varnugNA + globalsigma + sillbounded <= 1)
  } else { ### ! standard style
   
    if (anyeffect) {
      stopifnot(length(idx.error) == 1) ## idx
     
      globalsigma <- (modelinfo$sub[[idx.error]]$name =="$" &&
                      is.na(modelinfo$sub[[idx.error]]$param$var))
    } else {
      globalsigma <- (modelinfo$name=="$" && is.na(modelinfo$param$var))
    }
        
    if (globalsigma) {
      ## here SILL==VARIANCE and therefore, variance
      ## can be estimated by a formula without increasing the dimension
      if (anyeffect) {        
        eff <- rep(FALSE, nrow(minmax))
        eff[(csNAs[idx.error]+1) : csNAs[idx.error]] <- TRUE
      } else eff <- TRUE
      var.global <- which(eff & (ptype == VARPARAM | ptype == NUGGETVAR))[1]
      stopifnot(length(var.global) == 1)
      
      delete.idx <- c(delete.idx, var.global, which(minmax[,4]==1))
      trafo <- transform      
       if (is.null(trafo))
         transform <- function(x) {
           y <- numeric(length(x) + 1 + length(mixed.idx))
           y[-c(var.global, mixed.idx)] <- x
           y[c(var.global, mixed.idx)] <- 1.0
           y
         } 
       else
         transform <- function (x) {
           y <- numeric(length(x) + 1 + length(mixed.idx))
           y[-c(var.global, mixed.idx)] <- x
           y[c(var.global, mixed.idx)] <- 1.0
           trafo(y)
         }

       
    } else {
      delete.idx <- c(delete.idx, which(minmax[,4]==1))
    }
  }

  
  
  globalsigma <- globalsigma || varnugNA

  if (!is.null(users.lower))  {
    stopifnot(length(users.lower)==length(lower))
    idx <- !is.na(users.lower)
    lower[idx] <- users.lower[idx]
  }

  if (!is.null(users.upper))  {
    stopifnot(length(users.upper)==length(upper))
    idx <- !is.na(users.upper)
    upper[idx] <- users.upper[idx]
  }


  delete.idx <- c(delete.idx, mixed.idx)

  if (length(delete.idx)>0) {    
    upper <- upper[-delete.idx]
    lower <- lower[-delete.idx]
    autostart <- autostart[-delete.idx]
    autopar <- autopar[-delete.idx]
     varnames <- varnames[-delete.idx]
    SCALE.IDX <- SCALE.IDX[-delete.idx]
    ptype <- ptype[-delete.idx]
  }

  nvarWithoutBC <- length(lower)
  if (lambdaest) {
    varnames <- c(varnames, "BoxCox")
    autostart <- c(autostart, 1)
    autopar <- c(autopar, 1)
    lower <- c(lower, BC.lambdaLB)
    upper <- c(upper, BC.lambdaUB)
  }
  n.variab <- length(lower)

  if (any(lower > upper)) {
    Print(cbind(lower=lower, upper=upper))#
    stop("some lower bounds are greater than the upper bounds")
  }



######################################################################
###                                 Covariates                     ###
######################################################################
  if (PrintLevel>1) cat("\nCovariates...")


##############  Sinv, etc              #################
## 

  ##  ToDo !!! : check correct dimensions of X, etc

  Sinv <- list() ## coordinates
  Xges  <- list()
  
 ## startCoVar <- endCoVar <-  matrix(nr=sets, ncol=length(allcomponents))
  logdetbase <- 0

  ## all beta of mixed compents are fixed for each set
  ## and drawn independently among sets
  ## error part is independent for each repetition in/for each set



  if (anyeffect) {
    nCoVar <- endCoVar <- startCoVar <-
      matrix(NA, nrow=sets, ncol=length(allcomponents))
    
    cs <- 0
    s <- numeric(length(allcomponents))
    for (i in 1:sets) {   
      for (k in allcomponents) {
        s[k] <-
          if (effect[k]==RemainingError) 0
          else {
            if (length(modelinfo$sub[[k]]$param$X) > 0)
              ncol(modelinfo$sub[[k]]$param$X[[i]])
            else nrow(idx.na[[i]][[1]])
          }
      }
      
      nCoVar[i, ] <- s  
      
      endCoVar[i, ] <- cs + cumsum(s)          
      startCoVar[i, ] <- endCoVar[i, ] - s + 1 
      cs <- cs + sum(s)
    }
    nTotalComp <- apply(nCoVar, 2, sum)
    nCoVarSets <- apply(nCoVar, 1, sum)
    nCoVarAll <- sum(nCoVarSets)
   
    Ny <-
      sapply(modelinfo$sub,
             function(x) {
               if(x$name=="mixed" && !is.null(x$param$X)) sapply(x$param$X, nrow)
               else rep(0, sets)
             }) / vdim
    if (is.vector(Ny)) dim(Ny) <- c(1, length(Ny))   
  } 

  
  param <- as.double(transform(autostart)) ## only the 1s in mixed.idx needed  

   .C("PutValuesAtNA", param, PACKAGE="RandomFields", DUP=FALSE)

  for (i in 1:sets) {
    Xges[[i]] <- Sinv[[i]] <- list()
    if (anyeffect) {
      for (m in 1:ncol(idx.na[[i]][[1]])) {
        Xges[[i]][[m]] <- matrix(nrow=sum(idx.na[[i]][[1]][,m]), ncol=nCoVarAll)
        for (k in mmcomponents) {
          

          if (Ny[i,k] == 0) {
            stopifnot(nCoVar[i, k] == nrow(idx.na[[i]][[1]]))
          } else {
            if (Ny[i, k] != nrow(idx.na[[i]][[1]]))
              stop("length of data does not match X matrix")
          }
          Xges[[i]][[m]][, startCoVar[i, k] : endCoVar[i, k] ] <-
            (if (Ny[i,k] == 0) diag(nCoVar[i, k])
            else modelinfo$sub[[k]]$param$X[[i]]) [idx.na[[i]][[1]][,m], ]
         
        }
      }
    } else {
       nCoVarAll <- 0
       nCoVar <- matrix(nrow=1, ncol=1, 0)
   }

    for (k in allcomponents) {
      if (effect[k] > FixedEffect) {
        if (effect[k] < SpaceEffect) {       
          cov <- double(nCoVar[i, k]^2)
          
          .C("CovMatrixMLE", Xdistances[[i]], as.integer(is.dist),
             as.integer(lc[i]),
             as.integer(i), as.integer(k), 
             cov, PACKAGE="RandomFields", DUP=FALSE)
          dim(cov) <- rep(nCoVar[i, k], 2)
          
          cov <- chol(cov)
          logdetbase <- logdetbase + 2 * sum(log(diag(cov)))     
          Sinv[[i]][[k]] <- chol2inv(cov, LINPACK=TRUE)    
        } else {
          Sinv[[i]][[k]] <- if (effect[k] == SpaceEffect) NULL else list()
        }
      }
    }
  } # i in 1:sets



  eff <- efftmp <- rep(0, RemainingError + 1)
  for (k in mmcomponents) {
    if (effect[k] <= SpaceEffect) {
      eff[effect[k] + 1] <- eff[effect[k] + 1] + 1
    }
  }

  betanames <- character(nCoVarAll)
  for (i in 1:sets) {    
    for (k in mmcomponents) {
      if (effect[k] <= SpaceEffect) {
        efftmp[effect[k] + 1] <- efftmp[effect[k] + 1] + 1
        bn <- paste(EffectName[effect[k] + 1],
                    if (eff[effect[k] + 1] > 1) efftmp[effect[k] + 1], sep="")
        if (nCoVar[i, k]>1) bn <- paste(bn, 1:nCoVar[i, k], sep=".")
        betanames[startCoVar[i, k] : endCoVar[i, k]] <-
          if (sets == 1) bn else c(betanames, paste(i, bn, sep=""))
      }
    }
  }  
  

######################################################################
###                     Estimation part itself                     ###
######################################################################

  ## check optim.control 
  ## parscale will give the magnitude of the parameters to be estimated
  ##     passed to optim/optimise so that the optimiser estimates
  ##     values around 1

  parscale <- pmax(abs(lower), abs(upper)) / 10
  idx <- lower * upper > 0
  parscale[idx] <- sqrt(lower[idx] * upper[idx])
  stopifnot(all(is.finite(parscale)))
  if (length(optim.control)>0) {
    stopifnot(is.list(optim.control))
    forbidden.param <- c("parscale", "fnscale")
    ## fnscale=-1 turns the problem into a maximisation problem, see below
    if (any(!is.na(pmatch(names(optim.control), forbidden.param))))
      stop(paste(forbidden.param, collapse=" and "),
                 " may not be given in the optim.contol list")
  } else optim.control <- list()


###################  preparation  ################
  if (PrintLevel>1) cat("\npreparing fitting...")
  ## methods
  formals <- formals()
  allprimmeth <- c("autostart", "users.guess")
  lsq.orig.methods <- eval(formals$lsq.methods)
  nlsqinternal <- 3 ## cross checked after definition of weights below
  alllsqmeth <- c(lsq.orig.methods[-length(lsq.orig.methods)],
                  paste("internal", 1:nlsqinternal, sep=""))
  allmlemeth <- eval(formals$mle.methods)
  allcrossmeth <- eval(formals$cross.methods)
  allmethods <- c(allprimmeth, alllsqmeth, allmlemeth, allcrossmeth)

  ## how preceding methods have been considered ?
  ## note cm is used again at the very end when error checking
  cm <- cumsum(c(0, length(allprimmeth), length(alllsqmeth),
                     length(allmlemeth), length(allcrossmeth)))
  cm <- cbind(cm[-length(cm)] + 1, cm[-1])
  cm <- apply(cm, 1, function(x) x[1] : x[2])
  names(cm) <- c("prim", "lsq", "mle", "cross")

  methodprevto <-
    if (only.users) list(lsq="users.guess",mle="users.guess",cross="users.guess")
    else list(lsq=c(cm$prim),
              mle=c(cm$prim, cm$lsq),
              cross=c(cm$prim, cm$lsq, cm$cross)
              )

  
  ## index (start, end) to the various categories of
  ## information to be stored
  tblidx <- cumsum(c(0,
                     n.variab, # variables used in algorithm
                     length(lower), # their lower bounds
                     length(upper), # ... and upper bounds
                     ncovparam,  # param values to be estimated
                     rep(1, length(allmethods) - length(allprimmeth)),#method
                     ##                                                 score
                     nCoVarAll # coeff to estimated for covariates, i.e.
                     ##           mixed effects and trend parameters
                    ))


  tblidx <- rbind(tblidx[-length(tblidx)] + 1, tblidx[-1])
  idx <- tblidx[1, ] > tblidx[2, ]
  tblidx[, idx] <- 0
 
  dimnames(tblidx) <- list(c("start", "end"),
                           c("variab", "lower", "upper", "param", 
                             allmethods[-1:-length(allprimmeth)],
                             "covariab"
                             ##,  "lowbeta", "upbeta", only used for
                             ## cross-validation
                             ))
  maxtblidx <- max(tblidx)
  tblidx <- data.frame(tblidx)

  ## table of all information; col:various method; row:information to method

   tablenames <-
    c(
      if (n.variab > 0) {
        paste("v", if (is.null(ptype)) 1:n.variab else varnames, sep=":")
      },        
      if (n.variab > 0) { #paste("l", varnames, sep=":"),
        paste("lb", if (is.null(ptype)) 1:n.variab else varnames, sep=":")
      },
    #  if (nCoVarAll > 0) paste("lower", betanames, sep=":"),
      if (n.variab>0) { #paste("u", varnames, sep=":"),
        paste("ub", if (is.null(ptype)) 1:n.variab  else varnames, sep=":")
      }, #,
      if (nrow(minmax) > 0) {
        if (is.null(ptype) && FALSE) paste("p", 1:nrow(minmax), sep=":")
        else minmax.names
      },
      allmethods[-1:-length(allprimmeth)],
      betanames
      ## do not try to join the next two lines, since both
      ## varnames and betanames may contain nonsense if
      ## n.variab==0 and nCoVarAll==0, respectively

      ## if (nCoVarAll > 0) paste("upper", betanames, sep=":")
      )

#  print(list(tablenames, allmethods))

  
  param.table <- data.frame(matrix(NA, nrow=maxtblidx, ncol=length(allmethods),
                                   dimnames=list(tablenames, allmethods)))

  
#############################################################
## end preparation; remaining part is estimation  ###########
#############################################################

  MLELB <- LSQLB <- lower
  MLEUB <- LSQUB <- upper
  
##################################################
###############    PRIMITIVE METHODS   ###########
##################################################

  
  ##****************    autostart    *****************
  if (PrintLevel>1) cat("\nautostart...")
  M <- "autostart"
  default.param <- param.table[[M]][IDX("variab")] <- autostart
  param.table[[M]][IDX("param")] <- transform(autostart) 
  ## ****************    user's guess    *****************
  if (!is.null(users.guess)) {
    M <- "users.guess"
    if (length(users.guess) != nrow(minmax))
        stop("users.guess must contain all NA/NaN parameters")
    if (length(delete.idx) > 0) users.guess <- users.guess[-delete.idx]  
    if (any(idx <- users.guess < lower | users.guess > upper)) {
      m <- cbind(lower, users.guess, upper, idx)
      dimnames(m) <- list(rep("", length(lower)),
                          c("lower", "user", "upper", "outside bounds"))
      Print(m)
      stop("not all users.guesses within bounds\n change values of `lower' and `upper' or \nthose of the `lowerbound*.factor's and `upperbound*.factor's")
    }
    param.table[[M]][IDX("variab")] <- users.guess
    param.table[[M]][IDX("param")] <- transform(users.guess)
  }


 

  ##****************    autostart    *****************
  param.table[[M]][IDX("covariab")] <- get.covariates(autostart)
  ## ****************    user's guess    *****************

 # Print(autostart, users.guess)
  if (!is.null(users.guess))
    param.table[[M]][IDX("covariab")] <- get.covariates(users.guess)
 # Print("OK")


##################################################
###################  LSQ  ########################
  ## see above for the trafo definitions
  ##
  ## zuerst regression fit fuer variogram,
  ## dann schaetzung der Parameter, dann berechnung der covariates
  ## ToDo: kann verbessert werden durch einschluss der Trendschaetzung
  ##************   Empirical Variogram    ***********
  lsqMethods <- NULL
  ev <- list()

   
  if (stationary && !is.null(lsq.methods)) {
    if (PrintLevel>1) cat("\nempirical variogram...\n")
    sd <- vector("list", vdim)
    emp.vario <- vector("list", vdim)
    index.bv <- NULL
 
    for (j in 1:vdim) {
      for (i in 1:sets) {
        m <- 1
        if ((nrow(idx.na[[i]][[2]]) > 1) && PrintLevel>0)
          cat("only first column(s) used for empirical variogram\n")

        W <- werte[[i]][idx.na[[i]][[1]][, m], j ,idx.na[[i]][[2]][m, ]]
        if (nCoVarAll > 0) {
          regr <- lsfit(Xges[[i]][[m]], W, intercept=FALSE)
          W <- regr$residuals
        }
        
        if (length(nphi)==1)
          nphi <- c(0, nphi) # starting angle; lines per half circle
        if (length(ntheta)==1) ntheta <- c(0, ntheta) # see above
        if (length(ntime)==1) ntime <- c(ntime, 1) 
        ntime <- ntime * T[3] ## time endpoint; step
        
        bin <- if (length(bins)>1) bins else c(-1, seq(0, distance.factor *
                                                       maxdistances, len=bins+1))
        if (xgiven) {
          ev[[i]] <-
            EmpiricalVariogram(t(coord[[i]]), T=T, data=W, grid=FALSE, bin=bin,
                               phi=if ((spacedim>=2) && !isotropy) nphi,
                               theta=if ((spacedim>=3)&& !isotropy) ntheta,
                               deltaT=if (!is.null(T)) ntime
                               ) # 7.1.11 coord -> t(coord) !
         } else {
          n.bin <- vario2 <- vario <- rep(0, length(bin))
          stopifnot(length(lc) == 1)
          k <- 1
          for (g in 1:(lc[1]-1)) {
            # cat(g, "")
            for (h in (g+1):lc) {
              idx <- sum(Distances[[i]][k] > bin)
              n.bin[idx] <- n.bin[idx] + 2
              vario[idx] <- vario[idx] + mean((W[g] - W[h])^2)
              vario2[idx] <- vario2[idx] + mean((W[g] - W[h])^4)
              k <- k + 1
            }
          }
          vario <- vario / n.bin ## Achtung! n.bin ist bereits gedoppelt
          sdvario  <- sqrt(vario2 / n.bin / 2.0 - vario^2)
          sdvario[1] <- vario[1] <- 0
          n.bin[1] <- lc[1]
          centers <- 0.5 * (bin[-1] + bin[-length(bin)])
          centers[1] <- 0
          ev[[i]] <- list(centers=centers,
                          emp.vario=vario[-length(n.bin)],
                          sd=sdvario[-length(n.bin)],
                          n.bin=n.bin[-length(n.bin)])
        }
      }
        
      n.bin.raw <- sapply(ev, function(x) x$n.bin)
      n.bin <- rowSums(n.bin.raw)     
      
      sd[[j]] <- sqrt((sapply(ev, function(x) x$sd^2 * x$n.bin) %*%
                       rep(1, length(ev))) / n.bin)
      emp.vario[[j]] <- sapply(ev, function(x) x$emp.vario * x$n.bin) %*%
        rep(1, length(ev)) / n.bin
      index.bv <- cbind(index.bv, !is.na(emp.vario[[j]])) ## exclude bins without
      ##                                                    entry
    } # j in 1:vdim

    index.bv <- apply(index.bv, 1, all)
 
    if (sum(index.bv) <= 1)
      stop("not more than 1 value in empirical variogram that is not NA; check values of bins and distance.factor")

    

    sd <- lapply(sd, function(x) x[!is.nan(x)])
    sd <- sqrt(sum(sapply(sd, function(x) x^2)))
    ## bei vdim>1 werden nur die diag-elemente von binned.vario gefuellt:
    bvtext <- paste("emp.vario[[", 1:vdim, "]][index.bv]", collapse=", matrix(0, nrow=sum(index.bv), ncol=vdim) ,", sep="")
    binned.variogram <- eval(parse(text=paste("cbind(", bvtext,")", sep="")))
    binned.variogram <- matrix(t(binned.variogram))[ , , drop=TRUE]
    bin.centers <- as.matrix(ev[[1]]$centers)

    if (!is.null(ev[[1]]$phi)) {
      if (spacedim<2) stop("x dimension is less than two, but phi is given") 
      bin.centers <- cbind(as.vector(outer(bin.centers, cos(ev[[1]]$phi))),
                           as.vector(outer(bin.centers, sin(ev[[1]]$phi))))
    }
    if (!is.null(ev[[1]]$theta)) {
      if (spacedim<3)
        stop("x dimension is less than three, but theta is given") 
      if (ncol(bin.centers)==1) bin.centers <- cbind(bin.centers, 0)
      bin.centers <- cbind(as.vector(outer(bin.centers[, 1],
                                           cos(ev[[1]]$theta))),
                           as.vector(outer(bin.centers[, 2],
                                           cos(ev[[1]]$theta))),
                           rep(sin(ev[[1]]$theta), each=nrow(bin.centers)))
    } else {
      
    #  warning("must be ncol()")      
      if (ncol(bin.centers) < spacedim) { # dimension of bincenter vector
        ##                       smaller than dimension of location space      
        bin.centers <- 
          cbind(bin.centers, matrix(0, nrow=nrow(bin.centers),
                                    ncol=spacedim - ncol(bin.centers)
                                    ))
      }
    }
    if (!is.null(ev[[1]]$T)) {
      bin.centers <-
        cbind(matrix(rep(t(bin.centers), length(ev[[1]]$T)), byrow=TRUE,
                     ncol = ncol(bin.centers)),
              rep(ev[[1]]$T, each=nrow(bin.centers)))
    }

    if (isotropy) {
      bin.centers <- as.matrix(sqrt(apply(bin.centers^2, 1, sum)))
    }  
    bin.centers <- as.double(t(bin.centers[index.bv, ])) #
    ##  es muessen beim direkten C-aufruf die componenten der Punkte
    ##  hintereinander kommen (siehe auch variable coord, Xdistance). Deshalb t()
    
    evsd <- as.double(sd)
    evsd[is.na(evsd)] <- 0
    evsd[evsd==0] <- 10 * sum(evsd, na.rm=TRUE) ## == "infinity"

    bins             <- length(n.bin)
    binned.n         <- as.integer(n.bin)

    weights <- cbind(NA,                      # self
                     rep(1, bins),            # plain 
                     sqrt(binned.n),          # sqrt(#)
                     1 / evsd,                # sd^-1
                     sqrt(bins:1 * as.double(binned.n)), # internal1 # kann sonst
                     ##                   fehler verursachen, da integer overflow
                     bins:1,                  # internal2
                     sqrt(bins:1)             # internal3
                     )[index.bv, ]


 
    
    stopifnot(ncol(weights)==length(alllsqmeth))
    dimnames(weights) <- list(NULL, alllsqmeth)
    weights <- data.frame(weights)
    bins <- as.integer(sum(index.bv))
    rm(W)
    
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
    if (PrintLevel>1) cat("\nestimation part...")

    LSMIN <- Inf
    lsqMethods <- lsq.orig.methods[pmatch(lsq.methods, lsq.orig.methods)]
    if (!is.null(lsqMethods) &&
        any(is.na(lsqMethods))) stop("not all lsq.methods could be matched")
    if ("internal" %in% lsqMethods)
      lsqMethods <- c(lsqMethods, paste("internal", 1:nlsqinternal, sep=""))

    
     if (lambdaest) {
      ## koennte man prinzipiell besser machen...
      indices <- is.na(match(tablenames,c("v:BoxCox", "lb:BoxCox", "ub:BoxCox")))
      for (M in c(alllsqmeth)) {
        param.table[[M]][indices] <- 1
      #  param.table[[M]][!indices] <- 1
      }
    }

 
    for (M in c(alllsqmeth)) {
      if (!(M %in% lsqMethods)) next;
      if (PrintLevel>2) cat("\n", M) else cat(pch)

      param.table[[M]][IDX("variab")] <- default.param
      
      LSQsettings(M)
      LSMIN <- Inf ## must be before next "if (n.variab==0)"
      LSPARAM <- LSVARIAB <- NA 
      ##      if (n.variab == 0) {
 #       warning("trivial case may cause problems")
 #     } else {
      param.table[[M]][IDX("lower")] <- LSQLB
      param.table[[M]][IDX("upper")] <- LSQUB
      options(show.error.messages = show.error.message) ##

      
      if (n.variab == 0) {
        LStarget(param.table[IDX("variab"), methodprevto$lsq[1]])
      } else {
        if (n.variab == 1) {
          ##variab <-
          try(optimize(LStarget, lower = LSQLB, upper = LSQUB)$minimum,
              silent=silent)
        } else {
          min <- Inf
          for (i in methodprevto$lsq) { ## ! -- the parts that change if
            ##                               this part is copied for other methods
            if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {
              value <- LStarget(variab) ## !
              if (is.finite(value)) {
                param.table[tblidx[[M]][1], i] <- value
                if (value < min) {
                  min.variab <- variab
                  min <- value
                } else {
                  param.table[tblidx[[M]][1], i] <- NaN
                  next
                }
              }
            }
          }
          stopifnot(min.variab==LSVARIAB, min==LSMIN) ## check 

          lsq.optim.control <-
            c(optim.control, list(parscale=parscale, fnscale=min))          

          try(optim(LSVARIAB, LStarget, method ="L-BFGS-B", lower = LSQLB,
                    upper = LSQUB, control= lsq.optim.control)$par,
              silent=silent)
         } # n.variab > 1
      } # n.variab > 0
      options(show.error.messages = show.error.message)  
      ## side effect: minimum so far is in LSMIN and LSPARAM
      ## even if the algorithm finally fails
      if (is.finite(LSMIN)) {
        param.table[[M]][tblidx[[M]][1]] <- LSMIN
        param.table[[M]][IDX("variab")] <- LSVARIAB
        param.table[[M]][IDX("param")] <- LSPARAM
      } else {
        param.table[[M]] <- if (n.variab==0) NA else NaN
      }  
      param.table[[M]][IDX("covariab")] <- get.covariates(LSVARIAB)
    } # for M
  } # stationary


##################################################
### optional parameter grid for MLE and CROSS  ###

  
  idx <- IDX("variab")
  gridmax <- as.matrix(param.table[idx, cm$lsq])
  if (!any(is.finite(gridmax))) gridmax <- as.matrix(param.table[idx, ])
   
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
  
   MLEtarget <- NULL
   mleMethods <- (if (is.null(mle.methods)) NULL else
              allmlemeth[pmatch(mle.methods, allmlemeth)])
  if ("reml" %in% mleMethods && nCoVarAll == 0)
    mleMethods <- c("ml", "reml")

  ## lowerbound.scale.LS.factor <  lowerbound.scale.factor, usually
  ## LS optimisation should not run to a boundary (what often happens
  ## for the scale) since a boundary value is usually a bad initial
  ## value for MLE (heuristic statement). Therefore a small
  ## lowerbound.scale.LS.factor is used for LS optimisation.
  ## For MLE estimation we should include the true value of the scale;
  ## so the bounds must be larger. Here lower[SCALE] is corrected
  ## to be suitable for MLE estimation
  if (any(SCALE.IDX)) {
    MLELB[SCALE.IDX] <- MLELB[SCALE.IDX] *
      lowerbound.scale.LS.factor / lowerbound.scale.factor
  } else if (!is.null(ptype) &&
           any(idx <- ptype == DIAGPARAM | ptype == ANISOPARAM))  
    MLEUB[idx] <- MLEUB[idx] *
      lowerbound.scale.factor / lowerbound.scale.LS.factor

  ## fnscale <- -1 : maximisation
  for (M in c(allmlemeth)) {
    assign("LINEARPARTS2", rep(1, length(mixed.idx)), envir=ENVIR)
    if (!(M %in% mleMethods)) next;
    if (PrintLevel>2) cat("\n", M) else cat(pch)
    param.table[[M]][IDX("variab")] <- default.param
    if (M=="reml") {
      if (nCoVarAll == 0) { ## same as MLE
        param.table[[M]] <- param.table[["ml"]]
        param.table[[M]][tblidx[[M]][1]] <- param.table[[M]][tblidx[["ml"]][1]]
        next
      } else {
        stop("reml currently not programmed.")
      }
    }

    MLEsettings(M)
    MLEMAX <- -Inf ## must be before next "if (nMLEINDEX==0)"
    MLEVARIAB <- Inf
    MLEPARAM <- NA
    if (length(MLELB) == 0) {
       MLEtarget(NULL)
   } else {
      param.table[[M]][IDX("lower")] <- MLELB
      param.table[[M]][IDX("upper")] <- MLEUB
      options(show.error.messages = show.error.message) ##
      max <- -Inf
      for (i in methodprevto$mle) { ## ! -- the parts that change if
        ##                             this part is copied for other methods
        ## should mle be included when M=reml?
        ## same for lsq methods as well: should previous result be included?

        if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {

         
          value <- MLEtarget(variab) ## !
          if (is.finite(value)) {
            param.table[tblidx[[M]][1], i] <- value
            if (value > max) {
              max.variab <- variab
              max <- value
            }
          } else {
            param.table[tblidx[[M]][1], i] <- NaN
            next
          }
        }
      }
 
      mle.optim.control <-
          c(optim.control, list(parscale=parscale, fnscale=-max(abs(max), 0.1)))

      MLEINF <- FALSE
 
      try(
      optim(MLEVARIAB,
            MLEtarget, method="L-BFGS-B", lower = MLELB,
            upper=MLEUB, control=mle.optim.control)
          , silent=silent)

      if (MLEINF) {
        if (PrintLevel>2) Print("MLEINF", MLEVARIAB, MLEMAX) else cat("#")
        try(optim(MLEVARIAB,
                  MLEtarget, method="L-BFGS-B", lower = MLELB,
                  upper=MLEUB, control=mle.optim.control), silent=silent)
        if (PrintLevel>2) Print("MLEINF new", MLEVARIAB, MLEMAX)
      }
      
      
      options(show.error.messages = TRUE) ##
      mindistance <- pmax(minbounddistance, minboundreldist * abs(MLEVARIAB))
      
      onborderline <- 
        (abs(MLEVARIAB - MLELB) <
             pmax(mindistance,  ## absolute difference
                  minboundreldist * abs(MLELB) ## relative difference
                  )) |
         (abs(MLEVARIAB - MLEUB) <
          pmax(mindistance, minboundreldist * abs(MLEUB)))
    }

    if (PrintLevel>2) Print("mle first round", MLEVARIAB, MLEPARAM, MLEMAX)
    
    if (!is.finite(MLEMAX)) {
      if (PrintLevel>0) cat(M, "MLEtarget I failed.\n")
      param.table[[M]] <- MLEPARAM <- NaN
      variab <- MLELB ## to call for onborderline
      ml.residuals <- NA
    } else {
      param.table[[M]][tblidx[[M]][1]] <- MLEMAX
      param.table[[M]][IDX("variab")] <- MLEVARIAB
      param.table[[M]][IDX("param")] <- MLEPARAM
      param.table[[M]][IDX("covariab")] <- get.covariates(MLEVARIAB)     
      ml.residuals <- ML.RESIDUALS
       
      if (length(MLELB) > 0 && any(onborderline) && refine.onborder &&
          !only.users && n.variab > 1) {
        ## if the MLE result is close to the border, it usually means that
        ## the algorithm has failed, especially because of a bad starting
        ## value (least squares do not always give a good starting point,helas)
        ## so the brutal method:
        ## calculate the MLE values on a grid and start the optimization with
        ## the best grid point. Again, there is the believe that the
        ## least square give at least a hint what a good grid is
        MLEgridmin <- gridmin
        MLEgridmax <- gridmax
        
        if (any(is.na(MLEgridmin)) || any(is.na(MLEgridmax))) {
          if (PrintLevel > 1) {
            Print(cbind(MLELB, variab, MLEUB, onborderline),
                  MLEgridmin, MLEgridmax)
          }
          warning(paste(M, "converged to a boundary value -- ",
                        "better performance might be obtained",
                        "when allowing for more lsq.methods"))
        } else {
          if (PrintLevel>5) show(1, M, MLEMAX, MLEVARIAB) else cat(detailpch)
          MLEgridlength <-
            max(3, round(approximate.functioncalls^(1/n.variab)))
          ## grid is given by the extremes of the LS results
          ## so, therefore we should examine above at least 4 different sets
          ## of weights wichtig: gridmin/max basiert auf den reduzierten Variablen
          step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-2) # grid starts
          							# bit outside
          MLEgridmin <- pmax(MLEgridmin - step/2, MLELB)   # the extremes of LS
          MLEgridmax <- pmin(MLEgridmax + step/2, MLEUB)
          step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-1)

          startingvalues <- vector("list", len=length(step))
          for (i in 1:length(step)) {
            startingvalues[[i]] <- MLEgridmin[i] + step[i] * 0:(MLEgridlength-1)
          }

          startingvalues <- do.call("expand.grid", startingvalues)
          limit <- 10 * approximate.functioncalls
          if ((rn <- nrow(startingvalues)) > limit) {
            if (PrintLevel>1)
              cat("using only a random subset of the", rn, "grid points")
            rand <- runif(rn)
            startingvalues <-
              startingvalues[rand < quantile(rand, limit / rn), ]
            gc()
          }

          MLEMAX <- -Inf

          apply(startingvalues, 1, function(x) try(MLEtarget(x), silent=silent))

          if (PrintLevel>2) Print("mle grid search", MLEVARIAB, MLEPARAM, MLEMAX)
 
          
          ## side effect:Maximum is in MLEMAX!
          ##                             and optimal parameter is in MLEVARIAB
          if (PrintLevel>5) show(2, M, MLEMAX, MLEVARIAB)

          cat(detailpch)
          options(show.error.messages = show.error.message) ##
          try(optim(MLEVARIAB, MLEtarget, method ="L-BFGS-B",lower = MLELB,
                    upper = MLEUB, control=mle.optim.control)$par
              , silent=silent)
          options(show.error.messages = TRUE) ##
          if (!is.finite(MLEMAX) &&(PrintLevel>0))
            cat("MLtarget II failed.\n")
          ## do not check anymore whether there had been convergence or not.
          ## just take the best of the two strategies (initial value given by
          ## LS, initial value given by a grid), and be happy.
          if (PrintLevel>5) show(3, M, MLEMAX, MLEVARIAB) else 
          if (PrintLevel>2) Print("mle second round", MLEVARIAB, MLEPARAM, MLEMAX)
            
          if (is.finite(MLEMAX) && MLEMAX > param.table[[M]][tblidx[[M]][1]]) {
            param.table[[M]][tblidx[[M]][1]] <- MLEMAX
            param.table[[M]][IDX("variab")] <- MLEVARIAB
            param.table[[M]][IDX("param")] <- MLEPARAM
            param.table[[M]][IDX("covariab")] <- get.covariates(MLEVARIAB)
            ml.residuals <- ML.RESIDUALS
          }
        } # (is.na(MLEgridmin[1]))
      } # onborderline
    } # length(MLELB) != 0
  } ## M


########  estimation by cross validation  ########
  nCROSStotINDEX <- nCROSSINDEX <- n.variab
###############################################################
##  besser: lsqlower oder mlelower, dann andere Trafo, etc!  ## 
###############################################################
  if (FALSE) {
    werte <- as.matrix(werte)
    crosstrafo <- NULL;

    if (!missing(param)) {
      transform <- function(x) x
    }

    crosslower <- lower
    crossupper <- upper
    CROSSinvTrafo <- function(param) param

    CROSSINDEX <-  NA ## has to be deleted / rewritten

    if (standard.style) stop("standard style not allowed in cross validation")

    crossMethods <- (if (is.null(cross.methods)) NULL else
                     allcrossmeth[pmatch(cross.methods, allcrossmeth)])
    CROSS.lcrepet <- lc * repet
    cross.optim.control <-
      c(optim.control, list(parscale=list(parscale[CROSSINDEX]), fnscale=1))
    crossLBcovariates <-  crossUBcovariates <- NULL
    if (nCoVarAll > 0) {
      stopifnot(is.numeric(trend) && length(trend)==1)
      crossLBcovariates <- min(werte)
      crossUBcovariates <- max(werte)
      ## die neuen Grenzen sind 1. katastrophal schlecht; 2. werden
      ## sie nicht in die Tabelle uebernommen!!; 3. gibt es chaos mit
      ## den schnellen loesungen
      CROSSLB <- c(CROSSLB, crossLBcovariates)
      CROSSUB <- c(CROSSUB, crossUBcovariates)
      nCROSStotINDEX <- nCROSSINDEX + nCoVarAll # ??
      cross.optim.control$parscale <- ### auch nicht gut !!!
        list(cross.optim.control$parscale, rep(1, nCoVarAll))
    }
  }
  allcrossmeth <- NULL
  for (M in allcrossmeth) {
    if (!(M %in% crossMethods)) next;
    if (PrintLevel>2) cat("\n", M) else cat(pch)
     stopifnot(is.numeric(trend) && length(trend)==1)
    ## universal kriging not programmed yet
    crosssettings(M)    
    CROSSMIN <- Inf
    param.table[[M]] <- default.param ## in case the covariates are not estimated

    ## optimisation part
    RFparameters(Print=0)
    if (nCROSStotINDEX == 0) {
      crosstarget(numeric(0))
    } else {
      ixdCROSSINDEX <-  NaN ###
      if (length(ixdCROSSINDEX) > 0) {
        param.table[[M]][IDX("lower")][ixdCROSSINDEX] <- CROSSLB[1:nCROSSINDEX]
        param.table[[M]][IDX("upper")][ixdCROSSINDEX] <- CROSSUB[1:nCROSSINDEX]
      }
      if (nCoVarAll > 0) {
       ## param.table[[M]][IDX("lowbeta")] <- CROSSLB[nCROSSINDEX + 1:nCoVarAll]
       ## param.table[[M]][IDX("upbeta")] <- CROSSUB[nCROSSINDEX + 1:nCoVarAll]
      }
      options(show.error.messages = show.error.message) ##  
      if (nCROSStotINDEX==1) {
        variab <-
          try(optimize(crosstarget, lower = CROSSLB, upper = CROSSUB)$minimum,
              silent=silent)
      } else {
        min <- Inf
        for (i in methodprevto$cross) { ## ! -- the parts that change if
          ##                               this part is copied for other methods
          if (!is.na(param.table[1, i])) {
            variab <- CROSSinvTrafo(param.table[IDX("variab"), i])[CROSSINDEX]
            if (nCoVarAll > 0)
              variab <- c(variab, param.table[IDX("covariab"), i])
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
              silent=silent)
      } # nCROSStotINDEX > 1

      ## check onborderline
      variab <- CROSSMODEL$param[CROSSINDEX]
      if (nCoVarAll > 0) variab <- c(variab, CROSSMODEL$mean)

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
      param.table[[M]][tblidx[[M]][1]] <- CROSSMIN
      param.table[[M]][IDX("variab")] <- CROSSMODEL$param
      if (nCoVarAll > 0) {
        stopifnot(is.numeric(trend) && length(trend)==1)
        param.table[[M]][IDX("covariab")] <- CROSSMODEL$mean
      }
    } else {
      if (PrintLevel>0) cat(M, "target I failed.\n")
      param.table[[M]] <- NaN
    }
    if (nCROSStotINDEX > 0 && onborderline && refine.onborder) {
      CROSSgridmin <- c(CROSSinvTrafo(gridmin)[CROSSINDEX], crossLBcovariates)
      CROSSgridmax <- c(CROSSinvTrafo(gridmax)[CROSSINDEX], crossUBcovariates)
      if (any(is.na(CROSSgridmin)) || any(is.na(CROSSgridmax))) {
        warnung <- paste(M, "converged to a boundary value -- better performance might be",
                         "obtained when allowing for more lsq.methods")
      #  if (options(warn==2))
          Print(warnung)
        # else warning(warnung)
      } else {
        if (PrintLevel>5) show(1, M, CROSSMIN, CROSSMODEL$param)
        else cat(detailpch)
        CROSSgridlength <-
          max(3, round(approximate.functioncalls ^
                       (1/(nCROSStotINDEX + nCoVarAll))))
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
          if (PrintLevel>1)
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
          if (nCoVarAll > 0) {
            variab <- c(variab, param.table[idxCovar[1]:idxCovar[2], i])
          }
          options(show.error.messages = show.error.message) ##
          variab <-
            try(optim(variab, crosstarget, method ="L-BFGS-B",lower = CROSSLB,
                      upper = CROSSUB, control=cross.optim.control)$par,
                silent=silent)
          options(show.error.messages = TRUE) ##
          if (!is.numeric(variab) && (PrintLevel>0))
            cat("cross target II failed.\n")
        }
        if (PrintLevel>2) show(3, M, CROSSMIN, CROSSMODEL$par) 

        if (is.finite(CROSSMIN) && CROSSMIN < param.table[[M]][tblidx[[M]][1]]) {
          param.table[[M]][tblidx[[M]][1]] <- CROSSMIN
          param.table[[M]][IDX("variab")] <- CROSSMODEL$param
          if (nCoVarAll > 0) {
            stopifnot(is.numeric(trend) && length(trend)==1)
            param.table[[M]][IDX("covariab")] <- CROSSMODEL$mean
          }
        }
      } # (is.na(CROSSgridmin[1]))
    } # onborderline
  ##RFparameters(Print=save.RFparameters$PrintLevel)
    
    if (is.finite(param.table[[M]][tblidx[[M]][1]])) {
      param.table[[M]][IDX("covariab")] <-
        get.covariates(param.table[[M]][IDX("variab")])
    }
  } ## for M
  

######################################################################
###     calculate all target values for all optimal parameters     ###
######################################################################
  for (i in 1:length(allmethods)) if (!is.na(param.table[1, i])) {
    idx <- IDX("variab")
    idxCovar <- IDX("covariab")

    if (!lambdaest) {
      for (M in alllsqmeth) {
        cur <- param.table[tblidx[[M]][1], i]
        if (is.na(cur) && !is.nan(cur) && M %in% lsqMethods) {
          LSQsettings(M)
          param.table[tblidx[[M]][1], i] <- LStarget(param.table[idx, i])
        }
      }
    }

    for (M in allmlemeth) {
      cur <- param.table[tblidx[[M]][1], i]
      if (is.na(cur) && !is.nan(cur) && M %in% mleMethods) {
        MLEsettings(M)
        param.table[tblidx[[M]][1], i] <-
          MLEtarget(param.table[idx, i])
      }
    }
    
    for (M in allcrossmeth) {
      cur <- param.table[tblidx[[M]][1], i]
      if (is.na(cur) && !is.nan(cur) && M %in% crossMethods) {
        crosssettings(M)
        variab <- param.table[idx, i]
        if (nCoVarAll > 0) {
          variab <- c(variab, param.table[idxCovar, i])
        }
        param.table[tblidx[[M]][1], i] <- crosstarget(variab)
      }
    }
  }
  if (pch!="") cat("\n")



######################################################################
###                     error checks                               ###
######################################################################
  ## if rather close to nugget and nugget==fixed, do exception handling.
  ## This part is active only if
  ## scale.max.relative.factor < lowerbound.scale.LS.factor
  ## By default it is not, i.e., the following if-condition
  ## will "always" be FALSE.
  if (standard.style && !is.na(nugget)) {
    idx <- IDX("variab")
    alllsqscales <- param.table[idx, cm$lsq][SCALE, ]
    if (any(alllsqscales < mindistances/scale.max.relative.factor, na.rm=TRUE))
      warning(paste(sep="",
                    "Chosen model seems to be inappropriate!\n Probably a ",
                    if (nugget!=0.0) "larger ",
                    "nugget effect should be considered")
              )
  }


######################################################################
###                   format and return values                     ###
######################################################################


  ## if the covariance functions use natural scaling, just
  ## correct the final output by GNS$natscale
  ## (GNS$natscale==1 if no rescaling was used)
  ##
  ## currently natural scaling only for standard.style...


  if (use.naturalscaling && any(SCALE.IDX)) {
    idx <- IDX("param")
    for (i in 1:length(allmethods)) if (!is.na(param.table[1, i])) {
      param <- as.double(param.table[idx, i] + 0.0)
      .C("PutValuesAtNA", param, PACKAGE="RandomFields", DUP=FALSE)
      .C("expliciteDollarMLE", param, PACKAGE="RandomFields", DUP=FALSE)
      param.table[idx, i] <- param
    }
  }

 
  idx <- IDX("param")
  idxCovar <- IDX("covariab")
  idx.meth <- rep(FALSE, length(allmethods))
  res <- values.res <- list()
  
  for (i in 1:length(allmethods)) {
    M <- allmethods[i]
    if (idx.meth[i] <- !is.na(param.table[1, i]) || is.nan(param.table[1,i])
        || length(idx)==1 && idx==0) {
      .C("PutValuesAtNA", as.double(param.table[idx, i]),
         PACKAGE="RandomFields", DUP=FALSE)
      # modus: 1 : Modell wie gespeichert
      #        0 : Modell unter Annahme PracticalRange>0
      #        2 : natscale soweit wie moeglich zusammengezogen
      modus <- (old.param$Practical == 0) + (use.naturalscaling > 0)
      modelres <- GetModel(modelreg=ModelNr, modus=modus)
      if (modelres[[1]] == "|" && length(modelres[[1]]) <= 3) {#1 model nur
        modelres <- modelres[[3]]
      }
      
      if (trend.input) {
        mm <- which(sapply(modelres, function(x) {x[[1]]=="mixed"}))
        
        if (length(modelres) > 3) {
           sub <- modelres[-mm]
         } else {
           sub <- modelres[c(-1, -mm)][[1]]
         }   
         res[[M]] <-
          list(model=sub, trend=if (length(idxCovar)>0) param.table[idxCovar, i],
               residuals = if (M == "ml") ml.residuals,
               ml.value =  param.table[[M]][tblidx[["ml"]][1]]
               )
      } else {
        res[[M]] <- list(model=modelres,
                         residuals = if (M == "ml") ml.residuals,
                         ml.value =  param.table[[M]][tblidx[["ml"]][1]]
                         )
      }
    } # else  res[[M]] <- NA
  }
#  names(res) <- names(param.table)


  RFparameters(old.param)
  options(save.options)

  lower <- transform(lower)
  upper <- transform(upper)
  idx <- lower == upper
  lower[idx] = upper[idx] = NA

  .C("PutValuesAtNA", as.double(lower),
     PACKAGE="RandomFields", DUP=FALSE, NAOK=TRUE)  
  lower <- GetModel(modelreg=ModelNr, modus=1)
  .C("PutValuesAtNA", as.double(upper),
    PACKAGE="RandomFields", DUP=FALSE, NAOK=TRUE)  
  upper <- GetModel(modelreg=ModelNr, modus=1)

   
  return(c(list(ev = ev,
                table=as.matrix(param.table),
                lowerbounds=lower,
                upperbounds=upper,
                transfrom = transform,
                vario = "'$vario' is defunctioned. Use '$ml' instead of '$value$ml'!"
                ),
           res 
           ))
}

  
