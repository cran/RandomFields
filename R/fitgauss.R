### !!!!!!!!!!!! ACHTUNG !!!!!!!!!!!! TREND als cov-fct muss
### noch programmiert werden !!!

# RFsimulate:  Not implemented yet: If \code{model} is a formula or of class
#    \command{\dQuote{\link{RFformula}}},
#    the corresponding linear mixed model of the type
 #   \deqn{response = W*b + Z*u + e} is simulated

##   source("~/R/RF/RandomFields/R/MLES.R")

## Printlevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information


## jetzt nur noch global naturalscaling (ja / nen)
## spaeter eine unktion schreibbar, die den naturscaling umwandelt;
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


## NAs in data mit mixed model grundsaetzlich nicht kombinierbar !
## NAs in data mit trend (derzeit) nicht kombinierbar

## zentrale C -Schnittstellen
##    .C("PutValuesAtNA", Reg, param, PACKAGE="RandomFields")

## bins bei Distances automatisch


## bei repet sind die Trends/fixed effects gleich, es muessen aber die
## random effects unterschiedlich sein.
## bei list(data) werden auch trend/fixed effects unterschiedlich geschaetzt.


## Erweiterungen: Emilio's Bi-MLE, Covarianz-Matrix-INversion per fft oder
## per INLA, grosse Datensaetze spalten in kleinere "unabhaengige".


###################################
## !!! Mixed Model Equations !!! ##
###################################


list2RMmodelFit <- function(x, isRFsp=FALSE,
                            coords, gridTopology, data.RFparams, T) {
  stopifnot(is.list(x),
            all(c("model", "trend", "ml.value", "residuals") %in% names(x)))

  if (isRFsp) {
    stopifnot(!missing(coords) &&
              !missing(gridTopology) &&
              !missing(data.RFparams))

    ## convert residuals to RFsp class
    err <-
      try({
        lres <- length(x$residuals)
        if (lres > 0) {
          for (i in 1:lres) {
            gT <- if (length(gridTopology) < i) NULL else gridTopology[[i]]
            co <- if (length(coords)<i) NULL else coords[[i]]
            if (!is.null(x$residuals[[i]]))
                x$residuals[[i]]<-
                  conventional2RFspDataFrame(data=x$residuals[[i]],
                                             coords=co,
                                             gridTopology=gT,
                                             n=data.RFparams[[i]]$n,
                                             vdim=data.RFparams[[i]]$vdim,
                                             T = T,
                                             vdim_close_together=FALSE)
          }
        }
      }
          , silent=!TRUE)
    if(class(err)=="try-error")
      warning(paste("residuals could not be coerced to class'RFsp';",
                    err))
  }
    
  return(new(ZF_MODELEXT,
             list2RMmodel(x$model),
             trend = list2RMmodel(x$trend),
             variab = x$variab,
             param = x$param,
             covariab = x$covariab,
             likelihood = x$ml.value,
             AIC = x$AIC,
             AICc = x$AICc,
             BIC = x$BIC,
             residuals = x$residuals))
}
  


## used by RFratiotest, fitgauss, Crossvalidation
StandardizeData <- function(x, y, z, T, grid, data, distances, dimensions,
                            RFopt) {
  if (missing(grid)) grid <- NULL
  dist.given <- !is.null(distances)
  matrix.indep.of.x.assumed <- FALSE  
  rangex <- neu <- gridlist <- RFsp.coord <-  gridTopology <- data.RFparams <-
    NULL

  #Print(x, y, z)
  
  if (isSpObj(data)) data <- sp2RF(data)
  if (isRFsp <- is(data, "RFsp")){# ||(is.list(data) && is(data[[1]], "RFsp")))
    stopifnot(missing(x) || is.null(x),
              is.null(y), is.null(z), is.null(T),
              is.null(distances))
    if (!is.list(data)) data <- list(data)
    sets <- length(data)
 
    x <- T <- gridlist <- RFsp.coord <- gridTopology <- data.RFparams <-
      vector("list", sets)
    

    dimdata <- NULL
    for (i in 1:length(data)) {      
      gridtmp <- isGridded(data[[i]])
      compareGridBooleans(grid, gridtmp)
      data.RFparams[[i]] <- data[[i]]@.RFparams
      gridTopology[[i]] <- if (gridtmp) data[[i]]@grid else NULL
      RFsp.coord[[i]] <- if (!gridtmp) data[[i]]@coords else NULL
      dims <- if (gridtmp) data[[i]]@grid@cells.dim else nrow(data[[i]]@data)
      dims <- c(dims, data[[i]]@.RFparams$vdim)
      if (RFopt$general$vdim_close_together) dims <- rev(dims)
      dimdata <- rbind(dimdata, c(dims, data[[i]]@.RFparams$n))
      gridlist[[i]] <- gridtmp
      tmp <- rfspDataFrame2conventional(data[[i]])
      x[[i]] <- tmp$x
      if (!is.null(tmp$T)) T[[i]] <- tmp$T
      data[[i]] <- as.matrix(tmp$data)      
    }

    idx <- if (RFopt$general$vdim_close_together) 1 else length(dims) # vdim
    if (all(dimdata[, idx] == 1))
      dimdata <- dimdata[, -idx, drop=FALSE]
    if (all(dimdata[, ncol(dimdata)] == 1)) # repet
      dimdata <- dimdata[, -ncol(dimdata), drop=FALSE]
    distances <- NULL
  } else { # !isRFsp
    ## dimdata wird spaeter bestimmt
    if (dist.given) {
      stopifnot(missing(x) || is.null(x), is.null(y), is.null(z))
      if (!is.list(distances)) {
        distances <- list(distances)
        if (is.list(data))
          stop("if list of data is given then also for distances ")
        data <- list(as.matrix(data))
      } else if (!is.list(data)) {
        stop("if list of distances is given then also for data ")
        if (length(data) != length(distances))
          stop("length of distances does not match length of data")
      }      
      for (i in 1:length(distances)) {
        if (any(is.na(data)))
          stop("missing data are not allowed if distances are used.")
      }
      spdim <- tsdim <- as.integer(dimensions)
      stopifnot(missing(T)  || is.null(T))
      lcc <- sapply(distances, function(x) 0.5 * (1 + sqrt(1 + 8 * length(x))) )
      lc <- as.integer(lcc)
      if (any(lc != lcc)) stop("number of distances not of form k(k-1)/2")
      coord <- distances
      coord.units <- RFopt$coords$coordunits
    } else { ## distances not given
      distances <- NULL
      if (missing(x)) { ## dec 2012: matrix.indep.of.x.assumed
        if (!is.null(dnames <- dimnames(data)[[2]])) {
          if ((length(xi <- earth_coordinate_names(dnames)) == 2) ||
              (length(xi <- cartesian_coordinate_names(dnames)) > 0)) {
            x <- data[ , xi, drop=FALSE]
            data <- data[ , -xi, drop=FALSE]
          }
        }
        if (missing(x)) {
          matrix.indep.of.x.assumed <- TRUE
          x <- 1:nrow(as.matrix(data))
          x[1] <- 0 ## so no grid !
        }
      } #else {
      if (is.data.frame(x)) x <- as.matrix(x)
      if (is.data.frame(data) || is.vector(data)) data <- as.matrix(data)
      if (is.list(x)) {
        if (!is.null(y) & !is.list(y) |
            !is.null(z) & !is.list(z) |
            !is.null(T) & !is.list(T) |
            !is.list(data)) {
          stop("either all coordinates and data must be given by a list or none.")
        }
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
          stopifnot(!is.list(T))
          T <- list(T)
        }
        stopifnot(!is.list(data))
        data <- list(as.matrix(data))
      }
      ##}
    }
    sets <- length(data)
    dimdata <- matrix(nrow=sets, ncol=length(dim(data[[1]])))
    for (i in 1:sets) dimdata[i, ] <- dim(data[[i]])
    gridlist <- rep(list(grid), sets)

    #Print(neu, x, y, z, T, gridlist)
     # xxxx
  } # !isRFsp
  
  if (!dist.given) { ##   x coordinates, not distances
    coord <- neu <- list()
    for (i in 1:length(x)) {
      
      neu[[i]] <- CheckXT(x[[i]], y[[i]], z[[i]], T[[i]],
                       grid = gridlist[[i]],
                       length.data=length(data[[i]]),
                       printlevel = 0)

      if (i==1) {
        spdim <- as.integer(neu[[i]]$spatialdim)
        tsdim <- as.integer(spdim + !is.null(T[[i]]))
        if (!missing(dimensions) && !is.null(dimensions)) {          
          stop("dim should be given only when 'distances' are given")
        }
         rangex <- as.matrix(apply(neu[[i]]$x, 2, range))       
      } else if (spdim != neu[[i]]$spatialdim |
                 tsdim != as.integer(spdim + !is.null(T[[i]])))
        stop("dimensions of all sets must be the same")

      coord[[i]] <- neu[[i]]$x # stehend
      
      if (neu[[i]]$grid) {
        ## note: checkxt returns always gridtriple format !
 
        coord[[i]] <- apply(coord[[i]], 2, function(z) z[1] +  z[2] * (1:z[3]))
        args <- lapply(if (is.list(coord[[i]])) coord[[i]] else
                      apply(coord[[i]], 2, list),
                      function(z) c(z, recursive=TRUE))
        coord[[i]] <- as.matrix(do.call(base::expand.grid, args=args))
      }
      if (!is.null(T[[i]])) {
        tt <- T[[i]][1] + T[[i]][2] * (1:T[[i]][3])
        coord[[i]] <- cbind(coord[[i]][rep(1:nrow(coord[[i]]), length(tt)), ],
                        rep(tt, each=nrow(coord[[i]])))
      }
      coord[[i]] <- t(coord[[i]]) ## liegend
      
    }
    # x <- y <- z <- T <- NULL    
    lc <- sapply(coord, ncol)    
    coord.units <- neu[[1]]$coordunits
  }
  
  storage.mode(tsdim) <- "integer"
  storage.mode(spdim) <- "integer"

  return(list(coord=coord, data=data, dimdata=dimdata, isRFsp = isRFsp,
              RFsp.coord = RFsp.coord,
              distances = distances,
              neu = neu,
              dist.given=dist.given,
              gridTopology = gridTopology,
              data.RFparams = data.RFparams,
              lc=lc, spdim=spdim, tsdim=tsdim,
              rangex = rangex, coord.units=coord.units,
              matrix.indep.of.x.assumed = matrix.indep.of.x.assumed))
}



GetLowerUpper <- function(lower, upper, trafo, optim_var_elimination,
                          sillbounded, var.idx, nugget.idx, sill,
                          nonugget, varmin, varmax) {
  low <- trafo(lower)
  up <- trafo(upper)
  if (optim_var_elimination == "yes") {
    if (sillbounded) {
      low[c(var.idx, nugget.idx)] <- 0
      up[c(var.idx, nugget.idx)] <- sill
    } else {
      if (length(var.idx) == 1) {
        if (length(nugget.idx) == 1) {
          idx <- c(var.idx, nugget.idx)
          low[idx] <- 0
          up[idx] <- varmax[idx]
        } else {
          if (nonugget) {
            low[var.idx] <- 0
            up[var.idx] <- varmin[var.idx]
          }
        }
      }
    }
  }
  return(list(lower=low, upper=up))
}



vary.variables <- function(variab, lower, upper) {
  n.var <- length(lower)
  w.value <- runif(n.var, 3-0.5, 3+0.5) # weight
  n.modivar <- 5

  idx <- lower <= 0
  modilow <- modiupp <- numeric(n.var)
  modilow[idx] <- (lower[idx] + w.value[idx] * variab[idx]) / (w.value[idx] + 1)
  modilow[!idx] <-(lower[!idx]*variab[!idx]^w.value[!idx])^(1/(w.value[!idx]+1))
  modiupp[idx] <- (upper[idx] + w.value[idx] * variab[idx]) / (w.value[idx] + 1)
  modiupp[!idx] <-(upper[!idx]*variab[!idx]^w.value[!idx])^(1/(w.value[!idx]+1))
  
  modivar <- matrix(nrow=n.var, ncol=n.modivar)
  for (i in 1:n.var) {
    modivar[i,] <- seq(modilow[i], modiupp[i], length=n.modivar)
  }
  return(modivar)
}


vary.x <- function(rangex) {
  ## letzter Eintrag ist 0
  npt <- 5
  totdim <- ncol(rangex)
  x <- matrix(0,ncol=npt+1, nrow=totdim)
  for (i in 1:totdim) {
    x[i, 1:npt] <- runif(npt, rangex[1, i], rangex[2, i])
  }
  x
}


GetTrend <- function(coord, idx, model, vdim, param) {
  ## param from modelnfo !!
  ## coord hat Sequenz stehender Vektoren
  lc <- ncol(coord)
  if (!missing(param) &&
      all(names(param) %in% c("mean", "plane"))) {
    ret <- NULL
    if (length(param$mean)>0) {
      stopifnot(is.na(param$mean))
      ret <- matrix(0, nrow=length(idx), ncol=vdim)
      cum <- 0
      for (v in 1:vdim) {
        current <- length(idx <= lc)
        ret[cum + (1:current), v] <- 1
        idx <- idx[!current] - lc
        cum <- cum + current
      }
    }

    if (length(param$plane)>0) {
      stopifnot(all(is.na(param$plane)))
      ## Dimensionen passen hier nicht ! db muss length(idx) x (vdim * xdimOZ)
      ## sein, oae. Sonst passt es nicht mit coord[, idxLlx] zusammen.
      stop("plane not programmed yet")
      dB <- matrix(0, nrow=length(idx), ncol=vdim)
      cum <- 0
      for (v in 1:vdim) {
        idxLlc <- idx <= lc
        current <- length(idxLlc)
        dB[cum + (1:current), v] <- t(coord[, idxLlc, drop=FALSE])
        idx <- idx[!current] - lc
        cum <- cum + current
      } 
      ret <- cbind(ret, dB)
    }
    return(ret)
  } else {
    if (length(param$polydeg) > 0)
      stop("poly not programmed yet")
    if (length(param$arbitraryfct)>0)
      stop("arbitrary fct not programmed yet")
    stop("unknown parameter choice for trend")

    ## DO NOT DELETE
#    trendModel <-
 #     rfSplitTrendAndCov(model=model, xdimOZ????,
  #                       cathegories=list(det=c(FixedTrendEffect,
  #                                          FixedEffect)))
   # trendinfo <- .Ca ll("SetAndG ?? etModelInfo", model, dimXXXX,
    #                   dimXXX, as.integer(short), FALSE, HilfsReg, TRUE,
     #                  PACKAGE="RandomFields")

    
  }
}

GetValuesAtNA <- function(NAmodel, valuemodel, spdim, Time, shortnamelength,
                          skipchecks) {

  spdim <- as.integer(spdim)
  Time <- as.logical(Time)
  skipchecks <- as.logical(skipchecks)
  shortnamelength <- as.integer(shortnamelength)
  
  reg <- MODEL_SPLIT

  info <- neu <- list()
  models <- list(NAmodel, valuemodel)
  #Print(models)
  for (m in 1:length(models)) {
    
    info[[m]] <- .Call("SetAndGetModelInfo", reg, list("Cov", models[[m]]),
                  spdim, FALSE, # distances
                  TRUE, # in case it is a non-trans.inv model
                  Time, spdim,
                  shortnamelength,
                  FALSE, # allowforinteger
                  TRUE,  # excludetrend
                  PACKAGE="RandomFields")
    neu[[m]] <- GetModel(register=reg, modus=GETMODEL_DEL_MLE, spConform=FALSE,
                         do.notreturnparam=TRUE)
    
  }
 # Print(neu)

  NAs <- nrow(info[[1]]$minmax)
  ret <- .Call("Take2ndAtNaOf1st", reg, list("Cov", neu[[1]]),
               list("Cov", neu[[2]]), spdim, Time, spdim, 
               NAs, skipchecks, PACKAGE="RandomFields")
  
  return(ret)
}


ParScale <- function(optim.control, current, lower, upper) {
  if (!is.null(parscale <- optim.control$parscale)) {
    if (!is.numeric(parscale)) {
      stop("non numeric parscale not allowed yet")
##    parscale <- GetValuesAtNA(NAmodel=model, valuemodel=parscale, spdim=spdim,
#                                Time=time,
#                                shortnamelength=fit$short, 
#                                skipchecks=!is.null(trafo))
    } # else is.numeric: nothing to do
  } else parscale <- rep(NA, length(lower))

  if (!missing(current)) {
    idx <- is.finite(parscale)
    parscale[!idx] <- current[!idx]
    parscale[idx] <- sqrt(parscale[idx] * current[idx]) # geom mittel
  }
  
  if (any(idx <- !is.finite(parscale) | parscale == 0)) {
    parscale[idx] <- pmax(abs(lower[idx]), abs(upper[idx])) / 10 # siehe fit_scale_ratio unten
    idx <- idx & (lower * upper > 0)
    parscale[idx] <- sqrt(lower[idx] * upper[idx])
    stopifnot(all(is.finite(parscale)))
  }

#  Print("Parscale", lower, upper, parscale)
  
  return(parscale)
}


ModelSplitXT <- function(Reg, info.cov, trafo, variab,
                         lower, upper, rangex, modelinfo, model,
                         p.proj=1:length(variab), v.proj=1,
                         report.base) {
  
  vdim <- modelinfo$vdim
  tsdim <- modelinfo$tsdim
  ts.xdim <- modelinfo$ts.xdim
  xdimOZ <- modelinfo$xdimOZ
  is.dist.vector <- modelinfo$is.dist.vector


  if (ts.xdim == 1) return(list(p.proj=p.proj,
        x.proj = TRUE,
        v.proj=v.proj))
  lp <- length(variab[p.proj]) # not length(p.proj) sicherheitshalber
  statiso <- info.cov$trans.inv && info.cov$isotropic
   
  abs.tol <- 0

  truely.dep <- dep <- matrix(FALSE, nrow=lp, ncol=ts.xdim)
  varyx <- vary.x(rangex=rangex)
  varyy <- vary.x(rangex=rangex)

  modivar <- vary.variables(variab=variab, lower=lower, upper=upper)
    
  for (d in 1:ts.xdim) {
    x <- varyx
    x[-d, ] <- 0
    if (is.dist.vector) y <- NULL
    else {
      y <- varyy
      y[-d, ] <- 0
    }

    nonzeros <- integer(1)
    S <- double(ncol(x) * vdim^2)
    for (i in 1:lp) {

      for (m0 in 1:ncol(modivar)) {
        var <- modivar[, m0]
      
        for (m1 in 1:ncol(modivar)) {
          var[p.proj[i]] <- modivar[i, m1]

          .C("PutValuesAtNA", as.integer(Reg),
             as.double(trafo(var)), PACKAGE="RandomFields")          
          .Call("CovLoc", Reg, x, y, xdimOZ, ncol(x), S, 
                PACKAGE="RandomFields")
          dim(S) <- c(ncol(x), vdim, vdim)
         
          if (any(is.na(S)))
            stop("Model too complex to split it. Please set split=FALSE")

          ## Bei x_j=0, j!=i, gibt es Auswirkungen hinsichtlich
          ## des Parameters auf die Kov-Fkt?
          ## d.h. aendern sich die kov-Werte bei variablem Parameter,
          ## obwohl die x_j = 0 sind?
          if (m1==1) Sstore <- S[, v.proj, v.proj, drop=FALSE]
          else {
            difference <- S[,v.proj, v.proj, drop=FALSE] - Sstore
            if (any(abs(difference) > abs.tol)) {
              dep[i, d] <- TRUE

              ## ist die Kovarianzfunktion der Bauart C_1(x_1) + v C_2(x_2),
              ## und d=1, so bildet v eine Konstante bzg. C_1(x_1)        

              ## i.e.,
              ## note: s_1 C_1(x) + s_2 C_2(t)
              ## then s_2 dep(ends) on x, but not truely
              truely.dep[i, d] <- truely.dep[i, d] || 
                 !(all(apply(difference, 2:3, function(x) all(diff(x) == 0))))
            }
          }
        }
      }
    }

  }

  modelsplit <- NULL
  untreated_x <- rep(TRUE, ts.xdim)
  untrtd_p <- rep(TRUE, lp)
  
  for (d in 1:ts.xdim) {
    if (!untreated_x[d]) next
    rd <- dep[, d]
    if (!any(rd)) next
    same.class <- which(apply(dep == rd, 2, all))
    untreated_x[same.class] <- FALSE

    ## untreated_x wird false falls rd und der Rest keine Abhaengigkeiten
    ## mehr gemein haben
    untreated_x[d] <- !all(!rd | !dep[, -same.class, drop=FALSE])
 
    ## falls es Parameter gibt, die nur zur Koordinate d gehoeren,
    ## werden diese geschaetzt, zusammen mit allen anderen Paremetern,
    ## die zur Koordinate d (und zu anderen Koordinaten) gehoeren 
    trtd <- apply(rd & !dep[, -same.class, drop=FALSE], 1, all)
    if (any(trtd)) {
      untrtd_p <- untrtd_p & !trtd
      report <- "seperable"
      modelsplit[[length(modelsplit) + 1]] <-
        list(p.proj = p.proj[dep[, d]], 
             v.proj = v.proj,
             x.proj=as.integer(same.class),
             report = (if (missing(report.base)) report else
                       paste(report.base, report, sep=" : "))
             )
    }
  }

#  Print("A", modelsplit)
  
  # !truely.dep, die moegicherweise rausgeflogen waren
  nontruely <- apply(!truely.dep & dep, 1, any)
  if (any(nontruely)) {
    
    if (sum(nontruely) > 1) stop("more than 1 parameter is detected that is not truly dependent -- please use split=FALSE")
      ## Problem ist hier die nicht Eindeutigkeit bei
    ## C = v_1 C_1(x_1) + v_2 C_2(X_2) + v_2 C_3(x_3)
    ## da nun v_2 und v_3 als nontruely auftauchen und somit
    ## dass auf x_1 projezierte Modell nicht mehr identifizierbar ist.
    ## Letztendlich muessen die beiden (oder mehrere) zusammengefasst werden
    ## und in recursive.estimation, use.new.bounds adaequat beruecksichtigt
    ## werden 
    l <- length(modelsplit)
    modelsplit[[l + 1]] <- list(use.new.bounds=p.proj[nontruely])
    modelsplit[[l + 2]] <-
      list(p.proj = p.proj[nontruely], v.proj = v.proj,
           x.proj=which(apply(nontruely & dep, 2, any))
           ) ## oder vielleicht doch TRUE oder durchzaehlen fuer
    ##                     verschiedene Ebenen
    untrtd_p <- untrtd_p & !nontruely
  }

  ## Parameter die verwoben sind:
  if (any(untrtd_p & !nontruely)) {
    l <- length(modelsplit)
    report <- "simple space-time"
    modelsplit[[l + 1]] <- list(use.new.bounds = which(untrtd_p))
    modelsplit[[l + 2]] <-
      list(p.proj=which(untrtd_p), v.proj=v.proj,
           x.proj=which(apply(untrtd_p & dep, 2, any)),
           report = (if (missing(report.base)) report else
                     paste(report.base, report, sep=" : ")))
  }

  ## komplett, falls notwendig
  if ((any(nontruely) || any(untrtd_p)) && !all(untrtd_p) && !all(nontruely)) {
    l <- length(modelsplit)
    modelsplit[[l + 1]] <- list(use.new.bounds=p.proj)
    modelsplit[[l + 2]] <- list(p.proj=p.proj, x.proj=1:ts.xdim, v.proj=v.proj)
  }

  if (length(modelsplit) == 1) {
    modelsplit <- modelsplit[[1]]
    modelsplit$report <- NULL
  }

#  Print("B", modelsplit); dsffas

  return(modelsplit)
}


ModelSplit <- function(Reg, info.cov, trafo, variab,
                       lower, upper, rangex, modelinfo, model) {
  vdim <- modelinfo$vdim
  tsdim <- modelinfo$tsdim
  ts.xdim <- modelinfo$ts.xdim
  is.dist.vector <- modelinfo$is.dist.vector
  xdimOZ <- modelinfo$xdimOZ
  refined <- modelinfo$refined

  abs.tol <- 0
  restrictive <- FALSE

  lp <- length(variab)
  statiso <- info.cov$trans.inv && info.cov$isotropic
 # if (xor(is.dist, statiso)) stop("mismatch in ModelSplit -- contact author")

  if (vdim == 1) {

    if (statiso) return(NULL)
    return(ModelSplitXT(Reg=Reg, info.cov=info.cov, trafo=trafo,
                        variab=variab, lower=lower, upper=upper, rangex=rangex,
                        modelinfo=modelinfo, model=model))
  }
  
  
  overlap <- dep <- matrix(FALSE, nrow=lp, ncol=vdim)

  x <- vary.x(rangex=rangex)
  S <- double(ncol(x) * vdim^2)
  if (!is.dist.vector) y <- vary.x(rangex=rangex)

  modivar <- vary.variables(variab=variab, lower=lower, upper=upper)

  for (i in 1:lp) {

    for (m0 in 1:ncol(modivar)) {
      var <- modivar[, m0]
      
      for (m1 in 1:ncol(modivar)) {
        var[i] <- modivar[i, m1]

       .C("PutValuesAtNA", Reg, trafo(var), PACKAGE="RandomFields")
        .Call("CovLoc", Reg, x, if (is.dist.vector) NULL else y,
              xdimOZ, ncol(x), S, PACKAGE="RandomFields")
        dim(S) <- c(ncol(x), vdim, vdim)
        if (any(is.na(S)))
          stop("model too complex to split it. Please set split=FALSE")

         if (m1==1) Sstore <- S
        else {
          for (n in 1:vdim)
            dep[i, n] <-
              dep[i, n] | any(abs(S[,n,n] - Sstore[,n,n]) > abs.tol)
        }
      }
    }
  }

  modelsplit <- list()
  untreated_v <- logical(vdim)
  report <- "indep. multivariate"
  
  for (v in 1:vdim) {
    d1 <- dep[, v]
    if (!any(d1)) next
    overlap[, v] <- apply(dep[, -v, drop=FALSE] & d1, 1, any)
    untreated_v[v] <- if (restrictive) !any(overlap[, v])
             else (sum(overlap[, v]) < sum(d1))
    if (untreated_v[v]) {
      idx <- length(modelsplit) +1
       modelsplit[[idx]] <-
        ModelSplitXT(Reg=Reg, info.cov=info.cov, trafo=trafo, variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj=which(d1), v.proj = v, report.base=report)
      modelsplit[[idx]]$report <- report
    }
  }
 
  ## bislang oben nur auf univariat heruntergebrochen
  ## man koennte noch auf bivariate runterbrechen
  ## dies aber erst irgendwann;
  ## wuerde ich auch zu weiteren report-Ebenenen 2,3, etc. fuehren
  
  if (any(untreated_v)) {    
    overlapping <- apply(overlap[, untreated_v, drop=FALSE], 1, any)
    depending <- apply(dep[, untreated_v, drop=FALSE], 1, any)

    #Print(refined, any(depending),  any(overlapping))
    
    if (refined && any(depending) && any(overlapping)) {
      idx <- length(modelsplit)
      depend <- which(depending | overlapping)
      notok <- which(!depending & !overlapping)
      modelsplit[[idx + 1]] <- list(use.new.bounds=depend, fix=notok)
      report <- "simple multivariate"
      modelsplit[[idx + 2]] <-
        ModelSplitXT(Reg=Reg, info.cov=info.cov, trafo=trafo, variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj = depend, v.proj = 1:vdim, report.base=report)
      modelsplit[[idx+2]]$report <- report
    } else notok <- which(!depending | overlapping)
 
    if (length(notok) > 0) {
      idx <- length(modelsplit) + 1
      modelsplit[[idx]] <-
        ModelSplitXT(Reg=Reg, info.cov=info.cov, trafo=trafo, variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj = notok, v.proj = 1:vdim)
    }

    l <- length(modelsplit)
    modelsplit[[l+1]] <- list(use.new.bounds=1:lp)
    modelsplit[[l+2]] <-
      list(p.proj=1:lp,
           x.proj=if (ts.xdim==1 && tsdim > 1) TRUE else 1:ts.xdim,
           v.proj=1:vdim)
  }

  modelsplit[[1]]$all.p <- dep

  return(modelsplit) ## return(NULL) # 11.9.12
}



recurs.estim <- function(split, level, Reg, vdim, lc,
                         model, x, T, data,
                         lower, upper,
                         users.lower, users.upper,
                         guess,  
                         distances,
                         grid,
                         bc_lambda,
                         lsq.methods, mle.methods,
                         tsdim,
                         optim.control,
                         transform, trafo,
                         spConform,
                         practicalrange,
                         printlevel,
                         minmax,
                         sdvar
                         ) {
  M <- if (length(mle.methods) >= 1) mle.methods[length(mle.methods)]
       else lsq.methods[length(lsq.methods)]
  
#  Print(lsq.methods, mle.methods, M); kkkk
  
  w <- 0.5
  dist.given <- missing(x) || is.null(x)
  sets <- length(data)

  if (printlevel >= PL_FCTN_DETAILS) Print(split) #
  submodels <- NULL
  submodels_n <- 0
  fixed <- FALSE

  for (s in 1:length(split)) {
    sp <- split[[s]]
    if (!is.null(p <- sp$use.new.bounds)) {
      if (printlevel >= PL_STRUCTURE) {
        cat("    calculating new lower and upper bounds for the parameters ",
            paste(p, collapse=", "), sep="")
        cat("\n")
      }
       
      if (!is.null(sp$fix)) {
        fix.zero <- (minmax[sp$fix, 1] <= 0) & (minmax[sp$fix, 2] >=0)
        fix.one <- (minmax[sp$fix, 1] <= 1) & (minmax[sp$fix, 2] >= 1)
        if (!all(fix.zero | fix.one))
          stop("Some parameters could not be fixed. Set 'split_refined = FALSE")
        fix.one <- sp$fix[fix.one & !fix.zero] ### zero has priority
        fix.zero <- sp$fix[fix.zero]
        
        guess[fix.one] <- 1
        guess[fix.zero] <- 0
        
        fixed <- TRUE ## info for the next split[[]] that some variables
        ##               have been fixed. The values must be reported in the
        ##               intermediate results for AIC and logratio test
      }
      next
    }
    if (is.list(sp[[1]])) {
      res <- recurs.estim(split=sp, level=level+1, Reg=Reg,
                          vdim=vdim, lc=lc,
                          model=model, x=x, T=T, 
                          data= data,
                          lower=lower,
                          upper=upper,
                          users.lower=users.lower,
                          users.upper=users.upper,
                          guess= guess,  
                          distances=distances,
                          grid=grid,
                          bc_lambda=bc_lambda,
                          lsq.methods=lsq.methods,
                          mle.methods=mle.methods,
                          tsdim=tsdim,
                          optim.control=optim.control,
                          transform=transform,
                          trafo=trafo,
                          spConform = spConform,
                          practicalrange = practicalrange,
                          printlevel=printlevel,
                          minmax = minmax,
                          sdvar = sdvar
                          )
      if (!is.null(sp$report)) {
        submodels_n <- submodels_n + 1
        if (fixed) res$fixed <- list(zero=fix.zero, one=fix.one)
        sumodels[[submodels_n]] <- res
      }
    } else { # !is.list(sp[[1]])
      if (printlevel>=PL_RECURSIVE) {
        cat("splitting (",
            paste(rep(". ", level), collapse=""), format(s, width=2),
            ") : x-coord=", paste(sp$x, collapse=","),
            "; compon.=", paste(sp$v, collapse=","), sep="")
        if (printlevel>=PL_REC_DETAILS)
          cat( "; parameters=", paste(sp$p, collapse=", "), sep="")
        cat(" ")
      }

      guess <- pmax(lower, pmin(upper, guess, na.rm=TRUE), na.rm=TRUE)
          
      .C("PutValuesAtNAnoInit",
         Reg, trafo(guess), PACKAGE="RandomFields", NAOK=TRUE) # ok
      old.model <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE ,
                            do.notreturnparam=TRUE)
      
      guess[sp$p.proj] <- NA
      .C("PutValuesAtNAnoInit",
         Reg, trafo(guess), PACKAGE="RandomFields", NAOK=TRUE) # ok
      new.model <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE,
                            spConform=FALSE, do.notreturnparam=TRUE)

      if (!all(is.na(lower))) {
        .C("PutValuesAtNAnoInit",
           Reg, trafo(lower),PACKAGE="RandomFields", NAOK=TRUE) # ok
        lower.model <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE,
                                spConform=FALSE, do.notreturnparam=TRUE)
      } else lower.model <- NULL
      
     if (!all(is.na(upper))) {
       .C("PutValuesAtNAnoInit",
          Reg,  trafo(upper), PACKAGE="RandomFields", NAOK=TRUE) # ok
       upper.model <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE,
                               spConform=FALSE, do.notreturnparam=TRUE)
     } else upper.model <- NULL
  

 #     Print(new.model, lower.model, upper.model); 
      
      if ((new.vdim <- length(sp$v.proj)) < vdim) {
        m <- matrix(0, nrow=new.vdim, ncol=vdim)
        for (j in 1:new.vdim) m[j, sp$v.proj[j]] <- 1
        new.model   <- list("M", M=m, new.model)
        lower.model <- list("M", M=m, lower.model)
        upper.model <- list("M", M=m, upper.model)
        old.model   <- list("M", M=m, old.model)
        new.data <- list()
        for (i in 1:sets) {
          new.data[[i]] <- data[[i]][ , sp$v.proj, , drop=FALSE]
        }
      } else {
        new.data <- data
        vlen <- length(sp$v.proj)
        for (i in 1:length(new.data)) {
          dim(new.data[[i]]) <-
            c(lc[i], vlen, length(new.data[[i]]) / (lc[i] * vlen))          
        }
      }

      
      x.proj <- sp$x.proj
      ignored.x <- if (is.logical(x.proj)) integer(0) else (1:tsdim)[-x.proj]
      if (time <- !is.null(T[[1]])) {
        stopifnot(!is.logical(x.proj))
        ignore.T <- all(x.proj != tsdim)
        x.proj <- x.proj[x.proj != tsdim]
        ignored.x <- ignored.x[ignored.x != tsdim]
        lenT <- sapply(T, function(x) x[3])
      } else {
        lenT <- rep(1, length(new.data))
      }
  
      if (!dist.given && !is.logical(sp$x.proj) && length(sp$x.proj) < tsdim)  {
        if (grid) {
          new.x <- x
          for (i in 1:length(new.data)) {
            len <- new.x[[i]][3, ]

            stopifnot(prod(len) == lc[i])
            d <- dim(new.data[[i]])
            new.d <-  c(len, d[-1])
            dim(new.data[[i]]) <- new.d
                        
            new.data[[i]] <-
              aperm(new.data[[i]], c(x.proj, length(new.d) + (-1:0), ignored.x))
            repet <- d[3] * prod(len[ignored.x]) ## = repet * ignored
            dim(new.data[[i]]) <- c(prod(len[x.proj]), d[2], repet)
            new.x[[i]][, ignored.x] <- c(0, 0, 1)
          }
        } else { ## not grid
          dummy <- new.data
          new.x <- new.data <- list()
          xlen <- sapply(x, function(y) nrow(y))
          for (i in 1:length(dummy)) {
            d <- dim(dummy[[i]])
            dim(dummy[[i]]) <-
              c(xlen[i], lenT[i], length(dummy[[i]]) / (xlen[i] * lenT[i])) 
            xyz <- x[[i]]
            while(length(dummy[[i]]) > 0) {
              slot <- xyz[1, ignored.x]
              idx <- apply(t(xyz) == slot, 2, all)
              last <- length(new.x) + 1
              new.x[[last]] <- xyz[idx, , drop=FALSE]
              new.x[[last]][, ignored.x] <- 0
              new.data[[last]] <- dummy[[i]][idx, , ,drop=FALSE]
              dummy[[i]] <- dummy[[i]][-idx, , , drop=FALSE]
              xyz <- xyz[-idx, , drop=FALSE]
            }
          }
        }
      } else {
        new.x <- x
      }

      if (!is.null(transform)) {
        isna <- is.na(transform[[2]](guess))
        idx <- transform[[1]][isna]
        stopifnot(max(sp$p.proj) <= length(guess))
                
        f <- function(p) {
          q <- guess
          q[sp$p.proj] <- p
          z <- (transform[[2]](q))[isna]
          z     
        }
      
        
        new.transform <- list(idx, f)
      } else new.transform <- NULL

   #   file <- paste("whole", level, s, "rda", sep=".")
#     file <- paste("pars", level, s, "rda", sep=".")
 #     Print(file)
   
      general_spConform <- if (s<length(split)) FALSE else spConform

      if (TRUE) {
        res <-
          rffit.gauss(model=new.model,
                      x=new.x,
                      T=if (time && !ignore.T) T else NULL,
                      grid=grid,
                      data= new.data, 
                      lower= if (!all(is.na(lower))) lower.model,
                      upper= if (!all(is.na(upper))) upper.model,
                      bc_lambda=bc_lambda,
                      mle.methods=mle.methods,
                      lsq.methods=lsq.methods,                        
                      users.guess=old.model,  
                      distances=if (dist.given) distances,
                      dimensions = if (dist.given) tsdim,
                      optim.control=optim.control,
                      transform=new.transform,
                      recall = TRUE,
                      fit.split = FALSE,
                      general.practicalrange = practicalrange,
                      general.spConform = general_spConform,
                      sdvar = if (length(sp$v.proj) > 1)
                                  sdvar[sp$p.proj, sp$v.proj, drop=FALSE]
                      )

      

        if (!is.null(sp$report)) {
          submodels_n <- submodels_n + 1
          if (general_spConform) {
            res@p.proj <- as.integer(sp$p.proj)
            res@v.proj <- as.integer(sp$v.proj)
            res@x.proj <- sp$x.proj
            res@report <- sp$report
            res@true.tsdim <- as.integer(tsdim)
            res@true.vdim <- as.integer(vdim)
            if (fixed) res@fixed <- list(zero=fix.zero, one=fix.one)
          } else {
            res$p.proj <- as.integer(sp$p.proj)
            res$v.proj <- as.integer(sp$v.proj)
            res$x.proj <- sp$x.proj
            res$report <- sp$report
            res$true.tsdim <- as.integer(tsdim)
            res$true.vdim <- as.integer(vdim)
            if (fixed) res$fixed <- list(zero=fix.zero, one=fix.one)
          }
          submodels[[submodels_n]] <- res
        }
        
      } else {
        stop("forbidden area. contact author")
      }

      table <- if (general_spConform) res@table else res$table  
      np <- length(sp$p.proj)
      
  #    Print("end", lower, upper, table, M, np)
      lower[sp$p.proj] <- pmin(na.rm=TRUE, lower[sp$p.proj],
                               table[[M]][(np + 1) : (2 * np)])
      upper[sp$p.proj] <- pmax(na.rm=TRUE, upper[sp$p.proj],
                               table[[M]][(2 * np + 1) : (3 * np)])

      
      if (s==length(split)) {
        if (submodels_n >= 1) {
          if (general_spConform) res@submodels <- submodels
          else res$submodels <- submodels
        }
        return(res)
      }
    } # else, (is.list(sp[[1]])) 
 
    guess[sp$p.proj] <- res$table[[M]][1:length(sp$p.proj)]    
    fixed <- FALSE
     
  } # for s in split
  
  return(res)
} # fit$split



######################################################################


rffit.gauss <- function(model, x, y=NULL, z=NULL, T=NULL, grid, data, 
           lower=NULL, upper=NULL,
           bc_lambda, ## if missing then no BoxCox-Trafo
           mle.methods=MLMETHODS,
           lsq.methods=if (variogram.known) LSQMETHODS else NULL,
           ## "internal" : name should not be changed; should always be last
           ##              method!
           users.guess=NULL,  
           distances=NULL,
           dimensions, ## space time dim -- later a vector of grouped dimensions
           optim.control=NULL,
           transform=NULL,
           ##type = c("Gauss", "BrownResnick", "Smith", "Schlather",
           ##             "Poisson"),
           recall = FALSE,
           sdvar = NULL,
           ...) {
  ##         BoxCox ((x+c)^l - 1) / l; log(x+c); with c==1
  ##         bc_lambda : NA / value
  ##         bc_c: nur value
  ##       cross.methods=NULL,
  ##      cross.methods=c("cross.sq", "cross.abs", "cross.ign", "cross.crps"),

  ## ACHTUNG: durch rffit.gauss werden neu gesetzt:
  ##    practicalrange, optim_var_elimination
  ##


    #Print(x, y); ooo

  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  if (RFopt$general$modus_operandi == MODENAMES[normal + 1] &&
      RFopt$internal$warn_normal_mode) {
    RFoptions(internal.warn_normal_mode = FALSE)
    message("The modus_operandi='", MODENAMES[normal + 1], "' is save, but slow. If you like the MLE running\nfaster (at the price of being slightly less save) choose mode='easygoing' or\neven mode='", MODENAMES[sloppy + 1], "'.")
  }

  debug <- FALSE
  
 
  cross.methods <- NULL
  var.name <- "X"
  time.name <- "T"
 # orig.lower <- lower
 # orig.upper <- upper
  no.warning <- -1
#
  no.warning <- 2
  nonzeros <- integer(1)
  
 
#  i <- as.integer(.Machine$integer.max / spam.n[length(spam.n)])
#  if (i %% 2 == 0) i <- i-1
#  repeat {
#    if (any(i %% (3: as.integer(sqrt(i))) == 0))
#      cat(i, "not a prime\n") else { cat(i, "prime\n"); break;}
#    i <- i - 2
#  }



  
####################################################################
###                      Preludium                               ###
####################################################################

  save.options <- options()
  on.exit(options(save.options), add=TRUE)
  
  general <- RFopt$general
  pch <- general$pch
  printlevel <- general$printlevel
  if (printlevel < 1) pch <- ""
   
  fit <- RFopt$fit
  sill <- fit$sill
  optim_var_elimination <- orig.optim_var_elimination <-
    fit$optim_var_elimination
  solvesigma <- fit$solvesigma ## ?? noch was zu tun ? NA?
  bins <- fit$bins
  nphi <- fit$nphi
  ntheta <- fit$ntheta
  ntime <- fit$ntime
  use_spam <- fit$use_spam
  optimiser <- fit$optimiser

#  if (length(cross.methods) > 0) {
#    stop("not programmed yet")
#    if (is.na(optim_var_elimination)) optim_var_elimination <- FALSE
#    else if (!optim_var_elimination)
#      stop("if cross.methods then optim_var_elimination may not be used")
#  }

    
  ##  library(RandomFields, lib="~/TMP")
  if (general$practicalrange && fit$use_naturalscaling)
    stop("practicalrange must be FALSE if fit$use_naturalscaling=TRUE")
  if (fit$use_naturalscaling) RFoptions(general.practicalrange = 3)
 
  if (printlevel>=PL_STRUCTURE) cat("\nfunction defintions...\n")

  ##all the following save.* are used for debugging only
  silent <- printlevel <  PL_STRUCTURE   # optimize
  show.error.message <- TRUE # options

  detailpch <- if (pch=="") "" else '+'
  nuggetrange <- 1e-15 + 2 * RFopt$nugget$tol


### definitions from 

  EffectName <- c("DetTr", "Determ", "FixTr", "FixEff", "RanEff", "RanEff",
                  "XRnEff", "XRnEff", "SpcEff")

  Reg <- MODEL_MLE # GetModelRegister(if (fit$split) "mle" else "mlesplit") 
  TrendReg <- MODEL_MLETREND
  ##
 
  spam.min.n <- as.integer(400) # 400
  spam.tol <- 1e-8
  spam.min.p <- 0.8 # 0.8
  spam.n <- as.integer(1:500)
  spam.factor <- as.integer(4294967)# prime, so that n * factor is still integer
 
  ## Note: upper case variables are global variables !!

  
######################################################################
###                function definitions                            ###
######################################################################
    
  # Rcpp need gradient
  nloptr <- NULL ## just to avoid the warning of R CMD check
  if (optimiser != "optim") {
    # stopifnot(do.call(base::require, list(optimiser, quietly=silent)))
    stop("currently unavailable feature")
    optim.call <- optimiser
    if (optim.call == "minqa") optim.call <- "bobyqa"
    else if (optim.call =="pso") optim.call <- "psoptim"
    if (requireNamespace(optimiser, quietly = TRUE)) {
      optim.call <- do.call("::", list(optimiser, optim.call))
    } else {
      stop("to use '", optimiser, "' its package must be installed")
    }
  } else {
    optim.call <- optim
  }
  
  ## optim : standard
  ## optimx: slower, but better 
  ## soma  : extremely slow; sometimes better, sometimes worse
  ## nloptr: viele algorithmen, z.T. gut
  ## GenSA : extrem langsam
  ## minqa : gut, langsamer
  ## pso   : exxtrem langsam
  ## DEoptim: langsam, aber interessant; leider Dauerausgabe
    
  if (is.na(use_spam) || use_spam) {
    ip <- installed.packages()
    if ("spam" %in% ip) {
      use_spam <- TRUE
    } else {
      if (!is.na(use_spam)) stop("package 'spam' could not be loaded")
      use_spam <- FALSE
    } 
  }


  INVDIAGHESS <- function(par, fn, control=NULL, ...) {
    if (length(par) == 0) return(NULL)
    if (length(idx <- which("algorithm" == names(control))) > 0)
      control <- control[-idx];
    
    oH <- try(optimHess(par=par, fn=fn, control=control), silent=TRUE)
    if (class(oH) != "try-error") {
      zaehler <- 1
      while (zaehler <= 3) {
        e <- eigen(oH)
        values <- Re(e$values) ## some imaginary small values included in eigen
        vectors <- Re(e$vectors)
        if (all(values != 0)) {
          ## solve scheitert am Eigenwert 0
          ##    Print(vectors, values)
          invH <- vectors %*% (1/values * t(vectors))
          diagH <- -diag(invH)        
          if (any(diagH < 0)) {
            diagH[diagH<0] <- Inf
          }
          return(sqrt(diagH))
        }
        oH <- oH + rnorm(length(oH), 0, 1e-7)
      }
    }
    #Print(optimHess(par=par, fn=fn, control=control)) #
    stop("The Hessian matrix is strange and not invertable")
  }

  OPTIM <- function(par, fn, lower, upper, control, optimiser, silent) {
#    print(control)
    try(switch(optimiser,
               "optim" = {
               if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
               optim.call(par=par, fn=fn, lower=lower, upper=upper,
                          control=control, method ="L-BFGS-B")
               },
               "optimx" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=fn, lower=lower, upper=upper,
                            control=control, method ="L-BFGS-B")
               },
               "soma" = {
                 optim.call(function(...) - fn(...),
                            bounds=list(min=lower, max=upper),
                            options=list(), strategy="all2one")
               },
               "nloptr" = {
                 if (length(control$xtol_rel) == 0) control$xtol_rel <- 1e-4
                 optim.call(x0=par, eval_f=fn, lb=lower, ub=upper,
                            opts= control[-pmatch(c("parscale", "fnscale"),
                              names(control))] )
               },
               "GenSA" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=function(...) -fn(...),
                            lower=lower, upper=upper,
                            control=control[-pmatch(c("parscale", "fnscale"),
                              names(control))])
               },
               "minqa" =  {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=function(...) -fn(...),
                            lower=lower, upper=upper,
                            control=control[-pmatch(c("parscale", "fnscale"),
                              names(control))])
               },
               "pso" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=fn, lower=lower, upper=upper,
                            control=control[-pmatch(c("fnscale"),
                              names(control))])
               },
               "DEoptim" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 control <- control[-pmatch(c("parscale", "fnscale"),
                                            names(control))]
                 if (length(control)==0)
                   optim.call(fn=fn, lower=lower,upper=upper)
                 else
                   optim.call(fn=fn, lower=lower,upper=upper, control=control)
               }), silent=silent)
   }

  BoxCox <- function(x, lambda)
     if (abs(lambda) < 10^-20) log(x) else (x^lambda - 1) / lambda

  Solve <- function(S, halt=TRUE, LINPACK=FALSE, usespam=use_spam) {    
    n <- dim(S)[1]
    if (is.na(usespam)) {
     usespam <-
        (n > spam.min.n &&
         mean(S[(spam.n * spam.factor) %% (n^2)] < spam.tol) >= spam.min.p)
    }
    oldwarn <- options()$warn
    options(warn=no.warning)

    if (usespam) {
      S <- spam::as.spam(S)
      res <- try(spam::chol.spam(S, silent=silent))
      rm("S")
    } else {
      res <- try(chol(S, silent=silent))
    }

    options(warn=oldwarn)
    if (class(res) == "try-error") {      
       #Print(minmax, minmax[,1], model) ## S tandardSimpleModel abaendern !!!!
      #      Print(res, S, eigen(S)$value); xxxx
      if (halt) stop("matrix does not have full rank")
      else {
        assign("LOGDET", NA, envir=ENVIR)
        return(res)
      }
    }

    if (usespam) {
      dg <- as.numeric(spam::determinant(res, logarithm=TRUE))[1]
      res <- spam::solve.spam(res)
    } else {
      dg <- as.numeric(determinant(res, logarithm=TRUE))[1]
      res <- chol2inv(res) #,LINPACK=LINPACK)
    }
    assign("LOGDET", 2 * dg, envir=ENVIR)
    return(res)
  }

  LogDet <- function(S, usespam=use_spam) {
    n <- dim(S)[1]
    if (is.na(usespam)) {
      usespam <-
        (n > spam.min.n &&
         mean(S[(spam.n * spam.factor) %% (n^2)] < spam.tol) >= spam.min.p)
    }
    if (usespam) {
      S <- spam::as.spam(S)
      logdet<-2*as.numeric(spam::determinant(spam::chol.spam(S, silent=silent)))
      rm("S")
    } else {
      logdet <- 2 * sum(log(diag(chol(S, silent=silent))))
    }
    return(logdet)
  }

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
      if(dummy[[1]]!=ZF_SYMBOLS_PLUS) break
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
    if (printlevel>PL_FCTN_DETAILS) Print(LSMIN, format(variab, dig=20))#

    if (n.variab==0) return(NA)		##trivial case

    if (any((variab<LSQLB) | (variab>LSQUB))) {            
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 13.12.03 still happens ...
      if (printlevel>=PL_STRUCTURE)
        WarningMessage(variab, LSQLB, LSQUB, "LSQ")
      assign("BEYOND", BEYOND + 1, envir=ENVIR)
      penalty <- variab
      variab <- pmax(LSQLB, pmin(LSQUB, variab))
      penalty <- + sum(variab-penalty)^2
      LSMIN.save <- LSMIN
      assign("LSMIN", -Inf, envir=ENVIR)
      res <- LStarget(variab)
      assign("LSMIN", LSMIN.save, envir=ENVIR)
      return(res + penalty * (1 + abs(res)))
    }

    param <- as.double(trafo(variab[1:nvarWithoutbc]))
    
    .C("PutValuesAtNA", Reg, param, PACKAGE="RandomFields")
    
    model.values <- double(bins * vdim^2)
    
    .Call("VariogramIntern", Reg, bin.centers, bins, model.values,
       PACKAGE="RandomFields")

#    Print( RFgetModelInfo(reg=Reg, level=2, spC=F),
#          bin.centers, bins, model.values, length(bin.centers),length(model.values))
  
 
    if (any(!is.finite(model.values))) {
      if (printlevel>=PL_IMPORTANT) {
        message("LSQ missing values!")
        Print(format(variab, dig=20), format(param, dig=20)) #
        #print(cbind(bin.centers, as.matrix(model.values, ncol=vdim))) #
      }
      return(1E300)
    } 

    if (LSQ.SELF.WEIGHING) {
      ## weights = 1/ model.values^2
        gx <- binned.variogram / model.values
        gx <- gx[is.finite(gx)]
        if (length(gx) == 0) return(Inf)
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
    base <- 0
    ML.df <- ntotdata
    assign("MLEtarget", MLtarget, envir=ENVIR)

    assign("DO.REML", M =="reml" && anyFixedEffectt, envir=ENVIR)
    assign("DO.RML1", M =="rml1" && anyFixedEffectt, envir=ENVIR)

    if (DO.RML1) {
      ML.df <- ML.df - repet * rowSums(nCoVar[, FixedEffectts])
      rml.a <- rml.x <- rml.data <- vector("list", sets)
      idx <- NULL
      for (i in 1:sets) {
        for (k in mmcomponents)
          if (FixedEffectts[k]) idx <- c(idx, startXges[i, k] : endXges[i, k])
        
        rml.data[[i]] <- rml.a[[i]] <- rml.x[[i]] <-
          vector("list", length(idx.repet[[i]]))

        for (m in 1:len.rep[i]) {
          ortho <- eigen(tcrossprod(Xges[[i]][[m]][, idx, drop=FALSE]))
          ortho <-
            ortho$vectors[ , order(abs(ortho$value))[1:(ML.df[i]/repet[i])],
                          drop=FALSE]
          rml.x[[i]][[m]] <-
            crossprod(ortho, Xges[[i]][[m]][, -idx, drop=FALSE])
          if (!lambdaest) {
            w <- werte[[i]][idx.na[[i]][[m]], idx.na.rep[[i]][[m]], drop =FALSE]
            d <- dim(w)
            dim(w) <- c(d[1] * d[2], d[3])
            rml.data[[i]][[m]] <-  crossprod(ortho, w)                          
          }
          rml.a[[i]][[m]] <- ortho
        }          
      }      
      assign("RML.A", rml.a, envir=ENVIR)
      assign("RML.Xges", rml.x, envir=ENVIR)
      assign("RML.data", rml.data, envir=ENVIR)
      assign("MLEtarget", MLtarget, envir=ENVIR)        
    } else if (DO.REML) {       
      ML.df <-
        ML.df - nCoVar[, FixedEffectts]
      for (k in which(FixedEffectts)) {
        XtX <- 0
        for (i in 1:sets) {
          for (m in 1:len.rep[i]) {
            XtX <- XtX +
              crossprod(Xges[[i]][[m]][, startCoVar[i, k] : endCoVar[i, k]]) 
          }
          XtX <- try(chol(XtX))
          if (!is.numeric(XtX)) stop("X does not have full rank")
          base <- base - 2 * sum(log(diag(XtX)))
        }
      }
    } else if (M=="rml2") {
      stop("not programmed yet")
      ML.df <- ML.df - rowSums(nCoVar[, FixedEffectts]) 
      ## hier wirklich die riesengrosse Matrix erstellen wo wirklich nur
      ## rowSums(nCoVar[, FixedEffectts]) abgezogen wird und nicht
      ## das repet-fache ?!
      ## oder per Rechnung klein machen
      for (i in 1:sets) {
        for (k in mmcomponents) if (FixedEffectts[k])
          idx <- c(idx, startXges[i, k] : endXges[i, k])
        
        Xfix <- NULL
        for (m in 1:len.rep[i]) {
          Xfix <- rbind(Xfix,
                        Xges[[i]][[m]][, idx, drop=FALSE])                  
          Xnonfix <- rbind(Xnonfix, Xges[[i]][[m]][, -idx, drop=FALSE])
        }
        
        ortho <- eigen(tcrossprod(Xfix))
        rml.a[[i]] <-
          ortho$vectors[ , order(abs(ortho$value))[1:ML.df[i]], drop=FALSE]
        rml.x[[i]] <- crossprod(rml.a[[i]], Xnonfix)
        if (!lambdaest) rml.data[[i]] <- crossprod(rml.a[[i]], werte[[i]])
        
      }
    }

    ML.df <- sum(ML.df)
    if (globalsigma)  base <- base + ML.df * (1 - log(ML.df))
    
    dummy <- rowSums(nCoVar[, effect >= RandomEffect & effect <= SpVarEffect,
                            drop=FALSE])
    base <-  base + (sum(dummy) +  ML.df) * log(2 * pi)
    if (solvesigma) {
      dummy <- rowSums(nCoVar[, vareffect, drop=FALSE])
      base <- base + sum(dummy * (1 - log(dummy)))
    }
    
    
    assign("ML.df", ML.df, envir=ENVIR)
     assign("ML.base", base, envir=ENVIR)
  }
 

  linearpart <- function(sigma2, return.z=FALSE) {
    
    z <- list()
    reml.corr <- 0
   
    for (i in 1:sets) {      
      D <- matrix(0, ncol=nCoVarSets[i], nrow=nCoVarSets[i])
      rS <- 0
      for (m in 1:len.rep[i]) {      
        
        Spalte <- crossprod(if (DO.RML1) RML.Xges[[i]][[m]] else Xges[[i]][[m]],
                            Sinv[[i]][[idx.error[1]]][[m]]) # all k
        for (k in mmcomponents) {
          if (nCoVar[i, k] > 0 && (!FixedEffectts[k] || !DO.RML1)) {
            idx <- startXges[i, k] : endXges[i, k]
            
            D[, idx] <- D[, idx] + idx.repet[[i]][m] * Spalte %*%
              if (DO.RML1) RML.Xges[[i]][[m]][, idx, drop=FALSE]
              else Xges[[i]][[m]][, idx, drop=FALSE] 
          }
        }

         
        rS <- rS + Spalte %*% sumdata[[i]][[m]]
      }

      if (DO.REML) {
        idx <- NULL
        for (k in mmcomponents) if (FixedEffectts[k])
          idx <- c(idx, startXges[i, k] : endXges[i, k])
        
        reml.corr <- reml.corr + LogDet(D[idx, idx])      
      }

      j <- 1
      for (k in mmcomponents) {
        if (effect[k] >= RandomEffect && effect[k] <= SpVarEffect) {
          idx <- startXges[i, k] : endXges[i, k]
          if (effect[k] == RvarEffect) {
            D[idx, idx] <- D[idx, idx] + Sinv[[i]][[k]] / sigma2[j]
            j = j + 1
          } else {
            D[idx, idx] <- D[idx, idx] + Sinv[[i]][[k]]
          }
        }
      }

      if (DO.RML1) {
        idx <- NULL
        for (k in mmcomponents) if (FixedEffectts[k])
          idx <- c(idx, startXges[i, k] : endXges[i, k])
        D <- D[-idx, -idx]
      }
      
      z[[i]] <- try(solve(D, rS))
      if (!is.numeric(z[[i]]))
        stop("Design matrix shows linear dependence. Check the design matrix and \n  make sure that you do not use `trend' and the `mixed' model at the same time.")
    }
 
    assign("REML.CORRECTION", reml.corr, envir=ENVIR)
    
    if (return.z) {
      return(z)
    }

    cumsq <- j <- 0
    for (k in mmcomponents) {
      s2est <- 0
      if (effect[k] == RvarEffect) {
        idx <- startCoVar[i, k] : endCoVar[i, k]
        for (i in 1:sets) {
          zi <- z[idx, i]
          s2est <- s2est + crossprod(zi[idx], Sinv[[i]][[k]] %*% zi[idx])
        }
        cumsq <-  cumsq +  (sigma2[j] - s2est / nTotalComp[k])^2
        j <- j + 1
      }
    }

    if (cumsq < LINEARLARP) {
      assign("LINEARLARP", cumsq, envir=ENVIR)
      assign("LINEARLARPZ", z, envir=ENVIR)
      assign("LINEARLARPS2", sigma2, envir=ENVIR)
    }    
    return(cumsq)
  }

  MLtarget <- function(variab) {
 #   Print(variab, nvarWithoutbc)
    
    if (n.variab > 0) {
      variab <- variab + 0  ## unbedingt einfuegen, da bei R Fehler der Referenzierung !! 16.2.10
      
      if (printlevel>=PL_FCTN_DETAILS ) {
        Print(format(variab, dig=20)) #
        if (printlevel>=PL_FCTN_SUBDETAILS) Print(minmax) #
      }
      if (any((variab < MLELB) | (variab > MLEUB))) {
        ## for safety -- should not happen, older versions of the optimiser
        ## did not stick precisely to the given bounds
        ## 23.12.03 : still happens
        if (printlevel>=PL_STRUCTURE)
          WarningMessage(variab, MLELB, MLEUB, "MLE")   
        assign("BEYOND", BEYOND + 1, envir=ENVIR)
        penalty <- variab
        variab <- pmax(MLELB, pmin(MLEUB, variab)) 
        penalty <-  - sum(variab - penalty)^2 ## not the best ....
        MLEMAX.save <- MLEMAX
        assign("MLEMAX", Inf, envir=ENVIR)

#Print(penalty)        
        res <- MLtarget(variab)
        assign("MLEMAX", MLEMAX.save, envir=ENVIR)

        return(res + penalty * (1+ abs(res)))
      }

      param <- as.double(trafo(variab[1:nvarWithoutbc]))
#      if (printlevel>4)  Print(format(param, dig=20))
      .C("PutValuesAtNA", Reg, param, PACKAGE="RandomFields")
      options(show.error.messages = show.error.message)
    } else param <- NULL
    
    if (printlevel >= PL_FCTN_SUBDETAILS) {
      cat("\n\nAufruf von MLtarget\n===================\n")
      Print(RFgetModelInfo(register=Reg, level=14))#
    }

        
    logdet <- ML.base # Konstante; wird mit -0.5 multipliziert
    Werte <- list() ## nur bei BoxCox Schaetzung verwendet

    stopifnot(is.finite(logdet))

#    Print(logdet)
    
    for (i in 1:sets) {

      stopifnot(lc[i] * vdim < 10000)

      S <- double((lc[i] * vdim)^2)
      for (k in allcomponents) {

        if (effect[k] >= SpaceEffect) { ## ansonsten muss schon vorher gesetzt
          ##                              werden

          sp.eff <- effect[k] == SpaceEffect || effect[k] == SpVarEffect
          if (balanced[i] || sp.eff) {
            
            .C("setListElements", Reg, as.integer(i),
               as.integer(if (sp.eff) k else idx.error),
               as.integer(if (sp.eff) 1 else length(idx.error)),
               PACKAGE="RandomFields");

            dim(S) <- rep(lc[i] * vdim, 2)


            .Call("CovMatrixLoc", Reg, Xdistances[[i]], is.dist.vector,
               xdimOZ, as.integer(lc[i]), S, nonzeros,
                  PACKAGE="RandomFields")
            dim(S) <- rep(lc[i] * vdim, 2)
         
          }

#          Print(S, det(S), param); ffff

          ##  Print("xx", S); #    fffff
          #Print(RFgetModelInfo(Reg), effect, k, allcomponents, SpaceEffect)
          
          if (sp.eff) {
            Sinv[[i]][[k]] <<- Solve(S, LINPACK=TRUE)
            logdet <- logdet + LOGDET # side effect of Solve
          } else {  ## RemainingError/PrimitiveModel, genau 1x pro i aufgerufen
            for (m in 1:len.rep[i]) {
              if (balanced[i]) {
                Si <- S[idx.na[[i]][[m]], idx.na[[i]][[m]]]
              } else {
                nSidx <- len.idx.na[[i]][m]
                Si <- double(nSidx^2)

                .C("setListElements", Reg, as.integer(-i),
                   as.integer(if (sp.eff) k else idx.error),
                   as.integer(if (sp.eff) 1 else length(idx.error)),
                   PACKAGE="RandomFields");
               
                .Call("CovMatrixSelectedLoc", Reg, Xdistances[[i]],
                   is.dist.vector, xdimOZ, as.integer(lc[i]),
                   as.integer(idx.na[[i]][[m]] - 1),
                      as.integer(nSidx), Si, nonzeros,
                   PACKAGE="RandomFields")                
                dim(Si) <- rep(nSidx, 2)

                if (!all(eigen(Si)$value > 0))
                  Print("CovMatrixSelectedLoc", Reg, # Xdistances[[i]],
                   is.dist.vector, xdimOZ, as.integer(lc[i]),
                   as.integer(idx.na[[i]][[m]] - 1),
                      as.integer(nSidx), Si, nonzeros,
                   PACKAGE="RandomFields", eigen(Si)$value,
                      GetModel(Reg)) ## OK

               }

              if (DO.RML1) {
                Si <- crossprod(RML.A[[i]][[m]], Si) %*% RML.A[[i]][[m]]
              }

              # Print(Solve, Si)
            
              Sinv[[i]][[k]][[m]] <<- Solve(Si,  halt=FALSE, LINPACK=TRUE)
              
              if (class(Sinv[[i]][[k]][[m]]) == "try-error") {
                assign("MLEINF", TRUE, envir=ENVIR)
                if (printlevel>=PL_SUBIMPORTANT || is.na(MLEVARIAB)) {
                  cat("\nMLE: error in cholesky decomp. -- matrix pos def?")
                  if (printlevel>=PL_ERRORS) {
                    Print(variab, MLEVARIAB, S, Si, eigen(Si)$val)#
                  }
                }
                return(1E300)
              }
              
              logdet <- logdet + LOGDET * idx.repet[[i]][m]
              
              ## Si ist hier Si^{1/2}!! Zum Schluss wird mit -0.5 multipliziert
              
            }
          }
        }
      } # k, components
      S <- NULL
      
      Werte[[i]] <- list()
      for (m in 1:len.rep[i]) {
        if (lambdaest) {
          
          Werte[[i]][[m]] <-
            BoxCox(werte[[i]] [ idx.na[[i]][[m]], idx.na.rep[[i]][[m]],
                               drop = FALSE ], variab[n.variab])
          if (DO.RML1) 
            Werte[[i]][[m]] <- crossprod(RML.A[[i]][[m]],  Werte[[i]][[m]])
          for (k in which.deteff) {
            Werte[[i]][[m]] <- Werte[[i]][[m]] -
              Xges[[i]][[m]][, startXges[i, k] : endXges[i, k] ] 
          }
        } else {          
          Werte[[i]][[m]] <-
            if (DO.RML1) RML.data[[i]][[m]]
            else werte[[i]] [ idx.na[[i]][[m]], idx.na.rep[[i]][[m]],
                             drop = FALSE ]
        }
        
        sumdata[[i]][[m]] <<- as.vector(t(apply(Werte[[i]][[m]], 1, sum)))
        ##d <- dim(Werte[[i]][[m]])
                
        ##dim(Werte[[i]][[m]]) <- c(d[1] * d[2], d[3])
      } # for m
    } # for i (sets)

#    Print(Werte)


    assign("REML.CORRECTION", 0, envir=ENVIR)


    if (nCoVarAll > 0) {
      assign("LINEARLARP", Inf, envir=ENVIR)
      if (!solvesigma) {
        ## !solvesigma: alle Parameter der Kovarianzen auf
        ## globaler Ebene geschaetzt       
        linearpart(LINEARLARPS2)     
      } else  {
        OPTIM(LINEARLARPS2, linearpart,
              lower = LINEARLARPS2 / 10, upper = LINEARLARPS2 * 10,
              control= lsq.optim.control, optimiser=optimiser,silent=silent)
      }
    }
    
    quadratic <- 0
    quad.effect <- rep(0, max(1, length(mmcomponents)))
    MLtargetV <- list()
    for (i in 1:sets) {
      MLtargetV[[i]] <- list()
      for (m in 1:len.rep[i]) {       
        MLtargetV[[i]][[m]] <- Werte[[i]][[m]] 
        if (nCoVarAll>0)
          MLtargetV[[i]][[m]] <- MLtargetV[[i]][[m]] -
            as.vector(if (DO.RML1) RML.Xges[[i]][[m]] else Xges[[i]][[m]]
                      %*% LINEARLARPZ[[i]])
                                        # y-\sum A_i z_i
                
        quadratic <- quadratic +
          sum(MLtargetV[[i]][[m]] *
              (Sinv[[i]][[idx.error[1]]][[m]] %*% MLtargetV[[i]][[m]]))       
      }
      
      for (k in mmcomponents) if (effect[k] >= RandomEffect) {
        idx <- startCoVar[i, k] : endCoVar[i, k]
        quad.effect[k] <- quad.effect[k] +
          sum(LINEARLARPZ[[i]][idx] *
              (Sinv[[i]][[k]] %*% LINEARLARPZ[[i]][idx]))
      }
    }

    if (solvesigma) {
      nc <- rowSums(nCoVar[, vareffect])
      quad.effect[vareffect] <- log(quad.effect[vareffect] / nc) * nc
    }

    res <- -0.5 * (logdet + sum(quad.effect) + REML.CORRECTION +
                    if (globalsigma) ML.df * log(quadratic) else quadratic)
   
    if (is.na(res)) {
      Print(logdet, sum(quad.effect), REML.CORRECTION, globalsigma, ML.df,#
            quadratic,
            if (globalsigma) ML.df * log(quadratic) else quadratic)
      warning("NA results appeared")
    }
    
    if (printlevel >= PL_FCTN_SUBDETAILS) {
      Print(globalsigma, variab, param, var.global, ML.df, #
            quadratic / ML.df)
      readline(paste(res, "Bitte return druecken."))
    }

    if (res > MLEMAX) {
      if (varnugNA) { ## then globalsigma is true as well. So never change oder of
        ##               if's
        sill <- quadratic / ML.df ### check!!!
        param[var.global] <- sill * (1.0 - param[nugget.idx])        
        param[nugget.idx] <- sill * param[nugget.idx]
      } else {
        if (globalsigma)
          param[var.global] <- quadratic / ML.df  ### check!!!
      }

      assign("MLEMAX", res, envir=ENVIR)
      assign("ML.RESIDUALS", MLtargetV, envir=ENVIR)
      assign("MLEPARAM", param, envir=ENVIR)
      assign("MLEVARIAB", variab, envir=ENVIR)      
  
    }

    if (printlevel >= PL_FCTN_SUBDETAILS) Print(logdet, quadratic, res)#

    return(res)
  } # mltarget


  get.covariates <- function(Variab) {
    if (nCoVarAll == 0) return(NULL)
    max <- MLEMAX
    param <- MLEPARAM
    inf <- MLEINF
    variab <- MLEVARIAB 
    resid <- ML.RESIDUALS
    linpart <- LINEARLARP
    linz <- LINEARLARPZ
    lins2 <- LINEARLARPS2
    
    assign("MLEMAX", -Inf, envir=ENVIR)
    assign("LINEARLARP", Inf, envir=ENVIR)

    MLEsettings("ml")
    MLtarget(Variab)
         
    result <- unlist(LINEARLARPZ) 
    assign("MLEMAX", max, envir=ENVIR)
    assign("ML.RESIDUALS", resid, envir=ENVIR)
    assign("MLEPARAM", param, envir=ENVIR)
    assign("MLEVARIAB", variab, envir=ENVIR)
    assign("MLEINF", inf, envir=ENVIR)
    assign("LINEARLARP", linpart, envir=ENVIR)
    assign("LINEARLARPZ", linz, envir=ENVIR)
    assign("LINEARLARPS2", lins2, envir=ENVIR)   
    return(result)
  }

  


  ## to avoid warning on "no visible binding" we define the following
  ## variable that are used in the local functions:
  ENVIR <- environment()
  LSQ.SELF.WEIGHING <- LSQ.WEIGHTS <- LSQ.BINNEDSQUARE <- 
    DO.REML <- DO.REML1 <- RML.A <- RML.Xges <- RML.data <-
      REML.CORRECTION <- DO.RML1 <- 
        ML.base <- ML.df <- MLEtarget <-
          ML.RESIDUALS <- MLEMAX <- MLEINF <- MLEPARAM <- 
            CROSS.DIST <- CROSS.KRIGE <- CROSS.VAR <- CROSSMODEL <-
              LINEARLARP <- LINEARLARPS2 <- LINEARLARPZ  <-
                LOGDET <- NULL
  BEYOND <- 0
 
  
######################################################################
###              End of definitions of local functions             ###
######################################################################


######################################################################
###    Initial settings, coord part I (without model info)         ###
######################################################################

  time <- !is.null(T) && (!is.list(T) || !is.null(T[[1]]))

  if (printlevel>=PL_STRUCTURE) cat("\ninitial settings...\n")

  variab.units <- RFopt$coords$varunits


  #Print(x, y, z); pp

  Z <- StandardizeData(x=x, y=y, z=z, T=T, grid=grid, data=data,
                       distances=distances, dimensions=dimensions, RFopt=RFopt)
  coord <- Z$coord
  data <- Z$data
  dimdata <- Z$dimdata
  isRFsp <- Z$isRFsp
  dist.given <- Z$dist.given
  gridTopology <- Z$gridTopology
  data.RFparams <- Z$data.RFparams
  lc <- Z$lc
  spdim <- Z$spdim
  tsdim <- Z$tsdim
  rangex <- Z$rangex
  coord.units <- Z$coord.units
  matrix.indep.of.x.assumed <- Z$matrix.indep.of.x.assumed
  neu <- Z$neu
  RFsp.coord <- Z$RFsp.coord
  distances <- Z$distances
  

 # Print(data, neu, data.RFparams, sdvar)

  
#  Print(Z); kkk

  
  small.dataset <- max(lc) <= fit$smalldataset    
    
  ## check box-cox transformation 
  if (lambdaest <- !missing(bc_lambda) && is.na(bc_lambda)) {
    if (any(unlist(data) <= 0))   warning("negative data may cause errors")
    if (fit$bc_lambda_lb > 1 || fit$bc_lambda_ub<1)
      stop("bounds for Box-Cox lambda must satisfy lambda_lb <= 1 <= lambda_ub")
  }


################    analyses of orginal model        ###############
##### variables needed for analysis of trend, upper and lower input --
##### user cannot know what the internal represenatation is

  if (printlevel>=PL_STRUCTURE) cat("\nfirst analysis of model  ...\n")

  variable.names <- extractVarNames(model)
  ## gets response part of model, if model is a formula

  ##
  xdimOZ <- as.integer(tsdim - !is.null(T[[1]]))

  info.cov <- .Call("SetAndGetModelInfo", Reg, list("Cov", model),
                    spdim, FALSE, TRUE, time, xdimOZ,
                    as.integer(if (printlevel <= 4) fit$short
                               else min(4, fit$short)),
                    FALSE, TRUE,
                    PACKAGE="RandomFields")


  ## hier zum ersten mal model verwendet wichtig,
  ## da interne Darstellung abweichen kann. Z.B. dass ein optionaler Parameter
  ## auf einen Standardwert gesetzt wird
    model <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE, spConform=FALSE,
                    do.notreturnparam=TRUE) ##standardisiert
  #userdefined <- GetParameterModelUser(model)

  
  modelinfo <- RFgetModelInfo(register=Reg, level=2, spConform=FALSE)
  vdim <- modelinfo$vdim
 ## -- wichtig fuer GetValuesAtNA
    
  NAs <-  info.cov$NAs
  trans.inv <- info.cov$trans.inv ## note: only with respect to the
  ##              coordinates, mixed effect koennen andere wirkung haben
  isotropic <- info.cov$isotropic
  
  if (dist.given && (!trans.inv || !isotropic)) {
    #Print(dist.given, !trans.inv, !isotropic, recall)
    stop("only domain and isotropic models go along with distances")
  }
  minmax <- info.cov$minmax # 4 Spalten: 1:min, 2:max, 3:type, 4:is.nan, not na

  diag.idx <- which(minmax[, 3] == DIAGPARAM)
  if (length(diag.idx)>0) minmax[diag.idx[1], 1] <- fit$min_diag

  ncovparam <- nrow(minmax)
  ptype <- minmax[, 3]

  if (printlevel >= PL_RECURSIVE) print(minmax) #

  
  ## identify the effects:
  effect <- info.cov$effect
  idx.error <- effect >= PrimitiveModel
  mmcomponents <- which(!idx.error & effect > DeterministicEffect)
  
  idx.error <- which(idx.error)
  allcomponents <- c(mmcomponents, idx.error[1])
  if (length(idx.error) == 0)
    stop("there must be an error component in the model")
  if (matrix.indep.of.x.assumed && !info.cov$matrix.indep.of.x)
    stop("x-coordinates are neither given by 'x' nor by 'distances' nor by 'data',\n  but the model seem to require them")

  stopifnot(sum(NAs) == nrow(minmax))
   
  eff <- rep(FALSE, nrow(minmax))  ## lokale Hilfsvariable
  csNAs <- cumsum(c(0, NAs))
  for (k in 1:length(effect)) {
    if (effect[k] >= RandomEffect && effect[k] <= SpVarEffect && NAs[k] > 0)
      eff[(csNAs[k]+1) : csNAs[k]] <- TRUE
  }

  
  
  minmax[ptype == MIXEDVAR & eff, 1] <- fit$minmixedvar
  minmax[ptype == MIXEDVAR & eff, 2] <- fit$maxmixedvar
  rm("eff")


  FixedEffectts <- effect == FixedEffect | effect == FixedTrendEffect
  anymixedeffects <- any(effect >= FixedTrendEffect & effect <= SpVarEffect)
  anyRandomEffects <- any(effect >= RandomEffect & effect <= SpVarEffect)
 # variogram.known <- trans.inv && !info.cov$matrix.indep.of.x && !anymixedeffects
  variogram.known <- trans.inv && !anyRandomEffects

#  Print(variogram.known, trans.inv, !anyRandomEffects)

  DeterministicEffecs <- effect == DeterministicEffect | effect == DetTrendEffect
  trends <- effect == DetTrendEffect | effect == FixedTrendEffect
  which.deteff <- which(DeterministicEffecs)
  anyFixedEffectt <- any(FixedEffectts)
  vareffect <- which(effect == RvarEffect)
  if (length(vareffect) == 0) solvesigma <- FALSE
 
  is.dist.vector <- trans.inv && small.dataset || dist.given
  coord_system <- RFoptions()$coords$coord_system
  cartesian <- coord_system %in% c("auto", CARTESIAN_SYSTEMS)


  if (is.scalar.dist <- isotropic && is.dist.vector && cartesian) {
    xdimOZ <- as.integer(1)
  }
  ts.xdim <- xdimOZ + !is.null(T[[1]])

  
  #Print(is.dist.vector, trans.inv, small.dataset, dist.given)
######################################################################
## model specific settings, upper & lower
######################################################################
   if (printlevel>=PL_STRUCTURE) cat("\nupper & lower...\n")

#  xxxxxxxxxxxxxx
  
  trafoidx <- trafo <- NULL
  if (!is.null(transform)) {
    if (!recall) message("Note: if 'transform' is used, 'anisoT' within 'RMS' should be used instead of 'Aniso', if elemens of an anisotropy matrix are estimated")
    if (trafo <- is.list(transform) && length(transform)>0) {
      trafoidx <- transform[[1]]
      if (length(transform)==2 && nrow(minmax) == length(trafoidx) &&
          is.logical(trafoidx)) {        
        trafo <- transform[[2]]
      }
    }
  }

  if (!is.null(trafo) && !is.function(trafo)) {
    if (length(trafoidx) != nrow(minmax) || !is.numeric(trafoidx))
    .C("PutValuesAtNA", Reg, as.double((1:nrow(minmax)) * 1.11),
       PACKAGE="RandomFields", NAOK=TRUE) # ok
    model_with_NAs_replaced_by_ranks <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE)
    cat("transform must be a list of two elements:\n* A vector V of logicals whose length equals the number of NAs in the model.\n* A function taking a vector whose length equals sum(V). The output\n  is the original model where the NAs are replaced by real numbers.\n\nThe currently identified NAs are:\n")
    print(minmax[, c(1:2)]) #
    cat("\nwithin the following model (NA positions are given by n.nn) :\n")
    str(model_with_NAs_replaced_by_ranks, vec.len=20) #
    cat("\nHowever, 'RFfit' was called by\ntransform =")
    str(transform) #
    cat("\n")
    if (length(trafoidx) > nrow(minmax))
      cat("Note that the parameters of the trend are optimised analytically, hence they may not be considered, here.\n")
    return(NULL)
  }



#######################   upper,lower,user     ########################
  if (printlevel>=PL_STRUCTURE) cat("\nlower and upper ...\n")
  
  users.lower <- users.upper <- NULL
  if (!is.null(lower)) {
    users.lower <- try(GetValuesAtNA(NAmodel=model, valuemodel=lower,
                                     spdim=spdim,
                                     Time=time, shortnamelength=fit$short,
                                     skipchecks=!is.null(trafo)) )
    if (!is.numeric(users.lower)) {
      if (printlevel>=PL_IMPORTANT) Print(model, lower, users.lower)
      stop("'lower' does not match 'model'")
    }
  }
  if (!is.null(upper)) {
    users.upper <- try(GetValuesAtNA(NAmodel=model, valuemodel=upper,
                                     spdim=spdim,
                                     Time=time, shortnamelength=fit$short,
                                     skipchecks=!is.null(trafo)))
   if (!is.numeric(users.upper)) {
      if (printlevel>=PL_IMPORTANT) Print(model, upper, users.upper)
      stop("'upper' does not match 'model'")
    }
  }
 

#  Print(users.guess)
  
  if(!is.null(users.guess)) {
    Users.guess <- try(GetValuesAtNA(NAmodel=model, valuemodel=users.guess,
                                     spdim=spdim,
                                     Time=time, shortnamelength=fit$short,
                                     skipchecks=!is.null(trafo)))
    
    if (!is.numeric(Users.guess)) {
      if (printlevel>=PL_IMPORTANT) Print(model, users.guess, Users.guess)
      stop("'users.guess' does not match 'model'")
    }
    if (any(is.finite(Users.guess))) users.guess <- Users.guess
    else users.guess <- NULL
  }
 
  if (anyRandomEffects) {    
    stopifnot(model[[1]] %in% ZF_PLUSSELECT)
    model[[1]] <- ZF_SELECT[1]
  }

  if (is.dist.vector || anyRandomEffects) {
    info.cov <- .Call("SetAndGetModelInfo", Reg, list("Cov", model),
                      spdim, is.dist.vector, !is.dist.vector,
                      time, xdimOZ,
                      fit$short, FALSE, TRUE,
                      PACKAGE="RandomFields")

    modelinfo <- RFgetModelInfo(register=Reg, level=2, spConform=FALSE)
     
    model <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE, spConform=FALSE,
                      do.notreturnparam=TRUE)
  }


###########################      transform     #######################
## either given bu users.transform + users.min, users.max
## DIESER TEIL MUSS IMMER HINTER GetValuesAtNA STEHEN
  
  if (printlevel>=PL_STRUCTURE) cat("\ntransform ...\n")
  
  lower <- minmax[,1] 
  upper <- minmax[,2]
  delete.idx <- rep(FALSE, length(lower))
  if (is.null(trafo)) {
    if (any(minmax[,4]==1)) ## nan, not na
      stop("NaN only allowed if transform is given.")
    trafo <- function(x) x;
  } else { ## is function
    origNAs <- nrow(minmax)
    delete.idx <- !trafoidx
    if (any(minmax[,4]==1)) {
      if (any(minmax[delete.idx,4] != 1) || any(minmax[delete.idx,4] == 1)) 
        stop("NaNs do not match logical vector of transform")
    }
    optim_var_elimination <- "never"
    solvesigma <- FALSE
    try(z <- trafo(lower[!delete.idx]))
    if (!is.numeric(z) || !all(is.finite(z)) || !is.vector(z))
      stop("The transformation does not return a vector of finite numerical values.")
    if (length(z) != length(lower))
      stop("\n'transform' returns a vector of length ",
           length(z), ", but one of length ", origNAs, " is expected.",
           if (length(z) > origNAs) " Note that the parameters of the trend are optimised analytically, hence they may not be considered, here.",
           " Call 'RFfit' with `transform=list()' to get more information on the parameters of the model.\n\n")
  }
  

  ## Achtung which, da upper,lower etc um box-cox-Variable verlaengert
  ## werden koennten !
  SDVAR.IDX <- ptype == SDPARAM | ptype == VARPARAM | ptype == NUGGETVAR
  SIGN.VAR.IDX <- ptype == SIGNEDVARPARAM
  SIGN.SD.IDX <- ptype == SIGNEDSDPARAM
  ALL.SDVAR <- SDVAR.IDX | SIGN.VAR.IDX | SIGN.SD.IDX

  if (is.null(sdvar)) sdvar <- matrix(SDVAR.IDX, nrow=ncovparam, ncol=vdim)
  else stopifnot(all(rowSums(sdvar[SDVAR.IDX, ]) >= 1))
  SCALE.IDX <- ptype == SCALEPARAM  ## large capitals 
  var.idx <- which(ptype == VARPARAM)
  sd.idx <- which(ptype == SDPARAM)
  nugget.idx <- which(ptype == NUGGETVAR)
  MIXED.IDX <- which((ptype == MIXEDVAR) & (is.na(solvesigma) || solvesigma))
  mixed.idx <-  which(ptype == MIXEDVAR)
  solvesigma <- length(MIXED.IDX) > 0


#################################################################
##############     prepare constants in S, X,etc      ###########
#################################################################

##############         distances              #################

  ## note: the direct C call needs matrix where points are given column-wise
  ##       whereas the R function CovarianceFct need them row-wise,
  ##                   except for fctcall==CovarianceMatrix
  ## distances are used to calculate bounds for the scale parameter(s)
  ## to be eliminationated -- further it is used in MLtarget

###  30.3.06 folgenden code durch nachfolgenden ersetzt.
###     dd <- 0   
###     for (i in 1:tsdim) {
###      dx <- outer(coord[,i], coord[,i], "-")
###      distances[i, ] <- as.double(dx[lower.tri(dx)])
###      ## x <- c(1,6,0.1)
###      ## outer(x, x, "-")
###      ##     [,1] [,2] [,3]
###      ##[1,]  0.0 -5.0  0.9
###      ##[2,]  5.0  0.0  5.9
###      ##[3,] -0.9 -5.9  0.0
###      if (i<=spdim) dd <- dd + distances[i, ]^2 
###    }
###    dd <- sqrt(dd)
###    coord <- NULL
###    maxdistances <- max(dd)
###    mindistances <- min(dd[dd!=0])  ## notwendig falls !is.null(T)
###    rm("dd")
###
 if (printlevel>=PL_STRUCTURE) cat("\ndistances ...\n")
  sets <- length(lc)
  Xdistances <- vector("list", sets)

  n.sample.smallset <- fit$smalldataset / 2
  
  mindistances <- matrix(ncol=sets, nrow=2)

  if (is.scalar.dist) {
    for (i in 1:sets) {
     
      Xdistances[[i]] <-
        if (dist.given) distances[[i]] else as.vector(dist(t(coord[[i]]))) 
      idx <- Xdistances[[i]] < nuggetrange
      if (any(idx)) {
        if (!general$allowdistanceZero)
          stop("distance with value 0 identified -- use allowdistanceZero=T?")
        Xdistances[[i]][idx] <- nuggetrange
      }
      
      mindistances[, i] <- range(Xdistances[[i]][!idx])
    }
  } else if (dist.given) {
    stop("Only distances, not distance vectors, and corresponding models allowed")
  } else {
    for (i in 1:sets) {
      if (small.dataset && lc[i] <= n.sample.smallset) {
        size <- lc[i]
        coord.idx <- TRUE
      } else {
        size <- n.sample.smallset
        coord.idx <- sample.int(lc[i], size=size)
      }

      d <- .Call("vectordist", coord[[i]][, coord.idx, drop=FALSE], FALSE,
                 PACKAGE="RandomFields")
      
      absdist <-
        sqrt(colSums(if (is.null(T[[i]])) d^2 else 
                       d[1:spdim, d[ts.xdim, ]==0, drop=FALSE]^2))
      idx <- absdist < nuggetrange

      if (any(idx)) {
        if (!general$allowdistanceZero) 
          stop("distance with value 0 identified - use allowdistanceZero=TRUE?")
        coord[[i]] <- coord[[i]] *
          (1 + rnorm(length(coord[[i]]), 0, nuggetrange))
      }
      
      mindistances[, i] <- range(absdist[!idx])        
  
      if (printlevel>=PL_FCTN_SUBDETAILS) Print(trans.inv, isotropic)#
      if (is.dist.vector) {
        ##stopifnot(!isotropic)
        ## if (isotropic) { Xdistances[[i]] <- absdist } else
        Xdistances[[i]] <- d
      } else {
        Xdistances[[i]] <- coord[[i]]
      }
      storage.mode(Xdistances[[i]]) <- "double"
    }
  }
  
  maxdistances <- max(mindistances[2, ])
  mindistances <- min(mindistances[1, ])
  ## if small.dataset these values are not the true values. Can/Should
  ## they be improved??
 
  
##############         Coordinates & data    #################
  if (printlevel>=PL_STRUCTURE) cat("\ndistances and data...")
    ## note: the direct C call needs matrix where points are given column-wise
    ##       whereas the R function CovarianceFct need them row-wise,
    ##                   except for fctcall==CovarianceMatrix

  if (lambdaest && vdim > 1) {
    stop("lambda can be eliminationated only in the univariate case")
  }    

  ntotdata <- sapply(data, length)
  idx.na <- len.idx.na <- idx.na.rep <-  werte <- sumdata <-
    vector("list", length(coord))
  repet <- len.rep <- numeric(length(coord))
  idx.repet <-  vector("list", length(sets))

  balanced <- rep(TRUE, sets)


    ##   Print(ntotdata); ddddd

  #Print(isRFsp, dimdata, sets, data)


  sum.not.isna.data <- 0
  if (vdim > 1) {
    varlist <- list()
    for (k in 1:vdim) varlist[[k]] <- list()
  }
  for (i in 1:sets) {
    repet[i] <- length(data[[i]]) /  (lc[i] * vdim) ## repetitions

    # Print(repet, data, lc, vdim, i, sets)
    
    if (repet[i] != as.integer(repet[i]))
      stop("number of data does not match number of coordinates" )

    dim(data[[i]]) <- c(lc[i], vdim, repet[i])
    
    werte[[i]] <- array(data[[i]], c(lc[i] * vdim, repet[i])) # matrix only

#    stopifnot(all(abs(werte[[i]]) < 10))
    
    for (k in which.deteff) {
      ## for the deterministic effects:
      ## DeterministicEffecs <- effect == DeterministicEffect | effect == DetTrendEffect
      ## trends <- effect == DetTrendEffect | effect == FixedTrendEffect
      if (!trends[k]) { # i.e. == DeterministicEffect
        size <- nrow(modelinfo$sub[[k]]$param$X[[i]])
        if (size != 1) {  
          if (size != nrow(werte[[i]]))
            stop("matrix of data does not match deterministic effect")
        }
      }
    }


    if (!lambdaest) {
      if (!missing(bc_lambda)) werte[[i]] <- BoxCox(werte[[i]], bc_lambda)
      if (length(which.deteff) > 0) {
        deteffs <-
          rfSplitTrendAndCov(model=model, spatialdim=spdim,
                             xdimOZ = xdimOZ, Time = FALSE,
                             cathegories=list(det=c(DeterministicEffect,
                                                DetTrendEffect)))
        
        if (!FALSE)
        stop(FALSE)
## Achtung! Muss hier x,y,z,T aus coord basteln
#             
#         deteffs <- RFsimulate(model=deteffs$det, grid=FALSE,
#                               x=if(xlen>1) t(coord[[i]]) else coord[[i]][1,], 
#                               y=if(is.null(y)) NULL else coord[[i]][2,],
#                               z=if(is.null(z)) NULL else coord[[i]][3,],
#                               T=if(is.null(T[[i]])) NULL 
#                                 else coord[[i]][nrow(coord[[i]]),], 
#                               register = TrendReg, spConform=FALSE)
         
        werte[[i]] <- werte[[i]] - deteffs         
      
  #      for (k in which.deteff) {
  #        if (trends[k]) {
  #          if (length(modelinfo$sub[[k]]$mean) > 0)
  #            werte[[i]] <- werte[[i]] - modelinfo$sub[[k]]$mean
  #          if (length(modelinfo$sub[[k]]$plane) > 0)
  #            werte[[i]] <-
  #              werte[[i]] - as.vector(modelinfo$sub[[k]]$plane %*% coord[[i]])
  #          if (length(modelinfo$sub[[k]]$polydeg))
  #            stop("Poly not programmed yet")
  #          if (length(modelinfo$sub[[k]]$arbitraryfct))
  #            stop("Arbitary fct not programmed yet")
  #        } else {
  #          werte[[i]] <- werte[[i]] -
  #            as.vector(modelinfo$sub[[k]]$param$X[[i]])
  #        }
  #      }
      }
    }


    not.isna <- !is.na(werte[[i]])
    idx.na[[i]] <- idx.na.rep[[i]] <- sumdata[[i]] <- list()
    sum.not.isna.data <- sum.not.isna.data + sum(not.isna)
     
    if (general$na_rm_lines) {
      dim(not.isna) <- c(lc[i], vdim, repet[i])
      not.isna <- matrix(ncol=repet[i], byrow=TRUE,
                         rep(t(apply(not.isna, c(1, 3), all)), times=vdim))
    }

    
    if (any(colSums(not.isna) > fit$max_neighbours)) {

      stop("too many points -- currently composite likelihood is not programmed") # to do
      
      balanced[i] <- FALSE
    
      if (printlevel>=PL_IMPORTANT && fit$split)
        message("Too many locations to use standard eliminationation methods.\n",
            "Using approximative methods. However, it is *not* ensured\n",
            "that they will do well in detecting long memory dependencies.")
      if (is.scalar.dist) { ## distances
        j <- 1:lc
        composite <- list()
        li <- 0
        maxi <-  (fit$splitn_neighbours[2] - 1) / vdim 
        # mini <-  fit$splitn_neighbours[3] / vdim
        while (length(j) >= maxi) {
          distpos <- (j[1]-1) * (lc - j[1] / 2) + 1
          distidx <- distpos : (distpos + (lc-j[1]-1))
          
          locs <- (j[1]+1):lc
          locidx <- which(!is.na(pmatch(locs, j)))
          ##
          if (length(locidx) > maxi) {
            locidx <- locidx[order(Xdistances[[i]][distidx[locidx]])[1:maxi]]
          }
          locidx <- c(locs[locidx], j[1])    

          li <- li + 1
          composites[[li]] <- locidx
           j <- j[is.na(pmatch(j, locidx))]
        }
        if (length(j) > 0) {
          li <- li + 1
          composite[[li]] <- j
        }
        
        ## kann jetzt noch optimiert werden hinsichtlich schaetzqualitaet durch
        ## zusammenfassen
      } else if (is.dist.vector) {
        stop("distance vector not fully programmed yet!")
      } else {
        modelnr <- -1
        if (!is.null(users.guess)) {
          .C("PutValuesAtNA", Reg, users.guess,PACKAGE="RandomFields")
          modelnr <- Reg         
        }
        if (ntotdata <= fit$max_neighbours)
          composite <- list(1:ncol(Xdistances[[i]]))
        else 
          composite <- GetNeighbourhoods(modelnr,
                                 list(given=Xdistances[[i]], vdim=vdim),
                                 fit$splitfactor_neighbours,
                                 fit$max_neighbours,
                                 fit$splitn_neighbours,
                                 shared=FALSE)
      }
      old.not.na <- not.isna

#      Print(l)
      
      not.isna <- matrix(FALSE, nrow=lc[i] * vdim, repet[i] * length(composite))
      save.rep <- numeric(repet * length(composite))
      for (j in 1:length(composite)) {
        idx <- rep(FALSE, lc[i])
        idx[composite[[j]]] <- TRUE
        not.isna[idx, ( (j-1) * repet + 1 ) : (j * repet)] <-
          old.not.na[idx, ]
        save.rep[ ((j-1)*repet + 1) : (j * repet)] <- repet
      }
#      Print(save.rep)
    } else {
      # not: any(colSums(not.isna) > fit$max_neighbours)
      save.rep <- 1:repet[i]
    }
    
    ## idx.na[[i]] indiziert ob fuer feste coordinate & feste
    ## wiederholung irgendeine multivariate komponente NA ist
    ## idx.na.rep[[i]] holt sich, wo die wiederholungen sind
    repetitions <- list()
    j <- 1:ncol(not.isna)
    idx.na[[i]] <- list() # length(idx.na[[i]][[m]]) = coordinates * vdim 

    stopifnot(any(not.isna))
    
    
    while (length(j) > 0) {
      current <- not.isna[, j[1]]
      
      if (sum(current) < 2) {
        #
        stop("data seem to be rather clumsy.")
        j <- j[-1]
        next
      }
      rep.curr <- which(apply(not.isna==current, 2, all))    
      j <- j[is.na(pmatch(j, rep.curr))]
 
      idx.na[[i]][[length(idx.na[[i]]) + 1]] <- as.integer(which(current))     
      repetitions[[length(repetitions) + 1]] <- save.rep[rep.curr]

    }
    
    len.rep[i] <- length(repetitions)
    len.idx.na[[i]] <- sapply(idx.na[[i]], length)    
    idx.na.rep[[i]] <- repetitions   
    idx.repet[[i]] <- sapply(idx.na.rep[[i]], length)
    
    for (m in 1:len.rep[i]) {
      sumdata[[i]][[m]] <-
        as.vector(t(apply(werte[[i]][idx.na[[i]][[m]],
                                     idx.na.rep[[i]][[m]],
                                     drop =FALSE], 1, sum)))
    }

    if (vdim > 1) {
      dim(werte[[i]]) <- c(lc[i], vdim, repet[i])
      for (k in 1:vdim) varlist[[k]][[i]] <- werte[[i]][, k,]
      m <- apply(werte[[i]], 2, mean)
      s <- sqrt(apply(werte[[i]], 2, var))
      if (printlevel > 0 && !recall) {
        if (mean(abs(m)) != 0 && any(log(sd(m) / mean(abs(m))) > 1.5))
          message("Are the average values of the components rather different? If so, it might be\n worth thinking of standardising the values before calling RFfit.\n")
        else if (any(abs(log(s / s[1])) > 1.5))
          message("The standard deviations of the components are rather different. It might be\n better to standardise the components of the data before calling RFfit.\n")
      }
      dim(werte[[i]]) <- c(lc[i] * vdim, repet[i])
    }
  } # i in sets

  if (vdim > 1) {
    vardata <- sapply(varlist, function(x) var(unlist(x), na.rm=TRUE))
    remove("varlist")
  }  else {
    vardata <- var(unlist(werte), na.rm=TRUE)
  }
  varmax <- apply(sdvar, 1, function(x) if (any(x)) max(vardata[x]) else NA)
  varmin <- apply(sdvar, 1, function(x) if (any(x)) min(vardata[x]) else NA)
#  Print(varmax, varmin, sdvar, vardata)
  
 

  if (vdim>1 && printlevel>=PL_IMPORTANT && !recall)
    message("Due to the covariance model a ", vdim,
        "-variate random field is expected. Therefore, \nthe data matrix",
        " is assumed to consist of ", repet,
        " independent measurements for\neach point.",
        " Each realisation is given as the entries of ", vdim,
        " consecutive \ncolumns.")

 
##############      find upper and lower bounds      #################
  if (printlevel>=PL_STRUCTURE) cat("\nbounds...")
 
  txt <- "lower and upper are both lists or vectors of the same length or NULL"
  lengthmismatch <- "lengths of bound vectors do not match model"
  structuremismatch <- "structures of bounds do not match the model"

 
  varnames <- minmax.names <- attr(minmax, "dimnames")[[1]]
   
  ## autostart will give the starting values for LSQ
    ## appears if trafo is given. Then better do not search for
    ## automatic bounds

  autostart <- numeric(length(lower))
  neg <- lower <= 0
  autostart[neg] <-  0.5 * (lower[neg] + upper[neg])
  autostart[!neg] <- sqrt(lower[!neg]*upper[!neg])

  idx <- which(SDVAR.IDX)
#  Print(idx)
                

  if (length(idx) > 0) {
    ## lower bound of first model is treated differently!
    ## so the "main" model should be given first!             !!!!!
     
    
    ## lower[idx] <- 0
    ## first.idx <- nugget.idx
    ## if (is.null(first.idx)) first.idx <- var.idx
    ## if (is.null(first.idx)) first.idx <- sd.idx
     lower[idx] <- varmin[idx] / fit$lowerbound_var_factor / length(idx)
   #  lower[idx] <- 0 ## ??
   
    if (fit$lowerbound_var_factor == Inf && length(idx)>1) {
      idx2 <- idx[if (is.null(users.guess)) length(idx) else
                 which(users.guess == max(users.guess, na.rm=TRUE))[1] ]
      lower[idx2] <- varmin[idx] / 1e8
    }
    
    upper[idx] <- varmax[idx] * fit$upperbound_var_factor
    autostart[idx] <- ((length(idx)-1) * lower[idx] + upper[idx]) /length(idx)
    
    if (length(sd.idx) > 0) {
      lower[sd.idx] <- sqrt(lower[sd.idx])
      upper[sd.idx] <- sqrt(upper[sd.idx])
      autostart[sd.idx] <- sqrt(autostart[sd.idx])
    }       
  }
  
  
  if (any(SIGN.VAR.IDX)) {
    lower[SIGN.VAR.IDX] <-
      -(upper[SIGN.VAR.IDX]<- varmax[SIGN.VAR.IDX] * fit$upperbound_var_factor);
    autostart[SIGN.VAR.IDX] <- 0 
  }
  
  if (any(SIGN.SD.IDX)) {
    lower[SIGN.SD.IDX] <-
      -(upper[SIGN.SD.IDX]<-sqrt(varmax[SIGN.SD.IDX]*fit$upperbound_var_factor))
    autostart[SIGN.SD.IDX] <- 0 
  }
 
#
  
  lb.s.ls.f <- fit$lowerbound_scale_ls_factor
  up.s.f <- fit$upperbound_scale_factor
  if ("longitude" %in% neu[[1]]$coordunits) {
    if ("km" %in% neu[[1]]$new_coordunits) {
      lb.s.ls.f <- lb.s.ls.f / 40 # 40 = approx 7000 / 180
      up.s.f <- up.s.f * 40
    } else if ("miles" %in% neu[[1]]$new_coordunits) {
      lb.s.ls.f <- lb.s.ls.f / 20 # 20 = approx 7000 / 180
      up.s.f <- up.s.f * 20      
    } else stop("unknown coordinate transformation")
  }
  if (any(idx <- ptype == DIAGPARAM)) {
    lower[idx] <- 1 / (up.s.f * maxdistances)
    upper[idx] <- lb.s.ls.f / mindistances
    autostart[idx] <- 8 / (maxdistances + 7 * mindistances)
  }
  
  if (any(idx <- ptype == ANISOPARAM)) {
    if (is.null(trafo))
      warning("The algorithms RandomFields transpose the matrix Aniso to aniso -- this may cause problems when applying transform to the anisotropy parameters. To be safe, use only the parameter anisoT in RMfit.")
    lower[idx] <- -lb.s.ls.f / mindistances
    autostart[idx] <- 0
  }

  if (any(SCALE.IDX)) {
    idx <- which(SCALE.IDX)
    lower[idx] <- mindistances / lb.s.ls.f
    upper[idx] <- maxdistances * up.s.f
    autostart[idx] <- (maxdistances + 7 * mindistances) / 8      
   }
 
#Print(lower)
 #  Print(autostart, users.guess)
#  Print(lower, upper)

#  Print(lower)
  
###########################        split       #######################
  if (printlevel>=PL_STRUCTURE) cat("\nsplit...")
  if (fit$split && length(autostart)>3) {
    
    Range <- if (is.scalar.dist || dist.given) {
      as.matrix(c(0, diff(range(if (dist.given) distances else
                                as.double(neu[[1]]$x) ))))
    } else rangex

    
    new.param <-
      if (is.null(users.guess)) autostart else users.guess


### Achtung!! delete.idx darf davor nur fuer trafo gesetzt werden!!

    aux.reg <- MODEL_MLESPLIT

    stopifnot(spdim == tsdim - time)
    
    .Call("SetAndGetModelInfo", aux.reg, list("Cov", model),
          as.integer(spdim), 
          is.dist.vector, # 3.4.13: warum war dist auf FALSE gesetzt gewesen?
          !is.dist.vector, time,
          xdimOZ, # 3.4.13: gaendert von dim. ? Warum war nicht xdim ?
          ##      Kleiber-Datensatz funktioniert nur mit xdimOZ!?
          fit$short, FALSE, TRUE,
          PACKAGE="RandomFields")

     
   stopifnot(ncol(Range) == ts.xdim)
    keep <- !delete.idx

#    Print(lower, keep)
#     Print(lower[keep], lower)
#   print(minmax)

    split <- ModelSplit(Reg=aux.reg, info.cov=info.cov, trafo=trafo,
                        variab=new.param[keep],
                        lower=lower[keep], upper=upper[keep],
                        rangex = Range,
                        modelinfo=list(ts.xdim=ts.xdim, tsdim=tsdim,
                          xdimOZ = xdimOZ, vdim=vdim,
                          is.dist.vector=is.dist.vector,
                          refined = fit$split_refined),
                        model=model)

    if (printlevel>=PL_STRUCTURE) cat("\nsplitted...")  

    if (length(split) == 1)
      stop("split does not work in this case. Use split=FALSE")
    if (length(split) > 1) {
      coord <- werte <- Xdistances <- NULL
      if (printlevel>=PL_STRUCTURE) {
        cat("\n")
      }

      if (printlevel >= PL_RECURSIVE) Print(split) #

#      Print(keep, lower, upper)

      return(recurs.estim(split=split, level=0,  Reg=aux.reg,
                          vdim=vdim,
                          lc=lc,
                          model=model,
                          x=if (!dist.given) lapply(neu, function(x) x$x), 
                          T=if (!dist.given) lapply(neu, function(x) x$T),
                          data= data, ###
                          lower= rep(NA, sum(keep)),
                          upper= rep(NA, sum(keep)),
                          users.lower = {
                            if (is.null(users.lower)) rep(-Inf,sum(keep))
                            else users.lower[keep]
                          },
                          users.upper = {
                            if (is.null(users.upper)) rep(Inf, sum(keep))
                            else users.upper[keep]
                          },
                          guess=new.param[keep], # setzt default werte
                          distances=if (dist.given) distances,
                          grid=neu[[1]]$grid,
                          bc_lambda=bc_lambda,
                          lsq.methods = LSQMETHODS,
                          mle.methods=mle.methods,
                          tsdim = tsdim,
                          optim.control=optim.control,
                          transform=transform,
                          trafo=trafo,
                          spConform = general$spConform,
                          practicalrange = general$practicalrange,
                          printlevel=printlevel,
                          minmax=minmax,
                          sdvar = apply(split[[1]]$all.p, 2,
                            function(x) {x[!ALL.SDVAR] <- FALSE; x})
                          ## sdvar in spaltenrichtung vdim, in zeilenrichtung
                          ## die parameter
                          ))
    }
  }

  delete.idx <- which(delete.idx) ## !!
  .Call("Delete_y", Reg, PACKAGE="RandomFields")
 
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


  if (printlevel>=PL_STRUCTURE) cat("\nauto...")  
  
  autopar <- autostart
  autopar[autopar == 0] <- 1
  varnugNA <- ## is.na(nugget) and is.na(var)
    globalsigma <- ## nugget==0
      sillbounded <- ## !is.na(sill)
        FALSE

  ##Print(optim_var_elimination)
  if (optim_var_elimination == "try" || optim_var_elimination == "yes" ||
      (optim_var_elimination == "respect bounds" && is.null(users.lower) &&
       is.null(users.upper))) { ## useful for lsq and mle methods

    ssm <- StandardSimpleModel(model, tsdim=tsdim, aniso=FALSE, addnugget=FALSE)

    if (!is.character(ssm)) {
      optim_var_elimination <- "yes"
      l <- length(ssm)
      nonugget <- (!(ssm[[l]][[length(ssm[[l]])]][[1]] %in% ZF_NUGGET) ||
                   !is.finite(ssm[[l]]$var) && ssm[[l]]$var == 0.0)
      novariance <- (ssm[[2]][[length(ssm[[2]])]][[1]] %in% ZF_NUGGET ||
                     !is.finite(ssm[[2]]$var) && ssm[[2]]$var == 0.0)
    } 
  } else {
    ssm <- NULL
    nonugget <- novariance <- FALSE
  }

  if (length(var.idx) > 1 || sum(SCALE.IDX)>1 || sum(ptype == DIAGPARAM)>0 ||
      sum(ptype==ANISOPARAM)>0 || sum(ptype==INTEGERPARAM)>0 ||
      sum(ptype==TRENDPARAM)>1 || length(nugget.idx)>1 ||
      length(MIXED.IDX) > 0 || sum(ptype==MIXEDPARAM)>0) {
    optim_var_elimination <- "never"
  }

  var.global <- var.idx
  sillbounded <- !is.na(sill)

  if (sillbounded && optim_var_elimination != "yes")
    stop("sill can not be bounded for complicated models. Check also the value of 'optim_var_elimination' which should be 'try'.")

#  Print(optim_var_elimination, sillbounded)

  
  if (optim_var_elimination == "yes") {
    if (printlevel>=PL_STRUCTURE) cat("\nstandard style settings...")
  
    if (sillbounded) {
      ## only VARIANCE need to be optimised
      ## NUGGET = SILL - VARIANCE
      stopifnot(sill > 0)

      if (length(nugget.idx) != 1 && length(var.idx) != 1) {
        if (length(var.idx)==0 && length(nugget.idx)==0)
          stop("variance and nugget are given, so sill must be NA")
        else 
          stop("Sill fixed. Then variance and nugget must be both NA.")
      }

      autopar[c(nugget.idx, var.idx)] <-
        autostart[var.idx] <- sill / length(c(nugget.idx, var.idx))
      upper[nugget.idx] <- sill
      delete.idx <- c(delete.idx, var.idx)
      lower[nugget.idx] <- 0
      
      trafo <- function(x) {
        ## warum  mixed.idx  ?
        y <- numeric(length(x) + 1 + length(MIXED.IDX) )
        y[-c(var.idx, MIXED.IDX)] <- x
        y[var.idx] <- sill - y[nugget.idx]
        y[MIXED.IDX] <- 1
        y
      }
   
    } else { ## not sill bounded
      if (length(var.idx) == 1) {
        if (varnugNA <- length(nugget.idx) == 1) {
          if (!is.null(users.guess)) {
            users.guess[nugget.idx] <-
              users.guess[nugget.idx] / sum(users.guess[c(var.idx, nugget.idx)])
          }
          
          ## of interest currently 
          ## both NUGGET and VARIANCE have to be optimised
          ## instead optimise the SILL (which can be done by
          ## a formula) and eliminationate the part of NUGGET
          ## (which will be stored in NUGGET) !
          ##  consequences:
          ## * eliminationtation space increased only by 1, not 2
          ## * real nugget: NUGGET * SILL
          ## *     variance: (1-NUGGET) * SILL
          autostart[nugget.idx] <- 1 / 2 ## lassen, da nicht standard Alg.
          lower[nugget.idx] <- 0
          upper[nugget.idx] <- 1
          delete.idx <- c(delete.idx, var.idx)
          trafo <- function(x) {
            y <- numeric(length(x) + 1 + length(MIXED.IDX))
            y[-c(var.idx, MIXED.IDX)] <- x
            y[var.idx] <- 1.0 - y[nugget.idx]
            y[MIXED.IDX] <- 1
            y
          }
        } else { ## not sillbounded, is.na(variance), !is.na(nugget),
          ## but optim_var_elimination
          if (globalsigma <- nonugget) {
            ## here SILL==VARIANCE and therefore, variance
            ## can be eliminationated by a formula without increasing the dimension

            ## more than only the variance has to be eliminationated
            delete.idx <- c(delete.idx, var.idx)
            # if (length(lower) == 0) stop("trivial case not programmed yet")
            trafo <- function(x) {
              y <- numeric(length(x) + 1 + length(MIXED.IDX))
              y[-c(var.idx, MIXED.IDX)] <- x
              y[c(var.idx, MIXED.IDX)] <- 1.0 ## 30.12.09 ex: vardata 
              y
            }
          } else { ## not sillbounded, is.na(variance), nugget!=0
            l <- length(ssm)
            nugget <- ssm[[l]]$var
            lower[var.idx] <-
              pmax(lower[var.idx], 0)
            #(varmin[var.idx]-nugget) / fit$lowerbound_var_factor)
            if (lower[var.idx] < fit$lowerbound_sill) {
              if (printlevel>=PL_SUBIMPORTANT) {
                if (printlevel>=PL_SUBIMPORTANT)
                  cat("low.var=",lower[var.idx]," low_sill",fit$lowerbound_sill,
                      "\ eliminationated variance from data=", varmin[var.idx],
                      "nugget=", nugget, "\n")
                message("The given nugget effect might not be too small.")
              }
              lower[var.idx] <- fit$lowerbound_sill
            }
            autostart[var.idx] <- autopar[var.idx] <- lower[var.idx]
          }
        }
      } else { ## not sillbounded, !is.na(variance), optim_var_elimination
        if (novariance) {
          if (any(ptype %in% c(SCALEPARAM, ANYPARAM, CRITICALPARAM))) 
            stop("If variance=0, eliminationating scale or model parameters does not make sense")          
          lower[var.idx] <- 0
        }
        if (length(nugget.idx)==1) { ## and !is.na(param[VARIANCE])
          l <- length(ssm)
          variance <- ssm[[2]]$var
          lower[nugget.idx] <-
            pmax(0, (varmin[nugget.idx] - variance) / fit$lowerbound_var_factor)
          if (lower[nugget.idx] < fit$lowerbound_sill) {
            if (printlevel>PL_SUBIMPORTANT)
              cat("\nlower nugget bound=", lower[nugget.idx],
                  " < lower sill bound=", fit$lowerbound_sill,
                  " -- is the variance given correctly?\n",sep="")
            ## warning("Has param[VARIANCE] been given correctly?!")
            lower[nugget.idx] <- fit$lowerbound_sill
            autostart[nugget.idx] <- autopar[nugget.idx] <- lower[nugget.idx]
          }
  #        lower[nugget.idx] <- 0
        }
      } ##  else { ## not sillbounded, !is.na(variance)
    } ## else { ## not sill bounded

    stopifnot( varnugNA + globalsigma + sillbounded <= 1)
      
  } else { ### optim_var_elimination != "yes"
   
    globalsigma <- !anyRandomEffects &&
    ((modelinfo$name %in% DOLLAR && is.na(modelinfo$param$var)) ||
      (length(idx.error)==1 &&
         (modelinfo$sub[[idx.error]]$name %in% DOLLAR &&
          is.na(modelinfo$sub[[idx.error]]$param$var)
          ||
          modelinfo$sub[[idx.error]]$name %in% ZF_PLUS &&
          length(modelinfo$sub[[idx.error]]$sub) == 1 &&
          modelinfo$sub[[idx.error]]$sub[[1]]$name %in% DOLLAR &&
          is.na(modelinfo$sub[[idx.error]]$sub[[1]]$param$var) 
          )))
 
    globalsigma <- globalsigma && is.null(transform) &&
                   length(c(var.idx, nugget.idx)) == 1

    if (globalsigma) {
      delete.idx <- c(delete.idx, var.global)       
      trafo <- function(x) {
        y <- numeric(length(x) + 1 + length(MIXED.IDX))
        y[-c(var.global, MIXED.IDX)] <- x
        y[c(var.global, MIXED.IDX)] <- 1.0
        y
      } 
    }
  }
  
#  Print(lower, upper, "eXXnd"); kkkk
  
 
   RFoptions(fit.optim_var_elimination =# in case of a recursive call of fit.gauss
            if (optim_var_elimination == 'yes') 'yes' else 'never') 
  

#  Print(lower, users.lower, upper, users.upper, delete.idx)

  globalsigma <- globalsigma || varnugNA

  if (is.null(users.lower)) {
    users.lower <- rep(-Inf, length(lower))
  } else {
    idx <- !is.na(users.lower)
    lower[idx] <- users.lower[idx]
    users.lower[!idx] <- -Inf
  }
  if (is.null(users.upper)) {
    users.upper <- rep(Inf, length(upper))
  } else {
    idx <- !is.na(users.upper)
    upper[idx] <- users.upper[idx]
    users.upper[!idx] <- Inf
  }
  bounds <- apply(abs(minmax[, 5:6, drop=FALSE]), 1, max)
  if (solvesigma) delete.idx <- c(delete.idx, vareffect)



#  Print(lower, upper, users.lower, users.upper, "end"); 
  
    paramnames <- varnames
  if (length(delete.idx)>0) {    
    upper <- upper[-delete.idx]
    lower <- lower[-delete.idx]
    users.lower <- users.lower[-delete.idx]
    users.upper <- users.upper[-delete.idx]
    autostart <-autostart[-delete.idx]
    autopar <- autopar[-delete.idx]
    varnames <- varnames[-delete.idx]
    SCALE.IDX <- SCALE.IDX[-delete.idx]
    SDVAR.IDX <- SDVAR.IDX[-delete.idx]
    ptype <- ptype[-delete.idx]
    bounds <- bounds[-delete.idx]    
  }


#  Print(lower, users.lower, upper, users.upper, delete.idx, users.guess); kkk

  
  if (any(autostart<lower) || any(autostart>upper)) {
    if (printlevel >= PL_ERRORS)
      Print(cbind(lower, autostart, upper)) #, orig.lower,orig.upper)#
    autostart <- pmin(upper, pmax(lower, autostart))
  }

  nvarWithoutbc <- length(lower)
  if (lambdaest) {
    varnames <- c(varnames, "BoxCox")
    paramnames <- c(paramnames, "BoxCox")
    autostart <- c(autostart, 1)
    autopar <- c(autopar, 1)
    lower <- c(lower, fit$bc_lambda_lb)
    upper <- c(upper, fit$bc_lambda_ub)
    users.lower <- c(users.lower, fit$bc_lambda_lb)
    users.upper <- c(users.upper, fit$bc_lambda_ub)
    SDVAR.IDX <- c(SDVAR.IDX, rep(FALSE, length(fit$bc_lambda_lb)))
    SCALE.IDX <- c(SCALE.IDX, rep(FALSE, length(fit$bc_lambda_lb)))
  }
  n.variab <- length(lower)

  if (any(idx <- lower >= upper)) {
    lu <- cbind(lower=lower, upper=upper, idx)  
    stop(paste("Some lower bounds for the parameters are greater than ",
               "or equal to the upper bound\n",
               paste(collapse="\n ", dimnames(lu)[[1]], ":",
                     apply(lu, 1, function(x)
                           paste("\tlower=", x[1], ",\tupper=", x[2],
                                 if (x[3]) "  \tnot ok!", sep=""))
                     )
               ))
  }


  fill.in <-  trafo(autostart)
  if (length(MIXED.IDX) > 0) {  ## Aufruf .C(Cov*Loc wird vor optim benoetigt
    ## somit muessen die NA's (jedoch nur) in mixed.var gesetzt sein
    autostart[MIXED.IDX] <- 1.0
  }
  .C("PutValuesAtNA", Reg, trafo(autostart), PACKAGE="RandomFields")   

 
######################################################################
###                                 Covariates                     ###
######################################################################
  if (printlevel>=PL_STRUCTURE) cat("\nCovariates...")

  
##############  Sinv, etc              #################
## 

  ##  ToDo !!! : check correct dimensions of X, etc

 
 ## startCoVar <- endCoVar <-  matrix(nr=sets, ncol=length(allcomponents))
  logdetbase <- 0

  ## all beta of mixed compents are fixed for each set
  ## and drawn independently among sets
  ## error part is independent for each repetition in/for each set

 
  nCoVar <- endCoVar <- startCoVar <- startXges <- endXges <-
      matrix(0, nrow=sets, ncol=length(effect))
  
  if (anymixedeffects) {     
    cs <- 0
    s <- numeric(length(effect))
    for (i in 1:sets) {   
      for (k in allcomponents) {
        s[k] <-
          if (effect[k] >= PrimitiveModel) 0
          else {
            stopifnot(effect[k] != DetTrendEffect,
                      effect[k] != DeterministicEffect)
            
            if (effect[k] == FixedTrendEffect) {
              s0 <- length(modelinfo$sub[[k]]$param$mean) +
                length(modelinfo$sub[[k]]$param$plane) +
                  length(modelinfo$sub[[k]]$param$polydeg)
              if (s0 > 0) s0 else {
                stop("arbitrary fct not programmed yet")
              }
            } else {
              if (length(modelinfo$sub[[k]]$param$X) > 0)
                ncol(modelinfo$sub[[k]]$param$X[[i]])
              else lc[i] * vdim
            }
          }
      }
      
      nCoVar[i, ] <- s
      endXges[i, ] <- cumsum(s) 
      endCoVar[i, ] <- cs + endXges[i, ]   
      startXges[i, ] <- endXges[i, ] - s + 1
      startCoVar[i, ] <- endCoVar[i, ] - s + 1
      cs <- cs + sum(s)
    }
   
    Ny <-
      sapply(modelinfo$sub,
             function(x) {
               if((x$name %in% ZF_MIXED) && !is.null(x$param$X))
                 sapply(x$param$X, nrow)
               else if (x$name %in% ZF_TREND) rep(-1, sets)
               else rep(0, sets) # > SpVarEffect
             }) / vdim
    if (is.vector(Ny)) dim(Ny) <- c(1, length(Ny))   
  }

  nTotalComp <- apply(nCoVar, 2, sum)
  nCoVarSets <- apply(nCoVar, 1, sum)
  nCoVarAll <- as.integer(sum(nCoVarSets))


  Sinv <- list() ## coordinates
  Xges  <- list()
 
  for (i in 1:sets) {
    Xges[[i]] <- Sinv[[i]] <- list()
    if (anymixedeffects) {
      for (m in 1:len.rep[i]) {
        Xges[[i]][[m]] <- matrix(nrow=len.idx.na[[i]][m], ncol=nCoVarSets[i])
        for (k in mmcomponents) {
          if (Ny[i, k] >= 0) { ## -1==trend
            if (Ny[i, k] == 0) {
              stop("Unclear why the algorithm has entered this part")
              stopifnot(nCoVar[i, k] == lc[i] * vdim) # or: len.idx.na[i]??
            } else {
              if (Ny[i, k] != lc[i] * vdim) { 
                stop("length of data does not match X matrix")
              }
            }
          }
          
          Xges[[i]][[m]][, startXges[i, k] : endXges[i, k] ] <- 
            if (Ny[i,k] <= 0) { ## -1==trend
              if (trends[k]) {                
                GetTrend(coord=coord[[i]], idx=idx.na[[i]][[m]],
                         model=if(model[[1]] %in% ZF_PLUSSELECT)
                         model[[k]] else model,
                         vdim=vdim, param=modelinfo$sub[[k]]$param)
              } else {
                stop("unclear why the algorithm has entered this part")
                diag(nCoVar[i, k])
              }
            } else modelinfo$sub[[k]]$param$X[[i]][idx.na[[i]][[m]], ]
        }
      }
    }
    
    for (k in allcomponents) {
       if (effect[k] > FixedEffect) {
        if (effect[k] < SpaceEffect) {
          dummy <- modelinfo$sub[[k]]$sub[[1]]
          cov <- if (dummy$name %in% DOLLAR) {
            dummy$sub[[1]]$param$M[[i]]
          } else {
            dummy$param$M[[i]]
          }   
          Sinv[[i]][[k]] <- Solve(cov, LINPACK=TRUE)    
          logdetbase <- logdetbase + LOGDET
        } else {
          sp.eff <- effect[k] == SpaceEffect || effect[k] == SpVarEffect
          Sinv[[i]][[k]] <- if (sp.eff) NULL else list()
        }
      }
    }
  } # i in 1:sets

   
  eff <- efftmp <- rep(0, RemainingError + 1)
  for (k in mmcomponents) {
    if (effect[k] <= SpVarEffect) {
      eff[effect[k] + 1] <- eff[effect[k] + 1] + 1
    }
  }

  betanames <- character(nCoVarAll)
  for (i in 1:sets) {    
    for (k in mmcomponents) {
      if (effect[k] <= SpVarEffect) {
        efftmp[effect[k] + 1] <- efftmp[effect[k] + 1] + 1
        bn <- paste(EffectName[effect[k] + 1],
                    if (eff[effect[k] + 1] > 1) efftmp[effect[k] + 1], sep="")
        if (nCoVar[i, k]>1) bn <- paste(bn, 1:nCoVar[i, k], sep=".")
        betanames[startCoVar[i, k] : endCoVar[i, k]] <-
          if (sets == 1) bn else paste(i, bn, sep="")
      }
    }
  }
  rm("eff", "efftmp")



######################################################################
######################################################################
###                     Eliminationation part itself                     ###
######################################################################
######################################################################

  ## check optim.control 
  ## parscale will give the magnitude of the parameters to be eliminationated
  ##     passed to optim/optimise so that the optimiser eliminationates
  ##     values around 1 ##to do: irgendwo in paper die wichtigkeit beschreiben?

  parscale <- ParScale(optim.control, lower=lower, upper=upper)
  

 
  fit.fnscale <- optim.control$fnscale

  
  
  if (length(optim.control)>0) {
    opt.control <- optim.control
    stopifnot(is.list(opt.control))
    forbidden.param <- c("parscale", "fnscale", "algorithm")
    ## fnscale=-1 turns the problem into a maximisation problem, see below
    forbidden <- which(!is.na(pmatch(names(opt.control), forbidden.param)))
    forbidden.opt <- opt.control[forbidden]  
    if (length(forbidden) > 0)  opt.control <- opt.control[-forbidden]
  } else {
    opt.control <- list()
  }

  if (length(fit$algorithm) > 0 && fit$algorithm != "")
      opt.control$algorithm <- fit$algorithm
  if (length(optim.control)>0) {
    if (length(forbidden.opt$algorithm) > 0)
      opt.control$algorithm <- forbidden.opt$algorithm
  }


###################  preparation  ################
  if (printlevel>=PL_STRUCTURE) cat("\npreparing fitting...")
  ## methods
  formals <- formals()
  allprimmeth <- c("autostart", "users.guess")
  nlsqinternal <- 3 ## cross checked after definition of weights below
  alllsqmeth <- (if (variogram.known)
                 c(LSQMETHODS[-length(LSQMETHODS)],
                   paste("internal", 1:nlsqinternal, sep=""))
                 else NULL)
  allmlemeth <- eval(formals$mle.methods)

  
  # allcrossmeth <- eval(formals$cross.methods)
    allcrossmeth <- NULL
  allmethods <- c(allprimmeth, alllsqmeth, allmlemeth, allcrossmeth)

  ## how preceding methods have been considered ?
  ## note cm is used again at the very end when error checking
  cm <- cumsum(c(0, length(allprimmeth), length(alllsqmeth),
                     length(allmlemeth), length(allcrossmeth)))
  cm <- cbind(cm[-length(cm)] + 1, cm[-1])
  cm <- apply(cm, 1, function(x) x[1] : x[2])
  names(cm) <- c("prim", "lsq", "mle", "cross")

  methodprevto <-
    if (fit$only_users) {
      list(lsq="users.guess",mle="users.guess",cross="users.guess")
    } else list(lsq=c(cm$prim),
              mle=c(cm$prim, cm$lsq),
              cross=c(cm$prim, cm$lsq, cm$cross)
              )

  ## index (start, end) to the various categories of
  ## information to be stored
  IDX <- function(name) {idx <- tblidx[[name]]; idx[1]:idx[2]}
  tblidx <- cumsum(c(0,
                     n.variab, # variables used in algorithm
                     n.variab, # their lower bounds
                     n.variab, # ... and upper bounds
                     n.variab,  # sd of variabs                     
                     ncovparam,# param values to be eliminationated
                     rep(1, length(allmethods) - length(allprimmeth)),#method
                     ##                                                 score
                     1, 1, 1, ## AIC, AICc, BIC,
                     nCoVarAll, # coeff to eliminationated for covariates, i.e.
                     ##           mixed effects and trend parameters
                     ncovparam ## sd of params
                    ))

 
  tblidx <- rbind(tblidx[-length(tblidx)] + 1, tblidx[-1])
  idx <- tblidx[1, ] > tblidx[2, ]
  tblidx[, idx] <- 0
 
  dimnames(tblidx) <- list(c("start", "end"),                           
                           c("variab", "lower", "upper",
                             "sdvariab",
                             "param", 
                             allmethods[-1:-length(allprimmeth)],
                             "AIC", "AICc", "BIC",
                             "covariab",
                             "sdparam"
                             ##,  "lowbeta", "upbeta", only used for
                             ## cross-validation
                             ))
  maxtblidx <- max(tblidx)
  tblidx <- data.frame(tblidx)
  
   ## table of all information; col:various method; row:information to method
  
  tablenames <-
     c(if (n.variab > 0) paste("v", varnames, sep=":"),        
       if (n.variab > 0) paste("lb", varnames, sep=":"),
       if (n.variab > 0) paste("ub", varnames, sep=":"),
       if (n.variab > 0) paste("sdv", varnames, sep=":"),        
      minmax.names,
##       if (nrow(minmax) > 0) paste("p", 1:nrow(minmax), sep=":")
  ##                       else minmax.names,
      allmethods[-1:-length(allprimmeth)],
       "AIC", "AICc", "BIC",
      betanames,
      ## do not try to join the next two lines, since both
      ## varnames and betanames may contain nonsense if
      ## n.variab==0 and nCoVarAll==0, respectively
       paste("sd", minmax.names, sep=":")
       )
 
  param.table <- data.frame(matrix(NA, nrow=maxtblidx, ncol=length(allmethods),
                                   dimnames=list(tablenames, allmethods)))


  fit.fnscale <- if (is.null(fit.fnscale)) rep(NA, length(allmethods)) else
                  -abs(fit.fnscale)

#  Print(tblidx, tblidx, n.variab, varnames, tablenames, param.table)


#############################################################
## end preparation; remaining part is eliminationation  ###########
#############################################################

##################################################
###############    PRIMITIVE METHODS   ###########
##################################################
 
 
  ##****************    autostart    *****************
  if (printlevel>=PL_STRUCTURE) cat("\nautostart...")
  M <- "autostart"
  
  default.param <- param.table[[M]][IDX("variab")] <- autostart
  param.table[[M]][IDX("param")] <- trafo(autostart) 


  ## ****************    user's guess    *****************
  if (!is.null(users.guess)) {
    M <- "users.guess"
    if (length(delete.idx) > 0) users.guess <- users.guess[-delete.idx]
    idx <- users.guess < lower | users.guess > upper
    if (any(idx)) {
      if (recall) users.guess <- NULL
      else if (general$modus_operandi ==  MODENAMES[careless + 1]) {
        lower <- pmin(lower, users.guess)
        upper <- pmax(upper, users.guess)
      } else {
        m <- cbind(lower, users.guess, upper, idx)
        dimnames(m) <- list(rep("", length(lower)),
                            c("lower", "user", "upper", "outside"))
        cat("\n")
        print(m) ## nicht loeschen!
        
        stop("not all users.guesses within bounds\n change values of `lower' and `upper' or \nthose of the `lowerbound*.factor's and `upperbound*.factor's")
      }
    }

    if (!is.null(users.guess)) {
      param.table[[M]][IDX("variab")] <- users.guess
      param.table[[M]][IDX("param")] <- trafo(users.guess)
    }
  }



  MLELB <- LSQLB <- lower
  MLEUB <- LSQUB <- upper

##################################################
################### Empirical Variogram   ########################
  lsqMethods <- NULL
  ev <- list()

  if (variogram.known) {    
    if (printlevel>=PL_STRUCTURE) cat("\nempirical variogram...\n")
    sd <- vector("list", vdim)
    emp.vario <- vector("list", vdim)
    index.bv <- NULL

    cum <- rep(0, sets)
    for (j in 1:vdim) {
      for (i in 1:sets) {
        m <- 1
        if ((length(idx.na.rep[[i]]) > 1) && printlevel>PL_IMPORTANT)
          message("Only first column(s) used for empirical variogram.")
        
        idx <- idx.na[[i]][[m]][idx.na[[i]][[m]] > (j-1) * lc[i] &
                                idx.na[[i]][[m]] <= j * lc[i] ]

        lidx <- length(idx)
        W <- werte[[i]][idx, , drop=FALSE]
        if (nCoVarAll > 0) {
          oldwarn <- options()$warn
          options(warn=no.warning)
          X <- Xges[[i]][[m]][cum[i] + (1:lidx), , drop=FALSE]
          regr <- lsfit(X[, apply(X!=0, 2, any), drop=FALSE],
                        W, intercept=FALSE)
          options(warn=oldwarn)
          W <- regr$residuals
        }
        cum[i] <- cum[i] + lidx
        
        ##if (length(nphi)==1)
        ##  nphi <- c(0, nphi) # starting angle; lines per half circle
        ##if (length(ntheta)==1) ntheta <- c(0, ntheta) # see above
        if (length(ntime)==1) ntime <- c(ntime, 1)

        ntime <- ntime * T[[i]][2] ## time endpoint; step
        
        bin <-
          if (length(bins)>1) bins
          else c(-1, seq(0, fit$bin_dist_factor * maxdistances, len=bins+1))
        if (!dist.given) {
          ev[[i]] <-
            RFempiricalvariogram(t(coord[[i]][, idx - (j-1) * lc[i],
                                              drop=FALSE]),
                                 T=T[[i]], data=W, grid=FALSE,
                                 bin=bin,
                                 phi=if ((spdim>=2) && !isotropic) nphi,
                                 theta=if ((spdim>=3) && !isotropic) ntheta,
                                 deltaT=if (!is.null(T[[i]])) ntime,
                                 spConform=FALSE) # 7.1.11 coord -> t(coord) !
        } else {
          n.bin <- vario2 <- vario <- rep(0, length(bin))
          stopifnot(length(lc) == 1)
          k <- 1
          for (g in 1:(lc[1]-1)) {
            # cat(g, "")
            for (h in (g+1):lc) {
              idx <- sum(distances[[i]][k] > bin)
              n.bin[idx] <- n.bin[idx] + 2
              vario[idx] <- vario[idx] + mean((W[g] - W[h])^2, na.rm=TRUE)
              vario2[idx] <- vario2[idx] + mean((W[g] - W[h])^4, na.rm=TRUE)
              k <- k + 1
            }
          }
          vario <- vario / n.bin ## Achtung! n.bin ist bereits gedoppelt
          
          sdvario <-
            sqrt(pmax(0, vario2 / n.bin / 2.0 - vario^2)) ## numerische Fehler
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
                       rep(1, length(ev))) / n.bin) # j:vdim; ev contains the sets
      emp.vario[[j]] <- sapply(ev, function(x) x$emp.vario * x$n.bin) %*%
        rep(1, length(ev)) / n.bin
      index.bv <- cbind(index.bv, !is.na(emp.vario[[j]])) ## exclude bins without
      ##                                                    entry
    } # j in 1:vdim
   
    index.bv <- apply(index.bv, 1, all)
 
    if (sum(index.bv) <= 1)
      stop("not more than 1 value in empirical variogram that is not NA; check values of bins and bin_dist_factor")
    
     
    ## bei vdim>1 werden nur die diag-elemente von binned.vario gefuellt:
    bvtext <- paste("emp.vario[[", 1:vdim, "]][index.bv]", collapse=", matrix(0, nrow=sum(index.bv), ncol=vdim) ,", sep="")
    binned.variogram <- eval(parse(text=paste("cbind(", bvtext,")", sep="")))
    binned.variogram <- matrix(t(binned.variogram))[ , , drop=TRUE]

    stopifnot(sum(binned.variogram) >0)

    bin.centers <- as.matrix(ev[[1]]$centers)

    if (!is.null(ev[[1]]$phi)) {
      if (spdim<2) stop("x dimension is less than two, but phi is given") 
      bin.centers <- cbind(as.vector(outer(bin.centers, cos(ev[[1]]$phi))),
                           as.vector(outer(bin.centers, sin(ev[[1]]$phi))))
    }
    if (!is.null(ev[[1]]$theta)) {
      if (spdim<3)
        stop("x dimension is less than three, but theta is given") 
      if (ncol(bin.centers)==1) bin.centers <- cbind(bin.centers, 0)
      bin.centers <- cbind(as.vector(outer(bin.centers[, 1],
                                           cos(ev[[1]]$theta))),
                           as.vector(outer(bin.centers[, 2],
                                           cos(ev[[1]]$theta))),
                           rep(sin(ev[[1]]$theta), each=nrow(bin.centers)))
    } else {
      
      ##  warning("must be ncol()")      
      if (ncol(bin.centers) < spdim) { # dimension of bincenter vector
        ##                       smaller than dimension of location space      
        bin.centers <- 
          cbind(bin.centers, matrix(0, nrow=nrow(bin.centers),
                                    ncol=spdim - ncol(bin.centers)
                                    ))
      }
    }
    if (!is.null(ev[[1]]$T)) {
      bin.centers <-
        cbind(matrix(rep(t(bin.centers), length(ev[[1]]$T)), byrow=TRUE,
                     ncol = ncol(bin.centers)),
              rep(ev[[1]]$T, each=nrow(bin.centers)))
    }

    if (isotropic && cartesian) {
      bin.centers <- as.matrix(sqrt(apply(bin.centers^2, 1, sum)))
    }
    bin.ts.xdim <- ncol(bin.centers)
    bin.centers <- as.double(t(bin.centers[index.bv, ])) #
    ## es muessen beim direkten C-aufruf die componenten der Punkte
    ## hintereinander kommen (siehe auch variable coord, Xdistance). Deshalb t()
    ## aus den vdim sd's muss ein Gewicht gemacht werden
    evsd <- sapply(sd, function(x) x)^2    
    evsd <-
      if (is.matrix(evsd)) rowSums(evsd / rowSums(is.finite(evsd)), na.rm=TRUE)
      else evsd[!is.finite(evsd)] <- 0    
    evsd[evsd==0] <- 10 * sum(evsd, na.rm=TRUE) ## == "infinity"
    evsd <- as.double(evsd)

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
  } else { # not trans.inv
    if (!is.null(lsq.methods)) warning("submethods are allowed")
  }

 
##################################################
###################  AUTOSTART  ##################
  if (printlevel>=PL_STRUCTURE) cat("\ncovariab...")

  ## see above for the trafo definitions
  ##
  ## zuerst regression fit fuer variogram,
  ## dann schaetzung der Parameter, dann berechnung der covariates
  ## T.o.D.o.: kann verbessert werden durch einschluss der Trendschaetzung
  ## Xges auf RFsimulate basieren
  
#  if (anymixedeffects) {
#    for (i in 1:sets) {
#      for (m in 1:len.rep[i]) {
#        Xges[[i]][[m]] <- Xges[[i]][[m]] [idx.na[[i]][[m]], , drop=FALSE]
#      }
#    }
#  }

    ##****************    autostart    *****************
  
  MLEVARIAB <- autostart
  assign("LINEARLARPS2", rep(1, length(MIXED.IDX)), envir=ENVIR)
  param.table[[M]][IDX("covariab")] <- get.covariates(autostart)



  ## ****************    user's guess    *****************
  if (!is.null(users.guess)) {
    MLEVARIAB <- users.guess
    param.table[[M]][IDX("covariab")] <- get.covariates(users.guess)
  }



 
##################################################
#######################  LSQ  ####################
   ##***********   eliminationation part itself   **********     
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
  if (printlevel>=PL_STRUCTURE) cat("\neliminationation part...")
  if (trans.inv && !is.null(lsq.methods)) {
    
    #Print(is.dist.vector, bin.ts.xdim)

   .Call("SetAndGetModelInfo", Reg, list("Cov", model),
          spdim,
          ## die nachfolgenden beiden Festlegungen widersprechen sich,
          ## werden aber in unterschiedlichen Situationen benoetigt.

         ##
       #  is.dist.vector, # wo war dies notwendig??
         bin.ts.xdim==1, # fuer Bsp von Emilio notwendig!
         FALSE,
         time, bin.ts.xdim - time,
         fit$short, FALSE, TRUE,
         PACKAGE="RandomFields")

      
    LSMIN <- Inf
    lsqMethods <- LSQMETHODS[pmatch(lsq.methods, LSQMETHODS)]
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

   firstoptim <- TRUE
   for (M in c(lsqMethods[1], alllsqmeth)) {
      if (!(M %in% lsqMethods)) next;
      if (printlevel>=PL_STRUCTURE) cat("\n", M) else cat(pch)

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
        min <- Inf
        min.variab <- NULL
        for (i in methodprevto$lsq) { ## ! -- the parts that change if
          ##                               this part is copied for other methods

          ##   Print(M, param.table[IDX("variab"), i])
          if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {
             value <- LStarget(variab) ##
             ## Print(value, min.variab, LSVARIAB, min, LSMIN,is.finite(value))
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

        ##Print( methodprevto$lsq, param.table[IDX("variab"),methodprevto$lsq])
        ##Print(min.variab, LSVARIAB, min, LSMIN)
        
        stopifnot(min.variab==LSVARIAB, min==LSMIN) ## check
        
        fnscale <- if (is.null(fit.fnscale) || is.na(fit.fnscale[M]))
          min else fit.fnscale[M]
        #Print("lsq", fnscale)

        lsq.optim.control <-
          c(opt.control, list(parscale=parscale, fnscale=fnscale))
         
        OPTIM(LSVARIAB, LStarget, lower = LSQLB, upper = LSQUB,
                control=lsq.optim.control,optimiser=optimiser,silent=silent)
      
      } # n.variab > 0
      options(show.error.messages = show.error.message)  
      ## side effect: minimum so far is in LSMIN and LSPARAM
      ## even if the algorithm finally fails


 
      if (is.finite(LSMIN)) {
        param.table[[M]][tblidx[[M]][1]] <- LSMIN
        param.table[[M]][IDX("variab")] <- LSVARIAB
        param.table[[M]][IDX("param")] <- LSPARAM

        ps <- abs(LSVARIAB)
        zero <- ps == 0
        parscale[!zero] <- ps[!zero]
        #        Print(ps, zero, parscale)
       
        # Print(parscale)
        
      } else {
        param.table[[M]] <- if (n.variab==0) NA else NaN
      }

      param.table[[M]][IDX("covariab")] <- get.covariates(LSVARIAB)

      firstoptim <- FALSE
    } # for M

    .Call("SetAndGetModelInfo", Reg, list("Cov", model),
          spdim, is.dist.vector, !is.dist.vector, time, xdimOZ,
          fit$short, FALSE, TRUE,
          PACKAGE="RandomFields")
    
 } # trans.inv



  if (length(alllsqmeth) > 0) {
    ps <- matrix(NA, nrow=length(IDX("variab")), ncol=length(alllsqmeth))
    for (iM in 1:length(alllsqmeth)) {
      M <- alllsqmeth[iM]

#      Print(M, tblidx[[M]], param.table[[M]][tblidx[[M]][1]])
      
      if (!is.na(param.table[[M]][tblidx[[M]][1]])) {
        ps[ , iM] <- abs(param.table[[M]][IDX("variab")])
      }
    }
    ps <- apply(ps, 1, median, na.rm=TRUE)
    parscale <- ParScale(optim.control, ps, lower, upper)
  }


##################################################
### optional parameter grid for MLE and CROSS  ###
  if (printlevel>=PL_STRUCTURE) cat("\nmle param...")
  
  
  idx <- IDX("variab")
  gridmax <- as.matrix(param.table[idx, cm$lsq])
  if (!any(is.finite(gridmax))) gridmax <- param.table[idx, , drop=FALSE]
   
  gridmin <- apply(gridmax, 1, min, na.rm=TRUE)
  gridmax <- apply(gridmax, 1, max, na.rm=TRUE)

  gridbound <- lower
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  abase <- 0.25
  a <- is.na(gridmin[idx]) * (1-abase) + abase
  ## maybe there have not been any lsq eliminationate; then a=1
  gridmin[idx] <- (1-a) * gridmin[idx] + a * gridbound[idx]
  stopifnot(all(gridmin >= lower))
  gridbound <- upper
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  a <- is.na(gridmax[idx]) * (1-abase) + abase
  gridmax[idx] <- (1-a) * gridmax[idx] + a * gridbound[idx]
  stopifnot(all(gridmax <= upper))

  

##################################################
###################   MLE    #####################


  if (printlevel>=PL_STRUCTURE) cat("\nMLE XXX...")
  MLEtarget <- NULL
  mleMethods <- (if (is.null(mle.methods)) NULL else
                 allmlemeth[pmatch(mle.methods, allmlemeth)])
  if ("reml" %in% mleMethods && nCoVarAll == 0)
    mleMethods <- c("ml", "reml", "rml")

    
  .Call("SetAndGetModelInfo", Reg, list("Cov", model),
        spdim,
        ## die nachfolgenden beiden Festlegungen widersprechen sich,
        ## werden aber in unterschiedlichen Situationen benoetigt.
        is.dist.vector, # wo war dies notwendig??
        !is.dist.vector, time, xdimOZ - time,
        fit$short, FALSE, TRUE,
        PACKAGE="RandomFields")



  ## lowerbound_scale_ls_factor <  lowerbound_scale_factor, usually
  ## LS optimisation should not run to a boundary (what often happens
  ## for the scale) since a boundary value is usually a bad initial
  ## value for MLE (heuristic statement). Therefore a small
  ## lowerbound_scale_ls_factor is used for LS optimisation.
  ## For MLE eliminationation we should include the true value of the scale;
  ## so the bounds must be larger. Here lower[SCALE] is corrected
  ## to be suitable for MLE eliminationation

  #Print(MLELB, MLEUB, SCALE.IDX)
  if (FALSE) {
    ## 12.11.13 ausgeblendet, da die Grenzen zu eng werden koennen,
    ## so dass user's guess oder was in lsq gefunden wurde
    ## ausserhalb der neuen MLE bounds liegt
    if (any(SCALE.IDX)) {    
      MLELB[SCALE.IDX] <- MLELB[SCALE.IDX] *
        fit$lowerbound_scale_ls_factor / fit$lowerbound_scale_factor
    } else if (any(idx <- ptype == DIAGPARAM | ptype == ANISOPARAM))
      MLEUB[idx] <- MLEUB[idx] *
        fit$lowerbound_scale_factor / fit$lowerbound_scale_ls_factor
    MLELB[SCALE.IDX] <- pmin(users.lower[SCALE.IDX], MLELB[SCALE.IDX])
    MLEUB[SCALE.IDX] <- pmax(users.upper[SCALE.IDX], MLEUB[SCALE.IDX])
  }
  
  if (any(MLELB > MLEUB))
    stop("the users lower and upper bounds are too restricitve")

  ## fnscale <- -1 : maximisation
  for (M in c(allmlemeth)) {
 
    assign("LINEARLARPS2", rep(1, length(MIXED.IDX)), envir=ENVIR)
    if (!(M %in% mleMethods)) next;
    if (printlevel>=PL_STRUCTURE ) cat("\n", M) else cat(pch)
    param.table[[M]][IDX("variab")] <- default.param
    if (M!="ml" && !anyFixedEffectt) { 
      param.table[[M]] <- param.table[["ml"]]
      param.table[[M]][tblidx[[M]][1]] <- param.table[[M]][tblidx[["ml"]][1]]
      next
    }

#    print(param.table)
#    Print(MLEVARIAB, MLEPARAM, MLEMAX)

    MLEsettings(M)

#    Print(trafo, MLELB, MLEUB, fit$critical)
# Print(MLEtarget(c(0, 90, 0.6, 6.75))); hhhh
 

    
    MLEMAX <- -Inf ## must be before next "if (nMLEINDEX==0)"
    MLEVARIAB <- Inf ## nachfolgende for-schleife setzt MLEVARIAB
    MLEPARAM <- NA
    onborderline <- FALSE
    if (length(MLELB) == 0) { ## n.variab == 0
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

#        Print(i, param.table[IDX("variab"), i], MLEMAX, MLELB, MLEUB, LSQLB, LSQUB)

        if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {

          value <- MLEtarget(variab) ## !

#        Print(i, variab, value, MLEVARIAB)
          
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
   
      fnscale <-
        if (is.null(fit.fnscale) || is.na(fit.fnscale[M]))
          -max(abs(max), 0.1) else fit.fnscale[M]
      #

      mle.optim.control <-
        c(opt.control, list(parscale=parscale, fnscale=fnscale))

#      print(param.table)
#      Print(methodprevto$mle,variab,mle.optim.control, parscale, MLEVARIAB, MLELB, MLEUB)
#      Print(parscale, MLEVARIAB)
      stopifnot(length(parscale)==0 || length(parscale) == length(MLEVARIAB))
      
      MLEINF <- FALSE

      if (fit$critical < 2) {

#        Print("MLE", MLELB, MLEUB)
        OPTIM(MLEVARIAB, MLEtarget, lower = MLELB, upper=MLEUB,
              control=mle.optim.control, optimiser=optimiser, silent=silent)
        
        if (MLEINF) {
          if (printlevel>=PL_STRUCTURE )
            Print("MLEINF", MLEVARIAB, MLEMAX) else cat("#") #
 #         Print("MLE 2")
          OPTIM(MLEVARIAB, MLEtarget, lower = MLELB, upper=MLEUB,
                control=mle.optim.control,optimiser=optimiser,silent=silent)
          if (printlevel>=PL_STRUCTURE )
            Print("MLEINF new", MLEVARIAB, MLEMAX) #
        }
              
        options(show.error.messages = TRUE) ##
        mindistance <-
          pmax(fit$minbounddistance, fit$minboundreldist * abs(MLEVARIAB))
        
        onborderline <- 
          (abs(MLEVARIAB - MLELB) <
           pmax(mindistance,  ## absolute difference
                fit$minboundreldist * abs(MLELB) ## relative difference
                )) |
        (abs(MLEVARIAB - MLEUB) <
         pmax(mindistance, fit$minboundreldist * abs(MLEUB)))
      }
    } # length(MLELB) > 0
  
    if (printlevel>=PL_STRUCTURE )
      Print("mle first round", MLEVARIAB, MLEPARAM, MLEMAX) #
    
    if (!is.finite(MLEMAX)) {
      if (printlevel>=PL_IMPORTANT ) message(M, ": MLEtarget I failed.")
      param.table[[M]] <- MLEPARAM <- NaN
      variab <- MLELB ## to call for onborderline
      ml.residuals <- NA
    } else {
      param.table[[M]][tblidx[[M]][1]] <- MLEMAX
      param.table[[M]][IDX("variab")] <- MLEVARIAB
      
      param.table[[M]][IDX("param")] <- MLEPARAM
      param.table[[M]][IDX("covariab")] <- get.covariates(MLEVARIAB)     
      ml.residuals <- ML.RESIDUALS

      if (FALSE) Print(recall, fit$critical>=2 || any(onborderline),#
             !fit$only_users , fit$critical >= 0)
      
      if ((fit$critical>=2  || any(onborderline))
          && !fit$only_users && fit$critical >= 0) {
        ## if the MLE result is close to the border, it usually means that
        ## the algorithm has failed, especially because of a bad starting
        ## value (least squares do not always give a good starting point,helas)
        ## so the brutal method:
        ## calculate the MLE values on a grid and start the optimization with
        ## the best grid point. Again, there is the believe that the
        ## least square give at least a hint what a good grid is
        
        if (fit$critical == 0) {
          MLEgridmin <- gridmin
          MLEgridmax <- gridmax
          
          if (any(is.na(MLEgridmin)) || any(is.na(MLEgridmax))) {
            if (printlevel >= PL_SUBIMPORTANT ) {
              Print(cbind(MLELB, variab, MLEUB, onborderline), #
                    MLEgridmin, MLEgridmax)
            }
            warning(paste(M, "converged to a boundary value -- ",
                          "better performance might be obtained",
                          "when allowing for more lsq.methods"))
          } else {
            if (printlevel>=PL_FCTN_SUBDETAILS)
              show(1, M, MLEMAX, MLEVARIAB) else cat(detailpch)
            MLEgridlength <-
              max(3, round(fit$approximate_functioncalls^(1/n.variab)))
            ## grid is given by the extremes of the LS results
            ## so, therefore we should examine above at least 4 different sets
            ## of weights wichtig: gridmin/max basiert auf den reduzierten
            ## Variablen 
            step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-2) # grid starts
                                        # bit outside
            MLEgridmin <- pmax(MLEgridmin - runif(length(MLEgridmin)) * step/2,
                               MLELB)   # the extremes of LS
            MLEgridmax <- pmin(MLEgridmax + runif(length(MLEgridmax)) * step/2,
                               MLEUB)
            step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-1)
            
            startingvalues <- vector("list", length(step))
            for (i in 1:length(step)) {
              startingvalues[[i]] <-
                MLEgridmin[i] + step[i] * 0:(MLEgridlength-1)
            }
            
            startingvalues <- do.call(base::expand.grid, startingvalues)
            
            limit <- 10 * fit$approximate_functioncalls
            if ((rn <- nrow(startingvalues)) > limit) {              
              if (printlevel>=PL_STRUCTURE)
                cat("using only a random subset of the", rn, "grid points")
              rand <- runif(rn)
              startingvalues <-
                startingvalues[rand < quantile(rand, limit / rn), , drop=FALSE]
              gc()
            }
            
            MLEMAX <- -Inf
            
            apply(startingvalues, 1,
                  function(x) {
                  try(MLEtarget(x), silent=silent)
                  ## try-result:
                  ##     if (!is.numeric(z) && abs(z)<1e10) cat(paste(x, sep=","), "\n") else cat(".\n")
                  ##  if (MLEINF) stop("stop")
                })
            
            if (printlevel>=PL_FCTN_SUBDETAILS)
              Print("mle grid search", MLEVARIAB, MLEPARAM, MLEMAX) #
            
            ## side effect:Maximum is in MLEMAX!
            ##                             and optimal parameter is in MLEVARIAB
            if (printlevel>=PL_FCTN_SUBDETAILS)
              show(2, M, MLEMAX, MLEVARIAB)
            
            cat(detailpch)
            options(show.error.messages = show.error.message) ##


             
#            Print("critical")
            OPTIM(MLEVARIAB, MLEtarget, lower = MLELB, upper = MLEUB,
                  control=mle.optim.control, optimiser=optimiser,
                  silent=silent)
            options(show.error.messages = TRUE) ##
            if (!is.finite(MLEMAX) &&(printlevel>=PL_IMPORTANT))
              message("MLtarget II failed.\n")
            ## do not check anymore whether there had been convergence or not.
            ## just take the best of the two strategies (initial value given by
            ## LS, initial value given by a grid), and be happy.
            if (printlevel>=PL_FCTN_SUBDETAILS )
              show(3, M, MLEMAX, MLEVARIAB) else 
            if (printlevel>=PL_STRUCTURE )
              Print("mle second round", MLEVARIAB, MLEPARAM, MLEMAX) #
            
            if (is.finite(MLEMAX) && MLEMAX > param.table[[M]][tblidx[[M]][1]]){
              param.table[[M]][tblidx[[M]][1]] <- MLEMAX
              param.table[[M]][IDX("variab")] <- MLEVARIAB
              param.table[[M]][IDX("param")] <- MLEPARAM
              param.table[[M]][IDX("covariab")] <- get.covariates(MLEVARIAB)
              ml.residuals <- ML.RESIDUALS
            }
          } # (is.na(MLEgridmin[1]))
        } else { # fit$critical > 0  

          critical <- ptype == CRITICALPARAM;
          if (fit$critical>=3 || !any(critical)) critical <- rep(TRUE, n.variab)
          ncrit <- as.integer(fit$n_crit^(1/sum(critical)))
          if (!is.finite(ncrit)) ncrit <- 2
          if (ncrit^sum(critical) > 100 && printlevel>=PL_IMPORTANT)
            message("The optimisation may last a pretty long time!")
          
          if (ncrit > 1 || fit$critical>=2) {
            if (!is.null(transform)) {
              #Print(transform)
              stop("if 'transform' is given, 'critical' must be '0'")
            }
            w <- 0.5
            lowlist <- upplist <- list()
            for (i in 1:n.variab) {
              if (critical[i]) {               
                newparam <- seq(if (SDVAR.IDX[i]) 0 else MLELB[i], MLEUB[i],
                                len=ncrit+1)
                lowlist[[i]] <- newparam[-length(newparam)]
                upplist[[i]] <- newparam[-1]
              } else {
                lowlist[[i]] <- MLELB[i]
                upplist[[i]] <- MLEUB[i]
              }
            }
            lowlist <- as.matrix(do.call(base::expand.grid, lowlist))
            upplist <- as.matrix(do.call(base::expand.grid, upplist))
            
            orig.MLEVARIAB <- MLEVARIAB
            b.idx <- is.finite(bounds)

            stopifnot(length(lowlist) > 1)
          
            for (i in 1:nrow(lowlist)) {
              cat(detailpch)

              lowerupper <- GetLowerUpper(lowlist[i, ], upplist[i, ],
                                          trafo, optim_var_elimination, sillbounded,
                                          var.idx, nugget.idx, sill, nonugget,
                                          varmin=varmin, varmax=varmax)


              new.lower.vector <- lowerupper$lower
              new.upper.vector <- lowerupper$upper
              new.bounds <- TRUE
                        
              new.parscale <- guess <- fnscale <- NULL              

              model2 <- model
              if (anyRandomEffects) {    
                stopifnot(model2[[1]] %in% ZF_SELECT)
                model2[[1]] <- ZF_SYMBOLS_PLUS
                model2$subnr <- NULL
              }
   
              first.passage <- TRUE
              while (TRUE) {
                
                
                #Print("\n\n\n", recall, i, lowerupper, lowlist, upplist, trafo, sillbounded, var.idx, nugget.idx)
                    
                if (new.bounds) {
                  new.bounds <- FALSE
                  .C("PutValuesAtNA", Reg, as.double(new.lower.vector),
                     PACKAGE="RandomFields")
                  new.lower <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE,
                                        spConform=FALSE, do.notreturnparam=TRUE)

                  .C("PutValuesAtNA", Reg, as.double(new.upper.vector),
                     PACKAGE="RandomFields")
                  new.upper <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE,
                                        spConform=FALSE, do.notreturnparam=TRUE)
                }
               
                if (printlevel>=PL_STRUCTURE) cat("\nrecall rffit...\n")

             
                res <-
                  rffit.gauss(model=model2,
                              x=if (!dist.given) lapply(neu, function(x) x$x),
                              T=if (!dist.given) lapply(neu, function(x) x$T),
                              grid=grid,
                              data= data, ###
                              lower=new.lower,
                              upper=new.upper,
                              bc_lambda=bc_lambda,
                              mle.methods="ml",
                              lsq.methods=lsq.methods,
                              users.guess = guess,
                              distances=if (dist.given) distances,
                                        # lsq.methods=lsq.methods,
                              dimensions = dimensions,
                              optim.control=c(opt.control,
                                fnscale=list(fnscale),
                                parscale=list(new.parscale)),
                              transform=transform,
                              recall = TRUE,
                              general.pch = if (pch == "") "" else ":",
                              general.practicalrange = general$practicalrange,
                              general.spConform = FALSE,
                              fit.critical = -1,
                              fit.split =FALSE)
                
 
                guess <- res$table[[M]][IDX("param")]
                ## scale_ratio <- abs(log(abs(guess[-delete.idx]/new.parscale)))
                #Print(guess, new.lower.vector, users.lower, new.upper.vector,
                 #     users.upper, trafo, delete.idx)

                stopifnot(length(users.lower) + length(delete.idx)
                          == length(new.lower.vector))
                
                if (!all(guess >= new.lower.vector & guess <= new.upper.vector)
                    || !all(new.lower.vector[-delete.idx] > users.lower &
                            new.upper.vector[-delete.idx] < users.upper)) {
                  new.lower.vector <- pmin(new.lower.vector, guess)
                  new.upper.vector <- pmax(new.upper.vector, guess)
                  new.lower.vector[-delete.idx] <-
                    pmax(new.lower.vector[-delete.idx], users.lower)
                  new.upper.vector[-delete.idx] <-
                    pmin(new.upper.vector[-delete.idx], users.upper)
                  guess <- pmax(new.lower.vector, pmin(new.upper.vector, guess))
                  new.bounds <- TRUE
                }

                
                ml.value <- res$table[[M]][tblidx[[M]][1]]

                if (is.finite(ml.value)) {
                  if (ml.value > param.table[[M]][tblidx[[M]][1]]) {
                    if (printlevel > PL_RECURSIVE && !first.passage)
                      cat("parscale: table improved by ",
                          ml.value -  param.table[[M]][tblidx[[M]][1]], "\n") 
                    param.table[[M]][tblidx[[M]][1]] <- MLEMAX <- ml.value
                    for (j in c("variab", "param", "covariab"))
                      param.table[[M]][IDX(j)] <- res$table[[M]][IDX(j)]
                    ml.residuals <- res$ml$residuals
                  }
                }
                
                if (first.passage) {
                  old.ml.value <- ml.value
                } else {
                  if (printlevel > PL_REC_DETAILS) {
                    if (ml.value > old.ml.value) {                      
                      cat("parscale: mle improved by", ml.value-old.ml.value,
                          "\n")
                    }
                  }
                  break;
                }
                            
                ## urspruengliches Abbruchkriterium, das nicht gut fkt:
                ##value.ratio <- abs(log (abs(ml.value) / abs(MLEMAX)))
                ##outside <- (scale_ratio > fit$scale_ratio)  & MLEVARIAB != 0
                ##outside <- outside & (!b.idx | new.parscale >
                ##                      exp(fit$scale_ratio) * 
                ##                      (w * abs(guess) + (1-w) * bounds))
                ## if (!any(outside) && value.ratio <= fit$scale_ratio) {
               ##   break
               ## }

                
                zero <- new.parscale == 0
                new.parscale[zero] <-
                  pmax(abs(lower[zero]), abs(upper[zero])) / 10                
                new.parscale[b.idx] <-
                  w * new.parscale[b.idx] + (1-w) * bounds[b.idx]
                #new.parscale <- abs(guess[-delete.idx])

                #Print(new.parscale)
                restable <- as.matrix(res$table)
                used <- which(!apply(is.na(restable), 2, all)[-1:-2]) # ohne auto,user
                names.tbl <- names(tblidx)
                start <- which(names.tbl %in% alllsqmeth)              
                start <- if (length(start) == 0) 0 else min(start) - 1
                if (start>0) names.tbl <- names.tbl[-1:-start]
                 
                fnscale <- numeric(length(used))                
                for (j in 1:length(used)) {
                  fnscale[j] <- abs(restable[tblidx[[start + used[j]]][1],
                                              2 + used[j]])
                  if (names.tbl[used[j]] == "ml")
                    fnscale[j] <- -max(0.1, fnscale[j])
                }
                #Print("crit", fnscale)

                .C("PutValuesAtNA", Reg, guess,PACKAGE="RandomFields",
                   NAOK=TRUE) # ok
                guess <- GetModel(register=Reg, modus=GETMODEL_DEL_MLE,
                                  spConform=FALSE, do.notreturnparam=TRUE)
                first.passage <- FALSE
              } # while true
            } # for 1:nrow(lowlist)            
          } # ncrit  > 1          
        } # fit$critical > 0
      }
   
      if (!fit$only_users && (fit$reoptimise || any(SDVAR.IDX)) &&
          fit$critical >= 0 && fit$critical <= 1){
        if (pch!="") cat("$")
                
        new.parscale <- abs(revariab <- param.table[[M]][IDX("variab")])
        new.parscale[new.parscale == 0] <- 1e-3
        fnscale <- -max(0.1, abs(MLEMAX <- param.table[[M]][tblidx[[M]][1]]))
       # Print(new.parscale, fnscale)
      
        old.control<- if (exists("mle.optim.control")) mle.optim.control else NA
        mle.optim.control <- c(opt.control,
                               list(parscale=new.parscale, fnscale=fnscale))
        low <- MLELB
        
        if (any(SDVAR.IDX)) low[SDVAR.IDX] <- pmax(0, users.lower[SDVAR.IDX])

#        Print("after critical")
        OPTIM(revariab, MLEtarget, lower = low, upper=MLEUB,
              control=mle.optim.control, optimiser=optimiser, silent=silent)
        
        if (MLEMAX >  param.table[[M]][tblidx[[M]][1]]) {   
          if (printlevel > PL_SUBIMPORTANT)
            Print(old.control, fnscale, revariab, #
                  param.table[[M]][IDX("param")],
                  MLEMAX, MLEVARIAB,  MLEPARAM, MLELB, MLEUB)
          param.table[[M]][tblidx[[M]][1]] <- MLEMAX
          param.table[[M]][IDX("variab")] <- MLEVARIAB         
          param.table[[M]][IDX("param")] <- MLEPARAM
          param.table[[M]][IDX("covariab")] <- get.covariates(MLEVARIAB)

           
        }   
      }   
    } # is.finite(MLEMAX)
  } ## M

## nugget, var, s, nu

  ## Print(minmax, MLEtarget(c(0, 90.3, 0.6, 2.6^2)))


######################################################################
###     calculate all target values for all optimal parameters     ###
######################################################################
  for (i in 1:length(allmethods)) if (!is.na(param.table[1, i])) {
#    Print("calc", i, alllsqmeth, allmlemeth);
    
    idxCovar <- IDX("covariab")
    p <- param.table[IDX("variab"), i]

    if (!lambdaest) {
      for (M in alllsqmeth) {
        cur <- param.table[tblidx[[M]][1], i]
        if (is.na(cur) && !is.nan(cur) && M %in% lsqMethods) {
          LSQsettings(M)
          param.table[tblidx[[M]][1], i] <- LStarget(p)
        }
      }
    } 

    for (M in allmlemeth) {
      cur <- param.table[tblidx[[M]][1], i]
      #Print(cur)
      if (is.na(cur[1]) && !is.nan(cur[1]) && M %in% mleMethods) {
        MLEsettings(M)
        param.table[tblidx[[M]][1], i] <- MLEtarget(p)
      }   
      
      if (allmethods[i] == M) {    
        param.table[[M]][IDX("sdvariab")] <-
          INVDIAGHESS(param.table[[M]][IDX("variab")], MLEtarget,
                      control=mle.optim.control)
      }
    }
    
    for (M in allcrossmeth) {
      cur <- param.table[tblidx[[M]][1], i]
      if (is.na(cur) && !is.nan(cur) && M %in% NULL) { # crossMethods) {
        stop("not programmed ")
       ##  crosssettings(M) ## uncomment
        variab <- p
        if (nCoVarAll > 0) {
          variab <- c(variab, param.table[idxCovar, i])
        }
        param.table[tblidx[[M]][1], i] <- stop("") # crosstarget(variab)
      }
    }
  }
  if (pch!="" && !recall) cat("\n")



#  print(param.table)
#  Print(MLEtarget(c(0, 90.3, 0.6, 2.6^2)))
#  Print(MLEtarget(c(0.727, 39, 2.02, 6.5)), GetModel(Reg, modus=0))
  
  


######################################################################
###                     error checks                               ###
######################################################################
  ## if rather close to nugget and nugget==fixed, do exception handling.
  ## This part is active only if
  ## scale_max_relative.factor < lowerbound_scale_ls_factor
  ## By default it is not, i.e., the following if-condition
  ## will "always" be FALSE.
  
  if (optim_var_elimination == "yes" && length(nugget.idx) == 0 && any(SCALE.IDX)) {
    idx <- IDX("variab")
    alllsqscales <- param.table[idx, cm$lsq][SCALE.IDX, ]
   
    if (any(alllsqscales < mindistances/fit$scale_max_relative_factor,
            na.rm=TRUE))
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
  ## currently natural scaling only for optim_var_elimination...

  if (model[[1]] %in% ZF_SELECT) {
    model[[1]] <- ZF_SYMBOLS_PLUS
    model$subnr <- NULL
    info.cov <- .Call("SetAndGetModelInfo", Reg, list("Cov", model),
                      spdim, is.dist.vector, !is.dist.vector, time, xdimOZ,
                      fit$short, FALSE, TRUE,
                      PACKAGE="RandomFields")
  }

  
  lowerupper <- GetLowerUpper(lower, upper,
                              trafo, optim_var_elimination, sillbounded,
                              var.idx, nugget.idx, sill, nonugget,
                              varmin=varmin, varmax=varmax )
  lower <- lowerupper$lower
  upper <- lowerupper$upper
      idx <- lower == upper
  lower[idx] = upper[idx] = NA
  .C("PutValuesAtNA", Reg, as.double(lower), PACKAGE="RandomFields",
     NAOK=TRUE) # ok
  lower <- GetModel(register=Reg, modus=GETMODEL_SOLVE_MLE, replace.select=TRUE)
  .C("PutValuesAtNA", Reg, as.double(upper), PACKAGE="RandomFields",
     NAOK=TRUE) # ok
  upper <- GetModel(register=Reg, modus=GETMODEL_SOLVE_MLE, replace.select=TRUE)



  idx <- IDX("param") # henceforth

  ## die nachfolgenden Zeilen aendern die Tabellenwerte -- somit
  ## muessen diese Zeilen unmittelbar vor dem return stehen !!!
  if (fit$use_naturalscaling && any(SCALE.IDX)) {
    for (i in 1:length(allmethods)) {
      M <- allmethods[i]
      if (!is.na(param.table[[M]][1])) {
        param <- as.double(param.table[[M]][idx] + 0.0)
        .C("PutValuesAtNA", Reg, param, PACKAGE="RandomFields")
        .C("expliciteDollarMLE", Reg, param, PACKAGE="RandomFields")
        param.table[[M]][idx] <- param
      }
    }
  }

  
  idxCovar <- IDX("covariab")
  idx.meth <- rep(FALSE, length(allmethods))
  res <- values.res <- list()
  nparam <- as.integer(ncovparam + nCoVarAll)
  storage.mode(sum.not.isna.data) <- "integer"
  AICconst <-
    2 * nparam + 2 * nparam * (nparam + 1) / (sum.not.isna.data - nparam - 1)
  
  for (i in 1:length(allmethods)) {
    M <- allmethods[i]
    if (idx.meth[i] <- !is.na(param.table[1, i]) || is.nan(param.table[1,i])
        || length(idx)==1 && idx==0) {
      .C("PutValuesAtNA", Reg, as.double(param.table[idx, i]),
         PACKAGE="RandomFields")
      # modus: 1 : Modell wie gespeichert
      #        0 : Modell unter Annahme practicalrange>0
      #        2 : natscale soweit wie moeglich zusammengezogen
      modelres <- GetModel(register=Reg, modus=GETMODEL_SOLVE_MLE,
                           replace.select=TRUE)

     ## Print(M, modelres)
  
      ## trend werte in das Modell einfuegen fuer die fixed effects
      ## geht nur gut wenn sets==1
      if (length(idxCovar)>0) {
        regression <- param.table[idxCovar, i]
        #names(regression) <- betanames
        if (sets==1) {       
          if (keinplus <- !(modelres[[1]] %in% ZF_PLUSSELECT)) {
            stopifnot(length(effect) == 1)
            modelres <- list(ZF_SYMBOLS_PLUS, modelres)
          }
          i <- 1
        
          for (k in mmcomponents) {
            
            if (effect[k]==FixedTrendEffect) {
              nr <- startCoVar[i, k]
              if (length(modelres[[k+1]]$mean) > 0) {
                modelres[[k+1]]$mean=regression[nr]
                nr <- nr + 1
              }
              if (length(modelres[[k+1]]$plane) > 0) {
                modelres[[k+1]]$plane=regression[nr:(nr - 1 + tsdim)]
                nr <- nr + tsdim
              }
              if (length(modelres[[k+1]]$polydeg) > 0) {
                stop("not programmed yet")
              }
              if (length(modelres[[k+1]]$arbitraryfct) > 0) {
                stop("not programmed yet")
              }       
            } else if (effect[k]==FixedEffect) {
              
              modelres[[k+1]]$beta=regression[startCoVar[i, k] : endCoVar[i, k]]
            }
          }
          if (keinplus) modelres <- modelres[[2]]
        }
      }

      if (M == "ml")  {
        residu <- werte
        for (i in 1:length(residu)) {
          for (m in 1:len.rep[i]) {
            residu[[i]][idx.na[[i]][[m]], idx.na.rep[[i]][[m]]] <-
              ml.residuals[[i]][[m]]
          }

#          Print(dimdata, dimdata[i, ], residu[[i]])

          dim(residu[[i]]) <- dimdata[i, ] ## passt das so??
#          d <- dim(residu[[i]])
#          dim(residu[[i]]) <- c(d[1], d[2] * d[3])
        }
      }


      ml.value <- param.table[[M]][tblidx[["ml"]][1]]
      AIC <- 2 * nparam  - 2 * ml.value
      AICc <- AICconst - 2 * ml.value
      BIC <- log(sum.not.isna.data) * nparam - 2 * ml.value      
      param.table[[M]][tblidx[["AIC"]][1]] <- AIC
      param.table[[M]][tblidx[["AICc"]][1]] <- AICc
      param.table[[M]][tblidx[["BIC"]][1]] <- BIC

      if (nCoVarAll > 0) {
        c.table <- rbind(param.table[[M]][IDX("covariab")], NA)
        dimnames(c.table) <- list(NULL, betanames) 
      } else c.table <- NULL
      v.table <- rbind(param.table[[M]][IDX("variab")],
                       param.table[[M]][IDX("sdvariab")],
                       param.table[[M]][IDX("lower")],
                       param.table[[M]][IDX("upper")]
                       )

   
      
      dimnames(v.table) <- list(c("value", "sd", "lower", "upper"), varnames)
      p.table <- rbind(param.table[[M]][IDX("param")], NA)
#      Print(paramnames, varnames)
      dimnames(p.table) <- list(c("value", "sd"), paramnames)

#
       
      res[[M]] <- list(model=modelres,
                       trend=if (length(idxCovar)>0) regression,
                       variab = v.table,
                       param = p.table,
                       covariab = c.table,
                       ml.value = ml.value,
                       AIC = AIC,
                       AICc = AICc,
                       BIC = BIC,
                       residuals = if (M == "ml") residu
                       )
       
      class(res[[M]]) <- "RM_modelFit"
      
    } else {
      ## else  res[[M]] <- NA
      ## -- none --
    }
  }
  options(save.options)

                                        #  names(res) <- names(param.table)

#  Print("here")



  ## Print("here E", GetModel(Reg, modus=GETMODEL_SOLVE_MLE))
  
 #  Print("here E1")



  ## next lines must be the very last call as standards are overwritten
#  Print("Here")
  if (is.null(transform) && !varnugNA && !globalsigma && !sillbounded) {
#    Print(allmlemeth)
    
    trafo <- function(x) x
    sillbounded <- varnugNA <- FALSE
    optim_var_elimination <- "never"
    solvesigma <- FALSE
    
    nvarWithoutbc <- length(IDX("param"))
    MLELB <- rep(-Inf, nvarWithoutbc)
    MLEUB <- rep(Inf, nvarWithoutbc)    
      
    for (i in 1:length(allmlemeth)) {
      M <- allmlemeth[i]
      p <- param.table[[M]][IDX("param")]
      if (!is.na(p[1]) && M %in% mleMethods) {        
        MLEsettings(M)

        if (abs(cur - MLEtarget(p)) > 1e-10 && printlevel >= 1)
          message("likelihoods are calculated in two different ways leading",
                  " to the values ", cur, " and ",  MLEtarget(p),
                  ", which are not the same. Please contact author.\n")
       
        res[[M]]$param[2, 1:ncovparam] <- inv <- INVDIAGHESS(p, MLEtarget)
        param.table[[M]][IDX("sdparam")] <- inv       
      }
    }
  }

  ##  Print("end Here", res)
  if (printlevel >= 1 && BEYOND > 100)
    cat("\nNote: There are ",
        if (BEYOND > 1000)"very strong" else
        if (BEYOND > 250) "strong" else "some",
        " indications that the model might be overparametrised\nor that the bounds for the variables are too wide. Try narrower lower and\nupper bounds for the variables in the latter case. One of the critical\nparameters is 'lowerbound_var_factor' whose value (", fit$lowerbound_var_factor, ") xmight be reduced.\n", sep="")
  

  nice_modelinfo <- function(minmax) {
    minmax <- minmax[, -4, drop=FALSE]
    minmax <- as.data.frame(minmax)
    minmax$type <- TYPEOF_PARAM_NAMES[minmax$type + 1]
    minmax
  }
 
  if (!general$spConform) {
    Res <- c(list(ev = ev,
                  table=param.table,
                  n.variab = n.variab,
                  n.param = ncovparam,
                  n.covariates = nCoVarAll,
                  boxcox = lambdaest,
                  deleted = delete.idx,
                  lowerbounds=lower,
                  upperbounds=upper,
                  transform = transform,
                  number.of.data= sum.not.isna.data,
                  modelinfo = nice_modelinfo(minmax),               
                  number.of.parameters = nparam,
                  p.proj = integer(0),
                  v.proj = integer(0),
                  x.proj = TRUE,
                  fixed = NULL,
                  true.tsdim = tsdim,
                  true.vdim = vdim,
                  report = "",
                  submodels = NULL,
                  vario = "'$vario' is defunctioned. Use '$ml' instead of '$value$ml'!"
                  ),
             res)
    class(Res) <-  "RF_fit"
    return(Res)
  } else {
    ev2 <- lapply(ev, FUN=list2RFempVariog)
    res2 <- lapply(res, FUN=list2RMmodelFit, isRFsp=isRFsp, coord=RFsp.coord,
                   gridTopology=gridTopology, data.RFparams=data.RFparams,
                   Z$neu[[1]]$T)
    if (is.null(transform)) transform <- list()

#    str(res, 3)
 #   str(res2, 3)
   # cccc
    
    # res2 <- res2["autostart"] ## nur fuer test-Zwecke
 
    
    do.call.par <- c(list(Class = "RFfit",
                          Z=Z,
                          ev=ev2, 
                          table=param.table,
                          n.variab = n.variab,
                          n.param = ncovparam,
                          n.covariates = nCoVarAll,
                          boxcox = lambdaest,
                          deleted = delete.idx,
                          lowerbounds=list2RMmodel(lower),
                          upperbounds=list2RMmodel(upper),
                          transform=transform,                          
                          coord.units=coord.units,
                          variab.units=variab.units,
                          number.of.data = sum.not.isna.data,
                          number.of.parameters = nparam,
                          modelinfo = nice_modelinfo(minmax),
                          p.proj = integer(0),
                          v.proj = integer(0),
                          x.proj = TRUE,
                          fixed = NULL,
                          true.tsdim = tsdim,
                          true.vdim = vdim,
                          report = "",
                          submodels = NULL
                          ),
                     res2)
    
    fit2 <- do.call(methods::new, args=do.call.par)
                                        #   Print("here E2")
    return(fit2)
  }
}
