### functions used in several situations, but do not really
### belong to getNset.R
### in particular internal getNset functions staying on the R level


rfSplitTrendAndCov <-
  function(model, spatialdim, xdimOZ, Time,  forcedplus=FALSE,
           cathegories=list(trend=DetTrendEffect : SpVarEffect),
           remaining.name="cov") {
  ## splittet additives modell in trend-Teil und cov-Teil auf

    info <- .Call("SetAndGetModelInfo",
                 MODEL_SPLIT, list("Dummy", model),
                 as.integer(spatialdim),
                 FALSE,
                 TRUE, ## in case it is a non-domain model
                 as.logical(Time), 
                 as.integer(xdimOZ), 
                 as.integer(254), FALSE, TRUE,
                 PACKAGE="RandomFields")

    
  effect <- info$effect
  
  submodels <- list()
  for (j in 1:(length(cathegories)+1)) submodels[[j]]<-list()
  names(submodels) <- c(names(cathegories), remaining.name)

  if ( !(model[[1]] %in% ZF_PLUS)) model<- list(ZF_SYMBOLS_PLUS, model)        

  for (i in 2:length(model)) {
    j<-1
    while(j<=length(cathegories) && !(effect[i-1] %in% cathegories[[j]])) {
      j <- j+1
    }
    submodels[[j]][[length(submodels[[j]])+1]] <- model[[i]]
  }

  for (j in 1:(length(cathegories)+1)) {
    if (length(submodels[[j]]) > 1 ||
        (forcedplus && length(submodels[[j]])>0)) { # geaendert 2.11.11
      submodels[[j]] <- c(ZF_SYMBOLS_PLUS,  submodels[[j]])
    } else if (length(submodels[[j]])==1) {
      submodels[[j]] <- submodels[[j]][[1]]
    }
  }
  submodels[sapply(submodels, length)==0] <- NULL
  return(submodels)
}

StandardSimpleModel <- function(model, tsdim, aniso=TRUE, addnugget=TRUE, ...) {    ## bei mixed model ist eine zweite Ebene von '+' vorhanden
  ## einfache mixed model koennen aber wie Simple behandelt werden
  ## die zweite '+'-Ebene wird hier aufgeloest:
  ##
  ## Nun koennten noch Namen fuer die Untermodelle vergeben worden sein.
  ## Wenn dann "unlist" aufgerufen wird, werden die namen der Ober- und
  ## Untermodelle zu einem gemeinsamen Namen verbunden, der dann nicht
  ## mehr als als submodel-Name erkannt wird. Also hilft nur alle
  ## (potentiellen) Namen zu loeschen:
  if (model[[1]] %in% ZF_PLUS) {
    submodels <- sapply(model[-1], function(x) x[[1]])
    plus <- submodels %in% ZF_PLUS
    if (any(plus)) {
      plus <- which(plus) + 1
      subs <- lapply(model[plus], function(x) x[-1])
      names(subs) <- NULL
      subs <- unlist(subs, recursive=FALSE)
      names(subs) <- NULL
      trend <- model[-plus]
      names(trend) <- NULL      
      model <- c(trend, subs)
      ## Print(model,trend,subs); xxx
    }
  } 

 # Print("SSM X", model)

  model <- PrepareModel2(model, ...)
  storage.mode(tsdim) <- "integer"
 # userdefined <- GetParameterModelUser(model)
  
  vdim <- InitModel(MODEL_USER, list("Dummy", model), tsdim, NAOK=TRUE) # ok
  if (!is.numeric(vdim) || any(vdim > 1)) return("model not univariate")
  model <- GetModel(register=MODEL_USER, GETMODEL_DEL_NATSC) # , do.notreturnparam=TRUE)
  .C("DeleteKey", MODEL_USER)

  s <- rfSplitTrendAndCov(model=model, spatialdim=tsdim, xdimOZ=tsdim,
                          Time=FALSE, forcedplus=TRUE,
                          cathegories=list(simple=c(PrimitiveModel:SimpleModel),
                            trend=DetTrendEffect:FixedTrendEffect),
                          remaining.name="cov")
    
  if (!is.null(s$cov) || length(s$simple)>3)
    return("model is too complex to be simple")
  if (length(s$simple) <= 1) return("model is too primitive")
  s <- s$simple

  for (i in 2:length(s)) {
    if (! (s[[i]][[1]] %in% DOLLAR)) s[[i]] <- list(DOLLAR[1], var=1, s[[i]])
  }
  if (length(s)==2 && !(s[[2]][[length(s[[2]])]][[1]] %in% ZF_NUGGET)
              && addnugget){
    ## assuming nugget is missing
    s[[3]] <- list(DOLLAR[1], var=0, list(ZF_NUGGET[2]))
  }

  if (length(s) == 3) {  ## [[1]] is '+'
    if (s[[2]][[length(s[[2]])]][[1]] %in% ZF_NUGGET) {
      s[2:3] <- s[3:2] ## [[1]] is '+'
    }
    if (s[[2]][[length(s[[2]])]][[1]] %in% ZF_NUGGET)
      return("nugget model given twice")
    if (!(s[[3]][[length(s[[3]])]][[1]] %in% ZF_NUGGET))
      return("simple model must the sum of some domain and isotropic model and a nugget effect")
     if (length(c(s[[3]]$scale, s[[3]]$Aniso, s[[3]]$aniso, s[[3]]$proj)) > 0)
       return("nugget contains scale and/or anisotropic elements")
  }
  if (length(c(if (!aniso) s[[2]]$aniso, s[[2]]$Aniso, s[[2]]$proj)) > 0)
    return("model is a complicated scale structure")
       
  return(s)
}



earth_coordinate_names<- function(names) {
  ## Earth coordinates + possibly radius
  n <- substr(tolower(names), 1, 6)
  nc <- nchar(n)
  lon <- lat <- integer(length(n))
  for (i in 1:length(n)) {
    lon[i] <- substr("longit", 1, nc[i]) == n[i]
    lat[i] <- substr("latitu", 1, nc[i]) == n[i]
  }
  lonORlat <- lon | lat  
  earth <- all(nc[lonORlat] >= 2) && sum(lon==1) && sum(lat == 1)
  return(if (length(names)==2 | !earth) earth else         
         if ( (lo <- which(lon)) < (la <- which(lat)))
         which(lonORlat[]))
}

cartesian_coordinate_names <- function(names) {
  n <- substr(tolower(names), 1, 1)
  coords <-  c("T", "x", "y", "z")
  Txyz <- outer(n, coords, "==")
  cs <- colSums(Txyz)
  if (any(cs > 1) || sum(cs[1:2]) == 0 || any(diff(cs[-1]) > 0))
    return (integer(0))
  Txyz <- Txyz[, c(2:4, 1), drop=FALSE]
  ord <- apply(Txyz, 2, function(x) which(x > 0))
  ord <- order(unlist(ord))
  rs <- which(rowSums(Txyz) > 0)
  return(rs[ord])
}



data.columns <- function(data, xdim=0, force=FALSE, halt=TRUE) {
  #  Print("data.col", data)
  if (xdim>0 && xdim >= ncol(data)) stop("not enough columns in 'data'.")
  RFopt <- RFoptions()
  info <- RFoptions()$coords
  cn <- colnames(data)

  if (all(is.na(info$varnames))) {
    if (all(is.na(info$coordnames))) {
      if (is.null(cn)) {
        if (force) return(list(data=(xdim+1):ncol(data), x=1:xdim))
        if (halt)
          stop('colnames of data argument must contain "data" or "variable"')
        else return(NULL);
      }
      is.data <- (tolower(substr(cn, 1, 4)) == "data" |
                  tolower(substr(cn, 1, 4)) == "value" |
                  tolower(substr(cn, 1, 8)) == "variable")
      if (!any(is.data)) {
        if (force) return(list(data=(xdim+1):ncol(data), x=1:xdim))
        if (halt) stop('no colname starts with "data" or "variable"')
        else return(NULL);
      }
      is.data <- which(is.data)
      if (is.data[1] > 1) is.data <- is.data[1] : ncol(data)# coord am Anfang
      ##     dann wird Rest alles als data angenommen, egal welcher Name
    }
  } else {
    if (is.numeric(info$varnames)) {
      is.data <- rep(info$varnames, length.out=2)
      if (is.na(is.data[1])) is.data[1] <- 1
      if (is.na(is.data[2])) is.data[2] <- ncol(data)
      is.data <- is.data[1] : is.data[2]
    } else {
      if (RFopt$general$vdim_close_together)
        stop("'vdim_close_together' must be FALSE")
      l <- list()
      vdim <- length(info$varnames)
      for (v in 1:vdim)
        l[[v]] <-
          substring(cn, 1, nchar(info$varnames[v])) == info$varnames[v]
      repet <- sapply(l, sum)
      if (repet[1] == 0) stop("data names could not be detected") 
      if (any(repet != repet[1]))
        stop("detected repetitions are not equal for all components")
      m <- matrix(unlist(l), ncol=vdim)
      if (any(rowSums(m) > 1))
        stop("names of multivariate components not unique")
      is.data <- as.vector(t(apply(m, 2, which)))
    }
  }

  if (all(is.na(info$coordnames))) {
    is.x <- (1:ncol(data))[-is.data]
    if (xdim > 0) {
      if (length(is.x) < xdim)
        stop("not enough columns for coordinates found ")
      if (xdim < length(is.x) &&
          RFopt$general$printLevel >= PL_SUBIMPORTANT)
        message("column(s) '", paste(is.x[-1:-xdim], collapse=","),
                "' not used.\n")
      is.x <- is.x[1:xdim]
    }
  } else {
    if (is.numeric(info$coordnames)) {
      is.x <- rep(info$coordnames, length.out=2)
      if (is.na(is.x[1])) is.x[1] <- 1
      if (is.na(is.x[2])) is.x[2] <- ncol(data)
      is.x <- is.x[1] : is.x[2]
    } else {
      l <- list()
      len <- length(info$coordnames)
      for (i in 1:len)
        l[[v]] <-
          substring(cn, 1, nchar(info$coordnames[v])) == info$coordnames[v]
      is.x <- unlist(l)
      if (xdim > 0 && xdim != length(l))
        stop("expected dimension of coordinates does not match the found coordinates")
    }
     
    if (all(is.na(info$varnames))) {
      is.data <-  (1:ncol(data))[-is.x]
      if (length(is.data) == 0) stop("no columns for data found")
    } else {
     if (any(is.x %in% is.data))
       stop("column names and data names overlap.")
     if (length(is.x) + length(is.data) < ncol(data) &&
         RFopt$general$printLevel >= PL_SUBIMPORTANT)
       message("column(s) '",
               paste(1:ncol[c(-is.x, -is.data)], collapse=","),
               "' not used.\n")
    }
  }
  return(list(data=is.data, x=is.x) )
}

    
GetDataNames <- function(model, coords, locinfo, data) {
  if (!missing(data) && length(data) > 0) {
    cd <- try(CheckData(model=model, given=coords, data=data))      
    if (class(cd) != "try-error")
      return(list(coord.names=cd$coord.names, variab.names=cd$variab.names))
  }

  variab.names <- extractVarNames(model)
  coord.names <-
    if (!missing(coords) && is.matrix(coords)) colnames(coords) else NULL
                       # gets response part of model, if model is a formula

  Zeit <- locinfo$Zeit
  ts.xdim <- locinfo$spatialdim + locinfo$Zeit
  if (is.null(coord.names))
    coord.names <- paste("coords.x", 1:ts.xdim, sep="")
  if (Zeit) coord.names[ts.xdim] <- "coords.T1"

  return(list(coord.names=coord.names, variab.names=variab.names))
}



search.model.name <- function(cov, name, level) {
  if (length(name) == 0 || length(cov) ==0) return(cov);
  if (!is.na(pmatch(name[1], cov))) return(search.model.name(cov, name[-1], 1))

  for (i in 1:length(cov$submodels)) {
    found <- search.model.name(cov$submodels[[i]], name, 1)
    if (!is.null(found)) return(found)      
  }
  found <- search.model.name(cov$internal, name, 1)
  if (!is.null(found)) return(found)
  if (level == 0) stop("model name not found")
  return(NULL)
}




## to do: grid
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

  d <- dim(all$given$x)
  ts.xdim <- as.integer(d[1])
  n <- as.integer(d[2])

  if (model.nr >= 0) {
    natsc <- .C("MultiDimRange", as.integer(model.nr), natsc = double(ts.xdim),
                PACKAGE="RandomFields")$natsc
  } else natsc <- rep(1, ts.xdim)

 
  u <- numeric(ts.xdim)
  for (i in 1:ts.xdim) {
    u[i] <- length(unique(all$given$x[i,]))
  }

  if (all$given$grid) {
    stop("not programmed yet")
  } else {
    Range <- apply(all$given$x, 1, range)
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
      coord.idx <- floor((all$given$x - Range[1,]) / (len / parts))
      
      cumparts <- cumprod(parts)
      totparts <- as.integer(cumparts[length(cumparts)])
      Ccumparts <- as.integer(c(1, cumparts))
      cumparts <- Ccumparts[-length(Ccumparts)]
      
      ## elms.in.boxes <- integer(totparts)  
      ## neighbours <- integer(totparts)
      cumidx <- as.integer(colSums(coord.idx * cumparts))
      
      elms.in.boxes <- .Call("countelements", cumidx, n, totparts,
                             PACKAGE="RandomFields")
      
      neighbours <- .Call("countneighbours", ts.xdim, parts, locfactor, Ccumparts,
                     elms.in.boxes, PACKAGE="RandomFields")
      
      ## if there too many points within all the neighbours, then split
      ## into smaller boxes
      zaehler <- zaehler * 2
      
      ## image(neighbours, zlim=c(0:(prod(parts)-1)))
      if (!is.null(neighbours)) break;
    }
    
    l <- list()
    l[[1]] <- .Call("getelements", cumidx, ts.xdim, n, Ccumparts,
                    elms.in.boxes, PACKAGE="RandomFields")
    l1len <- sapply(l[[1]], length)
  }

  

  if (length(all$x) > 0) {
    if (all$grid) {
     stop("not programmed yet")
    } else {
      ## now calculate the boxes for the locations where we will interpolate
      i <- pmax(0, pmin(parts-1, floor((t(all$x) - Range[1,]) / (len / parts))))
      dim(i) <- rev(dim(all$x))
      i <- as.integer(colSums(i * cumparts))
      
      res.in.boxes <- .C("countelements", i, nrow(all$x), totparts,
                         PACKAGE="RandomFields")
      
      notzeros <- res.in.boxes > 0
      l[[3]] <- .Call("getelements", i, ts.xdim, as.integer(nrow(all$x)),
                      Ccumparts, res.in.boxes, PACKAGE="RandomFields")[notzeros]
      ## TO DO : idx[[3]] passt nicht, da sowohl fuer Daten
      ##         als auch coordinaten verwendet wird. Bei repet > 1
      ##         ist da ein Problem -- ueberpruefen ob repet=1
   

    }
  } else {
    notzeros <- TRUE
  }

  if (all$given$grid) {
    stop("not programmed yet")
  } else {
    ll <- .Call("getneighbours", ts.xdim, parts, locfactor, Ccumparts,
                neighbours, PACKAGE="RandomFields")[notzeros]
    less <- sapply(ll, function(x) sum(elms.in.boxes[x]) < minimum) | !shared
    ##                  if !shared then all(less)==TRUE

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
  }
  return(if (shared) l else lapply(l[[2]], function(x) unlist(l[[1]][x])))
}



FinishImputing <- function(data, all, n, spConform) {
  ## to do: grid
  
  #Print(data, all, tail(all$simu), spConform);
  if (is(data, "RFsp")) {
    if (spConform) {
      data@data[ , ] <- as.vector(all$simu)
      return(data)
    } else {
      values <- as.matrix(data@data)
      values[is.na(values)] <- all$simu
      
      #print(cbind(coordinates(data), values))
      return(cbind(coordinates(data), values))
    }
  } else { ## not RFsp
    #Print("for testing")
    coords <- all$fullgiven$x
    if (all$fullgiven$grid) {
      ## to do
      stop("not programmed yet")
    } else {
      ##  coords <- all$x
      colnames(coords) <- all$coord.names
      
      values <- data[, all$data.col]
      values[is.na(values)] <- all$simu
      
      if (!spConform)  return(cbind(coords, values))
      
      tmp.all <- conventional2RFspDataFrame(data=values, coords=coords,
                                            gridTopology=NULL,
                                            n=n, vdim=all$vdim,
                                            vdim_close_together=FALSE)
      all <- tmp.all@data         
      if (is(tmp.all, "RFspatialPointsDataFrame"))
        try(tmp.all <- as(tmp.all, "RFspatialGridDataFrame"), silent=TRUE)
      if (is(tmp.all, "RFpointsDataFrame"))
        try(tmp.all <- as(tmp.all, "RFgridDataFrame"), silent=TRUE)
    }
    return(tmp.all)
  }
}
