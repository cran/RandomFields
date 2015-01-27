


PrepareModel2 <- function(model, ...) {
  if (missing(model) || is.null(model)) stop("'model' must be given.")
  m <- parseModel(model, ...)
  
  #if (!is.list(model))
      
#  if (class(try(all.equal(parseModel(list2RMmodel(m)), m, check.attr=FALSE)))=="try-error"){
#    print("AlexWarning in PrepareModel2: ruecktrafo via list2RMmodel waere falsch oder nich dasselbe ... bitte model mir zuschicken")
#    print("doppelt hin und zurueck:")
#    print(try(all.equal(parseModel(list2RMmodel(parseModel(list2RMmodel(m)))), parseModel(list2RMmodel(m)), check.attr=FALSE)))
#  }

  if (notplus <- !(m[[1]] %in% ZF_PLUS)) m <- list(ZF_SYMBOLS_PLUS, m)
     
  for (i in 2:length(m)) {
    if ((m[[i]][[1]] %in% ZF_MIXED) && length(m[[i]]$X)==1 &&
        is.numeric(m[[i]]$X) && m[[i]]$X==1 && !is.null(m[[i]]$b)) {        
      m[[i]] <- list(ZF_TREND[2], mean=m[[i]]$b)
      if (RFoptions()$general$printlevel > PL_IMPORTANT)
        message(paste("The '1' in the mixed model definition has been replaced by '", ZF_TREND[1], "(mean=", m[[i]]$mean, ")'.", sep=""))
      }
  }

  if (notplus) m <- m[[2]]
  class(m) <- "RM_model"
  return(m)

#  if (class(model) != "formula") {
#    if (is.list(model)) return(model)
#    else stop("model of unknown form -- maybe you have used an obsolete definition. See ?RMmodel for the model definition")
#  }
#  return(listmodel)
}


PrepareModel <-  function(model, param, trend=NULL, 
                          nugget.remove=TRUE, method=NULL) {
  ## any of the users model definition (standard, nested, list) for the
  ## covariance function is transformed into a standard format, used
  ## especially in the c programs
  ##
  ## overwrites in some situation the simulation method for nugget.
  ## allows trend to be NA (or any other non finite value  -- is not checked!)
  ## trend has not been implemented yet!

  if (is(model, ZF_MODEL))
    stop("models of class ZF_MODEL cannot be combined with obsolete RandomFields functions")

  if (!is.null(method)) stop("to give method in PrepareModel is obsolete")
  
  if (!is.null(trend))      
     if (!is.numeric(trend) || length(trend)!=1)
        stop("in the obsolete setting, only constant mean can used")
 
  if (is.list(model) && is.character(model[[1]]) &&
      (is.null(names(model)) || names(model)[[1]]=="")) {
     if (!missing(param) && !is.null(param))
        stop("param cannot be given in the extended definition")
     if (is.null(trend)) return(model)
     trend <- list(ZF_TREND[2], mean=trend)
     if (model[[1]] %in% ZF_PLUS) return(c(model, list(trend)))
     else return(list(ZF_SYMBOLS_PLUS, model, trend))
  }
    
  printlevel <- RFoptions()$general$printlevel
  STOP <- function(txt) {
    if (printlevel>=PL_ERRORS) {
      cat("model: ")
      if (!missing.model) Print(model) else cat(" missing.\n") #
      cat("param: ")
      if (!missing.param) Print(param) else cat(" missing.\n") #
      cat("trend: ")
      Print(trend) # 
     }
    stop("(in PrepareModel) ", txt, call.=FALSE)
  }
   
  transform <- function(model) {
    if (!is.list(model)) {
      STOP("some elements of the model definition are not lists")
    }
    m <- list(DOLLAR[1], var=model$v)
    lm <- length(model) - 3 # var, scale/aniso, name
    if (!is.null(model$a)) m$aniso <- model$a else m$scale <- model$scale
##    model <- c(model, if (!is.null(model$a))
##               list(aniso=model$a) else list(scale=model$s)) ## ???

    
    if (!is.na(p <- pmatch("meth", names(model), duplicates.ok=TRUE))) {
      if (printlevel>=PL_ERRORS)  Print(p, model) #
      stop("method cannot be given with the model anymore. It must be given as a parameter to the function. See 'RFoptions' and 'RFsimulate'")
     }
  
    if (!is.null(model$me))
      stop("'mean' seems to be given within the inner model definitions"); 
    if (!is.character(model$m)) {
       stop("'model' was not given extacly once each odd number of list entries or additional unused list elements are given.")
    }
    m1 <- list(model$m)
    if (!is.null(model$k)) {
      lm <- lm - 1
      if (length(model$k) != 0)
        for (i in 1:length(model$k)) {
          eval(parse(text=paste("m1$k", i, " <- model$k[", i, "]", sep="")))
      }
    }
    if (lm != 0) {
      if (printlevel>=PL_ERRORS) Print(lm, model) #
      stop("some parameters do not fit")
    }
    m <- c(m, list(m1))

    return(m)
    
  } # end transform

  op.list <- c(ZF_SYMBOLS_PLUS, ZF_SYMBOLS_MULT)  ## if others use complex list definition !
  missing.model <- missing(model)
  missing.param <- missing(param) || is.null(param)

  if (full.model <- missing.param && is.null(model$param)) { ## full model
    if (RFoptions()$internal$warn_oldstyle)
      warning("the sequential list format is depreciated.")
    if (missing.model || (length(model)==0)) model <- list()
    else if (!is.list(model))
      STOP("if param is missing, model must be a list of lists (or a list in the extended notation)")
    if (is.null(trend) + is.null(model$mean) + is.null(model$trend)<2)
      STOP("trend/mean is given twice")
    if (!is.null(model$mean)) trend <- model$mean else
    if (!is.null(model$trend)) trend <- model$trend else trend <- NULL

    model$trend <- model$mean <- NULL
    ## the definition might be given at a deeper level as element
    ## $model of the list:
    if (is.list(model$model)) {
      if (!is.list(model$model[[1]]))
        STOP("if param is missing, the model$model must be a list of lists")
      model <- model$model
    }
    if (length(model)==0) { ## deterministic      
      return(if (is.null(trend)) NULL else list(ZF_TREND[2], mean=trend))
    }
    if (length(model) %% 2 !=1) STOP("list for model definition should be odd")
    if (length(model)==1)
      return(if (is.null(trend) ||
                 is.numeric(trend) && length(trend)==1 && !is.na(trend)&&trend==0)
             transform(model[[1]]) 
             else list(ZF_SYMBOLS_PLUS, transform(model[[1]]),
                       list(ZF_TREND[2], mean=trend)));

    op <- pmatch(c(model[seq(2, length(model), 2)], recursive=TRUE),
                 op.list, duplicates.ok=TRUE) - 1
    if (!all(is.finite(op))) STOP("operators are not all allowed; see the extended list definition for extensions")
    model <- model[seq(1, length(model), 2)]

    plus <- which(op==0)
    if (length(plus) == 0) {
      m <- list("*", lapply(model, transform))
    } else {
      plus <- c(0, plus, length(op)+1)
      m <- list(ZF_SYMBOLS_PLUS)
      for (i in 1:(length(plus) - 1)) {
        m[[i+1]] <-
          if (plus[i] + 1 == plus[i+1]) transform(model[[plus[i] + 1]])
          else list(ZF_SYMBOLS_MULT,
                    lapply(model[(plus[i] + 1) : plus[i+1]], transform))
      }
    }
   model <- m
  } else { ## standard definition or nested model
    if (missing.param) { ## a simple list of the model and the
      ##                    parameters is also possible
      if (is.null(param <- model$p)) STOP("is.null(model$param)")
      stopifnot(is.null(trend) || is.null(model$trend))
      if (is.null(trend)) trend <- model$trend
      if (!is.null(model$mean)) {
        if (!is.null(trend)) STOP("mean and trend given twice")
        trend <- model$mean
      }
      model <- model$model
    }
    stopifnot(is.character(model), length(model)==1)
    if (is.matrix(param)) { ## nested
      if (nrow(param) == 1)
        return(PrepareModel(model=model, param=c(param[1], 0, param[-1]),
                            trend=trend))
      name <- model
      model <- list(ZF_SYMBOLS_PLUS)#, method=method)
      for (i in 1:nrow(param)) {
        model <- c(model,
                   if (is.na(param[i, 2]) || param[i, 2] != 0)
                   list(list(DOLLAR[1], var=param[i, 1], scale=param[i, 2],
                             if (ncol(param) >2) list(name, k=param[i,-1:-2])
                             else list(name)))
                   else list(list(DOLLAR[1], var=param[i,1],
                                  list(ZF_NUGGET[2]))))
      }
    } else if (is.vector(param)) {  ## standard, simple way
      ## falls trend gegeben, dann ist param um 1 Komponente gekuerzt
      if (is.null(trend)) {
        trend <- param[1]
        param <- param[-1]
      } else message("It is assumed that no mean is given so that the first component of param is the variance")
      if (model == ZF_NUGGET[2]) {
        model <- transform(list(model=model, var=sum(param[1:2]), scale=1)) 
      } else {
        if  (length(param) > 3)
          model <- transform(list(model=model, var=param[1], scale=param[3],
                                  k=param[-1:-3]))
        else 
          model <- transform(list(model=model, var=param[1], scale=param[3]))
        if (is.na(param[2]) || param[2] != 0 || !nugget.remove) {# nugget
          model <- list(ZF_SYMBOLS_PLUS,
                        model,
                        transform(list(model=ZF_NUGGET[2], var=param[2], scale=1)))
        }
        ## if (!is.null(method)) model <- c(model, method=method) ## doppelt
      }
    } else stop("unknown format")  # end nested/standard definition
  }

  return(if (is.null(trend) ||
             is.numeric(trend) && length(trend)==1 &&  !is.na(trend) &&trend==0)
         return(model)
         else if (model[[1]] %in% ZF_PLUS)
                 c(model, list(list(ZF_TREND[2], mean=trend)))
         else list(ZF_SYMBOLS_PLUS, model, list(ZF_TREND[2], mean=trend)))
}


seq2grid <- function(x, name, grid, warn_ambiguous, gridtolerance) {
  xx <- matrix(nrow=3, ncol=length(x))
  step0 <- rep(FALSE, length(x))
  gridnotgiven <- missing(grid) || length(grid) == 0
  
  for (i in 1:length(x)) {
    if (length(x[[i]]) == 1) {
      xx[,i] <- c(x[[i]], 0, 1)
      next
    }
    step <- diff(x[[i]])
    if (step[1] == 0.0) {
      
      ok <- step0[i] <- all(step == 0.0)      
    } else {
      ok <- max(abs(step / step[1] - 1.0)) <= gridtolerance
    }

    if (!ok) {
      if (gridnotgiven) return(FALSE)
      if (!TRUE)
        Print(i, x[[i]][1:min(100, length(x[[i]]))], #
              step[1:min(100,length(step))],
              range(diff(step[1:min(100,length(step))])))
      stop("Different grid distances detected, but the grid must ",
           "have equal distances in each direction -- if sure that ",
           "it is a grid, increase the value of 'gridtolerance' which equals ",
           gridtolerance,".\n")
    }

    xx[,i] <- c(x[[i]][1], step[1], if (step0[i]) 1 else length(x[[i]]))
  }

  if (FALSE && gridnotgiven && warn_ambiguous && length(x) > 1) {
    RFoptions(internal.warn_ambiguous = FALSE)
    message("Ambiguous interpretation of coordinates. Better give 'grid=TRUE' explicitly. (This message appears only once per session.)")
  }

  if (any(step0)) {
    if (all(step0)) {
      if (gridnotgiven) return(FALSE)
      else stop("Within a grid, the coordinates must be distinguishable")
    } else {
      if (gridnotgiven && warn_ambiguous) {
        RFoptions(internal.warn_ambiguous = FALSE)
        warning("Interpretation as degenerated grid. Better give 'grid' explicitely. (This warning appears only once per session.)")
      }
    }
  }

  return(xx)
}

CheckXT <- function(x, y=NULL, z=NULL, T=NULL, grid, distances=NULL,
                    dim=NULL, # == spdim!
                    length.data,
                    y.ok = FALSE, 
                    printlevel = RFoptions()$general$printlevel){

  ## converts the given coordinates into standard formats
  ## (one for arbitrarily given locations and one for grid points)
  
  RFopt <- RFoptions()
  curunits <- RFopt$coords$coordunits
  newunits <-  RFopt$coords$new_coordunits
  coord_system <-  RFopt$coords$coord_system
  new_coord_system <-  RFopt$coords$new_coord_system
  
  if (!missing(distances) && length(distances) > 0) {
    stopifnot(is.matrix(distances) || (!missing(dim) && !is.null(dim)),
              (missing(grid) || length(grid) == 0),
              missing(x) || is.null(x),
              is.null(y),
              is.null(z),
              is.null(T))
    
    if (coord_system != new_coord_system && new_coord_system != "keep")
      stop("coordinate systems differ")

    if (class(distances) == "dist") {
      x <- as.vector(distances)
      l <- length(x)
    } else if (is.matrix(distances) || is.vector(distances)) {
      if (is.matrix(distances)) {        
        len <- nrow(distances)
        if (is.null(dim)) dim = ncol(distances)
        else if (dim != ncol(distances))
          stop("matrix of distances does not fit the given dimension")
      } else {
        len <- length(distances)
        if (is.null(dim))
          stop("dim is not given although 'distances' are used")
      }
      x <- distances
      l <- as.integer(1e-9 + 0.5 * (1 + sqrt(1 + 8 * len)))
      if (l * (l-1) / 2 != len) l <- NaN
    } else {
      stop("'distances' not of required format.")
    }

    return(list(x=as.matrix(x),
                T=NULL, 
                Zeit=FALSE,
                restotal = l,
                l = l,
                spatialdim=dim,
                grid = FALSE,
                distances = TRUE,
                coordunits = curunits,
                new_coordunits = newunits
                )
           )
  }

  stopifnot(!missing(x))
  if (is(x, "RFsp") || isSpObj(x)) {
    return(CheckXT(x=coordinates(x), y=y, z=z, T=T, grid=grid,
                   distances=distances, dim=dim, length.data=length.data,
                   y.ok=y.ok, printlevel=printlevel))
  }    
  if (is.raster(x)) x <- as(x, 'GridTopology')
 
  if ((missing(grid) || length(grid) == 0) && !missing(length.data)) {
    new <-  try(CheckXT(x=x, y=y, z=z, T=T, grid=TRUE, distances=distances,
                        dim=if (!missing(dim)) dim,
                        length.data = length.data, y.ok =y.ok,
                        printlevel = printlevel
                        ), silent=TRUE)
    if (grid <- (class(new) != "try-error")) {
      ratio <- length.data / new$restotal

      if (grid <- ratio == as.integer(ratio)) {
        if (printlevel>=PL_IMPORTANT && new$spatialdim > 1)
          message("Grid detected. If it is not a grid, set grid=FALSE.\n")
      }
    }
    return(if (grid) new else {
      CheckXT(x, y, z, T, grid=FALSE, distances,
              if (!missing(distances) && length(distances) > 0) dim=1,
              length.data = length.data,
              printlevel = printlevel) }
           )
  } # if (missing(grid) && !missing(length.data))


  gridtriple <- FALSE

  if (is.GridTopology <- is(x, "GridTopology")){
    x <- rbind(x@cellcentre.offset,
               x@cellsize,
               x@cells.dim)
    if ((missing(grid) || length(grid) == 0)) grid <- TRUE else stopifnot(grid)
    gridtriple <- TRUE
  }
  ##else {
  ##  is.GridTopology <- FALSE
  ##}

  
  if (is.data.frame(x)) {
    if (ncol(x)==1) x <- as.vector(x) else x <- as.matrix(x)
  }
  
  stopifnot(length(x) != 0)
#  stopifnot(all(unlist(lapply(as.list(x), FUN=function(li) is.numeric(li))))) ## wann benoetigt???

  stopifnot(is.numeric(x))# um RFsimulte(model, data) statt data=data abzufangen
  
  
#  stopifnot(all(is.finite(x)), all(is.finite(y)), all(is.finite(z))) ; s.u. unlist
 
     
  if (is.matrix(x)) {
    if (!is.numeric(x)) stop("x is not numeric.")
    if (!is.null(z)) stop("If x is a matrix, then z may not be given")
    if (!is.null(y)) {
      if (!y.ok) stop("If x is a matrix, then y may not be given")
      if (!is.null(T))
        stop("If x is a matrix and y is given, then T may not be given")
      if (!is.matrix(y) || ncol(y) != ncol(x) ||
          nrow(x)==3 && nrow(y)!=3 && ((missing(grid) || length(grid) == 0) ||
                                grid))
        stop("y does not match x (it must be a matrix)")
    }

    if (coord_system == COORD_SYS_NAMES[coord_auto + 1] && ncol(x) >= 2
        && ncol(x) <= 3 && !is.null(n <- dimnames(x)[[2]])) {
      if (any(idx <- earth_coordinate_names(n))) {
        if (length(idx) == 2 && !all(idx == 1:2))
          stop("earth coordinates not in order longitude/latitude")
        cur <- curunits[1]
        newunits <- RFopt$coords$new_coordunits
        curunits <- RFopt$coords$coordunits
        curunits[1:2] <- c("longitude", "latitude")
        if (newunits[1] == "") newunits[1] <-  UNITS_NAMES[units_km + 1]
        newunits[2:3] <- newunits[1]                
        if (RFopt$internal$warn_coordinates)
          message("\n\nNOTE: current units are ",
                  if (cur=="") "not given and" else paste("'", cur, "', but"),
                  " earth coordinates detected:\n",
                  "earth coordinates will be transformed into units of '",
                  newunits[1],
                  "'.\nIn particular, the values of all scale parameters of ",
                  "any model defined\nin R^3 (currently all models!) are ",
                  "understood in units of '", newunits[1],
                  "'.\nChange options 'coord_system' and/or 'units' if ",
                  "necessary.\n(This message appears only once per session.)\n")
        coord_system <- COORD_SYS_NAMES[earth + 1]
        RFoptions(coords.coord_system = coord_system,
                  coords.coordunits = curunits,
                  coords.new_coordunits = newunits,
                  internal.warn_coordinates=FALSE)

      } else {
         RFoptions(coords.coord_system =  COORD_SYS_NAMES[cartesian + 1])
      }
    }    
   
    spatialdim <- ncol(x)
    len <- nrow(x)
    if (spatialdim==1 && len != 3 && (missing(grid) || length(grid) == 0)) {
      if (length(x) <= 2) grid <- TRUE
      else {
        dx <- diff(x)
        grid <- max(abs(diff(dx))) < dx[1] * RFopt$general$gridtolerance
      }
    } # else {

    if ((missing(grid) || length(grid) == 0) &&
        any(apply(x, 2, function(z) (length(z) <= 2) || max(abs(diff(diff(z))))
                  > RFopt$general$gridtolerance))) {
      grid <- FALSE
    }

    if ((missing(grid) || length(grid) == 0) || !is.logical(grid)) {
      grid <- TRUE
      if (spatialdim > 1 && RFopt$internal$warn_ambiguous) {
        RFoptions(internal.warn_ambiguous = FALSE)
        warning("Ambiguous interpretation of the coordinates. Better give the logical parameter 'grid=TRUE' explicitely. (This warning appears only once per session.)")
      }
    }

    if (grid && !is.GridTopology) {
      if (gridtriple <- len==3) {
        if (printlevel >= PL_SUBIMPORTANT && RFopt$internal$warn_oldstyle) {
          message("x was interpreted as a gridtriple; the new gridtriple notation is:\n  1st row of x is interpreted as starting values of sequences,\n  2nd row as step,\n 3rd row as number of points (i.e. length),\n  in each of the ", ncol(x), " directions.")
        } 
      } else len <- rep(len, times=spatialdim)   # Alex 8.10.2011
    }

    if (grid && !gridtriple) {
      ## list with columns as list elements -- easier way to
      ## do it??
      x <- lapply(apply(x, 2, list), function(r) r[[1]])
      if (!is.null(y)) y <- lapply(apply(y, 2, list), function(r) r[[1]])
    }
  } else { ## x, y, z given separately
    if (is.null(y) && !is.null(z)) stop("y is not given, but z")
    xyzT <- list(x=if (!missing(x)) x, y=y, z=z, T=T)
    for (i in 1:4) {
      if (!is.null(xyzT[[i]]) && !is.numeric(xyzT[[i]])) {
        if (printlevel>PL_IMPORTANT) 
          message(names(xyzT)[i],
                  " not being numeric it is converted to numeric")
        assign(names(xyzT)[i], as.numeric(xyzT[[i]]))
      }
    }
    remove(xyzT)
    spatialdim <- 1 + (!is.null(y)) + (!is.null(z))
    if (spatialdim==1 && ((missing(grid) || length(grid) == 0) || !grid)) {
      ## ueberschreibt Einstellung des Nutzers im Falle d=1
      if (length(x) <= 2) newgrid <- TRUE
      else {
        dx <- diff(x)
        newgrid <- max(abs(diff(dx))) < dx[1] * RFopt$general$gridtolerance
      }
      if ((missing(grid) || length(grid) == 0)) grid <- newgrid
      else if (xor(newgrid, grid) && RFopt$internal$warn_on_grid) {
        RFoptions(internal.warn_on_grid = FALSE)
        message("coordinates", if (grid) " do not",
                " seem to be on a grid, but grid = ", grid)
      }
    }
    len <- c(length(x), length(y), length(z))[1:spatialdim]
    
    if (!(missing(grid) || length(grid) == 0) && !grid) { ## sicher nicht grid, ansonsten ausprobieren
      if (any(diff(len) != 0)) stop("some of x, y, z differ in length")
      x <- cbind(x, y, z)
      ## make a matrix out of the list
      len <- len[1]
    } else {
      if ((missing(grid) || length(grid) == 0) && any(len != len[1]))
        grid <- TRUE
      x <- list(x, y, z)[1:spatialdim]
    }
    y <- z <- NULL ## wichtig dass y = NULL ist, da unten die Abfrage
  }  ## end of x, y, z given separately 
  
  if (!all(is.finite(unlist(x)))) stop("coordinates are not all finite")


  if ((missing(grid) || length(grid) == 0) || grid) {
    if (gridtriple) {
      if (len != 3)
        stop("In case of simulating a grid with option gridtriple, exactly 3 numbers are needed for each direction")
      lr <- x[3,] # apply(x, 2, function(r) length(seq(r[1], r[2], r[3])))
      ##x[2,] <- x[1,] + (lr - 0.999) * x[3,] ## since own algorithm recalculates
      ##                               the sequence, this makes sure that
      ##                               I will certainly get the result of seq
      ##                               altough numerical errors may occurs
      restotal <- prod(x[3, ])
      if (!is.null(y) && !all(y[3,] == x[3,]))
        stop("the grids of x and y do not match ")        
    } else {     
      xx <- seq2grid(x, "x",  grid,
                     RFopt$internal$warn_ambiguous, RFopt$general$gridtolerance)
     if (!is.null(y)) {
        yy <- seq2grid(y, "y", grid,
                       RFopt$internal$warn_ambiguous,
                       RFopt$general$gridtolerance)
        if (xor(is.logical(xx), is.logical(yy)) ||
            (!is.logical(xx) && !all(yy[3,] == xx[3,])))
          stop("the grids for x and y do not match")      
      }
      if (missing(grid) || length(grid) == 0) grid <- !is.logical(xx)       
      if (grid) {
        x <- xx
        if (!is.null(y)) y <- yy
        restotal <- prod(len)
        len <- 3
      } else {
        x <- sapply(x, function(z) z)
        if (!is.null(y)) y <- sapply(y, function(z) z)
      }
    }
    if (grid && any(x[3, ] <= 0))
      stop(paste("step must be postive. Got as steps",
                 paste(x[3,], collapse=",")))
    ##if (len == 1) stop("Use grid=FALSE if only a single point is simulated")
  }
 
  if (!grid) { ## not grid
    restotal <- nrow(x)
    if (is.null(y)) {
      if (restotal < 200 && any(as.double(dist(x)) == 0)) {## 2000
        d <- as.matrix(dist(x))
        diag(d) <- 1
        idx <-  which(as.matrix(d) ==0)
        if (printlevel>PL_ERRORS)
          Print(x, dim(d), idx , cbind( 1 + ((idx-1)%% nrow(d)), #
                                       1 + as.integer((idx - 1)  / nrow(d))) )
        warning("locations are not distinguishable")
      }
      ## fuer hoehere Werte con total ist ueberpruefung nicht mehr praktikabel
    }
  }

  if (coord_system == "earth") {
    # if (ncol(x) > 4) stop("earth coordinates have maximal 3 components")
    global.units <- RFoptions()$coords$new_coordunits[1]
    if (global.units[1] == "") global.units <- "km"

    Raumdim <- if (grid) ncol(x) else nrow(x)
    new_is_cartesian <- new_coord_system %in% CARTESIAN_SYSTEMS
    if (new_is_cartesian) {
      if (!is.null(y))
        stop("cannot take 'y' as an argument if 'x' is a (internal) matrix")
      code <- switch(new_coord_system,
                     "cartesian" = CARTESIAN_COORD,
                     "gnomonic" = GNOMONIC_PROJ,
                     "orthographic" = ORTHOGRAPHIC_PROJ,
                     stop("unknown projection method")
                     )
      x <- RFfctn(RMtrafo(new=code), x, grid=grid, 
                  coords.new_coordunits=global.units,
                  coords.new_coord_system = "keep")
        
      if (new_coord_system == "cartesian") {
        Raumdim <- max(3, Raumdim)
        spatialdim <- Raumdim
      }
      dim(x) <- c(length(x) /Raumdim, Raumdim)
      #x <- t(x)

      ## never try to set the following lines outside the 'if (new_coord_system'
      ## as in case of ..="keep" none of the following lines should be set
      RFoptions(coords.coord_system = 
                if (new_is_cartesian) "cartesian" else new_coord_system)
      grid <- FALSE
    } else if (!(new_coord_system %in% c("keep", "sphere", "earth"))) {
      warning("unknown new coordinate system")
    }
  }

  if (Zeit <- !is.null(T)) {
    Ttriple <- length(T) == 3;
    if (length(T) <= 2) Tgrid <- TRUE
      else {
        dT <- diff(T)
        Tgrid <- max(abs(diff(dT))) < dT[1] * RFopt$general$gridtolerance
      }
    if (is.na(RFopt$general$Ttriple)) {
      if (Ttriple && Tgrid)
        stop("ambiguous definition of 'T'. Set RFoptions(Ttriple=TRUE) or ",
             "RFoptions(Ttriple=FALSE)")
      if (!Ttriple && !Tgrid) stop("'T' does not have a valid format")
    } else if (RFopt$general$Ttriple) {
      if (!Ttriple)
        stop("'T' is not given in triple format 'c(start, step, length)'")
      Tgrid <- FALSE
    } else {
      if (!Tgrid) stop("'T' does not define a grid")
      Ttriple <- FALSE
    }
    if (Tgrid)
      T <- as.vector(seq2grid(list(T), "T", Tgrid,
                              RFopt$internal$warn_ambiguous,
                              RFopt$general$gridtolerance))
    restotal <- restotal * T[3]
  }

  if (!missing(dim) && !is.null(dim) && spatialdim != dim) {
    stop("'dim' should be given only when 'distances' are given. Here, 'dim' contradicts the given coordinates.")
  }
  
  storage.mode(x) <- "double"

  if (is.null(y)) {
   return(list(x=x, T=T, Zeit=Zeit, restotal=restotal, l=len,
                spatialdim=spatialdim, grid=grid, distances=FALSE,
                coordunits = curunits,  new_coordunits = newunits))
  } else {
    storage.mode(y) <- "double"
    return(list(x=x, y=y, T=T, Zeit=Zeit, restotal=restotal, l=len,
                spatialdim=spatialdim, grid=grid, distances=FALSE,
                coordunits = curunits, new_coordunits = newunits))
  }
  
}



RFearth2cartesian <- function(coord, units=NULL, system = "cartesian",
                              grid=FALSE) {
  if (is.character(system)) system <- pmatch(system, ISONAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (is.null(units)) {
    global.units <- RFoptions()$coords$new_coordunits[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coord)) coord <- t(coord)
  res <- RFfctn(RMtrafo(new=system), coord, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  return(res)
}

RFearth2dist <- function(coord, units=NULL, system="cartesian",
                         grid=FALSE, ...) {
  if (is.character(system)) system <- pmatch(system, ISONAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (is.null(units)) {
    global.units <- RFoptions()$coords$new_coordunits[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coord)) coord <- t(coord)
  z <- RFfctn(RMtrafo(new=system), coord, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  return(dist(z, ...))
}



CheckData <- function(model, x, y, z, T, distances, grid, data,
                      given=NULL, allowFirstCols=TRUE, vdim, repet,
                      dim, ...) {
  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")
  if (missing(data)) stop("missing data")

  if (missing(model)) model <- list("null")
  else model <- PrepareModel2(model, ...)
 
  missing.given <- length(given) == 0
  missing.x <- missing(x) && missing.given 

  if (isSpObj(data)) data <- sp2RF(data)
  if (isRFsp <- is(data, "RFsp")) {
    if (length(given) != 0) {
      stop("the coordinates for the measured data are given ambigiously (by 'given' and within 'data')")
    }
    data <- selectDataAccordingFormula(data, model=model)
    if (!is.null(data@.RFparams)) {
      vdim <- data@.RFparams$vdim
      repet <- data@.RFparams$n
    } else {
      vdim <- repet <- NULL
    }
    variab.names <- colnames(data@data)

    if (grid.fulldata <- isGridded(data)) {
      data <- rfspDataFrame2conventional(data) ## zeitkritisch ! to do
      given <- list(x=data$x, T=data$T, grid=TRUE)
      data <- as.matrix(data$data)
    } else {
      given <- list(x=coordinates(data), grid=FALSE) ## stehend; darf nicht gaendert werden
      data <- as.matrix(data@data)
    }
  } else {
    if (!hasArg("vdim")) vdim <- NULL
    if (!hasArg("repet")) repet <- NULL

   # Print(allowFirstCols)
    if (missing.x) {
      if (allowFirstCols)
        data.col <- data.columns(data, dim, force=TRUE, halt=FALSE)
      else data.col <- data.columns(data)
      given <- list(x=data[, data.col$x, drop=FALSE], grid=FALSE)
      data <- data[ , data.col$data, drop=FALSE]
    }
  }
  
  if (!missing.given) {
    if (!missing(x))
      stop("x coordinates of the data given in two different ways")
    stopifnot(missing(y) || is.null(y), missing(z) || is.null(z),
              missing(T) || is.null(T))
    if (!is.list(given)) given <- list(x=given)
  } 

  if (length(given) == 0)
    given <- CheckXT(x = x, y = y, z = z, T = T, grid = grid, 
                   distances=distances, 
                   if (!missing(distances) && length(distances) > 0)
                     dim = if (hasArg("dim")) dim else 1
                   ##, length.data=length(data)
                   )
  else
    given <- do.call(CheckXT, given)
  
  if (missing(dim) || length(dim) == 0)
    dim <- as.integer(given$spatialdim + given$Zeit)
  else stopifnot(dim == given$spatialdim + given$Zeit)
      
  if (is.vector(data)) {
    dimdata <- c(length(data), 1)
  } else {
    if (!is.array(data)) {
      data <- as.matrix(data)
    }
    dimdata <-  base::dim(data) 
  }


  
  if (length(dim(data)) != 2) {
    ## Achtung! data verliert die Row-/col-Namen mit der
    ## naechsten Anweisung!!
    datanames <- dimnames(data)
    base::dim(data) <- c(dimdata[1], prod(dimdata[-1]))
  }
 
  
  vdimrepet <- as.integer(length(data) / given$restotal)
  if (vdimrepet * given$restotal != length(data))
    stop("dimension of data does not match coordinates")
  
  variab.names <- colnames(data)
  names <- GetDataNames(model=model, coords=given$x, locinfo=given)#ohne data!
  if (is.null(names$variab.names))
    names$variab.names <-
      if (class(variab.names) == "try-error") NULL else variab.names

  return(c(ts.xdim = dim,
           vdim =vdim,
           repet = repet,
           list(## xgr = if (given$grid) cbind(given$x, given$T) else NULL,
                fullgiven = given,
                coord.names=names$coord.names,
                variab.names=names$variab.names,            
                model = model, #userdefined = userdefined,
                data.col = if (missing.x && !isRFsp) data.col else NULL,
                fulldata = data),
           missing.x=missing.x
           ))
}
