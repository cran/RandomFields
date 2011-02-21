

PrepareModel <-  function(model, param, trend=NULL, method=NULL,
                          nugget.remove=TRUE) {
  ## any of the users model definition (standard, nested, list) for the
  ## covariance function is transformed into a standard format, used
  ## especially in the c programs
  ##
  ## overwrites in some situation the simulation method for nugget.
  ## allows trend to be NA (or any other non finite value  -- is not checked!)
  ## trend has not been implemented yet!
 
  if (is.list(model) && is.character(model[[1]]) &&
      (is.null(names(model)) || names(model)[[1]]=="")) {
    if (!is.null(method)) {
      if (!is.null(model$method)) stop("method given twice")
      model <- c(model[1], method = method, model[-1])
    }
    if(!missing(param) && !is.null(param))
      stop("param cannot be given in the extended definition")
    return(list(model=model, trend=trend))
  }
    
#  PrintLevel <- RFparametersrf()$Print
  PrintLevel <- 2
  STOP <- function(txt) {
    if (PrintLevel>1) {
      cat("model: ")
      if (!missing.model) Print(model) else cat(" missing.\n") #
      cat("param: ")
      if (!missing.param) Print(param) else cat(" missing.\n") #
      cat("trend: ")
      Print(trend) # 
      Print(method) #
    }
    stop("(in PrepareModel) ", txt, call.=FALSE)
  }
   
  transform <- function(model) {
    if (!is.list(model)) {
      STOP("some elements of the model definition are not lists")
    }
    m <- list("$", var=model$v)
    lm <- length(model) - 3 # var, scale/aniso, name
    if (!is.null(model$a)) m$aniso <- model$a else m$scale <- model$scale
    model <- c(model, if (!is.null(model$a))
               list(aniso=model$a) else list(scale=model$s)) ## ???
    if (is.na(p <- pmatch("meth", names(model), duplicates.ok=TRUE))) {
      if (!is.null(method)) m$method <- method
    } else {
      if (length(p) > 1) STOP("method may be given only once")
      lm <- lm - 1
      m$method <- model$meth
      model[[p]] <- NULL
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
      Print(lm, model) 
      stop("some parameters do not fit")
    }
    m <- c(m, list(m1))

    ##    Print(m)
    return(m)
    
  } # end transform

  op.list <- c("+", "*")  ## if others use complex list definition !
  missing.model <- missing(model)
  missing.param <- missing(param) || is.null(param)

  if (full.model <- missing.param && is.null(model$param)) { ## full model
    warning("the sequential list format is depreciated.")
    if (missing.model || (length(model)==0)) model <- list()
    else if (!is.list(model))
      STOP("if param is missing, model must be a list of lists (or a list in the extended notation)")
    if (is.null(trend) + is.null(model$mean) + is.null(model$trend)<2)
      STOP("trend/mean is given twice")
    if (!is.null(model$mean)) trend <- model$mean else
    if (!is.null(model$trend)) trend <- model$trend else trend <- NULL
    if (!is.null(model$meth)) {
      if (!is.null(method)) STOP("method is given twice")
      method <- model$meth
    }

    model$trend <- model$mean <- model$meth <- NULL
    ## the definition might be given at a deeper level as element
    ## $model of the list:
    if (is.list(model$model)) {
      if (!is.list(model$model[[1]]))
        STOP("if param is missing, the model$model must be a list of lists")
      model <- model$model
    }
    if (length(model)==0) ## deterministic 
      return(list(model=list(), trend = trend))
    if (length(model) %% 2 !=1) STOP("list for model definition should be odd")
    if (length(model)==1)
      return(list(model = transform(model[[1]]), trend=trend))

    op <- pmatch(c(model[seq(2, length(model), 2)], recursive=TRUE),
                 op.list, dup=TRUE) - 1
    if (!all(is.finite(op))) STOP("operators are not all allowed; see the extended list definition for extensions")
    model <- model[seq(1, length(model), 2)]

    plus <- which(op==0)
    if (length(plus) == 0) {
      m <- list("*", lapply(model, transform))
    } else {
      plus <- c(0, plus, length(op)+1)
      m <- list("+")
      for (i in 1:(length(plus) - 1)) {
        m[[i+1]] <-
          if (plus[i] + 1 == plus[i+1]) transform(model[[plus[i] + 1]])
          else list("*", lapply(model[(plus[i] + 1) : plus[i+1]], transform))
      }
    }
   model <- m     
  } else { ## standard definition or nested model
    if (missing.param) { ## a simple list of the model and the
      ##                    parameters is also possible
      if (is.null(param <- model$p)) STOP("is.null(model$param)")
      stopifnot(is.null(method) || is.null(model$meth))
      if (is.null(method))
        if (!is.null(model$meth)) method <- model$meth
      stopifnot(is.null(trend) || is.null(model$trend))
      if (is.null(trend)) trend <- model$trend
      if (!is.null(model$mean)) {
        if (!is.null(trend)) STOP("mean and trend given twice")
        trend <- model$mean
      }
      model <- model$model
    }
    if (length(method)>1) STOP("only one method should be given")
    stopifnot(is.character(model), length(model)==1)
    if (is.matrix(param)) { ## nested
      if (nrow(param) == 1)
        return(PrepareModel(model=model, param=c(param[1], 0, param[-1]),
                            trend=trend, method=method))
      name <- model
      model <- list("+")#, method=method)
      for (i in 1:nrow(param)) {
        model <- c(model,
                   if (is.na(param[i, 2]) || param[i, 2] != 0)
                   list(list("$", var=param[i, 1], scale=param[i, 2],
                             if (ncol(param) >2) list(name, k=param[i,-1:-2])
                             else list(name)))
                   else list(list("$", var=param[i,1], list("nugget"))))
      }
    } else if (is.vector(param)) {  ## standard, simple way
      ## falls trend gegeben, dann ist param um 1 Komponente gekuerzt
      if (is.null(trend)) {
        trend <- param[1]
        param <- param[-1]
      } else warning("it is assumed that no mean is given so that the first component of param is the variance")
      if (model == "nugget") {
        model <- transform(list(model=model, var=sum(param[1:2]), scale=1)) 
      } else {
        if  (length(param) > 3)
          model <- transform(list(model=model, var=param[1], scale=param[3],
                                  k=param[-1:-3], method=method))
        else 
          model <- transform(list(model=model, var=param[1], scale=param[3],
                                  method=method))
        if (is.na(param[2]) || param[2] != 0 || !nugget.remove) {# nugget
          model <- list("+",
                        model,
                        transform(list(model="nugget", var=param[2], scale=1)))
        }
        ## if (!is.null(method)) model <- c(model, method=method) ## doppelt
      }
    } else stop("unknown format")  # end nested/standard definition
  }
  return(list(model = model, trend=trend))
}




CheckXT <- function(x, y, z, T, grid, gridtriple){
 
  ## converts the given coordinates into standard formats
  ## (one for arbitrarily given locations and one for grid points)

  if (is.data.frame(x)) {
    if (ncol(x)==1) x <- as.vector(x) else x <- as.matrix(x)
  }
  stopifnot(length(x) != 0, #, !grid || is.logical(gridtriple),
            all(is.finite(x)), all(is.finite(y)), all(is.finite(z))
  )
  ## if (grid && !gridtriple) convert (x,y,z) to list(x,y,z)
  ## else to matrix cbind(x,y,z)
  if (missing(grid) && gridtriple) grid <- TRUE
  if (is.matrix(x)) {
    if (!is.null(y) || !is.null(z)) 
      stop("If x is a matrix, then y and z may not be given")
    spacedim <- ncol(x)
    len <- nrow(x)
    if (missing(grid) && spacedim==1) {
      dx <- diff(x)
      grid <- max(abs(diff(dx))) < dx[1] * 1e-12
    } else {
      stopifnot(is.logical(grid))
      if (missing(gridtriple)) gridtriple <- len==3
    } 
    if (grid && !gridtriple)
      ## list with columns as list elements -- easier way to do it??
      x <- lapply(apply(x, 2, list), function(r) r[[1]])
  } else { ## x, y, z given separately 
    if (is.null(y) && !is.null(z)) stop("y is not given, but z")
    spacedim <- 1 + (!is.null(y)) + (!is.null(z))
    if (missing(grid) && spacedim==1) {      
      dx <- diff(x)
      grid <- gridtriple || (max(abs(diff(dx))) < dx[1] * 1e-12)
    } else {
      stopifnot(is.logical(grid))
      if (missing(gridtriple)) gridtriple <- length(x) == 3
    }
    len <- c(length(x), length(y), length(z))[1:spacedim] 
    if (!grid || gridtriple) {
      if (any(diff(len) != 0))
        stop(if (gridtriple)
             "some of x, y, z are neither NULL nor have length 3" else
             "some of x, y, z differ in length")
      x <- cbind(x, y, z)
      ## make a matrix out of the list
      len <- len[1]
    } else {
      x <- list(x,y,z)[1:spacedim]
    }
    y <- z <- NULL
  }

  if (!all(is.finite(unlist(x)))) stop("coordinates are not all finite")
  if (grid) {
    if (gridtriple) {
      if (len != 3)
        stop("In case of simulating a grid with option gridtriple, exactly 3 numbers are needed for each direction")
      lr <- apply(x, 2, function(r) length(seq(r[1],r[2],r[3])))
      x[2,] <- x[1,] + (lr - 0.999) * x[3,] ## since own algorithm recalculates
      ##                               the sequence, this makes sure that
      ##                               I will certainly get the result of seq
      ##                               altough numerical errors may occurs
      total <- prod(lr)
    } else {
      eqdist <- function(x) {
        step <- diff(x)
        if (any(step == 0.0))
          stop("duplicated values detected: the definition of coordinates does not seem to define a grid of equidistant coordinates")
        if (max(abs(step / step[1] - 1.0)) > 1e-13) {
          diffx <- round(diff(diff(x)))
          Print(x[1:min(10000, length(x))], diffx[1:min(10000,length(diffx))], 14)
          stop("different grid distances detected, but the grid must have equal distances in each direction -- try gridtriple=TRUE that avoids numerical errors.")
        }
        return(c(x[1], x[length(x)]+0.001*step[1], step[1]))
      }
      total <- prod(len)
      x <- sapply(x, eqdist)
      len <- 3
    }
    if (any(x[3, ] <= 0)) {
      Print(x)    
      stop("step must be postive")
    }
    ##if (len == 1) stop("Use grid=FALSE if only a single point is simulated")
  } else {
    total <- nrow(x)
    if (total < 200 && any(as.double(dist(x)) == 0)) {## 2000
      d <- as.matrix(dist(x))
      diag(d) <- 1
      idx <-  which(as.matrix(d) ==0)
      Print(x, dim(d), idx , cbind( 1 + ((idx-1)%% nrow(d)),
                                   1 + as.integer((idx - 1)  / nrow(d))) )
      stop("locations must be distinguishable")
   }
    ## fuer hoehere Werte con total ist ueberpruefung nicht mehr praktikabel
 }
  
  if (Time <- !is.null(T)) {
    stopifnot(length(T)==3)
    lT <- length(seq(T[1],T[2],T[3]))
    T[2] <- T[1] + (lT - 0.999) * T[3]
    total <- total * lT
  }
  return(list(x=x, T=T, Time=Time, total=total, l=len, spacedim=spacedim,
              grid=grid))
}

