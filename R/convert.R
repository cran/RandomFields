
## since fitvario uses the returned ordering of the parameters by PrepareModel
## and convert.to.readable for the standard model definition, the ordering
## nugget/ main model should not be changed; 

PrepareModel <-  function(model, param, timespacedim, trend, method=NULL,
                          named=FALSE)
{
  
  ## any of the users model definition (standard, nested, list) for the
  ## covariance function is transformed into a standard format, used
  ## especially in the c programs
  ##
  ## overwrites in some situation the simulation method for nugget.
  ## allows trend to be NA (or any other non finite value  -- is not checked!)
  ## trend has not been implemented yet!
  
  PrintLevel <- RFparameters()$Print
  nugget.nr <- .C("GetModelNr", "nugget", as.integer(1), nr=integer(1),
                  PACKAGE="RandomFields")$nr
  STOP <- function(txt) {
    if (PrintLevel>1) {
      cat("model: ")
      if (!missing.model) str(model) else cat(" missing.\n")
      cat("param: ")
      if (!missing.param) str(param) else cat(" missing.\n")
      cat("timespacedim: ")
      if (!missing.timespacedim) str(timespacedim) else cat(" missing.\n")
      cat("trend: ")
      if (!missing.trend) str(trend) else cat(" missing.\n")
      cat("method: ")
      str(method)
    }
    stop(paste("(in PrepareModel)", txt), call.=FALSE)
  }

  
  op.list <- c("+", "*", "(")
  Mean <- 0
  if (missing.trend <- missing(trend)) trend <- NULL
  if (missing(method)) method <- NULL
  missing.model <- missing(model)
  if (missing.timespacedim <- missing(timespacedim)) timespacedim <- 1
  else timespacedim <- as.integer(timespacedim)
  missing.param <- missing(param) || is.null(param)
  if (length(timespacedim)==0) STOP("length(timespacedim)==0") 

  if (full.model <- missing.param && is.null(model$p)) { ## full model
    if (missing.model || (length(model)==0)) model <- list()
    else if (!is.list(model))
      STOP("if param is missing, model must be a list of lists")
    if (is.null(trend) + is.null(model$trend) + is.null(model$mean)<2)
      STOP("trend/mean is given twice")
    if (!is.null(model$mean)) Mean <- model$mean
    if (!is.null(model$trend)) trend <- model$trend
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
      return(list(covnr=integer(0), param=double(0), anisotropy=FALSE,
                  op=integer(0), mean=Mean, trend=trend, method=-1))
    if (length(model) %% 2 != 1) STOP("list for model definition should be odd")
    if (length(model)==1) op <- integer(0)
    else {
      op <- pmatch(c(model[seq(2, length(model), 2)], recursive=TRUE),
                   op.list, dup=TRUE) - 1
      if (!all(is.finite(op))) STOP("operators are not all known")
    }
    ## note: model is now reduced!
    model <- model[seq(1, length(model), 2)]
    if (sum(sapply(model, is.list))!=length(model)) 
      STOP("some elements of the model definition are not lists")
    
    ## extract the methods first since model and method start
    ## both with lower case m;
    ## model may be abbreviated by m, method by me
    methodX <- lapply(model, function(x) x$meth)
    model <- lapply(model, function(x) {
      if (!(is.na(p <- pmatch("meth", names(x), dup=TRUE)))) {
        if (length(p)!=1) STOP("method may be given only once")
        x[[p]] <- NULL
      }
      x
    })

    if (length(unlist(lapply(model, function(x) x$me))) > 0)
      STOP("'mean' seems to be given within the inner model definitions"); 
    if (!is.character(covnames <- sapply(model, function(x) x$m)))
      STOP("'model' was not given extacly once each odd number of list entries or additional unused list elements are given.")
    lm <-  as.integer(length(covnames))
    covnr <- .C("GetModelNr", covnames, lm, nr=integer(lm),
                PACKAGE="RandomFields")$nr
    names(covnr) <- covnames
        
    ## finalize the given methods
    explicite.method <- sapply(methodX, function(x) !is.null(x))
    if (length(method)>1 && length(method)!=lm ||
        length(method)<=1 && length(explicite.method)>0 &&
        length(explicite.method)!=lm)
      STOP("number of methods and number of models differ")
    if (any(explicite.method)) {
      methodX <- as.character(methodX)
      if (length(method)==0) method <- methodX
      else {
        method <- rep(method, len=lm)
        method[explicite.method] <- methodX
      }
    }

    k <- c(lapply(model, function(x) length(x$k)), recursive=TRUE)  ## used
    ##                                                 for checks later on 
    stopifnot(all(c(lapply(model, function(x) length(x$v)),
                    recursive=TRUE)==1))
    anull <- sapply(model, function(x) is.null(x$a))
    snull <- sapply(model, function(x) is.null(x$s))
    if (!all(anull) && !all(snull))
      STOP("do not mix isotropic definitions with anisotropic ones")
    scale.except <- integer(lm)
    aniso.except <- integer(lm)
    .C("GetModelException", covnr, lm, scale.except, aniso.except,
        PACKAGE="RandomFields", DUP=FALSE)
    if (anisotropy <- !all(anull)) {
      ## here, anisotropy parameter must be given also for nugget
      ## (think of higher dimensions, where nugget is only in one direction)
      if (any(anull & !aniso.except))
        STOP("anisotropy parameters are not given for all submodels")
      if (missing.timespacedim)
        timespacedim <- as.integer(sqrt(length(model[[which(!anull)[1]]]$a)))
      
      repl <- function(x) c(x, list(a=diag(rep(1.0, timespacedim))))
      model[anull] <- lapply(model[anull], repl)

      n.par <- 1 + timespacedim^2   
      param <- lapply(model, function(x) c(var=x$v, kappa=x$k,
                                           aniso=as.matrix(x$a), recursive=TRUE))
    } else {
      ## in case of isotropy, scale for nugget does not make sense
      ## so, if scale is not given (by default), add scale=1.0
      if (any(snull & !scale.except))
        STOP("scale parameter is not given for all submodels")
      model[snull] <- lapply(model[snull],
                             function(x) c(x, list(s=1.0)))
      n.par <- 2
      param <- lapply(model, function(x) c(var=x$v, kappa=x$k, scale=x$s))
    }
    ## as.character seems unnecessary, but internal representation does not seem
    ## to be character
  } else { ## standard definition or nested model
    if (missing.param) { ## a simple list of the model and the
      ##                    parameters is also possible
      if (is.null(param <- model$p)) STOP("is.null(model$param)")
      stopifnot(is.null(method) || is.null(model$meth))
      if (is.null(method)) method <- model$meth
      stopifnot(is.null(trend) || is.null(model$trend))
      if (is.null(trend)) trend <- model$trend
      if (!is.null(model$mean)) Mean <- model$mean
      model <- model$model
    }
    if (length(method)>1) STOP("only one method should be given")
    n.par <- 2
    anisotropy <- FALSE
    stopifnot(is.character(model), length(model)==1)
    if (is.matrix(param)) { ## nested models
      model <- rep(model, ncol(param))
      k <- rep(nrow(param) - 2, ncol(param))
      if (ncol(param)==1) op <- integer(0)
      else op <- rep(pmatch("+", op.list) - 1, ncol(param) - 1)
      param <- apply(param, 2, function(x) list(var=x[1], kappa=x[-1:-2],
                                                  scale=x[2]))
      for (i in 1:length(param)) if (param[[i]]$scale==0) {
        model[i] <- "nugget"
        param[[i]] <- c(nugget=param[[i]]$var, nugget.scale=1.0)
        k[i] <- 0
      }
      param <- lapply(param, function(x) unlist(x))
   } else if (is.vector(param)) {  ## standard, simple way
      ### falls trend gegeben, dann ist param um 1 Komponente gekuerzt
      if (is.null(trend)) {
        Mean <- param[1]
        param <- param[-1]
      } else if (Mean!=0) STOP("trend/mean given twice")
      k <- length(param) - 3
      
      if (is.na(param[2]) || (param[2]>0)) {## nugget
        param <- list(c(var=param[1], kappa=param[-1:-3], scale=param[3]),
                      c(nugget=param[2], nugget.scale=1.0))
        model <- c(model, "nugget") ## do not change ordering,
        ##                             since this ordering is used by other
        ##                             functions (fitvario)
        op <- pmatch("+", op.list) - 1
        ncov <- integer(2)
        k <- c(k, 0)
      } else {
        param <- list(c(var=param[1], kappa=param[-1:-3], scale=param[3]))
        op <- integer(0)
        ncov <- integer(1)
      }
    } else STOP("param does not have allowed format")
    if (Mean!=0 && !is.null(trend)) STOP("trend/mean are given twice")
    stopifnot(is.character(model)) 
    lm <-  as.integer(length(model))
    covnr <- .C("GetModelNr", model, lm, nr=integer(lm),
                 PACKAGE="RandomFields")$nr
    names(covnr) <- model
  } # end nested/standard definition

  ## common part:
  if (is.numeric(trend) && (length(trend)==1)) {
    Mean <- trend
    trend <- NULL
  } else {
    if (!is.na(Mean) && !is.numeric(Mean) || length(Mean)>1)
      stop("mean not a scalar")
  }
 
  if (any(covnr < 0)) {
    PrintModelList()
    STOP(paste("model(s) [",
               paste(paste("#", which(covnr<0), sep=""), collapse=", "),
               "] cannot be (uniquely) identified from the above list", sep=""))
  }

  kappas <- integer(lm)
  .C("GetNrParameters", as.integer(covnr), as.integer(lm),
     as.integer(timespacedim), kappas, PACKAGE="RandomFields", DUP=FALSE)

  if ( any(idx <- (sapply(param, length) != n.par + kappas) | (k!=kappas)) ) {
    STOP(paste("number of parameters in model #", which(idx), 
               " is not correct: ",
               kappas[idx], " genuine parameters and ",
               n.par, " parameters for variance and scale expected; got ",
               k[idx], " and ", sapply(param,length)[idx] - k[idx] ,
               "\n", sep=""))
  }
  
  if (is.null(method)) method <- -1
  else {
    if (any(method=="NULL")) STOP("either all methods must be specified or none")
    
    method <- .C("GetMethodNr", as.character(method),
                 as.integer(length(method)),
                 nr=integer(length(method)), PACKAGE="RandomFields")$nr
    if (any(method < 0)) {
      PrintMethodList()
      STOP(paste("given method(s) cannot be (uniquely) identified from the",
                 "above list: ", which(method<0))
           )
    }
    method <- rep(method, len=length(covnr))

    ## 30.8.05, now set in internal_InitSimulateRF:
    #if (!full.model)
    #  method[covnr==nugget.nr] <-
    #      .C("GetMethodNr", "nugget", as.integer(1), nr=integer(1),
    #        PACKAGE="RandomFields")$nr
  }
  
  param <- c(param,recursive=TRUE);
  if (!named) names(Mean) <- names(trend) <- names(covnr) <- names(param) <- NULL
  
  return(list(covnr=covnr,
              param=param,
              anisotropy=anisotropy,
              op=op,
              mean=Mean,
              trend=trend,
              method=method,
              timespacedim = timespacedim
         ))
}

convert.to.readable<- function(l, allowed=c("standard", "nested", "list")) {
  ## the inverse function for PrepareModel
  ## i.e. transform the coded covariance model into a user readable format
  ## takes always the simplest model out of the 'allowed' list
  assign(".methods", GetMethodNames())
  if (is.null(l$mean)) str(l)
  stopifnot(is.list(l),
            !is.null(l$covnr),
            !is.null(l$param),
            !is.null(l$anisotropy),
            !is.null(l$op),
            !is.null(l$mean), 
            !is.null(l$timespacedim)
            )
  PrintLevel <-  RFparameters()$Print
  op.list <- c("+", "*", "(")
  Allowed <- rep(FALSE, 3)
  Allowed[pmatch(allowed, c("standard","nested","list"))] <- TRUE 
  nugget.nr <- .C("GetModelNr", "nugget", as.integer(1), nr=integer(1),
                  PACKAGE="RandomFields")$nr
  ModelNames <- GetModelNames()
  
  method <- l$method
  names(method) <- names(l$param) <-  names(l$mean) <- names(l$trend) <- NULL
  if (method[1]==-1) method <- NULL
  stopifnot(all(method != -1))
  if (!is.null(method)) method <- .methods[1 + method]
  lc <- length(l$covnr)
  if (!is.null(l$trend)) {
    stopifnot(l$mean==0)
    l$mean <- NULL
  }
  dim.scale <- if (l$anisotropy) l$timespacedim^2 else 1
  k  <- c(0, cumsum(.C("GetNrParameters", as.integer(l$covnr), as.integer(lc),
                       as.integer(dim.scale),
                       k=integer(lc), PACKAGE="RandomFields", DUP=FALSE)$k +
                    1 + dim.scale))
  ol <- op.list[l$op + 1] == "+"
  ngg <- (l$covnr == nugget.nr) & c(ol,TRUE) & c(TRUE,ol)
  nuggets <- if (any(ngg)) l$param[k[c(ngg, FALSE)] + 1] else 0
  
  iso <- !l$anisotropy && (all(!is.na(nuggets)) || all(is.na(nuggets)))
  if (iso) {
    lp <- length(l$param)
    nugget <- sum(nuggets)
    if (all(ngg) && Allowed[1])
      return(list(model="nugget", param=as.double(c(l$mean, nugget, 0, 1)),
                  method=method, trend=l$trend)) # old.method 
 #   if ((lc==1) && Allowed[1])
 #      return(list(model=ModelNames[1 + l$covnr],
 #                 param=as.double(c(l$mean, l$param[1], 0, l$param[lp],
 #                   l$param[-c(1, lp)])),
 #                 method=method, trend=l$trend))
    pos <- which(!ngg)
    if ((length(pos)==1) && Allowed[1]) {  # old.method 
      if (!all(method[ngg] == "nugget") && PrintLevel>1)
        cat("nugget simulation methods overwritten.\n")  
      method <- method[!ngg]
      return(list(model=ModelNames[1 + l$covnr[pos]],
                  param=as.double(c(l$mean, l$param[k[pos]+1], nugget,
                    l$param[k[pos+1]],
                    if (k[pos+1]>=k[pos]+3)
                    l$param[(k[pos]+2):(k[pos+1]-1)])),
                  method=method, trend=l$trend)
             )
    }
  }
  
  if (iso && all(l$covnr[pos][1]==l$covnr[pos]) && Allowed[2]) {
    ## nested model
    param <- rbind(l$param[k[pos]+1], l$param[k[pos+1]]) # var and scale
    kappa <- NULL
    if (k[pos[1]+1]>=k[pos[1]]+3) ## i.e., there are additional parameters
      for (i in pos) 
        kappa <- cbind(kappa, l$param[(k[i]+2):(k[i+1]-1)])
    if (nugget==0) nugget <- NULL
    else {
      nugget <- c(nugget, rep(0, k[pos[1]+1]-k[pos[1]]-1))
      if (!all(method[ngg] == "nugget") && PrintLevel>1)
        cat("nugget simulation method overwritten.\n")  
      method <- method[!ngg]
    }
    method <- if (all(method==method[1])) method[1] else c(method, "nugget")
    ## do not change ordering of cbind(rbind( , ), nugget)
    ## since this ordering is used by other functions (mleRF)
    param <- cbind(rbind(param, kappa), nugget)
    attr(param, "dimnames") <- NULL 
    return(list(model=ModelNames[1 + l$covnr[!ngg][1]], param=param,
                method=method, mean=l$mean, trend=l$trend))
  }
  
  ## anisotropic
  model <- NULL
  if (!Allowed[3]) return(NULL)
  for (i in 1:lc) {
    if (l$anisotropy) {
      model <- c(model,list(list(model=ModelNames[1 + l$covnr[i]],
                                 var=l$param[k[i] + 1],
                                 kappas=if (k[i+1] >= k[i] + 2 + dim.scale)
                                    l$param[(k[i]+2):(k[i+1]-dim.scale)],
                                 aniso=l$param[(k[i+1]+1-dim.scale) : k[i+1]]
                                 )))
    } else {
      model <- c(model, list(list(model=ModelNames[1 + l$covnr[i]],
                                  var=l$param[k[i] + 1],
                                  kappas= if (k[i+1] >= k[i]+3)
                                     l$param[(k[i]+2): (k[i+1]-1)],
                                  scale=l$param[k[i+1]]
                                  )))
    }
    if (i<lc) model <- c(model, op=op.list[op=l$op[i]+1])
  }
  return(list(model=model, mean=l$mean, method=method, trend=l$trend))
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
  if (is.matrix(x)) {
    if (!is.null(y) || !is.null(z)) 
      stop("If x is a matrix, then y and z may not be given")
    spacedim <- ncol(x)
    l <- nrow(x)
    if (missing(grid) && spacedim==1) {## 1.10.03 from l==1
      dx <- diff(x)
      stopifnot(!gridtriple || l==3)
      grid <- max(abs(diff(dx))) < dx[1] * 1e-12
    } else stopifnot(is.logical(grid))
    if (grid && !gridtriple)
      ## list with columns as list elements -- easier way to do it??
      x <- lapply(apply(x, 2, list), function(r) r[[1]])
  } else { ## x, y, z given separately 
    if (is.null(y) && !is.null(z)) stop("y is not given, but z")
    spacedim <- 1 + (!is.null(y)) + (!is.null(z))
    if (missing(grid) && spacedim==1) {      
      dx <- diff(x)
      stopifnot(!gridtriple || length(x)==3)
      grid <- gridtriple || (max(abs(diff(dx))) < dx[1] * 1e-12)
    } else stopifnot(is.logical(grid))
    l <- c(length(x), length(y), length(z))[1:spacedim] 
    if (!grid || gridtriple) {
      if (any(diff(l) != 0))
        stop(if (gridtriple)
             "some of x, y, z are neither NULL nor have length 3" else
             "some of x, y, z differ in length")
      x <- cbind(x, y, z)
      ## make a matrix out of the list
      l <- l[1]
    } else {
      x <- list(x,y,z)[1:spacedim]
    }
    y <- z <- NULL
  }


  if (!all(is.finite(unlist(x)))) stop("coordinates are not all finite")
  if (grid) {
    if (gridtriple) {
      if (l != 3)
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
        if (any(step)==0)
          stop("duplicated values detected: the definition of coordinates does not seem to define a grid of equidistant coordinates")
        if (max(abs(step / step[1] - 1.0)) > 1e-13) {
          print(x[1:min(10000, length(x))])
          print(round(diff(diff(x)), 14))
          stop("different grid distances detected, but the grid must have equal distances in each direction -- try gridtriple=TRUE that avoids numerical errors.")
        }
        return(c(x[1], x[length(x)]+0.001*step[1], step[1]))
      }
      total <- prod(l)
      x <- sapply(x, eqdist)
      l <- 3
    }
    if (any(x[3, ] <= 0)) {
      print(x)    
      stop("step must be postive")
    }
    ##if (l == 1) stop("Use grid=FALSE if only a single point is simulated")
  } else {
    total <- nrow(x)
   if (total < 2000 && any(as.double(dist(x)) == 0))
     stop("locations must be distinguishable")
    ## fuer hoehere Werte con total ist ueberpruefung nicht mehr praktikabel
 }
  
  if (Time <- !is.null(T)) {
    stopifnot(length(T)==3)
    lT <- length(seq(T[1],T[2],T[3]))
    T[2] <- T[1] + (lT - 0.999) * T[3]
    total <- total * lT
  }
  return(list(x=x, T=T, Time=Time, total=total, l=l, spacedim=spacedim, grid=grid))
}

