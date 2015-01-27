

summary.RFcrossvalidate <- function(object, ...) {
  ## cf geoR
  if (is.null(object$error)) {
    #Print(names(object))
    res <- list()
    for (i in 1:length(object)) {
      res[[i]] <- summary.RFcrossvalidate(object[[i]])
    }
    names(res) <- names(object)
  } else {
    F <- function(object) {    
      paste(formatC("sd", flag="-", width=7), ":", sep="",
            formatC(object, digits=2, width=6, format = "f"))
    }
    res <- list()
    res$error <- rbind(summary(object$error), F(apply(object$error, 2, sd)))
    res$std.error <- rbind(summary(object$std.error), F(apply(object$std.error, 2, sd)))
  }
  class(res) <- "summary.RFcrossvalidate"  
  return(res)
}


print.RFcrossvalidate <- function(x, ...) {
#  str(x)
  print.summary.RFcrossvalidate(summary.RFcrossvalidate(x, ...))
}

print.summary.RFcrossvalidate <- function(x, ...) {
  ## cf geoR
  if (is.null(x$error)) {
    for (i in 1:length(x)) {
      cat("\n", names(x)[i], "\n")
      print.summary.RFcrossvalidate(x[[i]])
    }
    return(invisible(x))
  } else {
    colnames(x$error) <- paste("err:", colnames(x$error), sep="")
    colnames(x$std.error) <- paste("std:", colnames(x$std.error), sep="")
    res <- cbind(x$error, x$std.error)
    print(res, quote=FALSE) #
    return(invisible(res))
  }
}

RFcrossvalidate <- function(model, x, y=NULL, z=NULL, T=NULL, grid=NULL, data,
                            lower=NULL, upper=NULL,
                            bc_lambda, ## if missing then no BoxCox-Trafo
                            method="ml", # "reml", "rml1"),
                            users.guess=NULL,  
                            distances=NULL,
                            dim, ## space time dim -- later a vector of grouped dimensions
                            optim.control=NULL,
                            transform=NULL,
                            full = FALSE,
                            ...) {


  RFoptOld <- internal.rfoptions(seed=NA, ...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  refit <- RFopt$fit$cross_refit
  return.variance <- RFopt$krige$return.variance
  general <- RFopt$general
  printlevel <-general$printlevel
  details <- general$detailed_output
  if (general$modus_operandi == "neurotic") stop("crossvalidation is not a precise method")
  if (length(method) != 1) stop("exactly 1 method must be given")
   pch <- general$pch
  if (printlevel < 1) pch <- ""
 
  if (method %in%LSQMETHODS){
    sub.methods <- method
    methods <- NULL
  } else if (method %in% MLMETHODS) {
    sub.methods <- LSQMETHODS
    methods <- method
  } else stop("unknown method")
  
  res <- list()
  if (full) {
    ats <- approx_test_single(model, method=method)
    result <- ats$result
    models <- ats$fitted.model
  } else {
    if (class(model) == "RF_fit") {
      models <- list(model[[method]])
    } else if (class(model) == "RFfit") {
      models <- list(PrepareModel2(model[method]))
    } else {
      models <- list(model)
    }
  }

  
  for (i in 1:length(models)) {
    if (class(model) == "RFfit" || class(model) == "RF_fit") {
      fit <- models[[i]]
      if (details) fitted <- list(fit)
      refit <- FALSE
   #   if (missing(dim))
   #     dim <- if (class(model) == "RFfit") model@true.tsdim else model$true.dim
    }
    if (class(model) == "RFfit") {
      if (missing(data)) Z <- model@Z
      else {
        Z <- StandardizeData(x=x, y=y, z=z, T=T, grid=grid, data=data,
                             distances=distances, dimensions=dim,
                             RFopt=RFopt)
      }
    } else { # class(model) == "RF_fit" || just !refit
      if (!refit && class(model) != "RF_fit") {
        fit <- RFfit(model=model, x=x, y=y, z=z, T=T, grid=grid, data=data,
                     lower=lower, upper=upper, bc_lambda = bc_lambda,
                     methods=methods, sub.methods=sub.methods,
                     optim.control=optim.control, users.guess = users.guess,
                     distances = distances, dim=dim,
                     transform = transform, spConform=FALSE)
        if (details) fitted <- list(fit)
      }
      
      Z <- StandardizeData(x=x, y=y, z=z, T=T, grid=grid, data=data,
                           distances=distances, dimensions=dim,
                           RFopt=RFopt)
    }
   
    if (length(Z$data) > 1)
      stop("cross-valdation currently only works for single sets")

    predicted <- var <- error <- std.error <- Z$data[[1]]
    fitted <- list()
    len <- Z$dimdata[1, 1]
    Zcoord <- t(Z$coord[[1]])

    #oldRF <- RFoptions(); Print("A")
    
    for (nr in 1:len) {
#      Print(nr, len, RFoptions(), all.equal(oldRF, RFoptions()));     
       if (length(base::dim(Z$data[[1]])) == 2)
        Daten <- Z$data[[1]][-nr,   , drop=FALSE] ## nicht data -- inferiert!!
      else # dim == 3 assumed
        Daten <- Z$data[[1]][-nr, ,  , drop=FALSE]
      coords <- Zcoord[-nr, , drop=FALSE]
      if (printlevel >= 1 && pch!="") cat(pch)
      if (refit) {
        if (!is.null(Z$distances)) {
          stop("not programmed yet. Please contact author") # to do
        } else {
          fit <- RFfit(model=model,
                       x=coords, grid=grid, data=Daten,
                       lower=lower, upper=upper, bc_lambda = bc_lambda,
                       methods=methods, sub.methods=sub.methods,
                       optim.control=optim.control, users.guess = users.guess,
                       distances = distances, dim=dim,
                       transform = transform, spConform=FALSE)
          if (details) fitted[[nr]] <- fit
        }
      }
      
      interpol <- RFinterpolate(fit, x = Zcoord[nr, , drop=FALSE],
                                grid=FALSE,
                                data = cbind(coords, Daten),
                                 spConform = FALSE,
                                return_variance=TRUE, pch="")
      var[nr, ] <- as.vector(interpol$var)
      predicted[nr, ] <- as.vector(interpol$estim)
    }
    
    if (printlevel <= 1 && pch!="") cat("\n")

    error <- Z$data[[1]] - predicted
    std.error <- error / sqrt(var)
    
    res[[i]] <- list(data=Z$data, predicted=predicted, krige.var=var,
                error=error, std.error=std.error,
                p = 2 * pnorm(-abs(std.error)))
    if (details) res[[i]]$fitted <- fitted
    class(res[[i]]) <- "RFcrossvalidate"
  }
  

  if (length(res) == 1) res <- res[[1]]
  else {
     names(res) <- c(as.character(result$model1.name),
                     as.character(result$model2.name)[1])
     class(res) <- "RFcrossvalidate"
   }

  if (RFopt$general$returncall)
    attr(res, "call") <- as.character(deparse(match.call())) 
  attr(res, "coord_system") <- c(orig=RFopt$coords$coord_system,
                                 model=RFopt$coords$new_coord_system)
  res
}

