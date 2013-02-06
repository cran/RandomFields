
## it does not make sense to me at the moment that a space-time model
## for extremes is defined.

InitMaxStableRF <- function(x, y = NULL, z = NULL, grid, model, param,
                            maxstable,
                            method = NULL, register = 0, gridtriple = FALSE) {
  MaxStableList <- c("extremalGauss","BooleanFunction","BrownResnick")
  stopifnot(length(maxstable)==1)
  MaxStableNr <- pmatch(maxstable,MaxStableList)
  if (is.na(MaxStableNr)) stop("Unknown max-stable random field",
                               "\nPossible values for `maxstable': \"",
                               paste(MaxStableList,collapse="\", \""),
                               "\"",sep="")
  if (MaxStableNr==1) {
    if (is.null(method)) method <- "extremalGauss"
  } else if (MaxStableNr==2) { # maxstablenr ==2
    if (is.null(method)) method <- "max.MPP"
    else {
      if (!is.character(method) || (length(method)==0))
        stop("Method must be a string.")
      if (.C("GetMethodNr", as.character(method), nr = integer(1),
             PACKAGE="RandomFields")$nr
          !=
          .C("GetMethodNr",as.character("max.MPP"), nr = integer(1),
             PACKAGE="RandomFields")$nr) { 
        warning("Method does not match max-stable random field definition. Set to 'max.MPP'.")
      if (any(method[1] != method))
        warning("method='max.MPP' cannot be mixed with other methods")
        method <- "max.MPP"
      }
    }
  }
  
  return(InitSimulateRF(x=x, y=y, z=z, grid=grid, model=model, param=param,
                        trend=NULL, method=method, register=register,
                        gridtriple=gridtriple,
                        distribution="MaxStable")
          )
}

MaxStableRF <- function (x, y = NULL, z = NULL, grid, model, param, maxstable,
                         method = NULL, n = 1, register = 0,
                         gridtriple = FALSE,
                         ...) {
  old.param <- RFparameters(no.readonly=TRUE)
  RFpar <- list(...)
  if (length(RFpar)>0) RFparameters(RFpar)
  if (delete <- n>1 && !RFparameters()$Storing) RFparameters(Storing=TRUE)
  on.exit({if (delete) DeleteRegister(register);
           RFparameters(old.param);
         })
  error <- InitMaxStableRF(x=x, y=y, z=z, grid=grid, model=model, param=param,
                           maxstable=maxstable, method=method,
                           register=register, gridtriple=gridtriple)
  if (error > 0) stop("InitMaxStable: error ", error, " occured")

  
  return(if (n<1) NULL else DoSimulateRF(n=n, register=register, paired=FALSE))
}

