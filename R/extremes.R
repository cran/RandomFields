
## it does not make sense to me at the moment that a space-time model
## for extremes is defined.

"InitMaxStableRF" <-
function(x, y = NULL, z = NULL, grid, model, param, maxstable,
         method = NULL, register = 0, gridtriple = FALSE) 
{
  MaxStableList <- c("extremalGauss","BooleanFunction")
  stopifnot(length(maxstable)==1)
  MaxStableNr <- pmatch(maxstable,MaxStableList)
  if (is.na(MaxStableNr)) stop(paste("Unknown max-stable random field",
                                     "\nPossible values for `maxstable': \"",
                                     paste(MaxStableList,collapse="\", \""),
                                     "\"",sep=""))
    if (MaxStableNr==2) {
    if (is.null(method)) method <- "max.MPP"
    else{
      if (!is.character(method) || (length(method)==0))
        stop("Method must be a string.")
      if (.C("GetMethodNr", as.character(method),
             as.integer(1), nr = integer(1), PACKAGE="RandomFields")$nr
          !=
          .C("GetMethodNr",as.character("max.MPP"),
             as.integer(1), nr = integer(1), PACKAGE="RandomFields")$nr) { 
        warning("Method does not match max-stable random field definition. Set to `max.MPP'.")
      if (any(method[1] != method))
        warning("method='max.MPP' cannot be mixed with other methods")
        method <- "max.MPP"
      }
    }
  }
  return(InitSimulateRF(x=x, y=y, z=z, grid=grid, model=model, param=param,
                        method=method, register=register,
                        gridtriple=gridtriple,
                        distribution="MaxStable")
          )
}

"MaxStableRF" <-
function (x, y = NULL, z = NULL, grid, model, param,  maxstable,
          method = NULL, n = 1, register = 0, gridtriple = FALSE, ...) 
{
   old.param <- RFparameters(no.readonly=TRUE)
  if (n>1 && !old.param$Storing) {
    RFparameters(Storing=TRUE)
  }
  RFpar <- list(...)
  if (length(RFpar)>0) RFparameters(RFpar)
  on.exit({RFparameters(old.param);
           if (!old.param$Storing) DeleteRegister(register)})
  error <- InitMaxStableRF(x=x, y=y, z=z, grid=grid, model=model, param=param,
                           maxstable=maxstable, method=method,
                           register=register, gridtriple=gridtriple)
  if (error > 0) stop(paste("InitMaxStable: error", error, "occured"))

  str(GetRegisterInfo())
   
  return(DoSimulateRF(n=n, reg=register, paired=FALSE))
}








