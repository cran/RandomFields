"InitMaxStableRF" <-
function(x, y = NULL, z = NULL, grid, model, param, maxstable,
         method = NULL, register = 0, gridtriple = FALSE) 
{
  MaxStableList <- c("extremalGauss","BooleanFunction")
  MaxStableNr <- pmatch(maxstable,MaxStableList)
  if (is.na(MaxStableNr)) stop(paste("Unknown max-stable random field",
                                     "\nPossible values for `maxstable': \"",
                                     paste(MaxStableList,collapse="\", \""),
                                     "\"",sep=""))
  if (MaxStableNr==2) {
    if (is.null(method)) method <- "max.MPP"
    else{
      if (!is.character(method)) stop("Method must be a string.")
      if (.C("GetMethodNr",method, nr = integer(1))$nr !=
          .C("GetMethodNr","max.MPP", nr = integer(1))$nr) { 
        warning("Method does not match max-stable random field definition. Set to `max.MPP'.")  
        method <- "max.MPP"
      }
    }
  }
  return (InitSimulateRF(x=x, y=y, z=z, grid=grid, model=model, param=param,
                         method=method, register=register, gridtriple=gridtriple,
                         distribution="MaxStable")
          )
}

"MaxStableRF" <-
function (x, y = NULL, z = NULL, grid, model, param,  maxstable,
          method = NULL, n = 1, register = 0, gridtriple = FALSE) 
{
    if (InitSimulateRF(x=x, y=y, z=z, grid=grid, model=model, param=param,
                       method=method, register=register, gridtriple=gridtriple,
                       distribution="MaxStable") <= 0) {
        return(DoSimulateRF(n=n, reg=register))
    }
    else {
        return(NULL)
    }
}
