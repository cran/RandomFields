
setMethod("[", signature=c("RFsp"), def=brack)
setMethod("[", signature=c("RFspatialGridDataFrame"), def=brack)
setMethod("[", signature=c("RFspatialPointsDataFrame"), def=brack)

setMethod("[<-", signature=c("RFsp"),
          function(x, i, j, ..., value) {
            dots = list(...)
            if (length(dots)>0) warning("dots are ignored")
            if (missing(j)) 
              x@data[i] <- value
            else
              x@data[i,j] <- value
            return(x)
          })

trafo_pointsdata <- function(x, dim) {
  if (isgrid <- is(x, "RFgridDataFrame")) {
 #   Print(x)
    x <- as(x, "RFpointsDataFrame")  ## funktioniert nicht
    ##    Print(x)
  } else if ((is(x, "matrix") || is(x, "data.frame")) && !missing(dim)) {
    dc <- data.columns(x, xdim = dim, force=TRUE)
    x <- list(coords=x[, dc$x, drop=FALSE], data=x[, dc$data, drop=FALSE])
  } else {
    if (!is(x, "RFpointsDataFrame"))
      stop("method only for objects of class 'RFpointsDataFrame' and 'RFgridDataFrame'")
  }
             
  dummy <- dimnames(x@coords)[[2]][1]
  lab <- xylabs(if (is.null(dummy)) "" else dummy, "",
                units=x@.RFparams$coord.units)
  labdata <- names(x@data)
  colname <- colnames(x@data)
  if (isgrid) {
    return(list(coords=as.vector(x@coords),
                data=as.matrix(x@data),
                RFparams=x@.RFparams,
                lab=lab, labdata=labdata, colnames=colname))
  } else {
    ord <- order(x@coords)
    return(list(coords=x@coords[ord, ],
                data=as.matrix(x@data)[ord, , drop=FALSE],
                RFparams=x@.RFparams,
                lab=lab, labdata=labdata, colnames=colname))
  } 
}
