

AddUnits <- function(params) {
  ## see also empvario.R and fitgauss.R, if changed
  general <- RFoptions()$general
  return(c(params, list(coord.units=general$coord_units,
                        variab.units=general$variab_units)))
}

compareGridBooleans <- function(grid, gridtmp) {
  if (!missing(grid) && length(grid)>0 && grid!=gridtmp)
    message(paste("you specified grid=", as.character(grid),
                  " but isGridded(data)=", as.character(gridtmp),
                  ";  grid is set to ", as.character(gridtmp), sep=""))
}


## Generate Objects ########################################################

RFspatialGridDataFrame <- function(grid, data,
                                   proj4string = CRS(as.character(NA)),
                                   RFparams=list()) {
  grid <- convert2GridTopology(grid)
  tmp <- SpatialGridDataFrame(grid=grid, data=data, proj4string=proj4string)
  tmp <- as(tmp, "RFspatialGridDataFrame")
  tmp@.RFparams <- AddUnits(RFparams)
  validObject(tmp)
  return(tmp)
}

RFspatialPointsDataFrame <- function(coords, data, coords.nrs = numeric(0),
                                     proj4string = CRS(as.character(NA)), 
                                     match.ID = TRUE, bbox = NULL,
                                     coord.units = NULL,
                                     variab.units = NULL,
                                     RFparams=list()) {
  if (is.null(bbox)) {
    bbox <- t(apply(coords, 2, range))
    colnames(bbox) <- c("min", "max")
  }
  tmp <- SpatialPointsDataFrame(coords=coords, data=data, coords.nrs=coords.nrs,
                                proj4string=proj4string, 
                                match.ID=match.ID, bbox=bbox)
  tmp <- as(tmp, "RFspatialPointsDataFrame")
  tmp@.RFparams <- AddUnits(RFparams)

#  if (!is.null(dimnames(coords)))
#    dimnames(tmp@coords) <- dimnames(coords)
  validObject(tmp)
  return(tmp)
}

RFgridDataFrame <- function(data, grid,
                            RFparams=list()){
  grid <- convert2GridTopology(grid)
  data <- as.data.frame(data)
  return(new("RFgridDataFrame", data=data, grid=grid,
             .RFparams=AddUnits(RFparams)))
}

RFpointsDataFrame <- function(data=data.frame(NULL), coords=as.numeric(NULL),
                              RFparams=list()){
  data <- as.data.frame(data)
  if (is.null(dim(coords))) coords <- matrix(coords)
  return(new("RFpointsDataFrame", data=data, coords=coords,
             .RFparams=AddUnits(RFparams)))
}


convert2GridTopology <- function(grid){
  if (!is(grid, "GridTopology")) {
    if (is.null(dim(grid)))
      grid <- matrix(grid, ncol=1)
    stopifnot(nrow(grid)==3)
    grid <- GridTopology(cellcentre.offset=grid[1,],
                         cellsize=grid[2,],
                         cells.dim=grid[3,])
  }
  return(grid)
}


## Coerce Objects #########################################################

setAs("RFspatialGridDataFrame", "RFspatialPointsDataFrame",
      function(from, to) RFspatialPointsDataFrame(coords=coordinates(from),
                                                  data=from@data,
                                                  RFparams=from@.RFparams))
setAs("RFspatialPointsDataFrame", "RFspatialGridDataFrame",
      function(from, to) RFspatialGridDataFrame(grid=points2grid(from),
                                                data=from@data,
                                                RFparams=from@.RFparams))

setAs("RFpointsDataFrame", "RFgridDataFrame",
      function(from, to) {
        tmp.coord <- SpatialPoints(cbind(from@coords, 1:length(from@coords)))
        tmp.grid <- points2grid(tmp.coord)
        grid <- convert2GridTopology(tmp.grid[,1, drop=FALSE])
        RFgridDataFrame(data=from@data, grid=grid, RFparams=from@.RFparams)
      })
setAs("RFgridDataFrame", "RFpointsDataFrame", 
      function(from, to) {
        coords <- GridTopology2gridVectors(from@grid)[[1]]
        RFpointsDataFrame(data=from@data, coords=coords,
                          RFparams=from@.RFparams)
      })



brack <- function(x, i, j, ..., drop=FALSE) {
  dots = list(...)
  if (length(dots)>0) warning("dots are ignored")
  has.variance <- !is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance
  if (missing(j)) {
    x@data <- x@data[i]#, drop=drop]
    n <- x@.RFparams$n
    v <- x@.RFparams$vdim
    if (!is.numeric(i)) {
      if (is.logical(i)) {
        i <- which(i)
      } else {
        stopifnot(all(i %in% colnames(x@data)))
        i <- match(i, colnames(x@data))
      }
    }
    if (! (length(unique(table(i%%v, rep(0, length(i)))))==1) )
      stop(paste("for each variable selected, the same number of repetitions ",
                 "must be selected; you selected columns ",
                 paste(i, collapse=","), " but vdim=",v," and n=",n, sep=""))
    x@.RFparams$vdim <- v.new <- length(unique(i%%v))
    if (ret.has.var <- has.variance && any(i > n*v))
      x@.RFparams$has.variance <- ret.has.var
    x@.RFparams$n <- length(i) / v.new - ret.has.var
    
  }
  else
    x@data <- x@data[i,j]
  return(x)
}


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


## methods 'as.matrix' and 'cbind' as in the 'sp' package

as.matrix.RFspatialGridDataFrame = function(x, ..., byrow = FALSE) {
  if (ncol(x@data) > 1)
    warning("as.matrix.RFspatialGridDataFrame uses first column;
            \n pass subset or [] for other columns;
            \n only the first two space dimensions are used")
  matrix(x@data[[1]], x@grid@cells.dim[1], x@grid@cells.dim[2], byrow=byrow)
}

as.matrix.RFgridDataFrame = function(x, ..., byrow = FALSE) {
  if (ncol(x@data) > 1)
    warning("as.matrix.RFgridDataFrame uses first column;
            \n pass subset or [] for other columns;
            \n only the first two space dimensions are used")
  matrix(x@data[[1]], x@grid@cells.dim[1], 1, byrow=byrow)
}

cbind_RFsp = function(...) {  ##copied from sp package
  stop.ifnot.equal = function(a, b) {
    res = all.equal(a@grid, b@grid)
    if (!is.logical(res) || !res)
      stop("grid topology is not equal")
  }
  grds = list(...)
  ngrds = length(grds)
  if (ngrds < 1)
    stop("no arguments supplied")
  if (ngrds == 1)
    return(grds[[1]])
  ## verify matching topology:
  sapply(grds[2:ngrds], function(x) stop.ifnot.equal(x, grds[[1]]))
  gr = grds[[1]]
  gr@data = do.call("cbind", lapply(grds, function(x) x@data))
  ##for (i in 2:ngrds)
  ##	gr@data = cbind(gr@data, grds[[i]]@data)
  if (is(gr, "RFspatialGridDataFrame"))
    proj4string(gr) = CRS(proj4string(grds[[1]]))
  gr
}

cbind_RFspPoints = function(...) {  ##copied from sp package
  stop.ifnot.equal = function(a, b) {
    res = all.equal(a@coords, b@coords)
    if (!is.logical(res) || !res)
      stop("coords are not equal")
  }
  grds = list(...)
  ngrds = length(grds)
  if (ngrds < 1)
    stop("no arguments supplied")
  if (ngrds == 1)
    return(grds[[1]])
  ## verify matching topology:
  sapply(grds[2:ngrds], function(x) stop.ifnot.equal(x, grds[[1]]))
  gr = grds[[1]]
  gr@data = do.call("cbind", lapply(grds, function(x) x@data))
  ##for (i in 2:ngrds)
  ##	gr@data = cbind(gr@data, grds[[i]]@data)
  gr
}

cbind.RFspatialGridDataFrame <- function(...)
  cbind_RFsp(...)
cbind.RFgridDataFrame <- function(...)
  cbind_RFsp(...)

cbind.RFspatialPointsDataFrame <- function(...)
  cbind_RFspPoints(...)
cbind.RFpointsDataFrame <- function(...)
  cbind_RFspPoints(...)



## convert 'RFsp' objects to conventional format of 'RFsimulate',
## i.e. data is an array and x a matrix of coordinates or gridtriple defs.

spatialGridObject2conventional <- function(obj) {
  timespacedim <- length(obj@grid@cells.dim)
  data <- as.matrix(obj@data)

  has.variance <- !is.null(obj@.RFparams$has.variance) && obj@.RFparams$has.variance
  
  dim(data) <- NULL
  dim(data) <- c(obj@grid@cells.dim,
                 obj@.RFparams$vdim,
                 obj@.RFparams$n + has.variance)
  if (timespacedim > 1)
    data <- reflection(data, 2, drop=FALSE)
  ## re-ordering of 2nd space dimension since in sp objects, the 2nd dimension
  ## is in decreasing order
  
  perm <- c(timespacedim+1, 1:timespacedim, timespacedim+2) 
  data <- aperm(data, perm=perm)
  ## new order of dimensions: vdim, space-time-dims, n
  
  is.dim <- dim(data) != 1
  if (sum(is.dim) > 1)
    dim(data) <- dim(data)[is.dim] # drop dimensions length 1
  else
    dim(data) <- NULL
  
  x <- rbind(obj@grid@cellcentre.offset,
             obj@grid@cellsize,
             obj@grid@cells.dim)

  if (dimensions(obj)==1 ||
      !("coords.T1" %in% names(obj@grid@cellcentre.offset)))
    T <- NULL
  else {
    idxT1 <- which("coords.T1" == names(obj@grid@cellcentre.offset))
    T <- x[,  idxT1]
    x <- x[, -idxT1, drop=FALSE]
  }

  .RFparams <- obj@.RFparams
  
  return(list(data=data, x=x, T=T, .RFparams=.RFparams))
}

spatialPointsObject2conventional <- function(obj) {
  data <- as.matrix(obj@data)
  has.variance <- !is.null(obj@.RFparams$has.variance) && obj@.RFparams$has.variance
  dim(data) <- NULL
  dim(data) <- c(nrow(obj@data), obj@.RFparams$vdim, obj@.RFparams$n + has.variance)
  perm <- c(2,1,3)
  data <- aperm(data, perm=perm)
  is.dim <- dim(data) != 1
  if (sum(is.dim) > 1)
    dim(data) <- dim(data)[is.dim]
  else
    dim(data) <- NULL
  
  x <- obj@coords
  dimnames(x) <- NULL

  if (dimensions(obj)==1 || !("coords.T1" %in% colnames(obj@coords)))
    T <- NULL
  else {
    idxT1 <- which("coords.T1" == colnames(obj@coords))
    T <- points2grid(RFpointsDataFrame(coords=unique(x[, idxT1]),
                                       data=double(length(unique(x[, idxT1]))),
                                       RFparams=obj@.RFparams))
    dimdata <- dim(data)
    if (obj@.RFparams$vdim==1)
      dim(data) <- c(dimdata[1]/T@cells.dim, T@cells.dim, dimdata[-1])
    else
      dim(data) <- c(dimdata[1], dimdata[2]/T@cells.dim,
                     T@cells.dim, dimdata[-c(1,2)])
    x <- x[1:(nrow(x)/T@cells.dim), -idxT1, drop=FALSE]
    T <- c(T@cellcentre.offset, T@cellsize, T@cells.dim)
   }
  
  return(list(data=data, x=x, T=T, .RFparams=obj@.RFparams))
}

setGeneric(name = "RFspDataFrame2conventional", 
           def = function(obj)
           standardGeneric("RFspDataFrame2conventional") )

setMethod("RFspDataFrame2conventional", signature=c("RFspatialGridDataFrame"),
          function(obj)
          spatialGridObject2conventional(obj))
setMethod("RFspDataFrame2conventional", signature=c("RFgridDataFrame"),
          function(obj)
          spatialGridObject2conventional(obj))

setMethod("RFspDataFrame2conventional",
          signature=c("RFspatialPointsDataFrame"),
          function(obj)
          spatialPointsObject2conventional(obj))
setMethod("RFspDataFrame2conventional",
          signature=c("RFpointsDataFrame"),
          function(obj)
          spatialPointsObject2conventional(obj))


rfspDataFrame2dataArray_generic <- function(obj) {
  has.variance <- !is.null(obj@.RFparams$has.variance) && obj@.RFparams$has.variance
  arr <- unlist(obj@data)
  names(arr) <- NULL
  grid.length.v <- obj@grid@cells.dim
  dim(arr) <- c(grid.length.v,         # spatial and time dim
                obj@.RFparams$vdim,    # vdim
                obj@.RFparams$n + has.variance)       # repetitions
  timespacedim <- length(grid.length.v)
  if (timespacedim > 1)
    arr <- reflection(arr, 2, drop=FALSE)
    ## re-ordering of 2nd space dimension since in sp objects, the 2nd dimension
    ## is in decreasing order
  return(arr)
}      

setGeneric(name = "RFspDataFrame2dataArray", 
           def = function(obj) standardGeneric("RFspDataFrame2dataArray"))
setMethod("RFspDataFrame2dataArray", signature=c("RFspatialGridDataFrame"),
          function(obj) rfspDataFrame2dataArray_generic(obj))
          ## used in plot method
setMethod("RFspDataFrame2dataArray", signature=c("RFgridDataFrame"),
          function(obj) rfspDataFrame2dataArray_generic(obj))



##  convert GridTopology or 3-row matrix to list of vectors
setGeneric(name = "GridTopology2gridVectors", 
           def = function(grid)
           standardGeneric("GridTopology2gridVectors"))
setMethod("GridTopology2gridVectors",
          signature=c("GridTopology"),
          function(grid) {
            len <- length(grid@cells.dim)
            x <- list()
            for (i in 1:len)
              x[[i]] <- seq(from   = grid@cellcentre.offset[i],
                            to     = grid@cellcentre.offset[i] +
                                     ((grid@cells.dim[i]-1) * grid@cellsize[i]),
                            length = grid@cells.dim[i])
            return(x)
          })
setMethod("GridTopology2gridVectors",
          signature=c("matrix"),
          function(grid) {
            stopifnot(nrow(grid)==3)
            x <- list()
            for (i in 1:ncol(grid))
              x[[i]] <- seq(from  = grid[1,i],
                            to    = grid[1,i] + (grid[3,i]-1) * grid[2,i],
                            length= grid[3,i])
            return(x)
          })
##  convert GridTopology to 3-row matrix
as.matrix.GridTopology <- function(x, ...)
  return(rbind(x@cellcentre.offset, x@cellsize, x@cells.dim))



## extract or calculate! coordinates; 
                                                
setMethod(f = "coordinates", signature="RFpointsDataFrame",
          definition=function(obj) return(obj@coords) )
setMethod(f = "coordinates", signature="RFgridDataFrame",
          definition=function(obj) coordinates(as(obj, "RFpointsDataFrame")) )


extract.kriging.variance <- function(x){
  stopifnot(is(x, "RFsp"))
  has.variance <- !is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance
  if (!has.variance) return(NULL)
  return(x[(x@.RFparams$n*x@.RFparams$vdim + 1):(ncol(x@data))])
}

  
setGeneric(name = "variance", 
           def = function(obj) standardGeneric("variance"))
setMethod(f = "variance", signature="RFsp",
          definition=function(obj) extract.kriging.variance(obj))




## conventional 'RFsimulate' output to 'RFsp' class

print.RFpointsDataFrame = function(x, ..., digits = 6) {
  df = data.frame(coordinates=signif(coordinates(x), digits), x@data)
  row.names(df) = row.names(x@data)
  cat("Object of class RFpointsDataFrame\n")
  print(df, ...)
}
setMethod(f="show", signature="RFpointsDataFrame",
          definition=function(object) print.RFpointsDataFrame(object))


print.RFgridDataFrame = function(x, ...) {
	cat("Object of class RFgridDataFrame\n")
        cat("Grid topology:\n")
	print(data.frame(cellcentre.offset=x@grid@cellcentre.offset,
                         cellsize=x@grid@cellsize,
                         cells.dim=x@grid@cells.dim, row.names=""))
        cat("Points:\n")
        print(data.frame(coordinates=coordinates(x), x@data))
	cat("Data summary:\n")
        if (ncol(x@data) > 1)
          sobj = summary(x@data)
        else sobj = summary(x@data[[1]])
	#print(sobj)
	invisible(x)
}
setMethod(f="show", signature="RFgridDataFrame", 
	  definition=function(object) print.RFgridDataFrame(object))


print.RFspatialPointsDataFrame = function(x, ...) {
  cat("Object of class RFspatialPointsDataFrame\n")
  show(as(x, "SpatialPointsDataFrame"))
}
setMethod(f="show", signature="RFspatialPointsDataFrame", 
	  definition=function(object) print.RFspatialPointsDataFrame(object))


print.RFspatialGridDataFrame = function(x, ...) {
  cat("Object of class RFspatialGridDataFrame\n")
  show(as(x, "SpatialGridDataFrame"))
}
setMethod(f="show", signature="RFspatialGridDataFrame", 
	  definition=function(object) print.RFspatialGridDataFrame(object))


## extend methods for 'Spatial' to class 'RFsp'
if (!isGeneric("isGridded"))
  setGeneric(name = "isGridded", 
             def = function(obj)
             standardGeneric("isGridded") )

setMethod(f="isGridded", signature="RFsp", 
	  definition=function(obj)
          (is(obj, "RFgridDataFrame") || is(obj, "RFspatialGridDataFrame")))
setMethod(f="dimensions", signature="RFsp", 
	  definition=function(obj)
          if (is(obj, "RFgridDataFrame") || is(obj, "RFpointsDataFrame"))
          return(1)
          else (getMethod("dimensions", "Spatial")@.Data)(obj) )


summary.RFsp = function(object, ...) {
  if (is(object, "Spatial"))
    return((getMethod("summary", "Spatial")@.Data)(object))
  
  obj = list()
  obj[["class"]] = class(object)
            #obj[["bbox"]] = bbox(object)
            #obj[["is.projected"]] = is.projected(object)
            #obj[["proj4string"]] = object@proj4string@projargs
  if (is(object, "RFpointsDataFrame")) 
    obj[["npoints"]] = nrow(object@coords)
  if (is(object, "RFgridDataFrame")) {
    gr <- object@grid
    obj[["grid"]] =
      data.frame(cellcentre.offset = gr@cellcentre.offset, 
                 cellsize = gr@cellsize,
                 cells.dim = gr@cells.dim)
  }
  if ("data" %in% slotNames(object))
    if (nrow(object@data) > 0)
      if (ncol(object@data) > 1) 
        obj[["data"]] = summary(object@data)
      else obj[["data"]] = summary(object@data[[1]])
                                        #class(obj) = "summary.Spatial"
  return(obj)
}

          

## plot methods
## for 1-dim coordinates for grid and points
## for 2-or more dimensional coordinates for grid only

setMethod(f="plot", signature(x="RFgridDataFrame", y="missing"),
	  definition=function(x, y, nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            ...)
          plotRFgridDataFrame(x, nmax=nmax,
                              plot.variance=plot.variance, ...))


setMethod(f="plot", signature(x="RFpointsDataFrame", y="missing"),
	  definition=function(x, y, nmax = 6, sorted = FALSE,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            ...)
          plotRFpointsDataFrame(x, nmax=nmax, sorted=sorted,
                                plot.variance=plot.variance, ...))


setMethod(f="plot", signature(x="RFspatialGridDataFrame", y="missing"),
	  definition=function(
            x, MARGIN = c(1, 2), MARGIN.slices = NULL,
            n.slices = if (!missing(MARGIN.slices) && !is.null(MARGIN.slices))
                           10 else 1,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables, # = 1:vdim,
            zlim,
            legend=TRUE, ...)
          plotRFspatialGridDataFrame(x=x, MARGIN=MARGIN,
                                     MARGIN.slices=MARGIN.slices,
                                     n.slices=n.slices, nmax=nmax,
                                     plot.variance = plot.variance,
                                     select=select.variables,
                                     zlim=zlim,
                                     legend=legend,
                                       ...))
setMethod(f="plot",
          signature(x="RFspatialGridDataFrame", y="RFspatialPointsDataFrame"),
	  definition=function(
            x, y, MARGIN = c(1, 2), MARGIN.slices = NULL,
            n.slices = if (!missing(MARGIN.slices) && !is.null(MARGIN.slices))
                       10 else 1,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables,
            zlim,
            legend=TRUE, ...)
          plotRFspatialGridDataFrame(x=x, y=y, MARGIN=MARGIN,
                                      MARGIN.slices=MARGIN.slices,
                                      n.slices=n.slices, nmax=nmax,
                                      plot.variance = plot.variance,
                                      select=select.variables,
                                      zlim=zlim,
                                      legend=legend,
                                     ...))

setMethod(f="plot", signature(x="RFspatialPointsDataFrame", y="missing"),
	  definition=function(x, y, MARGIN = c(1, 2), nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables,
            zlim,
            legend=TRUE,...)
          plotRFspatialPointsDataFrame(x=x, MARGIN=MARGIN, nmax=nmax, 
                                       plot.variance = plot.variance,
                                       select=select.variables,
                                       zlim=zlim,
                                       legend=legend, ...))

setMethod(f="plot", signature(x="RFspatialPointsDataFrame",
            y="RFspatialPointsDataFrame"),
	  definition=function(x, y, MARGIN = c(1, 2), nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables,
            zlim,
            legend=TRUE, ...)
          plotRFspatialPointsDataFrame(x=x, y=y, MARGIN=MARGIN, nmax=nmax, 
                                       plot.variance = plot.variance,
                                       select=select.variables,
                                       zlim=zlim,
                                       legend=legend, ...))

errMsgNoPlotAvailable <- function(x, y)
  warning(paste("no plot method available for signature c(",
                class(x), ",", class(y), ")"))

setMethod(f="plot", signature(x="RFspatialPointsDataFrame",
            y="RFspatialGridDataFrame"),
	  definition=function(x, y) {
            errMsgNoPlotAvailable(x, y)
            return(invisible(NULL))
          })
 
plotRFgridDataFrame <- function(x, y, nmax, plot.variance, ...) {
  #if (!is(x, "RFgridDataFrame"))
  #  stop("method only for objects of class 'RFgridDataFrame'")
  x <- as(x, "RFpointsDataFrame")
  plotRFpointsDataFrame(x, y, nmax=nmax, sorted=TRUE,
                        plot.variance=plot.variance, ...)
}


plotRFpointsDataFrame <-
  function(x, y, nmax=6, sorted=FALSE,
           plot.variance=plot.variance, ...)
{
  if (!is(x, "RFpointsDataFrame"))
    stop("method only for objects of class 'RFpointsDataFrame'")

  graphics <- RFoptions()$graphics
  
  has.variance <- !is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance
  if (!has.variance) plot.variance <- FALSE
   
  x.coords <- as.vector(x@coords)
  if (!sorted){
    ord <- order(x.coords)
    x.coords <- x.coords[ord]
    x@data <- x@data[ord, , drop=FALSE]
  }
  
  n <- min(x@.RFparams$n, nmax) + plot.variance
  vdim <- x@.RFparams$vdim
  always.close <- n > 1 || graphics$always_close_screen
##  Print("A", always.close)
  
  nc <- ncol(x@data)
  if (nc < n*vdim)
    if (n==1) vdim <- nc else if (vdim==1) n <- nc else {
      stop("ncol(x@data) does not match 'x@.RFparams'; change 'x@.RFparams'")
    }

  dots <- list(...)
  dotnames <- names(dots)
  if ("bg" %in% dotnames) {
    par(bg=dots$bg)
    dots$bg <- NULL
  }

  dummy <- dimnames(x@coords)[[2]][1]
  lab <- xylabs(if (is.null(dummy)) "" else dummy, "", units=x@.RFparams$coord.units)
 
  if (!("xlab" %in% dotnames)) dots$xlab <- lab$x
  if (!("type" %in% dotnames)) dots$type <- "l"

  make.small.mar <- ("xlab" %in% dotnames &&
                     is.null(dots$xlab) && is.null(dots$ylab))
  
  ## variable names
  if (!is.null(names(x@data)) && all(nchar(names(x@data))>0))
    names.vdim <- unlist(lapply(strsplit(names(x@data)[1:vdim], ".n"),
                                FUN=function(li) li[[1]]))
  else {
    names.vdim <- paste("variable", 1:vdim)
    names(x@data) <- names.vdim
  }

  if (n>1){
    ylab.vec <- c(paste("realization ", 1:(n-plot.variance), sep=""),
                  if (plot.variance) "kriging variance")
  } else {
    ylab.vec <- if (vdim==1) colnames(x@data) else ""
  }

  if ("ylab" %in% dotnames) {
    if (!is.null(dots$ylab))
      ylab.vec[1:length(ylab.vec)] <- dots$ylab
    dots$ylab <- NULL
  }

  col <- 1:vdim
  if ("col" %in% dotnames) {
    if (!is.null(dots$col))
      col[1:length(col)] <- dots$col
    dots$col <- NULL
  }

  
  ArrangeDevice(graphics, c(1, n))
  split.screen(c(n,1))
  
#  if (always.close) {
#    close.screen(all.screens=TRUE)
#    par(mfrow=c(1,1))
#    split.screen(c(n,1))
#  }
                
  for (i in 1:n){
    screen(i)
    if (make.small.mar)
      par(oma=c(3,0,1,1)+.1, mar=c(0,3,0,0))
    else
      par(oma=c(4,0,1,1)+.1, mar=c(0,4,0,0))
    ylab <- ylab.vec[i]
    
    if (tmp.idx <- (plot.variance && i==n)){
      i <- x@.RFparams$n + plot.variance
    }

    do.call(plot,
            args=c(dots, list(
              x=x.coords, y=x@data[ , vdim*(i-1)+1],
              xaxt="n", yaxt="n", ylab=ylab, col=col[1]))
            )
    axis(2)
    if (tmp.idx) i <- n
    
    if(i==n){
      axis(1, outer=always.close)
      title(xlab=dots$xlab, outer=TRUE) # always.close) 
    }
    else axis(1, labels=FALSE)
    for (j in 1:vdim){
      if (j==1) next
      do.call(what="points", quote=TRUE,
              args=c(dots, list(
                x=x.coords, y=x@data[ , vdim*(i-1)+j], col=col[j]))
              )
      if (i==1) {
        if (!TRUE || vdim > 1) {
          legend("topright", col=col, lty=1, legend = c(names.vdim))
        }
      }
    }
  }
  if (always.close) close.screen(all.screens=TRUE)
}


### nur fuer Testzwecke
plotRFspatialGridDataFrameX <- function(x, y, MARGIN, MARGIN.slices,
                                        n.slices, nmax, plot.variance,
                                        select, zlim, legend,
                                        coord.units=NULL, variab.units=NULL,
                                        ...) {
  if (file.exists("/home/schlather/R/RF/svn/makefile")) {

    args <- list(x=x,
                 y=if (!missing(y)) y,
                 MARGIN=MARGIN, MARGIN.slices=MARGIN.slices,
                 n.slices=n.slices, nmax=nmax,
                 plot.variance = if (!missing(plot.variance)) plot.variance,
                 select = if (!missing(select)) select,
                 zlim = if (!missing(zlim)) zlim,
                 legend=if (!missing(legend)) legend,
                 coord.units=coord.units, variab.units, variab.units, ...)
    args <- args[!sapply(args, is.null)]
 #   Print(args)
    do.call("plotRFspatialGridDataFrame", args=args, envir=.GlobalEnv)
    
  } else {
    plotRFspatialGridDataFrame(x, y, MARGIN, MARGIN.slices,
                               n.slices, nmax, plot.variance,
                               select, zlim, legend,
                               coord.units=coord.units,
                               variab.units=variab.units, ...)
  }
}
