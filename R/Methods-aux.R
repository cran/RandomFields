
reflection <- function(data, orth, drop=FALSE)
  ##IMPORPANT NOTE! DO NOT CHANGE THE VARIABLE NAMES IN THIS SIGNATURE
  ## why ???
  ## since the variable data is pasted by its name
{
  d <- dim(data)
  return(do.call("[", c(list(data), rep(TRUE, orth-1), list(d[orth]:1),
                        rep(TRUE, length(d) - orth), drop=drop)))
}

AddUnits <- function(params) {
  ## see also empvario.R and fitgauss.R, if changed
  coords <- RFoptions()$general
  return(c(params, list(coord.units=coords$coordunits,
                        variab.units=coords$varunits)))
}

compareGridBooleans <- function(grid, gridtmp) {
  if (!missing(grid) && length(grid)>0 && grid!=gridtmp)
    message(paste("you specified grid=", as.character(grid),
                  " but isGridded(data)=", as.character(gridtmp),
                  ";  grid is set to ", as.character(gridtmp), sep=""))
}

isSpObj <- function(x)
  (is(x, "SpatialGridDataFrame") || is(x, "SpatialPointsDataFrame")) &&
  !is(x, "RFsp")


sp2RF <- function(sp, param=list(n=1, vdim=1)) {
  class(sp) <- paste("RF", tolower(substr(class(sp), 1, 1)),
                       substring(class(sp), 2),  sep="")
  sp@.RFparams <- AddUnits(param)
  validObject(sp)
  return(sp)
}

convert2GridTopology <- function(grid){
  if (!is(grid, "GridTopology")) {
    if (is.null(dim(grid)))
      grid <- matrix(grid, ncol=1)
    stopifnot(nrow(grid)==3)
    grid <- sp::GridTopology(cellcentre.offset=grid[1,],
                             cellsize=grid[2,],
                             cells.dim=grid[3,])
  }
  return(grid)
}
     


## Generate Objects ########################################################

RFspatialGridDataFrame <- function(grid, data,
                                   proj4string = sp::CRS(as.character(NA)),
                                   RFparams=list(n=1, vdim=1)) {
  grid <- convert2GridTopology(grid)
  tmp <- sp::SpatialGridDataFrame(grid=grid,
                                  data = if (is.data.frame(data)) data else
                                  data.frame(data),
                                  proj4string=proj4string)
  return(sp2RF(tmp, RFparams))
#  tmp <- as(tmp, "RFspatialGridDataFrame")
#  tmp@.RFparams <- AddUnits(RFparams)
#  validObject(tmp)
#  return(tmp)
}

RFspatialPointsDataFrame <- function(coords, data, coords.nrs = numeric(0),
                                     proj4string = sp::CRS(as.character(NA)), 
                                     match.ID = TRUE, bbox = NULL,
                                     coord.units = NULL,
                                     variab.units = NULL,
                                     RFparams=list(n=1, vdim=1)) {
  if (is.null(bbox)) {
    bbox <- t(apply(coords, 2, range))
    colnames(bbox) <- c("min", "max")
  }
  tmp <- sp::SpatialPointsDataFrame(coords=coords,
                                    data=if (is.data.frame(data)) data else
                                    data.frame(data),
                                    coords.nrs=coords.nrs,
                                    proj4string=proj4string, 
                                    match.ID=match.ID, bbox=bbox)
#  if (!is.null(dimnames(coords)))
#    dimnames(tmp@coords) <- dimnames(coords)
#
#
  return(sp2RF(tmp, RFparams))
#  str(tmp)
#  Print(class(tmp))
#  class(tmp) <-  "RFspatialPointsDataFrame"
# tmp <- as(tmp, "RFspatialPointsDataFrame")
  
#  tmp@.RFparams <- AddUnits(RFparams)
# Print(class(tmp))
#  str(tmp); lll
#   validObject(tmp)
#  return(tmp)
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




cbind_RFsp <- function(...) {  ##copied from sp package
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
  gr@data = do.call(base::cbind, lapply(grds, function(x) x@data))
  ##for (i in 2:ngrds)
  ##	gr@data = cbind(gr@data, grds[[i]]@data)
  if (is(gr, "RFspatialGridDataFrame"))
    sp::proj4string(gr) = sp::CRS(sp::proj4string(grds[[1]]))
  gr
}

cbind_RFspPoints <- function(...) {  ##copied from sp package
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
  gr@data = do.call(base::cbind, lapply(grds, function(x) x@data))
  ##for (i in 2:ngrds)
  ##	gr@data = cbind(gr@data, grds[[i]]@data)
  gr
}



extract.names <- function(names) {
  if (length(names) == 1) return(as.vector(names))
  nr <- strsplit(names[,1], ".")
  if (any(sapply(nr, length) != 2)) nr <- names[,1]
  else nr <- sapply(nr, function(x) x[1])

  nc <- strsplit(names[1,], ".")
  if (any(sapply(nc, length) != 2)) nc <- names[1,]
  else nc <- sapply(nc, function(x) x[1])

  return(list(nr, nc))
}



spatialGridObject2conventional <- function(obj, data.frame=FALSE) {
  timespacedim <- length(obj@grid@cells.dim)
  data <- as.matrix(obj@data)
  names <- colnames(data)
  
  has.variance <- !is.null(obj@.RFparams$has.variance) &&
    obj@.RFparams$has.variance
  dim(data) <- NULL
  vdimn <- c(obj@.RFparams$vdim, obj@.RFparams$n + has.variance)

  dim(data) <- c(obj@grid@cells.dim, vdimn)
#  Print(names)
  
  if (timespacedim > 1) data <- reflection(data, 2, drop=FALSE)
  ## re-ordering of 2nd space dimension since in sp objects, the 2nd dimension
  ## is in decreasing order


  if (data.frame) {
    dim(data) <- c(prod(obj@grid@cells.dim), prod(vdimn))
    colnames(data) <- names
    return(as.data.frame(data))
  }
  
  dim(names) <- vdimn
  vdim_close_together <- FALSE
  if (vdim_close_together) {
    perm <- c(timespacedim+1, 1:timespacedim, timespacedim+2) 
    data <- aperm(data, perm=perm)
    names <- aperm(names, perm[-1]) ### ?????
  }
  ## new order of dimensions: vdim, space-time-dims, n

  is.dim <- dim(data) != 1
  if (sum(is.dim) > 1) {
    dim(data) <- dim(data)[is.dim] # drop dimensions length 1
    l <- list()
    l[length(obj@grid@cells.dim) + (1:2)] <- extract.names(names)
    dimnames(data) <- l[is.dim]
  } else {
    dim(data) <- NULL
    #names(data) <- names
  }

  x <- rbind(obj@grid@cellcentre.offset,
             obj@grid@cellsize,
             obj@grid@cells.dim)

 # Print(obj, "TTTT", is(obj, "RFsp"))
  
  if (dimensions(obj)==1 ||
      !("coords.T1" %in% names(obj@grid@cellcentre.offset)))
    T <- NULL
  else {
    idxT1 <- which("coords.T1" == names(obj@grid@cellcentre.offset))
    T <- x[,  idxT1]
    x <- x[, -idxT1, drop=FALSE]
  }

  .RFparams <- obj@.RFparams
  
  return(list(data=data, x=x, T=T, .RFparams=.RFparams, .names=names))
}

spatialPointsObject2conventional <- function(obj) {
  data <- as.matrix(obj@data)
  Enames <- names <- colnames(data)
  
  has.variance <-
    !is.null(obj@.RFparams$has.variance) && obj@.RFparams$has.variance
  dim(data) <- NULL
  vdimn <- c(obj@.RFparams$vdim, obj@.RFparams$n + has.variance)
  dim(data) <- c(nrow(obj@data), obj@.RFparams$vdim, vdimn)
  
  dim(Enames) <- vdimn
  Enames <- extract.names(Enames)
  vdim_close_together <- FALSE
  if (vdim_close_together) {
    perm <- c(2,1,3)
    data <- aperm(data, perm=perm)
    Enames <- aperm(Enames, perm[-1]) ### ?????
  }

  x <- obj@coords
  dimnames(x) <- NULL
  if (dimensions(obj)==1 || !("coords.T1" %in% colnames(obj@coords))) {
    T <- NULL
    is.dim <- dim(data) != 1
    if (sum(is.dim) > 1) {    
      dim(data) <- dim(data)[is.dim] # drop dimensions length 1
      dimnames(data) <- c(list(NULL), Enames)[is.dim]
    } else {
      dim(data) <- NULL
      ##names(data) <- names
    }  
  } else {
    idxT1 <- which("coords.T1" == colnames(obj@coords))
    T <- sp::points2grid(RFpointsDataFrame(coords=unique(x[, idxT1]),
                                           data=double(length(unique(x[,idxT1]))),
                                           RFparams=obj@.RFparams))
    dimdata <- dim(data)
    if (obj@.RFparams$vdim==1) {
      dim(data) <- c(dimdata[1]/T@cells.dim, T@cells.dim, dimdata[-1:-2])
      dimnames(data) <- list(NULL,
                             paste("T", 1:T@cells.dim, sep=""),
                             Enames[[2]])
    } else {
      dim(data) <- c(dimdata[1], dimdata[2]/T@cells.dim,
                     T@cells.dim, dimdata[-1])
      dimnames(data) <- list(NULL,
                             paste("T", 1:T@cells.dim, sep=""),
                             Enames[[1]], Enames[[2]])
     
    }
    x <- x[1:(nrow(x)/T@cells.dim), -idxT1, drop=FALSE]
    T <- c(T@cellcentre.offset, T@cellsize, T@cells.dim)
  }
  return(list(data=data, x=x, T=T, .RFparams=obj@.RFparams))
}


rfspDataFrame2conventional <- function(obj) {
  if (is(obj, "RFspatialPointsDataFrame") || is(obj, "RFpointsDataFrame"))
    return(spatialPointsObject2conventional(obj))
  else if (is(obj, "RFspatialGridDataFrame") || is(obj, "RFgridDataFrame"))
    return(spatialGridObject2conventional(obj))
  else stop("unknown class in 'RFspDataFrame2conventional'")
}
