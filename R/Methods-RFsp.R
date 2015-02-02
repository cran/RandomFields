
## Coerce Objects #########################################################

as.data.frame.RFpointsDataFrame <-
  as.data.frame.RFspatialPointsDataFrame <- function(x, ...) {
  #str(x); kkkk
  cbind(x@data, x@coords)
}
as.data.frame.RFgridDataFrame <-
  as.data.frame.RFspatialGridDataFrame <- function(x, ...) 
  spatialGridObject2conventional(x, TRUE)

setAs("RFspatialPointsDataFrame", "data.frame",
      function(from, to) from@data
)
setAs("RFspatialGridDataFrame", "data.frame",
      function(from, to) spatialGridObject2conventional(from, TRUE)$data
)


setAs("SpatialPointsDataFrame", "RFspatialPointsDataFrame",
      function(from, to) sp2RF(from))
setAs("SpatialGridDataFrame", "RFspatialGridDataFrame",
      function(from, to) sp2RF(from))

setAs("RFspatialGridDataFrame", "RFspatialPointsDataFrame",
      function(from, to) RFspatialPointsDataFrame(coords=coordinates(from),
                                                  data=from@data,
                                                  RFparams=from@.RFparams))
setAs("RFspatialPointsDataFrame", "RFspatialGridDataFrame",
      function(from, to) RFspatialGridDataFrame(grid=sp::points2grid(from),
                                                data=from@data,
                                                RFparams=from@.RFparams))

setAs("RFpointsDataFrame", "RFgridDataFrame",
      function(from, to) {
        tmp.coord <- sp::SpatialPoints(cbind(from@coords, 1:length(from@coords)))
        tmp.grid <- sp::points2grid(tmp.coord)
        grid <- convert2GridTopology(tmp.grid[,1, drop=FALSE])
        RFgridDataFrame(data=from@data, grid=grid, RFparams=from@.RFparams)
      })
setAs("RFgridDataFrame", "RFpointsDataFrame", 
      function(from, to) {
        coords <- GridTopology2gridVectors(from@grid)[[1]]
        RFpointsDataFrame(data=from@data, coords=coords,
                          RFparams=from@.RFparams)
      })



## methods 'as.matrix' and 'cbind' as in the 'sp' package

as.matrix.RFgridDataFrame <- 
 as.matrix.RFspatialGridDataFrame <- function(x, ...) {
  z <- as.array.RFgridDataFrame(x, ...) 
  if (length(dim(z)) > 2)
    stop("the data set cannot be turned into a matrix. Try 'as.array'")
  z
}

as.array.RFgridDataFrame <-
  as.array.RFspatialGridDataFrame <- function(x, ...) {
    z <- as.matrix(x@data)    
    dim <- c(x@grid@cells.dim, x@.RFparams$vdim, x@.RFparams$n)
    dim <- dim[dim > 1]
    if (length(dim) == 1) return(as.vector(z))
    dim(z) <- dim
    if ((l <- length(x@grid@cells.dim)) > 1) {
      idx <- as.list(rep(TRUE, length(dim)))
      idx[[2]] <- (dim[2]) : 1
      #Print(idx, c(list(z), drop=FALSE, idx))
      do.call("[", c(list(z), drop=FALSE, idx))
    } else z
  }


as.vector.RFgridDataFrame <-
  as.vector.RFspatialGridDataFrame <- function(x, ...) {
    as.vector(as.array.RFgridDataFrame(x))
  }

as.matrix.RFpointsDataFrame <-
  as.matrix.RFspatialPointsDataFrame <- function(x, ...) {
  z <- as.array.RFspatialPointsDataFrame(x, ...) 
  if (length(dim(z)) > 2)
    stop("the data set cannot be turned into a matrix. Try 'as.array'")
  z
}
as.array.RFpointsDataFrame <-
  as.array.RFspatialPointsDataFrame <- function(x, ...) {
  z <- as.matrix(x@data)
  dim <- c(x@.RFparams$vdim, x@.RFparams$n)
  if (length(x@.RFparams$T) == 3) dim <- c(x@.RFparams$T[3], dim)
  dim <- c(length(z) / prod(dim), dim)
  dim <- dim[dim > 1]
  if (length(dim) == 1) return(as.vector(z))
  dim(z) <- dim
  z
}
as.vector.RFpointsDataFrame <-
  as.vector.RFspatialPointsDataFrame <- function(x, ...) {
  as.vector(as.array.RFspatialPointsDataFrame(x))
}



#setAs("RFgridDataFrame", "matrix", gridtomatrix)
#setAs("RFspatialGridDataFrame", "matrix", gridtomatrix)


cbind.RFspatialGridDataFrame <- function(...)
  cbind_RFsp(...)
cbind.RFgridDataFrame <- function(...)
  cbind_RFsp(...)

cbind.RFspatialPointsDataFrame <- function(...)
  cbind_RFspPoints(...)
cbind.RFpointsDataFrame <- function(...)
  cbind_RFspPoints(...)

range.RFspatialGridDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)
range.RFgridDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)
range.RFspatialPointsDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)
range.RFpointsDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)


hist.RFspatialGridDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)
hist.RFgridDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)
hist.RFspatialPointsDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)
hist.RFpointsDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)


## convert 'RFsp' objects to conventional format of 'RFsimulate',
## i.e. data is an array and x a matrix of coordinates or gridtriple defs.

setGeneric(name = "RFspDataFrame2conventional", 
           def = function(obj) standardGeneric("RFspDataFrame2conventional"))
setMethod("RFspDataFrame2conventional",
          signature=c("RFspatialGridDataFrame"),
          function(obj) spatialGridObject2conventional(obj))
setMethod("RFspDataFrame2conventional", signature=c("RFgridDataFrame"),
          function(obj) spatialGridObject2conventional(obj))
setMethod("RFspDataFrame2conventional",
          signature=c("RFspatialPointsDataFrame"),
          function(obj) spatialPointsObject2conventional(obj))
setMethod("RFspDataFrame2conventional", signature=c("RFpointsDataFrame"),
          function(obj) spatialPointsObject2conventional(obj))


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
           def = function(grid) standardGeneric("GridTopology2gridVectors"))
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

summary.RFpointsDataFrame <- function(object, digits = 6, ...) {
  df = data.frame(coordinates=signif(coordinates(object), digits), object@data)
  row.names(df) = row.names(object@data)
  class(df) <- "summary.RFpointsDataFrame"
  df
}
print.summary.RFpointsDataFrame <- function(x, ...)
  print.data.frame(x, ...)#

print.RFpointsDataFrame <- function(x, ...) {
  sx <- summary.RFpointsDataFrame(x, ...)
  cat("Object of class RFpointsDataFrame\n")
  print.summary.RFpointsDataFrame(sx) #
  invisible(sx)
}

setMethod(f="show", signature="RFpointsDataFrame",
          definition=function(object) print.RFpointsDataFrame(object))

summary.RFgridDataFrame <- function(object, ...) {
  if (ncol(object@data) > 1) summary(object@data) else summary(object@data[[1]])
} 
print.RFgridDataFrame <- function(x, ...) {
  cat("Object of class RFgridDataFrame\n")
  cat("Grid topology:\n")
  print.data.frame(data.frame(cellcentre.offset=x@grid@cellcentre.offset, #
                   cellsize=x@grid@cellsize,
                   cells.dim=x@grid@cells.dim, row.names=""))
  cat("Points:\n")
  print.data.frame(data.frame(coordinates=coordinates(x), x@data)) #
  cat("Data summary:\n")
  summary.RFgridDataFrame(x, ...)
  invisible(x)
}
setMethod(f="show", signature="RFgridDataFrame", 
	  definition=function(object) print.RFgridDataFrame(object))

summary.RFspatialPointsDataFrame <- function(object, ...) summary(object@data)
print.RFspatialPointsDataFrame <- function(x, ...){
  if (!hasArg("silent") || !list(...)$silent) {
    cat("Object of class 'RFspatialPointsDataFrame'\n")
    str(x@data, no.list=TRUE, give.attr=FALSE)#
  }
  invisible(x@data)
}
setMethod(f="show", signature="RFspatialPointsDataFrame", 
	  definition=function(object) print.RFspatialPointsDataFrame(object))


summary.RFspatialGridDataFrame <- function(object, ...) summary(object@data)
print.RFspatialGridDataFrame <- function(x,...) {
 if (!hasArg("silent") || !list(...)$silent) {
    cat("Object of class 'RFspatialGridDataFrame'\n")
    utils::str(x@data, no.list=TRUE, give.attr=FALSE)#
  }
  invisible(x@data)
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
	  definition=function(obj) {
            (getMethod("dimensions", "Spatial")@.Data)(obj)
          })
          
setMethod(f="dimensions", signature="RFspatialDataFrame", 
	  definition=function(obj) {
            if (is(obj, "RFdataFrame")) return(1)
            else (getMethod("dimensions", "Spatial")@.Data)(obj)
          })
setMethod(f="dimensions", signature="RFdataFrame", 
	  definition=function(obj) return(1))


summary.RFsp <- function(object, ...) {
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
  class(obj) = "summary.RFsp"
  return(obj)
}

print.summary.RFsp <- function(x, ...){
  x <- summary.RFsp(x)
  str(x, give.attr=FALSE) #
  invisible(x)
}

print.RFsp <- function(x, ...) print.summary.RFsp(summary.RFsp(x, ...))#

## plot methods
## for 1-dim coordinates for grid and points
## for 2-or more dimensional coordinates for grid only

#for (.x in c("RFgridDataFrame", "RFpointsDataFrame"))
#  for (.y in c("missing", "RFgridDataFrame", "RFpointsDataFrame"))
#  setMethod(f="plot", signature(x=.x, y=.y),
#            definition=function(x, y, nmax = 6,
#              plot.variance = (!is.null(x@.RFparams$has.variance) &&
#                               x@.RFparams$has.variance),
#              ...)
#            plotRFdataFrame(x, y=y, nmax=nmax, plot.variance=plot.variance, ...)
#            )


setMethod(f="plot", signature(x="RFdataFrame", y="missing"),
          definition=function(x, y, nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),legend=TRUE, 
            ...)
          plotRFdataFrame(x, y=y, nmax=nmax, plot.variance=plot.variance,
                          legend=legend, ...)
          )

setMethod(f="plot", signature(x="RFdataFrame", y="RFdataFrame"),
          definition=function(x, y, nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),legend=TRUE, 
            ...)
          plotRFdataFrame(x, y=y, nmax=nmax, plot.variance=plot.variance,
                          legend=legend, ...)
          )

setMethod(f="plot", signature(x="RFdataFrame", y="matrix"),
          definition=function(x, y, nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),legend=TRUE, 
            ...)
          plotRFdataFrame(x, y=y, nmax=nmax, plot.variance=plot.variance,
                          legend=legend, ...)
          )


setMethod(f="plot", signature(x="RFdataFrame", y="data.frame"),
          definition=function(x, y, nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),legend=TRUE, 
            ...)
          plotRFdataFrame(x, y=y, nmax=nmax, plot.variance=plot.variance,
                          legend=legend, ...)
          )



setMethod(f="plot", signature(x="RFspatialDataFrame", y="missing"),
	  definition=function(
            x, y, MARGIN = c(1, 2),
            MARGIN.slices = NULL,
            n.slices= if (is.null(MARGIN.slices)) 1 else 10,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables, # = 1:vdim,
            zlim,
            legend=TRUE,  MARGIN.movie = NULL, ...)
          plotRFspatialDataFrame(x=x, y=y, MARGIN=MARGIN,
                                 MARGIN.slices=MARGIN.slices,
                                 n.slices=n.slices, nmax=nmax,
                                 plot.variance = plot.variance,
                                 select=select.variables,
                                 zlim=zlim,
                                 legend=legend,
                                 MARGIN.movie = MARGIN.movie,
                                 ...))

setMethod(f="plot",
          signature(x="RFspatialDataFrame", y="RFspatialGridDataFrame"),
	  definition=function(
            x, y, MARGIN = c(1, 2),
            MARGIN.slices = NULL,
            n.slices= if (is.null(MARGIN.slices)) 1 else 10,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables, # = 1:vdim,
            zlim,
            legend=TRUE,
            MARGIN.movie = NULL,
            ...)
          plotRFspatialDataFrame(x=x, y=y, MARGIN=MARGIN,
                                 MARGIN.slices=MARGIN.slices,
                                 n.slices=n.slices, nmax=nmax,
                                 plot.variance = plot.variance,
                                 select=select.variables,
                                 zlim=zlim,
                                 legend=legend,
                                 MARGIN.movie = MARGIN.movie,
                                 ...))

setMethod(f="plot",
          signature(x="RFspatialDataFrame", y="RFspatialPointsDataFrame"),
	  definition=function(
            x, y, MARGIN = c(1, 2),
            MARGIN.slices = NULL,
            n.slices= if (is.null(MARGIN.slices)) 1 else 10,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables, # = 1:vdim,
            zlim,
            legend=TRUE,
            MARGIN.movie = NULL,
            ...)
          plotRFspatialDataFrame(x=x, y=y, MARGIN=MARGIN,
                                 MARGIN.slices=MARGIN.slices,
                                 n.slices=n.slices, nmax=nmax,
                                 plot.variance = plot.variance,
                                 select=select.variables,
                                 zlim=zlim,
                                 legend=legend,
                                 MARGIN.movie = MARGIN.movie,
                                 ...))
setMethod(f="plot",
          signature(x="RFspatialDataFrame", y="matrix"),
	  definition=function(
            x, y, MARGIN = c(1, 2),
            MARGIN.slices = NULL,
            n.slices= if (is.null(MARGIN.slices)) 1 else 10,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables, # = 1:vdim,
            zlim,
            legend=TRUE,
            MARGIN.movie = NULL,
            ...)
          plotRFspatialDataFrame(x=x, y=y, MARGIN=MARGIN,
                                 MARGIN.slices=MARGIN.slices,
                                 n.slices=n.slices, nmax=nmax,
                                 plot.variance = plot.variance,
                                 select=select.variables,
                                 zlim=zlim,
                                 legend=legend,
                                 MARGIN.movie = MARGIN.movie,
                                 ...))

setMethod(f="plot",
          signature(x="RFspatialDataFrame", y="data.frame"),
	  definition=function(
            x, y, MARGIN = c(1, 2),
            MARGIN.slices = NULL,
            n.slices= if (is.null(MARGIN.slices)) 1 else 10,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables, # = 1:vdim,
            zlim,
            legend=TRUE,
            MARGIN.movie = NULL,
            ...)
          plotRFspatialDataFrame(x=x, y=y, MARGIN=MARGIN,
                                 MARGIN.slices=MARGIN.slices,
                                 n.slices=n.slices, nmax=nmax,
                                 plot.variance = plot.variance,
                                 select=select.variables,
                                 zlim=zlim,
                                 legend=legend,
                                 MARGIN.movie = MARGIN.movie,
                                 ...))

setMethod(f="persp",
          signature(x="RFspatialGridDataFrame"),
	  definition=function(
            x, y, MARGIN = c(1, 2),
            MARGIN.slices = NULL,
            n.slices= if (is.null(MARGIN.slices)) 1 else 10,
            nmax = 6,
            plot.variance = (!is.null(x@.RFparams$has.variance) &&
                             x@.RFparams$has.variance),
            select.variables, # = 1:vdim,
            zlim,
            legend=TRUE, MARGIN.movie = NULL, ...)
          plotRFspatialDataFrame(x=x, y=y, MARGIN=MARGIN,
                                 MARGIN.slices=MARGIN.slices,
                                 n.slices=n.slices, nmax=nmax,
                                 plot.variance = plot.variance,
                                 select=select.variables,
                                 zlim=zlim,
                                 legend=legend,
                                 MARGIN.movie = MARGIN.movie,
                                 ..., plotmethod="persp"))

contour.RFspatialGridDataFrame <-
  function( x, y, zlim,
           MARGIN = c(1, 2),
           MARGIN.slices = NULL,
           n.slices = if (is.null(MARGIN.slices)) 1 else 10,
           nmax = 6,
           plot.variance = (!is.null(x@.RFparams$has.variance) &&
                            x@.RFparams$has.variance),
           select.variables, # = 1:vdim,
            legend=TRUE,  MARGIN.movie = NULL, ...)
  plotRFspatialDataFrame(x=x, y=y, MARGIN=MARGIN,
                         MARGIN.slices=MARGIN.slices,
                         n.slices=n.slices, nmax=nmax,
                         plot.variance = plot.variance,
                         select=select.variables,
                         zlim=zlim,
                         legend=legend,MARGIN.movie = MARGIN.movie,
                         ..., plotmethod="contour")


errMsgNoPlotAvailable <- function(x, y)
  warning(paste("no plot method available for signature c(",
                class(x), ",", class(y), ")"))

setMethod(f="plot", signature(x="RFspatialPointsDataFrame",
            y="RFspatialGridDataFrame"),
	  definition=function(x, y) {
            errMsgNoPlotAvailable(x, y)
            return(invisible(NULL))
          })
 



#summary.RandomFieldsReturn <- function(object, ...) {
#  class(object) <- "summary.RandomFieldsReturn"
#  object
#}
#print.RandomFieldsReturn <- function(x, give.attr=FALSE, ...) {
#  print.summary.RandomFieldsReturn(summary.RandomFieldsReturn(x, ...))#
#}
#
#print.summary.RandomFieldsReturn <- function(x, give.attr=FALSE, ...) {
#  #str(x)
#  class(x) <- "double"
#  str(x, no.list=TRUE, give.attr=FALSE) #
#  if (!give.attr) {
#    cs <- attr(x, "coord_system")
#    if (!(cs[1] %in% c("auto", "cartesian")) ||
#        (cs[2] != "keep" && cs[1] != cs[2])) {
#      cat(' - attr(*, "coord_system")=')
#      str(cs, no.list=TRUE, give.attr=FALSE) #
#    }
#  }
#  # class(x) <- "RandomFieldsReturn"
#  invisible(x)
#}

