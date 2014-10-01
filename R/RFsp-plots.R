

## use only via calls of 'eval(body())' from 'plotRFspatialGridDataFrame' or
## 'plotRFspatialPointsDataFrame'

default.image.par <- function(data.range, var.range, legend=TRUE) {
  data.col <- colorspace::heat_hcl(12, c. =c(80, 30), l = c(30, 90), power = c(1/5, 1.5))
  var.col <- colorspace::rainbow_hcl(12, c = 50, l = 70)

  list(data=list(dot.name="col", default.col=data.col,
         pch=16, cex=1, range=data.range),
       var=list(dot.name="var.col", default.col=var.col,
          pch=16, cex=1, range=var.range),
       #mar = c(2,2,1,0),
       #mar.leg = c(0, 2.2, 2, 0.2),
       legend = legend,
       lower.leg = if (legend) 0.85 else 1,
       arrows = list(reduction = 1.5, nx.vectors = 20, leg.pos=c(1, 0.7)),
       text.col="blue"
     )
}

my.arrows <- function(xy, z, r, thinning, col, nrow) {
  half <- as.integer(thinning / 2)
  thinned <- c(rep(FALSE, half), TRUE, rep(FALSE, thinning - half))
  if (!missing(nrow) && !is.null(nrow)) {
    thinned <- as.vector(outer(rep(thinned, length = nrow),
                               rep(thinned, length = nrow(xy) / nrow), "&"))
  }

 # Print(thinned, z, xy, nrow, thinning); xxx
  
  arrows(x0=xy[thinned, 1] - r/2*z[thinned,1],
         y0=xy[thinned, 2] - r/2*z[thinned,2],
         x1=xy[thinned, 1] + r/2*z[thinned,1],
         y1=xy[thinned, 2] + r/2*z[thinned,2], length=0.03, col=col)
}

prepareplotRFsp <- function(x, vdim, select, plot.var,
                            data.range, var.range,
                            MARGIN, n,
                            n.slices, plot.legend, zlim,
                            ...) {

#  Print(vdim, select)
  
  if (vdim == 1 && !identical(select, vdim))
    stop("the given select.variables does not match the data")
   
  timespacedim <-
    if (is(x, "RFspatialGridDataFrame")) length(x@grid@cellsize)
    else ncol(x@coords)
  if (!(length(MARGIN)==2))
    stop("MARGIN must have length 2")
  if (!all(MARGIN %in% 1:timespacedim))
    stop("chosen MARGINS are out of bounds")

  if (!missing(zlim)){
    if (is.character(zlim)) stopifnot(zlim=="joint")
    mychk <- function(zlim) stopifnot((is.null(dim(zlim)) && length(zlim)==2) ||
                                      (is.matrix(zlim) && nrow(zlim)==2))
    if (is.numeric(zlim)) mychk(zlim)
    if (is.list(zlim)) {
      stopifnot(names(zlim) %in% c("data", "var"))
      lapply(zlim, mychk)
    }
  }
  
 
  coord.units=x@.RFparams$coord.units
  variab.units=x@.RFparams$variab.units;

  graphics <- RFoptions()$graphics
 
  image.par <- default.image.par(data.range, var.range, legend=plot.legend) 
  names.rep <- c(paste("realization", 1:(n-plot.var), sep=" "),
                 "kriging variance")
  names.vdim <- 
    if (!is.null(names(x@data)) && all(nchar(names(x@data))>0)) {
      if (is.list(select)) {
        u <- unlist(lapply(strsplit(names(x@data), ".n"),
                           FUN=function(li) li[[1]]))
        lapply(select, function(indices)
               if (length(indices)==3) {
                 if (FALSE) warning("first component interpreted as scalar, the other two as vector")
                  paste(u[indices[1]], paste(u[indices[-1]], collapse="/"),
                        sep=" and ")
               } else if (length(indices)>3) {
                 stop("currently, only two-dimensional vectors can be plotted")
               } else paste(u[indices], collapse="/"))
      } else unlist(lapply(strsplit(names(x@data)[unlist(select)], ".n"),
                           FUN=function(li) li[[1]]))
    } else {
      paste("variable", select)
    }
  
  names.coords <-
    if (.hasSlot(x, "coords"))
      dimnames(x@coords)[[2]]
    else
      names(x@grid@cellcentre.offset)
     
   dots <- list(...)
  names.graphics <- names(graphics)
  for (i in 1:length(graphics)) {
    dots[[names.graphics[i]]] <- NULL
  }  
  dotnames <- names(dots)

  if (bgInDots <- "bg" %in% dotnames) {
    bgdots <- dots$bg
    dots$bg <- NULL
  }

  lab <- xylabs(names.coords[MARGIN[1]], names.coords[MARGIN[2]],
                units=coord.units )
                 
  
  if (!("xlab" %in% dotnames)) dots$xlab <- lab$x
  if (!("ylab" %in% dotnames)) dots$ylab <- lab$y
  if (!("pch" %in% dotnames)) dots$pch=image.par$data$pch
  if (!("cex" %in% dotnames)) dots$cex=image.par$data$cex


  ## colours, legend
  for (i in c("data", "var")) {
    if (i=="var" && !plot.var) next
    if (is.null(colour <- dots[[ image.par[[i]]$dot.name ]]))
      colour <- image.par[[i]]$default.col
    image.par[[i]]$col <- if (is.list(colour)) colour else list(colour)
    dots[[ image.par[[i]]$dot.name ]] <- NULL
    lencol <- length(colour)
    image.par[[i]]$range <-
      apply(image.par[[i]]$range, 2, function(x) {
        if (is.logical(all.equal(x[1], x[2])))
          warning("range of data is a single value")
        stopifnot(all(is.finite(x)))         
        r <- range(outer(x, 1 + 0.05 * c(-1,1), "*"))
        if (identical(r[1], r[2])) r[2] <- r[2]+1
        return(r)
      } )

    if (!missing(zlim)) {
      idx <- unlist(lapply(as.list(select), 
                           FUN=function(x) if (length(x) != 2) x[1] else NA))
      idx <- idx[!is.na(idx)]
      if (is.character(zlim) && zlim=="joint") {
        image.par[[i]]$range[,idx] <- c(min(image.par[[i]]$range[,idx]),
                                        max(image.par[[i]]$range[,idx]))
      } else {
        zz <- if (is.list(zlim)) zlim[[i]] else zlim
        if (!is.null(zz))
          image.par[[i]]$range[,idx] <- matrix(zz, nrow=2, ncol=length(idx))
      }
    }

    image.par[[i]]$z.legend <-
      as.matrix(apply(image.par[[i]]$range, 2,
                      function(x) seq(x[1], x[2], length=lencol)))
    image.par[[i]]$breaks <- 
      as.matrix(apply(image.par[[i]]$range, 2,
                      function(x) seq(x[1]-abs((x[2] - x[1])*1e-3),
                                      x[2]+abs((x[2]-x[1])*1e-3), len=lencol+1)
                      ))
  }


  ### Splitting & Legend Plotting
  if (vdim>1) # note: one of n and n.slices is always equal to 1
    split.main <- c(n * n.slices, length(select))
  else {
    if (n * n.slices > 1)
      split.main <- c(ceiling(n * n.slices/2), 2) else {
        split.main <- c(1,1)
      }
  }
  
  ArrangeDevice(graphics, figs=split.main) ## NIE par() o.ae. vor ArrangeDevice !!!!
  
  if (length(dev.list()) > 0 &&
      any(par()$mfcol != c(1,1))) par(mfcol=c(1,1)) ## first figure appears
  par(cex=dots$cex) ## NIE par() o.ae. vor ArrangeDevice !!!!
  if (bgInDots) par(bgdots)


  
  xlab.given <-  1 - as.integer("xlab" %in% dotnames && (is.null(dots$xlab) ||
                                                         dots$xlab==""))
  #if (always.close <- any(split.main > 1) || graphics$always_close_screen){
  mar <- c(1,1,0,0)
  image.par$data$mar.leg <- mar #+ c(0,0,2,0)
  image.par$var$mar.leg <- mar + c(0,0,0,0)
  oma.top <- 2*plot.legend + if (is.null(dots$main)) 0 else 2# + xlab.given
  oma.left <- 2 * (1 + xlab.given)
  oma.bottom <- oma.left + 0*(plot.legend && plot.var)
  oma <- c(oma.bottom, oma.left, oma.top, xlab.given) + 0.2
  #} else {
  #  mar <- c(rep(2 * (1 + xlab.given), 2), 0, 0)
  #  image.par$mar.leg <- c(1, mar[2], 2, 0) 
  #  oma <- rep(0.2, 4) + (!is.null(dots$main)) * c(0,0,2,0)
  #}
  figs <- c(length(select), prod(split.main) / length(select))
  
  
#  scr <- split.screen(rbind(c(0,1,if (plot.var) 1-image.par$lower.leg else 0,
#                              image.par$lower.leg),
#                            if (plot.legend) c(0, 1, image.par$lower.leg, 1),
#                            if (plot.var && plot.legend)
#                            c(0,1,0, 1-image.par$lower.leg)
#                            ))

 scr <- split.screen(rbind(
    c(0,1, 0,
      if (plot.var) 1-2*(1-image.par$lower.leg) else image.par$lower.leg),
    if (plot.legend) c(0, 1, image.par$lower.leg, 1),
    if (plot.var && plot.legend) c(0,1,1-2*(1-image.par$lower.leg),
                                   image.par$lower.leg)
    ))

  scr.main <- matrix(nrow = figs[2],
                     split.screen(split.main, scr[1]),
                     byrow=TRUE) # statt nr=split.main[1]
  scr.legends <- integer(0)
  
  if (plot.legend) {
    SCR <- scr
    for (i in c("data", "var")) {
      if (i == "var" && !plot.var) next
      SCR <- SCR[-1]
      screen(SCR[1])

      scr.leg <- split.screen(figs=c(1, length(select)))
      scr.legends <- c(scr.legends, scr.leg)
      
     for (jx in 1:length(scr.leg)) {
        j <- if (is.list(select)) select[[jx]] else select[jx]
        screen(scr.leg[jx])
        par(oma=oma, mar=image.par[[i]]$mar.leg)
        if (!is.list(select) || length(select[[jx]]) != 2) {
          #Print(i,j,image.par)
          lab <- xylabs("", "", units=coord.units)
          image(x=image.par[[i]]$z.legend[, j[1]], y=c(0,1),
                z=matrix(ncol=2, rep(image.par[[i]]$z.legend[, j[1]], 2)),
                axes=FALSE,
                xlab=lab$x, ylab=lab$y,
                col=image.par[[i]]$col[[1 +(jx-1) %% length(image.par[[i]]$col)]]
              )
          axis(3, mgp=if (plot.var) c(3,0,0), hadj=if (plot.var) -0.5 else NA)#1 + 2 * (i =="data"))
          box()
        }
      }
    }
  }
  
  #close.screen(scr)

  return(c(image.par,
           dots=list(dots),
           list(names.coords = names.coords, names.rep = names.rep,
                names.vdim = names.vdim, 
                mar=mar, oma=oma, scr.main=scr.main, scr.legends=scr.legends,
                scr = scr,
               split.main=split.main,
               #always.close=always.close, 
                grPrintlevel = graphics$grPrintlevel)))
}





PlotTitle <- function(x, main) {
  p <- x@.RFparams
  if (!is.null(p$krige.method))  main <- paste(p$krige.method, "\n", main)
  par(mar=rep(0, 4), new=TRUE)
  plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE)
  text(0, 0.5, labels=main, adj=1, xpd=NA, col="blue", cex=0.8) # links davon
}




plotRFspatialDataFrame <-  
  function(x, y,
           MARGIN, #=c(1,2),      # which dimensions are to be plotted
           MARGIN.slices, #=NULL,
              # in which dimension a sequence of slices is to be plotted
           n.slices, #=if (!is.null(MARGIN.slices)) 6 else 1,
           nmax, #=6,             # max number of repetitins plotted
           plot.variance, # = !is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance,
           select, #1:vdim,
           zlim, # default: missing,
           legend, #=TRUE,
           MARGIN.movie,
          ..., plotmethod="image") {

  x.grid <- is(x, "RFspatialGridDataFrame")
  do.slices <- !is.null(MARGIN.slices)
  do.movie <- !is.null(MARGIN.movie)
  if (length(MARGIN.slices) > 1) stop("MARGIN.slices must be a scalar.")
  if (length(MARGIN.movie) > 1) stop("MARGIN.movie must be a scalar.")
  if (!x.grid) {
    if (is(x, "RFspatialPointsDataFrame")) {
      if (do.slices || n.slices[length(n.slices)] != 1)
        stop("'MARGIN.slices' must be 'NULL' and 'n.slices' must be 1.")
    } else {
      stop("method only for objects of class 'RFspatialPointsDataFrame' and 'RFspatialGridDataFrame'")
    }
  }
  NEWMARGIN <- MARGIN

  has.variance <-
    !is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance
  if (!has.variance) plot.variance <- FALSE
  
  if (x.grid) {
    conventional <- rfspDataFrame2conventional(x)
    x.grid.vectors <-
      GridTopology2gridVectors(cbind(conventional$x,conventional$T))
    
    ## array with  dims  (space-time-dims, vdim, n)  AND drop=FALSE!! 
    timespacedim <-  genuine.timespacedim <- length(x.grid.vectors)
    if (timespacedim!=length(x@grid@cellsize))
      stop("sollte nicht auftauchen: programming error in plotRFspatialGridDataFrame, timespacedim wrong ... (AM)")
    
    if (do.slices){
      if (!(MARGIN.slices <= timespacedim))
        stop("chosen MARGIN.slices out of bounds")
      if (MARGIN.slices %in% MARGIN)
        stop("MARGIN.slices must be different from MARGIN")
    }
    
    if (length(n.slices)!=1 && length(n.slices)!=3)
      stop("n.slices must be an integer of length 1 or 3")

    data.arr <- RFspDataFrame2dataArray(x)
    vdim <- dim(data.arr)[timespacedim+1]
    n.orig <- dim(data.arr)[timespacedim+2] ## including the kriging variance as one repet
    n.ohne.var <- n.orig - has.variance
    n <- min(n.ohne.var, nmax) + plot.variance

    ## want to have at least 3 space-time dims, if only 2,
    ## generate artificial dim
    dimdata <- dim(data.arr)  # new dims


    if (timespacedim <= 3){
      if (timespacedim == 3)
        dim(data.arr) <- c(dimdata[1:3], 1, dimdata[4:5])
      else if (timespacedim==2)
        dim(data.arr) <- c(dimdata[1:2], 1, 1, dimdata[3:4])
      else stop("dimension too small: dim=", timespacedim)
      timespacedim <- 4
    }
   if (!all(MARGIN %in% 1:(length(dimdata) - 2))) stop("MARGIN out of range.")
    
    if (any(MARGIN.slices %in% MARGIN.movie))
      stop("MARGIN.slices and MARGIN.movie are not disjoint.")
    if (length(MARGIN.slices) < 1)
      MARGIN.slices <-
        (1:(max(MARGIN, MARGIN.movie) + 1))[-c(MARGIN, MARGIN.movie)][1]
    if (length(MARGIN.movie) < 1)
      MARGIN.movie <-
        (1:(max(MARGIN, MARGIN.slices) + 1))[-c(MARGIN, MARGIN.slices)][1]
    
    dimdata <- dim(data.arr)
    vdimrep <- dimdata[(-1:0) + length(dimdata)]
    if (!all(c(MARGIN.movie, MARGIN.slices) %in%
             1:(length(dimdata) - 2))) stop("MARGINs out of range.")


    xx <- x.grid.vectors[[MARGIN[1]]]
    xy <- x.grid.vectors[[MARGIN[2]]]

    ## reduce the data array by taking a section with respect
    ## to the dimensions not covered by MARGIN*
    mar.vec <- c(MARGIN, MARGIN.slices, MARGIN.movie)
    ind <- as.list(rep(1, length(dimdata) - 2))
    ind[mar.vec] <- TRUE
    data.arr <- do.call("[", c(list(data.arr), ind, TRUE, TRUE, drop=FALSE))
    dim(data.arr) <- c(dimdata[sort(mar.vec)], vdimrep)

    ## re-oredered dimensions such that MARGIN and MARGIN.slices
    ## and MARGIN.movie are the first 4 dimensions
    perm.tmp <- c(mar.vec, (-1:0) + length(dimdata))
    data.arr <- aperm(data.arr, perm.tmp)
    dimdata <- dim(data.arr)
    NEWMARGIN <- 1:2
    MARGIN.slices <- 3
    MARGIN.movie <- 4

    if (do.slices) {
      if (n != 1) {
        n <- 1
        message("only first realization is shown")
      }
      if (plot.variance){
        plot.variance <- FALSE
        message("plot.variance was set to FALSE")
      }
      mar.len <- dimdata[MARGIN.slices]
      if (n.slices[length(n.slices)] > mar.len)
        n.slices[length(n.slices)] <- mar.len
      slices.ind <-
        if (length(n.slices) == 1) seq(1, mar.len, length=n.slices)
        else seq(n.slices[1], n.slices[2], length=n.slices[3])
      slices.ind <- unique(round(slices.ind))
      slices.ind <- slices.ind[slices.ind >= 1 & slices.ind <= mar.len]
      ## which indices in the slices dimension are taken
    } else {
      slices.ind <- 1
      ## the first dimension which is not in MARGIN
    }  
    n.slices <- length(slices.ind)
           
   ## ersten n.orig-1 sind wiederholungen, die letzte 'Spalte'
    ## ist die Varianz falls existent
    data.idx <- 1 : (n.ohne.var*vdim)    
    
    all.i <- as.matrix(expand.grid(1:n.slices, 1:n)[2:1]) ## i, ii
    coords <- as.matrix(expand.grid(xx, xy))
    m.range <- if (do.movie) 1:dimdata[MARGIN.movie] else 1

#    Print(m.range, do.movie, dimdata, MARGIN.movie); kk
    
  } else { ## not grid
    vdim <- x@.RFparams$vdim
    n <- min(x@.RFparams$n, nmax) + plot.variance
    nc <- ncol(x@data)
    if (nc < n*vdim)
    if (n==1) vdim <- nc else if (vdim==1) n <- nc else {
      stop("ncol(x@data) does not match 'x@.RFparams'; change 'x@.RFparams'")
    }
    data.idx <- 1:(x@.RFparams$n*vdim)
    all.i <- cbind(1:n, 1)
    genuine.timespacedim <- ncol(x@coords)
    coords <- x@coords[, MARGIN]
    m.range <- 1
  }

  if (!(missing.y <- missing(y))) {
    if (do.slices)
      stop("'y' and 'MARGIN.slices' may not be given at the same time")
    #
    if (is(y,  'RFspatialGridDataFrame')) {
      y.coords <- as(y, "RFspatialPointsDataFrame")@coords
      y.data <- y@data
    } else if (is(y, "matrix") || is(y, "data.frame")) {
      dc <- data.columns(y, xdim = dimensions(x), force=TRUE)
      y.coords <- y[, dc$x, drop=FALSE]
      y.data <- y[, dc$data, drop=FALSE]
    } else  {
      y.coords <- y@coords
      y.data <- y@data
    }     
  }

    

  data.range <-
    apply(as.matrix(1:vdim), 1,  # statt data.idx  # pro vdim eine legend
                    function(z) 
                    range(if (missing.y)
                          x@data[z + vdim * 0:(n-plot.variance-1)] else {
                            idx <- z + vdim * 0:(n-plot.variance-1)
                           # Print(idx, ncol(y.data), idx[idx <= ncol(y.data)])
                            c(#x@data[idx],
                              y.data[, idx[idx <= ncol(y.data)]])
                          },
                          na.rm=TRUE))

 # Print(data.range, vdim * 0:(n-plot.variance-1))
 
  var.range <- if (plot.variance) sapply(x@data[-data.idx],
                                         range, na.rm=TRUE) else NULL

  if (missing(select)) select <- 1:vdim
  image.par <- prepareplotRFsp(x=x, vdim=vdim, select=select,
                               data.range = data.range, var.range=var.range,
                               plot.var=plot.variance, MARGIN=NEWMARGIN,
                               n=n, n.slices=n.slices, plot.legend=legend,
                               zlim=zlim,
                              ...)

 image.par$names.vdim <-
    add.units(image.par$names.vdim, x@.RFparams$variab.units)
  

  if (x.grid) {
    nx.vectors <- min(length(xx), image.par$arrow$nx.vectors)
    thinning <- as.integer( (length(xx)-1) / nx.vectors)
  } else {
    nx.vectors <- min(nrow(coords), image.par$arrow$nx.vectors^2)
    thinning <- as.integer( (nrow(coords)-1) / nx.vectors^2)
  }

  ## split the left part according to 'split.main' for different vdims and
  ## repetitions

 # Print(dim(data.arr))
  
  for (m in m.range) {  
    for (jx in 1:length(select)) {       
      j <- if (is.list(select)) select[[jx]] else select[jx]
      for (ix in 1:nrow(all.i)) {

       # Print(m.range, length(select), all.i)
        
        i <- all.i[ix, ]
        dots <- dots.with.main.lab <- image.par$dots
        main <- dots$main
        dots$main <- NULL
        lab <- xylabs("", "", units=x@.RFparams$coord.units)
        dots$xlab <- lab$x
        dots$ylab <- lab$y
        if (do.plot.var <- (plot.variance && i[1]==n)){
          k <- if (x.grid) n.orig else x@.RFparams$n + 1;
          dv <- "var"
        } else {
          k <- i[1]
          dv <- "data"
        }
        
        screen(image.par$scr.main[ix, jx])
        par(mar=image.par$mar, oma=image.par$oma)
        
        col <-
          image.par[[dv]]$col[[1 + (j[1]-1) %% length(image.par[[dv]]$col) ]]
        breaks <- image.par[[dv]]$breaks[, j[1]]

        genuine.image <- (length(j) == 1 || length(j)==3)
        if (x.grid) {
          dots$type <- NULL
          dots$col <- if (genuine.image) col else par()$bg

          plot.return <- do.call(plotmethod,
                                 args=c(dots, list(
                                   x=xx, y=xy, z=data.arr[,,i[2], m, j[1], k],
                                   zlim = image.par[[dv]]$range[, j[1]],
                                   axes=plotmethod == "persp"))
                                 )
        } else {
          idx <- if (n==1) j else if (vdim==1) k else (k-1)*vdim+j
          dots$col <- if (genuine.image)
            col[ cut(x@data[,idx[1]], breaks=breaks) ] else par()$bg
          
          do.call(graphics::plot,
                  args=c(dots, list(x=coords[, 1], y=coords[, 2],
                    axes=FALSE)))
          box()
        }
        
        if (n.slices > 1)
          legend("bottomright", bty="n",
                 legend=paste(image.par$names.coords[MARGIN.slices], "=",
                   x.grid.vectors[[MARGIN.slices]][slices.ind[i[2]]]))
        
        if (do.plot.arrows <- length(j) >= 2 && !do.plot.var) {
          jj <-  if (length(j) == 3) j[-1] else jj <- j
          rx <- range(coords[, 1])
          ry <- range(coords[, 2])
          col.arrow <- if (length(image.par[["data"]]$col) >= jj[1] &&
                           length(image.par[["data"]]$col[[jj[1]]]) == 1)
            image.par[["data"]]$col[[jj[1]]] else "black"
          
          if (ix == 1) {
            factor <- image.par$arrow$reduction *
              sqrt(diff(rx) * diff(ry) / max(x@data[jj[1]]^2 +
                                             x@data[jj[2]]^2)) / nx.vectors
          }
          my.arrows(coords, x@data[jj], r = factor, thinning = thinning,
                    col = col.arrow, 
                    nrow = if (x.grid) length(xx))
        }
        
        if (!do.plot.var && !missing.y && (length(j) == 1 || (length(j)==3))) {
          idx <- if (n==1) j else if (vdim==1) i[1] else (i[1]-1)*vdim+j
          if (ncol(y.data) < idx) idx <- 1
          if (plotmethod == "persp") {
                                        # theta = 30, phi = 30, expand = 0.5, 
            xy <- trans3d(y.coords[, MARGIN[1]], y.coords[, MARGIN[2]],
                          data[ , idx], pmat=plot.return)
            points(xy, pch=16, col="black")          
          } else {
            col2 <- col[ cut(y.data[ , idx], breaks=breaks) ]
            dots2 <- dots
            dots2[c("type", "pch", "lty", "col", "bg", "cex", "lwd")] <- NULL
            addpoints <- function(pch, col, cex) {
              do.call(graphics::plot.xy,
                      args=c(dots2,
                        list(xy=xy.coords(y.coords[, MARGIN[1]],
                               y.coords[, MARGIN[2]]),
                             type="p", pch=pch, lty=1, col=col, bg=NA, cex=cex,
                             lwd=1)))
            }
            if (plotmethod=="image") addpoints(15, "darkgray", dots$cex*2)
            addpoints(dots$pch, col2, dots$cex)
          }
        }
        
        if (ix==1 ||
            ((image.par$split.main[1] != nrow(all.i)) &&
             (ix <= image.par$split.main[2]))) { # nrow(all.i) || ) #!image.par$always.close ||
          axis(1, outer=TRUE)#image.par$always.close)
        }
        if (jx==1 &&
            ((image.par$split.main[2] == length(select)) ||
             ((ix-1) %% image.par$split.main[2] == 0))) # !image.par$always.close || 
          axis(2, outer=TRUE)#image.par$always.close)
        ##if (!image.par$always.close){
        ##  dots2 <- dots.with.main.lab
        ##  dots2$xlab <- dots2$ylab <- ""
        ##  do.call(graphics::title, args=c(dots2, list(outer=TRUE, line=NA)))
        ##}
        ## if (!image.par$always.close) {
        ##   if (x.grid) do.call(graphics::title, args=c(dots, list(outer=FALSE, line=NA)))
        ##   else do.call(graphics::title, args=c(dots, line=2))  # line=-1 (ersetzt (AM))
        ##   ##if (i==1) axis(1, outer=image.par$always.close)
        ##   ##if (j==1) axis(2, outer=image.par$always.close)
        ## }
        
        if (all(i==1) && (image.par$grPrintlevel > 1 || vdim>1)) {
          mtext(text = image.par$names.vdim[jx],  # names(x)[j[1]]
                side=3, line=-1,
                col = image.par$text.col, cex=dots$cex)
          ##legend("topleft", bty="n", legend=c("", image.par$names.vdim[jx]),
          ##      text.col=image.par$text.col)
        }
        if (n>1 && jx==1){
          mtext(text = image.par$names.rep[ix], side=3, line=-2, cex=dots$cex)
          ##legend.pos <- "topright" #if (vdim==1) "topright" else "left"
          ##legend(legend.pos, bty="n", legend=c("", image.par$names.rep[ix]))
        }
        
        if (do.plot.arrows && ix == 1) {
          ## do not merge with ix==1 above !!
          ## reason: screen is changed here and going back
          ## is not a good idea
          
          if (image.par$legend) {
            screen(image.par$scr.legends[jx])
            do.call(graphics::plot, args=c(dots, list(x=Inf, y=Inf,
                                      xlim=rx, ylim=c(0,1), axes=FALSE)))
            len <- max(pretty(diff(rx) / image.par$arrows$nx.vectors/2 /factor))
            x.arrow<- cbind(mean(rx), image.par$arrows$leg.pos[1+genuine.image])
            my.arrows(x.arrow, cbind(len, 0), r = factor,
                      thinning=0, col=col.arrow)
            text(x.arrow, pos=1, labels = len)
          }
        }            
      }
    }
  
    
                                        #if (image.par$always.close) {
                                        #if (x.grid){
    dots.with.main.lab$type <- NULL # woher kommt dieses type??
    ##      falls nicht auf NULL gesetzt gibt es einen Fehler. Warum?
    do.call(graphics::title,
            args=list(main=dots.with.main.lab$main,
              outer=TRUE, line=image.par$oma[3]-1.5))
    dots.with.main.lab$main <- NULL
    do.call(graphics::title, args=c(dots.with.main.lab, list(outer=TRUE, line=NA)))
    ##} else {
    ##  do.call(graphics::title, args=c(dots, list(outer=TRUE)))      
    ##}
    

    if (image.par$grPrintlevel > 0) {
      screen(image.par$scr[1 + image.par$legend], new=FALSE)
      PlotTitle(x, if (is.null(main)) "" else main)
    }
  }

  
  
  graphics <- RFoptions()$graphics
  if (graphics$always_close_screen){
    close.screen(all.screens=TRUE)
  } else {
    close.screen(c(image.par$scr.main, image.par$scr.legends, image.par$scr))
  }

  return(invisible())
}



  
  





