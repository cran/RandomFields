regression <- function(x, y, main, scr,
                       value=function(regr) regr$coeff[2], variable="var",
                       mode=c("nographics", "plot", "interactive"),
                       regr.select=c("coefficients"),
                       averaging=FALSE, cex.main=0.85,
                       col.passive = "grey",
                       col.active = "green",
                       col.main.passive = "black",
                       col.main.active = "red",
                       col.points.chosen = "lightblue",
                       col.abline = "yellow",
                       col.abline.chosen = "darkblue",
                       col.smooth = "red",
                       ...) {
  bg <- par()$bg
  par(bg="white")
  mode <- match.arg(mode)
  if (is.array(y)) {
    stopifnot(length(x)==dim(y)[1])
    y <- as.vector(y)
  } else stopifnot(length(y) %% length(x) == 0)
  if (length(y)>length(x)) {
    if (averaging) y <- apply(matrix(y, nrow=length(x)), 1, mean)
    else x <- rep(x, len=length(y))
  }
  regr <- lm(y ~ x)[regr.select]
  val <- value(regr)
  names(val) <- NULL
  x.u <- y.u <- regr.u <- val.u <- NULL
  ## curmain must be set since within repeat loop curmain is not set if
  ## the user immediately breaks the loop
  curmain <- main <- paste(main,": ",variable,"=", format(val, dig=3), sep="")
  sm <- ksmooth(x, y, n.p=100, kernel="box", bandwidth=0.5)[c("x", "y")]
  if (mode!="nographics") { 
    screen(scr)
    par(mar=c(4.1, 4.1, 3, 0.1))
    plot(x, y,
         col=if (mode=="plot") col.passive else col.active, main=main,
         col.main=if (mode=="plot")  col.main.passive else col.main.active,
         cex.main=cex.main, ...)
    abline(regr, col=abline)
    lines(sm$x, sm$y, col=col.smooth)
    if (mode=="interactive") {
      repeat {
        if (length(loc <- Locator(2))==0) break;
        r <- range(loc$x)
        x.u <- x[idx <- x >= r[1] & x <= r[2]]
        y.u <- y[idx]
        if (diff(r)==0 || length(x.u) <= 1) next;
        regr.u <- lm(y.u ~ x.u)[regr.select]
        val.u <- value(regr.u)
        names(val.u) <- NULL
        curmain <- paste(main, ", ",variable,".u=", format(val.u, dig=3), sep="")
        screen(scr)
        plot(x, y, col=col.active, main=curmain, col.main=col.main.active,
             cex.main=cex.main, ...)
        points(x.u, y.u, col=col.points.chosen, ...)
        abline(regr, col=col.abline)
        lines(sm$x, sm$y, col=col.smooth)
        abline(regr.u, col=col.abline.chosen)
      }
      screen(scr)
      plot(x, y, col=col.passive, main=curmain, col.main=col.main.passive,
           cex.main=cex.main, ...)
      if (!is.null(x.u)) points(x.u, y.u, col=col.points.chosen, ...)
      abline(regr, col=col.abline)
      lines(sm$x, sm$y, col=col.smooth)
      if (!is.null(x.u)) abline(regr.u, col=col.abline.chosen)
    }
  }
  par(bg=bg)
  return(list(regr=regr, val=val,
         regr.u=regr.u, val.u=val.u, x.u=x.u, y.u=y.u,
         sm=sm))
}

  
hurst <-  function(x, y = NULL, z = NULL, data,
                   gridtriple = FALSE, sort=TRUE,
                   block.sequ = unique(round(exp(seq(log(3000), log(dim[1]),
                     len=min(100, dim[1]))))),
                   fft.m = c(1, min(1000, (fft.len - 1) / 10)),
                   fft.max.length = Inf, ## longer ts are cut down
                   method=c("dfa", "fft", "var"),
                   mode=c("plot", "interactive"),
                   pch=16, cex=0.2, cex.main=0.85,
                   PrintLevel=RFparameters()$Print,
                   height=3.5,
                   ...
                   ) {
  l.method <- eval(formals()$method)
  pch <- rep(pch, len=length(l.method))
  cex <- rep(cex, len=length(l.method))
  T <- NULL # if (!is.null(T)) stop("time not programmed yet")
  grid <- TRUE # if (!grid) stop("only grids are possible")
  
  method <- l.method[pmatch(method, l.method)]
  do.dfa <- any(method=="dfa")
  do.fft <- any(method=="fft")
  do.var <- any(method=="var")

  l.block.sequ <- l.dfa.var <- dfa <-
    l.varmeth.sequ <- l.varmeth.var <- varmeth <- 
      l.lambda <- l.I.lambda <- fft <- NULL
  
  modes <- c("nographics", "plot", "interactive")
  if (any(is.na(mode <- modes[pmatch(mode, modes)])))
    stop("unknown values of `mode'")
  
  if (missing(x)) {
    gridtriple <- TRUE
    ## stopifnot(grid) 
    x <- if (is.array(data)) rbind(1, dim(data), 1) else c(1, length(data), 1)
  }

  ct <- CheckXT(x=x, y=y, z=z, T=T, grid=grid, gridtriple=gridtriple)
  dim <- apply(cbind(ct$x, ct$T), 2,
               function(x) length(seq(x[1], x[2], x[3])))

  if (PrintLevel>2) cat("(formatting) ")
  if (ncol(ct$x)>1) {
    if (!is.array(data)) data <- array(as.double(data), dim=dim)
    if (sort) { ## so that the first dimension gives the side with
      ##           the greatest length (in points)
      ord <- order(dim, decreasing=TRUE)
      if (any(diff(ord)<0)) {
        data <- aperm(data, ord)
        ct$x <- ct$x[, ord, drop=FALSE]
      }
    }
    gc()
  } else {
    data <- as.matrix(data)
    dim <- c(dim, 1)
  }
  repet <- prod(dim[-1])

  ## periodogram, code taken from jean-francois coeurjolly
  if (do.fft) {
    if (PrintLevel>2) cat("periodogramm")
    fft.len <- min(dim[1], fft.max.length)
    fft.m <- round(fft.m)
    stopifnot(diff(fft.m)>0, all(fft.m>0), all(fft.m<=fft.len))
    l.I.lambda <- double(repet * (fft.m[2] - fft.m[1] + 1));
    .C("periodogram",
       data,
       as.integer(dim[1]), 
       as.integer(repet),## Produkt der anderen Dimensionen
       as.integer(fft.m),## Ausschnitt aus Fourier-Trafo aus Stueck nachf. Laenge
       as.integer(fft.len),## Reihe zerhackt in Stuecke dieser Laenge 
       as.integer(fft.len / 2), ## WOSO(?)-Sch\"aetzer
       l.I.lambda,
       PACKAGE="RandomFields", DUP=FALSE)
    l.lambda <-  log((2 * pi * (fft.m[1]:fft.m[2])) / fft.len)
  } 


  ## detrended fluctuation analysis
  if (do.dfa || do.var) {
    if (PrintLevel>2) {
      if (do.dfa) cat("detrended fluctuation; ")
      if (do.var) cat("aggregated variation; ")
    }
    stopifnot(all(diff(block.sequ)>0))
    l.block.sequ <- log(block.sequ)
    dfa.len <- length(block.sequ)
    l.dfa.var <- double(dfa.len * repet)
    l.varmeth.var <- double(dfa.len * repet)
    ## wird data zerstoert ?!
    ## log already returned!
    .C("detrendedfluc", as.double(data), as.integer(dim[1]), as.integer(repet),
       as.integer(block.sequ), as.integer(dfa.len),
       l.dfa.var, l.varmeth.var, PACKAGE="RandomFields", DUP=FALSE)
    ## 1:dfa.len since data could be a matrix; l.block.sequ has length dfa.len
    ## and is then shorter than l.varmeth.var!
    varmeth.idx <- is.finite(l.varmeth.var[1:dfa.len])
    l.varmeth.var <- l.varmeth.var[varmeth.idx]
    l.varmeth.sequ <- l.block.sequ[varmeth.idx]
  }

  gc() 
  if (PrintLevel>2) cat("\n")

  if (any(mode=="plot" | mode=="interactive"))
  {
    plots <- do.dfa + do.fft + do.var
    get(getOption("device"))(height=height, width=height * plots)
    par(bg="white")
    screens <- seq(0, 1, len=plots+1)
    screens <- split.screen(figs=cbind(screens[-plots-1], screens[-1], 0, 1))
    on.exit(close.screen(screens))
  }

  ## calculation and plot are separated into blocks since calculation time
  ## is sometimes huge. So here, the user might wait for any result quite
  ## a long time, but is than bothered only once.
  for (m in mode) {
    scr <- 0
    if (do.dfa) {
      scr <- scr + 1
      dfa <- regression(l.block.sequ, l.dfa.var, variable="H", pch=pch[1],
                        main="detrended fluctuation analysis", scr=scr,
                        cex=cex[1], value=function(regr) regr$coeff[2]/2, mode=m,
                        cex.main=cex.main, ...)
    }
    if (do.fft) {
      scr <- scr + 1
      fft <- regression(l.lambda, l.I.lambda, main="periodogram", pch=pch[2],
                        scr=scr,  value=function(regr) 0.5 * (1-regr$coeff[2]),
                        averaging = repet > 1,
                        mode=m, cex=cex[2], variable="H",
                        cex.main=cex.main, ...)
    }
    if (do.var) {
      scr <- scr + 1
      varmeth <- regression(l.varmeth.sequ, l.varmeth.var, variable="H",
                            pch=pch[3],
                            main="aggregated variation", scr=scr, cex=cex[3],
                            value=function(regr) 1 + 0.5 * (regr$coeff[2]-2),
                            mode=m, cex.main=cex.main,...)    
    }
  }
  
  if (PrintLevel>1) {
    cat("#################### Hurst #################### \n")
    cat(c(dfa.H=dfa$val, varmeth.H=varmeth$val, fft.H=fft$val))
    cat("---- by interactively defined regression interval:\n")
    cat(c(dfa.H=dfa$val.u, varmeth.H=varmeth$val.u, fft.H=fft$val.u))
    #cat("---- beta (Cauchy model) ----\n")
    #print(2 * (1 - c(dfa.H=dfa$val, fft.H=fft$val, varmeth.H=varmeth$val)))
    #print(2 * (1 - c(dfa.H=dfa$val.u, fft.H=fft$val.u,varmeth.H=varmeth$val.u)))
    cat("############################################### \n")
  }
  if (any(mode=="plot" | mode=="interactive")) close.screen(screens)
   
  return(list(dfa=list(x=l.block.sequ, y=l.dfa.var, regr=dfa$regr,
                sm=dfa$sm,
                x.u=dfa$x.u, y.u=dfa$y.u,  regr.u=dfa$regr.u,
                H=dfa$val, H.u=dfa$val.u),
              varmeth=list(x=l.varmeth.sequ, y=l.varmeth.var, regr=varmeth$regr,
                sm=varmeth$sm,
                x.u=varmeth$x.u, y.u=varmeth$y.u, regr.u=varmeth$regr.u,
                H=varmeth$val, H.u=varmeth$val.u),
              fft=list(x=l.lambda, y=l.I.lambda, regr=fft$regr,
                sm=fft$sm,
                x.u=fft$x.u, y.u=fft$y.u, regr.u=fft$regr.u,
                H = fft$val, H.u = fft$val.u)
              )
         )
}


fractal.dim <-
  function(x, y = NULL, z = NULL, data,
           grid=TRUE, gridtriple = FALSE,
           bin= seq(min(ct$x[3, ]) / 2, 
             min((ct$x[2,]-ct$x[1,]) / 4, vario.n * min(ct$x[3,]) + 1),
             min(ct$x[3,])),
           vario.n=5,
           sort=TRUE,
           #box.sequ=unique(round(exp(seq(log(1),
           #  log(min(dim - 1, 50)), len=100)))),
           #box.enlarge.y=1,
           #range.sequ=unique(round(exp(seq(log(1),
           #  log(min(dim - 1, 50)), len=100)))),
           fft.m = c(65, 86), ## in % of range of l.lambda
           fft.max.length=Inf,
           fft.max.regr=150000,
           fft.shift = 50, # in %; 50:WOSA; 100: no overlapping
           method=c("variogram", "fft"),# "box","range", not correctly implement.
           mode=c("plot", "interactive"),
           pch=16, cex=0.2, cex.main=0.85,
           PrintLevel = RFparameters()$Print,
           height=3.5,
           ...) {
  l.method <- eval(formals()$method)
  pch <- rep(pch, len=length(l.method))
  cex <- rep(cex, len=length(l.method))
  T <- NULL # if (!is.null(T)) stop("time not programmed yet")

  method <- l.method[pmatch(method, l.method)]
  do.vario <- any(method=="variogram")
  do.box <- any(method=="box")     # physicists box counting method
  do.range <- any(method=="range") # Dubuc et al.
  do.fft <- any(method=="fft")

  ev <- l.midbins <- l.binvario <- vario <-
    Ml.box.sequ <- l.box.count <- box <-
      Ml.range.sequ <- l.range.count <- rnge <-
        l.lambda <- l.I.lambda <- fft <- NULL

  modes <- c("nographics", "plot", "interactive")
  if (any(is.na(mode <- modes[pmatch(mode, modes)])))
    stop("unknown values of `mode'")
  
  if (missing(x)) {
    gridtriple <- TRUE
    stopifnot(grid)
    x <- if (is.array(data)) rbind(1, dim(data), 1) else c(1, length(data), 1)
  }
  
  ct <- CheckXT(x=x, y=y, z=z, T=T, grid=grid, gridtriple=gridtriple)
  
  ## variogram method (grid or arbitray locations)
  if (do.vario) {
    if (PrintLevel>2) cat("variogram; ")
    ev <- EmpiricalVariogram(x=x, y=y, z=z, T=T, data=data, grid=grid,
                             bin=bin, gridtriple=gridtriple)
    idx <-  is.finite(l.binvario <- log(ev$emp))
    l.midbins <- log((ev$centers[idx])[1:vario.n])
    l.binvario <- (l.binvario[idx])[1:vario.n]
  }

  if (grid) {
    dim <- apply(cbind(ct$x, ct$T), 2,
                 function(x) length(seq(x[1], x[2], x[3])))
    #box.sequ  ## bind box.sequ NOW -- dim will be changed later on!

    if (PrintLevel>2) cat("(formatting) ")
    if (ncol(ct$x)>1) {
      if (!is.array(data)) data <- array(as.double(data), dim=dim)
      if (sort) { ## see fractal
        ord <- order(dim, decreasing=TRUE)
        if (any(diff(ord)<0)) {
          data <- aperm(data, ord)
          ct$x <- ct$x[, ord, drop=FALSE]
        }
      }
    } else {
      data <- as.matrix(data)
      dim <- c(dim, 1)
    }
    repet <- prod(dim[-1])
    gc()
    
    if (do.range) {
      if (PrintLevel>2) cat("ranges")
      lrs <- length(range.sequ) ## range.sequ is bound right here!
      l.range.count <- double(lrs * repet)
      ## logarithm is already taken within minmax
      .C("minmax", as.double(data), as.integer(dim[1]),
         as.integer(repet), as.integer(range.sequ), as.integer(lrs),
         l.range.count, PACKAGE="RandomFields", DUP=FALSE, NAOK=TRUE)
      box.length.correction <- 0 ## might be set differently
      ##                            for testing or development
      Ml.range.sequ <- -log(range.sequ + box.length.correction)
    }
    
    if (do.box) {
      stop("this point cannot be reached")
      if (PrintLevel>2) cat("box counting; ")
      x <- ct$x[,1] / ct$x[3,1] ## alles auf integer
      factor <- diff(x[c(1,2)]) / diff(range(data)) * box.enlarge.y
      ## note: data changes its shape! -- should be irrelevant to any procedure!
      data <- matrix(data, nrow=dim[1])
      l.box.count <- double(length(box.sequ) * repet)
      .C("boxcounting",
         as.double(rbind(data[1,], data, data[nrow(data),])),
         as.integer(dim[1]),
         as.integer(repet), as.double(factor),
         as.integer(box.sequ), as.integer(length(box.sequ)),
         l.box.count, PACKAGE="RandomFields", DUP=FALSE)
      gc()
      
      box.length.correction <- 0 ## might be set differently
      ##                            for testing or development
      Ml.box.sequ <- -log(box.sequ + box.length.correction)
      l.box.count <- log(l.box.count)
    }

    if (do.fft) {
      ## periodogram, code taken from jean-francois coeurjolly and WOSA
      if (PrintLevel>2) cat("periodogramm; ")
      fft.len <- min(dim[1],  fft.max.length) ## the virtual length
     
      ## fft.m is given in percent [in log-scale]; it is now rewriten
      ## in absolute indices for the fft vector
      stopifnot(all(fft.m>=0), all(fft.m<=100))
      ## 0.5 * fft.len not fft.len since only first half should be considered
      ## in any case! (Since the fourier transform returns a symmetric vector.)
      spann <- log(2 * pi / c(0.5 * fft.len, 1))
      fft.m <- round(2 * pi / exp(spann[1] + diff(spann) * (1 - fft.m / 100)))
     
      l.lambda <- log((2 * pi * (fft.m[1]:fft.m[2])) / fft.len)
      l.I.lambda <- double(repet * (fft.m[2] - fft.m[1] + 1));
      .C("periodogram",
         data,
         as.integer(dim[1]),
         as.integer(repet),## Produkt der anderen Dimensionen
         as.integer(fft.m),## Ausschnitt aus Fourier-Tr aus Stueck nachf. Laenge
         as.integer(fft.len),## Reihe zerhackt in Stuecke dieser Laenge 
         as.integer(fft.len / 100 * fft.shift), ## shift (WOSA-Sch\"aetzer)
         l.I.lambda,
         PACKAGE="RandomFields", DUP=FALSE)
    }
  } else {
    if (do.box || do.range || do.fft) {
      if (PrintLevel>0) {
        cat("\nThe following methods are available only for grids:\n")
        if (do.box)   cat("  * box counting\n")
        if (do.range) cat("  * max-min method\n")
        if (do.fft)   cat("  * periodogram/WOSA\n")
      }
      do.box <- do.range <- do.fft <- FALSE
    }
  }
  if (PrintLevel>2) cat("\n")


  if (any(mode=="plot" | mode=="interactive")) {
    plots <- do.vario + do.box + do.range + do.fft
    get(getOption("device"))(height=height, width = height * min(3.4, plots))
    par(bg="white")
    screens <- seq(0, 1, len=plots+1)
    screens <- split.screen(figs=cbind(screens[-plots-1], screens[-1], 0, 1))
    on.exit(close.screen(screens))
  }
  
  for (m in mode) {
    scr <- 0
    if (do.vario) {
      scr <- scr + 1
      vario <-
        regression(l.midbins, l.binvario, variable="D", pch=pch[1],
                   main="variogram method", scr=scr, cex=cex[1],
                   value =function(regr) ncol(ct$x) + 1 - regr$coeff[2]/2,
                   mode=m, cex.main=cex.main, ...)
    }
    if (do.box) {
      scr <- scr + 1
      box <- regression(Ml.box.sequ, l.box.count,
                        variable="D", main="box counting", scr=scr, pch=pch[2],
                        value =function(regr) ncol(ct$x) - 1 + regr$coeff[2],
                        mode=m, cex=cex[2], cex.main=cex.main, ...)
    }
    if (do.range) {
      scr <- scr + 1
      rnge <- regression(Ml.range.sequ, l.range.count, variable="D",
                         pch=pch[3],
                         main="max-min method", scr=scr, mode=m, cex=cex[3],
                         value = function(regr) ncol(ct$x) - 1 + regr$coeff[2],
                         cex.main=cex.main,...)
    }
    if (do.fft) {
      scr <- scr + 1
      if (length(l.lambda)>fft.max.regr && PrintLevel>0)
        cat("too large data set; no fft regression fit")
      else
        fft <- regression(l.lambda, l.I.lambda, main="periodogram", pch=pch[4],
                          scr=scr,
                          value=function(regr) 2.5 + 0.5 * regr$coeff[2],
                          averaging=fft.len!=dim[1],
                          mode=m, variable="D", cex=cex[4], cex.main=cex.main,
                          ...)
    }
  }

  
  if (PrintLevel>1) {
    cat("##################### fractal D ############### \n")
    cat(c(D.vario=vario$val,# D.box=box$val, D.range=rnge$val,
          D.fft=fft$val))
    cat("---- by interactively defined regression interval:\n")
    cat(c(D.vario=vario$val.u, # D.box=box$val.u, D.range=rnge$val.u,
            D.fft=fft$val.u))
    #cat("----alpha (Cauchy)----\n")
    #DIM <- 1
    #print(2 * (DIM + 1 - c(D.vario=vario$val, D.box=box$val, D.range=rnge$val)))
    #print(2 * (DIM + 1 - c(D.vario=vario$val.u, D.box=box$val.u,
    #                       D.range=rnge$val.u))) 
    cat("############################################### \n")
  }
  if (any(mode=="plot" | mode=="interactive")) close.screen(screens)
 
  return(list(vario=list(ev=ev, vario.n=vario.n,
                x=l.midbins, y=l.binvario, regr=vario$regr,
                sm=vario$sm,
                x.u=vario$x.u, y.u=vario$y.u, regr.u=vario$regr.u,
                D = vario$val, D.u = vario$val.u),
#              box = list(x=Ml.box.sequ, y=l.box.count, regr=box$regr,
#                sm=box$sm,
#                x.u=box$x.u, y.u=box$y.u, regr.u=box$regr.u,
#                D = box$val, D.u = box$val.u),
#              range=list(x=Ml.range.sequ, y=l.range.count, regr=rnge$regr,
#                sm=rnge$sm,
#                x.u=rnge$x.u, y.u=rnge$y.u, regr.u=rnge$regr.u,
#                D = rnge$val, D.u = rnge$val.u),
              fft=list(x=l.lambda, y=l.I.lambda, regr=fft$regr,
                sm=fft$sm,
                x.u=fft$x.u, y.u=fft$y.u, regr.u=fft$regr.u,
                D = fft$val, D.u = fft$val.u)
              )
         )
}








