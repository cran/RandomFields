


.RandomFields.env <- new.env()

sleep.milli <- function(milli) {
  .C("sleepMilli", as.integer(milli))
  invisible(NULL)
}

Dev <- function(on, dev, ps=NULL, cur.cex=TRUE, paper="special",
                width=5, height=5, quiet=FALSE, innerwidth, innerheight,
                mai, horizontal=FALSE, ...){
  
   if (!missing(innerwidth) || !missing(innerheight)) {
    stopifnot(!missing(innerwidth), !missing(innerheight), !missing(mai))
    height <- innerheight + sum(mai[c(1,3)])
    width <- innerwidth + sum(mai[c(2,4)])
  }
  
  ## function to handle output device:
  ##   on: T=output device is activated; F=device will be closed
  ##   dev: device number or TRUE= postscript or FALSE=pdf
  ##        or character (then the name is interpreted as function name)
  ##   ps : postscript file name; only needed when dev is logical
  ##
  ##  Dev takes over the par.options of previous plots
  if (on) {
    par.orig <- if (!is.null(dev.list())) par(no.readonly=TRUE) else NULL
    if (exists(".dev.orig", envir=.RandomFields.env)) {
      warning("Dev has been still open (.dev.orig exists). Closing.")
      if (!is.null(try(Dev(FALSE, get(".dev.orig",
                                      envir=.RandomFields.env)$dev))))  
        rm(".dev.orig", envir=.RandomFields.env)
    }
    if ((cur.cex <- cur.cex && !is.null(par.orig)) && !is.null(dev.list())) {
      par.orig <- par(no.readonly=TRUE)
    }

     
    par.orig$new <- FALSE
    par.orig$fin <- par.orig$mai <- par.orig$pin <-
      par.orig$plt <- NULL
    
    devPrev <- dev.cur()
    if (is.logical(dev) || is.character(dev)) {
      keep <- FALSE
      if (is.null(ps)) stop("no name for the postscript file is given")
      else {
        if (is.logical(dev)) {
          ext <- "eps"
          exts <- c("ps", "eps", "pdf")
          if (any(l <- (splt <- rev(strsplit(ps,"\\.")[[1]]))[1] == exts) &&
              length(splt)>1){
            dev <- !l[3]
            ps  <- paste(rev(splt[-1]), collapse=".")
            ext <- exts[l][1]
          }
          if (dev) {
            fn <- paste(ps, ext, sep=".")
            if (!file.create(fn)) stop("The file ", fn," cannot be created")
            postscript(fn, paper=paper, width=width, height=height,
                       horizontal = horizontal, ...)
          } else {
            fn <- paste(ps,".pdf",sep="")
            pdf(fn, width=width, height=height, ...)
          }
        } else { # character
          fn <-  paste(ps,".", dev, sep="")
          if (!file.create(fn)) stop("The file ", fn, " cannot be created")
         txt <- paste(dev,"('", fn, "',width=width,height=height,...)", sep="")
          eval(parse(text=txt))
        }
        if (!quiet) cat("creating", fn, "\n")
       }
      if (!missing(mai)) par(mai=mai)
    } else {
      if (dev %in% dev.list()) {
        dev.set(dev)
      } else {
        stopifnot(is.finite(height+width))
        ScreenDevice(height=height, width=width)
      }
      keep <- dev < 3
    }

    if (cur.cex) par(par.orig) # uncommented 12.8.04 + nach unten
    if (exists("abc")) return()
    if (cur.cex && FALSE) { ## komisches Verhalten !! wenn die beiden Befehle
      ##              zusammengefasst werden (gekippte eps in Latex)
      ##              29.5.05
#      par(par.orig[39]) # $mfg
#      par(par.orig[-39])# uncommented 12. 8.04 + nach unten
    }
    assign(".dev.orig", 
           list(dev.prev=devPrev, dev.cur=dev.cur(), keep=keep),
           envir=.RandomFields.env)
  } else { # off
    if (!exists(".dev.orig", envir=.RandomFields.env)) stop("Dev is not open")
    if (dev.cur() != get(".dev.orig", envir=.RandomFields.env)$dev.cur) {
      warning("Dev is not the currently active device")
      dev.set(get(".dev.orig", envir=.RandomFields.env)$dev.cur)
    }
    if (length(dev.list())>0)
      if (get(".dev.orig", envir=.RandomFields.env)$keep) par(new=FALSE)
      else dev.off() ## obsolete
    if ((devPrev <- get(".dev.orig", envir=.RandomFields.env)$dev.prev) != 1)
      dev.set(devPrev)   
    rm(".dev.orig", envir=.RandomFields.env) 
  }
  invisible(NULL)
}
	
hostname<-function(){.C("hostname", h=paste(seq(0,0,l=100), collapse=""),
                        as.integer(100), PACKAGE="RandomFields")$h}

pid <- function() {.C("pid", i=integer(1), PACKAGE="RandomFields")$i}


FileExists <- function(file, printlevel=RFoptions()$general$printlevel) {
    ## for parallel simulation studies: the same data output file should not
  ## be created twice. So:
  ## 1. if file exists then assume another process has done the work already
  ## 2. if file.lock existss then assume another process is doing the work
  ## 3.a. otherwise create file.lock to show other processes that the process
  ##      will do the work
  ## 3.b. check if another process has started with the same work at the same
  ##      time it may happen that in case of simulatenous creation of file.lock
  ##      no process will do the work...(then the lock file will rest.)
  lock.ext <- ".lock";
  if (file.exists(file)) { #1.
    if (printlevel>=PL_ERRORS ) cat("'", file, "' already exists.\n");
    return(1)
  } else { 
    LockFile <- paste(file, lock.ext, sep="")
    if (file.exists(LockFile)) { #2.
      if (printlevel>=PL_ERRORS ) cat("'",file,"' is locked.\n");
      return(2);
    }
    PID <- pid();
    write(file=LockFile,c(PID,hostname()),ncolumns=2,append=TRUE); #3.a.
    Pid <- matrix(scan(LockFile,what=character(0), quiet=TRUE),nrow=2)
    if ((sum(Pid[1,]==PID)!=1) || (sum(Pid[1,]>PID)>0)){ #3.b.
      if (printlevel>PL_ERRORS )
        cat("Lock file of '", file, "' is knocked out.\n");
      return(3);
    }
  }
  return(0);
}

LockRemove <- function(file) {
  ## removes auxiliary files created by FileExists
  lock.ext <- ".lock";
  file.remove(paste(file, lock.ext, sep=""))
}


plotWithCircles <- function(data, factor=1.0,
                            xlim=range(data[,1])+c(-maxr,maxr),
                            ylim=range(data[,2])+c(-maxr,maxr),
                            col=1, fill=0, ...) {
  ## marked point process: presents positive values of data as radii of circles
  CIRLLE.X <- cos(seq(0,2*pi,l=20))
  CIRLLE.Y <- sin(seq(0,2*pi,l=20))
  circle <- function(x,r) { polygon(x[1]+ r* CIRLLE.X,x[2]+ r* CIRLLE.Y,
                                    col=fill, border=col) }
  ##r <- x$NormedData - min(x$NormedData) +1
  ##r <- r/max(r)/nrow(x$coord) * diff(xlim) * diff(ylim) * 2.5;
  maxr <- max(data[,3])
  plot(Inf, Inf, xlim=xlim, ylim=ylim, xlab="", ylab="",...)
  for (i in 1:nrow(data)) { circle(data[i,c(1,2)], factor*data[i,3]) }
}


Print <- function(..., digits=6, empty.lines=2) { #
  max.elements <- 99
  l <- list(...)
  n <- as.character(match.call())[-1]
  cat(paste(rep("\n", empty.lines), collapse="")) #
  for (i in 1:length(l)) {
    cat(n[i]) #
    if (!is.list(l[[i]]) && is.vector(l[[i]])) {
      L <- length(l[[i]])
      if (L==0) cat(" = <zero>")#
      else {
        cat(" [", L, "] = ", sep="")
        cat(if (is.numeric(l[[i]]))
            round(l[[i]][1:min(L , max.elements)], digits=digits)#
            else l[[i]][1:min(L , max.elements)]) #
        if (max.elements < L) cat(" ...")
      }
    } else {
       if (is.list(l[[i]])) {
        cat(" =  ") #
        str(l[[i]], digits.d=digits) #
      } else {
        cat(" =")
        if (length(l[[i]]) <= 100 && FALSE) {
          print(if (is.numeric(l[[i]])) round(l[[i]], digits=digits)#
                else l[[i]])
        } else {
          if (length(l[[i]]) > 1 && !is.vector(l[[i]]) && !is.matrix(l[[i]])
              && !is.array(l[[i]])) cat("\n")
          str(l[[i]]) #
        }
      }
    }
    cat("\n")
  }
}


vectordist <- function(x, diag=FALSE) {
  storage.mode(x) <- "double"
  res <- .Call("vectordist", t(x), diag)
  dimnames(res) <- list(dimnames(x)[[2]], NULL)
  return(t(res));
}

xylabs <- function(x, y, T=NULL, units=NULL) {
  if (is.null(units)) units <- RFoptions()$coords$coordunits
  xlab <- if (is.null(x)) NULL
          else if (units[1]=="") x else paste(x, " [", units[1], "]", sep="")
  ylab <- if (is.null(y)) NULL
          else if (units[2]=="") y else paste(y, " [", units[2], "]", sep="")
  Tlab <- if (is.null(T)) NULL
          else if (units[3]=="") T else paste(T, " [", units[3], "]", sep="")
  return (list(xlab=xlab, ylab=ylab, Tlab=Tlab))
}

add.units <- function(x,  units=NULL) {
    if (is.null(x)) return(NULL)
  if (is.null(units)) units <- RFoptions()$coords$varunits
  return(ifelse(units=="", x, paste(x, " [", units, "]", sep="")))
}


ScreenDevice <- function(height, width) {
  if (height > 0  && width > 0) {
    GD <- getOption("device")

    ispdf <- is.logical(all.equal(GD, "pdf"))
    isjpg <- is.logical(all.equal(GD, "jpeg"))
    if (ispdf) {
      GD <- pdf
    } else if (isjpg) {
      GD <- jpeg
    } else {
      ispdf <- is.function(GD) && is.logical(all.equal(args(pdf), args(GD)))
      isjpg <- is.function(GD) && is.logical(all.equal(args(jpeg), args(GD)))
    }
   
#   Print(args(pdf), args(GD), all.equal(args(pdf), args(GD)),
  #        all.equal(args(jpeg), args(GD)),
  #        is.function(GD), RFoptions()$graphics)
    
    if (ispdf || isjpg) {
      
      ##     Print(ispdf, isjpg)
      ##     Print(ispdf, RFoptions()$graphics        )
      
      graphics <- RFoptions()$graphics        
      if (graphics$filenumber == 999)
        stop("number of pictures exceeds max. number of pictures")        
      file <- graphics$file
      ext <- c(".pdf", ".jpg")[isjpg + 1]
      
       if (file != "") {
        if (!graphics$onefile) {
         file <- paste(file,
                        formatC(format="d", flag="0", graphics$filenumber,
                                width=3),
                        ext, sep="")
         RFoptions(graphics.filenumber = graphics$filenumber + 1)
        } else file <- paste(file, ext, sep="")
        ##            Print(file, dev.list())
         if (ispdf) GD(height=height, width=width, onefile=FALSE, file=file)
        else if (isjpg) GD(height=height, width=width,
                           file=file, units="in", res=graphics$resolution)
        else stop("device not recognized by 'RandomFields'")
       return()
      } else {
        if (isjpg) stop("'file' in 'RFoptions' not set.")
      }
    }
     
    args <- names(as.list(args(GD)))
    if (all(c("width", "height") %in% args) &&
        ( !any(c("file", "filename") %in% args)) || ispdf) {
      GD(height=height, width=width)        
      ##      Print("OK", height, width, RFoptions()$graphics, par()$cex)
     return()
    }
    
   if (RFoptions()$internal$warn_aspect_ratio) {
      RFoptions(warn_aspect_ratio = FALSE)
      cat("The graphical device does not seem to be a standard screen device. Hence the\naspect ratio might not be correct. (This message appears only one per session.)")
    }
  }
}


ArrangeDevice <- function(graphics, figs, dh=2.8, h.outer=1.2,
                          dw = 2.5, w.outer=0.7) {
  if (is.na(graphics$always_open_screen)) {
      open <- interactive()
      RFoptions(graphics.always_open_screen = open)
  } else open <- graphics$always_open_screen

  if (graphics$always_close_screen) {
    close.screen(all.screens=TRUE)
    if (is.finite(graphics$height) && graphics$height>0) {
      if (length(dev.list()) > 0) dev.off() ## OK
    }
  }
  
  H <-  graphics$height
  if (is.finite(H) && H>0) {
    H <- H * pmin(1, graphics$increase_upto[1] / figs[1],
                  graphics$increase_upto[2] / figs[2])
    DH <- H * dh / (dh + h.outer)
    HO <- H - DH
    curH <- figs[1] * DH + HO
    W <- H * (dw + w.outer) / (dh + h.outer)
    DW <- W * dw / (dw + w.outer)
    WO <- W - DW
    curW <-  figs[2] * DW + WO
    if (open) ScreenDevice(height=curH, width=curW)
    return(c(curH, curW)) 
  } else {
    return(rep(NA, 2))
  }
}


#restore_par <- function(oldpar) {
#  do.call(graphics::par, oldpar)
#  graphics::par(cex = oldpar$cex) ## muss extra aufgerufen werden. noch mehr?
#}
