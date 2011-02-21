if (FALSE) {
  library(RandomFields, lib="~/TMP"); source("~/R/RF/RandomFields/R/ShowModels.R");source("~/R/RF/RandomFields/R/rf.R");source("~/R/RF/RandomFields/R/auxiliary.R");source("~/R/RF/RandomFields/R/evalpar.R"); options(locatorBell=FALSE); RFparameters(Print=1); x <- seq(1,10,0.1); user.action("start"); ShowModels(x=x); getactionlist()

source("~/R/RF/RandomFields/R/ShowModels.R");source("~/R/RF/RandomFields/R/rf.R");source("~/R/RF/RandomFields/R/auxiliary.R");source("~/R/RF/RandomFields/R/evalpar.R"); user.action("replay"); ShowModels(x=x)
}

# putactionlist(getactionlist()[1:(length(getactionlist())-1)])


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

    par.orig$new <- par.orig$fin <- par.orig$mai <- par.orig$pin <-
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
      if (dev %in% dev.list()) dev.set(dev)
      else {
        stopifnot(is.finite(height+width))
        do.call(getOption("device"), list(height=height, width=width))
      }
      keep <- dev < 3
    }

    if (cur.cex) par(par.orig) # uncommented 12.8.04 + nach unten
    if (cur.cex && FALSE) { ## komisches Verhalten !! wenn die beiden Befehle
      ##              zusammengefasst werden (gekippte eps in Latex)
      ##              29.5.05
      par(par.orig[39]) # $mfg
      par(par.orig[-39])# uncommented 12. 8.04 + nach unten
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
      else dev.off()
    if ((devPrev <- get(".dev.orig", envir=.RandomFields.env)$dev.prev) != 1)
      dev.set(devPrev)
    rm(".dev.orig", envir=.RandomFields.env) 
  }
  invisible(NULL)
}
	
hostname<-function(){.C("hostname", h=paste(seq(0,0,l=100), collapse=""),
                        as.integer(100), PACKAGE="RandomFields")$h}

pid <- function() {.C("pid", i=integer(1), PACKAGE="RandomFields")$i}


FileExists <- function(file, PrintLevel=RFparameters()$Print) {
  ## for parallel simulation studies: the same data output file should not
  ## be created twice. So:
  ## 1. if file exists then assume another process has done the work already
  ## 2. if file.lock existss then assume another process is doing the work
  ## 3.a. otherwise create file.lock to show other processes that the process
  ##      is going to do the work
  ## 3.b. check if another process has started with the same work at the same
  ##      time it may happen that in case of simulatenous creation of file.lock
  ##      no process will do the work...(then the lock file will rest.)
  lock.ext <- ".lock";
  if (file.exists(file)) { #1.
    if (PrintLevel>2) cat("'", file, "' already exists.\n");
    return(1)
  } else { 
    LockFile <- paste(file, lock.ext, sep="")
    if (file.exists(LockFile)) { #2.
      if (PrintLevel>2) cat("'",file,"' is locked.\n");
      return(2);
    }
    PID<-pid();
    ## to do : dir<-strsplit(file,"/") ## check if directory exists !!
    ##         otherwise write will fail.
    write(file=LockFile,c(PID,hostname()),ncolumns=2,append=TRUE); #3.a.
    Pid <- matrix(scan(LockFile,what=character(0), quiet=TRUE),nrow=2)
    if ((sum(Pid[1,]==PID)!=1) || (sum(Pid[1,]>PID)>0)){ #3.b.
      if (PrintLevel>2) cat("Lock file of '", file, "' is knocked out.\n");
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
  CIRCLE.X <- cos(seq(0,2*pi,l=20))
  CIRCLE.Y <- sin(seq(0,2*pi,l=20))
  circle <- function(x,r) { polygon(x[1]+ r* CIRCLE.X,x[2]+ r* CIRCLE.Y,
                                    col=fill, border=col) }
  ##r <- x$NormedData - min(x$NormedData) +1
  ##r <- r/max(r)/nrow(x$coord) * diff(xlim) * diff(ylim) * 2.5;
  maxr <- max(data[,3])
  plot(Inf, Inf, xlim=xlim, ylim=ylim, xlab="", ylab="",...)
  for (i in 1:nrow(data)) { circle(data[i,c(1,2)], factor*data[i,3]) }
}


Print <- function(...) {
  l <- list(...)
  n <- as.character(match.call())[-1]
  cat("\n\n")
  for (i in 1:length(l)) {
    if (!is.list(l[[i]]) && is.vector(l[[i]])) cat("\n",n[i], "=", l[[i]])
    else {
      cat("\n", n[i], "= ")
      if (is.list(l[[i]])) {
        cat("  ")
        str(l[[i]])
      } else {
        cat("\n")
        print(l[[i]])
      }
    }
  }
  cat("\n")
}
