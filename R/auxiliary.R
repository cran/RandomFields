
Dev <- function(on, dev, ps=NULL, cur.cex=TRUE, paper="special",
                width=5, height=5, quiet=FALSE, ...){
  ## function to handle output device:
  ##   on: T=output device is activated; F=device will be closed
  ##   dev: device number or TRUE= postscript or FALSE=pdf
  ##        or character (then the name is interpreted as function name)
  ##   ps : postscript file name; only needed when dev is logical
  ##
  ##  Dev takes over the par.options of previous plots
  if (on) {
    par.orig <- if (!is.null(dev.list())) par(no.readonly=TRUE) else NULL
    if (exists(".dev.orig", envir=.GlobalEnv)) {
      warning("Dev has been still open (.dev.orig exists). Closing.")
      if (!is.null(try(Dev(FALSE, .dev.orig$dev))))
        rm(.dev.orig, envir=.GlobalEnv)
    }
    if ((cur.cex <- cur.cex && !is.null(par.orig)) && !is.null(dev.list()))
      par.orig <- par(no.readonly=TRUE)
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
          if (any(l <- (splt <- rev(strsplit(ps,"\\.")[[1]]))[1] == exts)){
            dev <- !l[3]
            ps  <- paste(rev(splt[-1]), collapse=".")
            ext <- exts[l][1]
          }
          if (dev) {
            fn <- paste(ps,ext,sep=".")
            postscript(fn, horizontal=FALSE, paper=paper, width=width,
                       height=height, ...)
          } else {
            fn <- paste(ps,".pdf",sep="")
            pdf(fn, horizontal=FALSE, width=width, height=height, ...)
          }
        } else { # character
          fn <-  paste(ps,".", dev, sep="")
          txt <- paste(dev,"('", fn, "',width=width,height=height,...)", sep="")
          eval(parse(text=txt))
        }
        if (!quiet) cat("creating", fn, "\n")
       # if (cur.cex) par(par.orig)
      }
    } else {
      if (dev %in% dev.list()) dev.set(dev)
      else {
        # print(c(height, width))
        stopifnot(is.finite(height+width))
        get(getOption("device"))(height=height, width=width)
      }
      keep <- dev < 3
    }
    assign(".dev.orig", 
           list(dev.prev=devPrev, dev.cur=dev.cur(), keep=keep),
           envir=.GlobalEnv)
  } else { # off
    if (!exists(".dev.orig", envir=.GlobalEnv)) stop("Dev is not open")
    if (dev.cur() != .dev.orig$dev.cur) {
      warning("Dev is not the currently active device")
      dev.set(.dev.orig$dev.cur)
    }
    if (length(dev.list())>0)
      if (.dev.orig$keep) par(new=FALSE) else dev.off()
    dev.set(.dev.orig$dev.prev)
    rm(.dev.orig, envir=.GlobalEnv) 
  }
  invisible(NULL)
};
	
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
    if (PrintLevel>1) cat("'", file, "' already exists.\n");
    return(1)
  } else { 
    LockFile <- paste(file, lock.ext, sep="")
    if (file.exists(LockFile)) { #2.
      if (PrintLevel>1) cat("'",file,"' is locked.\n");
      return(2);
    }
    PID<-pid();
    ## to do : dir<-strsplit(file,"/") ## check if directory exists !!
    ##         otherwise write will fail.
    write(file=LockFile,c(PID,hostname()),ncolumns=2,append=TRUE); #3.a.
    Pid <- matrix(scan(LockFile,what=character(0)),nrow=2)
    if ((sum(Pid[1,]==PID)!=1) || (sum(Pid[1,]>PID)>0)){ #3.b.
      if (PrintLevel>1) cat("Lock file of '", file, "' is knocked out.\n");
      return(3);
    }
  }
  return(0);
}

LockRemove <- function(file) {
  ## removes auxiliary files created by FileExists
  lock.ext <- ".lock";
  system(paste("rm -f ",file,lock.ext,sep=""))
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

