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

useraction <- function(action=c("none", "start.register",
                          "continue.register",
                          "replay", "endless",
                          "delete"), sleep, wait, PrintLevel,
                       actionlist) {
  ## define how the system should behave -- shall it register
  ## the mouse clicks and readlines or replay it?
  ##
  ## none* : neither registrating nor replaying
  ## start.register : start register
  ## continue.register* : if registration was interrupted by "none"
  ## replay* : replay once
  ## endless* : replay again and again
  ## delete: delete what have been stored
  ## * : only these values are used globally
  ##
  ## wait : regulates the blinking speed of the crosses (and of readline)
  ## sleep : sleeping time before the blinking starts
  action <- match.arg(action)
  if (.Platform$OS.type != "unix" && action %in% c("replay", "endless")) {
    warning("replay does not work correctly on non-unix systems")
  } 
  assign(".action.mismatch", "user input replay mismatch",
         envir=.RandomFields.env)
  assign(".action.mode", action, envir=.RandomFields.env)
  Set <- function(name, value, default)
    if (!missing(value)) assign(name, value, envir=.RandomFields.env)
    else if (!exists(name, envir=.RandomFields.env))
      assign(name, default, envir=.RandomFields.env)
  Set(".action.sleep", sleep, 2.5)
  Set(".action.wait", wait, 0.1)
  Set(".action.print", PrintLevel, 0)
  Set(".action.list", actionlist, list())
  Set(".action.n", default=Inf)
  switch(action,
         none = {},
         start.register = {
           assign(".action.list", list(), envir=.RandomFields.env);
           assign(".action.mode", "continue.register", envir=.RandomFields.env)
         },
         continue.register = {},
         replay = {
           assign(".action.n", 1, envir=.RandomFields.env)
         },
         endless = {
           assign(".action.n", 1, envir=.RandomFields.env)
         },
         delete = {
           assign(".action.list", list(), envir=.RandomFields.env)
           assign(".action.mode", "none", envir=.RandomFields.env)
         }
         )
  invisible(NULL)
}

getactionlist <- function() {
  if (exists(".action.list", envir=.RandomFields.env))
    get(".action.list", envir=.RandomFields.env)
  else stop("action list not initialised")
}

getactions <- function() {
  if (!exists(".action.list", envir=.RandomFields.env))
    stop("action list not initialised")
  list(mismatch=get(".action.mismatch", envir=.RandomFields.env),
       sleep=get(".action.sleep", envir=.RandomFields.env),
       wait=get(".action.wait", envir=.RandomFields.env),
       mode=get(".action.mode", envir=.RandomFields.env),
       actions=get(".action.list", envir=.RandomFields.env),
       n=get(".action.n", envir=.RandomFields.env),
       print=get(".action.print", envir=.RandomFields.env)
       )
}

putactions <- function(l) {
  formerstate <- getactions()
  assign(".action.mismatch", l$mismatch, envir=.RandomFields.env)# only an error
  ##                                                               message
  assign(".action.sleep", l$sleep, envir=.RandomFields.env)
  assign(".action.wait", l$wait, envir=.RandomFields.env)
  assign(".action.list", l$actions, envir=.RandomFields.env) # list, each 
  ##        component is itself a list containing a string (readline) or a list
  ##        of points (locator) and, maybe, a tracing information, see also the
  ##        function userinput 
  assign(".action.n", l$n, envir=.RandomFields.env)
  assign(".action.mode", l$mode, envir=.RandomFields.env)
  assign(".action.print", l$print, envir=.RandomFields.env)
  invisible(formerstate)
}

userinput <- function(fct, info=NULL, prompt, n,  type="n", pch=par()$pch,
                      cex=par()$cex, ...) {
  ## auxiliary function called by Locator and Readline for all types modes:
  ## none, register, replay, ...
  ##
  ## fct : currently only locator and readline are allowed
  ## info: intended for the documentation of the .action.list.
  ##        Background: .action.list gets quickly a very long list.
  ##                    So editing the list manually (you will do it at
  ##                    one stage!) without landmarks will not be too funny
  ##        So, if info!=NULL, its value is stored in list.
  ## prompt : usual parameter of the function readline
  ## n, type, pch, cex,... : usual parameters of the function `locator'
  
  mode <- if (exists(".action.mode", envir=.RandomFields.env))
    get(".action.mode", envir=.RandomFields.env) else "none"
  if (mode %in% c("none", "continue.register")) {
    l <- switch(fct,
                locator = locator(n=n, type=type, pch=pch, cex=cex, ...),
                readline = readline(prompt=prompt)
                )
    if (mode=="continue.register")
      assign(".action.list",
             append(get(".action.list", envir=.RandomFields.env),
                    list(c(l, if (is.null(info)) list() else list(info=info)))),
                    envir=.RandomFields.env)
    return(l)
  } else { ## endless or replay
    action.n <- get(".action.n", envir=.RandomFields.env)
    action.list <- get(".action.list", envir=.RandomFields.env)
    action.wait <- get(".action.wait", envir=.RandomFields.env) * 1000
    action.sleep <- get(".action.sleep", envir=.RandomFields.env) * 1000
    
    if (action.n <= length(action.list)) {
      l <- action.list[[action.n]]
      if (get(".action.print", envir=.RandomFields.env) > 1) {
        cat("action #", action.n, sep="")
        if (!is.null(l$info)) cat(" (", l$info, ")", sep="")
        cat("\n")
      }
      assign(".action.n",
             if (action.n==length(action.list) && mode=="endless") 1
             else action.n + 1, envir=.RandomFields.env)
            
      flash <- function(fct, x, y, wait=action.wait,
                        col=rep(c("white", "red"), 4), ...) {
        for (c in col) {
          fct(x, y, col=c, ...)
          sleep.milli(wait)
        }
      }
      
      switch(fct,
             locator = {
               l <- l[1:2]
               sleep.milli(action.sleep)
               if (is.null(l$x)) {
                 for (col in rep(c( "black", "white", "red"), 3)) {
                   text(mean(axTicks(1)), mean(axTicks(2)), "QUIT", col=col,
                        cex=2) 
                   sleep.milli(action.wait)
                 }
                 sleep.milli(action.sleep)
                 return(NULL)
               }
#               if (!is.list(l)) stop(.action.mismatch)
               switch (type,
                       n = {
                         for (i in 1:length(l$x))
                           flash(points, l$x[i], l$y[i], pch="X", cex=2 * cex,
                                 type="p",...)
                       },
                       p = {
                         for (i in 1:length(l$x))
                           flash(points, l$x[i], l$y[i], pch=pch, cex=2 * cex,
                                 type="p", ...)
                       },
                       l = {
                         if (length(l$x)==1)
                           flash(points, l$x, l$y, type="p", ...) 
                         else for (i in 2:length(l$x))
                           flash(lines, l$x[(i-1):i], l$y[(i-1):i],
                                 wait=action.wait / 3, ...)
                       },
                       o =  {
                         flash(points, l$x, l$y, pch=pch, type="p") 
                         if (length(l$x)>1)
                           for (i in 2:length(l$x)) {
                             flash(points, l$x[i], l$y[i], type="p",
                                   wait=action.wait / 5, ...)
                             flash(lines, l$x[(i-1):i], l$y[(i-1):i],
                                   wait=action.wait / 5, ...)
                           }
                       })
             },
             readline = {
               if (!is.character(l[[1]]))
                 stop(get(".action.mismatch", envir=.RandomFields.env))
               l <- l[[1]]
               wait <- 3 * action.wait
               cat(prompt, "")
               for (i in 1:nchar(l)) {
                 sleep.milli(wait)
                 cat(substring(l[1], i, i))
               }
               sleep.milli(wait)
               cat("\n")
             }
             )
      return(l)
    } else {
      cat("warning: end of action list; switching to terminal mode action=`none'\n")
      assign(".action.mode", "none", envir=.RandomFields.env)
      userinput(fct, n=n, prompt=prompt, type=type, pch=pch, cex=cex, ...)
    }
  }
}

Locator <- function(n, type="n", info=NULL, ...)
  userinput("locator", info=info, n=n, type=type, ...)

Readline <- function(prompt="", info=NULL)
  userinput("readline", info=info, prompt=prompt)
  

Dev <- function(on, dev, ps=NULL, cur.cex=TRUE, paper="special",
                width=5, height=5, quiet=FALSE, innerwidth, innerheight,
                mai, ...){
  
  # print(ls(all=TRUE, envir=.RandomFields.env))
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
        rm(.dev.orig, envir=.RandomFields.env)
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
            if (!file.create(fn)) stop(paste("The file", fn,"cannot be created"))
            postscript(fn, horizontal=FALSE, paper=paper, width=width,
                       height=height, ...)
          } else {
            fn <- paste(ps,".pdf",sep="")
            pdf(fn, horizontal=FALSE, width=width, height=height, ...)
          }
        } else { # character
          fn <-  paste(ps,".", dev, sep="")
          if (!file.create(fn)) stop(paste("The file", fn, "cannot be created"))
         txt <- paste(dev,"('", fn, "',width=width,height=height,...)", sep="")
          eval(parse(text=txt))
        }
        if (!quiet) cat("creating", fn, "\n")
       }
      if (!missing(mai)) par(mai=mai)
    } else {
      if (dev %in% dev.list()) dev.set(dev)
      else {
        # print(c(height, width))
        stopifnot(is.finite(height+width))
        get(getOption("device"))(height=height, width=width)
      }
      keep <- dev < 3
    }
    # if (cur.cex) par(par.orig) # uncommented 12.8.04 + nach unten
    if (cur.cex) { ## komisches Verhalten !! wenn die beiden Befehle
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
    rm(.dev.orig, envir=.RandomFields.env) 
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

