

.RandomFields.env <- new.env()

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
  
