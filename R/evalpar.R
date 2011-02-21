

eval.parameters <- function(variable, entry, update, simulate, dev,
                            create=TRUE, 
                            col.rect="red", col.bg="blue", col.sep="grey",
                            col.left="red", col.mid="white", col.right="white",
                            col.line="red", col.txt="black",
                            cex=0.8, cex.i=cex, sep="-----",
                            ...) {  
  ## since the variables are created within the function with names
  ## given by the user, all local variables used in this functions start with
  ## points
  ##
  ## the first variable in ... is the list for the variables in the menue
  ## additional parameters are parameters of the function 'simulate'
  ## usually theses additional parameters control the behaviour of simulate,
  ## e.g. partial simulation, instead of complete simulaton.
  ## the parameter values given are the standard parameters; these standard
  ## parameters can be locally be overwritten by special menue point definitions
  ##
  ## eval.parameters returns the first variable in the '...' list
  
  entry <- rev(entry)
  .length <- length(entry) + 1 
  .bg <- par()$bg
  par(bg="white")
  show.parameters <- function() {
    ## the function plots the menu
    screen(dev)
    par(mar=c(0,0.2,0,0))
    plot(Inf, Inf,
         xlim=c(0,1.115), ylim=c(0.8, length(entry)+0.5), 
         axes=FALSE, xlab="", ylab="")
    ## note! if x is changed then also at !e$delta
    for (.j in 1:length(entry)) {
      .e <- entry[[.j]]
      ## menue points may appear according to conditions
      if (!is.null(.e$cond) && !eval(parse(text=paste(variable, "$", .e$cond))))
        next
      if (is.null(.e$var)) {
        ## only text is shown
        .t.cex <- if (is.null(.e$cex)) cex else .e$cex       
        .lab <- if (is.expression(.e$name)) .e$name else paste(sep,.e$name,sep)
        .txt.x <- 0.5
        .adj.x <- 0.5
        if (strwidth(.lab, "figure", cex=.t.cex)>1) {
          .txt.x <- 0
          .adj.x <- 0
        }
        text(.txt.x, .j, lab=.lab,
             col=if (is.null(.e$col)) col.sep else .e$col,
             cex=.t.cex, adj=c(.adj.x, 0.5))
        next
      }
      ## a rectangle for the free input is shown
      rect(1.1, .j, 1.15, .j + 0.25, col=col.rect, xpd=NA)
      if (is.null(.e$val)) {
        ## menue point for input of a string
        text(0.5, .j+0.25,
             lab=paste(.e$name, ":", eval(parse(text=paste(variable,
                                                  "$",.e$var)))),
             adj=c(0.5, 0), cex=cex)
        rect(0, .j, 1, .j + 0.25, col=col.rect)
        next
      }
      .col <- if (is.null(.e$col)) col.bg else .e$col
      rect(0, .j, 1, .j + 0.25, col=.col) ## a menue bar
      ## either incremental, or absolute vaue, or logical, or discrete
      if (is.logical(.e$val)) {
        ## logical menue bar
        ## programmer may define the value of the variable being
        ## the not-value of the shown value. To manage this,
        ## xor is used in the next line
        .f <- xor(!.e$val, eval(parse(text=paste(variable, "$",.e$var))))
        if (is.expression(.e$name))
          .lab <- eval(parse(text=paste("expression(", as.character(.e$name),
                              "==",.f,")")))        
        else  .lab <- paste(.e$name,"=", .f)
        text(0.5, .j+0.25, lab=.lab, adj=c(0.5, 0), cex=cex)
        text(x=0.5, y=.j + 0.125, lab="negation", cex=cex.i,
             col=col.mid, adj=c(0.5,0.5))
        next
      }
      if (is.function(.e$val)) {
        ## either incremental or absolute value meanu bar
        .x <- c(0,0.25,0.5,0.75,1)
        .f <- eval(parse(text=paste(variable, "$", .e$var)))
        if (is.expression(.e$name))
          .lab <- eval(parse(text=paste("expression(", as.character(.e$name),
                              "==", format(.f, digits=3),")")))        
        else  .lab <- paste(.e$name,"=", format(.f, digits=3))
        text(0.5, .j + 0.25, adj=c(0.5, 0), lab=.lab, col=col.txt, cex=cex)
        .lab <- character(length(.x))
        .lab <- formatC(.ff <- .e$val(.x), digits=3, width=1)
        if (.e$delta) {
          ## incremental
          .lab[.ff==0] <- paste("", .lab[.ff==0], "+", sep="") ## "-" !!!??
          .lab[.ff>0] <- paste("+",.lab[.ff>0], sep="")
        }
        text(x=.x, y=.j, adj=c(0.5, 1), lab=.lab, col=col.txt, cex=cex.i)
        if (!.e$delta)
          ## absolute values: a vertical bar is shown at the value
          lines(x=rep(.f,2)/.ff[5], y=c(.j, .j + 0.25), col=col.line, lwd=2)
        next
      }
      if (is.character(.e$val)) {
        ## discrete variables: show current value on top,
        ## previous and next values in menue bar
        .evar <- eval(parse(text=paste(variable, "$",.e$var)))
        .left <- if (.evar>1) substr(.e$val[.evar - 1], 1, 13) else ""
        .right <-
          if (.evar<length(.e$val)) substr(.e$val[.evar + 1], 1, 13) else ""
        
        text(0.5, .j + 0.25, adj=c(0.5, 0), col=col.txt, cex=cex,
            lab=paste(.e$name,":",
              .e$val[eval(parse(text=paste(variable, "$",.e$var)))]))
        text(x=c(0.0, 0.5), y=.j + 0.125, adj=c(0,0.5),
             col=c(col.left,col.right), lab=c(.left,.right), cex=cex.i)
        next
      }
      stop("unknown entry type")
    }
  }

  .namen <- names(list(...))
  eval(parse(text=paste(.namen,"<-","list(...)[[",1:length(.namen),"]]")))
  
  ## eventuell ""h"" herausnehmen da sonst staendig kopiert
  ## hier: so programmiert, dass h irgendwo sein, kann. Weiter unten
  ## wird gefordert, dass h der erste zusaetzliche Parameter ist !!!!
  .namen.orig <- .namen
  .name.main <- is.na(pmatch(.namen.orig, variable))
  if (sum(!.name.main)==1) .namen.orig <- .namen.orig[.name.main]

  ## store standard values of the additional parametes given
  ## these values might be overwritten if special menue points
  ## are chosen (e.g. partial simulation instead of complete one)
  if (length(.namen.orig)>0) {
    eval(parse(text=paste(".orig.",.namen.orig,"<- list(...)[[",
                 (1:length(.namen))[.name.main],"]]", sep="")))
    .set.default <- paste(.namen.orig, "<- .orig.", .namen.orig, sep="")
  } else .set.default <- ""
  
  for (j in 1:length(entry)) {
    ## check that names are not vector of characters
    .e <- entry[[j]]    
    if (!is.character(.e$name) && !is.expression(.e$name)) {
      stop("'name' of entry ", length(entry) + 1 - j,
                 " is not a character of length 1")
    }
    if (!is.character(.e$var) && !is.null(.e$var)) {
      stop("'var' of entry ",  length(entry) + 1 - j,
                 " is not a character of length 1")
    }
  }   

  if (create) {
    ## if for some of the menu points values are not given yet
    ## some initial values are created, if create==TRUE,
    ## if there are missing values and create==FALSE, there will
    ## be an error later on -- anywhere.
    for (j in 1:length(entry)) {
      .e <- entry[[j]]      
      if (is.null(.e$var)) next
      if (is.null(eval(parse(text=paste(variable,"$",.e$var))))) {
        cat("Creating `",variable,"$",.e$var,"' ...\n", sep="")
        if (is.logical(.e$val)) {
          eval(parse(text=paste(variable, "$",.e$var,"<- TRUE")))
        } else
        if (is.function(.e$val)) {
          eval(parse(text=paste(variable, "$", .e$var, "<- .e$val(.5)")))     
        } else
        if (is.character(.e$val)) eval(parse(text=paste(variable, "$", .e$var,
                                               "<- 1")))
        else if (is.null(.e$val)) eval(parse(text=paste(variable, "$", .e$var,
                                               "<- ''"))) 
      }
    }
  }
  
  show.parameters()
  .simu.txt <- paste(.namen[1], "<- simulate(",
                    paste(paste(.namen,"=",.namen), collapse=","), ")")

  ## .history contains last changed value; it is erased if simulate is called
  .zaehler <- 0

  eval(parse(text=paste(variable, "$.history <- list()")))
  while (!is.null(.loc <- Locator(1, info=variable))){
    .j <- as.integer(.loc$y + 0.2)
    if ((.j>=1) && (.j<=length(entry)) && ((.loc$y+0.2) %% 1<0.8) &&
        (.loc$x>=-0.1) && (.loc$x<=1.15)){
      .e <- entry[[.j]]      
      if (!is.null(.e$cond) && !eval(parse(text=paste("(",variable,
                                             "$",.e$cond,")"))))
        next
      if (is.null(.e$var)) {
        if (!is.null(.e$val)) {
          .zaehler <- .zaehler + 1
          eval(parse(text=paste(variable,
                       "$.history[[.zaehler]] <- list(.length - .j, .e$val)")))
          ## call function simulate or a function defined in $val 
          if (.e$val=="simulate") {
            ## to do: define "undo" and "reundo" 
            eval(parse(text=.set.default))
            eval(parse(text=paste(names(.e), "<- .e$", names(.e))))
            eval(parse(text=.simu.txt), envir=environment())
          } else {
            eval(parse(text=.e$val))
          }
        }
        show.parameters()
        next
      }
      .noupdate <- FALSE
      .evar <- eval(parse(text=paste(variable, "$", .e$var)))
      if ((.loc$x>1.1) || is.null(.e$val)) {
        ## free input on xterm
        .newvalue <-
          Readline(prompt=paste(.e$name, " (", .evar , ") : ", sep=""),
                   info=paste(variable, .e$name))
        if (.newvalue=="") next
        if (.newvalue=="exit immediately") stop("user forced exit")
        else {
          if (is.logical(.e$val)) .newvalue <- as.logical(.newvalue) else
          if (is.function(.e$val)) .newvalue <- as.numeric(.newvalue) else
          if (is.character(.e$val)) .newvalue <- as.integer(.newvalue)
          eval(parse(text=paste(variable, "$", .e$var, "<- .newvalue")))
        }
      } else {
        if (is.character(.e$val)) {
          .newvalue <- .evar
          if (.loc$x < 0.5) {
            if (.newvalue > 1) .newvalue <- .newvalue - 1
          } else {
            if (.newvalue < length(.e$val)) .newvalue <- .newvalue + 1
          }
          eval(parse(text=paste(variable, "$", .e$var, "<- .newvalue")))
        } else
        if (is.logical(.e$val)) {
          .newvalue <- eval(parse(text=paste(variable, "$", .e$var, "<- !",
                                    variable, "$", .e$var)))
        } else
        if (is.function(.e$val)) {
         .newvalue <- eval(parse(text=paste(variable, "$", .e$var,
                                   "<-.e$val(.loc$x,",variable, "$",.e$var,")")))
        }
      } # else
      .zaehler <- .zaehler + 1
      eval(parse(text=paste(variable, "$.history[[.zaehler]] <- ",
                   "list(.length -.j, .evar, .newvalue)")))
      if ((is.null(.e$update) && update && !.noupdate) ||
          (!is.null(.e$update) && .e$update)) {
        screen(dev)
        eval(parse(text=.set.default))
        eval(parse(text=paste(names(.e), "<- .e$", names(.e))))
        eval(parse(text=.simu.txt), envir=environment())
      }
      show.parameters()
    }
  }
  
  screen(dev)
  plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
  par(bg=.bg)
  ## new6.4.04. orig:      return(eval(parse(text=names(list(...))[1])))
  return(eval(parse(text=strsplit(strsplit(variable, "\\$")[[1]][1],
                      "\\[")[[1]][1])))
}


