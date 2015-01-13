
## meta-function that generates functions like RMplus, RMwhittle, which
## the user will use to generate explicit covariance models, i.e. objects
## of class 'RMmodels'

#   rfGenerateModels(TRUE)
param.text.fct <- function(catheg, names, havedistr=TRUE){
  x <- paste("if (hasArg(", names, ") && !is.null(subst <- substitute(", names,
             "))) {\n    u <- try(is.numeric(",
             names,
             ") || is.logical(", names,
             ") || is.language(", names,
             ")\n\t || is.list(", names,
             ") || is(", names,
             ", class2='RMmodel'), silent=TRUE)\n",
             "    if (is.logical(u) && u) ", catheg, "[['", names,
             "']] <- ", names,
             "\n    else if (substr(deparse(subst), 1, 1)=='R') ", catheg,
             "[['", names,  "']] <- ", names,
             "\n    else ", sep="")
  if (havedistr) {
    paste(x, catheg, "[['", names,
          "']] <- do.call('", ZF_DISTR[1], "', list(subst))\n  }", sep="")
  } else {
    paste(x, "stop('random parameter not allowed')\n  }");
  }
}

rfGenerateModels <- function(assigning,
                             RFpath = "~/R/RF/svn/RandomFields",
                             RMmodels.file = paste(RFpath, "RandomFields/R/RMmodels.R",
                               sep="/")
                             ) {
  
  # if file already exists, remove it.
  if (assigning && file.exists(RMmodels.file))
    file.remove(RMmodels.file)
#  Print(RMmodels.file, file.exists(RMmodels.file)) #
#  Print(RMmodels.file)
  
  write(file = RMmodels.file, append = TRUE,
        "\n## This file is created automatically by 'rfGenerateModels'.\n\n")

  ## defined constants
  diminf <- 999999
  
  # define empty strings
  empty <- paste(rep(" ", MAXCHAR), collapse="")
  empty2 <- paste(rep(" ", MAXCHAR), collapse="")
  # inialized attribute parameter
  
  nr <- GetCurrentNrOfModels(TRUE)
  va <- nr * MAXVARIANTS

  # get attribute parameter
  A <- .C("GetAttr", nr=integer(va), type=integer(va), operator=integer(va),
          monotone=integer(va), finiterange=integer(va),
          simpleArguments=integer(va), internal=integer(va),
          domains=integer(va), isos=integer(va),
          maxdim=integer(va), vdim=integer(va),
          includevariants= as.integer(TRUE),
          va = integer(1),
          PACKAGE="RandomFields")
  va <- A$va
  
#  idx <- .C("GetModelList", idx=integer(nr^2), as.integer(TRUE),
#            PACKAGE="RandomFields")$idx
 # dim(idx) <- c(nr, nr)
#  Print(A)

  i <- 1
  while (i <= va) {
    step <- 1
    ## sequential steps for each model
    
    if (A$internal[i]) {
        cat(i, "internal", .C("GetModelName",as.integer(A$nr[i]),
                 name=empty, nick=empty2, PACKAGE="RandomFields")$name,"\n")
        i <- i + 1
      next
    }
    # get model name
    ret <- .C("GetModelName",as.integer(A$nr[i]),
              name=empty, nick=empty2, PACKAGE="RandomFields")
    nick <- ret$nick
    #cat(i, nick, "\n");

    type <- A$type[i]
    iso <- A$isos[i]
    while (i + step <= va  &&
           nick ==  .C("GetModelName",as.integer(A$nr[i + step]),
               name=empty, nick=empty2, PACKAGE="RandomFields")$nick) {
      cat("...variant added\n")
      type <- c(type, A$type[i + step])
      iso <- c(iso, A$isos[i + step])
      step <- step + 1
    }

   
    finiterange <- as.logical(A$finiterange)
    finiterange[A$finiterange < 0] <- NA
 
    ## get names of submodels
    subname.info<- .Call("GetSubNames", as.integer(A$nr[i]), PACKAGE="RandomFields")
    subnames <- subname.info[[1]]
    subintern <- subname.info[[2]]
    subnames.notintern <- subnames[!subintern]
     
    # get names of  parameters
    paramnames <- .Call("GetParameterNames", as.integer(A$nr[i]),
                        PACKAGE="RandomFields")
    internal <- which(paramnames == INTERNAL_PARAM)
    if (length(internal) > 0) paramnames <- paramnames[-internal]   
    elmnt <- which(paramnames == "element")
    if (length(elmnt) > 0)  {
      stopifnot(length(elmnt) == 1)
      paramnames <- c(paramnames[-elmnt], "element")
    }
    par.intern <- paramnames %in% subnames
 #   Print(subnames, paramnames)
    
    if (any(par.intern)) stop(nick, ": subnames (",
                              paste(subnames, collapse=", "),
                              ") and parameter names (",
                              paste(paramnames, collapse=", "),
                              ") match.")
   
    ex.anysub <- length(subnames)>0
    ex.sub <- length(subnames.notintern)>0
     
    ex.par <- length(paramnames)>0
    ex.std <- ((nick != ZF_DOLLAR[1] && any(isNegDef(type))) ||
               nick == "RMball"
               || nick == ZF_PLUS[1] || nick[1] == ZF_MULT[1]) &&
                 nick != "RMtrafo"
    
    
    cat(i, "\t", nick, ",\t\told name ", ret$name, "\t", ex.std, "\t",
        paste(type, collapse="/"), "\n", sep="")
    
    if(nick == ZF_DOLLAR[1]){ 
      text.fct.head <-
        paste(nick, " <- function(phi, var, scale, Aniso, proj, anisoT)")
    } else {
      text.fct.head <-
        paste(nick,
              " <- function(",
              if (ex.sub) {
                paste(paste(subnames.notintern, collapse=", "), sep="")
              },
              if (ex.sub && (ex.par || ex.std)) ", ",
              if (ex.par) {
                paste(paste(paramnames, collapse=", "), sep="")
              },
              if (ex.par && ex.std) ", ",
              if (ex.std){
                "var, scale, Aniso, proj"
              }, 
              ")",
              sep="")
    }

     
    if (ex.par) {
      par.body <- param.text.fct("par.model", paramnames,
                                 any(isNegDef(type)) || any(type==ShapeType))
      idx <- paramnames == 'envir'
      if (any(idx))
        par.body[idx] <-
          "par.model[['envir']] <- if (hasArg(envir)) envir else new.env()"
    } else par.body <- NULL
 
    text.fct.body <-
      paste("{\n  ",
            "cl <- match.call()",
            "\n  ",
            "submodels <- par.general <- par.model <- list() \n  ",
            ## get submodels
            if (ex.anysub) {
              paste("if (hasArg(", subnames, ")) submodels[['", subnames,
                    "']] <- ", subnames, sep="", collapse="\n  ")
            },
            if (ex.anysub) "\n  ",
            "\n  ",
            ## get model specific parameter
            if (ex.par) paste(par.body, collapse="\n  "),
            if (ex.par) "\n  ",
            ## get general model parameter
            if (ex.std) {
              paste(param.text.fct("par.general",
                                   c("var", "scale", "Aniso", "proj")),
                    collapse="\n  ")
            },
                      "\n  ",
             # create RMmodel object
            "model <- new('", ZF_MODEL, "', ",
            "call = cl, ",
            "name = ", "'", nick, "'", ", \n  \t\t",
            "submodels = submodels, ",   "\n  \t\t",
            "par.model = par.model, ",
            "par.general = par.general)",

            "\n  ",
            "return(model)\n}\n",
            sep=""
            )

    text.fct <- paste(text.fct.head, text.fct.body)
    
    # assign class 'RMmodelgenerator' (ZF_MODEL_FACTORY) and attributes like stationarity
    # to the function:

     text.assign.class <-
      paste(nick, " <- new('", ZF_MODEL_FACTORY,                  "',", "\n\t",
         ".Data = ",        nick,                                  ",", "\n\t",
         "type = ", "c('", paste(TYPENAMES[type+1], collapse="', '"), "'),",
            "\n\t",
         "isotropy = ", "c('", paste(ISONAMES[iso+1], collapse="', '"), "'),",
            "\n\t",
         "domain = ",       "'", DOMAIN_NAMES[A$domains[i]+1],    "',", "\n\t",
         "operator = ",     as.logical(A$operator[i]),             ",", "\n\t",
         "monotone = ",    "'", MONOTONE_NAMES[A$monotone[i] + 1 - MISMATCH],
                                                                  "',", "\n\t",
         "finiterange = ",  finiterange[i],          ",", "\n\t",
         "simpleArguments = ",  as.logical(A$simpleArguments[i]),  ",", "\n\t",
         "maxdim = ", if(A$maxdim[i]>diminf) Inf else A$maxdim[i], ",", "\n\t",
         "vdim = ",         A$vdim[i],                                  "\n\t",
         ")",
         sep="")
 
    text <- paste(text.fct, "\n", text.assign.class, "\n\n\n", sep="")
  
    if (assigning) {
      #sink(file = RMmodels.file, append = TRUE, type='output')
      write(file = RMmodels.file, append = TRUE, text)
      #cat(text)
      #sink()
      #unlink(RMmodels.file)
    }
    i <- i + step 
  }  ## matches for (i in 1:nr) {
 


  # if help page to the function does not exist, throw warning
  if (length(as.character(help(nick))) == 0) {
    if (file.exists("/home/schlather/R/RF/RandomFields/R/rf.RXX")||
        file.exists("do.not.rm.this.file")) {
      if (!any(nick == c("list of exceptions"))) {
        warn <- paste("Warning: help page for '", nick,"' does not exist.",
                      sep="")
        cat(warn, "\n")
      } 
    }
  }
  invisible()
}


kind <- function(Zeilen, i, start, cont="", ignore=" ", endofname=" ") {
 # Print(Zeilen[i], substr(Zeilen[i], 1, nchar(start)) == start, substr(Zeilen[i], 1, nchar(start)) , start)
  if (substr(Zeilen[i], 1, nchar(start)) == start) {
    s <- substring(Zeilen[i], nchar(start) + 1)
 #   Print(s)
    j <- 2
    while (j <= nchar(s) && substr(s, j, j) != endofname) j <- j + 1
    stopifnot(j <= nchar(s)) 
    name <- substr(s, 1, j -1)
    u <- strsplit(substring(s, j + 1), "//")[[1]]
    kommentar <- paste(u[-1], collapse = " ")
    RC <- nchar(strsplit(kommentar, "RC")[[1]][1]) < nchar(kommentar)
    value <- u <- u[1]
    i <- i + 1
    if (any(cont != "")) {
      repeat {
        j <- nchar(u)
        if (ignore != "") while (substr(u, j, j) %in% ignore) j <- j - 1
        if (substr(u, j, j) %in% cont) {
          u <- strsplit(Zeilen[i], "//")[[1]][1]
          value <- paste(value, u)
          i <- i + 1
        } else break
      }
    }
    res <- list(name=name, value=value, RC=RC, i=i)
 #   Print(res)
    return(res)
  } else return(NULL)
}

clean <- function(x) paste(strsplit(paste(strsplit(x, "\t")[[1]], collapse=""), " ")[[1]], collapse="")

CC <- function(x) {
  if (is.numeric(x)) {    
    Real <- TRUE
    Integer <- x == as.integer(x)
    y <- x
  } else {
    stopifnot(is.character(x))
    if (substr(x, 1, 1) =='"') return(x)
    warn <- options()$warn
    options(warn = -1)
    y <- try(as.numeric(x), silent=TRUE)
#    Print(x, y)  
    options(warn = warn)
    if (Real <- !class(y) == "try-error") {
      Integer <- nchar(strsplit(x, "\\.")[[1]][1]) == nchar(x)
    } else Real <- Integer <- FALSE
    y <- paste(strsplit(paste(strsplit(x, "\t")[[1]], collapse=""), " ")[[1]],
               collapse="")
  }
  if (Real) {
    y <- paste("as.", if (Integer) "integer" else "double", "(", y, ")", sep="")
  }
  return(y)
}

rfGenerateConstants <-
  function(RFpath = "~/R/RF/svn/RandomFields",
           RCauto.file = paste(RFpath, "RandomFields/R/RCauto.R", sep="/")
           ) {
  s <- scan(paste(RFpath, "RandomFields/src/AutoRandomFields.h", sep="/"),
            what=character(), sep="\n", blank.lines.skip=FALSE, skip=2)
  #for (i in 1:length(s)) cat(s[i], "\n")
  i <- 1
  typedef <- character(0)
  write(file = RCauto.file,
        "# This file is created automatically by 'rfGenerateConstants'")
  while (i <= length(s)) {    
    if (!is.null(k <- kind(s, i, "#define", "\\", ""))) {
      value <- k$value
      if (length(typedef) > 0)
        for (j in 1:length(typedef)) {
          #Print(value, typedef[j])
          value <- paste(strsplit(value, typedef[j])[[1]], collapse=" ")
         # Print(value)
        }
      v <- clean(k$name)
      w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
      line <- paste(w, v , "\t<-", CC(value))
      write(file = RCauto.file, append = TRUE, line)
    } else if (!is.null(k <- kind(s, i, "typedef enum", c(",", "{"), " "))) {
      value <- strsplit(k$value, "\\{")[[1]]
      typedef <- c(typedef,
                   paste("\\(", strsplit(value[1], ";")[[1]][1] ,"\\)", sep=""))
      value <- strsplit(strsplit(value[2], "\\}")[[1]][1], ",")[[1]]
      for (j in 1:length(value)) {
        v <- clean(value[j])
        w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
        line <- paste(w, v , "\t<-", CC(j - 1))
        write(file = RCauto.file, append = TRUE, line)
      }
      write(file = RCauto.file, append = TRUE, "")
    } else if (!is.null(k <- kind(s, i, "typedef", "\\", ""))) {
      typedef <- c(typedef,
                   paste("\\(", strsplit(k$value, ";")[[1]][1] ,"\\)", sep=""))
    } else if (!is.null(k <- kind(s, i, "extern const", ",", " "))) {
      ## ignored
    } else write(file = RCauto.file, append = TRUE, "")
    i <- if (is.null(k)) i+1 else k$i
  }

  
  s <- scan(paste(RFpath, "RandomFields/src/AutoRandomFields.cc", sep="/"),
            what=character(), sep="\n", blank.lines.skip=FALSE, skip=1)
  #Print(tail(s))
  i <- 1
  while (i <= length(s)) {
    if (!is.null(k <- kind(s, i, "  *", c(",", "{"), " ", endofname="["))) {
      #Print(k)
      value <- strsplit((k$value), "\\{")[[1]]
      value <- strsplit(value[2], "\\}")[[1]][1]
      v <- clean(k$name)
      w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
      line <- paste(w, v, "<-\nc(", value, ")")
     #cat(line, "\n")
      write(file = RCauto.file, append = TRUE, line)
      write(file = RCauto.file, append = TRUE, "")
   } else write(file = RCauto.file, append = TRUE, "")
    i <- if (is.null(k)) i+1 else k$i
  }
 
  return(NULL)
}


rfGenerateTest <- function(files = NULL,
                           RFpath = "~/R/RF/svn/RandomFields/RandomFields") {
   if (length(files) == 0) return()
   for (f in 1:length(files)) {
    cat("creating ", f, ".R", sep="")
    s <- scan(paste(RFpath, "/man/", files[f], ".Rd", sep=""),
              what=character(), sep="\n", blank.lines.skip=FALSE, skip=2)
    i <- 1
    while (i <= length(s) && substr(s[i], 1, 8) != "\\dontrun") i <- i + 1
    i <- i + 1
    if (i <= length(s))
      for (j in i:length(s)) {
        if (substr(s[j], 1, 1) == "%" ) s[j] <- ""
        if (substr(s[j], 1, 1) %in% c("\\", "}") ) {
          j <- j -1
          break;
        }
      }
    if (j >= i) {
      out <- paste(RFpath, "/tests/", files[f], ".R", sep="")
      write(file = out, append = FALSE, "if (RFoptions()$internal$do_tests){")
      write(file = out, append = TRUE, s[i:j])
      write(file = out, append = TRUE, "}")
    }
  }

   cat(out, "created")
}
