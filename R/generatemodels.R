
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
                             RMmodels.file = paste(RFpath, "R/RMmodels.R",
                               sep="/")
                             ) {
 
  # if file already exists, remove it.
  if (assigning && file.exists(RMmodels.file))
    file.remove(RMmodels.file)
  Print(RMmodels.file, file.exists(RMmodels.file)) #
  
  write(file = RMmodels.file, append = TRUE,
        "\n## This file is created automatically by 'rfGenerateModels'.\n\n")

  ## defined constants
  diminf <- 999999
  
  ### get covariance model information from c
  assign(".p", GetrfParameters(TRUE))
  # define empty strings
  empty <- paste(rep(" ", .p$covmaxchar), collapse="")
  empty2 <- paste(rep(" ", .p$covmaxchar), collapse="")
  # inialized attribute parameter
  nr <- .p$covnr
  # get attribute parameter
  A <- .C("GetAttr", type=integer(nr), operator=integer(nr),
          monotone=integer(nr), finiterange=integer(nr),
          simpleArguments=integer(nr), internal=integer(nr),
          domains=integer(nr), isos=integer(nr),
          maxdim=integer(nr), vdim=integer(nr), PACKAGE="RandomFields")
  #
  
  idx <- .C("GetModelList", idx=integer(nr^2), as.integer(TRUE),
            PACKAGE="RandomFields")$idx
  dim(idx) <- c(nr, nr)

  for (i in 1:nr) {
  # sequential steps for each model
    
    if (A$internal[i]) {
        #Print("internal", .C("GetModelName",as.integer(i-1),
        #         name=empty, nick=empty2, PACKAGE="RandomFields"))
      next
    }
    # get model name
    ret <- .C("GetModelName",as.integer(i-1),
              name=empty, nick=empty2, PACKAGE="RandomFields")
    nick <- ret$nick

    ## get names of submodels
    subname.info<- .Call("GetSubNames", as.integer(i-1), PACKAGE="RandomFields")
    subnames <- subname.info[[1]]
    subintern <- subname.info[[2]]
    subnames.notintern <- subnames[!subintern]
     
    # get names of  parameters
    paramnames <- .Call("GetParameterNames", as.integer(i-1),
                        PACKAGE="RandomFields")
    elmnt <- which(paramnames == "element")
    if (length(elmnt) > 0)  {
      stopifnot(length(elmnt) == 1)
      paramnames <- c(paramnames[-elmnt], "element")
    }
    par.intern <- paramnames %in% subnames
    if (any(par.intern)) stop(nick, ": subnames (",
                              paste(subnames, collapse=", "),
                              ") and parameter names (",
                              paste(paramnames, collapse=", "),
                              ") match.")
    
    ex.anysub <- length(subnames)>0
    ex.sub <- length(subnames.notintern)>0
     
    ex.par <- length(paramnames)>0
    ex.std <- ((nick != ZF_DOLLAR[1] && isNegDef(A$type[i])) ||
               nick == "RMball"
               || nick == ZF_PLUS[1] || nick[1] == ZF_MULT[1])
    
    cat(i, "\t", nick, ",\t\told name ", ret$name, "\t", ex.std, "\t",
        A$type[i], "\n", sep="")
    
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
                                 isNegDef(A$type[i]) || A$type[i]==ShapeType)
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
         "type = ",         "'", RC_TYPE[A$type[i]+1],            "',", "\n\t",
         "domain = ",       "'", RC_DOMAIN[A$domains[i]+1],       "',", "\n\t",
         "isotropy = ",     "'", RC_ISOTROPY[A$isos[i]+1],        "',", "\n\t",
         "operator = ",     as.logical(A$operator[i]),             ",", "\n\t",
         "monotone = ",    "'", RC_MONOTONE[A$monotone[i] + MON_MISMATCH],
                                                                  "',", "\n\t",
         "finiterange = ",  as.logical(A$finiterange[i]),          ",", "\n\t",
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


