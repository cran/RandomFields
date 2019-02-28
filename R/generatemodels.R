
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  


## meta-function that generates functions like RMplus, RMwhittle, which
## the user will use to generate explicit covariance models, i.e. objects
## of class 'RMmodels'

param.text.fct <- function(catheg, nick, names, havedistr=TRUE, Const=NULL,
                           ismath=FALSE){
  ifHasArg <- paste("  if (hasArg('", names,
                    "') && !is.null(subst <- substitute(", names,
                    "))) \n", sep="")
  if (ismath && any(idx <- names == "new" & Const == NN2)) {
    ifHasArg[idx] <- "  if (!(hasArg('new') && !is.null(subst <- substitute(new)))) new <- UNREDUCED\n"
  }
  
  x <- paste(ifHasArg, "\t", catheg, "[['", names, "']] <- ", sep="")
  for (i in 1:length(names)) {
    if (#!ismath && ## deleted Feb 2017
       names[i] == "proj")
      x[i] <- paste(x[i], "CheckMixed(proj, subst, PROJECTION_NAMES)", sep="")
## Version brasilien    
##       x[i] <- paste(x[i], "CheckProj(proj, subst)", sep="")
   else if (length(Const) > 0 && Const[i] == MixedInputType) 
       x[i] <- paste(x[i], "CheckMixed(", names[i], ", subst, ",
                     toupper(nick), "_",  toupper(names[i]),
                     # ", ", havedistr,
                     ")", sep="")
    else if (length(Const) > 0 && Const[i] == CharInputType) 
       x[i] <- paste(x[i], "CheckChar(", names[i], ", subst, ",
                     toupper(nick), "_",  toupper(names[i]),
                     ## ", ",  havedistr,
                     ", FALSE",
                     ")", sep="")
    else if (length(Const) > 0 && Const[i] >= NN1)
      x[i] <- paste(x[i], "CheckChar(", names[i], ", subst, ",
                    NAMES_OF_NAMES[Const[i] - NN1 + 1], ", ",
                    havedistr, ")", sep="")
   else if (ismath)
       x[i] <- paste(x[i], "CheckMaths(", names[i], ", subst, ",
                    havedistr, ")", sep="")    
    else x[i] <- paste(x[i], "CheckArg(", names[i], ", subst, ",
                       havedistr, ")", sep="")
  }
  x
}


rfGenerateModels <- function(package="RandomFields", assigning,
                             RFpath = "~/svn/RandomFields/RandomFields",
                             RMmodels.file = paste(RFpath,
						   "R/RMmodels.R",
						   sep="/"),
			     PL = RFoptions()$basic$printlevel
			     ) {

  # if file already exists, remove it.
  if (assigning && file.exists(RMmodels.file))
    file.remove(RMmodels.file)
  
   write(file = RMmodels.file, append = TRUE,
        "\n## This file has been created automatically by 'rfGenerateModels'.\n\n")

  ## defined constants
  diminf <- 999999
  
  # define empty strings
  empty <- paste(rep(" ", MAXCHAR_RF), collapse="")
  empty2 <- paste(rep(" ", MAXCHAR_RF), collapse="")
  # inialized attribute parameter
  
  nr <- GetCurrentNrOfModels()
  vn <- nr * MAXVARIANTS

  # get attribute parameter
  A <- .C(C_GetAttr, nr=integer(vn), type=integer(vn), operator=integer(vn),
          monotone=integer(vn), finiterange=integer(vn),
          simpleArguments=integer(vn), internal=integer(vn),
          domains=integer(vn), isos=integer(vn),
          maxdim=integer(vn), vdim=integer(vn),
          includevariants= as.integer(TRUE),
          paramtype = integer(vn * MAXPARAM),
          va = integer(1) )
  va <- A$va
  dim(A$paramtype) <- c(MAXPARAM, vn)

  i <- 1
  while (i <= va) {
    step <- 1
    ## sequential steps for each model

    if (A$internal[i] && A$internal[i] != INTERN_SHOW) {
       if (PL > 1)
	cat(i, "internal", .C(C_GetModelName,as.integer(A$nr[i]),
			      name=empty, nick=empty2)$name,"\n")      
      i <- i + 1
      next
    }

    domains <-  A$domains[i]
    if (domains == PREVMODEL_D) domains <- c(XONLY, KERNEL)
    
    # get model name
    ret <- .C(C_GetModelName,as.integer(A$nr[i]),
              name=empty, nick=empty2)
    internal.nick <- nick <- ret$nick
    if (A$internal[i]) internal.nick <- paste("i", internal.nick, sep="")
    ismath <- substr(nick, 1, 2) == "R."
    
    #if (PL > 1)cat(i, nick, "\n");

    type <- A$type[i]
    iso <- A$isos[i]
    while (i + step <= va  &&
           nick ==  .C(C_GetModelName,as.integer(A$nr[i + step]),
               name=empty, nick=empty2)$nick) {
      if (PL > 1)cat("...variant added\n")
      type <- c(type, A$type[i + step])
      iso <- c(iso, A$isos[i + step])
      step <- step + 1
    }

   
    finiterange <- as.logical(A$finiterange)
    finiterange[A$finiterange < 0] <- NA
 
    ## get names of submodels
    subname.info<- .Call(C_GetSubNames, as.integer(A$nr[i]))
    subnames <- subname.info[[1]]
    subintern <- subname.info[[2]]
    subnames.notintern <- subnames[!subintern]
     
    # get names of  parameters
    paramnames <- .Call(C_GetParameterNames, as.integer(A$nr[i]) )
    internal <- which(paramnames == INTERNAL_PARAM)
    if (length(internal) > 0) paramnames <- paramnames[-internal]
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
    ex.std <- ((nick != DOLLAR[2] && any(isVariogram(type))) ||
               nick %in% c("RMball", "RMsum", "RMconstant",
                           "RMfixcov", "RMcovariate")
               || nick == RM_PLUS[1] || nick[1] == RM_MULT[1]) &&
      !(nick %in% c("RMtrafo",  "RMsine")) && !ismath
    
    

    std.variables <-
      if (nick %in% c("RMid")) NULL
      else if (nick %in% c("RMfixcov")) c("var", "proj")
      else if (nick %in%  c("RMconstant", "RMcovariate")) "var"
      else if (nick == "RMnugget") c("var", "Aniso", "proj")
      else c("var", "scale", "Aniso", "proj")

   # cat(nick, ":", std.variables, "\n")
    
  #   if (PL > 1) cat( std.variables, ex.std)
   #stopifnot(i < 50)
    
    
    if (PL > 1)cat(i, "\t", internal.nick, ",\t",
        paste(std.variables, collapse=", "), "\t",
        ex.std, "\t",
        paste(DOMAIN_NAMES[domains+1], collapse="; "), "\t",
        paste(type, collapse="/"), "\n", sep="")
    
    if(nick == DOLLAR[2]){ 
      text.fct.head <-
        paste(nick, " <- function(phi, var, scale, Aniso, proj, anisoT)")

      ##Print(type, isVariogram(type), ShapeType, paramnames, A$paramtype[,i])
      
    } else {
      text.fct.head <-
        paste(internal.nick,
              " <- function(",
              if (ex.sub) {
                paste(paste(subnames.notintern, collapse=", "), sep="")
              },
              if (ex.sub && (ex.par || ex.std)) ", ",
              if (ex.par) {
                paste(paste(paramnames, collapse=", "), sep="")
              },
              if (ex.par && ex.std) ", ",
              if (ex.std) paste(std.variables, collapse =", "), 
              ")",
              sep="")
    }

     
    if (ex.par) {
      par.body <- param.text.fct(catheg="par.model", nick=nick,
                                 names=paramnames,
                                 havedistr=any(isVariogram(type))
                                 || any(type==ShapeType)
#                                 || any(type==RandomOrShapeType)
                                ,
                                 Const=A$paramtype[1:length(paramnames), i],
                                 ismath=ismath)
      if (any(idx <- paramnames == 'envir'))
        warning(internal.nick, ": envir not internal")
#      if (any(idx))
#        par.body[idx] <-
 #         "par.model[['envir']] <- if (hasArg(envir)) envir else new.env()"
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
            "\n",
            ## get model specific parameter
            if (ex.par) paste(par.body, collapse="\n"),
            if (ex.par) "\n  ",
            ## get general model parameter
            if (ex.std) {
              paste(param.text.fct(catheg="par.general", nick=nick,
                                   names=std.variables),
                    collapse="\n  ")
            },
            "\n  ",
             # create RMmodel object
            "model <- methods::new('", CLASS_CLIST, "', ",
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
    
    # assign class CLASS_RM (CLASS_RM) and attributes like stationarity
    # to the function:

     text.assign.class <-
      paste(internal.nick, " <- new(CLASS_RM, \n\t",
         ".Data = ",        internal.nick,                                  ",", "\n\t",
         "type = ", "c('", paste(TYPE_NAMES[type+1], collapse="', '"), "'),",
            "\n\t",
         "isotropy = ", "c('", paste(ISO_NAMES[iso+1], collapse="', '"), "'),",
            "\n\t",
         "domain = ", "c('", paste(DOMAIN_NAMES[domains+1], collapse="', '"),   "'),", "\n\t",
         "operator = ",     as.logical(A$operator[i]),             ",", "\n\t",
         "monotone = ",    "'", MONOTONE_NAMES[A$monotone[i] + 1 - MON_UNSET],
                                                                  "',", "\n\t",
         "finiterange = ",  finiterange[i],          ",", "\n\t",
         "simpleArguments = ",  as.logical(A$simpleArguments[i]),  ",", "\n\t",
         "maxdim = ", if(A$maxdim[i]>diminf) Inf else A$maxdim[i], ",", "\n\t",
         "vdim = ",         A$vdim[i],                                  "\n\t",
         ")",
         sep="")
 
    text <- paste(text.fct, "\n", text.assign.class, "\n\n\n", sep="")

    if (internal.nick == "RMwhittle") {
      if (PL > 1)cat(text.assign.class)
     # stop("KKKK")
    }
  
    if (assigning) {
      #sink(file = RMmodels.file, append = TRUE, type='output')
      write(file = RMmodels.file, append = TRUE, text)
      #if (PL > 1) cat(text)
      #sink()
      #unlink(RMmodels.file)
    }
    i <- i + step    
  }  ## matches for (i in 1:nr)


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

clean_value <- function(s) {
  s <- clean_name(s)
  if (s == "NA_INTEGER") s <- "as.integer(NA)"
  return(s)
}

kind <- function(Zeilen, i, start, cont="", ignore=" ",
		 endofname=" ", stops=NULL) {
  if (substr(Zeilen[i], 1, nchar(start)) == start) {
   
    s <- substring(Zeilen[i], nchar(start) + 1)
    j <- 2
    while (j <= nchar(s) && substr(s, j, j) != endofname) j <- j + 1
    stopifnot(j <= nchar(s))
    name <- substr(s, 1, j -1)
    u <- strsplit(substring(s, j + 1), "//")[[1]]
    kommentar <- paste(u[-1], collapse = " ")
    RC <- nchar(kommentar) > 0 &&
      nchar(strsplit(kommentar, "RC")[[1]][1]) < nchar(kommentar)
    value <- u <- clean_value(u[1])
    i <- i + 1
    RCX <- NULL
    if (any(cont != "")) {
       repeat {
 	if (i > length(Zeilen) ||
	    (length(stops) > 0 && length(strsplit(u, stops)[[1]]) > 1)) break;
	## cat(i, " >>", Zeilen[i], "<<\n", sep="")
        j <- nchar(u)
	if (ignore != "") while (substr(u, j, j) %in% ignore) {
	  j <- j - 1
	}
	##	print(u);	print(j) ;print(cont)
	##print(c(substr(u, j,j), cont, !(substr(u, j,j) %in% cont) , j))
        if (!(substr(u, j,j) %in% cont) && j > 0) break
	while (Zeilen[i] == "") i <- i + 1
	u <- strsplit(Zeilen[i], "//")[[1]]
	kommentar <- paste(u[-1], collapse = " ")
	u <- u[1]
	RCX <- c(RCX, nchar(kommentar) > 0 &&
		 nchar(strsplit(kommentar, "RC")[[1]][1]) < nchar(kommentar))
 	value <- paste(value, clean_value(u[1]))
 	## cat("u=", u, "\n")
	i <- i + 1
      }
    }
    res <- list(name=name, value=value,
		RC=if (length(RCX) == 0 || all(!RCX)) RC else RCX,
		i=i)
   ## cat("value=", value, "\n\n")
    return(res)
  } else return(NULL)
}

clean_name <- function(x) {
##  str(x)
  coll <- paste(strsplit(x, "\t")[[1]], collapse="")
  i <- 1
  while(i <= nchar(coll)) {
##    Print(i, coll, substr(coll, i, i))
    if (substr(coll, i, i) %in% c("\"", "'")) {
      repeat {
	i <- i + 1
	if (substr(coll, i, i) %in% c("\"", "'")) {
	  i <- i + 1
	  break;
	}
      }
    }
    if (substr(coll, i, i) == " ") {
      coll <- paste0(substr(coll, 1, i-1), substr(coll, i+1, nchar(coll)))
    } else i <- i + 1
    ## x <- paste(strsplit(coll, " ")[[1]],  collapse=""); Print(x); xx
   }
##  Print(coll)
  coll
}


CC <- function(x, envir) {
  if (is.numeric(x)) {
    Real <- TRUE
    Integer <- x == as.integer(x)
    y <- x
  } else {
    stopifnot(is.character(x))
    if (substr(x, 1, 1) =='"')  #'"')
	return(x)
    warn <- options()$warn
    options(warn = -1)
    y <- try(as.numeric(eval(parse(text=x), envir=envir)), silent=TRUE)
    options(warn = warn)
    if (Real <- !is(y, "try-error") && !is.na(y)) {
      Integer <- nchar(strsplit(x, "\\.")[[1]][1]) == nchar(x) &&
                 abs(y) <= .Machine$integer.max
##      cat(Integer, "", y,"", .Machine$integer.max, "",abs(y) <= .Machine$integer.max)
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
  function(package="RandomFields", aux.package = "RandomFieldsUtils",
	   RFpath = paste0("~/svn/",package, "/", package),
           RCauto.file = paste(RFpath, "R/aaa_auto.R", sep="/"),
	   header.source =
	   c(if (length(aux.package) > 0) paste0("../../", aux.package,"/",
				 aux.package, "/src/Auto", aux.package, ".h"),
	     paste0("src/Auto",package,".h")),
	   c.source = paste0("src/Auto", package, ".cc")) {
	   
    write(file = RCauto.file,
	  "# This file has been created automatically by 'rfGenerateConstants'")
    envir <- new.env()
    for (s in header.source) {
      write(file = RCauto.file, append = TRUE, paste("\n\n ## from ", s))
      rfGenerateConstantsHeader(RFpath = RFpath, RCauto.file=RCauto.file,
				header.source = s, envir=envir) 
    }

    cat("** header files done\n");
    
    for (s in c.source) {
      write(file = RCauto.file, append = TRUE, paste("\n\n ## from ", s))
      rfGenerateConstantsC(package = package,
	RFpath = RFpath, RCauto.file=RCauto.file,
	c.source = s)
      
    }
  }
    
rfGenerateConstantsHeader <- function(RFpath, RCauto.file, header.source,
				      envir) {
  s <- scan(paste(RFpath, header.source, sep="/"),
            what=character(), sep="\n", blank.lines.skip=FALSE, skip=2)
  #if (PL > 1)for (i in 1:length(s)) cat(s[i], "\n")
  i <- 1
  typedef <- character(0) ## nur fuer alte defs a la "typedef int dom_type;"
  while (i <= length(s)) {
    if (!is.null(k <- kind(s, i, start="#define", cont="\\", ignore=""))) {
     ## if (i>130 && i<140)
      value <- clean_value(k$value)
      if (length(typedef) > 0)
        for (j in 1:length(typedef)) {
          value <- paste(strsplit(value, typedef[j])[[1]], collapse=" ")
        }
      v <- clean_name(k$name)
      ##      Print(i,  k, v);     stopifnot(i < 70)
      
      if (v != "") {
	w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
	line <- paste(w, v , "\t<-", CC(value, envir=envir))
	eval(parse(text=line), envir=envir)
	write(file = RCauto.file, append = TRUE, line)
      }
    } else if (!is.null(k <- kind(s, i, start ="typedef enum",
				  cont=c(",", "{"), ignore=" ", stops=";"))) {
					#
      ##Print(i,  k, v);      stopifnot(i < 72)

      value <- strsplit(clean_value(k$value), "\\{")[[1]]
      if (value[1] != "")
	typedef <- c(typedef, paste("\\(", strsplit(value[1], ";")[[1]][1] ,
				    "\\)", sep=""))
 
      value <- strsplit(strsplit(value[2], "\\}")[[1]][1], ",")[[1]]
      
      zaehler <- 0
      for (j in 1:length(value)) {
        v <- clean_value(value[j])
##	Print(v, value)
	if (v != "") {
	  w <- if ((length(k$RC) == 1 && k$RC) || (length(k$RC) > 1) && k$RC[j])
		 paste("RC_", v, " <-", sep="") else ""
	  v <- strsplit(v, "=")[[1]]
	  if (length(v) > 1) {
	    stopifnot(length(v) == 2)
	    nr <- v[2]
	    v <- v[1]
	  } else {
	    nr <- zaehler
	    zaehler <- zaehler + 1
	  }
	  line <- paste(w, v , "\t<-", CC(nr, envir=envir))
	  eval(parse(text=line), envir=envir)
	  write(file = RCauto.file, append = TRUE, line)
	}
      }
      write(file = RCauto.file, append = TRUE, "")
    } else if (!is.null(k <- kind(s, i, start="typedef", cont="\\",
				  ignore="", stops=";"))) {
     if (value[1] != "")
       typedef <- c(typedef,
		    paste("\\(", strsplit(k$value, ";")[[1]][1] ,"\\)", sep=""))
    } else if (!is.null(k <- kind(s, i, start="extern const", cont=",",
				  ignore=" ", stops=";"))) {
      ## ignored
    } else write(file = RCauto.file, append = TRUE, "")
    i <- if (is.null(k)) i+1 else k$i
  }
}
  
rfGenerateConstantsC <- function(RFpath, RCauto.file, c.source, package="none") {
  s <- scan(paste(RFpath, c.source, sep="/"),
            what=character(), sep="\n", blank.lines.skip=FALSE, skip=1)
  i <- 1
  nl <- TRUE
  while (i <= length(s)) {
     if (!is.null(k <- kind(s, i, "  *", c(",", "{"), " ", endofname="[",
			    stops="}"))) {
       ##      Print(k); stopifnot(k$i != 55)
       value <- strsplit((clean_value(k$value)), "\\{")[[1]]
       value <- strsplit(value[2], "\\}")[[1]][1]
       v <- clean_name(k$name)
       if (v != "") {
	 w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
	 line <- paste(w, v, "<-\nc(", value, ")")
	 write(file = RCauto.file, append = TRUE, line)
	 write(file = RCauto.file, append = TRUE, "")
	 nl <- TRUE
       }
     } else if (nl) {
       write(file = RCauto.file, append = TRUE, "")
       nl <- FALSE
     }
     i <- if (is.null(k)) i+1 else k$i
  }

 
  if (package == "RandomFields") {

    define_char <- function(name, value) {
      if (is(value, "try-error")) value <- NULL
      write(file = RCauto.file, append = TRUE,
	    paste("\n", name, " <- c('", sep="",
		  paste(value, collapse="', '"),
		  "')")
	    )
    }
    define_num <- function(name, value) {
      if (is(value, "try-error")) value <- NULL
      write(file = RCauto.file, append = TRUE,
	    paste("\n", name, " <- c(", sep="",
		  paste(value, collapse=", "),
		  ")")
	    )
    }
    
    envir <- as.environment("package:RandomFields")
    all <- ls(envir=envir)
    genuine <- all[substr(all, 1, 2) %in% c("iR", "RM", "R.", "RP" ,"RF", "RR")]
 
    define_char("list2RMmodel_Names", genuine) ##, RM_TREND
    define_char("list2RMmodel_oldNames",
		try(RFgetModelNames(newnames=FALSE))) # , RM_INTERNALMIXED

 
    names1 <- c("RMwhittle",
                RFgetModelNames(type=TYPE_NAMES[c(TcfType, PosDefType) + 1],
                                isotropy=ISO_NAMES[ISOTROPIC + 1],
                                operator=FALSE,
                                group.by=NULL,
                                valid.in.dim = 1,#if (sim_only1dim)1 else 2,
                                simpleArguments = TRUE,
                                vdim=1))

    do.not.include <-
      c("RMnugget", # macht kaum sinn
        "RMwendland","RMcardinalsine","RMpoweredexp", # aliase
        "RMparswmX", "RMtent", # convenience models
        "RMconstant", ## macht aerger
        "RMlsfbm", ## nur fuer |x| < 1 definiert
        "RMdagum", ## internal parameter
        "RMgneiting" ## integer parameter
        )
    names1 <- sort(names1[!(names1  %in% do.not.include)])
    define_char("rfgui1_Names", names1)
    

    names2 <- c("RMwhittle",
                RFgetModelNames(type=TYPE_NAMES[c(TcfType, PosDefType) + 1],
                                isotropy=ISO_NAMES[ISOTROPIC + 1],
                                operator=FALSE,
                                group.by=NULL,
                                valid.in.dim = 2,#if (sim_only1dim)1 else 2,
                                simpleArguments = TRUE,
                                vdim=1))
    names2 <- sort(names2[!(names2  %in% do.not.include)])
     define_char("rfgui2_Names", names2)              
    }
 
  return(NULL)
}


rfGenerateTest <- function(package = "RandomFields",
			   files = NULL,
			   RFpath = paste0("~/svn/", package, "/", package)) {
   start.after <-  "\\dontrun"
   end.before <- c("\\", "}")
   initial.text <- "if (RFoptions()$internal$do_tests){"
   final.text <- "}"
   comment <- "%"

   if (length(files) == 0) return()
   ncomment <- nchar(comment)[1]
   nendbefore <- nchar(end.before)[1]
   for (f in 1:length(files)) {
    cat("creating ", files[f], ".R\n", sep="")
    s <- scan(paste(RFpath, "/man/", files[f], ".Rd", sep=""),
              what=character(), sep="\n", blank.lines.skip=FALSE, skip=2)
    i <- 1
    while (i <= length(s) && substr(s[i], 1,
               nchar(start.after)) != start.after) i <- i + 1 
    i <- i + 1
    if (i <= length(s))
      for (j in i:length(s)) {
        if (substr(s[j], 1, ncomment) %in% comment) s[j] <- ""
        if (substr(s[j], 1, nendbefore) %in% end.before) {
          j <- j -1
          break;
        }
      }
    if (j >= i) {
      out <- paste(RFpath, "/tests/", files[f], ".R", sep="")
      write(file = out, append = FALSE, initial.text)
      write(file = out, append = TRUE, s[i:j])
      write(file = out, append = TRUE, final.text)
    }
  }

  return(NULL)
}


.R.fmax <- function(...) {
  a <- list(...)
  if (length(a) == 1) if (is(a[[1]], CLASS_CLIST)) a[[1]] else R.c(a[[1]])
  else R.fmax(if (is(a[[1]], CLASS_CLIST)) a[[1]] else R.c(a[[1]]),
	      do.call(".R.fmax", a[-1]))
}


.R.fmin <- function(...) {
  a <- list(...)
  if (length(a) == 1) if (is(a[[1]], CLASS_CLIST)) a[[1]] else R.c(a[[1]])
  else R.fmin(if (is(a[[1]], CLASS_CLIST)) a[[1]] else R.c(a[[1]]),
	      do.call(".R.fmin", a[-1]))
}



rfGenerateMaths <- function(package = "RandomFields",
			    files = "/usr/include/tgmath.h",
			    do.cfile = FALSE,
                            ## copy also in ../private/lit
                            Cfile = "QMath", Rfile = "RQmodels",
                            RFpath = paste0("~/svn/",package,"/", package)) {
  
  prefix <- "R."
  start.after <-  "/* Unary functions"
  end.before <- c("#define carg")
  initial.text <- "if (RFoptions()$internal$do_tests){"
  final.text <- "}"
  comment <- c("/*", "#i", "#e")
  nendbefore <- nchar(end.before)[1]

  if (FALSE) {
  
  stopifnot(length(files) > 0)
  if (do.cfile) {
     cfile <-  paste(RFpath, "/src/", Cfile, ".cc", sep="")
    write(file = cfile,
	  c("// This file has been created automatically by 'rfGenerateMaths'",
	  "#include <math.h>",
	  "#include \"RF.h\"",
	  "#include \"primitive.h\""
	  ))
  }
  
  manfile <- paste(RFpath, "/man/", Cfile, ".Rd", sep="")
  write(file = manfile, 
        scan(paste(manfile, 0, sep="."),
             what=character(), sep="\n", blank.lines.skip=FALSE, skip=0))  

  Rfile <- paste(RFpath, "/R/", Rfile, ".R", sep="")
  write(file = Rfile, 
          "# This file has been created automatically by 'rfGenerateMaths'")
  
  usage <- include <- list()
  for (f in 1:length(files)) {
    cat("scanning ", f, files[f],"\n")
    s <- scan(files[f], what=character(), sep="\n", blank.lines.skip=FALSE,
              skip=0)
    i <- 1
    while (i <= length(s) && substr(s[i], 1, nchar(start.after)) !=start.after){
      i <- i + 1
    }
    i <- i + 1
 
    usage[[f]] <- include[[f]] <- rep("", length(s))
    while (i <= length(s)) {

      if (substr(s[i], 1, nendbefore) %in% end.before) break
      if (!is.null(k <- kind(s, i, start="#define", cont="\\", ignore=" ",
			     endofname=")"))) {
        x <- strsplit(k$name, "\\(")
        name <- clean_name(x[[1]][1])
	if (name == "") { i <- i + 1; next }
	if (!(name %in% c("frexp", "ldexp", "remquo", "scalbn", "scalbln",
                          #"ilogb",
			  "fma"))) { 
          args <- length(strsplit(x[[1]][2], ",")[[1]])
                                        #cat(name, args, k$name, "\n")
	  if (do.cfile) {
	    write(file = cfile, append = TRUE,
		  paste("void Math", name,
			"(double *x, cov_model *cov, double *v){",
			"\nMATH_DEFAULT\n", sep=""))
	    write(file = cfile, append = TRUE,
		  if (name %in% c("cos", "sin", "tan")) {
		    paste("*v = ", name,
                        "(GLOBAL.coords.anglemode == radians",
                        " ? w[0] : w[0] * piD180);", sep="")
		  } else if (name %in% c("acos", "asin", "atan", "atan2")) {
		    paste("*v = ", name, "(w[0]",
			  if (args == 2) ", w[1]", ");\n",
			  "if (GLOBAL.coords.anglemode == degree) *v/=piD180);",
			  sep="")
		  } else {
		    paste("*v = ", name, "(w[0]", if (args == 2) ", w[1]", ");",
			  sep="")
		  })
	    write(file = cfile, append = TRUE, "\n}\n\n")
          
 	    include[[f]][i] <-
	    paste(
              'IncludeModel(".', name, '", MathDefType, 0, 0, ', args,
			    ', NULL, XONLY,\n\t',## auch fuer kernel schreiben??
			    'PREVMODEL_I,checkMath,rangeMath, PREF_TREND,\n\t',
			    'false, SCALAR, PREVMODEL_DEP, false, false); \n',
	      'nickname("', name, '");\n',
	      'kappanames("a", REALSXP',
			  if (args == 2) ', "b", REALSXP', ');\n',
	      'addCov(Math', name, ', NULL, NULL);\n',
	      'AddVariant(TrendType, PREVMODEL_I);\n', sep="")
	  }
	    
          write(file = manfile, append = TRUE,
                paste("\\alias{", prefix, name, "}", sep=""))
          usage[[f]][i] <- paste(prefix, name, "(a", if (args == 2) ", b",
                                 ")", sep="")

          if (name %in% c(
"asin",
"atan",
"atan2",
"cos",
"sin",
"tan",
"acosh",
"asinh",
"atanh",
"cosh",
"sinh",
"tanh",
"exp",
"log",
"expm1",
"log1p",
#"logb",
"exp2",
"log2",
"pow",
"sqrt",
"hypot",
"cbrt",
"ceil",
"fabs", # abs
"floor",
"fmod",
"round",
"trunc",
"erf",
"erfc",
"gamma", # gamma
"lgamma",
#"nearbyint",
#"lrint",
#"llrint",
#"lround",
#"llround",
#"copysign",
#"rint",
#"nextafter",
#"nexttoward",
#"remainder",
#"fdim",
"fmax",
"fmin")) {
	    Rname <- switch(name,
			    "fabs" = "abs",
			    "pow"  = "^",
			    "fmin" = "min",
			    "fmax" = "max",
			    "fmod" = "%%",
			    "ceil" = "ceiling",
			    name)
##	    if (name == "fabs") "abs" else
####           if (name == "tgamma") "gamma" else
##           if (name == "pow") "^" else
##           if (name == "fmin") "min" else
##           if (name == "fmax") "max" else name
           
	    write(file = manfile, append = TRUE,
		  paste("\\alias{", if (Rname=="%%") "\\%\\%" else Rname, "}",
			sep=""))
	    args <- switch(Rname,
	            "atan2" = "y, x",
                    "min"   = "...",
		    "max"   = "...",
		    "round" = "x, digits=0",
		    "^"   = "x, y",
		    "%%"  = "x, y",
	#	    "logb" = "x, base=exp(1)",
                    "x")
	    innerargs <- switch(Rname,
				"round" = "x",
				args)
	    CLASS_CLIST_ANY <- c("c(CLASS_CLIST,'ANY')", "c('ANY',CLASS_CLIST)")
	    signature <- switch(Rname,
				"round" = "c(CLASS_CLIST, 'missing')",
				# "logb" = "c(CLASS_CLIST, 'missing')",
				"^"=CLASS_CLIST_ANY,
				"%%"=CLASS_CLIST_ANY,
				"atan2"=CLASS_CLIST_ANY,
				"CLASS_CLIST")

	    if (args != "...")
	      usage[[f]][i] <- paste(Rname, "(", args, ")\n",
				     usage[[f]][i], sep="")
	   
	    write(file = Rfile, append=TRUE,
		  if (exists(Rname, envir=parent.env(parent.env(.GlobalEnv)))) {
		    if (args != "...")
		      paste(sep="", "setMethod(\"", Rname,
			    "\", signature = ", signature,
			    ", definition=function(",
			    args, ") ", prefix, name, 
			    "(", innerargs, "))")
		    else {
		      if (FALSE)
	              paste(Rname,
			    " <- function(", args, ") if (any(sapply(list(",
			    args, "), function(x) is(x, \"", CLASS_CLIST,
			    "\")))) .", prefix, name, 
			    "(", args, ") else ", "base::", Rname,"(",
			    args, ")", sep="")
		    }
		 } else paste(sep="", name, " <- ", prefix, name)
		 )
	  }
        } # !name in
      } # !is.null k
      i <- if (is.null(k)) i+1 else k$i
    } ## while i <= length
  } # for 1:files
    

  if (do.cfile)
    write(file = cfile, append=TRUE, "void includeStandardMath() {")
  write(file = manfile, append=TRUE,
        scan(paste(manfile, 1, sep="."),
             what=character(), sep="\n", blank.lines.skip=FALSE, skip=0))  
  for (f in 1:length(files)) {
    if (do.cfile)
      write(file = cfile, append=TRUE, include[[f]][include[[f]] != ""]);
    write(file = manfile, append=TRUE, usage[[f]][usage[[f]] != ""]);
  }
  if (do.cfile)
    write(file = cfile, append=TRUE, "}")
  write(file = manfile, append=TRUE,
        scan(paste(manfile, 2, sep="."),
             what=character(), sep="\n", blank.lines.skip=FALSE, skip=0))

} ## NO QMath.* generated anymore
    
  return(NULL)
 }


