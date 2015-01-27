########################################################################
## Hier sind Funktionen, die vermutlich kaum ein Anwender interessiert
## und wohl eher zum internen Gebrauch, wenn auch im NAMESPACE stehend;
##
## Ausnahme: rfGenerate* stehen in generatemodels.R
########################################################################

PrintModelList <-function (operators=FALSE, internal=FALSE,
                           newstyle=TRUE) {
   stopifnot(internal <=2, internal >=0, operators<=1, operators>=0)
    .C("PrintModelList", as.integer(internal), as.integer(operators),
       as.integer(newstyle),
       PACKAGE="RandomFields")
    invisible()
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


checkExamples <- function(exclude=NULL, include=1:length(.fct.list),
                          ask=FALSE, echo=TRUE, halt=FALSE, ignore.all=FALSE,
                          path="RandomFields", package="RandomFields",
                          read.rd.files=TRUE) {
  .exclude <- exclude
  .ask <- ask
  .echo <- echo
  .halt <- halt
  .ignore.all <- ignore.all
  .path <- path
  .package <- package
  useDynLib <- importClassesFrom <- import <-
    importFrom <- exportClasses <-
      importMethodsFrom <- exportMethods <- S3method <- function(...) NULL
  .env <- new.env()

  ##  if (!require(.package, character.only=TRUE)) return(NA)
  exportPattern <- function(p) {
    all.pattern <- p %in% c("^[^\\.]", "^[^.]", ".") | get("all.pattern", .env)
    if (!.ignore.all) assign("all.pattern", all.pattern, .env)
    if (all.pattern) return(NULL)
    stopifnot(nchar(p)==2, substr(p,1,1)=="^")
    assign("p", c(get("p", .env), substring(p, 2)), .env)
  }
  export <- function(...) {
    ## code from 'rm'
    dots <- match.call(expand.dots = FALSE)$...
    z <-deparse(substitute(...))
    if (length(dots) && !all(sapply(dots, function(x) is.symbol(x) || 
                                    is.character(x)))) 
      stop("... must contain names or character strings")
    z <- sapply(dots, as.character)
    assign("export", c(get("export", .env), z), .env)
  }
  assign("export", NULL, .env)
  assign("all.pattern", FALSE, .env)
  assign("p", NULL, .env)
  source(paste(.path, "NAMESPACE", sep="/"), local=TRUE)
  if (is.logical(read.rd.files) && !read.rd.files) {
    .package.env <- parent.env(.GlobalEnv)
    while (attr(.package.env, "name") != paste("package:", .package, sep="")) {
      .package.env <- parent.env(.package.env)
    }
    .orig.fct.list <- ls(envir=.package.env)
    .ok <- (get("all.pattern", .env) |
            substr(.orig.fct.list, 1, 1) %in% get("p", .env) | 
            .orig.fct.list %in% get("export", .env))
    .fct.list <- .orig.fct.list[.ok]
  } else {
    .path <- read.rd.files
    if (is.logical(.path))
      .path <- "/home/schlather/svn/RandomFields/RandomFields/man"
    .files <- dir(.path, pattern="d$")
    .fct.list <- character(length(.files))
    for (i in 1:length(.files)) {
      cat(i, .files[i], "\n")
      #if (i == 152) {cat("jumped\n"); next}
      .content <- scan(paste(.path, .files[i], sep="/") , what=character(),
                       quiet=TRUE)
      .content <- strsplit(.content, "alias\\{")
      .content <- .content[which(lapply(.content, length) > 1)][[1]][2]
      .fct.list[i] <-
        strsplit(strsplit(.content,"\\}")[[1]][1], ",")[[1]][1]
    }
  }
  Print(.fct.list, length(.fct.list)) #
  .include <- include
  RFoptions(graphics.always_close_screen = TRUE) 
  .RFopt <- RFoptions()
  .not_working_no <- .not_working <- NULL
  .included.fl <- .fct.list[.include]
  .not.isna <- !is.na(.included.fl)
  .include <- .include[.not.isna]
  .included.fl <- .included.fl[.not.isna]
  Print(.included.fl, get("export", .env), get("p", .env)); #return(NA)
  .max.fct.list <- max(.included.fl)
  for (.idx in .include) {
    if (.idx %in% .exclude) next
    cat("\n\n\n\n\n", .idx, " ", .package, ":", .fct.list[.idx],
        " (total=", length(.fct.list), ") \n", sep="")
    RFoptions(LIST=.RFopt)
    .C("ResetWarnings", as.integer(FALSE))
    if (.echo) cat(.idx, "")
    .tryok <- TRUE
    if (.halt) {
      do.call(utils::example, list(.fct.list[.idx], ask=.ask, echo=.echo))
    } else {
       if (is(try(do.call(utils::example, list(.fct.list[.idx], ask=.ask,
                                         echo=.echo))), "try-error")) {
        .not_working_no<- c(.not_working_no, .idx)
        .not_working <- c(.not_working, .fct.list[.idx])
        .tryok <- FALSE
      }
    }
    RFoptions(storing = FALSE);
    cat("****** '", .fct.list[.idx], "' (", .idx, ") done. ******\n")
    if (.tryok && !is.na(RFoptions()$general$seed)) {
      Print(.not_working, paste(.not_working_no, collapse=", ")) #
      stop("seed not NA: ", .fct.list[.idx])
    }
  }
  Print(.not_working, paste(.not_working_no, collapse=", ")) #
  .ret <- list(.not_working, .not_working_no)
  names(.ret) <- c(.package, "")
  return(.ret)
}


FinalizeExample <- function() {
  if (!interactive()) {
    close.screen(all.screens = TRUE)    
    n <- length(dev.list()) - 2 ## otherwise R CMD check complains
    ##                             for missing graphic device
    if (n > 0) {
      for (i in 1:n) {      
        dev.off() ## OK
      }
    }
  }
  RFoptions(seed = NA)
}



showManpages <- function(path="/home/schlather/svn/RandomFields/RandomFields/man") {
  files <- dir(path)
  result <- character(length(files))
  for (i in 1:length(files)) {
    cat("\n\n\n\n\n")
    system(paste("grep \"\\examples\" -A 10000 ", path, "/", files[i], sep=""))
    result[i] <- readline(files[i])
  }
  result <- cbind(files, result)
  result <- result[result[, 2]!= "", ]
  result[result[, 2]=="t", 2] <- "To Do"
  result[result[, 2]=="s", 2] <- "dontshow"
  return(result)
}
# showManpages()

compareVersions <- function(path = "~/svn/RandomFields/cran") {
  l <- list(list(subpath="RandomFields/R", pattern="*R"),
            list(subpath="RandomFields/src", pattern="*\\.h"),
            list(subpath = "RandomFields/src", pattern = "*\\.cc"))
  
  for (i in 1:length(l)) {
    subpath <- l[[i]]$subpath
    pattern <- l[[i]]$pattern
    files <- dir(path=subpath, pattern=pattern)
    
    for (f in files) {
      cmd <- paste("diff ",  path, "/", subpath, "/", f, " ",
                   subpath,"/", f, sep="")
      Print(cmd)  #
      system(cmd)
    }
  }
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


#restore_par <- function(oldpar) {
#  do.call(graphics::par, oldpar)
#  graphics::par(cex = oldpar$cex) ## muss extra aufgerufen werden. noch mehr?
#}

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




reverse_dependencies_with_maintainers <-
  function(packages, which = c("Depends", "Imports", "LinkingTo"),
           recursive = FALSE) {
    ## function taken from CRAN developer website. 
    repos <- getOption("repos")["CRAN"]
    if (substr(repos, 1, 1) == "@") repos <- "http://cran.r-project.org"
    Print(repos) #
    contrib.url(repos, "source") # trigger chooseCRANmirror() if required
    description <- sprintf("%s/web/packages/packages.rds", repos)
    con <- if(substring(description, 1L, 7L) == "file://")
      file(description, "rb")
    else
      url(description, "rb")
    on.exit(close(con))
    db <- readRDS(gzcon(con))
    rownames(db) <- NULL
    
    rdepends <- tools::package_dependencies(packages, db, which,
                                            recursive = recursive,
                                            reverse = TRUE)
    rdepends <- sort(unique(unlist(rdepends)))
    pos <- match(rdepends, db[, "Package"], nomatch = 0L)
    
    db[pos, c("Package", "Version", "Maintainer")]
  }

Dependencies <- function(install = all(pkgs == all.pkgs),
                         check=TRUE, pkgs = all.pkgs, dir = "Dependencies") {
  all <- reverse_dependencies_with_maintainers("RandomFields", which="all")
  all.pkgs <- all[, 1]
  PKGS <- paste(all[,1], "_", all[,2], ".tar.gz", sep="")    
  
  ## getOption("repos")["CRAN"]
  URL <- "http://cran.r-project.org/src/contrib/"
  if (install) {
    system(paste("mkdir ", dir))
    system(paste("rm ", dir, "/*tar.gz*", sep=""))
    for (i in 1:length(pkgs)) {
      cat("PACKAGE:", PKGS[i], ":", i, "out of ", length(pkgs),"\n")
      x <- system(paste("(cd ", dir, "; wget ", URL, PKGS[i], ")", sep=""))
      if (x != 0) stop(PKGS[i], "not downloadable")
    ## extended version see RandomFields V 3.0.51 or earlier     
    }
  }
 
  if (check) {
    tools::check_packages_in_dir(dir=dir)
    return(NULL)

    ## old:
    for (i in 1:length(pkgs)) {
      command <- paste("(cd ", dir, "; R CMD check --as-cran", PKGS[i],")")
      Print(command) #
      x <- system(command)
      if (x != 0) stop(PKGS[i], "failed")
    }
  }
}
# R Under development (unstable) (2014-12-09 r67142) -- "Unsuffered Consequences"

#Dependencies()
