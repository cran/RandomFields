

CheckMixed <- function(arg, subst, names) {
  ## never allows for distributions
  if (is.character(arg)) {
    if (length(arg) != 1) stop("'", deparse(substitute(arg)),
			       "' must be a single string")      
    arg <- -pmatch(toupper(arg), toupper(names))
    if (is.na(arg)) stop("'", deparse(substitute(arg)),
			 "': unknown character sequence")
    return(arg)
  }
  u <- try(is.numeric(arg) || is.logical(arg) || is.language(arg)
           || is.list(arg) || is(arg, class2=CLASS_CLIST), silent=TRUE)
  if (is.logical(u) && u) arg
  else stop("'",  deparse(substitute(arg)),
	    "': submodels and random parameter not allowed")
}



CheckProj <- function(arg, subst) {
  if (is.character(arg)) {
    if (length(arg) != 1) stop("'proj' must be a single string")
    ## anyway most of the following part is written as if p was vector
    p <- -pmatch(arg, PROJECTION_NAMES)
    if (any(is.na(p))) {
      p <- pmatch(arg, RFoptions()$coord$coordnames)
      if (length(p) != length(arg))
	stop("projection definition '", arg,  "' is odd")
      if (any(is.na(p))) { 
	p <- pmatch(arg, c("x", "y", "z", "t", "X", "Y", "Z", "T"))
	if (!any(is.na(p))) { 
	  p <- (p - 1) %% 4 + 1
	  p[p == 4] <-
	  pmatch("time", PROJECTION_NAMES) -1 -length(PROJECTION_NAMES)
	} else {
	  p <- arg ## must be a single character string
	  while (!(substr(p, 1, 1) %in% 0:9)) p <- substring(p, 2)
	  warn <- options()$warn
	  options(warn = 0)            
	  p <- try(as.vector(p), silent=TRUE)
	  options(warn=warn)
	  if (is(p, "try-error"))
	    stop("'\"", arg,
		 "\"' is interpretated as a projection defintion but the character string is not recognized")
	}
        }
    }
    return(p)
  }
  u <- try(is.numeric(arg) || is.logical(arg) || is.language(arg)
	   || is.list(arg) || is(arg, class2='RMmodel'), silent=TRUE)
  if (is.logical(u) && u) return(arg)
  else if (substr(deparse(subst), 1, 1)=='R') arg
  else  do.call('RRdistr', list(subst))
}



CheckMaths <- function(arg, subst, distr) {
  u <- try(is.numeric(arg) || is.logical(arg) || is.language(arg)
           || is.list(arg) || is(arg, class2='RMmodel'), silent=TRUE)
  if (is.logical(u) && u) arg
  else if (is.character(arg)) do.call('R.p', list(arg))
  else if (substr(deparse(subst), 1, 1)=='R') arg
  else if (distr) do.call('RRdistr', list(subst))
  else stop('random parameter not allowed')
}


CheckArg <- function(arg, subst, distr) {
  u <- try(is.numeric(arg) || is.logical(arg) || is.language(arg)
           || is.list(arg) || is(arg, class2=CLASS_CLIST), silent=TRUE)
  
  ##Print(subst);  Print(distr);  Print(u)
  if (is.logical(u) && u) arg
  else if (substr(deparse(subst), 1, 1)=='R') arg
  else if (distr) do.call('RRdistr', list(subst))
  else stop('random parameter not allowed')
}

CheckChar <- function(arg, subst, names, distr) {
  if (is.character(arg)) {
    if (length(arg) != 1) stop("'", deparse(substitute(arg)),
			       "' must be a single string")   
    arg <- pmatch(arg, names) - 1
    if (any(is.na(arg))) stop("'", deparse(substitute(arg)),
			      "': unknown string value")
    return(arg)
  }
  u <- try(is.numeric(arg) || is.logical(arg) || is.language(arg)
           || is.list(arg) || is(arg, class2=CLASS_CLIST), silent=TRUE)
  if (is.logical(u) && u) arg
  else if (substr(deparse(subst), 1, 1)=='R') arg
  else if (distr) do.call('RRdistr', list(subst))
  else stop('random parameter not allowed')
}


copyProp <- function(what, from) {
  return(new(CLASS_RM,
             .Data = what,
             type = from["type"],
             isotropy = from["isotropy"],
             domain = from["domain"],
             operator = from["operator"],
             monotone = from["monotone"],
             finiterange = from["finiterange"],
             simpleArguments = from["simpleArguments"],
             maxdim = from["maxdim"],
             vdim = from["vdim"]
             ))
}
