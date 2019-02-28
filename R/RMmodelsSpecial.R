
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


RRdistr <- function(name, nrow, ncol,  ## ddistr, pdistr, qdistr, rdistr,
		    envir,  ...) {
  if (!missing(name)) {
    u <- try(isRMmodel(name), silent=TRUE)
    if (is.logical(u) && u) return(name)
  }

  cl <- match.call()  
  par.general <- submodels <- par.model <- list()

  if (FALSE)
  if (!missing(ddistr) && length(ddistr) != 0) {
    stopifnot(!missing(ddistr), !missing(pdistr), !missing(qdistr),
	      !missing(rdistr), !missing(envir))
    par.model <- c(list(name=name, nrow=nrow, ncol=ncol,
			ddistr=ddistr, pdistr=pdistr, qdistr=qdistr,
			rdistr=rdistr, envir=envir),
		   list(...))
    model <- new(CLASS_CLIST, call = cl, name = RM_DISTR[1],
		 submodels = submodels, 
		 par.model = par.model, par.general = par.general)
    return(model) 
  }

  ll <- as.list(substitute(name))
  name <- ll[[1]]
  if (is.character(name)) ll <- list(...)
  else {
    name <- as.character(name)
    ll <- ll[-1]
  }

  if (length(ll) > 0) {
    par.names <- names(ll)
    if (length(par.names) == 0 || any(par.names == ""))
      stop("In distribution families\n  (i) all parameters must be named;\n (ii) in more complicated models all distributions taken from R must be\n      encapsulated by 'RRdistr';\n(iii) use 'RRloc' to modifiy scale and location")
    num <- sapply(ll, function(x) {
      is.numeric(x) || is.symbol(x) || {u <- try(eval(x), silent=TRUE);
                                        !(is(u, "try-error"))}
    })
    if (!all(num)) {
      subs <- ll[!num]
      n <- names(subs)
      for (i in 1:length(subs)) {
        if (!is.language(subs[[i]]))
          stop("type of parameter (function, constant) cannot be determined")
        par.model[[n[i]]] <-
          if (substr(deparse(subs[[i]]), 1, 1)=='R') eval(subs[[i]]) else
              do.call(RM_DISTR[1], list(subs[[i]]))
      }
    }
    if (any(num)) {
      ll <- ll[num]
      n <- names(ll)
      for (i in 1:length(ll)) par.model[n[i]] <- eval(ll[[i]])
    }
  }
  
  var <- c('x', 'q', 'p', 'n')
  fctnames <-  c('d', 'p', 'q', 'r')
  pm <- list()
  pm[['name']] <- name
  if (hasArg(ncol)) pm[['ncol']] <- ncol
  if (hasArg(nrow)) pm[['nrow']] <- nrow
  for (ii in 1:length(fctnames)) {
    i <- fctnames[ii]
    iname <- paste(i, name, sep="")
    if (!exists(iname)) {
      stop("'", name,
     "' is considered as a name for a distribution family,\n  but the function '",
           iname, "' is not visible.")
    }
    pm[[paste(i, "distr", sep="")]] <-
      eval(parse(text=paste("quote(", iname, "(",
                   var[ii], "=", var[ii],
                   if (length(ll) > 0 && length(par.names)>0)
                      paste(",",
                            paste(par.names, "=", par.names, collapse=", ")),
                   "))", sep="")))
  }
  pm[['envir']] <- if (hasArg(envir)) envir else new.env()
  par.model <- c(pm, par.model)
    
  model <- new(CLASS_CLIST, call = cl, name = RM_DISTR[1], submodels = submodels, 
               par.model = par.model, par.general = par.general)
  return(model) 
}

RRdistr <- new(CLASS_RM,
               .Data = RRdistr,
               type = TYPE_NAMES[RandomType + 1],
               domain = DOMAIN_NAMES[DOMAIN_MISMATCH + 1],
               isotropy = ISO_NAMES[ISO_MISMATCH + 1],
               operator = FALSE,
               monotone =  MONOTONE_NAMES[NOT_MONOTONE + 1 - MON_UNSET],
               finiterange = FALSE,
               simpleArguments = FALSE,
               maxdim = Inf,
               vdim = 1
               )


GetSymbols <- function(ll) {
  idx <- sapply(ll, function(x) is.symbol(x) || !is.language(x))
  symbols <- as.character(ll[idx])
  if (!all(idx)) 
    for (i in which(!idx)) {
      symbols <- c(symbols, GetSymbols(as.list(ll[[i]])))
    }
  return(symbols)
}
       

RMuser <- function(type, domain, isotropy, vdim, beta,
                   coordnames = c("x", "y", "z", "T"),
                   fctn, fst, snd, envir,
                   var, scale, Aniso, proj) {
	cl <- match.call()
	submodels <- par.general <- par.model <- list() 
	
	if (!hasArg(type)) type <- TYPE_NAMES[ShapeType + 1]
        if (is.numeric(type)) par.model[['type']] <- type
        else if (is.character(type))
          par.model[['type']] <- pmatch(type, TYPE_NAMES) - 1

        
        if (any(par.model[['type']] < ProcessType))
          message("It is likely that the function you are defining is already available in 'RandomFields', or hasn't got the claimed property ('",
                  TYPE_NAMES[1+par.model[['type']]],
                  "'). (If you are not sure whether your function is positive/negative definite, please contact schlather@math.uni-mannheim.de.)\nUsing predefined functions leads to (much!) shorter computing times (up to a factor 100).\nSee ?RM for an overview over the implemented models. Further,\nsome simulation methods do not work at all for user defined functions.")
        else if (any(par.model[['type']] == TrendType))
          message("Please make sure that the defined function is not available in 'RandomFields'.\nUsing predefined functions leads to (much!) shorter computing times (up to a factor 100).\nSee ?RMmodelsTrend  for an overview over the implemented models. Further,\nsome simulation methods do not work at all for user defined functions.");

	if (hasArg(domain)) {
	  if (is.numeric(domain)) par.model[['domain']] <- domain
	  else if (is.character(domain))
                   par.model[['domain']] <- pmatch(domain, DOMAIN_NAMES) - 1
	  else stop("wrong value for 'domain'")
	}
	if (hasArg(isotropy)) {
	  if (is.numeric(isotropy)) par.model[['isotropy']] <- isotropy
	  else if (is.character(isotropy))
                   par.model[['isotropy']] <- pmatch(isotropy, ISO_NAMES) - 1
	  else stop("wrong value for 'isotropy'")
	}

	if (hasArg(vdim)) {
	  if (is.numeric(vdim)) par.model[['vdim']] <- vdim
	  else stop("wrong value for 'vdim'")
	}
	if (hasArg(beta)) {
	  if (is.numeric(beta) || is.language(beta) || 
	      is.list(beta))
	     par.model[['beta']] <- beta
	  else if (substr(deparse(substitute(beta)), 1, 1)=='R')
	    submodels[['beta']] <- beta
	  else submodels[['beta']] <- RRdistr(beta)
	}
 	if (hasArg(fctn)) {
          f <- substitute(fctn)
          par <- coordnames %in% GetSymbols(as.list(as.list(f)[-1]))
          par.model[['fctn']] <- f
 	} else stop("'fctn' not given")
       
	if (hasArg(fst)) {
          f <- substitute(fst)
          if (any(xor(par,
                      coordnames %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'fst' do not match the ones in 'fctn'")
          par.model[['fst']] <- f
	} else if (hasArg(snd)) stop("'fst' not given")
        
	if (hasArg(snd)) {         
	  f <- substitute(snd)          
          if (any(xor(par,
                      coordnames %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'snd' do not match the ones in 'fctn'")
	  par.model[['snd']] <- f
	}

        par.model[['coordnames']] <- which(par)    
	par.model[['envir']] <- if (hasArg(envir)) envir else new.env()
	par.general[['var']] <-if (hasArg(var)) var else RM_DEFAULT
	par.general[['scale']] <-if (hasArg(scale)) scale else RM_DEFAULT
	par.general[['Aniso']] <-if (hasArg(Aniso)) Aniso else RM_DEFAULT
	par.general[['proj']] <-if (hasArg(proj)) proj else RM_DEFAULT

	model <- new(CLASS_CLIST, call = cl, name = RM_USER[1], 
			submodels = submodels, 
			par.model = par.model, par.general = par.general)
	return(model) 
 } 

RMuser <- new(CLASS_RM,
              .Data = RMuser,
              type = TYPE_NAMES[PosDefType + 1],
              domain = DOMAIN_NAMES[PREVMODEL_D + 1],
              isotropy = ISO_NAMES[PREVMODEL_I + 1],
              operator = FALSE,
              monotone =  MONOTONE_NAMES[NOT_MONOTONE + 1 - MON_UNSET], # [MON_PARAMETER]
              finiterange = TRUE,
              simpleArguments = FALSE,
              maxdim = Inf,
              vdim = -1
              )



RMdeclare <- function(...) {
  PL <- RFoptions(GETOPTIONS="basic")$printlevel
     
  CL <- match.call()
  cl <- as.list(CL)
  model.name <- as.character(cl[[1]])
  cl <- cl[-1]
  symbol <- sapply(cl, is.symbol)
  par.names <- as.character(names(cl))
  par.names[symbol] <- as.character(cl[symbol])
  len <- length(par.names)
  if (len == 0 || sum(par.names != "") != len ||
      any(is.na(par.names)))
    stop("'", model.name,
         "' takes only scalar named variables whose value will be used in 'params' within an 'RF*' call.")
  if (any(sapply(cl, length) != 1) ) stop("'", model.name, "' takes only scalar variables.")
  if (len != length(unique(par.names)))
    stop("'", model.name, "' takes only distinct variable names")

  par.values <- numeric(len)
  for (i in 1:len) {    
    x <- try(base::...elt(i), silent=TRUE)
    if (is(x, "try-error") || !is.numeric(x) || length(x) != 1) {
      x <- NA
      if (PL > PL_IMPORTANT) 
        message("   '", par.names[i], "' set to ", x, ". Consider setting '", par.names[i],
                "' explicitely in the argument 'params'.\n       If this has been done already, check whether '", par.names[i], "' is a scalar.")
    } else {
      if (symbol[i] && !is.na(x) && PL > PL_IMPORTANT)
        message("   '", par.names[i], "' set to ", x, ".")
    }
    par.values[i] <- x    
  }
 
  par.model <- as.list(par.values)
  names(par.model) <- par.names 
 # Print(par.model, par.values, par.names)
  
  model <- new(CLASS_CLIST, call = CL, name = model.name, submodels = list(), 
               par.model = par.model, par.general = list())

#  print(model);
  
  return(model) 
}

RMdeclare <- new('RMmodelgenerator',
               .Data = RMdeclare,
               type = TYPE_NAMES[TrendType + 1],
               domain = DOMAIN_NAMES[XONLY+ 1],
               isotropy = ISO_NAMES[PREVMODEL_I + 1],
               operator = FALSE,
               monotone =  MONOTONE_NAMES[NOT_MONOTONE + 1],
               finiterange = FALSE,
               simpleArguments = TRUE,
               maxdim = Inf,
               vdim = -2
               )



RMcovariate <- function(formula=NULL, data, x, y=NULL, z=NULL, T=NULL, grid,
                        raw, norm, addNA, factor) {
  if (!missing(factor)) {
    if (!missing(addNA) && addNA)
      stop("'addNA' and 'factor' may not be given at the same time.")
    isna <- is.na(factor)
    if (any(xor(isna[1], isna))) stop("If 'factor' has NAs then all of the values must be NAs")
  }
  if (missing(data)) {
    if ("formula" %in% class(formula))
      stop("data argument 'data' has not been given")
    data <- formula
    formula <- NULL
  }

  N <- "data.frame"
  if (is.matrix(data)) {
    if (!is.null(colnames(data)) && !is.null(formula)) data <- as.data.frame(data)
  } else if (is.factor(data)) {
    if (!is.null(formula)) stop("if 'data' is a factor, a formula may not be given")
    data <- data.frame(own.factor = data)
    N <- "factor"
  }

  if (is.data.frame(data)) {
    if (is.null(formula)) 
      formula <- eval(parse(text=paste("~", paste(names(data), collapse="+"))))
    data <- model.matrix(object=formula, data=data)
    if (RFoptions()$basic$printlevel > 0)
      message("Intercept added to a model component that is a ", N, ".")
  }
  
  Call <- iRMcovariate
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    ans <- Call(norm=norm, data=data, raw=raw, addNA=addNA)
  } else {
    new <- C_UnifyXT(x=x, y=y, z=z, T=T, grid=grid, printlevel=0)
    ans <- Call(norm=norm, data=data, x=new, raw=raw, addNA=addNA)
  }
  ans
}
RMcovariate <- copyProp(RMcovariate, iRMcovariate)
  
 
RMfixcov <- function(M, x, y=NULL, z=NULL, T=NULL, grid, var, proj, raw, norm) {
  Call <- iRMfixcov
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    Call(norm=norm, M=M, raw=raw, var=var, proj=proj)
  } else {
    new <- C_UnifyXT(x, y, z, T, grid, printlevel=0)
    Call(norm=norm, M=M, x=new, raw=raw, var=var, proj=proj)
  }
}
RMfixcov <- copyProp(RMfixcov, iRMfixcov)


RMcov <- function(gamma, x, y=NULL, z=NULL, T=NULL, grid, a,
                  var, scale, Aniso, proj, raw, norm) {
  Call <- iRMcov
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    x <- list(1) ## origin
  } else if (is.character(x)) {
    if (!is.null(y) || !is.null(z) || !is.null(T) || !missing(grid))
      stop("If x is a character, the arguments y,z,T,grid may not be given.")
    x <- pmatch(x, RMCOV_X)
    if (is.na(x)) stop("unknown choice of 'x'")
    x <- list(as.double(x))
  } else x <- C_UnifyXT(x, y, z, T, grid, printlevel=0)
  Call(gamma=gamma, x = x, a=a, scale=scale, Aniso=Aniso, proj=proj, var=var)
}
RMcov <- copyProp(RMcov, iRMcov)


RMpolynome <- function(degree, dim, value=NA,
		       coordnames = c("x", "y", "z", "T"),
                       proj=1:4) {
  if (degree < 0 || degree > 5) stop("the degree is out of range")
  if (dim < 0  || dim > 4) stop("the dimension is out of range")
  x <- as.matrix(do.call("expand.grid", rep(list(0:degree), dim)))
  sums <- rowSums(x)
  y <- NULL
  for (i in 0:degree) {
    idx <- sums == i
    y <- rbind(y, x[idx, ])
  }
  n <- nrow(y)
#  y <- as.vector(y)
  z <- paste(rep(paste(" ", coordnames[1:dim], sep=""), each=n),
             ifelse(y>1, "^", ""),
             ifelse(y>1, y, ""), sep="")
  z[ y == 0] <- ""
  dim(z) <- dim(y)
  z <- apply(z, 1, paste, collapse="", sep="")
  m <- length(z) - length(value)
  if (m > 0) value <- c(value, rep(NA, m)) else
  if (m < 0) value <- value[1:length(z)]
  cat( paste( value, z, collapse = " + ", sep=""), "\n" )
  
  z <- paste(rep(paste("R.p(", proj[1:dim], ")", sep=""), each=n),
             ifelse(y>1, "^", ""),
             ifelse(y>1, y, ""), sep="")
  z[ y == 0 ] <- ""
  if (length(z) > 100) stop("maximum is ", MAXSUB, "^2 terms")
  dim(z) <- dim(y)
  z <- apply(z, 1, function(x)  paste(x[x!=""] , collapse="*", sep=""))
  value <- paste(R_CONST, "(", value, ")", sep="")
  z <- paste(value , z, sep="*")
  z[1] <- value[1]
  idx <- as.integer(length(z) / MAXSUB) * MAXSUB  
  if (idx > 0) {
    zz <- z[ 1:idx ]
    dim(zz) <- c(MAXSUB, length(zz) / MAXSUB)
    zz <- apply(zz, 2, function(x) paste("RMplus(", paste(x, collapse=", "), ")" ))
  } else zz <- NULL
  if (idx < length(z)) {
#    Print(zz, idx, idx + 1 == length(z))
    zz[length(zz) + 1] <-  if (idx + 1 == length(z)) z[length(z)] else 
       paste("RMplus(", paste(z[(idx+1) : length(z)], collapse=", "), ")" )
      
#    Print("A", idx, zz, (idx + 1) : length(z))
  }

  if (length(zz) > 1)
    zz <- paste("RMplus(", paste(zz, collapse=", "), ")")

 #  Print( zz)
  ##invisible
  return(eval(parse(text = zz)))
}
RMpolynome <- copyProp(RMpolynome, R.c)






xRMranef <- function(formula=NULL, data, x, y=NULL, z=NULL, T=NULL, grid,
                    var, scale, Aniso, proj, raw, norm) {
  if (hasArg("data") || is.numeric(formula) || is.data.frame(data) || 
      is(x, "RFsp") || isSpObj(x)) {
   formula <- RMcovariate(formula=formula, data=data, x=x, y=y, z=z, T=T,
                           grid=grid, raw=raw, norm=norm, addNA = TRUE)
  } else {
    if (!is(formula, "RMmodel"))
      stop("'formula' is not a 'RMmodel' as expected")
     if (!missing(data) || !missing(x) || !is.null(y) ||
         !is.null(z) || !missing(grid) || !missing(raw) || !missing(norm))
       stop("If 'formula' is an 'RMmodel' then only 'var', 'scale', 'Aniso', and 'proj' might be given") 
  }
  if (missing(var)) {
    if (RFoptions()$basic$printlevel > 0)
      message("Note that if 'var' is not given in 'RMranef', 'var' is set to 'NA' i.e., the variance is estimated'.")
    var <- NA
  }
  #RMraneffct(formula, var, scale, Aniso, proj)
}

  
XXXRMprod <- function(phi, var, scale, Aniso, proj) {
  #RMraneffct(phi, var, scale, Aniso, proj)
}
