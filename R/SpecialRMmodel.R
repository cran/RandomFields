
RMstrokorbMono <- function(phi) stop("Please use 'RMm2r' instead")
RMstrokorbBall <- function(phi) stop("Please use 'RMm3b' instead")
RMstrokorbPoly <- function(phi) stop("Please use 'RMmps' instead")

RMcoord <- function(C0, coord, dist)
{
	cl <- match.call()

	submodels <- list() 
	if(hasArg(C0)) submodels[['C0']] <- C0

	par.model <- list() 
	if(hasArg(coord)) par.model[['coord']] <- coord
	if(hasArg(dist)) par.model[['dist']] <- dist

	par.general <- list()
	par.general[['var']] <- ZF_DEFAULT_STRING
	par.general[['scale']] <- ZF_DEFAULT_STRING
	par.general[['Aniso']] <- ZF_DEFAULT_STRING
	par.general[['proj']]  <- ZF_DEFAULT_STRING

	model <- new(ZF_MODEL, call = cl, name = ZF_COORD,
	submodels = submodels,
	par.model = par.model, par.general = par.general)

	return(model) 
} 

RMcoord <- new('RMmodelgenerator',
               .Data = RMcoord,
               type = RC_TYPE[OtherType + 1],
               domain = RC_DOMAIN[PREVMODELD + 1],
               isotropy = RC_ISOTROPY[RC_CARTESIAN_COORD + 1],
               operator = TRUE,
               monotone = RC_MONOTONE[NOTMONOTONE],
               finiterange = TRUE,
               maxdim = Inf,
               vdim = -2
               )


internalRMmixed <- function(X, beta, cov, coord, dist, element)
{
	cl <- match.call()

	submodels <- list() 

	par.model <- list() 
	if(hasArg(element)) par.model[['element']] <- element
	if(hasArg(X)) par.model[['X']] <- X
	if(hasArg(beta)) par.model[['beta']] <- beta
	if(hasArg(coord)) par.model[['coord']] <- coord
	if(hasArg(dist)) par.model[['dist']] <- dist
	if(hasArg(cov)) submodels[['cov']] <- cov #par.model[['cov']] <- cov

	par.general <- list()
	par.general[['var']] <- ZF_DEFAULT_STRING
	par.general[['scale']] <- ZF_DEFAULT_STRING
	par.general[['Aniso']] <- ZF_DEFAULT_STRING
	par.general[['proj']]  <- ZF_DEFAULT_STRING

	model <- new(ZF_MODEL, call = cl, name = ZF_MIXED[1],
	submodels = submodels,
	par.model = par.model, par.general = par.general)

	return(model) 
} 

internalRMmixed <- new('RMmodelgenerator',
                       .Data = internalRMmixed,
                       type = RC_TYPE[OtherType + 1],
                       domain = RC_DOMAIN[PREVMODELD + 1],
                       isotropy = RC_ISOTROPY[RC_CARTESIAN_COORD + 1],
                       operator = TRUE,
                       monotone =  RC_MONOTONE[NOTMONOTONE],
                       finiterange = TRUE,
                       maxdim = Inf,
                       vdim = -2
                       )


RRdistr <- function(fct, nrow, ncol, envir) {
  if (!missing(fct)) {
    u <- try(isModel(fct), silent=TRUE)
    if (is.logical(u) && u)
      return(fct)
  }

  cl <- match.call()
  
  par.general <- submodels <- par.model <- list()

 # Print(substitute(fct), rate);  xxxx
  
  ll <- as.list(substitute(fct))
  name <- as.character(ll[[1]])
  ll <- ll[-1]
  if (length(ll) > 0) {
    par.names <- names(ll)
    if (length(par.names) == 0 || any(par.names == ""))
      stop("In distribution families, all parameters must be named.")
   # print(par.names)
    num <- sapply(ll, function(x) {
      is.numeric(x) || is.symbol(x) || {u <- try(eval(x), silent=TRUE);
                                        !(class(u)=="try-error")}
    })
    if (!all(num)) {
      subs <- ll[!num]
      n <- names(subs)
      for (i in 1:length(subs)) {
      #  Print(subs[[i]], n[[i]], deparse(subs[[i]]), substr(deparse(subs[[i]]), 1, 1))
   #     ffff
        if (!is.language(subs[[i]]))
          stop("type of parameter (function, constant) cannot be determined")
        par.model[[n[i]]] <-
          if (substr(deparse(subs[[i]]), 1, 1)=='R') eval(subs[[i]]) else
              do.call(ZF_DISTR[1], list(subs[[i]]))
      }
   # Print(submodels); xxx
    }
    if (any(num)) {
      ll <- ll[num]
      n <- names(ll)
      for (i in 1:length(ll)) par.model[n[i]] <- eval(ll[[i]])
    }
    #Print(submodels, par.model, ll, num, eval(par.model[[1]]))
  }
  
  var <- c('x', 'q', 'p', 'n')
  fctnames <-  c('d', 'p', 'q', 'r')
  for (ii in 1:length(fctnames)) {
    i <- fctnames[ii]
    par.model[[paste(i, "distr", sep="")]] <-
      eval(parse(text=paste("quote(", i, name, "(",
                   var[ii], "=", var[ii],
                   if (length(ll) > 0 && length(par.names)>0)
                      paste(",",
                            paste(par.names, "=", par.names, collapse=", ")),
                   "))", sep="")))
  }
  
  if (hasArg(ncol)) par.model[['ncol']] <- ncol
  if (hasArg(nrow)) par.model[['nrow']] <- nrow
  par.model[['envir']] <-  if (hasArg(envir)) envir else new.env()
    
  model <- new(ZF_MODEL, call = cl, name = ZF_DISTR[1], submodels = submodels, 
               par.model = par.model, par.general = par.general)
  return(model) 
}

RRdistr <- new('RMmodelgenerator',
               .Data = RRdistr,
               type = RC_TYPE[.RandomType + 1],
               domain = RC_DOMAIN[PREVMODELD + 1],
               isotropy = RC_ISOTROPY[PREVMODELI + 1],
               operator = FALSE,
               monotone =  RC_MONOTONE[NOTMONOTONE],
               finiterange = FALSE,
               maxdim = Inf,
               vdim = 1
               )




GetSymbols <- function(ll) {
  idx <- sapply(ll, function(x) is.symbol(x) || !is.language(x))
  symbols <- as.character(ll[idx])
  #Print(ll, idx)
  if (!all(idx)) 
    for (i in which(!idx)) {
      symbols <- c(symbols, GetSymbols(as.list(ll[[i]])))
    }
  return(symbols)
}
       

RMuser <- function(type, domain, isotropy, vdim, beta,
                   variab.names = c("x", "y", "z", "T"),
                   fctn, fst, snd, envir, var, scale, Aniso, proj) {
	cl <- match.call()
	submodels <- par.general <- par.model <- list() 
	
	if (!hasArg(type)) type <- RC_TYPE[ShapeType + 1]
        if (is.numeric(type)) par.model[['type']] <- type
        else if (is.character(type))
          par.model[['type']] <- pmatch(type, RC_TYPE) - 1
        if (type < ProcessType)
          message("It is likely that the defined function is already available in 'RandomFields'. Using predefined functions leads to (much) shorter computing times. See ?RMmodels for an overview over the implemented models. Further, some simulation methods do not work at all for user defined functions.")
        else if (type == TrendType)
          message("Please make sure that the defined function is not available in 'RandomFields'. Using predefined functions leads to (much) shorter computing times. Further, some simulation methods do not work at all for user defined functions.");

	if (hasArg(domain)) {
	  if (is.numeric(domain)) par.model[['domain']] <- domain
	  else if (is.character(domain))
                   par.model[['domain']] <- pmatch(domain, RC_DOMAIN) - 1
	  else stop("wrong value for 'domain'")
	}
	if (hasArg(isotropy)) {
	  if (is.numeric(isotropy)) par.model[['isotropy']] <- isotropy
	  else if (is.character(isotropy))
                   par.model[['isotropy']] <- pmatch(isotropy, RC_ISOTROPY) - 1
	  else stop("wrong value for 'isotropy'")
	}

        
   #     Print(par.model, type, domain, isotropy, RC_TYPE, RC_DOMAIN, RC_ISOTROPY); xxx        
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
          par <- variab.names %in% GetSymbols(as.list(as.list(f)[-1]))
          par.model[['fctn']] <- f
 	} else stop("'fctn' not given")
       
	if (hasArg(fst)) {
          f <- substitute(fst)
          if (any(xor(par,
                      variab.names %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'fst' do not match the ones in 'fctn'")
          par.model[['fst']] <- f
	} else if (hasArg(snd)) stop("'fst' not given")
        
	if (hasArg(snd)) {         
	  f <- substitute(snd)          
          if (any(xor(par,
                      variab.names %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'snd' do not match the ones in 'fctn'")
	  par.model[['snd']] <- f
	}

        ##      Print(par.names, par);        xxx
        par.model[['variab.names']] <- which(par)    
	par.model[['envir']] <- if (hasArg(envir)) envir else new.env()
	par.general[['var']] <-if (hasArg(var)) var else ZF_DEFAULT_STRING
	par.general[['scale']] <-if (hasArg(scale)) scale else ZF_DEFAULT_STRING
	par.general[['Aniso']] <-if (hasArg(Aniso)) Aniso else ZF_DEFAULT_STRING
	par.general[['proj']] <-if (hasArg(proj)) proj else ZF_DEFAULT_STRING

	model <- new(ZF_MODEL, call = cl, name = ZF_USER[1], 
			submodels = submodels, 
			par.model = par.model, par.general = par.general)
	return(model) 
 } 

RMuser <- new('RMmodelgenerator',
              .Data = RMuser,
              type = RC_TYPE[PosDefType + 1],
              domain = RC_DOMAIN[PREVMODELD + 1],
              isotropy = RC_ISOTROPY[PREVMODELI + 1],
              operator = FALSE,
              monotone =  RC_MONOTONE[NOTMONOTONE], # [MON_PARAMETER]
              finiterange = TRUE,
              maxdim = Inf,
              vdim = -1
              )


f <- substitute(exp(-x-y+z(), zz=a))
