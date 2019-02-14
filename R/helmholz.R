

xRMderivative <- function(phi, partial, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('partial') && !is.null(subst <- substitute(partial))) 
	par.model[['partial']] <- CheckArg(partial, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new(CLASS_CLIST, call = cl, name = 'RMderivative', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

xRMderivative <- new(CLASS_RM,
	.Data = xRMderivative,
	type = c('positive definite'),
	isotropy = c('zero-space-isotropic'),
	domain = 'single variable',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



xRMhelmholtz <- function(phi, component, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('component') && !is.null(subst <- substitute(component))) 
	par.model[['component']] <- CheckArg(component, subst, TRUE)
  if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.model[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new(CLASS_CLIST, call = cl, name = 'RMhelmholtz', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

xRMhelmholtz <- new(CLASS_RM,
	.Data = xRMhelmholtz,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = 'single variable',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)
