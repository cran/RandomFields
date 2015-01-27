
RMcauchytbm <- function(alpha, beta, gamma, var, scale, Aniso, proj) {
  return(RMtbm(fulldim=gamma,
               RMgencauchy(alpha, beta, var, scale, Aniso, proj)))
}
  
RMcardinalsine <- function(var, scale, Aniso, proj) {
  return(RMwave(var, scale, Aniso, proj))
}
  
RMgneitingdiff <- function(nu, taper.scale, scale, var, Aniso, proj){
  return(RMmult(RMgengneiting(kappa=3, mu=1.5, scale=taper.scale) *
                RMwhittle(nu=nu, scale=scale),
                var=var, Aniso=Aniso, proj=proj))
}

RMparswmX <- function(nudiag, rho, var, scale, Aniso, proj) {
  return(RMschur(M=rho, RMparswm(nudiag, var, scale, Aniso, proj)))
}

RMpoweredexp <- function(alpha, var, scale, Aniso, proj) {
  return(RMstable(alpha, var, scale, Aniso, proj))
}

RMtent <- function(var, scale, Aniso, proj) {
  return(RMaskey(alpha=1.0, var, scale, Aniso, proj))
}

RMwendland <- function(kappa, mu, var, scale, Aniso, proj) {
  return(RMgengneiting(kappa, mu, var, scale, Aniso, proj))
}

RMcovariate <- function(c, x, y=NULL, z=NULL, T=NULL, grid, factor, var) {
  if (missing(x) ) {
    if (!is.null(y) || !is.null(z) || !is.null(T) || !missing(grid))
      stop("y, z, T, grid  may only be given if 'x' is given")
    if (missing(var)) RMcovariateIntern(c=c, factor=factor)
    else RMS(var=var,  RMcovariateIntern(c=c, factor=factor))
  } else {
    x <- CheckXT(x, y, z, T, grid, printlevel=0)
      if (missing(var)) RMcovariateIntern(c=c, factor=factor, x=x$x,
                                          if (!is.null(T)) T=x$T, grid=x$grid)
      else RMS(var=var, RMcovariateIntern(c=c, factor=factor, x=x$x,
                   if (!is.null(T)) T=x$T, grid=x$grid))
  }
}
