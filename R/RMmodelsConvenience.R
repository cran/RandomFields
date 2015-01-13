
RMcauchytbm <- function(alpha, beta, gamma, var, scale, Aniso, proj) {
  return(RMtbm(fulldim=gamma,
               RMgencauchy(alpha, beta, var, scale, Aniso, proj)))
}
  
RMcardinalsine <- function(var, scale, Aniso, proj) {
  return(RMwave(var, scale, Aniso, proj))
}
  
RMgneitingdiff <- function(nu, taper.scale, scale, var, Aniso, proj){
  return(RMmult(RMgengneiting(kappa=3, mu=1.5, scale=taper.scale),
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
