
RMcauchytbm <- function(alpha, beta, gamma, var, scale, Aniso, proj) {
  return(RMtbm(fulldim=gamma,
               RMgencauchy(alpha, beta, var, scale, Aniso, proj)))
}
  
RMcardinalsine <- function(var, scale, Aniso, proj) {
  return(RMwave(var, scale, Aniso, proj))
}
  
RMgneitingdiff <- function(nu, taper.scale, scale, var, Aniso, proj){
  return(RMmult(RMgengneiting(kappa=3, mu=1.5,  scale=taper.scale),
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

RFearth2cartesian <- function(coord, units=NULL) {
  if (is.null(units)) {
    global.units <- RFoptions()$coords$new_coord_units[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coord)) coord <- t(coord)
  return(t(RFfctn(RMtrafo(RC_CARTESIAN_COORD), coord, grid=FALSE,
                  coords.new_coord_units=units,
                  coords.coordinate_system="earth")))
}

RFearth2dist <- function(coord, units=NULL, ...) {
  if (is.null(units)) {
    global.units <- RFoptions()$coords$new_coord_units[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coord)) coord <- t(coord)
  z <- t(RFfctn(RMtrafo(RC_CARTESIAN_COORD), coord, grid=FALSE,
                coords.new_coord_units=units,
                coords.coordinate_system="earth"))
  return(dist(z, ...))
}
