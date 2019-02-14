# accessing 'RMmodel' and RMmodelgenerator via '['-operator
# e.g. RMwhittle["domain"]

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





RMcauchytbm <- function(alpha, beta, gamma, var, scale, Aniso, proj) {
  return(RMtbm(fulldim=gamma,
               RMgencauchy(alpha, beta, var, scale, Aniso, proj)))
}
RMcauchytbm <- new(CLASS_RM, 
	.Data = RMcauchytbm,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)
  

RMgneitingdiff <- function(nu, taper.scale, scale, var, Aniso, proj){
  return(RMmult(RMgengneiting(kappa=3, mu=1.5, scale=taper.scale) *
                RMwhittle(nu=nu, scale=scale),
                var=var, Aniso=Aniso, proj=proj))
}
RMgneitingdiff <- new(CLASS_RM, 
	.Data = RMgneitingdiff,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)

RMparswmX <- function(nudiag, rho, var, scale, Aniso, proj) {
  return(RMschur(M=rho, RMparswm(nudiag, var, scale, Aniso, proj)))
}
RMparswmX <- copyProp(RMparswmX, RMparswm)

RMpoweredexp <- function(alpha, var, scale, Aniso, proj) {
  return(RMstable(alpha, var, scale, Aniso, proj))
}
RMpoweredexp <- copyProp(RMpoweredexp, RMstable)

RMtent <- function(var, scale, Aniso, proj) {
  return(RMaskey(alpha=1.0, var, scale, Aniso, proj))
}
RMgneitingdiff <- new(CLASS_RM, 
	.Data = RMgneitingdiff,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)

R.lon <- copyProp(function() R.p(proj=1, new="spherical system"), R.p)
R.lat <- copyProp(function() R.p(proj=2, new="spherical system"), R.p)

RMhandcock <- function(nu, notinvnu, var, scale, Aniso, proj){
  RMS(scale = 1/sqrt(2), RMmatern(nu, notinvnu, var, scale, Aniso, proj))
}
RMhandcock <- copyProp(RMhandcock, RMmatern)


RMcardinalsine <- RMwave
RMwendland <- RMgengneiting

RMchoquet <- function(b) stop("not implemented yet")
