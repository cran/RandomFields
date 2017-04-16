## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2017 -- 2017 Martin Schlather
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


CHOICE <- c("low", "medium", "high")
TAILCHOICE <- c("compact", "exponential", "power")
UNIVARIATE <- 1

## make sure that the results can be reobtained in the future
## -> use explicite the method for simiulation, not RPgauss

## get peoples suggestion on nu(x) and also other
## non-stationary cov fcts

## can covariates be involved?

## non-Gaussian marginal distribution: empirically based; Box-Cox based?!

## Gneiting, whittle, Cauchy, ex(?!) multivariate verallg. von cauchy?

xRFdatasets <- function(nonstationarity=CHOICE,
                       trendnonstationarity=CHOICE,
                       anisotropy=CHOICE,
                       differentiability=CHOICE,
                       tail = TAILCHOICE,
                       grid= c(FALSE, TRUE),
                       
                       locations = CHOICE, ## or number
                       spacedimension=1:3,
                       multivariate = UNIVARIATE,
                       time=c(FALSE, TRUE),

                       holdout_points,

                       covariates,

                       trend,
                       
                       n = 1,                       
                       seed = 0,
                       marginal= "Gaussian"
                      ) {
  res <- array(dim=c(length(nonstationarity),
                 length(trendnonstationarity),
                 length(anisotropy),
                 length(differentiability),
                 length(tail),
                 length(grid),
                 n
                 ))



 # Bsplines

}
                       
