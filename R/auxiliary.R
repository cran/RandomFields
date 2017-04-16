
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


.RandomFields.env <- new.env()



vectordist <- function(x, diag=FALSE) {
  storage.mode(x) <- "double"
  res <- .Call(C_vectordist, t(x), diag)
  dimnames(res) <- list(dimnames(x)[[2]], NULL)
  return(t(res));
}



my.legend <- function(lu.x, lu.y, zlim, col, cex=1, ...) {
  ## uses already the legend code of R-1.3.0
  cn <- length(col)
  if (cn < 43) {
    col <- rep(col, each=ceiling(43 / cn))
    cn <- length(col)
  }
  filler <- vector("character", length=(cn-3)/2)
  legend(lu.x, lu.y, y.intersp=0.03, x.intersp=0.1, 
         legend=c(format(zlim[2], dig=2), filler,
             format(mean(zlim), dig=2), filler,
             format(zlim[1], dig=2)),
         lty=1, col=rev(col),cex=cex, ...)
}
