


.RandomFields.env <- new.env()



vectordist <- function(x, diag=FALSE) {
  storage.mode(x) <- "double"
  res <- .Call("vectordist", t(x), diag)
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
