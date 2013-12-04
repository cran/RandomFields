### here general function to fit maxstable processes

fit.maxstable.scale <- function(...) {
  stop("de Haan estimator not programmed yet") ## joined paper with estimating trend
  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

}


fit.extremal.gauss <- function(...) {
 stop("estimation of Schlather model not programmed yet") ## joined paper with estimating trend
  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]


}


fit.smith <- function(...) {
  stop("estimation of Smith model not programmed yet") ## joined paper with estimating trend
  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

}

