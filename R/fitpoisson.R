### here function to fit Poisson fields

fit.poisson <- function(...) {
  stop("estimation of Poisson fields not programmed yet")

  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  
}
