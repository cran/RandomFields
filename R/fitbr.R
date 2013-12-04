### here function to fit Brown-Resnick processes

fit.br <- function(...) {
  stop("estimation of Brown-Resnick fields not programmed yet")
  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
}
