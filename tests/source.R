
library(RandomFields,
        lib=if(TRUE && file.exists("/home/schlather/TMP/RandomFields")) "~/TMP")

.path <- "/home/schlather/R/RF/RandomFields/R/"
if (EXTENDED.TESTING <- file.exists(paste(.path, "rf.R", sep=""))) {
  EXTENDED.TESTING <- FALSE
  Source <- function(x) {
    cat(x, "...\n")
    x <- paste(.path, x, sep="")
    source(x)
  }
  Source("D.H.R")
  Source("MLE.R")
  Source("ShowModels.R")
  Source("auxiliary.R")
  Source("convert.R")
  Source("empvario.R")
  Source("evalpar.R")
  Source("extremes.R")
  Source("getNset.R")
  Source("modelling.R")
  Source("rf.R")
} 
