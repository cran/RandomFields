
library(RandomFields, lib=if (file.exists("/home/schlather/TMP/RandomFields"))
                       "~/TMP")

q()

.path <- "/home/schlather/article/R/NEW.RF/RandomFields/R/"
if (file.exists(paste(.path, "rf.R", sep=""))) {
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
} else q()
