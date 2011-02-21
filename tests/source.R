
# .x <- dir(); for (.i in .x) source(.i)

library(RandomFields,
        lib=if (TRUE && file.exists("~/TMP/RandomFields") &&
          file.exists("/home/schlather/TMP/RandomFields")) "~/TMP")

.path <- "~/R/RF/RandomFields/R/"
.path2 <- "/home/schlather/R/RF/RandomFields/R/"
if (EXTENDED.TESTING <- file.exists(paste(.path, "rf.R", sep="")) &&
    file.exists(paste(.path2, "rf.R", sep=""))) {
#  EXTENDED.TESTING <- FALSE
  Source <- function(x) {
    cat(x, "...\n")
    x <- paste(.path, x, sep="")
    source(x)
  }
  Source("auxiliary.R")
  Source("convert.R")
  Source("D.H.R")
  Source("empvario.R")
#  Source("evalpar.R")
  Source("extremes.R")
  Source("getNset.R")
  Source("MLE.R")
  Source("modelling.R")
  Source("rf.R")
  Source("ShowModels.R")
} 

# .f <- dir(); for (.i in .f) {cat("\n\n\n\n\n\n\n", .i);source(.i)}
