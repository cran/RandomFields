
# .x <- dir(); for (.i in .x) source(.i)

library(RandomFields,
        lib=if (TRUE && file.exists("~/R/OLD") &&
          file.exists("/home/schlather/R/OLD/RandomFields")) "~/R/OLD")

.path <- "~/R/old.randomfields/RandomFields/R/"
.path2 <- "/home/schlather/R/old.randomfields/RandomFields/R/"
if (EXTENDED.TESTING <- file.exists(paste(.path, "rf.R", sep="")) &&
    file.exists(paste(.path2, "rf.R", sep=""))) {
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

# .f <- dir(); for (.i in .f) {cat("\n\n\n\n\n\n\n", .i);source(.i)}
