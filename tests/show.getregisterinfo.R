
if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")
PrintModelList()
print(GetMethodNames())

x <- seq(0, 1, 0.1)
y <- seq(2, 4, 0.2)

RFparameters(Storing=TRUE, Print=1, CE.force=TRUE, TBMCE.force=TRUE)
grid <- TRUE
nug <- 0.5

for (method in GetMethodNames()) {
  cat("\n\n", method)
  if (interactive()) readline("-- press return")
  try(GaussRF(x, y, model="exp", param=c(0.1, 1.1, nug, 1.7 - 0.7),
              grid=grid, method=method))
  str(GetRegisterInfo(0, TRUE))
}


for (method in c("coins")) {
  if (interactive()) readline(paste("\n\n", method, "-- press return"))
  GaussRF(x, y, model="circular", param=c(0.1, 1.1, nug, 1.7),
          grid=grid, method=method)
  str(GetRegisterInfo(0, TRUE))
}


