if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

DeleteAllRegisters()
RFparameters(Print=5)
y <- x <- seq(0,10, 0.05)
model <- list("*",
              list("$", var=1, aniso=matrix(nr=2, c(1,0)),
                   list("fractgauss", k=1)),
              list("$",var=1, aniso=matrix(nc=2, c(0,0,0,1)),
                   list("whittle", nu=0.6))
                   )


z <- GaussRF(x, y, model=model, grid=TRUE)
str(GetRegisterInfo(0))
print(range(z))
image(x, y, z)


PrintModelList();
