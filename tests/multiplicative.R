source("source.R")

DeleteAllRegisters()
RFparameters(Print=5)
y <- x <- seq(0,10, 0.05)
model <- list(list(model="fractgauss", kappa=1, var=1, aniso=c(1,0,0,0)),
              "*",
              list(model="whittle", kappa=0.6, var=1, aniso=c(0,0,0,1)))
z <- GaussRF(x, y, model=model, grid=TRUE)
str(GetRegisterInfo(0))
print(range(z))
image(x, y, z)


PrintModelList();
