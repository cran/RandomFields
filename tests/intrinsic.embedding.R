# source("intrinsic.embedding.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")
PrintModelList()
print(GetMethodNames())


RFparameters(Storing=TRUE, Print=8, CE.force=TRUE, TBMCE.force=TRUE)
meth <- "intr"

######################################################################
## 2 D
######################################################################

## definition of the realisation
x <- seq(0, 2, len=12)
y <- seq(0, 1, len=15)
grid.size <- c(length(x), length(y))
model <- list(list(model="fractalB", var=1.1, kappa=1, aniso=c(1, 0, 0, 1)))# num [1:2, 1:2]  2.07e-05  1.66e+00 -5.52e-01  6.90e-06
model <- list(list(model="fractalB", var=1.1, kappa=1, aniso=c(1, 1, 0, 3)))# num [1:2, 1:2]  2.07e-05  1.66e+00 -5.52e-01  6.90e-06

# determination of the (minimal) size of the torus
DeleteRegister()
RFparameters(Storing=TRUE, Print=6, CE.force=TRUE, TBMCE.force=TRUE)
InitGaussRF(x, y, model=model, grid=TRUE, method=meth)
ce.info <- GetRegisterInfo()$method[[1]]$mem$new$method[[1]]$mem
blocks <- ceiling(ce.info$size / grid.size)
size <- blocks * grid.size
RFparameters(Print=1)
str(GetRegisterInfo())

DeleteRegister()
RFparameters(Print=1)
n <- if (interactive()) 200000 else 2;
rf <- GaussRF(x, y, model=model, grid=TRUE, method=meth, n=n,
             local.dependent=TRUE, pch="")
0.5 * var(rf[1,length(y),] - rf[length(x), 1 , ])
rm("rf")
Variogram(matrix(c(2, -1), ncol=2), model=model) ## true value

# simulation and plot of the torus
DeleteRegister()
rf <- GaussRF(x, y, model=model, grid=TRUE, method=meth, n=prod(blocks),
             local.dependent=TRUE, local.mmin=size)
hei <- 8
do.call(getOption("device"),
        list(height=hei, wid=hei / blocks[2] / diff(range(y)) *
             blocks[1] * diff(range(x))))

close.screen(close.screen())
split.screen(rev(blocks[1:2]))
k <- 0
time <- 1
for (j in 1:blocks[2]) {
  for (i in 1:blocks[1]) {
    k <- k + 1
    screen(k)
    par(mar=rep(1, 4) * 0.02)
    image(rf[,,(blocks[2]-j) * blocks[1]  + i], zlim=c(-3, 3), axes=FALSE)
  }
}


######################################################################
## 3 D
######################################################################

# definition of the realisation
x <- seq(0, 2, len=12)
y <- seq(0, 1, len=15)
z <- seq(0, 0.5, len=10)

grid.size <- c(length(x), length(y), length(z))
model <- list(list(model="fractalB", var=1.1, kappa=1,
                   aniso=c(2,1,0.5, -1,1,2, -1,3,-1)))

# determination of the (minimal) size of the torus
DeleteRegister()
#RFparameters(local.mmin=c(60, 72, 80))
InitGaussRF(x, y, z=z, model=model, grid=TRUE, method=meth)
ce.info <- GetRegisterInfo()$method[[1]]$mem$new$method[[1]]$mem
blocks <- ceiling(ce.info$size / grid.size)
size <- blocks * grid.size
# str(GetRegisterInfo(0, TRUE))

DeleteRegister()
RFparameters(Print=1)
n <- if (interactive()) 10000 else 2;
rf <- GaussRF(x, y, z, model=model, grid=TRUE, method=meth, n=n,
             local.dependent=TRUE, pch="")
cat(0.5 * var(rf[1,length(y), 1, ] - rf[length(x), 1 , length(z), ]),
    Variogram(matrix(c(2, -1, 0.5), ncol=3), model=model), "\n") ## true value
rm("rf")

