#source("RFtest.model.R")

if (file.exists("source.R")) source("source.R")

x <- function(...) {
  str(PrepareModel(...))
  cat("--------------------------------\n")
  str(m <- convert.to.readable(PrepareModel(...)))
  invisible(m)
}
#x <- function(...) str(PrepareModel(...))

x(model=list(list(model="whi",kappa=5,var=2,s=4),"+",
    list(model="whi",kappa=1,var=6,s=1),"*",
    list(model="nugg",var=3,s=0)
   ),
  ti=1, me="ci")


# stop("")
x(model="nugg",param=c(0,2,3,4),ti=1)
x(model="nugg",param=c(0,0,3,4),ti=1)
x(model="nugg",param=c(0,3,0,4),ti=1)
x(model="nugg",param=c(0,3,0,4),ti=1,me="ci")

x(model="sphe",param=c(1,2,0,4),ti=1)
x(model="sphe",param=c(1,2,3,4),ti=1)
x(model="sphe",param=c(1,2,3,4),ti=1, me="ci")

x(model="sphe",param=cbind(c(2,4),c(8,1)),ti=1, me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0)),ti=1, me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0)),ti=1, trend=1, me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0)),ti=1, trend=c(1,2),me="ci")#not simple
x(model="sphe",param=cbind(c(2,4),c(3,0),c(7,0)),ti=1, me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0),c(8,1)),ti=1, me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0),c(8,1)),ti=1, trend=1)

x(model="whi",param=c(1,2,3,4,5),ti=1)
x(model="whi",param=c(1,2,0,4,5),ti=1)
x(model="whi",param=c(1,2,3,4,5),ti=1, me="ci")
x(model="whi",param=c(1,2,0,4,5),ti=1, me="ci")

x(model="whi",param=cbind(c(2,4,5),c(8,1,1)),ti=1, me="ci")
x(model="whi",param=cbind(c(2,4,5),c(3,0,0)),ti=1, me="ci", trend=9)
x(model="whi",param=cbind(c(2,4,5),c(3,0,0),c(7,0,0)),ti=1, me="ci")
x(model="whi",param=cbind(c(2,4,5),c(3,0,0),c(8,1,1)),ti=1, me="ci")

# list to old, vector
x(model=list(list(model="whi",kappa=5,var=2,s=4)),ti=1, me="ci")
x(model=list(list(model="whi",kappa=5,var=2,s=4),"+",
    list(model="whi",kappa=1,var=3,s=0)
    ),
  ti=1, me="ci")
## note difference to the next result !!!!!!!!!
xx <- convert.to.readable(PrepareModel(model=
                                   list(list(model="whi",kappa=5,var=2,s=4),
                                        "+",
                                        list(model="whi",kappa=1,var=3,s=0)),
                                   ti=1, me="ci"))
x(model=xx$mo, pa=xx$pa, ti=1, me=xx$me)


# list to nested
x(model=list(list(model="whi",kappa=5,var=2,s=4),"+",
    list(model="nugg",var=3,s=0),"+",
    list(model="nugg",var=7,s=0)
    ),
  ti=1, me="ci")

x(model=list(list(model="whi",kappa=5,var=2,s=4),"+",
    list(model="whi",kappa=1,var=6,s=1),"+",
    list(model="nugg",var=3,s=0)
   ),
  ti=1, me="ci")

x(model=list(list(model="whi",kappa=5,var=2,s=4),"+",
    list(model="whi",kappa=1,var=6,s=1),"+",
    list(model="nugg",var=3)
   ),
  ti=1, me="ci")


# multiplicative model
x(model=list(list(model="whi",kappa=5,var=2,s=4),"+",
    list(model="whi",kappa=9,var=6,s=1),"*",
    list(model="nugg",var=3,s=0)
   ),
  ti=1, me="ci")

## anisotropy
x(model=list(list(model="whi",kappa=5,var=2,a=1:4),"+",
    list(model="whi",kappa=1,var=6,a=2:5),"+",
    list(model="nugg",var=3,a=3:6)
   ),
  ti=2, me="ci")
