#    source("RFtest.model.R")

if (EXTENDED.TESTING <- file.exists("source.R")) { source("source.R")
} else if (file.exists(f <- "~/R/RF/RandomFields/tests/source.R")) source(f) 

x <- function(...) {
  print(unlist(list(...)))
  str(PrepareModel(...))
  cat("--------------------------------\n")
  str(m <- PrepareModel(...))
  invisible(m)
}
#x <- function(...) str(PrepareModel(...))

# stop("")
x(model="nugg",param=c(0,2,3,4))
x(model="nugg",param=c(0,0,3,4))
x(model="nugg",param=c(0,3,0,4))
x(model="nugg",param=c(0,3,0,4),me="ci")

x(model="sphe",param=c(1,2,0,4))
x(model="sphe",param=c(1,2,3,4))
x(model="sphe",param=c(1,2,3,4), me="ci")

x(model="sphe",param=cbind(c(2,4),c(8,1)), me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0)), me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0)), trend=1, me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0)), trend=c(1,2),me="ci")#not simple
x(model="sphe",param=cbind(c(2,4),c(3,0),c(7,0)), me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0),c(8,1)), me="ci")
x(model="sphe",param=cbind(c(2,4),c(3,0),c(8,1)), trend=1)

x(model="whi",param=c(1,2,3,4,5))
x(model="whi",param=c(1,2,0,4,5))
x(model="whi",param=c(1,2,3,4,5), me="ci")
x(model="whi",param=c(1,2,0,4,5), me="ci")

x(model="whi",param=cbind(c(2,4,5),c(8,1,1)), me="ci")
x(model="whi",param=cbind(c(2,4,5),c(3,0,0)), me="ci", trend=9)
x(model="whi",param=cbind(c(2,4,5),c(3,0,0),c(7,0,0)), me="ci")
x(model="whi",param=cbind(c(2,4,5),c(3,0,0),c(8,1,1)), me="ci")

# list to old, vector
x(model=list(list(model="whi",k=5,var=2,s=4)), me="ci")
x(model=list(list(model="whi",k=5,var=2,s=4),"+",
    list(model="whi",k=1,var=3,s=0)
    ),
   me="ci")
## note difference to the next result !!!!!!!!!
xx <- PrepareModel(model= list(list(model="whi",k=5,var=2,s=4),
                     "+",
                     list(model="whi",k=1,var=3,s=0)),
                   me="ci")
x(model=xx$mo, pa=xx$pa,  me=xx$me)


# list to nested
x(model=list(list(model="whi",k=5,var=2,s=4),"+",
    list(model="nugg",var=3,s=0),"+",
    list(model="nugg",var=7,s=0)
    ),
  me="ci")

x(model=list(list(model="whi",k=5,var=2,s=4),"+",
    list(model="whi",k=1,var=6,s=1),"+",
    list(model="nugg",var=3,s=0)
   ),
  me="ci")

x(model=list(list(model="whi",k=5,var=2,s=4),"+",
    list(model="whi",k=1,var=6,s=1),"+",
    list(model="nugg",var=3,s=0)
   ),
  me="ci")


# multiplicative model
x(model=list(list(model="whi",k=5,var=2,s=4),"+",
    list(model="whi",k=9,var=6,s=1),"*",
    list(model="nugg",var=3,s=0)
   ),
  me="ci")

## anisotropy
x(model=list(list(model="whi",k=5,var=2,a=1:4),"+",
    list(model="whi",k=1,var=6,a=2:5),"+",
    list(model="nugg",var=3,a=3:6)
   ),
  me="ci")

# geht auch
x(model=list(list(model="whi",k=5,var=2,s=4),"+",
    list(model="whi",k=1,var=6,s=1),"*",
    list(model="nugg",var=3,s=0)
   ),
   me="ci")

