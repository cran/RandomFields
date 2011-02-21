# source("RFtest.Covariance.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

epsilon <-  0.000000001

tf <- function(x, l, i=1) {
  x <- as.matrix(x) 
  if (!is.null(l[[i]]$s)) {
    return(sqrt(rowSums(x^2))/l[[i]]$s)
  } else {
    return(sqrt(rowSums((x %*% matrix(l[[i]]$a, nrow=ncol(x)))^2)))
  }
}

tfst <- function(x, m) {
  if (is.matrix(x)) {
    if (is.null(m$a)) {
     stop("tfst -- no $aniso"); ERR
    }
    else {
      x <- x %*% matrix(m$a, nrow=ncol(x))
      return(cbind(sqrt(rowSums(x[,-ncol(x),drop=FALSE]^2)),
                    abs(x[,ncol(x)])))
    }
  } else {
    stop("tfst -- iso"); ERR
  }
}

st.fct <- function(x, m) {
  eval(parse(text=paste(letters[1:6], "<- m$u1$k[", 1:6,"];")))
  y <- abs(tfst(x, m))
  switch(e,
         psi <- sqrt(y[,ncol(y)]^c +1)^d,
         psi <- sqrt((d^(-1) * y[,ncol(y)]^c + 1) / (y[,ncol(y)]^c + 1)),
         psi <- sqrt(-log(y[,ncol(y)]^c + 1/d)/log(d))
         )
  if (ncol(y)==2) {
    cov <- c("stable", "whittle", "cauchy")[b]
    phi <- CovarianceFct(y[,1] / psi, model=cov, param=c(0, m$v, 0, 1, a))
  } else {
    stopifnot(ncol(y)==1); ERR
    phi <- 1
  }
  return(psi^(-f) * phi)
}

aniso <- function(dim) matrix(runif(dim^2), ncol=dim)
id <- function(dim) diag(dim)

x0 <- list(c(1:4), cbind(1:4, c(3,2,5,6)))
x1 <- list(cbind(0:4, c(0,3,2,5,6)), cbind(0:4, c(0,3,2,5,6), c(0,2:5)))


models <-
  list(list(list("$", var=3, scale=5, list(model="exp")),
            function(x, l) { exp(-tf(x,list(l),1)) * l$v},
            x0),       
       list(list("+", list("$", var=7, scale =5, u1=list("exp")),
                 list("$", var=3, scale =5, u1=list("gencau", a=2, b=3))),
            function(x, l) {
               exp(-tf(x,l,2)) * l[[2]]$v +
                (1+tf(x,l,3)^l[[3]]$u1$a)^(-l[[3]]$u1$b/l[[3]]$u1$a) * l[[3]]$v
            },
            x0),
        list(list("*", list("$", var=7, scale =5, u1=list("exp")),
                 list("$", var=3, scale =5, u1=list("gencau", a = 2,b=3))),
            function(x, l) {
               exp(-tf(x,l,2)) * l[[2]]$v *
                (1+tf(x,l,3)^l[[3]]$u1$a)^(-l[[3]]$u1$b/l[[3]]$u1$a) * l[[3]]$v
            },
            x0),
        list(list("*", list("$", var=7, aniso=aniso, u1=list("exp")),
                 list("$", var=3, aniso=aniso, u1=list("gencau", a=2,b=3))),
            function(x, l) {
               exp(-tf(x,l,2)) * l[[2]]$v *
                (1+tf(x,l,3)^l[[3]]$u1$a)^(-l[[3]]$u1$b/l[[3]]$u1$a) * l[[3]]$v
            },
            x0),
       list(list("$", var=7, aniso=aniso,
                 u1=list("nsst", delta= 6, phi=list("stable", 1.1), psi =list("gencau", alpha= 1.5, beta = .5))), #kappa=c(1.1, 1, 1.5,0.7, 1, 3)
            function(x,l) { x <-  tfst(x, l);  l$var * ( d<-(1 + Variogram( x[,ncol(x)], model=l$u1$psi)) ) ^
                            (-l$u1$d/2) * CovarianceFct( x[,1:(ncol(x)-1)]/sqrt(d), model=l$u1$phi) },
            x1),       
       list(list("+",
                 list("$", var=7, aniso=aniso,
                 u1=list("nsst", delta= 6, phi=list("stable", 1.1), psi =list("gencau", alpha= 1.5, beta = .5))),
                 list("$", var=5, aniso=aniso,
                      u1=list("nsst", delta= 6, phi=list("matern", nu=1.1),
                        psi =list("fractalB", alpha= 1.2)))),                      
            function(x,l){
              x1 <-  tfst(x, l[[2]]);
              x2 <-  tfst(x, l[[3]]);
              l$var * ( d<-(1 + Variogram( x1[,ncol(x1)], model=l[[2]]$u1$psi)) ) ^ (-l[[2]]$u1$d/2) *
                CovarianceFct( x1[,1:(ncol(x1)-1)]/sqrt(d), model=l[[2]]$u1$phi) +
              l$var * ( d<-(1 + Variogram( x2[,ncol(x2)], model=l[[3]]$u1$psi)) ) ^ (-l[[3]]$u1$d/2) *
                CovarianceFct( x2[,1:(ncol(x2)-1)]/sqrt(d), model=l[[3]]$u1$phi)  
            },
            x1))



  
 #      list(list("$", var=7, aniso=aniso,
 #         u1=list("nsst", kappa=c(1.1, 2, 1.5,0.7, 2, 3))),
 #           function(x,l) st.fct(x, l),
 #           x1),
 #      list(list("$",  var=7, aniso=aniso,
 #                u1=list("nsst", kappa=c(1.1, 3, 1.5,0.7, 3, 3))),
 #           function(x,l) st.fct(x, l),
 #           x1),
 #
 #      list(list("+",
 #                list("*",
 #                     list("$", var=7, aniso=aniso,
 #                         u1=list("nsst", kappa=c(1.3, 3, 1.2,0.5, 2, 4))),
 #                     list("$", var=5, aniso=aniso,
 #                          u1=list("nsst", kappa=c(1.2, 1, 1.5,0.7, 2, 3)))),
 #                list("$", var=5, aniso=aniso,
 #                     u1=list("nsst", kappa=c(1.2, 1, 1.5,0.7, 1, 3)))
 #                ),
 #           function(x,l) st.fct(x,l[[2]][[2]]) * st.fct(x,l[[2]][[3]]) + st.fct(x,l[[3]]),
 #           x1)       



{
for (model in models) {
  cat("\n",paste(c(lapply(model[[1]], f <- function(l) if(!is.function(l)) l[[1]])),"..")) #, recursive=TRUE
  for (x in model[[3]]) {
    f <- model[[2]]
    M <- model[[1]]
    for (i in 1:length(model[[1]])){
      if (is.function(M$aniso))
        M$aniso <- M$aniso(ncol(as.matrix(x)))
      if (is.list(M[[i]]) && is.function(M[[i]]$aniso))
        M[[i]]$aniso <- M[[i]]$aniso(ncol(as.matrix(x)))
      if (is.list(M[[i]]))
        for (j in 1:length(M[[i]]))
          if (is.list(M[[i]][[j]]) && is.function(M[[i]][[j]]$aniso))
            M[[i]][[j]]$aniso <- M[[i]][[j]]$aniso(ncol(as.matrix(x)))
    }
    #print(M)
    if ( any(abs(CovarianceFct(as.matrix(x), model=M) - f(x, M)) > f(x, M) * epsilon)){
      #print(model)
      print("X")
      print(x)
      print("C code")
      print(CovarianceFct(x, model=M))
      print("f eval")
      print(f(x, M))
      print("-----------")
      print(abs(CovarianceFct(x, model=M) - f(x, M)))
      print("X-----------")
      print(f(x, M) * epsilon)
      print("differences"); #ERR
    } else cat(". OK")
  }
}
cat("\n")
}
