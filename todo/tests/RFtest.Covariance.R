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
  eval(parse(text=paste(letters[1:6], "<- m$k[", 1:6,"];")))
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
  list(
       list(list(list(model="exp", var=3, scale =5)),
            function(x, l) { exp(-tf(x,l,1)) * l[[1]]$v},
            x0),
       list(list(list(model="exp", var=7, scale =5),"+",
                 list(model="gencau", var=3, scale =5, k=c(2,3))),
            function(x, l) {
               exp(-tf(x,l,1)) * l[[1]]$v +
                (1+tf(x,l,3)^l[[3]]$k[1])^(-l[[3]]$k[2]/l[[3]]$k[1]) * l[[3]]$v
            },
            x0),
        list(list(list(model="exp", var=7, scale =5),"*",
                 list(model="gencau", var=3, scale =5, k=c(2,3))),
            function(x, l) {
               exp(-tf(x,l,1)) * l[[1]]$v *
                (1+tf(x,l,3)^l[[3]]$k[1])^(-l[[3]]$k[2]/l[[3]]$k[1]) * l[[3]]$v
            },
            x0),
        list(list(list(model="exp", var=7, aniso=aniso),"*",
                 list(model="gencau", var=3, aniso=aniso, k=c(2,3))),
            function(x, l) {
               exp(-tf(x,l,1)) * l[[1]]$v *
                (1+tf(x,l,3)^l[[3]]$k[1])^(-l[[3]]$k[2]/l[[3]]$k[1]) * l[[3]]$v
            },
            x0),
       list(list(list(model="nsst", var=7, aniso=aniso,
                      kappa=c(1.1, 1, 1.5,0.7, 1, 3))),
            function(x,l) st.fct(x, l[[1]]),
            x1),
       list(list(list(model="nsst", var=7, aniso=aniso,
                      kappa=c(1.1, 2, 1.5,0.7, 2, 3))),
            function(x,l) st.fct(x, l[[1]]),
            x1),
       list(list(list(model="nsst", var=7, aniso=aniso,
                      kappa=c(1.1, 3, 1.5,0.7, 3, 3))),
            function(x,l) st.fct(x, l[[1]]),
            x1),
      list(list(list(model="nsst", var=7, aniso=aniso,
                      kappa=c(1.3, 3, 1.2,0.5, 2, 4)),"+",
                list(model="nsst", var=5, aniso=aniso,
                     kappa=c(1.2, 1, 1.5,0.7, 2, 3))),
            function(x,l) st.fct(x, l[[1]]) +  st.fct(x, l[[3]]),
            x1),
       list(list(list(model="nsst", var=7, aniso=aniso,
                      kappa=c(1.3, 3, 1.2,0.5, 2, 4)),"*",
                list(model="nsst", var=5, aniso=aniso,
                     kappa=c(1.2, 1, 1.5,0.7, 2, 3)), "+",
                list(model="nsst", var=5, aniso=aniso,
                     kappa=c(1.2, 1, 1.5,0.7, 1, 3))),
            function(x,l)st.fct(x,l[[1]]) *  st.fct(x,l[[3]]) + st.fct(x,l[[5]]),
            x1)
       )
  
{
for (model in models) {
  cat("\n",paste(c(lapply(model[[1]], function(l) l[[1]]), recursive=TRUE)),"..")
  for (x in model[[3]]) {
    f <- model[[2]]
    M <- model[[1]]
    for (i in 1:length(model[[1]]))
      if (is.list(M[[i]]) && is.function(M[[i]]$aniso))
        M[[i]]$aniso <- M[[i]]$aniso(ncol(as.matrix(x)))
    #print(M)
    if ( any(abs(CovarianceFct(as.matrix(x), M) - f(x, M)) > f(x, M) * epsilon)){
      #print(model)
      print("X")
      print(x)
      print("C code")
      print(CovarianceFct(x, M))
      print("f eval")
      print(f(x, M))
      print("-----------")
      print(abs(CovarianceFct(x, M) - f(x, M)))
      print("X-----------")
      print(f(x, M) * epsilon)
      stop("differences"); ERR
    } else cat(". OK")
  }
}
cat("\n")
}
