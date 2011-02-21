model <- list("a", list("y", 6, 7, c(0, NA, NA)), list(5, 6, NA,8 ))


getNApos <- function(l, n=0) {
#  print(l)
  if (is.list(l)) {
    ret <- NULL
    for (i in 1:length(l)) {
      nas <- getNApos(l[[i]], n)
      n <- nas[[2]]
      ret <- c(ret, nas[[1]])
    }
#   str(list("r", list(ret, n)))
    return(list(ret, n))
  } else {
    x <- which(is.na(l))
#    print(c(length(l), n))
#    str(list("x", list(x + n, length(l) + n)))
    return(list(x + n, length(l) + n))
  }
}


putpos <- function(l, pos, what, n=0) {
#  print("neu")
#  str(l)
  if (is.list(l)) {
    for (i in 1:length(l)) {
#      print(i)
#      str(l)
      lneu <- putpos(l[[i]], pos, what, n)
      l[[i]] <- lneu[[1]]
      n <- lneu[[2]]
    }
#    str(list("r", list(l, n)))
    return(list(l, n))
  } else {
    x <- n + 1:length(l)
    idx <- x %in% pos 
    widx <- pos %in%  x
    if (any(idx)) {
      l[idx] <- what[widx]
    }
#    str(list("x", list(l, n + length(l))))
    return(list(l, n + length(l)))
  }
}

getpos <- function(l, pos, n=0) {
#  print("neu")
#  str(l)
  if (is.list(l)) {
    ret <- NULL
    for (i in 1:length(l)) {
#      print(i)
#      str(l)
      lneu <- getpos(l[[i]], pos, n)
      ret <- c(ret, lneu[[1]])
      n <- lneu[[2]]
    }
#    str(list("r", list(l, n)))
    return(list(ret, n))
  } else {
    x <- n + 1:length(l)
    idx <- x %in% pos
    ret <- l[idx]
#    str(list("x", list(l, n + length(l))))
    return(list(ret, n + length(l)))
  }
}


cat(rep("\n", 15))
p <- getNApos(model)[[1]]
(pp <- putpos(model, p, c(23, 24, 25)))
(pp2 <- getpos(pp[[1]], p))
