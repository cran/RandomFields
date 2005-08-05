
f <- expression(exp(-x^a))
Df <- D(f, "x")
D2f <- D(Df, "x")

rmin <- function(a=1, s=1) {
  x <- 0
  g00 <- eval(f) 
  x <- s
  g0 <- eval(f)
  g1 <- eval(Df)
  g2 <- eval(D2f)
  
  nenner <- (0.5 * g2 + g00 -g0)
  if  return(NA)
  if (nenner < 0) return(Inf)
  return((g2 -g1) / nenner  -1)
}

a <- seq(1.5, 1.95, 0.05)
s <- seq(0.0001, 4.5, 0.1)
z <- matrix(nrow=length(s), ncol=length(a))
for (i in 1:length(a))
  for (j in 1:length(s))
    z[j, i] = pmax(1, rmin(a[i], s[j]))
image(s, a, z, zlim=c(-2, 2))








df <- D(expression((1 + x^(2 *a))^(-b/a)), "x")
ddf <- D(df, "x")
a <- 1.2
b <- 2.9
x <- seq(0.00001, 5, 0.001)
z <- eval(ddf)
t <- x^2
zz <- t^(a-1) / (1+t^a)^(b/a+2) * b * (-(2 * a - 1) + t^a  * (1+2*b))
plot(z, zz, col=rainbow(length(x)))


z <- eval(df)
zz <- 2 * x  * (- t^(a-1) / (1 + t^a)^(b/a+1))
plot(z, zz)
