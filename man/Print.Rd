\name{Print}
\alias{Print}
\title{Nice print function returning also the names automatically}
\description{
  prints variable names and the values
}
\usage{
Print(...)
}
\arguments{
  \item{...}{any object that can be \command{print}-ed}
}
\value{
  prints the names and the values; for vectors \command{cat}
  is used and for lists \command{str}
}
\examples{
  a <- 4
  b <- list(c=5, g=7)
  m <- matrix(1:4, nc=2)
  Print(a, b, m)
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\keyword{print}
