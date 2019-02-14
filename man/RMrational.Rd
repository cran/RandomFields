\name{RMrational}
\alias{RMrational}
\title{Rational function}
\description{
  Defines a simple rational function.
    
  \deqn{f(h) = (a_1 + a_2 z(h)) / (1 + z(h))}
  where
  \deqn{z(h) = h^\top A A^\top h}
}
 

\usage{
RMrational(A, a) 
}
\arguments{
 \item{A}{a \eqn{d \times d} matrix}
 \item{a}{a vector of one or two components; the second component has
   default value zero.}
}

\value{
 \command{\link{RMrational}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}. 
}

\me
\seealso{
 \command{\link{RMmodel}},
  \link{S10}.
 }


\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
# see S10
\dontshow{FinalizeExample()}
}
