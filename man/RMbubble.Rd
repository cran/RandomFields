\name{RMbubble}
\alias{RMbubble}
\title{Bubble model for arbitrary areas of scales}
\description{
  A model that allows for arbitray areas of scale applied to an isotropic
  model, i.e.
  \deqn{
    C(x,y) = \phi(\|x -y \| / s) 
  }
  as long as \eqn{s_x = s_y = s}. Here,
  \eqn{s_x} is the scaling at location \eqn{x},

  The cross-correlations between areas of
  different scales are given through a modified distance
  \eqn{d}. Let \eqn{z_{s}} be a finite
  subset of \eqn{R^d} depending on the scale \eqn{s}.
  Let \eqn{w_u} be a weight for an auxiliary point \eqn{u\in z_{s}}
  with \eqn{\sum_{u \in z_s} w_u = 1}. 
  Let \eqn{\tau_x = s_x^{-2}}. Then
  \deqn{d^2(x, y) = \min\{\tau(x), \tau(y)\} \|x - y\|^2
    + \sum_{\xi \in_{span(\tau(x), \tau(y))}}
      \sum_{u \in z_{\xi^{-0.5}}} w_u
    \|x - u\|^2 \Delta \xi
  }
  Here, \eqn{span(\tau(x), \tau(y))} is the finite set of values
  \eqn{s^{-2}} that are realized on the locations of interest
  and \eqn{\Delta \xi} is the difference of two
  realized and ordered values of the scaling \eqn{s}.
  
}
\usage{
RMbubble(phi, scaling, z, weight, minscale, barycentre,
         var, scale, Aniso, proj)}
\arguments{
 \item{phi}{isotropic submodel}
 \item{scaling}{model that gives the non-stationary scaling \eqn{s_x}}
 \item{z}{
   matrix of the union of all \eqn{z_s}.
   The number of rows equals the dimension of the field.
   If not given, the locations with non-vanishing gradient are taken.
 }
 \item{weight}{vector of weights \eqn{w}
   whose length equals the number of columns of \code{z}.
   The points given by \code{z} might be weighted.
 }
 \item{minscale}{vector for partioning \eqn{z} into classes \eqn{z_s}.
   Its length equals the number of columns of \code{z}. The vector
   values must be descending. See details.
   If not given then \eqn{z_s=}\code{z} for all \eqn{s}.
   Else see details.
 }
 \item{barycentre}{logical. If \code{FALSE} and \code{z} is not given,
   the reference locations are those with non-vashing gradient. If
   \code{TRUE} then, for each realized value of the scale, the
   barycentre of the corresponding reference locations is used instead
   of the reference locations themselves. This leads to higher
   correlations, but also to highly non-stationary cross-correlation
   between the areas of different scale.

   The argument has no effect when \code{z} is given.

   Default: \code{FALSE}.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 \command{\link{RMbubble}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}

\details{
   \code{minscale} gives the minimal scale \eqn{s} value above which
   the corresponding points \code{z} define the set \eqn{z_s}.
   The validity of the set \eqn{z_s} ends with the next lower value
   given.

    Let \code{minscale = (10, 10, 10, 7, 7, 7, 0.5)}. Then for some
   \eqn{d}-dimensional vectors \eqn{z_1,\ldots, z_7} we have

   \deqn{z_s = \{ z_1, z_2, z_3 \}, s \ge 10}
   \deqn{z_s = \{ z_4, z_5, z_5 \},  7 \ge s < 10}
   \deqn{z_s = \{ z_7 \}, s \ge 0.5}

   Note that, in this case, all realized scaling values must be
   \eqn{\ge 0.5}. Note further, that the weights for the subset must
   sum up to one, i.e.
   \deqn{w_1+w_2 +w_3=w_4 + w_5 + w_6 = w_7 = 1.}
}

\note{This model is defined only for grids.}

\references{
  \itemize{
    \item Bonat, W.H. , Ribeiro, P. Jr. and Schlather, M. (2019)
    Modelling non-stationarity in scale. In preparation.
  }
}


\me

\seealso{
  \command{\link{RMSadvanced}},
  \command{\link{RMblend}},
  \command{\link{RMscale}}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- seq(0,1, if (interactive()) 0.02 else 0.5)
d <- sqrt(rowSums(as.matrix(expand.grid(x-0.5, x-0.5))^2))
d <- matrix(d < 0.25, nc=length(x))
image(d)

scale <- RMcovariate(data=as.double(d) * 2 + 0.5, raw=TRUE)

## two models:
## the frist uses the standard approach for determining the
##           reference point z, which is based on gradients
## the second takes the centre of the ball
model1 <- RMbubble(RMexp(), scaling=scale)
model2 <- RMbubble(RMexp(), scaling=scale, z=c(0.5, 0.5))
model3 <- RMbubble(RMexp(), scaling=scale, barycentre=TRUE) # approx. of model2

## model2 has slightly higher correlations than model1:
C1 <- RFcovmatrix(model1, x, x)
C2 <- RFcovmatrix(model2, x, x)
C3 <- RFcovmatrix(model3, x, x) 
print(range(C2 - C1))
dev.new(); hist(C2 - C1)
print(range(C3 - C2)) # only small differences to C2
print(mean(C3 - C2))
dev.new(); hist(C3 - C2)

plot(z1 <- RFsimulate(model1, x, x))
plot(z2 <- RFsimulate(model2, x, x))
plot(z3 <- RFsimulate(model3, x, x)) # only tiny differences to z2


## in the following we compare the standard bubble model with
## the models RMblend, RMscale and RMS (so, model2 above
## performs even better)
biwm <- RMbiwm(nudiag=c(0.5, 0.5), nured=1, rhored=1, cdiag=c(1, 1), 
                s=c(0.5, 2.5, 0.5))
blend <- RMblend(multi=biwm, blend=RMcovariate(data = as.double(d), raw=TRUE))
plot(zblend <- RFsimulate(blend, x, x)) ## takes a while ...
Cblend <- RFcovmatrix(blend, x, x)

Mscale <- RMscale(RMexp(), scaling = scale, penalty=RMid() / 2)
plot(zscale <- RFsimulate(Mscale, x, x))
Cscale <- RFcovmatrix(Mscale, x, x)

Mscale2 <- RMscale(RMexp(), scaling = scale, penalty=RMid() / 20000)
plot(zscale2 <- RFsimulate(Mscale2, x, x))
Cscale2 <- RFcovmatrix(Mscale2, x, x)

S <- RMexp(scale = scale)
plot(zS <- RFsimulate(S, x, x))
CS <- RFcovmatrix(S, x, x)

print(range(C1 - CS))
print(range(C1 - Cscale))
print(range(C1 - Cscale2))
print(range(C1 - Cblend))
dev.new(); hist(C1-CS)     ## C1 is better
dev.new(); hist(C1-Cscale) ## C1 is better
dev.new(); hist(C1-Cscale2) ## both are equally good. Maybe C1 slightly better
dev.new(); hist(C1-Cblend) ## C1 is better




\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
