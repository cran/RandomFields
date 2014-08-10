
\name{RFsimulate.more.examples}
\alias{RFsimulate.more.examples}
\title{Further Examples for the Simulation of Random Fields}

\description{
 This man page will give a collection of basic examples for the use of
 \code{\link{RFsimulate}}.

 For other kind of random fields (binary, max-stable, etc.) or
 more sophisticated approches see \link{RFsimulateAdvanced}.

 See  \link{RFsimulate.sophisticated.examples} for further examples.
}


\seealso{
 \command{\link{RFsimulate}},
 \command{\link{RFsimulateAdvanced}}
 \link{RFsimulate.sophisticated.examples}
}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

\dontshow{\dontrun{
#############################################################
## ##
## Basic examples: Unconditional simulation ## 
## ##
#############################################################

#############################################################
## Example 1: Location formats / Plot method ##
#############################################################

RFgetModelNames(type="CovType") ## the complete list covariance models
## and variogram models
## our choice is RMstable

## define the model:
model <- RMtrend(mean=0.5) + # mean
RMstable(alpha=1, var=4, scale=10) +
# see help("RMstable")
# for additional arguments
RMnugget(var=1) # nugget

## define the locations:
from <- 0
to <- 20
step <- 1 ## nicer, but also time consuming if step <- 0.1
x.seq <- seq(from, to, step) 
y.seq <- seq(from, to, step)

## simulate and get image of output:
simulated <- RFsimulate(model, x=x.seq, y=y.seq)
plot(simulated)
## for comparison only:
image(x.seq, y.seq, RFspDataFrame2conventional(simulated)$data)

#############################################################
## ... using seq(from, by, length.out) notation
len <- 21 ## more generally len <- (to-from)/step+1
x.short <- y.short <- c(from, step, len)
short.matrix <- cbind(x.short, y.short)
simulated <- RFsimulate(model = model, x=short.matrix)
plot(simulated) 

#############################################################
## ... using GridTopology for the locations:
len <- 21 ## more generally len <- (to-from)/step+1
gridtop <- GridTopology(c(from,from), c(step,step), c(len,len))
simulated <- RFsimulate(model = model, x=gridtop)
plot(simulated) 

############################################################
## arbitrary points
x <- runif(100, max=20) 
y <- runif(100, max=20) # 100 points in 2 dimensional space
simulated <- RFsimulate(model = model, x=x, y=y)
plot(simulated)

#############################################################
## using the 1-dimensional plot routine
## simulate 1-dimensional random field first
x.seq <- seq(from, to, step) # grid
simulated <- RFsimulate(model = model, x=x.seq)
plot(simulated) 


#############################################################
## Example 2: Simulation several realizations at once      ##
## Access to simulated data                                ##
#############################################################

model <- RMstable(alpha=1.5)
step <- 1 ## nicer, but also time consuming if step <- 0.1
x.seq <- seq(0, 20, step) 
y.seq <- seq(0, 20, step) # grid
simulated <- RFsimulate(model, x=x.seq, y=y.seq,
n=4) # 4 realizations at once
plot(simulated) 
summary(simulated@data$variable1.n2)
# summary of simulated univariate data of 2nd realization


#############################################################
## Example 3: simulating with trend                        ##
#############################################################

##s function 
model <- RMexp(var=0.3) +
RMtrend(arbitraryfct = function(x,y) 2*sin(x)*cos(y))
x.seq <- y.seq <- seq(-5,5,0.1)
simulated <- RFsimulate(model, x=x.seq,y=y.seq)
plot(simulated)

#############################################################
## with linear trend surface: 3x-y 

model1 <- RMexp() + RMtrend(plane=c(3,-1), fctcoeff=1)
# or equivalently:
model2 <- RMexp() + RMtrend(arbitraryfct=function(x,y) 3*x-y)

simulated <- RFsimulate(model1, x=x.seq, y=y.seq)
persp(x.seq,y.seq, RFspDataFrame2conventional(simulated)$data,
phi=30, theta=-3)


#############################################################
## Example 4: Brownian motion (using Stein's method) ##
#############################################################

# Brownian motion (1 dimensional)
alpha <- 1 # in [0,2)
x.seq <- seq(0, 10, 0.001)
simulated <- RFsimulate(model = RMfbm(alpha=alpha),
x=x.seq) 
plot(simulated) 

 
#############################################################
## Example 5: Models that depend on "submodels"; ##
## Combining models / Operators on models ##
#############################################################
 
RFgetModelNames(operator=TRUE) ## list all models
## that are an operator;
## our choice is RMcoxisham
D <- as.matrix(1) # a 1x1-correlation matrix for RMcoxisham
submodel <- RMwhittle(nu=0.3)  # submodel on which the operator
##                                model RMcoxisham will be applied
model <- RMcoxisham(submodel, mu=0, D=D, beta=2)
x.seq <- y.seq <- seq(-10,10,0.1)
simulated <- RFsimulate(model = model,
x=x.seq, y=y.seq) 
plot(simulated)

#############################################################
## further nesting of models is possible: 

#model2 <- RMcoxisham(RMintexp(RMdewijsian(alpha=1)),
# mu=0, D=D, beta=2)
#simulated <- RFsimulate(model = model2,
# x=x.seq, y=y.seq) 
#plot(simulated)

#############################################################
## addition of random fields using RMplus

model1 <- RMexp(var=5) + RMwhittle(nu=1, var=5)
## Alternatively, the common variance argument var=5
## can be included in the RMplus model:
model2 <- RMplus(C0 = RMexp(), C1 = RMwhittle(nu=1), var=5)
x.seq <- y.seq <- seq(-10,10,0.5)

simulated1 <- RFsimulate(model = model1, x=x.seq, y=y.seq) 
simulated2 <- RFsimulate(model = model2, x=x.seq, y=y.seq)
# compare (should give the same
plot(simulated1)
plot(simulated2)
sum(abs(simulated1@data - simulated2@data)) # should be numerically zero
 

#############################################################
## Example 6: A bivariate random field ##
#############################################################

RFgetModelNames(vdim=2) ## list all bivariate models
## our choice is RMbiWM
model <- RMbiwm(nudiag=c(1.3, 0.7), nured=2.5,
s=c(1, 1, 1), cdiag=c(0.7, 0.8), rhored=1)
x.seq <- y.seq <- seq(-10,10,0.5)
simulated <- RFsimulate(model = model, x=x.seq, y=y.seq) 
plot(simulated)



#############################################################
## ##
## Basic examples: Conditional simulation ## 
## ##
#############################################################

#############################################################
## Example 7: ways to pass given data ##
#############################################################

# simulate given locations and corresponding data
# (simulate measurements)
x <- runif(n=100, min=-1, max=1)
y <- runif(n=100, min=-1, max=1)
data <- RFsimulate(model = RMexp(), x=x, y=y, grid=FALSE)
# locations for conditional simulation
x.seq.cond <- y.seq.cond <- seq(-1.5,1.5,length=100)

#############################################################
## pass given data and locations as SP4-object
## given.data is an SP4-object

cond.simulated <- RFsimulate(RMexp(), x=x.seq.cond, y=y.seq.cond, data=data)

#############################################################
## or equivalently: pass given data and locations as a matrix

cond.simulated.2 <- RFsimulate(RMexp(), x=x.seq.cond, y=y.seq.cond, data = cbind(x, y, data@data))

all.equal(cond.simulated, cond.simulated.2) ## TRUE
plot(cond.simulated, data)
 

#############################################################
## multiple realizations

# simulate corresponding data twice (2 measurements)
given.data.2realize <- RFsimulate(model = RMwhittle(nu=1.1),
x=x.seq.cond, y=y.seq.cond,
data=data, n=2)
plot(given.data.2realize, data)
 
#############################################################
## simulation not on a grid

x.cond <- runif(1000, -2, 2)
y.cond <- runif(1000, -2, 2)
cond.simulated <- RFsimulate(model=RMwhittle(nu=1.4),
x=x.cond, y=y.cond, grid=FALSE, data=data)
plot(cond.simulated, data)


#############################################################
## Example 8: allow measurement errors ##
#############################################################

# err.model specifies the error model, typically RMnugget
cond.simulated <- RFsimulate(model=RMexp(),
x=x.seq.cond, y=y.seq.cond,
data=data,
err.model=RMnugget(var=1))
plot(cond.simulated, data)

}}
\dontshow{FinalizeExample()}

}