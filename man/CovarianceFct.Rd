\name{CovarianceFct}
\alias{Variogram}
\alias{PrintModelList}
\alias{GetModelNames}
\title{Covariance And Variogram Models}
\description{
  \code{CovarianceFct} returns the values of an isotropic covariance function
  
  \code{Variogram} returns the values of an isotropic variogram model

  \code{PrintModelList} prints the list of currently implemented models

  \code{GetModelNames} returns a list of currently implemented models
 }
\usage{
CovarianceFct(x,model,param,dim=1)
Variogram(x,model,param,dim=1)
PrintModelList()
GetModelNames()
}
\arguments{
  \item{x}{a vector of distances at which the covariance function
    or variogram should be evaluated}
  \item{model}{character; name of the covariance function or
    variogram model; see below,
    or type \code{PrintModelList()} for all
    options}
  \item{param}{vector of parameters; 
    \code{param=c(NA,variance,nugget,scale,...)}, in this order;\cr
    The dots \code{...} stand for additional parameters of the
    model
  }
  \item{dim}{dimension of the space in which the model is applied}
}
\details{
  The first component of \emph{param} is reserved for the \code{mean}
  of a random field and thus ignored in the evaluation of the covariance
  function or variogram.  The parameters mean, variance, nugget, and scale
  must be given in this order; additional
  parameters have to be supplied in case of a parametrised class of
  models (e.g. \code{hyperbolic}, see below),
  in the order \eqn{\kappa_1}{a}, \eqn{\kappa_2}{b}, \eqn{\kappa_3}{c}.
 
  The implemented models are in standard notation of a covariance
  function (variance 1,
  nugget 0, scale=1) and for positive real arguments \eqn{x}{x}:
  
  \itemize{
    \item \code{bessel}
    \deqn{C(x)=2^a \Gamma(a+1)x^{-a} J_a(x)}{C(x)=
      2^a Gamma(a+1)x^(-a) J_a(x)}
    The parameter \eqn{\kappa}{a} is greater than or equal to
    \eqn{\frac{d-2}2}{(d-2)/2}, where \eqn{d}{d} is the
    dimension of the random field.
    
    \item \code{cauchy}
    \deqn{C(x)=\left(1+x^2\right)^{-\kappa}}{C(x)=(1+x^2)^(-a)}
    The parameter \eqn{\kappa}{a} is positive. 
    The model possesses two generalisations, the \code{gencauchy}
    model and the \code{hyperbolic} model.
    
    \item \code{cauchytbm}
    \deqn{C(x)=\left(1+\left(1-\frac{\kappa_2}{\kappa_3}
      \right)x^{\kappa_1}\right)
      \left(1+x^{\kappa_1}\right)^{-\frac{\kappa_2}{\kappa_1}-1}}{
      C(x)=
      (1+(1-b/c)x^a)(1+x^a)^(-b/a-1)}
    The parameter \eqn{\kappa_1}{a} is in (0,2],  \eqn{\kappa_2}{b}
    is positive, and \eqn{\kappa_3}{c} is an integer. 
    The model is valid for dimensions \eqn{d\le\kappa_3}{d<=c}.
    \cr
    It allows for simulating random fields where 
    fractal dimension and Hurst coefficient can be chosen
    independently. 
    It has negative correlations for \eqn{\kappa_2>1}{b>c} and large
    \eqn{x}{x}.
    
    \item \code{circular}
    \deqn{C(x)=
      \left(1-\frac 2\pi
      \left(x \sqrt{1-x^2} +
      \arcsin(x)\right)\right)
      1_{[0,1]}(x)}{
      C(x)=1-2/pi*(x sqrt(1-x^2)+asin(x))   if
      0<=x<=1, 0 otherwise}
    This isotropic covariance function is valid only for dimensions
    less than or equal to 2.
    
    \item \code{cone}\cr
    This model is used only for methods based on marked point processes
    (see \code{\link{RFMethods}}); it is defined only in two dimensions.
    The corresponding (boolean)
    function is a truncated cone with socle. The base has radius
    \eqn{\frac12}{1/2}. The model has three parameters, \eqn{\kappa_1}{a},
    \eqn{\kappa_2}{b}, and \eqn{\kappa_3}{c}:\cr
    \eqn{kappa_1}{a} gives the radius of the top circle of the cone, given
    as part of the socle radius; \eqn{kappa_1\in[0,1)}{a in [0,1)}.\cr
    \eqn{kappa_2}{b} gives the height of the socle.\cr
    \eqn{kappa_3}{c} gives the height of the truncated cone.\cr
    
    \item \code{cubic}
    \deqn{C(x)=(1- 7x^2+8.75x^3-3.5x^5+0.75 x^7)1_{[0,1]}(x)}{C(x)=
      1- 7 x^2 + 8.75 x^3 - 3.5 x^5 + 0.75 x^7  if 0<=x<=1,
      0 otherwise}
    This model is valid only for dimensions less than or equal to 3. 
    It is a 2 times differentiable covariance functions with compact
    support. %(See Chiles&Delfiner, 1998)
    
    \item \code{exponential}
    \deqn{C(x)=e^{-x},  \quad x\ge0}{C(x)=exp(-x)}
    This model is a special case of the \code{whittlematern} model
    (for \eqn{\kappa=\frac12}{a=1/2} there)
    and the \code{stable} class (for \eqn{\kappa=1}{a=1}).
    
    %%   \item\code{expPLUScirc}
    %%   \deqn{\kappa_1 e^{-x} +
    %%     (1-\kappa_1)
    %%     \left(1-\frac2\pi
    %%     \left(\frac xb\sqrt{1-(x/\kappa_2)^2}+\arcsin(x/\kappa_2)\right)
    %%     \right) 1_{[0,1]}(x)}{
    %%     a exp(-x) +
    %%     (1-a)(1-2/pi(x/b*sqrt(1-x^2/b^2)+asin(x/b)))   if
    %%     0<=x<=1, 0 otherwise}
    %%   The parameter \eqn{\kappa_1}{a} is in \eqn{[0,1]} and
    %%   \eqn{\kappa_2}{b} is positive. 
    %%   This isotropic covariance function is valid only for dimensions less
    %%   than or equal to 2.
   
    \item \code{gaussian}
    \deqn{C(x)=e^{-x^2}}{C(x)=exp(-x^2)}
    This model is a special case of the \code{stable} class
    (for \eqn{\kappa=2}{a=2} there).
    Note that the corresponding function for the random coins
    method (cf. the methods based on marked point processes in
    \code{\link{RFMethods}}) is
    \deqn{e^{- 2 x^2}.}{exp(-2 x^2).}
    See \code{gneiting} for an alternative model that does not have
    the disadvantages of the Gaussian model.
    
    \item \code{gencauchy} (generalised \code{cauchy})\cr
    \deqn{C(x)=(1+x^{\kappa_1})^{-\kappa_2/\kappa_1}}{C(x)=
      \left(1+x^a\right)^(-b/a)}
    The parameter \eqn{\kappa_1}{a} is in (0,2], and \eqn{\kappa_2}{b}
    is positive.
    \cr
    This model allows for simulating random fields where
    fractal dimension and Hurst coefficient can be chosen
    independently.

    \item \code{gengneiting} (generalised \code{gneiting}\cr
    if \eqn{\kappa_1=1}{a=1} then
    \deqn{C(x)=\left(1+(\kappa_2+1)x\right) * (1-x)^{\kappa_2+1}
      1_{[0,1]}(x)}{C(x)=[1 + (b+1) * x] * (1-x)^(b+1) if 0<=x<=1,
      0 otherwise}
    if \eqn{\kappa_1=2}{a=2} then
    \deqn{C(x)=\left(1+(\kappa_2+2)x+\left((\kappa_2+2)^2-1\right)x^2/3\right)
      (1-x)^{\kappa_2+2}  1_{[0,1]}(x)}{C(x)=
      [1 + (b+2) * x + ((b+2)^2-1) * x^2 / 3] * (1-x)^(b+2)
      if 0<=x<=1, 0 otherwise}
    if \eqn{\kappa_1=3}{a=3} then
    \deqn{C(x)=\left(1+(\kappa_2+3)x+\left(2(\kappa_2+3)^2-3\right)x^2/5
      +\left((\kappa_2+3)^2-4\right)(\kappa_2+3)x^3/15\right)(1-x)^{\kappa_2+3}
      1_{[0,1]}(x)}{C(x)=[1 + (b+3) * x +  (2 * (b+3)^2 - 3) * x^2 / 5
      + ((b+3)^2 - 4) * (b+3) * x^3 / 15] * (1-x)^(b+3)
      if 0<=x<=1, 0 otherwise}
    The parameter \eqn{kappa_1}{a} is a positive integer; here only the
    cases \eqn{kappa_1=1, 2, 3}{a=1, 2, 3} are implemented.
    The parameter \eqn{kappa_2}{b} is greater than or equal to
    \eqn{(d + 2\kappa_1 +1)/2}{(d + 2a +1)/2} where \eqn{d} is the
    dimension of the random field.
    
    
    \item \code{gneiting}
    \deqn{C(x)=\left(1 + 8 sx + 25 (sx)^2 + 32
      (sx)^3\right)(1-sx)^8 1_{[0,1]}(sx)}{C(x)=
      (1 + 8 s x + 25 s^2 x^2 + 32
      s^3 x^3)*(1-s x)^8   if 0<=x<=1, 0 otherwise}
    where
    \eqn{s=\frac{10\sqrt2}{47}\approx 0.3}{s = 10 sqrt(2) / 47 ~= 0.3}.
    This isotropic covariance function is valid only for dimensions less
    than or equal to 3. 
    It is a 6 times differentiable covariance functions with compact
    support.\cr
    It is an alternative to the \code{gaussian} model since
    its graph is visually hardly distinguishable from the graph of
    the Gaussian model, but possesses neither the mathematical and nor the
    numerical disadvantages of the Gaussian model.\cr
    This model is a special case of \code{gengneiting} (for
    \eqn{kappa_1=3}{a=3} and \eqn{kappa_2=5}{b=5} there).
    
    \item \code{gneitingdiff}
    \deqn{C(x)=\left(1 + 8 \frac x{\kappa_2}
      + 25 \frac {x^2}{\kappa_2^2}
      + 32 \frac {x^3}{\kappa_2^3}\right)
      \left(1-\frac{x}{\kappa_2}\right)^8 
      \;\frac{2^{1-\kappa_1}}{\Gamma(\kappa_1)}
      \,x^{\kappa_1} K_{\kappa_1}(x)1_{[0,\kappa_2]}(x)}{
      C(x)=(1 + 8 x/b + 25 (x/b)^2 + 32
      (x/b)^3)*(1-x/b)^8 * 2^{1-a} Gamma(a)^{-1} x^a K_a(x)   if
      0<=x<=b, 0 otherwise}
    This isotropic covariance function is valid only for dimensions less
    than or equal to 3. 
    The parameters \eqn{\kappa_1}{a} and \eqn{\kappa_2}{b} are
    positive.\cr
    This class of models with compact support
    allows for smooth parametrisation of the differentiability up to
    order 6.     
    
    \item \code{holeeffect}
    \deqn{C(x)=e^{-\kappa x} \cos(x),  \quad x\ge0}{C(x)=
      exp(-a x) cos(x) }
    This model is valid
    for dimension 1 iff \eqn{\kappa\ge1}{a>=0},
    for dimension 2 iff \eqn{\kappa\ge1}{a>=1},
    and for dimension 3 iff \eqn{\kappa\ge \sqrt{3}}{a >= sqrt(3)}.
    
    \item \code{hyperbolic}
    \deqn{C(x)=\frac{1}{\kappa_3^{\kappa_2}
	K_{\kappa_2}(\kappa_1 \kappa_3)}
      \left(\kappa_3^2 +x^2\right)^{{\kappa_2}/2}
      K_{{\kappa_2}}\left(
      \kappa_1 \left(\kappa_3^2 + x^2\right)^{1/2}\right),  \quad
      x>0}{C(x)=
      c^(-b) (K_b(a*c))^(-1) * (c^2 +x^2)^(0.5 b) *
      K_b(sqrt(a(c^2 + x^2)))}
    The parameters are such that\cr
    \eqn{\kappa_3\ge0}{c>=0},  \eqn{\kappa_1>0}{a>0 }  and 
    \eqn{\kappa_2>0,\quad}{b>0,   }
    or\cr
    \eqn{\kappa_3>0}{c>0 },  \eqn{\kappa_1>0}{a>0 }  and 
    \eqn{\kappa_2=0,\quad}{b=0,   }
    or\cr
    \eqn{\kappa_3>0}{c>0 },  \eqn{\kappa_1\ge0}{a>=0}, and 
    \eqn{\kappa_2<0}{b<0}.\cr
    Note that this class is over-parametrised; always one
    of the three parameters 
    \eqn{\kappa_1}{a}, \eqn{\kappa_3}{c}, and scale
    can be eliminiated in the formula. Therefore, one of these
    parameters should be kept fixed in any simulation study.
    \cr
    The model contains as special cases the \code{whittlematern}
    model and the \code{cauchy} model, for 
    \eqn{\kappa_3=0}{c=0} and \eqn{\kappa_1=0}{a=0}, respectively.
    
    \item matern\cr
    See \code{whittlematern}.
    
    \item \code{nugget}
    \deqn{C(x)=1_{\{0\}}(x)}{1(x==0)}
    Here \code{param[2]}, the variance,  must be zero!
    
    \item \code{pentamodel}
    \deqn{C(x)= \left(1 - \frac{22}3 x^2 +33 x^4 -
      \frac{77}2 x^5 + \frac{33}2
      x^7 -\frac{11}2 x^9 + \frac 56 x^{11}
      \right)1_{[0,1]}(x)}{C(x)=
      1 - 22/3 x^2 +33 x^4
      - 77/2 x^5 + 33/2 x^7 - 11/2 x^9 + 5/6 x^11  if 0<=x<=1,  
      0 otherwise}
    valid only for dimensions less than or equal to 3. 
    This is a 4 times differentiable covariance functions with compact
    support.
    %(See Chiles&Delfiner, 1998)
    
    \item \code{power}
    \deqn{C(x)= (1-x)^\kappa 1_{[0,1]}(x)}{C(x)=
      (1-x)^a   if 0<=x<=1, 0 otherwise}
    This covariance function is valid for dimension \eqn{d}{d} if
    \eqn{\kappa\ge\frac{d+1}2}{a >= (d+1)/2}. 
    For \eqn{\kappa=1}{a=1} we get the well-known triangle (or tent)
    model, which is valid on the real line, only.
    % proposition 3.8 in phd thesis tilmann gneiting
    % Golubov, Zastavnyi
    
    \item powered exponential\cr
    See \code{stable}.
    
    \item \code{qexponential}
    \deqn{C(x)=\frac{2 e^{-x}-\kappa e^{-2x}}{2-\kappa}}{
      C(x) = (2 exp(-x)-a exp(-2x))/(2-a)}
    The parameter \eqn{\kappa}{a} takes values in \eqn{[0,1]}{[0,1]}.   

    % \item rational quadratic model\cr
    %  See \code{cauchy} for \eqn{\kappa=1}{a=1}.
    % (Cressie)
    
    \item \code{spherical}
    \deqn{C(x)=\left(1- \frac32 x+\frac 12 x^3\right)
      1_{[0,1]}(x)}{C(x)=
      1 - 1.5 x + 0.5 x^3  if 0<=x<=1, 0 otherwise}
    This isotropic covariance function is valid only for dimensions
    less than or equal to 3.
    
    \item \code{stable}
    \deqn{C(x)=\exp\left(-x^\kappa\right)}{C(x)=exp(-x^a)}
    The parameter \eqn{\kappa}{a} is in \eqn{[0,2]}{[0,2]}.
    See \code{exponential} and \code{gaussian} for special cases.
    
    \item symmetric stable\cr
    See \code{stable}.
    
    \item tent model\cr
    See \code{power}.
    
    \item triangle\cr
    See \code{power}.
    
    \item \code{wave}
    \deqn{C(x)=\frac{\sin x}x, \quad x>0}{C(x)=sin(x)/x if x>0}
    This isotropic covariance function is valid only for dimensions less
    than or equal to 3.
    It is a special case of the \code{bessel} model
    (for \eqn{\kappa=2}{a=3}).
    
    \item \code{whittlematern}
    \deqn{C(x)=2^{1-\kappa} \Gamma(\kappa)^{-1} x^\kappa
      K_\kappa(x)}{C(x)=2^{1-a} Gamma(a)^{-1} x^a K_a(x),
    }
    The parameter \eqn{\kappa}{a} is positive.
    \cr
    This is the model of choice if the smoothness of a random field is to
    be parametrised. It is a special case of the
    \code{hyperbolic} model (for \eqn{\kappa_3=0}{c=0} there).
  }
  
  Let \eqn{\rm cov}{cov} be a model given in standard notation. 
  Then the covariance model
  applied with arbitrary variance,
  nugget, and scale equals
  \deqn{\rm \qquad nugget + variance * \rm cov( (\cdot)/ scale).
  }{nugget + variance * cov( (.)/scale).}
  
  For a given covariance function \eqn{\rm cov}{cov} the variogram
  \eqn{\gamma}{gamma} equals 
  \deqn{\gamma(x) = {\rm cov}(0) - {\rm cov}(x).}{
    gamma(x) = cov(0) - cov(x).}
  
  Note that the value of the covariance function or variogram
  depends also on
  \code{\link{RFparameters}()$PracticalRange}.  If the latter is
  \code{TRUE} then the covariance function is internally
  rescaled such that cov\eqn{(1)\approx 0.05}{(1)~=0.05} for standard
  parameters (\code{scale==1}).
  
  Some models allow certain parameter combinations only for certain
  dimensions.  As any model in \eqn{d}{d} dimensions is also valid in 1
  dimension, the default in \code{CovarianceFct} and \code{Variogram}
  is \code{dim=1}.
}
\value{
  \code{CovarianceFct} returns a vector of values of the covariance function.
  
  \code{Variogram} returns a vector of values of the variogram model.
  
  \code{PrintModelList} prints a table of the currently implemented covariance
  functions and the matching methods.
  \code{PrintModelList} returns \code{NULL}.

  \code{GetModelNames} returns a list of implemented models
}
\references{
  Overviews:
  \itemize{
    \item Chiles, J.-P. and Delfiner, P. (1999)
    \emph{Geostatistics. Modeling Spatial Uncertainty.}
    New York: Wiley.
    
    \item  Gneiting, T. and Schlather, M. (2001)
    Statistical modeling with covariance functions.
    \emph{In preparation.}
    
    \item Schlather, M. (1999) \emph{An introduction to positive definite
      functions and to unconditional simulation of random fields.}
    Technical report ST 99-10, Dept. of Maths and Statistics,
    Lancaster University.

    \item Schlather, M. (2001) Models for stationary max-stable
      random fields. \emph{Submitted.}
    
    \item Yaglom, A.M. (1987) \emph{Correlation Theory of Stationary and
      Related Random Functions I, Basic Results.}
    New York: Springer.
    
    \item Wackernagel, H. (1998) \emph{Multivariate Geostatistics.} Berlin:
    Springer, 2nd edition.
  }

  Cauchy models, generalisations and extensions
  \itemize{
    \item Gneiting, T. and Schlather, M. (2001)
    Stochastic models which separate fractal dimension and Hurst effect.
    \emph{Submitted.}
  }
 
  Gneiting's models
  \itemize{
    \item Gneiting, T. (1999)
    Correlation functions for atmospheric data analysis
    \emph{Q. J. Roy. Meteor. Soc., Part A} \bold{125}, 2449-2464. 
  }
   
  Holeeffect model
  \itemize{
    \item Zastavnyi, V.P. (1993)
    Positive definite functions depending on a norm,
    \emph{Russian Acad. Sci. Dokl. Math.} \bold{46}, 112-114. 
  }

  Hyperbolic model
  \itemize{
    \item Shkarofsky, I.P. (1968) Generalized turbulence space-correlation and
    wave-number spectrum-function pairs. \emph{Can. J. Phys.} \bold{46},
    2133--2153.
  }
  
  Power model
  \itemize{
    \item Golubov, B.I. (1981) On Abel-Poisson type and Riesz means,
    \emph{Analysis Mathematica} \bold{7}, 161-184.
    
    \item  Zastavnyi, V.P. (2000) On positive definiteness of some
    functions, \emph{J. Multiv. Analys.} \bold{73}, 55-81.
  }
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
\seealso{\code{\link{EmpiricalVariogram}},
  \code{\link{RandomFields}},
  \code{\link{RFparameters}},
  \code{\link{ShowModels}}.}

\examples{
 PrintModelList()
 CovarianceFct(0:100, "bessel", c(NA,2,1,5,0.5))
}
\keyword{spatial}






