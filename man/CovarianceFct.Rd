\name{CovarianceFct}
\alias{CovarianceFct}
\alias{Covariance}
\alias{Variogram}
\alias{CovMatrix}
\title{Basic Covariance And Variogram Models}
\description{
  \code{CovarianceFct} returns the values of a covariance function;
  see \command{\link{Covariance}} for sophisticated models
  
  \code{Variogram} returns the values of a variogram model
}
\usage{
Covariance(x, y=NULL, model, param=NULL, dim=ifelse(is.matrix(x),ncol(x),1),
               Distances, fctcall=c("Cov", "Variogram", "CovMatrix"))
CovarianceFct(...)
CovMatrix(...)

Variogram(x, model, param, dim=ifelse(is.matrix(x),ncol(x),1))
}
\arguments{
  \item{x}{vector or \eqn{(n \times \code{dim})}{(n x
      \code{dim})}-matrix.  In particular,
    if the model is isotropic or \code{dim=1} then \code{x}
    is a vector.}
  \item{y}{second vector or matrix in case of non-stationary covariance
    functions}
  \item{model}{
    for basic models, \code{model} is one of the names given in the
    Details.
    \cr
%    For advanced users we refer to \link{modeldef}.
  }
  \item{param}{
    The simplest form of \code{param} is the vector 
    \code{param=c(mean,variance,nugget,scale,...)}, in this order;\cr
    The dots \code{...} stand for additional parameters of the
    model, e.g. the smoothing parameter in the \code{whittle}
    model.
    Within this function \code{mean} is not interpreted and can
    take an arbitrary value.
      \cr
    %For advanced users we refer to \link{modeldef}.  
  }
  \item{dim}{dimension of the space in which the model is applied}
  \item{Distances}{for covariance matrices, the lower triangular part
    of the distance matrix can be given instead of the values \code{x}
    themselves}
  \item{fctcall}{this parameter should not be changed by the user}
  \item{...}{The function \code{CovarianceFct} is identical to the function
    \code{Covariance}.
  }
}
\details{
  \bold{Here, only the basic, isotropic models are listed;
    see \command{\link{sophisticated}}  models for nonisotropic and hyper models.
  }
  
  See \command{\link{GetModel}} for commands in R to get information
  about implemented models and currently used ones.
  
  The implemented models are in standard notation for a
  covariance function (variance 1, nugget 0, scale 1) and for positive
  real arguments \eqn{h}: 
  
  \itemize{
    \item \code{+} see `\link{sophisticated}'
    \item \code{*} see `\link{sophisticated}'
    \item \code{$} see `\link{sophisticated}'
    \item \code{ave1} see `\link{sophisticated}'
    \item \code{ave2} see `\link{sophisticated}'
    
    \item \code{bessel}
    \deqn{C(h)=2^\nu \Gamma(\nu+1)h^{-\nu} J_\nu(h)}
    The parameter \eqn{\nu} is greater than or equal to
    \eqn{\frac{d-2}2}{(d-2)/2}, where \eqn{d}{d} is the
    dimension of the random field.

    \item Brownian motion\cr
    see \code{fractalB} 
    
    \item cardinal sine\cr
    see \code{wave}
    
    \item \code{cauchy (normal scale mixture)}
    \deqn{C(h)=\left(1+h^2\right)^{-\beta}}{C(h)=(1+h^2)^(-\beta)}
    The parameter \eqn{\beta} is positive. 
    The model possesses two generalisations, the \code{gencauchy}
    model and the \code{hyperbolic} model.
    See also \code{nonstatcauchy} in \command{\link{Covariance}}.
    
    \item \code{cauchytbm}
    \deqn{C(h)= (1+(1-\beta/\gamma)h^\alpha)(1+h^\alpha)^(-\beta/\alpha-1)}
    The parameter \eqn{\alpha} is in (0,2] and  \eqn{\beta}
    is positive.
    The model is valid for dimensions \eqn{d\le\gamma};
    this has been shown for integer \eqn{\gamma}, but the
    package allows real values of \eqn{\gamma}.
    \cr
    It allows for simulating random fields where 
    fractal dimension and Hurst coefficient can be chosen
    independently. 
    It has negative correlations for \eqn{\beta>\gamma} and large
    \eqn{h}{h}.

    This model is equivalent to the model
    \code{list("tbm3", n=gamma, list("gencauchy", alpha=alpha,
      beta=beta))}
    
    
    \item \code{circular}
    \deqn{C(h)=
      \left(1-\frac 2\pi
      \left(h \sqrt{1-h^2} +
      \arcsin(h)\right)\right)
      1_{[0,1]}(h)}{
      C(h)=1-2/pi*(h sqrt(1-h^2)+asin(h))   if
      0<=h<=1, 0 otherwise}
    This isotropic covariance function is valid only for dimensions
    less than or equal to 2.

   
    \item \code{cone}\cr
    This model is used only for methods based on marked point processes
    (see \command{\link{RFMethods}}); it is defined only in two dimensions.
    The corresponding (boolean)
    function is a truncated cone with socle. The base has radius
    \eqn{\frac12}{1/2}. The model has three parameters, \eqn{r},
    \eqn{s}, and \eqn{h}:\cr
    \eqn{r} gives the radius of the top circle of the cone, given
    as part of the socle radius; \eqn{r \in [0,1)}.\cr
    \eqn{s} gives the height of the socle.\cr
    \eqn{h} gives the height of the truncated cone.\cr
       
    \item \code{coxisham} see \link{sophisticated}.
    
    \item \code{cutoff} see \link{sophisticated}.

    
    \item \code{cubic}
    \deqn{C(h)=(1- 7h^2+8.75h^3-3.5h^5+0.75 h^7)1_{[0,1]}(h)}{C(h)=
      1- 7 h^2 + 8.75 h^3 - 3.5 h^5 + 0.75 h^7  if 0<=h<=1,
      0 otherwise}
    This model is valid only for dimensions less than or equal to 3. 
    It is a 2 times differentiable covariance functions with compact
    support. %(See Chiles&Delfiner, 1998)
        
    \item \code{dampedcosine} (hole effect model)
    \deqn{C(h)= e^{-\lambda h} \cos(h),   \quad h\ge0}
    This model is valid
    for dimension 1 iff \eqn{\lambda \ge 1},
    for dimension 2 iff \eqn{\lambda \ge 1},
    and for dimension 3 iff \eqn{\lambda \ge \sqrt{3}}.

    \item \code{DeWijsian}
    \deqn{\gamma(h) = \log(\|h\|^\alpha + 1)}
    generalised version of the DeWijsian model with \eqn{\alpha \in (0,2]}

 
    \item \code{EAxxA} and  see `\link{sophisticated}'

    \item \code{EtAxxA} and  see `\link{sophisticated}'
    
    
    \item \code{exponential (normal scale mixture)}
    \deqn{C(h)=e^{-h},  \quad h\ge0}{C(h)=exp(-h)}
    This model is a special case of the \code{whittle} model
    (for \eqn{\nu=\frac12}{nu = 1/2} there)
    and the \code{stable} class (for \eqn{\alpha = 1}).
    
     
    \item \code{FD}
    \deqn{C(k) = \frac{(-1)^k \Gamma(1-a/2)^2}{\Gamma(1-a/2+k)
	\Gamma(1-a/2-k), 
	\qquad k \in {\bf N}}
    }{C(k) = {(-1)^k \Gamma(1-a/2)^2} / {\Gamma(1-a/2+k)
	\Gamma(1-a/2-k), 
	k in  N}}
    and linearly interpolated otherwise.
    Here, \eqn{\Gamma} is the Gamma function and \eqn{a \in [-1, 1)}. 
    The model is defined in 1 dimension only.
    
    Remark: the fractionally differenced process
    stems from time series modelling
    where the grid locations are multiples
    of the scale parameter.  
    
    \item \code{fractalB} (fractal Brownian motion)
    \deqn{gamma(h) = h^\alpha}   
    Here, \eqn{\alpha \in (0,2]}.
    (Implemented for up to three dimensions). See also \code{genB}.

    \item \code{fractgauss}
    \deqn{C(h) = 0.5 (|h+1|^{\alpha} - 2|h|^{\alpha} +
      |h-1|^{\alpha})}
    This model is the covariance function for the fractional Gaussian noise
    with Hurst parameter \eqn{H=\alpha /2}, \eqn{\alpha \in (0,2]}.
    In particular, the model is valid only
    in one dimension.
	
    \item \code{gauss (normal scale mixture)}
    \deqn{C(h)=e^{-h^2}}{C(h)=exp(-h^2)}
    This model is a special case of the \code{stable} class
    (for \eqn{\kappa=2}{a=2} there).
    Note that the corresponding function for the random coins
    method (cf. the methods based on marked point processes in
    \command{\link{RFMethods}}) is
    \deqn{e^{- 2 h^2}.}{exp(-2 h^2).}
    See \code{gneiting} for an alternative model that does not have
    the disadvantages of the Gaussian model.

    
    \item \code{genB} (generalised fractal Brownian motion)
    \deqn{\gamma(h) = (h^{\alpha}+1)^{\delta} -1 }
    Here, \eqn{\alpha \in(0,2]} and
    \eqn{\delta \in (0,1)}.
    (Implemented for up to three dimensions). See also \code{fractalB}.

    
    \item \code{gencauchy} (generalised \code{cauchy}; normal scale mixture)\cr
    \deqn{C(h)=
      \left(1+h^\alpha\right)^(-\beta/\alpha)}
    The parameter \eqn{\alpha} is in (0,2], and \eqn{\beta}
    is positive.
    \cr
    This model allows for simulating random fields where
    fractal dimension and Hurst coefficient can be chosen
    independently.

    \item \code{gengneiting} (generalised \code{gneiting})\cr
    If \eqn{n=1} then
    \deqn{C(h)=\left(1+(\alpha+1)h\right) * (1-h)^{\alpha+1}
      1_{[0,1]}(h)}
    If \eqn{n=2} then
    \deqn{C(h)=\left(1+(\alpha+2)h+\left((\alpha+2)^2-1\right)h^2/3\right)
      (1-h)^{\alpha+2}  1_{[0,1]}(h)}
    If \eqn{n=3} then
    \deqn{C(h)=\left(1+(\alpha+3)h+\left(2(\alpha+3)^2-3\right)h^2/5
      +\left((\alpha+3)^2-4\right)(\alpha+3)h^3/15\right)(1-h)^{\alpha+3}
      1_{[0,1]}(h)}
    The parameter \eqn{n} is a positive integer; here only the
    cases \eqn{n=1, 2, 3} are implemented.
    The parameter \eqn{\alpha} is greater than or equal to
    \eqn{(d + 2n +1)/2} where \eqn{d} is the
    dimension of the random field.

    % the differentiability is ??

    
    \item \code{gneiting}
    \deqn{C(h)=\left(1 + 8 sh + 25 (sh)^2 + 32
      (sh)^3\right)(1-sh)^8 1_{[0,1]}(sh)}{C(h)=
      (1 + 8 s h + 25 s^2 h^2 + 32
      s^3 h^3)*(1-s h)^8   if 0 <= h <= 1/s, 0 otherwise}
    where
    \eqn{s=0.301187465825}.
    This isotropic covariance function is valid only for dimensions less
    than or equal to 3. 
    It is a 6 times differentiable covariance functions with compact
    support.\cr
    It is an alternative to the \code{gaussian} model since
    its graph is visually hardly distinguishable from the graph of
    the Gaussian model, but possesses neither the mathematical and nor the
    numerical disadvantages of the Gaussian model.\cr
    This model is a special case of \code{gengneiting} (for
    \eqn{n=3} and \eqn{\alpha=5} there).
    Note that, in the original work by Gneiting (1999), 
    \eqn{s=\frac{10\sqrt2}{47}\approx 0.3008965}{s = 10 sqrt(2) / 47 ~=
      .3008965}, a numerical value slightly deviating from the
    optimal one.
    

    \item gneitingdiff is obsolete, see the last example in
    \link{Sophisticated} for a user's definition of \code{gneitingdiff}.
    
    \deqn{C(h)=( 1 + 8 h \alpha^{-1}
      + 25 h^2\alpha^{-2}
      + 32 h^3 \alpha^{-3} )
      ( 1-h \alpha^{-1} )^8 
      2^{1-\nu} (\Gamma(\nu))^{-1}
      h^{\nu} K_{\nu}(h) 1_{[0,\alpha]}(h)}
    This isotropic covariance function is valid only for dimensions less
    than or equal to 3. 
    The parameters \eqn{\nu} and \eqn{\alpha} are
    positive.\cr
    This class of models with compact support
    allows for smooth parametrisation of the differentiability up to
    order 6.     
       
    \item \code{hyperbolic (normal scale mixture)}
    \deqn{C(h)= \delta^{-\lambda}
	(K_{\lambda}(\nu \delta))^{-1}
      ( \delta^2 + h^2 )^{\lambda/2}
      K_{\lambda}(
      \nu [ \delta^2 + h^2 ]^{1/2} )}
    The parameters are such that\cr
    \eqn{\delta\ge0},  \eqn{\nu>0}  and
    \eqn{\lambda>0,\quad}{\lambda>0, }
    or\cr
    \eqn{\delta>0},  \eqn{\nu>0}  and  \eqn{\lambda=0,\quad}{\lambda>0, }
    or\cr
    \eqn{\delta>0},  \eqn{\nu\ge0}, and   \eqn{\lambda<0}.\cr
    Note that this class is over-parametrised; always one
    of the three parameters 
    \eqn{\nu}, \eqn{\delta}, and scale
    can be eliminated in the formula. Therefore, one of these
    parameters should be kept fixed in any simulation study.
    \cr
    The model contains as special cases the \code{whittle}
    model and the \code{cauchy} model, for 
    \eqn{\delta=0} and \eqn{\nu=0}, respectively.

    See also \code{nonstathyperbolic} in \command{\link{Covariance}}.

    
    \item \code{iacocesare} (non-separabel space time model)
    \deqn{C(h, t)=(1+\|h\|^\nu+|t|^\lambda)^{-\delta}}
    The parameters \eqn{\nu} and \eqn{\lambda} take values
    in \eqn{[1,2]}; the parameters \eqn{\delta} must be greater
    than or equal to half the space-time dimension.
    
   
    \item J-Bessel\cr
    see \code{bessel}

    \item K-Bessel\cr
    see \code{whittle} and \code{matern}
    
    \item linear with sill\cr
    See \code{power} (\code{a=1} there).

    \item \code{lgd1} (local-global distinguisher)
    \deqn{C(h)=
      1-\frac\beta{\alpha+\beta}|h|^{\alpha}, |h|\le 1 \qquad \hbox{and} \qquad
      \frac\alpha{\alpha+\beta}|h|^{-\beta}, |h|> 1
    }{C(h)=
      1-\beta^{-1}{\alpha+\beta}|h|^{\alpha}, |h|\le 1    and   
      \alpha^{-1}{\alpha+\beta}|h|^{-\beta}, |h|> 1
    }
    Here \eqn{\beta>0} and \eqn{\alpha} is in
    \eqn{(0,(3 - d)/2]} for dimension \eqn{d=1,2}.
    The random field has fractal dimension
    \eqn{d + 1 - \alpha/2}
    and Hurst coefficient \eqn{1 -\beta/2} for
    \eqn{\beta \in (0,1]}


%    \item \code{matern}\cr
%    equals \eqn{W_\nu(\sqrt{2\nu}x)} where \eqn{W_\nu} is the
    
    \item matern (normal scale mixture)\cr
    \deqn{C(x)=W_a(x) = 2^{1-\nu} \Gamma(\nu)^{-1}
       (\sqrt{2 \nu} x)^\nu K_\nu(\sqrt{2 \nu}x)}
     The parameter \eqn{\nu} is positive.
    \cr
    This is the model of choice if the smoothness of a random field is to
    be parametrised: if \eqn{\nu > m} then the
    graph is \eqn{m} times differentiable.

    In contrast to the \code{whittle} model
    this model separates the effects of the scaling parameter and the
    shape parameter. For \eqn{\nu=0.5} we get the exponential
    model; for \eqn{\nu=\infty} we get \eqn{C(x) = \exp(0.5 x^2)}.
    
    The model \eqn{C(x \sqrt{2})} equals the Handcock-Wallis (1994)
    parameterisation.

    The model allows further to replace \eqn{nu} by \eqn{1/\nu},
    setting the second parameter \code{invnu=TRUE}.

    See also \code{whittle}, and
    \code{nonstatwhittle} in \command{\link{Covariance}}.


    \item \code{M} and  see `\link{sophisticated}'
    \item \code{mastein} see `\link{sophisticated}'
    \item \code{mixed} see `\link{sophisticated}'
    
   
    \item \code{nugget}
    \deqn{C(h)=1_{\{0\}}(h)}{C(h)=1_{0}(h)}
    If the model is used in \code{param}-definition mode,
    either \code{param[2]}, the \code{variance},
    or \code{param[3]}, the \code{nugget}, must be zero.
    If the model is used in the list-definition mode,
    the anisotropy matrix must be given in an anisotropic
    context, but not
    the scale parameter in an isotropic context.
    See also \command{\link{sophisticated}}.
     
 
    \item \code{penta}
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
    \deqn{C(x)= (1-x)^a 1_{[0,1]}(x)}
    This covariance function is valid for dimension \eqn{d}{d} if
    \eqn{a \ge (d+1)/2}. 
    For \eqn{\kappa=1}{a=1} we get the well-known triangle (or tent)
    model, which is valid on the real line, only.
    % proposition 3.8 in phd thesis tilmann gneiting
    % Golubov, Zastavnyi

    
    \item powered exponential\cr
    See \code{stable}.

    
    \item \code{qexponential}
    \deqn{C(x)= ( 2 e^{-x} - \alpha e^{-2x} ) / ( 2 - \alpha )}
    The parameter \eqn{\alpha} takes values in \eqn{[0,1]}.   
    
    
    % \item rational quadratic model\cr
    %  See \code{cauchy} for \eqn{\kappa=1}{a=1}.
    % (Cressie)
    \item \code{rational} and  see `\link{sophisticated}'

    
    \item \code{spherical}
    \deqn{C(x)=\left(1- 1.5 x+0.5 x^3\right)
      1_{[0,1]}(x)}
    %{C(x)=
%      1 - 1.5 x + 0.5 x^3  if 0<=x<=1, 0 otherwise}
    This isotropic covariance function is valid only for dimensions
    less than or equal to 3.

    
    \item \code{stable}
    \deqn{C(x)=\exp\left(-x^\alpha\right)}
    The parameter \eqn{\alpha} is in \eqn{(0,2]}.
    See \code{exponential} and \code{gaussian} for special cases.

    \item \code{Stein} and  see `\link{sophisticated}'
    \item \code{steinst1} and  see `\link{sophisticated}'
   
    \item symmetric stable\cr
    See \code{stable}.
    
    \item \code{tbm2} and  see `\link{sophisticated}'
    \item \code{tbm3} and  see `\link{sophisticated}'

    \item tent model\cr
    See \code{power}.
    
    \item triangle\cr
    See \code{power}.
    
    \item \code{wave}
    \deqn{C(x)=\frac{\sin x}x, \quad x>0 \qquad \hbox{and } C(0)=1}{
      C(x)=sin(x)/x if x>0 and C(0)=1}
    This isotropic covariance function is valid only for dimensions less
    than or equal to 3.
    It is a special case of the \code{bessel} model
    (for \eqn{\kappa}{a}\eqn{=0.5}).
    
    \item \code{whittle (normal scale mixture)}
    \deqn{C(x)=W_\nu(x) = 2^{1-\nu} \Gamma(\nu)^{-1} x^\nu
      K_\nu(x)}
    The parameter \eqn{\nu} is positive.
    \cr
    This is the model of choice if the smoothness of a random field is to
    be parametrised: if \eqn{\nu > m} then the
    graph is \eqn{m} times differentiable.

    The model is a special case of the
    \code{hyperbolic} model (for \eqn{\nu_3=0}{c=0} there).

    See also \code{nonstWM} in \link{sophisticated}.
   }
  
   Let \eqn{\code{cov}} be a model given in standard notation. 
   Then the covariance model
   applied with arbitrary variance and scale equals
  \deqn{\code{variance} * \code{cov}( (\cdot)/ \code{scale}).
  }{variance * cov( (.)/scale).}
  The parameters can be passed by the vector \code{param},
  \code{param=c(mean, variance, nugget, scale, ...)}.
  Here \sQuote{...} stands for additional parameters such as \eqn{\nu}
  in the \code{whittle} model.
  In case a model has several parameters, as in \code{hyperbolic},
  the parameters must be given in the sequence they are explained
  aboved. However, it is strongly recommended to use the \code{list}
  notation explained in  \code{sophisticated}. The \code{list}
  definition available in \pkg{RandomFields} V 1.x, is depreciated!
  
  For a given covariance function \eqn{cov}{cov} the variogram
  \eqn{\gamma} equals 
  \deqn{\gamma(x) = cov(0) -  cov(x).}
  
  Note:
  \itemize{
    \item  The value of the covariance function or variogram
  depends also on
  \command{\link{RFparameters}}\code{()$PracticalRange}.  If the latter is
  \code{TRUE} and the covariance model is isotropic
  then the covariance function is internally
  rescaled such that cov\eqn{(1)\approx 0.05}{(1)~=0.05} for standard
  parameters (\code{scale=1}).

  \item   Some models allow certain parameter combinations only for certain
    dimensions.  As any model valid in \eqn{d}{d} dimensions is also valid in 1
    dimension, the default in \code{CovarianceFct} and \code{Variogram}
    is \code{dim=1}.
  }
}
\value{
  \code{CovarianceFct} returns a vector of values of the covariance
  function.
  
  \code{Variogram} returns a vector of values of the variogram model.

  \code{CovMatrix} return a covariance matrix. Here a matrix of
  of coordinates (\code{x}) or a vector or a matrix of \code{Distances}
  is expected.
  \code{CovMatrix} allows also for variogram models. Then negative of
  variogram matrix is returned.
}
\references{
  Overviews:
  \itemize{
    \item Chiles, J.-P. and Delfiner, P. (1999)
    \emph{Geostatistics. Modeling Spatial Uncertainty.}
    New York: Wiley.
    
    \item  Gneiting, T. and Schlather, M. (2004)
    Statistical modeling with covariance functions.
    \emph{In preparation.}
    
    \item Schlather, M. (1999) \emph{An introduction to positive definite
      functions and to unconditional simulation of random fields.}
    Technical report ST 99-10, Dept. of Maths and Statistics,
    Lancaster University.

    \item Schlather, M. (2002) Models for stationary max-stable
      random fields. \emph{Extremes} \bold{5}, 33-44.
    
    \item Yaglom, A.M. (1987) \emph{Correlation Theory of Stationary and
      Related Random Functions I, Basic Results.}
    New York: Springer.
    
    \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
    Springer, 3nd edition.
  }

  Cauchy models, generalisations and extensions
  \itemize{
    \item Gneiting, T. and Schlather, M. (2004)
    Stochastic models which separate fractal dimension and Hurst effect.
    \emph{SIAM review} \bold{46}, 269-282.% see also lgd
  }

  Dagum model
  \itemize{
    \item Porcu, E., Zini, A. and Pini, R. (2007)
    \emph{Stats. Probab. Lett.}, \bold{77}, 83-89. 

    \item Berg, C., Mateu, J. and Porcu, E. (2008)
    \emph{Bernoulli}, \bold{14}, 1134-1149.
  }

  Generalised fractal Brownian motion
   \itemize{
   \item Gneiting, T. (2002) Nonseparable, stationary covariance
    functions for space-time data, \emph{JASA} \bold{97}, 590-600.
  }
  
  Gneiting's models
  \itemize{
    \item Gneiting, T. (1999)
    Correlation functions for atmospheric data analysis.
    \emph{Q. J. Roy. Meteor. Soc., Part A} \bold{125}, 2449-2464. 
  }
   
  Holeeffect model
  \itemize{
    \item Zastavnyi, V.P. (1993)
    Positive definite functions depending on a norm.
    \emph{Russian Acad. Sci. Dokl. Math.} \bold{46}, 112-114. 
  }

  Hyperbolic model
  \itemize{
    \item Shkarofsky, I.P. (1968) Generalized turbulence space-correlation and
    wave-number spectrum-function pairs. \emph{Can. J. Phys.} \bold{46},
    2133-2153.
  }
  
  fractalB
  \itemize{
    \item  Stein, M.L. (2002)
    Fast and exact simulation of fractional Brownian surfaces.
    \emph{J. Comput. Graph. Statist.} \bold{11}, 587-599.
  }
  
  
  genB
  \itemize{
    \item Schlather, M. (2010)
    On some covariance models based on normal scale mixtures.
    \emph{Bernoulli}, \bold{16}, 780-797.
  }

  
  lgd
  \itemize{
    \item Gneiting, T. and Schlather, M. (2004)
    Stochastic models which separate fractal dimension and Hurst effect.
    \emph{SIAM review} % see also cauchy
  }
   
  Power model
  \itemize{
    \item Golubov, B.I. (1981) On Abel-Poisson type and Riesz means,
    \emph{Analysis Mathematica} \bold{7}, 161-184.
    
    \item  Zastavnyi, V.P. (2000) On positive definiteness of some
    functions, \emph{J. Multiv. Analys.} \bold{73}, 55-81.
  }

}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}
}
\seealso{
  \command{\link{sophisticated}},
  \command{\link{EmpiricalVariogram}},
  \command{\link{GetModel}},
  \command{\link{GetPracticalRange}},
  \command{\link{parameter.range}},
  \code{\link{RandomFields}},
  \command{\link{RFparameters}},
  \command{\link{ShowModels}}.}

\examples{
% library(RandomFields, lib="~/TMP")

PrintModelList()
x <- 0:100

## the following five model definitions are the same!
## 
## (1) very traditional form
(cv <- CovarianceFct(x, model="bessel", param=c(NA,2,1,5,0.5)))
plot(x, cv)


## (2) above model in the very general list definition
model <- list("+",
              list("$", var=2, scale=5, list("bessel", 0.5)),
              list("nugget"))
cv <- CovarianceFct(x, model=model)
points(x, cv, col="red", pch=20) ## no differnce to first


% source("/home/schlather/R/RF/RandomFields/tests/source.R")
## (3) nested model definition
## this kind of definiton models is depreciated from Version 2.0 on
cv <- CovarianceFct(x, model="bessel",
                  param=rbind(c(2, 5, 0.5), c(1, 0, 0)))
points(x, cv, col="blue", pch=20, cex=0.5) 
 

## (4) anisotropic notation
 model <- list("+",
               list("$", var=2, aniso=as.matrix(0.2),
                    list("bessel", nu=0.5)
                   ),
               list("nugget")
              ) 
cv <- CovarianceFct(as.matrix(x), model=model)
points(x, cv, col="green", pch=4) 
 


## Depreciated list defintions in Version 1.x
## this way of defining a model still works, but
## is not supported anymore
## (isotropic version)
model <- list(list(model="bessel", var=2, kappa=0.5, scale=5),
              "+",
              list(model="nugget", var=1, scale=1))
cv <- CovarianceFct(x, model=model)
points(x, cv, col="black", pch=5) 

}
\keyword{spatial}






