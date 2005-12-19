\name{CovarianceFct}
\alias{CovarianceFct}
\alias{Variogram}
\alias{PrintModelList}
\alias{GetModelList}
\alias{GetModelNames}
\title{Covariance And Variogram Models}
\description{
  \code{CovarianceFct} returns the values of a covariance function
  
  \code{Variogram} returns the values of a variogram model

  \code{PrintModelList} prints the list of currently implemented models
  including the corresponding simulation methods

  \code{GetModelList} returns a matrix of currently implemented models
  and their simulation methods

  \code{GetModelNames} returns a list of currently implemented models
 }
\usage{
CovarianceFct(x, model, param, dim=ifelse(is.matrix(x),ncol(x),1),
              fctcall="Covariance")

Variogram(x, model, param, dim=ifelse(is.matrix(x),ncol(x),1))

PrintModelList()

GetModelList(abbr=TRUE)

GetModelNames()
}
\arguments{
  \item{x}{vector or \eqn{(n \times \code{dim})}{(n x
      \code{dim})}-matrix.  In particular,
    if the model is isotropic or \code{dim=1} then \code{x}
    is a vector.}

  \item{model}{character or list;
    if character then name of the covariance function or
    variogram model - see below, or type \code{PrintModelList()} for
    all options; see Details for the definition of the model by a list. 
  }
  \item{param}{vector or matrix of parameters or missing, see Details
    and Examples; 
    The simplest form is that \code{param} is vector of the form
    \code{param=c(NA,variance,nugget,scale,...)}, in this order;\cr
    The dots \code{...} stand for additional parameters of the
    model. 
  }
  \item{dim}{dimension of the space in which the model is applied}
  \item{fctcall}{this parameter should not be changed by the user}
  \item{abbr}{logical or numerical. If \code{TRUE} the names for the methods are
    abbreviated. If numerical, \code{abbr} gives the number of letters.} 
}
\details{  
  The implemented models are in standard notation for a
  covariance function (variance 1, nugget 0, scale 1) and for positive
  real arguments \eqn{x} (and \eqn{t}): 
  
  \itemize{
    \item \code{bessel}
    \deqn{C(x)=2^a \Gamma(a+1)x^{-a} J_a(x)}{C(x)=
      2^a Gamma(a+1)x^(-a) J_a(x)}
    The parameter \eqn{\kappa}{a} is greater than or equal to
    \eqn{\frac{d-2}2}{(d-2)/2}, where \eqn{d}{d} is the
    dimension of the random field.

    \item Brownian motion\cr
    see \code{fractalB} 
    
    \item cardinal sine\cr
    see \code{wave}
    
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
    The parameter \eqn{\kappa_1}{a} is in (0,2] and  \eqn{\kappa_2}{b}
    is positive.
    The model is valid for dimensions \eqn{d\le\kappa_3}{d<=c};
    this has been shown for integer \eqn{\kappa_3}{c}, but the
    package allows real values of \eqn{\kappa_3}{c}.
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

    \item \code{constant}\cr
    Identically constant. Any scale parameter is ignored.
    
    \item \code{cone}\cr
    This model is used only for methods based on marked point processes
    (see \command{\link{RFMethods}}); it is defined only in two dimensions.
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
    
    \item \code{cutoff} (hypermodel, see also below)
    \deqn{C(x)=\phi(xB), 0\le xA \le \kappa_2}{C(x) = phi(xB),
      0 <= xA <= b}
    
    \deqn{C(x)=b_1 (r^{\kappa_3} - xB)^{2 \kappa_3}, \kappa_2\le xA \le
      \kappa_2\kappa_3}{C(x) = b_0 (r - xB)^{2 c},   b <= xA <= bc}
    
    \deqn{C(x)=0,  \kappa_2\kappa_3\le xA}{C(x) = 0,   bc <= xA}
   
    The cutoff model is a functional of the covariance function
    \eqn{\phi}{phi}. Since the model itself is indifferent for
    scale or anisotropy parameters, the latter must be given only
    for the submodels.
    See below for general comments on hypermodels.
    
    The first parameter, \eqn{\kappa_1}{a}, gives the number of
    subsequent models that build \eqn{\phi}{phi};
    \eqn{\kappa_2\ge0}{b>0}, \eqn{\kappa_3>0}{c>0}.
    The parameters \eqn{r} and \eqn{b_0}
    are chosen internally such that \eqn{C} is a smooth function.
    The parameters \eqn{A} and \eqn{B} are the inverse scale parameters for
    the hypermodel and submodel, respectively. Note that \code{cutoff} seldemly
    works, if \eqn{A} and \eqn{B} are not identical.

     The algorithm that checks the given parameters knows
     only about some few necessary conditions.
     Hence it is not ensured that
     the cutoff-model is a valid covariance function for any
    choice of phi and the parameters. 
    
    For certain models \eqn{\phi}{phi}, i.e. \code{stable},
    \code{whittle} and \code{gencauchy}, some sufficient conditions
    are known.
    
    \item \code{dampedcosine}
    \deqn{C(x)=e^{-\kappa x} \cos(x),  \quad x\ge0}{C(x)=
      exp(-a x) cos(x) }
    This model is valid
    for dimension 1 iff \eqn{\kappa\ge1}{a>=0},
    for dimension 2 iff \eqn{\kappa\ge1}{a>=1},
    and for dimension 3 iff \eqn{\kappa\ge \sqrt{3}}{a >= sqrt(3)}.

    \item \code{exponential}
    \deqn{C(x)=e^{-x},  \quad x\ge0}{C(x)=exp(-x)}
    This model is a special case of the \code{whittlematern} model
    (for \eqn{\kappa=\frac12}{a=1/2} there)
    and the \code{stable} class (for \eqn{\kappa=1}{a=1}).
    
    \item \code{fractalB} (fractal Brownian motion)
      \deqn{\gamma(x) = x^\kappa}{gamma(x) = x^a}   
      The parameter \eqn{\kappa}{a} is in \eqn{(0,2]}.
      (Implemented for up to three dimensions)   
      
    \item \code{FD}
    \deqn{C(k) = \frac{(-1)^k \Gamma(1-\kappa/2)^2}{\Gamma(1-\kappa/2+k)
	\Gamma(1-\kappa/2-k), 
      \qquad k \in {\bf N}}}{C(k) = (-1)^k
      G(1-a/2)^2/(G(1-a/2+k)G(1-a/2-k)) for integer k}
    and linearly interpolated otherwise.
    Here, \eqn{\Gamma}{G} is the Gamma function.
    The parameter \eqn{\kappa}{a} is in
    \eqn{[-1, 1)}. 
    The model is defined in 1 dimension only.
    
    Remark: the fractionally differenced process
    stems from time series modelling
    where the grid locations are multiples
    of the scale parameter.  
    
    \item \code{fractgauss}
    \deqn{C(x) = 0.5 (|x+1|^{\kappa_1} - 2|x|^{\kappa_1} +
      |x-1|^{\kappa_1})}{C(x) = 0.5 (|x+1|^a - 2|x|^a + |x-1|^a)}
      This model is the covariance function for the fractional Gaussian noise
      with Hurst parameter \eqn{H=\kappa_1 /2}{H =a/2}, \eqn{\kappa_1 \in
	(0,2]}{a in (0,2]. In particular, the model is valid only
	in one dimension.
      }
	
    \item \code{gauss}
    \deqn{C(x)=e^{-x^2}}{C(x)=exp(-x^2)}
    This model is a special case of the \code{stable} class
    (for \eqn{\kappa=2}{a=2} there).
    Note that the corresponding function for the random coins
    method (cf. the methods based on marked point processes in
    \command{\link{RFMethods}}) is
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

    \item \code{gengneiting} (generalised \code{gneiting})\cr
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
    The parameter \eqn{\kappa_1}{a} is a positive integer; here only the
    cases \eqn{\kappa_1=1, 2, 3}{a=1, 2, 3} are implemented.
    The parameter \eqn{\kappa_2}{b} is greater than or equal to
    \eqn{(d + 2\kappa_1 +1)/2}{(d + 2a +1)/2} where \eqn{d} is the
    dimension of the random field.

    % the differentiability is ??
    
    \item \code{gneiting}
    \deqn{C(x)=\left(1 + 8 sx + 25 (sx)^2 + 32
      (sx)^3\right)(1-sx)^8 1_{[0,1]}(sx)}{C(x)=
      (1 + 8 s x + 25 s^2 x^2 + 32
      s^3 x^3)*(1-s x)^8   if 0<=s x<=1, 0 otherwise}
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
    \eqn{\kappa_1=3}{a=3} and \eqn{\kappa_2=5}{b=5} there).
    Note that, in the original work by Gneiting (1999), 
    \eqn{s=\frac{10\sqrt2}{47}\approx 0.3}{s = 10 sqrt(2) / 47 ~=
      .3008965}, a numerical value slightly deviating from the
    optimal one.
    

    \item gneitingdiff is obsolete, see example below
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
       
    \item \code{hyperbolic}
    \deqn{C(x)=\frac{1}{\kappa_3^{\kappa_2}
	K_{\kappa_2}(\kappa_1 \kappa_3)}
      \left(\kappa_3^2 +x^2\right)^{{\kappa_2}/2}
      K_{{\kappa_2}}\left(
      \kappa_1 \left(\kappa_3^2 + x^2\right)^{1/2}\right),  \quad
      x>0}{C(x)=
      c^(-b) (K_b(a*c))^(-1) * (c^2 +x^2)^(0.5 b) *
      K_b(a sqrt(c^2 + x^2))}
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
    can be eliminated in the formula. Therefore, one of these
    parameters should be kept fixed in any simulation study.
    \cr
    The model contains as special cases the \code{whittlematern}
    model and the \code{cauchy} model, for 
    \eqn{\kappa_3=0}{c=0} and \eqn{\kappa_1=0}{a=0}, respectively.

    \item \code{iacocesare}
    \deqn{C(x, t)=(1+\|x\|^{\kappa_1}+|t|^{\kappa_2})^{-\kappa_3}}{
      C(x, t) = (1+||x||^a+|t|^b)^{-c}}
    The parameters \eqn{\kappa_1}{a} and \eqn{\kappa_2}{b} take values
    in \eqn{[1,2]}; the parameters \eqn{\kappa_3}{c} must be greater
    than or equal to half the space-time dimension.
    
    \item J-Bessel\cr
    see \code{bessel}

    \item K-Bessel\cr
    see \code{whittlematern}
    
    \item linear with sill\cr
    See \code{power} (\code{a=1} there).

    \item \code{lgd1} (local-global distinguisher)
    \deqn{C(x)=
      1-\frac\beta{\alpha+\beta}|x|^{\alpha}, |x|\le 1 \qquad \hbox{and} \qquad
    \frac\alpha{\alpha+\beta}|x|^{-\beta}, |x|> 1
    }{
      C(x)=
      1- b(a+b)^{-1}|x|^a for |x|\le 1   and  
    a(a+b)^{-1}|x|^-b for |x|> 1
    }
    Here \eqn{\beta>0}{b>0} and \eqn{\alpha}{a} is in
    \eqn{(0, \frac12 (3 - d)]}{(0, 1.5-d/2]} for dimension \eqn{d=1,2}.
    The random field has fractal dimension
    \eqn{d + 1 - \frac\alpha2}{d + 1 - a/2}
    and Hurst coefficient \eqn{1 - \frac\beta2}{1 - b/2} for
    \eqn{\beta\in(0,1]}{b in (0,1]}

    \item \code{mastein} (hypermodel for non-separabel space time modelling)
    \deqn{C(x, t)=\frac{\Gamma(\kappa_2 + \gamma(t))\Gamma(\kappa 2 +
	\kappa_3)}{ \Gamma(\kappa_2 +
	\gamma(t) + \kappa_3) \Gamma(\kappa_2)}
      W_{\kappa_2 + \gamma(t)}(\|x - Vt\|)}{
      C(x, t)=Gamma(b + g(t)) Gamma(b + c) / [Gamma(b + g(t) + c)
      Gamma(b)]
      W_{b + g(t)}(\|x - Vt\|)}

    \eqn{\Gamma}{Gamma} is the Gamma function; \eqn{\gamma}{g}
    is a variogram; \eqn{W} is the Whittle-Matern model.
    The first parameter, \eqn{\kappa_1}{a}, gives the number of
    subsequent models that build \eqn{\gamma}{g}. Here, the
    names of covariance models can also be used; the algorithm
    chooses the corresponding variograms then.
    The parameter \eqn{\kappa}{b} is the smoothness parameter
    of the Whittle-Matern model (for \eqn{t=0}) and must be positive.
    Finally, \eqn{c} must be greater than or equal to half the
    dimension of \eqn{x}.
    Instead of the velocity parameter \eqn{V}, the anisotropy matrix
    for the hyper model is chosen appropriately. Note that
    the anisotropy matrix must be such that \eqn(x,t) is transformed
    into a purely spatial vector, i.e. the entries in
    last column of the matrix are all naught.
    On the other hand, all entries of the
    anisotropy matrices in the submodels that build \eqn{\gamma}{gamma}
    is naught except the very last, purely temporal one.
   
    
    \item matern\cr
    See \code{whittlematern}.
    
   \item \code{nsst} (Non-Separable Space-Time model)
    \deqn{C(x,t)= \psi(t)^{-\kappa_6} \phi(x / \psi(t))}{C(x,t)=
      psi(t)^{-f} \phi(x / psi(t))}
    This model is used for space-time modelling where the spatial
    component is isotropic.\cr
    \eqn{\phi} is the \code{stable} model if \eqn{\kappa_2=1}{b=1};\cr
    \eqn{\phi} is the \code{whittlematern} model if \eqn{\kappa_2=2}{b=2};\cr
    \eqn{\phi} is the \code{cauchy} model if \eqn{\kappa_2=3}{b=3};\cr
    Here, \eqn{kappa_1}{a} is the respective parameter for the model.

    The function \eqn{\psi}{psi} satisfies\cr
    \eqn{\psi^2(t) = (t^{\kappa_3}+1)^{\kappa_4}}{psi^2(t) =
      (t^c+1)^d} if \eqn{\kappa_5=1}{e=1}\cr
    \eqn{\psi^2(t) = \frac{\kappa_4^{-1}t^{\kappa_3}+1}{t^{\kappa_3}+1}
    }{psi^2(t) = (d^{-1} t^c+1)/(t^c+1)} if \eqn{\kappa_5=2}{e=2}\cr
    \eqn{\psi^2(t) = -\log(t^{\kappa_3}+\kappa_4^{-1})/
	\log\kappa_4}{psi^2(t) = -\log(t^c+1/d)/log d}
      if \eqn{\kappa_5=3}{e=3}\cr
    The parameter \eqn{\kappa_6}{f} must be greater than or equal to
    the genuine spatial dimension of the field. Furthermore, \eqn{\kappa_3\in
      (0,2]}{c in (0,2]}  and \eqn{\kappa_4\in
      (0,1)}{d in (0,1)}.
     The spatial dimension must be \code{>=1}.

    \item \code{nsst2}
    \deqn{C(x,t)= \psi(t)^{-\kappa_7} \phi(x /\psi(t))}{C(x,t)=
      psi(t)^{-g} \phi(x / psi(t))}
    This model is used for space-time modelling where the spatial
    component is isotropic.  Here\cr
    \eqn{\phi} is the \code{gencauchy} model
    if \eqn{\kappa_3=1}{c=1}.\cr
    The parameters \eqn{kappa_1}{a} and \eqn{kappa_2}{b}
    are the respective parameters for the model.
    The function \eqn{\psi}{psi} satisfies\cr
    \eqn{\psi^2(t) = (t^{\kappa_4}+1)^{\kappa_5}}{psi^2(t) =
      (t^d+1)^e} if \eqn{\kappa_6=1}{f=1}\cr
    \eqn{\psi^2(t) = \frac{\kappa_5^{-1}t^{\kappa_4}+1}{t^{\kappa_4}+1}
    }{psi^2(t) = (e^{-1} t^d+1)/(t^d+1)} if \eqn{\kappa_6=2}{f=2}\cr
    \eqn{\psi^2(t) = -\log(t^{\kappa_4}+\kappa_5^{-1})/
	\log\kappa_5}{psi^2(t) = -\log(t^d+1/e)/log e}
      if \eqn{\kappa_6=3}{f=3}\cr
    The parameter \eqn{\kappa_7}{g} must be greater than or equal to
    the genuine spatial dimension of the field. Furthermore, \eqn{\kappa_4\in
      (0,2]}{d in (0,2]}  and \eqn{\kappa_5\in
      (0,1]}{e in (0,1]}.
    Necessarily, \code{dim>=2}.
    The spatial dimension must be \code{>=1}.

    \item \code{nugget}
    \deqn{C(x)=1_{\{0\}}(x)}{1(x==0)}
    If the model is used in \code{param}-definition mode,
    either \code{param[2]}, the \code{variance},
    or \code{param[3]}, the \code{nugget}, must be zero.
    If the model is used in the list-definition mode,
    the anisotropy matrix must be given in an anisotropic
    context, but not
    the scale parameter in an isotropic context.
    
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
    The parameter \eqn{\kappa}{a} is in \eqn{(0,2]}{(0,2]}.
    See \code{exponential} and \code{gaussian} for special cases.

    
    \item \code{Stein} (hypermodel, see also below)
    \deqn{C(x)=a_0 + a_2 (xB)^2 + \phi(xB), 0\le x \le \kappa_2}{C(x) = a_0 + a_2
      (xB)^2 + phi(xB),   0 <= xA <= b}
    
    \deqn{C(x)=b_1 (\kappa_3 - xB)^3/(xB), \kappa_2\le xA \le
      \kappa_2\kappa_3}{C(x) = b_0 (c - xB)^3/(xB),   b <= xA <= bc}
    
    \deqn{C(x)=0,  \kappa_2\kappa_3\le x}{C(x) = 0,   bc <= xA}
   
    The Stein model is a functional of the covariance function
    \eqn{\phi}{phi}. Since the model itself is indifferent for
    scale or anisotropy parameters, the latter must be given only
    for the submodels.
    See below for  general comments on hypermodels.
    
    The first parameter, \eqn{\kappa_1}{a}, gives the number of
    subsequent models that build \eqn{\phi}{phi};
    \eqn{\kappa_2\ge0}{b>0}, \eqn{\kappa_3\ge1}{c>=1}.
    The parameters \eqn{a_0}, \eqn{a_2} and \eqn{b_0}
    are chosen internally such that \eqn{C} becomes a smooth function.
    The parameters \eqn{A} and \eqn{B} are the inverse scale parameters for
    the hypermodel and submodel, respectively. Note that \code{Stein} seldemly
    works, if \eqn{A} and \eqn{B} are not identical.

    The algorithm that checks the given parameters knows
    only about some few necessary conditions.
    Hence it is not ensured that
    the Stein-model is a valid covariance function for any
    choice of phi and the parameters. 
    
    For certain models \eqn{\phi}{phi}, i.e. \code{stable},
    \code{whittle}, \code{gencauchy}, and the variogram
    model \code{fractalB}
    some sufficient conditions are known.
    
    
    \item \code{steinst1} (non-separabel space time model)
    \deqn{C(x, t) = W_{\kappa_1}(y) -
      \frac{\langle x, z \rangle t}{(\kappa_1 - 1)(2\kappa_1
	+ \kappa_2)} W_{\kappa_1 -1}(y)}{
      C(x, t) = W_a(y) - <x, z> t W_{a-1}(y) / [(a - 1)(2a + b)]}

    Here, \eqn{W_{\kappa_1}}{W_a} is the Whittle-Matern model with
    smoothness parameter \eqn{\kappa_1}{a};
    \eqn{\kappa_2}{b} is greater than or equal to the
    space-time dimension \sQuote{dim};
    \eqn{y=\|(x,t)\|}{y = ||(x,t)||}.
    The components of \eqn{z} are given by \eqn{\kappa_3,
      \ldots \kappa_{1+\sQuote{dim}}}{c, d, ...}; the norm of \eqn{z}
    must less than or equal to 1.
    
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
    (for \eqn{\kappa}{a}\eqn{=0.5}).
    
    \item \code{whittlematern}
    \deqn{C(x)=W_a(x) = 2^{1-\kappa} \Gamma(\kappa)^{-1} x^\kappa
      K_\kappa(x)}{C(x)=W_a(x) =2^{1-a} Gamma(a)^{-1} x^a K_a(x),
    }
    The parameter \eqn{\kappa}{a} is positive.
    \cr
    This is the model of choice if the smoothness of a random field is to
    be parametrised: if \eqn{\kappa\ge}{a>=}\eqn{(2m+1)/2} then the
    graph is \eqn{m} times differentiable.

    The model is a special case of the
    \code{hyperbolic} model (for \eqn{\kappa_3=0}{c=0} there).
    
   }
  
  Let \eqn{\rm cov}{cov} be a model given in standard notation. 
  Then the covariance model
  applied with arbitrary variance and scale equals
  \deqn{\rm \qquad variance * \rm cov( (\cdot)/ scale).
  }{variance * cov( (.)/scale).}

  For a given covariance function \eqn{\rm cov}{cov} the variogram
  \eqn{\gamma}{gamma} equals 
  \deqn{\gamma(x) = {\rm cov}(0) - {\rm cov}(x).}{
    gamma(x) = cov(0) - cov(x).}
  
  Note that the value of the covariance function or variogram
  depends also on
  \command{\link{RFparameters}}\code{()$PracticalRange}.  If the latter is
  \code{TRUE} and the covariance model is isotropic
  then the covariance function is internally
  rescaled such that cov\eqn{(1)\approx 0.05}{(1)~=0.05} for standard
  parameters (\code{scale==1}).

  The model and the parameters can be specified by three different
  forms; the first \sQuote{standard} form allows for the specification of the
  covariance model as given above for an isotropic random
  field.  The second form defines isotropic nested models using
  matrices.  The third
  form allows for defining anisotropic and/or space-time models
  using lists;
  here any basic models can arbitrarily be combined by multiplication and
  summation.
  \itemize{
    \item \code{model} is a string; \code{param} is a vector of the form
    \code{param=c(mean,variance,nugget,scale,...)}. (These components
    might be given separately or bound to a simple list passed to
    \code{model}.)
    
    The first component of \emph{param} is reserved for the \code{mean}
    of a random field and thus ignored in the evaluation of the covariance
    function or variogram.  The parameters mean, variance, nugget, and scale
    must be given in this order; additional
    parameters have to be supplied in case of a parametrised class of
    models (e.g. \code{hyperbolic}, see below),
    in the order \eqn{\kappa_1}{a}, \eqn{\kappa_2}{b}, \eqn{\kappa_3}{c}.
    
    Let \eqn{\rm cov}{cov} be a model given in standard notation. 
    Then the covariance model
    applied with arbitrary variance,
    nugget, and scale equals
    \deqn{\rm \qquad nugget + variance * \rm cov( (\cdot)/ scale).
    }{nugget + variance * cov( (.)/scale).}
     
    Some models allow certain parameter combinations only for certain
    dimensions.  As any model valid in \eqn{d}{d} dimensions is also valid in 1
    dimension, the default in \code{CovarianceFct} and \code{Variogram}
    is \code{dim=1}.

    \item \code{model} is a string; \code{param} is a matrix with columns
    of the form \code{c(variance, scale, ...)}.

    Except that the entries for the \code{mean} and the \code{nugget}
    are missing all explanations given above also apply here. 
    Each column defines a summand of the nested model. A nugget effect
    is indicated by \code{scale=0}; possibly additional parameters
    are ignored.
    
    \item \code{model} is a list as specified below; \code{param} is
    missing.

    \code{model = list(l.1, OP.1, l.2, OP.2, ..., l.n)}
    where \eqn{n} is at most 10 (except \code{cutoff embedding} is
    used, see \link{RFMethods}). The lists \code{l.i}
    are all either of the  form \code{l.i =
      list(model=,var=,kappas=,scale=,method=)} 
    or of the form  \code{l.i = list(model=,var=,kappas=,aniso=,method=)}.
    \code{model} is a string; \code{var} gives the variance;
    \code{scale} is a scalar whereas \code{aniso} is a \eqn{d \times
      d}{d x d} matrix, which is multiplied from the right to the
    \eqn{(n\times d)}{n x d} matrix of points;
    at the transformed points the values of the (isotropic)
    random field (with scale 1) are
    calculated. The dimension \eqn{d} of matrix must match the
    number of columns of \code{x}. The models given by \code{l.i}
    can be combined by \code{OP.i="+"} or \code{OP.i="*"}. 
    \code{method} is ignored here; it can be set in
    \command{\link{GaussRF}}.

    \item Hypermodels\cr
    hypermodels are functions or functionals of covariance functions
    or variograms. The first parameter is always the number of the following
    covariance models included.
    
%    The hypermodel inherits the anisotropy parameters (or the scale
%    parameter) from the first submodel. The given anisotropy parameters
%    are ignored.
    
    \bold{Important!} Hyper models are in an experimental stage:\cr
    (i) the (current) algorithm does not allow for a complete
    check whether the parameters for a hypermodel are well chosen.
    So, only use parameter combinations for which you are sure they lead to a
    positive definite function.
    \cr
    (ii) behaviour and parameters may change in future version!
  }
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
    
    \item Wackernagel, H. (1998) \emph{Multivariate Geostatistics.} Berlin:
    Springer, 2nd edition.
  }

  Cauchy models, generalisations and extensions
  \itemize{
    \item Gneiting, T. and Schlather, M. (2004)
    Stochastic models which separate fractal dimension and Hurst effect.
    \emph{SIAM review} \bold{46}, 269-282.% see also lgd
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

  Iaco-Cesare model
 \itemize{
   \item de Cesare, L., Myers, D.E., and Posa, D. (2002)
   FORTRAN programs for space-time modeling. \emph{Computers \&
     Geosciences} \bold{28}, 205-212.
   \item  de Iaco, S.. Myers, D.E., and Posa, D. (2002)
   Nonseparable space-time covariance models: some parameteric
   families. \emph{Math. Geol.} \bold{34}, 23-42.
  }
  
  lgd
  \itemize{
    \item Gneiting, T. and Schlather, M. (2004)
    Stochastic models which separate fractal dimension and Hurst effect.
    \emph{SIAM review} % see also cauchy
  }

  Ma-Stein model
  \itemize{
    \item Ma, C. (2003)
    Spatio-temporal covariance functions generated by mixtures.
    \emph{Math. Geol.}, \bold{34}, 965-975.
    \item Stein, M.L. (2005) Space-time covariance functions. \emph{JASA},
    \bold{100}, 310-321.
  }
  
  nsst
  \itemize{
    \item Gneiting, T. (2001) Nonseparable, stationary covariance
    functions for space-time data, \emph{JASA} \bold{97}, 590-600.
    
    \item Gneiting, T. and  Schlather, M. (2001)
    Space-time covariance models.
    In El-Shaarawi, A.H. and Piegorsch, W.W.:
    \emph{The Encyclopedia of Environmetrics.} Chichester: Wiley.
    }

  Power model
  \itemize{
    \item Golubov, B.I. (1981) On Abel-Poisson type and Riesz means,
    \emph{Analysis Mathematica} \bold{7}, 161-184.
    
    \item  Zastavnyi, V.P. (2000) On positive definiteness of some
    functions, \emph{J. Multiv. Analys.} \bold{73}, 55-81.
  }

  fractalB
  \itemize{
    \item  Stein, M.L. (2002)
    Fast and exact simulation of fractional Brownian surfaces.
    {\emph J. Comput. Graph. Statist.} \bold{11}, 587-599.
  }
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html};
  
  Yindeng Jiang \email{jiangyindeng@gmail.com} (circulant embedding
  methods \sQuote{cutoff} and \sQuote{intrinsic})}
\seealso{\command{\link{EmpiricalVariogram}},
  \command{\link{GetPracticalRange}},
  \command{\link{parameter.range}},
  \code{\link{RandomFields}},
  \command{\link{RFparameters}},
  \command{\link{ShowModels}}.}

\examples{
 PrintModelList()
 x <- 0:100

 # the following five model definitions are the same!
 ## (1) very traditional form
 (cv <- CovarianceFct(x, model="bessel", c(NA,2,1,5,0.5)))

 ## (2) traditional form in list notation
 model <- list(model="bessel", param=c(NA,2,1,5,0.5))
 cv - CovarianceFct(x, model=model)

 ## (3) nested model definition
 cv - CovarianceFct(x, model="bessel",
                    param=cbind(c(2, 5, 0.5), c(1, 0, 0)))

 #### most general notation in form of lists
 ## (4) isotropic notation 
 model <- list(list(model="bessel", var=2, kappa=0.5, scale=5),
               "+",
               list(model="nugget", var=1))
 cv - CovarianceFct(x, model=model)
              
 ## (5) anisotropic notation
 model <- list(list(model="bessel", var=2, kappa=0.5, aniso=0.2), 
               "+",
               list(model="nugget", var=1, aniso=1))
 cv - CovarianceFct(as.matrix(x), model=model)

 # The model gneitingdiff was defined in RandomFields v1.0.
 # This isotropic covariance function is valid for dimensions less
 # than or equal to 3 and has two positive parameters.
 # It is a class of models with compact support that allows for
 # smooth parametrisation of the differentiability up to order 6.     
 # The former model `gneitingdiff' must now be coded as
 gneitingdiff <-  function(p){
   list(list(m="gneiting", v=p[2], s=p[6]*p[4]),
       "*",
       list(m="whittle", k=p[5], v=1.0, s=p[4]),
       "+",
       list(m="nugget", v=p[3]))
 }
 # and then 
 param <- c(NA, runif(5,max=10)) 
 CovarianceFct(0:100,gneitingdiff(param))
 ## instead of formerly CovarianceFct(x,"gneitingdiff",param)

 # definition of a hypermodel is more complex
 model <- list(list(model="mastein", var=1,
                    aniso=c(1, -0.5, 0, 0), kappa=c(1, 0.5, 1.5)),
              "(",
              list(model="exp", var=1, aniso=c(0, 0, ,0, 1)))
 CovarianceFct(cbind(0:10, 0:10), model=model)

}
\keyword{spatial}






