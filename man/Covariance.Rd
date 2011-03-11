\name{Sophisticated Models}
\alias{sophisticated}
\alias{Sophisticated}
\title{Sophicated Covariance And Variogram Models}
\description{
  \code{Covariance} returns the values of complex stationary and
  nonstationary covariance functions;
   see \command{\link{CovarianceFct}} for basic isotropic models
 }
\details{
  \bold{Here only the non-isotropic and hyper models are listed;
    see \command{\link{CovarianceFct}} for basic isotropic models.}
  
  The implemented models are in standard notation for a
  covariance function (variance 1, nugget 0, scale 1) and for positive
  real arguments \eqn{h} (and \eqn{t}) for the stationary models or parts: 
  
  \itemize{
    \item \code{+}\cr
    Operator that adds up at most 10 submodels
    
    \item \code{*}\cr
    Operator that multiplies at most 10 submodels

    \item \code{$}\cr
    \deqn{C(x, y) = v C(x / s, y / s)}
    \deqn{C(x, y) = v C(x a, y a)}
    \deqn{C(x, y) = v C(A x, A y)}
    \deqn{C(x, y) = v C(p x, p y)}
    Operator that modifies the the variance (\eqn{v=}\code{var}) and the
    coordinates or distances by 
    \itemize{
      \item the scale (\eqn{s=}\code{scale}) or
      \item  the anisotropy matrix \eqn{a=}\code{anisoT} multiplied from the
      right or
      \item the anisotropy matrix \eqn{A} multiplied from the left or
      \item \eqn{p=}\code{proj} on a lower dimensional space along
      the coordinate axis
    }
    The parameter \code{scale} is positive, \code{aniso} and \code{A}
    are matrices, and \code{proj} is a vector indices with between 1 and
    the dimension of \eqn{x}.
    Note, at most one of the parameters,
    \code{anisoT}, \code{A}, \code{proj} may be given at the same time.

    The operator \code{$} has 1 submodel.
    If the dimension of the field is 1 or \code{aniso} is not given, the
    operator allows for derivatives.

    \item \code{ave1} % bernoulli 2010 paper example 13
    \deqn{C(h, u) = |E + 2Ahh^tA|^{-1/2}
      \phi(\sqrt(\|h\|^2/ 2 + (z^th + u)^2 (1 - 2h^tA (E+2Ahh^tA)^{-1} Ah)))}{
      C(h, u) = |E + 2 A h h^t A|^{-1/2}
      phi(sqrt[|h|^2/2 + (z^t h + u)^2 (1 - 2h^t A (E + 2 A h h^t A)^{-1} A h)])}
    where \eqn{E} is the identity matrix.
    
    \eqn{A} is a symmetric positive definite \eqn{(d-1) \times (d-1)}{(d-1)
      x (d-1)} and \eqn{z} is a \eqn{d-1} dimensional vector.
    The function \eqn{\phi} is normal mixture model, e.g.
    \code{whittle} model, see \command{\link{CovarianceFct}}
    and \command{\link{PrintModelList}}().

    
    \item \code{ave2} (nonstationary)\cr
    Here \eqn{C(h) = C_0(h, 0)} where \eqn{C_0} is the
    \code{ave1} model.
    
    
    \item \code{biWM} (bivariate model)\cr
   
    \deqn{   C_{ij}(h) = c_{ij} W_{\nu_{ij}} ( h / s_{ij}) }

    where \eqn{W_nu} is the \code{whittle} model and \eqn{i,j=1,2}.
    For (i=j) the  constants  \eqn{\nu_{ii}, s_{ii}, c_{ii} > 0}.
    For the offdiagonal elements with have \eqn{C_{12} = C_{21}},
    \eqn{s_{12}=s_{21}  > 0},
    \eqn{\nu_{12} =\nu_{21} = 0.5 (\nu_{11} + \nu{22}) /
      \nu_{red}}
    for some constant \eqn{\nu_{red} \in (0,1]}.
    The scalar \eqn{c_{12} =c_{21} = \rho_{red} \sqrt{f m c_{11} c_{22}}}
    where
    \deqn{f = \Gamma(\nu_{11} + d/2) * \Gamma(\nu_{22} + d/2) /
      \Gamma(\nu_{11}) / \Gamma(\nu_{22}) *
      (\Gamma(\nu_{12}) / \Gamma(\nu_{12}+d/2))^2 *
      ( s_{12}^{2*\nu_{12}} / s_{11}^{\nu_{11}} / s_{22}^{\nu_{22}} / )^2}
    and \eqn{\Gamma} is the Gamma function and \eqn{d} is the dimension
    of the space. The constant \eqn{m} is 
    the infimum of the function \eqn{g} on \eqn{[0,\infty)},
    \deqn{
      g(t) = (1/s_{12}^2 +t^2)^{2\nu_{12} + d}
      (1/s_{11}^2 + t^2)^{-\nu_{11}-d/2}
      (1/s_{22}^2 + t^2)^{-\nu_{22}-d/2}     
    }
    see the reference below for details on the infimum.

    The model now has the parameters

    
    \code{nu} \eqn{= (nu_{11}, nu_{22})}\cr
    \code{nured12} \eqn{=\nu_{red}}\cr
    \code{s} \eqn{= (s_{11}, s_{22})}\cr
    \code{s12} \eqn{= s_{12} = s_{21}}\
    \code{c} \eqn{= (c_{11}, c_{22})}\cr
    \code{rhored} \eqn{=\rho_{red}}
   See also \code{parsbiWM}.
  


   \item \code{constant}\cr
    This model is designes for the use in \command{\link{fitvario}}
    as a part of a linear model definition.
    Its only parameter is a covariance matrix of appropriate size
    to match the number of (non-repeated) observations or
    the number of columns of parameters \code{X} in model
    \command{mixed}, see \link{sophisticated}.
    
	
 
%    \item \code{coxisham} (non-separabel space time model)
%    \deqn{C(h,t) = (\det M_t)^{-1/2} \ehp(- (h-z t)^\top M_t (h- z t) /
%      det M_t)
%    }{C(h,t) = (det M)^{-1/2} exp(- (h-z t)^T M (h- z t) / det M)
%    }
%    where
%    \deqn{
%      M_t = \left( \begin{array}{cc} 1+t^2 & -\kappa_3 t^2 \\ -\kappa_3 t^2 &
%      1+t^2 \end{array}\right)
%    }{
%      M = rbind(c(1+t^2, -c t^2), c(-c t^2, 1+t^2))
%    }
%    The model is valid in 2 spatial dimensions.
%    The parameter \eqn{z} is passed as
%    \eqn{z=(\kappa_1, \kappa_2)}{z=(a, b)}.
%    The parameter \eqn{\kappa_3}{c} is in \eqn{[-1,1]}.


    \item \code{coxisham}
    \deqn{
      C(h, u) = |E + u^\beta D|^{-1/2} \phi([
      (h - u \mu)^t (E + u^\beta D)^{-1} (h - u\mu)]^{1/2})
    }
    Here $\eqn{mu} is vector;
    \eqn{E} is the identity matrix and \eqn{D} is a correlation 
    matrix with \eqn{|D| >0}. Currently implementation is done only for
    \eqn{d=2}.  The parameter \eqn{\beta}
    is in \eqn{(0,2]} and equals 2 by default.


    \item \code{curlfree} (multivariate)
    \deqn{
      ( - \nabla_x \nabla_x^T ) C_0(x, t)
    }
    \eqn{C_0} is a univariate
    covariance model that is motion invariant and at least twice differentiable
    in the first component. The operator is applied to the first
    component only. The model returns the potential field in the first
    component, the corresponding
    curlfree field and field of sources and sinks in the last component.
    The above formula for the covariance function only gives the part
    for the curlfree field.  The complete matrix-valued
    correlation function, including all
    components, is more complicated.

    \eqn{C_0} is either a spatiotemporal model (then \eqn{t} is the
    time component) or it is an isotropic model. Then, the first \eqn{Dspace}
    coordinates are considered as \eqn{x} coordinates and the remaining
    ones as \eqn{t} coordinates. By default, \eqn{Dspace} equals the
    dimension of the field (and \eqn{t} is identically \eqn{0}).

    See also the models \code{divfree} and \code{vector}.
   

    
    \item \code{cutoff} 
    \deqn{C(h)=\phi(h),                      0\le h \le d}
    \deqn{C(h) = b_0 ((dr)^a - h^a)^{2 a},   d \le h \le dr}
    \deqn{C(h) = 0,                           dr\le h}
   
    The cutoff model is a functional of the covariance function
    \eqn{\phi}{phi}.
    
    Here, \eqn{d>0} should be the diameter of the domain on which
    simulation is done . The parameter \eqn{a>0} has been shown to be
    optimal for  \eqn{a = 1/2} or \eqn{a =1}.
    
    The parameters \eqn{r} and \eqn{b_0}
    are chosen internally such that \eqn{C} is a smooth function.

    NOTE: The algorithm that checks the given parameters knows
    only about some few necessary conditions.
    Hence it is not ensured that
    the cutoff-model is a valid covariance function for any
    choice of phi and the parameters. 
    
    For certain models \eqn{\phi}{phi}, i.e. \code{stable},
    \code{whittle} and \code{gencauchy}, some sufficient conditions
    are known.

    

    \item \code{delayeffect} (bivariate)
    \deqn{ C_{11}(h) = C_{22}(h) = C_0(h) \qquad  C_{12}(h) = C_0(h +
    r), C_{21}(h) = C_0(-h + r)}{
      C_{11}(h) = C_{22}(h) = C_0(h),  C_{12}(h) = C_0(h + r), C_{21}(h) = C_0(-h + r)}
    Here \eqn{r} is a vector of the dimension of the random field,
    and \eqn{C_0} is a translation invariant, univariate covariance model.


     \item \code{divfree} (multivariate)
    \deqn{
     ( - \Delta E + \nabla \nabla^T ) C_0(x, t)
    }
    \eqn{C_0} is a univariate
    covariance model that is motion invariant and at least twice differentiable
    in the first component. The operator is applied to the first
    component only. The model returns the potential field in the first
    component, the corresponding
    divfree field and the field of curl strength in the last component.
    The above formula for the covariance function only gives the part
    for the divfree field.  The complete matrix-valued
    correlation function, including all
    components, is more  complicated.

    \eqn{C_0} is either a spatiotemporal model (then \eqn{t} is the
    time component) or it is an isotropic model. Then, the first \eqn{Dspace}
    coordinates are considered as \eqn{x} coordinates and the remaining
    ones as \eqn{t} coordinates. By default, \eqn{Dspace} equals the
    dimension of the field (and \eqn{t} is identically \eqn{0}).

    See also the models \code{curlfree} and \code{vector}.


    \item \code{EtAxxA} (auxiliary function)
    \deqn{
      S(x) = E + R^t A^t x x^t A R, \qquad x \in R^3
    }{
      S(x) = E + R^t A^t x x^t A R,  x \in R^3
    }
    where \eqn{E} and \eqn{A} are arbitrary \eqn{3 \times 3}{3 x 3} matrices
    and \eqn{R} is a rotation matrix,
    \deqn{
      R =\left(
      \begin{array}{lll}
      \cos (\alpha x_3) & -\sin(\alpha x_3)& 0
      \\
       \sin(\alpha x_3)& \cos (\alpha x_3) &  0
      \\
      0 & 0 & 1
      \end{array}
      \right)     
    }{R = matrix(c(cos(alpha x_3), -sin(alpha x_3), 0
             sin(alpha x_3),  cos(alpha x_3), 0,
             0,                0,             1), nc=3, by=TRUE)
	   }
   
    This is not a covariance function, but can be used as a submodel for
    certain classes of non-stationary covariance functions.

    
    \item \code{Exp}
    \deqn{ C(h) = \exp(-\gamma(h)) }
    where \eqn{\gamma} is a valid variogram.
    If a stationary covariance model \eqn{C}is given in stead of \eqn{\gamma},
    this is automatically turned into a variogram model, i.e.
    \eqn{ C(h) = \exp(-C(0)+C(h)) }.

    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    \item \code{M}
    \deqn{
      C(h) = M^t \phi(h) M
    }
    Here \eqn{phi} is a \eqn{k}-variate variogram or covariance,
    and \eqn{M} is any \eqn{m \times k}{m x k} matrix.

    \item \code{ma1}
    \deqn{
      C(h) = (\theta / (1 - (1-\theta) * C_0(h)))^\alpha
    }
    Here, \eqn{C_0} is any correlation function,
    \eqn{\alpha \in (0,\infty)} and \eqn{\theta \in (0,1)}.
    

    \item \code{ma2}
    \deqn{
      C(h) = (1 - exp(-\gamma(h)))/\gamma(h)
    }
    Here \eqn{\gamma} is a variogram model.

    \item \code{mastein}
    \deqn{C(h, t)=\frac{\Gamma(\nu + \gamma(t))\Gamma(\nu + \delta)}{
	\Gamma(\nu + 
	\gamma(t) + \delta) \Gamma(\nu)}
      W_{\nu + \gamma(t)}(\|h - Vt\|)
    }{C(h, t)= [ \Gamma(\nu + \gamma(t))\Gamma(\nu + \delta) ] / [
	\Gamma(\nu + 
	\gamma(t) + \delta) \Gamma(\nu) ]
      W_{\nu + \gamma(t)}(\|h - Vt\|)}
    
    \eqn{\Gamma} is the Gamma function; \eqn{\gamma(t)}
    is a variogram on the real axis; \eqn{W} is the Whittle-Matern model.
    Here, the
    names of covariance models can also be used; the algorithm
    chooses the corresponding variograms then.
    The parameter \eqn{\nu} is the smoothness parameter
    of the Whittle-Matern model (for \eqn{t=0}) and must be positive.
    Finally, \eqn{\delta} must be greater than or equal to half the
    dimension of \eqn{h}.
    Instead of the velocity parameter \eqn{V} in original model
    description, a preceeding anisotropy matrix is chosen appropriately:
    
    \deqn{\left( \begin{array}{cc} A & -V \\ 0 &
      1\end{array}\right)}{
      matrix(c(A,  -V,
               0,   1), nr=2, by=TRUE)
    }
    
    A is a spatial transformation matrix.
    (I.e. (x,t) is multiplied from left on the above matrix and
    the first elements of the obtained vector are intepreted as
    new spatial components and only these components are used to form
    the argument in the Whittle-Matern function.)
    The last component in the new coordinates is the time which is
    passed to \eqn{\gamma}{g}. (Velocity is assumed to be zero in
    the new coordinates.)
    
    %Note that
    %the anisotropy matrix must be such that \eqn{(x,t)} is transformed
    %into a purely spatial vector, i.e. the entries in
    %last column of the matrix are all naught.
    %On the other hand, all entries of the
    %anisotropy matrices in the submodels that build \eqn{\gamma}{gamma}
    %is naught except the very last, purely temporal one.

    Note, that for numerical reasons, \eqn{\nu+\gamma+d} may not exceed
    the value 80.0. If exceeded the algorithm fails.
    

    \item \code{mixed}
    This model is designed for the use in \command{\link{fitvario}}
    to build up linear regression models with fixed effects, mixed
    effects, including geoadditive parts.
    
    The model has two parameters. The first, \code{X}
    is a  matrix of independent variables.
    The second, \code{b}, is a vector of regression coefficients.
    Furthermore a submodel, \code{covb}, may give the covariance structure
    for \code{b}.

    Let \code{n} the number of (non-repeated) observations.
    The following combinations are allowed:
    \itemize{
      \item only \code{X} is given. Then \code{X} is a scalar
      or a vector of length \code{n}, and \code{X} defines a known mean.
      \item \code{X} and \code{b} are given. Then \code{X} is a
      \eqn{(n \times m)}{n x m} matrix
      where \eqn{m} is the length of the vector \eqn{b}. Then a fixed
      effect is defined.
      \item \code{X} and \code{covb} are given.
      \itemize{
	\item if \code{covb} is the model \cite{constant}, 
	then we have a random model (maybe with preceeding model
	\code{$}).
	\item if \code{covb} is any other model then we have
	a geoadditive part
      }
    }

    The data in the fitvario may contain \code{NA}s, but not \code{X}.
    
    \item \code{mqam} (multivariate quasi-arithmetic mean)
    \deqn{C_{ij}(h) = \rho_{ij} \phi(\theta \phi^{-1}(C_i(h)) + (1 - \theta) \phi^{-1}(C_j(h)))}
    where \eqn{\phi} is a completely monotone function
    and \eqn{C_i} are suitable covariance functions.
    
    The submodel \eqn{\phi} is given (by name) as first submodel.
    Since \eqn{\phi} is completely monotone if and only if
    \eqn{\phi(\| . \|^2)}{\phi( \| \cdot \|^2)} is a valid
    covariance function for all dimensions, e.g. \code{stable},
    \code{gauss}, \code{exponential},  \eqn{\phi} is given
    by the name of the corresponding covariance function \eqn{C},
    i.e. \eqn{phi( . ) = C(sqrt( . ))}{\phi(\cdot) = C(\sqrt(\cdot))}.

    Warning: \code{RandomFields} cannot check whether the combination
    of \eqn{\phi} and \eqn{C_i} is valid.
  

    \item \code{natsc}
    \deqn{    C(h) = C_0(h / s)  }
    Where \eqn{C_0} is any stationary and isotropic model. The parameter
    \code{s} is chosen by \code{natsc} such that the practical range
    (or the mathematical range, if finite) is 1.

  
   \item \code{nonstWM}
    \deqn{C(x, y)=\Gamma(\mu) \Gamma(\nu(x))^{-1/2} \Gamma(\nu(y))^{-1/2}
      W_{\mu} (\|x-y\|)}
    \deqn{= 2^{1-\mu}  \Gamma(\nu(x))^{-1/2} \Gamma(\nu(y))^{-1/2}
      \|x-y\|^\mu
      K_\nu(\|x-y\|)}
    where  \eqn{\mu = [\nu(x) + \nu(y)]/2} and
    \eqn{\nu} is a positive function.
    If \eqn{\nu} is a scalar use the variable \code{nu}.
    If \eqn{\nu} is a function, use the submodel \code{Nu}.
    Note that for \code{Nu} the usual list structure applies
    and only the defined covriance models can be used.

    
    \item \code{nsst} (Non-Separable Space-Time model)
    \deqn{C(h,u)= (\psi(u)+1)^{-\delta/2} \phi(h /\sqrt(\psi(u) +1))}
    The parameter \eqn{\delta} must be greater than or equal to
    the spatial dimension of the field. 
    \eqn{\phi} is normal mixture model and
    \eqn{\psi} is a variogram.
    
    This model is used for space-time modelling where the spatial
    component is isotropic.\cr
    

    \item \code{nugget} (multivariat model)
    \deqn{C(h)= {\rm diag}(1,\ldots,1) 1_{\{0\}}(h)}{diag(1,\ldots,1) 1(h==0)}
    The components of the multivariate vector are always independent.
    The models adapts the multivariate dimension to the calling model.
   

    \item \code{parsbiWM} (bivariate model)\cr
     
    \deqn{   C_{ij}(h) = c_{ij} W_{\nu_{ij}} (h / s) }

    where \eqn{W_nu} is the \code{whittle} model and \eqn{i,j=1,2}.
    For (i=j) the  constants  \eqn{\nu_{ii}, c_{ii} \ge 0} and \eqn{s>0}.
    For the offdiagonal elements with have \eqn{C_{12} = C_{21}}.
    Furthermore,  \eqn{\nu_{12} =\nu_{21} = 0.5 (\nu_{11} + \nu{22})}
    and the scalar \eqn{c_{12} =c_{21} = \rho_{red} \sqrt{f m c_{11} c_{22}}}
    where
    \deqn{f = \Gamma(\nu_{11} + d/2) * \Gamma(\nu_{22} + d/2) /
      \Gamma(\nu_{11}) / \Gamma(\nu_{22}) *
      (\Gamma(\nu_{12}) / \Gamma(\nu_{12}+d/2))^2}
    and \eqn{\Gamma} is the Gamma function and \eqn{d} is the dimension
    of the space. The constant \eqn{m} is 
    the infimum of the function \eqn{g} on \eqn{[0,\infty)}, 
    \deqn{
      g(t) = (1/s_{12}^2 +t^2)^{2\nu_{12} + d}
      (1/s_{11}^2 + t^2)^{-\nu_{11}-d/2}
      (1/s_{22}^2 + t^2)^{-\nu_{22}-d/2}
    }
    see the reference below for details on the infimum.

    The model now has the parameters

    \code{nu} \eqn{= (\nu_{11}, \nu_{22})}\cr
    \code{s} \eqn{= (s_{11}, s_{22})}\cr
    \code{s12} \eqn{= s_{12} = s_{21}}\cr
     \code{c} \eqn{= (c_{11}, c_{22})}\cr
    \code{rhored} \eqn{=\rho_{red}}
  
    See also \code{biWM}.

    
    \item \code{Pow}
    \deqn{ \gamma(h) = (\gamma_0(h))^\alpha}
    or 
    \deqn{ C(h) = C_0(0) - [C_0(0) - C_0(h)]^\alpha}
    where \eqn{\gamma_0} is a valid variogram or
    \eqn{C_0} is a valid covariance function,
    and \eqn{\alpha \in [0, 1]}.
    
    
    \item \code{qam} (Quasi-arithmetic mean)
    \deqn{C(h) = \phi(\sum_i \theta_i \phi^{-1}(C_i(h)))}
    where \eqn{\phi} is a completely monotone function
    and \eqn{C_i} are suitable covariance functions.
    
    The submodel \eqn{\phi} is given (by name) as first submodel.
    Since \eqn{\phi} is completely monotone if and only if
    \eqn{\phi(\| . \|^2)}{\phi( \| \cdot \|^2)} is a valid
    covariance function for all dimensions, e.g. \code{stable},
    \code{gauss}, \code{exponential}, \eqn{\phi} is given
    by the name of the corresponding covariance function \eqn{C},
    i.e. \eqn{phi( . ) = C(sqrt( . ))}{\phi(\cdot) = C(\sqrt(\cdot))}.

    Warning: \code{RandomFields} cannot check whether the combination
    of \eqn{\phi} and \eqn{C_i} is valid.
 
    
   \item \code{rational} (auxiliary)
    \deqn{
      S(x) = (a_0 + a_1 * x^t A A^t x)/ (1 + x^t A A^t x)
    }
    where is some \eqn{d\times d}{d x d} matrix
    and \eqn{a = (a_0, a_1)} is a 2-dimensional vector.

 
    \item \code{Rotat}(auxiliary function)
    \deqn{
      S^t(x) = x^t R, \qquad x \in R^3
    }{
      S^t(x) = x^t R, x \in R^3
    }
    where 
    and \eqn{R} is a rotation matrix,
    \deqn{
      R =\left(
      \begin{array}{lll}
      \cos (\alpha x_3) & -\sin(\alpha x_3)& 0
      \\
       \sin(\alpha x_3)& \cos \alpha x_3 &  0
      \\
      0 & 0 & 1
      \end{array}
      \right)     
    }{R = matrix(c(cos(alpha x_3), -sin(alpha x_3), 0
             sin(alpha x_3),  cos(alpha x_3), 0,
             0,                0,             1), nc=3, by=TRUE)
	   }
	  
    This is not a covariance function, but can be used a submodel for
    certain classes of non-stationary covariance functions.
    
    
    \item \code{Stein} 
    \deqn{C(h)=a_0 + a_2 (h)^2 + \phi(h), 0\le h \le D}
    
    \deqn{C(h)=b_0 (rD - h)^3/(h),   r \le h \le rD}
    
    \deqn{C(h) = 0,   rD \le h}
   
    The Stein model is a functional of the covariance function
    \eqn{\phi}{phi}. 
    
     Here, \eqn{D>0} should be the diameter of the domain on which
    simulation is done,\eqn{r\ge1}.
    The parameters \eqn{a_0}, \eqn{a_2} and \eqn{b_0}
    are chosen internally such that \eqn{C} becomes a smooth function.
    
   % Note that the scale parameter for the submodel must be \code{double(0)}.

    NOTE: The algorithm that checks the given parameters knows
    only about some few necessary conditions.
    Hence it is not ensured that
    the Stein-model is a valid covariance function for any
    choice of phi and the parameters. 
    
    For certain models \eqn{\phi}{phi}, i.e. \code{stable},
    \code{whittle}, \code{gencauchy}, and the variogram
    model \code{fractalB}
    some sufficient conditions are known.

 
    \item \code{steinst1} (non-separabel space time model)
    \deqn{C(h, t) = W_{\nu}(y) -
      \frac{\langle h, z \rangle t}{(\nu - 1)(2\nu
	+ d)} W_{\nu -1}(y)
    }{
      C(h, t) = W_\nu(y) - <h, z> t W_{\nu-1}(y) / [ (\nu-1) (2\nu + d
    +1) ]}

    Here, \eqn{W_{\nu}} is the Whittle-Matern model with
    smoothness parameter \eqn{\nu};
    \eqn{y=\|(h,t)\|}{y = ||(h,t)||}.
    \eqn{z} is a vector whose norm 
    must less than or equal to 1.

    \item \code{stp}
    \deqn{
      C(x,y) = |S_x|^{1/4} |S_y|^{1/4} |A|^{-1/2}
      \phi(Q(x,y)^{1/2})
    }
    where
    \deqn{
      Q(x,y) = c^2 - m^2 + h^t (S_x + 2(m + c)M) A^{-1} (A_y + 2
      (m-c)M)h,
    }
    \deqn{
      c = -z^t h + \xi_2(x) - \xi_2(y),
    }
    \deqn{
      A = S_x + S_y + 4 M h h^t M      
    }
    \deqn{
      m = h^t M h
      .}
    \deqn{
      h = H(x) - H(y)
    }
    The parameters are
    \itemize{
      \item \eqn{S_x} (strictly) positive definite matrices for \eqn{x
	\in R^d}
      \item \eqn{M} an arbitrary \eqn{d \times d}{d x d} matrix
      \item \eqn{z \in R^d}  arbitrary
      \item \eqn{H} arbitrary d-variate function on \eqn{R^d}
      \item \eqn{\xi} arbitrary univariate function on \eqn{R^d}
      \item \eqn{\phi} a normal mixture model
    }
    The model allows for mimicking cyclonic behaviour.
  
   

    \item \code{tbm2}
    \deqn{
      C(h) = \frac{d}{dh}\int_0^h \frac{u\phi(u)}{\sqrt{h^2 - u^2}} d u
    }{
      C(h) = d/dh int_0^h [ u phi(u) ] / [ sqrt{h^2 - u^2} ] d u
    }
    for some stationary and isotropic covariance \eqn{\phi}
    that is valid in at least 2 dimensions.

    This operator is currently only designed for internal use!

    \item \code{tbm3}
    \deqn{
      C(h) = \phi(h) + h \phi'(h) / n
    }{
      C(h) = phi(h) + h phi'(h) / n
    }
 
    which, for \code{n=1} reduced to the standard TBM operator
    \deqn{
      C(h) =\frac {d}{d h} h \phi(h)
    }{
      C(h) = d/dh [ h phi(h) ]
    }
    for some stationary and isotropic covariance \eqn{\phi}
    that is valid in at least \eqn{n+2} dimensions.
    \code{n} should be an integer.
    
    This operator is currently only designed for internal use!
    
    \item \code{vector} (multivariate)
    \deqn{
      ( -0.5 * (a + 1) \Delta E + a \nabla \nabla^T ) C_0(x, t)
    }
    \eqn{C_0} is a univariate
    covariance model that is motion invariant and at least twice differentiable
    in the first component. The operator is applied to the first
    component only. The parameter \eqn{a} is in \eqn{[-1, 1]}. If \eqn{a=-1}
    then the field is curl free; if \eqn{a=1} then the field is
    divergence free.

    \eqn{C_0} is either a spatiotemporal model (then \eqn{t} is the
    time component) or it is an isotropic model. Then, the first \eqn{Dspace}
    coordinates are considered as \eqn{x} coordinates and the remaining
    ones as \eqn{t} coordinates. By default, \eqn{Dspace} equals the
    dimension of the field (and \eqn{t} is identically \eqn{0}).

    See also the models \code{divfree} and \code{curlfree}
  

     
  
  }
  
  See \command{\link{CovarianceFct}} for comments on the use of a
  covariance model.
  
  However, for the above sophicated models,
  the following differences should be considered:
  \itemize{
    \item
    \command{\link{RFparameters}}\code{()$PracticalRange}
    is usually not defined for the above models
    \item
    only the list notation 
    can be used, but not the simple model definitions
    with \code{model="name"} and
    \code{param=c(mean, variance, nugget, scale,...)}.
    \item
    the use of \command{Covariance} is obligatory if the model is
    non-stationary.
    \item
    the anisotropy matrix belonging to a hypermodel 
    is applied first to the coordinates before any
    call of the submodels.
  }
  
  To use the above models, a new, very flexible, straight forward
  list notation is needed. Background of this notation is that
  we have \sQuote{primitives}, i.e. functions that are positive
  definite. And we have \sQuote{operators}, i.e. functionals that
  make out of given variograms, covariance functions etc. new
  models. Examples are \code{"+"}, \code{"*"}, or Gneiting's
  \code{"nsst"}. Consequently, we need also an operator, called
  \code{"$"}, that changes the variance and the scale.

  E.g. a standard exponential model (variance=1, scale=1, nugget=0)
  is now simply written as
 
  \deqn{\hbox{list("exponential")}}{list("exponential")}

  (And no \code{param} must be given!)
  
  Further, a standard exponential model with a nugget effect,
  nugget variance 3, is now written as
  \cr
 
  list("+",\cr
         list("exponential"),\cr
         list("$", var=3, list("nugget"))\cr
  )\cr

  Here, only the relevant parameters need to be given; the missing
  parameters get standard values whenever standard values exist,
  e.g. variance equals 1 if not given.
  Further, the parameters can (and must) be called by names, which makes
  complex models much more readable.
  Submodels, as \code{list("exponential")} in the second example above,
  can (but need not) be called by name.
}


\value{
  \code{CovarianceFct} and \code{Covariance}
  return a vector of values of the covariance function.
}

\references{
  Overviews:
  \itemize{
    \item see reference list in \code{\link{CovarianceFct}}
  }

  ave1, ave2
  \itemize{
    \item Schlather, M. (2010)
    On some covariance models based on normal scale mixtures.
    \emph{Bernoulli}, \bold{16}, 780-797. (Example 13)
  }

  biWM, parsbiWM
   \itemize{
    \item Gneiting, T., Kleiber, W., Schlather, M. (2011)
    Matern covariance functions for multivariate random fields
    \emph{JASA}
  }

 coxisham
 \itemize{
   \item Cox, D.R., Isham, V.S. (1988)
   A simple spatial-temporal model of rainfall.
   \emph{Proc. R. Soc. Lond. A}, \bold{415}, 317-328.

    \item Schlather, M. (2010)
    On some covariance models based on normal scale mixtures.
    \emph{Bernoulli}, \bold{16}, 780-797.
  }

  curlfree
  \itemize{
    \item see vector
  }
  
  cutoff
  \itemize{
    \item Gneiting, T., Sevecikova, H, Percival, D.B., Schlather M.,
    Jiang Y. (2006) Fast and Exact Simulation of Large {G}aussian
    Lattice Systems in {$R^2$}: Exploring the Limits.
    \emph{J. Comput. Graph. Stat.} \bold{15}, 483-501.
    \item Stein, M.
  }


  delayeffect
  \itemize{
    \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
     Springer, 3nd edition.
   }
   

   divfree
   \itemize{
   \item see vector
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


  vector
  \itemize{
    \item Fuselier, E.J. (2006) \emph{Refined Error Estimates for
      Matrix-Valued Radial Basis Functions}
    PhD thesis. Texas A&M University
    \item Scheuerer, M. and Schlather, M. (2011)
    Covariance Models for Random Vector Fields
    \emph{Submitted}
  }


  Ma-Stein model
  \itemize{
    \item Ma, C. (2003)
    Spatio-temporal covariance functions generated by mixtures.
    \emph{Math. Geol.}, \bold{34}, 965-975.
    \item Stein, M.L. (2005) Space-time covariance functions. \emph{JASA},
    \bold{100}, 310-321.
  }
  
  
  ma1/ma2
  \itemize{
    \item
  }

  mixed
  \itemize{
    \item Ober, U., Erbe, M., Porcu, E., Schlather, M. and Simianer, H.
    (2011) Kernel-Based Best Linear Unbiased Prediction with Genomic
    Data.
    \emph{Submitted.}
  }
  
  nonstWM/hyperbolic/cauchy
  \itemize{
    \item Stein, M. (2005)
    Nonstationary Spatial Covariance Functions.
    Tech. Rep., 2005
  }

  nsst
  \itemize{
    \item Gneiting, T. (1997) Normal scale mixtures and dual probability
    densitites, \emph{J. Stat. Comput. Simul.} \bold{59}, 375-384.
    
    \item Gneiting, T. (2002) Nonseparable, stationary covariance
    functions for space-time data, \emph{JASA} \bold{97}, 590-600.
    
    \item Gneiting, T. and  Schlather, M. (2001)
    Space-time covariance models.
    In El-Shaarawi, A.H. and Piegorsch, W.W.:
    \emph{The Encyclopedia of Environmetrics.} Chichester: Wiley.

    \item Zastavnyi, V. and Porcu, E. (2011)
    Caracterization theorems for the Gneiting class space-time
    covariances.
     \emph{Bernoulli}, \bold{??}.

     \item Schlather, M. (2010)
     On some covariance models based on normal scale mixtures.
     \emph{Bernoulli}, \bold{16}, 780-797.
  }

  Quasi-arithmetic means (qam, mqam)
  \itemize{
    \item Porcu, E., Mateu, J. & Cchristakos, G. (2007) Quasi-arithmetic
    means of covariance functions with potential applications to
    space-time data.  Submitted to Journal of Multivariate Analysis.
    \item
  }
  
  Paciorek-Stein (steinst1)
  \itemize{
    \item Stein, M. (2005)
    Nonstationary Spatial Covariance Functions.
    Tech. Rep., 2005
    \item Paciorek, C. (2003)
    \emph{Nonstationary Gaussian Processes for Regression and Spatial
      Modelling},
    Carnegie Mellon University, Department of Statistics, PhD thesis.
  }


  Stein
  \itemize{
    \item Stein, M.
  }
  
  stp
  \itemize{
    \item Schlather, M. (2008)
    On some covariance models based on normal scale mixtures.
    \emph{Submitted}
  }
  

  tbm
  \itemize{
    \item Gneiting, T. (1999)
    On the derivatives of radial positive definite function.
    \emph{J. Math. Anal. Appl}, \bold{236}, 86-99

    \item
    Matheron, G. (1973).
    The intrinsic random functions and their applications.
    \emph{Adv . Appl. Probab.}, \bold{5}, 439-468.
  }
    

}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}
}
\seealso{
  \command{\link{CovarianceFct}},
  \command{\link{EmpiricalVariogram}},
  \command{\link{GetPracticalRange}},
  \command{\link{parameter.range}},
  \code{\link{RandomFields}},
  \command{\link{RFparameters}},
  \command{\link{ShowModels}}.
}

\examples{
%  library(RandomFields,lib="~/TMP");  source("/home/schlather/R/RF/RandomFields/tests/source.R")
% RFparameters(Print=9)
PrintModelList(op=TRUE)

## the subsequent model can be used to model rainfall...
y <- x <- seq(0, 10, len=25) # better 256 -- but will take a while 
T <- c(0, 10, 1) # better 0.1
col <- c(topo.colors(300)[1:100], cm.colors(300)[c((1:50) * 2, 101:150)])

model <- list("coxisham", mu=c(1, 1), D=matrix(nr=2, c(1, 0.5, 0.5, 1)),
              list("whittle", nu=1)
              )

system.time(z <- GaussRF(x, y, T=T, grid =TRUE, spectral.lines=1500,
                       model = model))

zlim <- range(z)
time <- dim(z)[3]
for (i in 1:time) {
  Print(i)
  sleep.milli(100)
  image(x, y, z[, , i], add=i>1, col=col, zlim=zlim)
}


####################################################
####################################################

 # the following five model definitions are the same!
 ## (1) very traditional form
 (cv <- CovarianceFct(x, model="bessel", param=c(NA, 2 , 1, 5, 0.5)))

 ## (2) traditional form in list notation
 model <- list(model="bessel", param=c(NA, 2, 1, 5, 0.5))
 cv - CovarianceFct(x, model=model)

 ## (3) nested model definition
 cv - CovarianceFct(x, model="bessel",
                    param=rbind(c(2, 5, 0.5), c(1, 0, 0)))

 #### most general notation in form of lists
 ## (4) isotropic notation 
 model <- list("+",
               list("$", var=2, scale=5, list("bessel", 0.5)),
               list("nugget"))
 cv - CovarianceFct(x, model=model)
              
 ## (5) anisotropic notation
 model <- list("+",
               list("$", var=2, aniso=0.2, list("bessel", 0.5)),
               list("nugget"))
 cv - CovarianceFct(as.matrix(x), model=model)




####################################################
####################################################

 # The model gneitingdiff was defined in RandomFields v1.0.
 # This isotropic covariance function is valid for dimensions less
 # than or equal to 3 and has two positive parameters.
 # It is a class of models with compact support that allows for
 # smooth parametrisation of the differentiability up to order 6.     
 # The former model `gneitingdiff' should now be coded as

 gneitingdiff <- function(p){
    list("+",
         list("$", var=p[3], list("nugget")),
         list("$", scale=p[4],
              list("*", 
                   list("$", var=p[2], scale=p[6], list("gneiting")),
                   list("whittle", nu=p[5])
                  )
              )
         )
 }

 # and then 
 param <- c(NA, runif(5, max=10)) 
 CovarianceFct(0:100, model=gneitingdiff(param))
 ## instead of formerly CovarianceFct(x,"gneitingdiff",param)
}


\keyword{spatial}






