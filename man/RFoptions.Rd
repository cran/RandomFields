\name{RFoptions}
\alias{RFoptions}
\title{Setting control arguments}
\description{
 \command{\link{RFoptions}} sets and returns control arguments for the analysis
 and the simulation of random fields
}
\usage{
 RFoptions(..., no.readonly = TRUE)
}
\arguments{
 \item{...}{arguments in \code{tag = value} form, or a list of tagged
 values.}
 \item{no.readonly}{If \command{\link{RFoptions}} is called without
   argument then all arguments are returned in a list. If
   \code{no.readonly=TRUE} then only rewritable arguments are returned.
 }
}
\details{
  The subsections below comment on\cr
  \bold{1. \code{general}: General options}\cr
  \bold{2. \code{br}: Options for Brown-Resnick
    Fields}\cr
  \bold{3. \code{circulant}: Options for circulant embedding methods
    \command{\link{RPcirculant}}}\cr
  \bold{4. \code{coords}: Options for coordinates and units, see
    \link{coordinate systems}}\cr
  \bold{5. \code{direct}: Options for simulating by simple matrix decomposition}\cr
  \bold{6. \code{distr}: Options for distributions, in particular \command{\link{RRrectangular}}}\cr
  \bold{7. \code{empvario}: Options for calculating the empirical variogram}\cr
  \bold{8. \code{fit}: Options for \command{\link{RFfit}},
    \command{\link{RFratiotest}}, and \command{\link{RFcrossvalidate}}}\cr
  \bold{9. \code{gauss}: Options for simulating Gaussian random fields}\cr
  \bold{10. \code{graphics}: Options for graphical output}\cr
  \bold{11. \code{gui}: Options for \command{\link{RFgui}}}\cr
  \bold{12. \code{hyper}: Options for simulating hyperplane tessellations}\cr
  \bold{13. \code{krige}: Options for Kriging}\cr
  \bold{14. \code{maxstable}: Options for simulating max-stable random fields}\cr
  \bold{15. \code{mpp}: Options for the random coins (shot noise) methods}\cr
  \bold{16. \code{nugget}: Options for the nugget effect}\cr
  \bold{17. \code{registers}: Register numbers}\cr
  \bold{18. \code{sequ}: Options for the sequential method}\cr
  \bold{19. \code{solve}: Options for solving linear systems}\cr
  \bold{20. \code{special}: Options for some special methods}\cr
  \bold{21. \code{spectral}: Options for the spectral (turning bands) method}\cr
  \bold{22. \code{tbm}: Options for the turning bands method}\cr
  \bold{23. \code{internal}: Internal}\cr
  \bold{General comments}\cr\cr\cr
  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \bold{16. Options for RFloglikelihood}\cr
% 
% "auto", "full", "composite", "selection"




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{1. General options}
\describe{
  % only works for RFsimulate, not for RFsimulate
  \item{\code{allowdistanceZero}}{boolean. Only used in
    \command{\link{RFinterpolate}} and in \command{\link{RFfit}}.
    If true, then
    multiple observations or identical locations
    are allowed within a single data set.
    In this case, the coordinates are slightly scattered, so that
    the points have some tiny distances.
    
    Default: \code{FALSE}. 
  }
  
  
  \item{\code{cPrintlevel}}{
    \code{cPrintlevel} is automatically set to \code{printlevel}
    when \code{printlevel} is changed.
    Standard users will never use a value higher than 3.
    
    0 : no messages\cr
    1 : messages and warnings when the user's input looks odd\cr
    2 : messages (and internal errors) documenting the choice of the
    simulation method\cr
    3 : further user relevant informations\cr
    4 : information on recursive function calls\cr
    5 : function flow information of central functions \cr
    6 : errors that are internally treated\cr
    7 : details on building up the covariance structure\cr
    8 : details on taking the square root of the covariance matrix\cr
    9 : details on intermediate calculations\cr
    10 : further details on intermediate calculations\cr
    
    Note that \code{printlevel} works
    on the R level whereas \code{cPrintlevel} works on the C level.
    
    Default: 1 \cr % [also  do].\cr
  }
  
  \item{\code{detailed_output}}{logical.
    if \code{TRUE} some function, e.g. \code{\link{RFcrossvalidate}}
    will return additional information.
  }
  
  \item{\code{every}}{integer.
    if greater than zero, then every \code{every}th iteration is
    printed if simulated by TBM or random coin method. The value zero
    means that nothing is printed.
    
    Default: \code{0} % [do].
  }
  
  
  \item{\code{exactness}}{logical or NA. Currently only used when simulating
    Gaussian random fields.
    \itemize{
      \item \code{TRUE}: \command{\link{RPcoins}},
      \command{\link{RPhyperplane}}, \command{\link{RPsequential}},
      \command{\link{RPspectral}} and \command{\link{RPtbm}} and
      \emph{approximative} circulant embedding  are excluded.
      If the circulant embedding method is considered as badly
      behaved, then the matrix decomposition methods are preferred.
      
      \item \code{FALSE}: all the methods are allowed.
      If the circulant embedding method is
      considered as badly behaved or the number of points to be
      simulated is large, the turning bands methods are
      rather preferred.
      
      \item \code{NA}: Similar to \code{FALSE}, but
      some inexact algorithms get less preference.
    }
    
    Default: \code{NA} .
  }
  
  \item{\code{expected_number_simu}}{positive integer which is usally set
    internally as the value of the argument \code{n} in
    \command{\link{RFsimulate}}. The argument \code{expected_number_simu}
    should be set only by an advanced users and only if 
    \command{\link{RFsimulate}} will be called with argument \code{n} alone.
  }
  
  
  \item{\code{gridtolerance}}{
    used in \command{\link{RFsimulate}} to see if the coordinates build a
    grid for x, y, z, T-values. This argument is also used
    in case of conditional
    simulation where the data locations might ly on
    a grid.
    
    Default: \code{1e-6}
  }

  \item{asList}{logical. Lists of arguments are treated slightly
    different from non-lists. If \code{asList=FALSE} they are treated the
    same way as non-lists. This options being set to \code{FALSE} after
    calling \command{RFoptions} it should be set as first element of a list.
    
    Default: \code{TRUE}
  }
  
  \item{\code{modus_operandi}}{character. One of the values
    \code{"careless"}, \code{"sloppy"}, \code{"easygoing"},
    \code{"normal"}, \code{"precise"}, \code{"pedantic"},
    \code{"neurotic"} .
    \bold{This argument is in an experimental stage and its definition
      and effects will change very likely in near future.}%to do
    This argument sets a lot of argument at once related to estimation
    and simulation. \code{"careless"} prefers rather fast algorithms,
    but the results
    might be very rough approximations. By way of contrast,
    \code{"neurotic"} will try very
    hard to return exact result at the cost of hugh computing times.
    
    Default: \code{"normal"}
  }
  
  
  \item{\code{na_rm_lines}}{
    logical. If \code{TRUE} then a line of the data that contains a
    \code{NA} value is deleted. Otherwise it is tried to deal with the
    \code{NA} value at higher costs of computing time.
    
    Default: \code{FALSE}.
  }
  
  \item{\code{output}}{character.
    one of the values \code{"sp"} (if and only if
    \code{spConform=TRUE}), 
    \code{"RandomFields"} (if and only if \code{spConform=FALSE}),
    \code{"geoR"}.

    The output mode \code{geoR} currently adds some attributes such as
    the call of the function.

    NOTE: \code{output} is in an experimental stage, whose effects might
    change in future. Currently, \code{output} changes the values of
    \code{reportcoord}, \code{returncall} and \code{spConform}.
  }
  
  \item{\code{pch}}{character.
    \command{\link{RFfit}}: shown before evaluating any method;
    if \code{pch!=""} then one or two
    additional steps in the MLE methods are
    marked by \dQuote{+} and \dQuote{#}.
    
    Simulation: 
    
    The character is printed after each
    performed simulation if more than one simulation is performed at
    once. If \code{pch='!'} then an absolute
    counter is shown instead of the character.
    If \code{pch='\%'} then a
    counter of percentages is shown instead of the character.
    Note that also \sQuote{\eqn{\mbox{\textasciicircum}}{^}H}s are printed in
    the last two cases, 
    which may have undesirable interactions with some few other R
    functions, e.g. \command{\link[utils]{Sweave}}.
    
    Default: \code{'*'}. % [do]. 
  }
 
 
 

  \item{\code{practicalrange}}{logical or integer.
    If not \code{FALSE} the range of primitive
    covariance functions is
    adjusted so that cov(1) is zero for models with finite range.
    (Operators are too complex to be adjusted; for anisotropic
    covariance the practical range is not well defined.)
    
    The value of cov(1) is about 0.05 (for \code{scale=1})
    for models without range. See \command{\link{RMmodel}} or type
    \cr
    \code{\link{RFgetModelNames}(type="positive definite",
      domain="single variable", isotropy="isotropic", operator=FALSE, vdim=1)}
    \cr
    for the list of primitive models.
    \itemize{
      \item \code{FALSE} : the practical range ajustment is not used.
      \item \code{TRUE} : \code{practicalrange} is applicable only if
      the value is known exactly, or, at least, can be approximated by
      a closed formula.
      \item \code{2} : if the practical range is not known exactly it
      is approximated numerically.	
    }
    
    Default: \code{FALSE} .
  }
  
  \item{\code{printlevel}}{If \code{printlevel}\eqn{\le0}{<=0}
    there is not any output on the screen. The
    higher the number the more tracing information is given. 
    Standard users will never use a value higher than 3.

    0 : no messages\cr
    1 : important (error) messages and warnings\cr
    2 : less important messages\cr
    3 : details, but still for the user\cr
    4 : recursive call tracing (only used within \command{\link{RFfit}})\cr
    5 : function flow information of large functions\cr
    6 : errors that are internally treated\cr
    7 : details on intermediate calculations\cr
    8 : further details on intermediate calculations\cr
    
    Default: 1 %[also do].\cr
  }

  \item{reportcoord}{character.
    Current values are \code{"always"}, \code{"important"},
    \code{"warn"}, \code{"never"},

    Both \code{"warn"} and \code{"important"} have any effect only
    if the coordinate system is changed internally. In this case
    \code{"warn"} yields a displayed warning message whereas
    \code{"important"} adds an attribute to the result as in the
    case \code{"always"}.
    
    If \code{"always"} or \code{"important"}
    the reports are added as attribute to the results.
    Note that in this case the class of the result may change
    (e.g. from \code{"numeric"} to \code{"atomic"}).

    Default: \code{"warn"}
  }

  \item{returncall}{logical. If \code{TRUE}
    then the call is returned as an attribute

    Default: \code{TRUE}
  }
   
  \item{seed}{integer. If \code{NULL} or \code{NA}
    \command{\link{set.seed}} is \bold{not} called.
    Otherwise, \code{\link[base]{set.seed}(seed)} is set
    before simulations are performed, e.g. by
    \command{\link{RFsimulate}} or \command{\link{RFdistr}}.
    
    If the argument is set locally, i.e., within a function,
    it has the usual local effect. If it is set globally, i.e. by
    \command{RFoptions} the \code{seed} is fixed
    for \bold{all subsequent} calls.
    
    If the number of simulations \code{n} is greater than one
    and if \code{RFoptions(seed=seed)} is set, the \eqn{i}th
    simulation is started with the seed \sQuote{\code{seed}\eqn{+i-1}}.
    %   The function \code{set.seed} should not be used in case \code{n}
    %    is greater than 1.
    %
    %Vgle!   
    %set.seed(5)
    %RFsimulate(RPschlather(RMmatern(nu=2), xi=1, mu=1, s=1), x, grid=F, n=5)@data
    %set.seed(5)
    %RFsimulate(RPschlather(RMmatern(nu=2.01), xi=1, mu=1, s=1), x,grid=F,n=5)@data
    %RFoptions(cPr=3, seed=5)
    %RFsimulate(RPschlather(RMmatern(nu=2), xi=1, mu=1, s=1), x, grid=F, n=5)@data
    %RFsimulate(RPschlather(RMmatern(nu=2.01), xi=1, mu=1, s=1), x, grid=F,n=5)@data

    
    Note also that \command{\link{RFratiotest}} has its own argument
    \code{seed} with a slightly different meaning.
  }

  \item{\code{set}}{integer.
    Certain models (e.g. \command{\link{RMfixcov}} and
    \command{\link{RMcovariate}}) 
      allow for lists as arguments.
    \code{set} selects a certain list element.
    If necessary the list is recycled.
  }
  

  \item{\code{spConform}}{logical.
    \code{spConform=TRUE} might be used by
    a standard user as this allows the comfortable use of \command{plot},
    for instance, while \code{spConform=FALSE} is \bold{much} faster and
    and consumes \bold{much less memory}, hence might
    be used by programmers or advanced users.
    
    Details: if \code{spConform=TRUE} then \command{\link{RFsimulate}} and
    many other functions 
    return an \code{sp}-object (which is an S4 object). Otherwise, matrices
    or lists are 
    returned as defined in RandomFields 2.0, see the manuals for the
    specific functions. Frequently, the latter have now a class attribute
    to make the output nicer.
    
    Note: for large data sets (to be generated),
    \code{spConform=TRUE} should \bold{not} be used.

    See also \code{output}.
    
    Default: \code{TRUE} %[do].
  }
  

  \item{\code{skipchecks}}{logical.
    If \code{TRUE}, several checks whether the given parameter values
    and the dimension are within the allowed range is skipped.
    Do not change the value of this variable except you really
    know what you do.
    
    Default: \code{FALSE} $ %[also do]. 
  }
  
  \item{\code{storing}}{Logical.
    If \code{FALSE} then the intermediate results are
    destroyed after the simulation of the random field(s)
    or if an error had occured.
    If \code{storing=TRUE}, then
    additional simulations can be performed by calling
    \command{\link{RFsimulate}} with at most the argument \code{n}.
    This call can then be much faster, but the a rather large
    amount of memory could be kept.
    
    When \code{storing} turned from \code{TRUE} to \code{FALSE} by
    global call then all registers are deleted.
    Advanced:
    With \code{\link{RFoptions}(storing=list(FALSE, register,
      model_register))}
    single registers can be deleted.
    
    Default: \code{FALSE} %[do]. 
  }
  
  \item{\code{Ttriple}}{Logical or \code{NA}.
    If \code{TRUE}, then triple for the time argument \code{T} is
    expected, containing start, step (by), length.
    If \code{FALSE} a sequence on a grid is expected.
    If \code{NA} then the decision is automatic, but will lead to an
    error if ambiguous.
  }
  
  \item{\code{vdim_close_together}}{logical. Used especially in functions that
    create covariance matrices. If the model is multivariate, then two
    ways of ordering the matrix exist. To consider first all variables at
    a certain location (\code{vdim_close_together=TRUE}) or to consider first
    all locations keeping the variable fixed
    (\code{vdim_close_together=FALSE}).
    Note that several simulation methods rely on the value \code{FALSE},
    so that these methods will not work anymore if
    \code{vdim_close_together=TRUE}. 
    
    Default: \code{FALSE}.
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{2. Options for Brown-Resnick Fields}
\describe{
  \item{\code{corr_factorr}}{ to do
  }
  \item{\code{deltaAM}}{ to do
  }
  
  \item{\code{maxtrendmem}}{
    integer; the maximal number of trends for shifted locations that may
    be stored at the same time when simulating BR processes via
    \code{RPbrshifted}; if \code{maxtrendmem} is large, multiple trend
    evaluation may be avoided. 
    
    Default: \code{1e8} .
  }
  \item{\code{meshsize}}{
    positive; width of the grid on which the shape functions in the M3
    representation of BR processes are simulated; only used for
    simulation of BR processes via \code{RPbrmixed}.
    
    Default:  \code{0.1} .
  }

  \item{\code{optim_mixed}}{\code{0, 1, 2}; only used for simulation of BR
    processes via \code{RPbrmixed}.\cr
    If \code{optim_mixed=0}, the arguments
    \code{\link[=RPbrmixed]{lambda}} and 
    \code{\link[=RPbrmixed]{areamat}} of \command{\link{RPbrshifted}}
    are used for the simulation.\cr
    If \code{optim_mixed=1}, \code{\link[=RPbrmixed]{lambda}} is estimated for
    \code{\link[=RPbrmixed]{areamat=1}}.\cr
    If \code{optim_mixed=2}, \code{\link[=RPbrmixed]{areamat}} is optimized and
    \code{\link[=RPbrmixed]{lambda}} is estimated.
    
    Default: \code{1} .
  }
  \item{\code{optim_mixed_maxpoints}}{
    positive integer; only used for simulation of BR processes via
    \code{RPbrmixed} with \code{optim_mixed>0}. Maximal number of Poisson
    points used for the optimization of
    \code{\link[=RPbrmixed]{areamat}} and the estimation 
    of \code{\link[=RPbrmixed]{lambda}}.
    
    Default: \code{10000} .
  }
  \item{\code{optim_mixed_tol}}{
    value in \eqn{[0,1]}; only used for simulation of BR processes via
    \code{RPbrmixed} with \code{optim_mixed=2}. In this case,
    \code{\link[=RPbrmixed]{areamat}} is optimized under the constraint that the
    probability of drawing the shape function incorrectly is bounded by
    \code{optim_mixed_tol} (cf. Oesting et al., 2012).
    
    Default: \code{0.01} .
  }
  
  \item{\code{variobound}}{
    positive; the shape functions in the mixed moving maxima
    representation are cut off where the variogram belonging
    to \code{phi} exceeds \code{variobound}.
    
    Default: \code{8.0} .
  }
  
  \item{\code{vertnumber}}{
    positive integer; for an efficient simulation of the shape functions
    in the M3 representation of BR processes, the component \eqn{E} from
    of the domain \eqn{[x_0, \infty] \times E}{[x_0, Inf] x E} of the
    underlying Poisson point process is sub-dividedinto cubes
    (cf. Oesting et al., 2012); \code{vertical}  is the number of
    vertical breaks of \eqn{E}; only used for simulation of BR processes
    via \code{RPbrmixed} with \code{optim_mixed=2}.
    
    Default: \code{7} .
  }
   
}


 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{3. \code{circulant}: Options for circulant embedding methods, cf. \command{\link{RPcirculant}}}\cr
   These options influence the standard circulant embedding
   method, cutoff circulant embedding  intrinsic circulant embedding.
   It can also influence \command{\link{RPtbm}} if the line is simulated
   with any circulant embedding method.
   \cr
 \describe{
 \item{\code{approx_maxgrid}}{See \command{\link{RPcirculant}}}
 \item{\code{approx_step}}{See \command{\link{RPcirculant}}}
 \item{\code{dependent}}{See \command{\link{RPcirculant}}}
 \item{\code{force}}{See \command{\link{RPcirculant}}}
 \item{\code{maxGB}}{See \command{\link{RPcirculant}}}
 \item{\code{maxmem}}{See \command{\link{RPcirculant}}}
 \item{\code{mmin}}{See \command{\link{RPcirculant}}}
 \item{\code{strategy}}{See \command{\link{RPcirculant}}}
 \item{\code{tolIm}}{See \command{\link{RPcirculant}}}
 \item{\code{tolRe}}{See \command{\link{RPcirculant}}}
 \item{\code{trials}}{See \command{\link{RPcirculant}}}
 \item{\code{useprimes}}{See \command{\link{RPcirculant}}}
 }
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{4. \code{coords}: Options for coordinates and units}
\describe{
  \item{\code{coord_system}}{character. See \link{coordinate systems}
  }
  
  \item{\code{coordunits}}{ See \link{coordinate systems}
  }
 
  \item{\code{coordnames}}{See \link{coordinate systems}
  }

  
  \item{\code{new_coord_system}}{See \link{coordinate systems}
  }

  \item{\code{new_coordunits}}{See \link{coordinate systems}
  }

  \item{\code{polar_coord} }{See \link{coordinate systems}
  }

  \item{\code{varnames}}{See \link{coordinate systems}
  }  

  \item{\code{varunits}}{See \link{coordinate systems}
  }
    
  \item{\code{xyz_notation}}{See \link{coordinate systems}
  }
 
  \item{\code{zenit}}{See \link{coordinate systems}
  }
  
}





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{5. \code{direct}: Options for simulating by simple matrix decomposition}
\describe{
  \item{\code{max_variab}}{See \command{\link{RPdirect}}}
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{6. \code{distr}: Options for distributions, in particular \command{\link{RRrectangular}}}
\describe{
  \item{\code{innermin}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.
    The minimal length of the interval where the Taylor expansion shall
    be valid.
   
    Default: \code{1e-20} .
  }

  \item{\code{maxit}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.

    The number of iterative steps where the
    the constant of the Taylor development is increased,
    to find an upper bound for the given function.
    
    Default: \code{20} .
  }

  \item{\code{maxsteps}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.

    \code{maxsteps} is usually the number of steps in the middle part of
    the approximation. From this value and the length between
    the determined endpoints for  the approximation at the origin and in
    the tail, the step length is calculated. If the step length is less
    than \code{minsteplen} the number of steps is reduced.
    
    Default: \code{1000} .
  }

  
  \item{\code{mcmc_n}}{
    In case of the use of MCMC it leaves out \eqn{n-1}
    member of the Markov chain bevor the \eqn{n} member
    is returned. See also maxsteps.
    
    Default: \code{15} .
  }

  
  \item{\code{minsteplen}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.
    The minimal step length
    for the middle part of approximation, which is a step function,
    
    Default: \code{0} (i.e. not used as a criterion.)
  }
  
  \item{\code{outermax}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.
    The largest possible endpoint for the middle part that
    approximates the function by a step function. See also \code{innermax}.

    Default: 20.
  }
 
  \item{\code{parts}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.

    \code{parts} determines the number of tests that are performed to
    check whether a proposed power function is an upper bound for
    the given function, at the origin and the tail.
    
    Default: \code{8} .
  }

  \item{\code{repetitions}}{
    Minimal number of realisations to determine a quantity of the
    distribution by MCMC. E.g. to determine the integral value \eqn{c}
    in the paper of Oesting, Schlather, Zhou.
    
    Default: 1000.
  }
  
  \item{\code{safety}}{ 
     Default value to simulate from the
     \command{\link{RRrectangular}} distribution.

     First, at the origin, the first power function of the Taylor
     expansion is taken as potential upper function.
     The constant of the power function are increased by factor
     \eqn{1 + }\code{safety} and the exponent of the function
     similarly decreased. A number of test evaluations
     is performed to check whether this modified function is indeed
     a upper bound. If not, the considered interval at the origin
     is reduced iteratively, the constants of the power function
     further increased and the exponent decreased.
     If \code{maxit} iteration have been performed without success,
     the search for an upper bound fails.
     The search at the origin also fails if the interval around
     the origin has become less than \code{innermin}.
     

     Similar procedure is performed for the tail.
     
   
    Default: \code{0.08} .
  }
  
 
  
}


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{7. \code{empvario}: Options for calculating the empirical variogram}
\describe{ 
  \item{\code{fft}}{
    Logical. Determines whether FFT should be used for data on a grid
    Default: \code{TRUE}.
  }
  \item{\code{phi0}}{
    numeric. In case of anisotropic fields directional cones are
    considered. The argument \code{phi0} determines the starting angle.
    
    Default: \code{0}.
  }
  \item{\code{pseudovariogram}}{
    logical. Only in the multivariate case. Whether the
    pseudovariogram or the crossvariogram should be calculated.
    
    Default: \code{FALSE}.
  }
  \item{\code{theta0}}{
    numeric. In case of anisotropic fields directional cones are
    considered. The argument \code{theta0} determines one of the
    boundaries, hence all boundaries for a given fixed number of cones.
    The argument \code{theta0} determines the starting value of the
    second anglue in polar coordinate representation in 3 dimensions.
     
    Default: \code{0}.
  }
  \item{\code{tol0}}{
    numeric. Estimated values of the empirical variogram
    below \code{tol0} times the grid step in the third dimension
    are considered to be zero. Hence the respective values are set
    to zero.
    
    Default: \code{1e-13}.
  }
 }
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{8. \code{fit}: Options for \command{\link{RFfit}},
   \command{\link{RFratiotest}}, and \command{\link{RFcrossvalidate}}}

\describe{
  \item{\code{algorithm}}{
    See \link{RFfitOptimiser}.

    Default: \code{NULL}
  }
  
  \item{\code{approximate_functioncalls}}{
    In case the parameter vector is too close to the given
    bounds, the ML target function is evaluated on a grid
    to get a new initial value for the ML estimation. 
    The number of points of the grid is approximately
    \code{approximate_functioncalls}.
 
    Default: \code{50}
  }

  \item{\code{boxcox_lb}}{
    lower bound for the Box-Cox transformation
    
    Default: \code{-10}.
  }
  \item{\code{boxcox_ub}}{
    upper bound for the Box-Cox transformation
    
    Default: \code{10}.
  }

  \item{\code{bin_dist_factor}}{
    numeric. The empirical variogram is calculated up the distance
    \code{bin_dist_factor} times (maximum distance among any pair of locations)
    
    Default: \code{0.5}.
  }
  
  \item{\code{bins}}{vector of explicit boundaries for the bins or the
    number of bins for the empirical variogram (used in the
    LSQ target function, which is described at the beginning
    of the Details).
    Note that for anisotropic models, the value of \code{bins} might
    be enlarged.
    
    Default: \code{20}.
  }

 \item{\code{critical}}{logical or signed integer.
    
    If \code{critical=FALSE} and if the result of
    any maximum likelihood method 
    is on a borderline, then the optimisation is redone
    in a modified way (which takes about double extra time)
    
    If \code{critical=TRUE} and if the result of
    any maximum likelihood method
    is on a borderline, then a kind of profile likelihood
    optimization is done (which takes about 10 times extra time)
    
    If \code{critical>=2} then a kind of profile likelihood
    optimization is always done (which takes about \code{n_crit}
    times extra time) for an automatically chosen selection
    of the model parameters.
    
    If \code{critical>=3} then a kind of profile likelihood
    optimization is always done (which takes about \code{n_crit}
    times extra time) for all the parameters.
    
    If \code{critical<0} then none of the refined methods
    are performed.
    
    Default: \code{TRUE}.
  }

  \item{\code{cross_refit}}{logical. 
    For each of the subset of the cross-validation method
    the parameters have to be fitted to the given model.
    If \code{cross_refit} is \code{TRUE}, this is done, but takes a huge
    amount of time. If \code{FALSE}, the model is fitted only once to
    the data and the value at
    each point is predicted with the same model given
    the values of the other points.
    
    Default: \code{FALSE}.
  }

  \item{estimate_variance}{
    see \command{\link{RFlikelihood}}.
  }

  \item{factr, factr_recall}{
    See the argument \code{control} in \link[stats]{optim}.
    \code{factr_recall} is used for intermediate calculations.
  }
  
  \item{likelihood}{character -- not programmed yet.
    types of likelihood are \code{"auto"}, \code{"full"},
    \code{"composite"}, \code{"tesselation"}; 
    
    Default: \code{"auto"}    
  }

  \item{\code{lowerbound_scale_factor}}{
    The lower bound for the scale is determined as
    
    (minimum distance between different pairs of points) /\cr
     \code{lowerbound_scale_factor}.
    
    Default: \code{3}.
  }

  \item{\code{lowerbound_scale_ls_factor}}{ For the LSQ target
    function a different lower bound 
    for the scale is used. It is determined as
    
    (minimum distance between different pairs of points) / \cr
    \code{lowerbound_scale_ls_factor}.
    
    Default: \code{5}.
  }

 % \item{\code{lowerbound_sill}}{absolute lower bound for variance
%    and nugget. See \code{lowerbound_var_factor}.
%    
%    Default: \code{1E-10}.
%  }
  
  \item{\code{lowerbound_var_factor}}{
    The lower bound for the nugget and the variance is determined
    as var(\code{data}) / \code{lowerbound_var_factor}.
    If a standard model definition is given and
    either the nugget or the variance is fixed,
    the parameter to be estimated
    must also be greater than \code{lowerbound_sill}.
    
    Default: \code{10000}.
  }

  \item{\code{maxmixedvar}}{
    upper bound for variance in a mixed model;
    so, the covariance model for mixed model part might
    be calibrated appropriately
  }

  \item{\code{max_neighbours}}{integer.
    Maximum number of locations (with depending values)
    that are allowed.
    
    Default: \code{5000}.
  }
  
  \item{\code{minbounddistance}}{
    If any value of the parameter vector
    returned from the ML estimation
    is closer than \code{minbounddistance}
    to any of the bounds or if any value
    has a relative distance smaller than
    \code{minboundreldist}, then it is assumed that
    the MLE algorithm has dropped into a local minimum,
    and it will be continued with evaluating the
    ML target function on a grid, cf. the beginning paragraphs
    of the Details.
    
    Default: \code{0.001}.
  }
  
  \item{\code{minboundreldist}}{relative distance to the bounds
    below which a part of the algorithm is considered as
    having failed. See \code{minbounddistance}.
    
    Default: \code{0.02}.
  }

  \item{\code{min_diag}}{
    Minimal value of any estimated diagonal matrix element.    
    
    Default: \code{1e-7}.
  }

  \item{\code{n_crit}}{integer.
    The approximate profiles that are considered.
    
    Default: \code{10}.
  }

    \item{\code{nphi}}{scalar or vector of 2 components.
    If it is a vector then the first component gives the first angle
    of the xy plane
    and the second one gives the number of directions on the half circle.
    If scalar then the first angle is assumed to be zero.
    Note that a good estimation of the variogramm by LSQ with a
    anisotropic model a large value for \code{ntheta} might be needed
    (about 20).
    
    Default: \code{1}. 
  }
  \item{\code{ntheta}}{scalar or vector of 2 components.
    If it is a vector then the first component gives the first angle
    in the third direction
    and the second one gives the number of directions on the half circle.
    If scalar then the first angle is assumed to be zero.
    
    Note that a good estimation of the variogramm by LSQ with a
    anisotropic model a large value for \code{ntheta} might be needed
    (about 20). 
    
    Default: \code{1}. 
  }
  \item{\code{ntime}}{scalar or vector of 2 components.
    if \code{ntimes} is a vector, then the first component are the
    maximum time distance (in units of the grid length \code{T[3]}) and the
    second component gives the step size (in units of the grid length
    \code{T[3]}). If scalar then the step size is assumed to 1 (in units
    of the grid length \code{T[3]}).
    
    Default: \code{20}. 
  }

  \item{\code{only_users}}{boolean.
    If true then only \code{users_guess} is used as a
    starting point for the fitting algorithms
    
    Default: \code{FALSE}. 
  }

  \item{\code{optimiser}}{
    See \link{RFfitOptimiser}.
    
    Default: \code{"optim"}.
  }
  
% \item{\code{optim_var_elimination}}{This argument takes the values
%    \code{'never'}, \code{'respect bound'}, \code{'try'},
%    \code{'yes'}, and should only be 
%    set by the advanced user. Background of this option is that
%    a global variance can optimized analytically.
%    
%    The meaning of the values is as follows.
%    \itemize{
%      \item{\code{'never'}}{
%	A global variance is never tried to be eliminated
%      }
%      \item{\code{'respect bound'}}{
%	A global variance is eliminated if such a variance
%	is detected and the user did not indicate bounds
%	for the parameters.
%      }
%      \item{\code{'try'}}{
%	A global variance is eliminated if such a variance
%	is detected.
%      }
%      \item{\code{'yes'}}{
%	A global variance is tried to be eliminated altough
%	the algorithm did not find an indication. Here, the full
%	responsibility is left to the user. (This option might
%	make sense if \code{transform} is given.)
%	This option is only overwritten when it does not make sense,
%	e.g. no variance is estimated.
%      }
%    }
%    Default: \code{'respect bound'}.
%  }

  \item{pgtol, pgtol_recall}{
    See the argument \code{control} in \link[stats]{optim}.
    \code{pgtol_recall} is used for intermediate calculations.
  }
  
  \item{\code{refine_onborder}}{logical.
    If \code{TRUE} and an estimated parameter of the model
    is close to the boundary, a second search for the optimum
    is started.

    Default: \code{TRUE}
    
  }
 
  \item{\code{minmixedvar}}{
    lower bound for variance in a mixed model;
    so, the covariance model for mixed model part might
    be calibrated appropriately

    Default:  1/1000
  }
  
%  \item{\code{solvesigma}}{Logical. -- experimental stage!
%    If a mixed effect part is present where the variance
%    has to be estimated, then this variance parameter is solved
%    iteratively within the profile likelihood function, if
%    \code{solvesigma=TRUE}.This makes sense
%    if the number of independent variables is very small.
%    If \code{solvesigma=FALSE} then the variance parameter is
%    treated as any other parameter to be estimated.
%    Default: \code{FALSE}.
%  }
  
  \item{\code{ratiotest_approx}}{logical.
    if \code{TRUE} the approximative formula that twice the
    difference of the likelihoods follow about a \eqn{\chi^2}
    distribution is used. The parameter of freedom equals
    the number of parameters to be estimated for the covariance
    function, including those for the covariates.

    Default: \code{TRUE}
  }
  \item{\code{reoptimise}}{logical.
    If \code{TRUE && !only_users} then at a very last step,
    the optimisation is redone with currently best parameters
    and likelihood as scale parameter for \command{\link{optim}}.
    
    Default: \code{TRUE}. 
  }
  

  \item{\code{scale_max_relative_factor}}{ If the initial scale
    value for the ML estimation 
    obtained by the LSQ target function is
    less than
    \eqn{(minimum distance
      between different pairs of points) / }
    \code{scale_max_relative_factor}
 
    a warning is given that probably a nugget effect
    is present. 
    Note: if \code{scale_max_relative_factor} is greater
    than \code{lowerbound_scale_ls_factor} then
    no warning is given as
    the scale has the lower bound \eqn{(minimum distance
      between different pairs of points) / }
    \code{lowerbound_scale_ls_factor}.
    
    Default: \code{1000}
  }

  \item{\code{scale_ratio}}{
    \command{\link{RFfit}} uses \code{parscale} and \code{fnscale}
    in the calls of \command{\link{optim}}. As these arguments should
    have the magnitude of the estimated values, \command{\link{RFfit}}
    checks this by calculating the absolute log ratios.
    If they are larger than \code{scale_ratio},
    \code{parscale} and \code{fnscale} are reset and the optimisation
    is redone. 
    
    Default: \code{0.1}. 
  }
   
  \item{\code{shortnamelength}}{
    The names of the variables in the returned table are
    abbreviated by taking the first \code{shortnamelength}
    letters.
    
    Default: \code{4}. 
  }
 
%  \item{\code{sill}}{ currently not maintained anymore.    
%    Additionally to estimating \code{nugget} and \code{variance}
%    separately, they may also be estimated together under the
%    condition that \code{nugget} + \code{variance} = \code{sill}.
%    For the latter a finite value for \code{sill} has to be supplied,
%    and \code{nugget} and \code{variance} are set to \code{NA}.
%    
%    \code{sill} is only used for the standard model. 
%    
%    Default: \code{NA}. 
%  }

  \item{\code{smalldataset}}{
    If the number of locations is considered as small, then some more data
    are kept in the storage to accelerate the estimation algorithm.
    
    Default: \code{2000}.
  }
  
  \item{\code{split}}{integer.
    If the number of parameters to be numerically optimised is larger
    than or equal to \code{split} then \command{\link{RFfit}} checks whether a
    space-time covariance model or a multivariate covariance model
    can be split into components, so that certain parameters
    can be estimated separately.
    
    Default: \code{4}. 
  }
   
  \item{\code{cliquesize}}{integer.
    \command{\link{RFfit}} tries to split the data set
    into parts of size splitn_neighbours[2] or less, but never more than 
    \code{splitn_neighbours[3]} and never less than
    splitn_neighbours[1].     
    
    Default: \code{c(200, 1000, 3000)}.
  }
  
  \item{\code{splitfactor_neighbours}}{
     The total number of neighbouring boxes in each direction
     \eqn{1 + 2\code{splitfactor}}, including the current box itself.
    
    Default: \code{2}.
  }
  \item{\code{split_refined}}{logical.
    If \code{TRUE} then also submodels are fitted if splitted.
    This takes more time, but \command{\link[=anova.RF_fit]{anova}} and
    \command{\link{RFratiotest}}, for instance,
    will give additional information.

    Default: \code{TRUE}.
  }
  
  \item{\code{upperbound_scale_factor}}{
    The upper bound for the scale is determined
    as

    \code{upperbound_scale_factor} * (maximum distance
      between all pairs of points).
    
    Default: \code{3}.
  }
  \item{\code{upperbound_var_factor}}{ The upper bound for the
    variance and the nugget is determined 
    as \code{upperbound_var_factor} * var(\code{data})
    
    Default: \code{10}.
  }

  \item{\code{use_naturalscaling}}{
    logical. Only used if model is given in standard (simple) way.
    If \code{TRUE} then \emph{internally}, rescaled
    covariance functions will be used for which
    cov(1)\eqn{\approx}{~=}0.05.
    \code{use_naturalscaling} has the advantage that \code{scale}
    and the form parameters of the model get \sQuote{orthogonal},
    but \code{use_naturalscaling} does not work for all models.
    
    Note that this argument does not influence
    the output of \command{\link{RFfit}}: the parameter vector
    returned by \command{\link{RFfit}} refers
    \emph{always} to the standard covariance model as given in
    \command{\link{RMmodel}}. (In contrast to \code{practicalrange}
    in \command{\link{RFoptions}}.)\cr
    Advantages if \code{use_naturalscaling=TRUE}:
    
    \itemize{
      \item \code{scale} and the shape parameter of a parameterised
      covariance model can be estimated better if they are estimated
      simultaneously.
      \item The estimated bounds calculated by means of
      \code{upperbound_scale_factor} and \code{lowerbound_scale_factor},
      etc. might be more realistic.
      \item in case of anisotropic models, the inverse of the elements
      of the anisotropy matrix should be in the above bounds.
    }
    Disadvantages if \code{use_naturalscaling=TRUE}:
    \itemize{
      \item For some covariance models with additional parameters, the
      rescaling factor has to be determined numerically.
      Then, more time is needed to perform \command{\link{RFfit}}.
    }
    
    Default: \code{TRUE}.
  }

}


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{9. \code{gauss}: Options for simulating Gaussian random fields}
\describe{
  \item{\code{approx_zero}}{
    Value below which a correlation is considered to be essentially zero.
    This argument is used to determine the practical range of covariance
    function with non-compact support.
    
    Default: \code{0.05}
  }
 
 \item{\code{boxcox}}{
   real vector of one or two components.
   If the first component is \code{Inf} then no transformation
   is performed. Otherwise the BoxCox transformation is performed.
   Note that Box Cox only works in a Gaussian framework.
   Note further that either \code{boxcox} or \code{loggauss}
   may be given.
  
   Default \code{c(Inf, 0)}
  }


  \item{\code{direct_bestvar}}{integer.
    When searching for an appropriate simuation method
    the matrix decomposition method (\code{method="direct"})
    is preferred if the number of variables is less than or equal to 
    \code{direct_bestvariables}.
    
    Default is \code{1200}.
  }

  \item{\code{loggauss}}{
    logical. Whether a log-Gauss random fields should be returned.
    See also \code{boxcox} for a generalisation.
  }

  \item{\code{paired}}{
    (\dQuote{Antithetic pairs}.)
    Logical. If \code{TRUE} then the second half of the
    simulations is logical. If \code{TRUE} then the second half of the
    simulations is obtained by
    only changing the signs of all the standard Gaussian random variables,
    on which the first half of the 
    simulations is based. Default is \code{FALSE}.
  }
  
  \item{\code{stationary_only}}{
    See \command{\link{RPgauss}}
  }
 
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{10. \code{graphics}: Options for graphical output}
\describe{
  \item{\code{always_close_screen}}{logical.
    If \code{FALSE} the current device is kept as it is.
    Otherwise the action depends on the value of
    \code{height}: if \code{height} is not positive
    then \command{\link[graphics]{close.screen}}
    is performed on the current device.
    Else the current device is closed.
    
    Default: \code{FALSE}.
  }
  \item{\code{always_open_screen}}{logical.
    if \code{TRUE} a new graphical window is opened for every
    \command{\link[graphics]{plot}} if a standard graphical output is
    used. If \code{NA} then the value is set
    to \code{\link[base]{interactive}()}. %ok
    
    Default: \code{TRUE}.
  }
  \item{\code{grPrintlevel}}{
    integer values 0, 1, 2. The higher the more text is shown in the
    plot.
    
    Default: \code{1}.
  }
  
  \item{\code{height}}{real number.
    if \code{height} is greater than zero then it gives the height
    of a single figure in a plot create by \pkg{RandomFields};
    always a new window is opened. See also \code{always_close_screen}.
    
    If plots with multiple figures are
    shown, the height and width of the plot
    will be increased by a factor up the
    ones given by \code{increase_upto}.
       
    Default: \code{6}.
  }

  \code{pdffile}{ argument \code{file} in \command{\link{pdf}}
    If \code{""} then no internal naming is performed.

    Default: \code{""}.
  }

  \item{\code{jpegfile}}{ argument \code{file} in \command{\link{jpeg}}
    If \code{""} then no internal naming is performed.

    Default: \code{""}.
  }

  \item{\code{filenumber}}{integer. Starting number of the file if
    \code{onefile=FALSE}. It is set to 0 whenever \code{pdffile} is
    changed and \code{onefile=FALSE}. 

    Default 0.
  }
  
  \item{\code{onefile}}{about the behaviour of argument \code{onefile} in
    \command{\link{pdf}}

    Default: \code{FALSE}.
  }

  
  \item{\code{increase_upto}}{
    See \code{height}.
    
    Default: \code{c(3,4)}.
  }
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{11. \code{gui}: Options for c\command{\link{RFgui}}}
\describe{
  \item{\code{alwaysSimulate}}{
    logical. If \code{TRUE} then a new random field is simulated
    whenever a parameter is changed. Otherwise only the covariance
    function or the variogram is re-plotted; simulations are performed
    only when the correponding button is pressed.
    
    Default: \code{TRUE}.
  }
  \item{\code{simu_method}}{
    \code{"RPcirculant"},
    \code{"RPcutoff"},
    \code{"RPintrinsic"},
    \code{"RPtbm"}, 
    \code{"RPspectral"},
    \code{"RPdirect"},
    \code{"RPsequential"},
    \code{"RPaverage"},
    \code{"RPnugget"},
    \code{"RPcoins"},
    \code{"RPhyperplane"},
    \code{"RPspecific"},
    \code{"any method"}.
    
				      
    Default: \code{"RPcirculant"}.
  }
  \item{\code{size}}{vector of 2 components.
    Grid size of the simulated stochastic processes.
    The two components of the vector correspond to one-dimensional and
    two-dimensional processes, respectively.
    
    Default: \code{c(1024, 64)}.
  }
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{12. \code{hyper}: Options for simulating hyperplane tessellations}
\describe{
  \item{\code{mar_distr}}{integer.
    This argument should not be changed yet.
    
    It  codes the marginal distribution used in the
    simulation:
    
    0 : uniform distribution\cr
    1 : Frechet distribution with form argument \code{mar.param}\cr
    2 : Bernoulli distribution (Binomial with \eqn{n=1}) with
    argument \code{mar.param}
    
    
    Default: \code{0} .
  }

  \item{\code{mar_param}}{Argument used for the marginal
    distribution. The argument should not be changed yet.
    
    Default: \code{NA} .
  }

  \item{\code{maxlines}}{integer.
    Maximum number of allowed lines.
    
    Default: \code{1000} .
  }
  
  \item{\code{superpos}}{integer.
    number of superposed hyperplane tessellations.
    
    Default: \code{300} .
  }
}



 
\bold{13. \code{krige}: Options for Kriging}
\describe{
  
  \item{\code{cholesky_R}}{ obsolete }

  \item{\code{fillall}}{
    logical value for imputing.
    If true all the components are estimated whether they are
    \code{NA} or not.
    
    Default: \code{TRUE}.
  }
  
  \item{\code{locmaxn}}{
    Kriging is conditions on maximal \code{locmaxn} points.
    If the data contain more points, neighbourhood kriging is performed.
    
    Default: \code{8000}. 
  }
  
  \item{\code{locsplitfactor}}{
    In case of neighbourhood kriging, the area is split into small
    boxes. The complete neighbourhood contains (2 *
    \code{locsplitfactor} +1) boxes in each direction.
    
    
    Default: \code{2}. 
  }
  
%  \item{\code{locsplitn}}{vector of 3 components.
%    A box should contain no more than \code{locsplitn[1]}
%    points, but never less than \code{locsplitn[2]}. If
%    a box had originally less than \code{locsplitn[2]} points,
%    then the box is increased until at least \code{locsplitn[3]}
%    points are in the box.    
   \item{\code{locsplitn}}{vector of 3 components.
    A box should contain no more than \code{locsplitn[3]}
    points, but never less than \code{locsplitn[1]}. If
    a box had originally less than \code{locsplitn[1]} points,
    then the box is increased until at least \code{locsplitn[2]}
    points are in the box.
   
    Default: \code{c(200, 1000, 5000)}. 
  }
  
  \item{\code{method}}{ obsolete }
  
  \item{\code{return.variance}}{logical.
    If \code{FALSE} the kriged field is
    returned. If \code{TRUE} a list of two elements, \code{estim} and
    \code{var}, i.e. the kriged field and the kriging variances,
    is returned.
    
    Default: \code{FALSE}. 
  }
  
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{14. \code{maxstable}: Options for simulating max-stable random fields}
\describe{
  \item{\code{check_every}}{
    integer. In order to get a precise simulation result, by definition,
    the maximum must be taken, for each shape function, over alle
    locations of interest. Clearly, small values will not play a role.
    To this end, the global minimum has to be determined.
    The calculation of the global minimum is expensive and therefor
    should not be  done too frequently. On the other hand,
    rare updates increases the computing times for taking the maximum
    over a single shape functions. Here, after every \code{check_every}
    considered shape function, the global minimum is calculated.
    It is expected that a good choice for  \code{check_every} is in
    in the interval \eqn{[10, 100]}.
     
    (For ease and for concerns of efficiency, the more adequate, local
    minimum is not considered.)
    
    Default: 30 .
  }
    
  \item{\code{density_ratio}}{
    value in \eqn{[0,1]}. This argument is considered only
    if \code{flat=-1} and the simulation is performed on a grid.
    Then, the ratio between the highest and the lowest value is
    calculated within the convex hull of the grid. If the
    value is less than \code{density_ratio} then the grid points
    are considered separately. Else the density is considered to be
    constant in the convex hull of the grid.
    
    Default: 0.0.
  }
  
  \item{\code{eps_zhou}}{positive real number, which 
    gives the aimed relative precision.
    E.g. if \code{eps_zhou=0.01} then the first 2 digits should be
    correct.
   
    Default: 0.01
  }

  \item{\code{flat}}{\code{-1, FALSE, TRUE}.
    The argument is considered only if the simulation is performed on a
    grid. 
    If \code{flat} is logical, then the density
    is considered to \code{flat} in the convex hull of the grid.
    If \code{flat=-1} the choice is done automatically.
    
    Default: -1 .
  }

  \item{\code{max_gauss}}{
    The simulation of the max-stable process based on random fields uses
    a stopping rule that necessarily needs a finite upper endpoint
    of the marginal distribution of the random field.
    In the case of
    \command{\link[=RPbrownresnick]{Brown-Resnick processes}},
    \command{\link[=RPschlather]{extremal Gaussian fields}},
    and
    \command{\link[=RPopitz]{extremal t fields}},
    the upper endpoint is approximated by \code{standardmax}. 
    
    Default: \code{3.0} .
  }
  
  \item{\code{max_n_zhou}}{positive integer.
    The overall constant \eqn{c} in the paper of
    Oesting, Schlather, Zhou (2014) has to be determined
    by MCMC, if the shape functions are random.
    
    The two arguments, \code{min_n_zhou} and \code{max_n_zhou},
    give the minimal and the maximal
    number of simulations that are performed. To economize
    computer time the values of \eqn{c} is partially estimated 
    when the shape functions are simulated. If the number
    of shape functions is larger than the number of simulations
    given by \code{eps_zhou} then
    no further simulation is performed to determine \eqn{c}.
   
    Default:  1000 and 10000000, respectively.
  }

  \item{\code{maxpoints}}{
    positive integer; the maximal number of Poisson points to be simulated
    for one realization of the max-stable random field. This option will
    not be considered for most of the users.
     
    Default: \code{2e9} .
  }

  \item{\code{mcmc_zhou}}{positive integer.
    In case of random shape functions, an MCMC step is required.
    \code{mcmc_zhou}-1 equals the number of members of the MCMC chain
    that are left out before the next value of the chain is returned.
    
    Default: 20
  }
  
  \item{\code{min_n_zhou}}{see \code{max_n_zhou}}

  \item{\code{xi}}{
    Extreme value index.    Default: \code{2e9} .
  }

  \item{\code{mcmc_zhou}}{positive integer.
    In case of random shape functions, an MCMC step is required.
    \code{mcmc_zhou}-1 equals the number of members of the MCMC chain
    that are left out before the next value of the chain is returned.
    
    Default: 20
  }
  
  \item{\code{min_n_zhou}}{see \code{max_n_zhou}}

  \item{\code{xi}}{
    Extreme value index.
    
    Default: \code{1.0} .
  }
}
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{15. \code{mpp}: Options for the random coins (shot noise) methods}
\describe{
\item{\code{about_zero}}{
    In certain cases (\link{Coins},\link{RMtruncsupport}),
    functions are assumed to zero if the value is less than \code{about_zero}.
    
    Default: \code{0.001} .
  }
  \item{\code{n_estim_E}}{integer. Number of draws from the
    distribution of the scale to estimate the mean of the distribution.
    This is used only if the mean of the scale distribution
    is not explicitely given.
    Default: \code{50000} .
  }
  
    
 
  \item{\code{scatter_size}, \code{scatter_max}}{
    Used in function \code{RMscatter} that calculates
    \eqn{\sum_{i=1}^n f(x + h_i)} for some function \eqn{f} and
    for some distances \eqn{h_i}.
      
    %    see \code{size} and \code{max} in \code{\link{RMscatter}}.
    Real valued and integer valued, respectively, or \code{NA}.
    Let \eqn{\varepsilon=}\code{about_zero}, \eqn{s=}\code{scatter_size} and \eqn{m=}\code{scatter_max}.
    We distinguish 4 cases:
    \itemize{
      \item \code{scatter_size > 0} and \code{scatter_max >= 0}\cr
      Here, \eqn{n} equals \eqn{(2m)^d}.
      and \eqn{h_i \in M = \{ (k s, \ldots, k s),\ldots, (m s, \ldots, m
	s)\}}
      with \eqn{k=-m}.
      
      \item \code{scatter_size > 0} and \code{scatter_max < 0}\cr
       same as the previous case, but \eqn{m} is chosen such that
      \eqn{f(k_i e_i s_i) \approx \varepsilon}, \eqn{-k_i\in N},
      \eqn{i=1,\ldots,d} and 
      \eqn{f(m_i e_i s_i) \approx \varepsilon}, \eqn{m \in N}.
      
      \item \code{scatter_size <= 0} and \code{scatter_max >= 0}\cr
      This option is possible only for grids.
      Here \eqn{h_i} runs on the given grid \eqn{i=1,\ldots,d},
      but at most \code{scatter_max} steps.
      
      \item \code{scatter_size <= 0} and \code{scatter_max < 0}\cr
      this option is possible only for grids.
      Here, \eqn{h_i} runs over the whole grid.
    }
  }
  \item{\code{shape_power}}{
    Shape functions are powered by \code{shape_power} before used as
    intensity function for the point process.
    
    Default: \code{2.0}.
  }
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{16. \code{nugget}: Options for the nugget effect}\cr
Simulating a nugget effect is per se trivial.
However, it gets complicated
and best methods (including \code{direct} and \code{circulant
  embedding}!) fail if zonal anisotropies are considered,
where sets of points have to be identified that belong to the
same subspace of eigenvalue 0 of the anisotropy matrix.
\describe{
  \item{\code{tol}}{
    The nugget tolerance influences two different kind of models
    \itemize{
      \item \command{\link{RPnugget}}
      \item \command{\link{R.is}}
    }
    See there for more information.
  }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{17. \code{registers}: Register numbers}\cr
Model for different purposes are or can be stored at
different places. They are called registers and have non-negative
numbers up to 21 (currently).
The user can use the registers 0..9.
\describe{
  \item{\code{register}}{number in 0:9; place where intermediate calculation
    for random field simulation are stored; 
    the number refers to 10 internal registers 0..9.
    
    Changing the register number only makes sense, when
    two different random fields, say, are to be simulated
    alternatingly, several times in a row. Then the
    simlulation speed can be increased if several registers
    are used, \code{storing=TRUE} and \command{\link{RFsimulate}}
    is used with the only argument \code{n}.
     
    Default: \code{0} %[also do].
  } 
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{18. \code{sequ}: Options for the sequential method}
\describe{
  \item{\code{back_steps}}{See \command{\link{RPsequential}}}
  \item{\code{initial}}{See \command{\link{RPsequential}}}
  \item{\code{max_variables}}{See \command{\link{RPsequential}}}
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{19. \code{solve}: Options for solving linear systems}

\describe{ % to do
  \item{\code{max_chol}}{integer. Maximum number of rows of a matrix in
    a Cholesky decomposition

    Default: \eqn{8192}
  }
  \item{\code{max_svn}}{integer. Maximum number of rows of a matrix in
    a svd decomposition

    Default: \eqn{6555}
  }

   \item{\code{matrix_methods}}{
    vector of at most 3 integers that gives the sequence of methods
    in order to inverse a matrix or to calculate its square root: 

    0 : Choleskey decomposition \cr
    1 : SVD\cr
    2 : spam (sparse matrix algorithm)\cr
    3 : OR\cr
    4 : LU\cr
    5 : <none>

    Note that if \code{use_spam} is not \code{false} the algorithm
    checks whether a sparse matrix algorithm should be used and which is
    then tried first.

    Values larger than 4 are used internally:\cr
    5 : no further method available\cr
    6 : not initialised\cr
    7 : diagonal matrix found\cr
    
    Default: \code{5}.
  }
  \item{\code{spam_factor}}{
    integer. See argument \code{spam_sample_n}.

    Default: 4294967
  }

  
  \item{\code{spam_min_n}}{
    integer. Has the matrix 
    
    Default: 400
  }
  \item{\code{spam_min_p}}{
    number in \eqn{(0,1)} giving the proportion of
    zero about which an sparse matrix algorithm is used.

    Default: \code{0.8}
  }
  \item{\code{spam_pivot}}{
    integer. Pivoting algorithm for sparse matrices:
    0:none; 1:MMD, 2:RCM

    See package \code{spam} for details.

    Default: 1
  }
  \item{\code{spam_sample_n}}{
    Whether a matrix is sparse or not is tested by a
    \sQuote{random} sample of size \code{spam_sample_n};
    The selection of the sample is iteratively
    obtained by multiplying the index by \code{spam_factor}
    modulo the size of the matrix.

    Default: 500. 
  }
  \item{\code{spam_tol}}{
    largest absolute value being considered as zero.
    Default: \code{DBL_EPSILON}
  }
  
 \item{\code{svdtol}}{
    When the svd decomposition is used for calculating the square root of
    a matrix then the absolute componentwise difference between
    this matrix and the square of the square root must be less
    than \code{svdtol}. No check is performed if
    \code{svdtol} is negative.
 
    When the svd decomposition is used for calculating the inverse of
    a matrix then a diagonal value is set to zero if it is less than
    \code{svdtol}.
    
    Default: \code{1e-8}
  } \item{\code{use_spam}}{
    Should the package \code{spam} (sparse matrices)
    be used for matrix calculations?
    If \code{TRUE} \pkg{spam} is always used. If \code{FALSE},
    it is never used. If \code{NA} its use is determined by
    the size and the sparsity of the matrix.
    
    Default: \code{NA}. 
  }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{20. \code{special}: Options for specific methods}
\describe{
  \item{\code{multicopies}}{Only used by \command{\link{RMmult}}.
    The covariance functions are multiplied if the corresponding
    independent random fields are multiplied. To get
    an approximative Gaussian random fields with a multiplicative
    covariance functions the average over \code{multicopies}
    products of random fields is calculated.
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{21. \code{spectral}: Options for the spectral (turning bands) method}
\describe{
  \item{\code{ergodic}}{
    In case of an additive model and \code{ergodic=FALSE},
    the additive component are chosen proportional to their
    variance. In total \code{lines} are simulated. If
    \code{ergodic=TRUE}, the components are simulated
    separately and then added.
    
    Default: \code{FALSE}.
  }

 \item{\code{prop_factor}}{see \command{\link{RPspectral}}
    
  }

    
  \item{\code{sigma}}{see \command{\link{RPspectral}}
  }

  \item{\code{sp_grid}}{   see \command{\link{RPspectral}}
  }
 
  \item{\code{sp_lines}}{
    see \command{\link{RPspectral}}
  }
 
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{22. \code{tbm}: Options for the turning bands method}
\describe{
  \item{\code{center}}{Scalar or vector.
    If not \code{NA}, the \code{center} is used as the center of
    the turning bands for \code{TBM2} and \code{TBM3}.
    Otherwise the center is determined
    automatically such that the line length is minimal.
    See also \code{points} and the examples below.
    
    Default: \code{NA} .
  }

  \item{\code{fulldim}}{
    positiv integer. The dimension of the space into which the
    simulated field is embedded. So, the value \code{fulldim}
    must be at least the dimension of the field.

    Default: 3.
  }
  
  \item{\code{grid}}{Logical. 
    The angle of the lines is random if
    \code{grid=FALSE}, 
    and \eqn{k\pi/}{k*pi/}\code{lines}
    for \eqn{k}{k} in \code{1:lines},
    otherwise.
    
    This option is used by both \command{\link{RPspectral}}
    and \command{\link{RPtbm}}, the latter only when the dimension is 2.
    
    Default: \code{TRUE} .
  }
  
  \item{\code{layers}}{
    Logical or integer. If \code{TRUE} then the turning layers are used whenever
    a time component is given.
    If \code{NA} the turning layers are used only when the
    traditional TBM is not applicable.
    If \code{FALSE} then turning layers may never be used.
    
    Default: \code{TRUE} .
  }

  \item{\code{lines}}{
    Number of lines used.
    
    Default: \code{60} .
  }
  
  \item{\code{linesimustep}}{
    If \code{linesimustep} is positive the grid on the line has lag
    \code{linesimustep}. 
    See also \code{linesimufactor}.
    
    Default: \code{0.0} .
  }

 \item{\code{linesimufactor}}{ \code{linesimufactor} or
    \code{linesimustep} must be non-negative; if
    \code{linesimustep}
    is positive then \code{linesimufactor} is ignored.
    If both
    arguments are naught then \code{points} is used (and must be
    positive).
    The grid on the line is \code{linesimufactor}-times
    finer than the smallest distance. 
    See also \code{linesimustep}.
    
    Default: \code{2.0} .
  }

  \item{\code{points}}{integer. If greater than 0,
    \code{points} gives the number of points simulated on the TBM
    line, hence 
    must be greater than the minimal number of points given by
    the size of the simulated field and the two paramters
    \code{TBMx.linesimufactor} and \code{TBMx.linesimustep}.
    If \code{points} is not positive the number of points is
    determined automatically.
    The use of \code{center} and \code{points} is highlighted
    in an example below.
    
    Default: \code{0}.
  }

  \item{\code{reduceddim}}{
    if positiv integer, then the value itself. If negativ, then
    the value is substracted from fulldim.

    Default: \code{-2}.
  }
  
}
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{23. \code{internal}: Internal options mostly for warnings and
  messages}

All these options should not be changed by the user unless
he/she really known what he/she is doing.

Most of the options below change their value in a session
without the user's notice.


\describe{

  \item{\code{do_tests}}{
    Internal variable. Do not use it.
    Default: \code{FALSE}.
  }

  \item{\code{examples_reduced}}{non-negative integer.
    If positve, then the design of any simulation in \pkg{RandomFields}
    is internally reduced in size (roughly downto the given value in each
    direction). Warnings report this behaviour.
    This option is necessary to run the examples of \pkg{RandomFields}
    under the time constraint of CRAN.
    }
 
  \item{\code{stored.init}}{internally used logical argument.
    This option is closely related to 
    \code{storing} which controls whether intermediate calculations
    should be stored to have faster repeated simulations.

    This user option is internally overwritten if the user calls several
    simulations at once. This current value is stored in \code{stored.init}.
  
    Default: \code{FALSE}.
  }
  
  \item{\code{warn_ambiguous}}{internally used logical argument.
    Usually, the argument \code{grid} in \command{\link{RFsimulate}},
    for instance, can or should be given. If not given,
    the system takes a default definition.
    Additionally a message is displayed in this case if
    \code{ambiguous=TRUE}.
    
    Default: \code{FALSE}.
  }
  \item{\code{warn_aspect_ratio}}{internally used logical argument.
    if \code{TRUE} then a warning is given not a standard graphical
    device is used and the package plots try to keep a certain aspect
    ratio. 
    
    Default: \code{TRUE}
  }
      
  \item{\code{warn_colour_palette}}{internally used logical argument.
    If none of the packages \pkg{RColorBrewer} and \pkg{colorspace}
    are available and graphics are displayed, a message is displayed.

    Default: \code{TRUE}.
  }
  \item{\code{warn_constant}}{
    The definition of \command{\link{RMconstant}} has changed.
    A warning is displayed if the command is used. \code{warn_constant}
    will become obsolete in future versions.
    
    Default: \code{TRUE}.
  }
  
  \item{\code{warn_coordinates}}{internally used logical argument.
    If \code{TRUE} then a transformation from earth coordinates to
    cartesian coordinates is reported.
    
    Default: \code{TRUE}.
  }
 
  \item{\code{warn_missing_zenit}}{ % to do
    Only for Earth systems: a missing zenit is frequently a cause
    for errors that are difficult to understand. Therefore, in such
    cases an additional warning message is displayed.
    
    Default: \code{TRUE}
  }
  
  \item{\code{warn_newAniso}}{ obsolete.\cr
    internally used logical argument.
    If \code{newAniso=TRUE} and the argument \code{Aniso} is used in the model
    definition, then a message is displayed that the matrix \code{Aniso}
    is multiplied from the right by \eqn{x}, where up to Version 2.0
    the argument \code{aniso} was available which was multiplied from
    the left by \eqn{x}.
     
    Default: \code{TRUE}.
  }
  
%  \item{\code{warn_new_definitions}}{internally used logical argument.
%    If \code{warn_new_defintions=TRUE} then a warning is returned when 
%    models are used whose definition has changed recently.
%
%Default: \code{TRUE}.
%  }

  \item{\code{warn_newstyle}}{internally used logical argument.
    If \code{TRUE} a message is displayed the by the argument
    \code{spConform=FALSE} oldstyle return values are obtained instead
    of S4 objects.
    
    Default: \code{TRUE}.
  }
  
  \item{\code{warn_normal_mode}}{internally used logical argument.
    if \code{TRUE} then the function \command{\link{RFfit}}
    displays the message that other values for the option
    \code{modus_operandi} are available.
    
     Default: \code{TRUE}.
   }
   
  \item{\code{warn_oldstyle}}{internally used logical argument.
    If \code{TRUE} a warning is given if an obsolete function
    from Version 2 is used.

    Default: \code{TRUE}.
  }

  \item{\code{warn_on_grid}}{internally used logical argument.
    If a (one-dimensional) grid is given, but the argument
    \code{grid=FALSE}, e.g. in \code{RFsimulate}, this contraction is
    reported if \code{warn_on_grid=TRUE}
    
    
    Default: \code{TRUE}.
  }
  \item{\code{warn_scale}}{internally used logical argument.
    If \code{warn_scale=TRUE} then a scale less than 10 [km] is reported
    if earth coordinates are transformed to cartesian coordinates.
    
    Default: \code{TRUE}.
  }
  
  \item{\code{warn_var}}{
    In some cases, \pkg{RandomFields} cannot detect whether the
    variance is non-negative. If \code{TRUE} then a warning is displayed
    in such a case.
    Default: \code{TRUE}.
  }
 
}




%\bold{General comments}\cr
%Most of the above arguments determine the basic settings of a
%simulation, e.g. \code{matrix_method} (which chooses the method to
%calculate a square
%root of a positive definite matrix). The values of
%such arguments are read by
%\command{\link{RFsimulate}} and stored in an internal register
%and used in preliminary, deterministic calculations (e.g. calculation
%of the square of a covariance matrix)
%Some few other arguments like \code{lines} (which determines the number of
%i.i.d. processes to be simulated on the line)
%are only relevant when generating
%random numbers.
%These arguments are read only when actually simulations are performed, and
%are marked by \dQuote{[do]}.

%The difference of these two classes of arguments might become
%important only if\code{general.storing = TRUE} and
%\command{\link{RFsimulate}} is called with the only argument
%\code{n}.

}

 
 
\value{
 \code{NULL} if any argument is given, and the full list of
 arguments, otherwise.

 if \code{no.readonly=FALSE} then additionally,
 a list called \code{readonly} is included containing
 \cr
 * \code{covmaxchar}: the maximum length of a model name
 * \code{covnr}: number of currently implemented
 variogram/covariance models (-1 means that none of the functions like
 \command{\link{RFsimulate}}, \command{\link{RFfit}}
 , etc., have been called yet.)
 \cr
 * \code{distrmaxchar}: max. name length for a distribution
 \cr
 * \code{distrnr}: number of currently implemented
 distributions
 \cr
 * \code{maxdim}: maximum number of dimensions for a
 random field
 \cr
 * \code{maxmodels}: maximum number of elementary models in
 definition of a complex covariance model
 \cr
 * \code{methodmaxchar}: max. name length for methods
 \cr
 * \code{methodnr}: number of currently implemented simulation methods
}

\references{
  \itemize{
    \item General
    \itemize{
      \item
      Schlather, M. (1999) \emph{An introduction to positive definite
	functions and to unconditional simulation of random fields.}
      Technical report ST 99-10, Dept. of Maths and Statistics,
      Lancaster University.

      \item Schlather, M. (2011) Construction of covariance functions and
      unconditional simulation of random fields. In Porcu, E., Montero, J.M.
      and Schlather, M., \emph{Space-Time Processes and Challenges Related
	to Environmental Problems.} New York: Springer.
      % \item Schlather, M. (2002) Models for stationary max-stable
    }

    \item rectangular distribution; \code{eps_zhou}
    \itemize{
      \item
      Oesting, M., Schlather, M. and Zhou, C. (2013) On the Normalized
      Spectral Representation of Max-Stable Processes on a compact set. 
      \emph{arXiv},  \bold{1310.1813}
    }

    \item \code{shape_power}
    \itemize{
      \item
      Ballani, F. and Schlather, M. (2015) In preparation.
    }
  }
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\seealso{\command{\link{RFsimulate}},
  \link{RFoptionsAdvanced},
 \code{\link[=RandomFields-package]{RandomFields}},
 and \command{\link{RFgetMethodNames}}.}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
RFoptions()


############################################################
##                                                        ## 
## use of exactness                                       ##
##                                                        ##
############################################################
x <- seq(0, 1, 1/30)
model <- RMgauss()

for (exactness in c(NA, FALSE, TRUE)) { 
  readline(paste("\n\nexactness: `", exactness, "'; press return"))
  z <- RFsimulate(model, x, x, exactness=exactness,
                  stationary_only=NA, storing=TRUE)
  print(RFgetModelInfo(which="internal")$internal$name)
}

\dontshow{FinalizeExample()}

}

\keyword{spatial}

