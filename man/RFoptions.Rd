\name{RFoptions}
\alias{RFoptions}
\title{Setting control parameters}
\description{
 \command{\link{RFoptions}} sets and returns control parameters for the analysis
 and the simulation of random fields
}
\usage{
 RFoptions(..., no.readonly = TRUE)
}
\arguments{
 \item{...}{arguments in \code{tag = value} form, or a list of tagged
 values.}
 \item{no.readonly}{If \command{\link{RFoptions}} is called without
 parameter then all parameters are returned in a list. If
 \code{no.readonly=TRUE} then only rewritable parameters are returned.
 }
}
\details{
 The subsections below comment on options for\cr
 \bold{1. General options}\cr
 \bold{2. Gaussian random field simulation (\command{\link{RFsimulate}})}\cr
 \bold{3. Kriging (\command{\link{RFinterpolate}})}\cr
 \bold{4. Circulant embedding methods \command{\link{RPcirculant}},
 \command{\link{RPcutoff}}, \command{\link{RPintrinsic}}}\cr
 \bold{5. Matrix decomposition \command{\link{RPdirect}}}\cr
 \bold{6. Nugget effect \command{\link{RPnugget}}}\cr
 \bold{7. Sequential method \command{\link{RPsequential}}}\cr
 \bold{8. Spectral (turning bands) method \command{\link{RPspectral}}}\cr
 \bold{9. Turning band \command{\link{RPtbm}}}\cr
 \bold{10. Random coins (shot noise) \command{\link{RPcoins}},
 \command{\link{RPaverage}}}\cr
 \bold{11. Hyperplane tessellations \command{\link{RPhyperplane}}}\cr
 \bold{12. Max-stable random fields (\command{\link{RFsimulate}})}\cr
 \bold{13. \command{\link{BRmethods}}}\cr
 \bold{14. \command{\link{RRrectangular}}}\cr
 \bold{15. \command{\link{RFfit}}, \command{\link{RFratiotest}}}\cr
 \bold{16. \command{\link{RFempiricalvariogram}}}\cr
 \bold{17. \command{\link{RFgui}}}\cr
 \bold{18. Graphics}\cr
 \bold{19. Warnings and messages}\cr
 \bold{General comments}\cr\cr\cr
 
 
 \bold{1. General options}\cr
 \describe{
 \item{\code{modus_operandi}}{character. One of the values
   \code{"careless"}, \code{"sloppy"}, \code{"easygoing"},
   \code{"normal"}, \code{"precise"}, \code{"pedantic"},
   \code{"neurotic"} .
   \bold{This parameter is in an experimental stage and its definition
     and effects will change very likely in near future.}%to do
   This parameter sets a lot of parameter at once related to estimation
   and simulation. \code{"careless"} prefers rather fast algorithms,
   but the results
   might be very rough approximations. By way of contrast,
   \code{"neurotic"} will try very
   hard to return exact result at the cost of hugh computing times.

   Default: \code{"normal"}
 }
 \item{\code{printlevel}}{If \code{printlevel}\eqn{\le0}{<=0}
 there is not any output on the screen. The
 higher the number the more tracing information is given. 
 
 0 : no message\cr
 1 : important (error) messages and warnings\cr
 2 : less important messages\cr
 3 : recursive call tracing\cr
 4 : details of recursive call tracing\cr
 5 : function flow information of large functions\cr
 6 : errors that are internally treated\cr
 7 : details on intermediate calculations\cr
 8 : further details on intermediate calculations\cr
 
 Default: 1 [also do].\cr
 }
 
 \item{\code{storing}}{Logical.
 If \code{FALSE} then the intermediate results are
 destroyed after the simulation of the random field(s)
 or if an error had occured.
 If \code{storing=TRUE}, then
 additional simulations can be performed by calling
 \command{\link{RFsimulate}} with at most the parameter \code{n}.
 This call can then be much faster, but the a rather large
 amount of memory could be kept.
 
 When \code{storing} turned from \code{TRUE} to \code{FALSE} by
 global call then all registers are deleted.
 Advanced:
 With \code{\link{RFoptions}(storing=list(FALSE, register,
	model_register))}
 single registers can be deleted.
 
 Default: \code{FALSE} [do]. 
 }
 
 \item{\code{skipchecks}}{logical.
 If \code{TRUE}, several checks whether the given parameter values
 and the dimension are within the allowed range is skipped.
 Do not change the value of this variable except you really
 know what you do.
 
 Default: \code{FALSE} [also do]. 
 }
 
 \item{\code{every}}{integer.
 if greater than zero, then every \code{every}th iteration is
 printed if simulated by TBM or random coin method. The value zero
 means that nothing is printed.
 
 Default: \code{0} [do].
 }
 
 \item{\code{register}}{0:9; place where intermediate calculation
 for random field simulation are stored; 
 the number refers to 10 internal registers 0..9.
 
 Changing the register number only makes sense, when
 two different random fields, say, are to simulated in
 alternatingly, several times in a row. Then the
 simlulation speed can be increased if several registers
 are used and \code{storing=TRUE}.
 
 Default: \code{0} [also do].
 }
 
 \item{\code{interpolregister}}{Register number for kriging

 Default: \code{9} [also do].
 }
 
 \item{\code{condregister}}{Register number for conditional simulation

 Default: \code{8} [also do].
 }
 
 \item{\code{errregister}}{Additional register for conditional
 simulation if error term is envolved.

 Default: \code{7} [also do].
 }
 
 \item{\code{guiregister}}{Register number for RFgui.

 Default: \code{0} [also do].
 }
 
 \item{\code{gridtolerance}}{
 used in \command{\link{RFsimulate}} to see if the coordinates build a
 grid for x, y, z, T-values. This parameter is also used
 in case of conditional
 simulation where the data locations might ly on
 a grid.
 
 Default: \code{1e-6}
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

 Default: \code{'*'} [do]. 
 }
 
 \item{\code{practicalrange}}{logical or integer.
 If not \code{FALSE} the range of primitive
 covariance functions is
 adjusted so that cov(1) is zero for models with finite range.
 (Operators are too complex to be adjusted; for anisotropic
 covariance the practical range is not well defined.)
 
 The value of cov(1) is about 0.05 (for \code{scale=1})
 for models without range. See \command{\link{RMmodel}} or type
 \code{\link{RFgetModelNames}(type="positive definite",
   domain="single variable", isotropy="isotropic", operator=FALSE, vdim=1)}
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
 
 \item{\code{spConform}}{logical.
 If \code{TRUE} then \command{\link{RFsimulate}} and many other functions
 return an \code{sp}-object. Otherwise matrices or lists are
 returned as defined in RandomFields 2.0, see the manuals for the
 specific functions.

 Default: \code{TRUE} [do].
 }
 
 \item{\code{cPrintlevel}}{\code{printlevel} above works
 on the R level whereas \code{cPrintlevel} works on the C level.
 \code{cPrintlevel} is automatically set when \code{printlevel}
 is changed, namely as \code{min(printlevel, 5)}. 
 
 0 : no message\cr
 1 : important (error) messages and warnings\cr
 2 : less important messages\cr
 3 : recursive call tracing\cr
 4 : details of recursive call tracing\cr
 5 : function flow information of large functions\cr
 6 : errors that are internally treated\cr
 7 : details on building up the covariance structure\cr
 8 : details on taking the square root of the covariance matrix\cr
 9 : details on intermediate calculations\cr
 10 : further details on intermediate calculations\cr
 
 Default: 1 [also  do].\cr
 }

 \item{\code{exactness}}{logical or NA.
 \itemize{
	\item \code{TRUE}: add.MPP, hyperplanes and
	all turning bands methods are excluded.
	If the circulant embedding method is considered as badly
	behaved, then the matrix decomposition methods are preferred.
	\item \code{FALSE}: if the circulant embedding method is
	considered as badly behaved or the number of points to be
	simulated is large, the turning bands methods are
	rather preferred.
	\item \code{NA}: approximative non-exact methods are excluded,
	i.e. TBM2 if the Abel transform of the covariance function
	cannot be given explicitely.
 }
 
 Default: \code{NA} .
 }

 \item{\code{matrix_inversion}}{
   vector of integers gives the methods (and their sequence)
   to inverse a matrix: 

     0 : Choleskey decomposition \cr
     1 : QR decomposition\cr
     2 : SVD\cr
   
   This option is currently used only in  \command{\link{RFinterpolate}}.% to
   % do  
   
   Default: \code{c(0, 2)}.
 }
 
 \item{\code{matrix_tolerance}}{
   In case of SVD, it is checked whether the eigenvalue is less than
   \code{matrix_tolerance}. Then the eigenvalue is not inverted, but
   the inverse is set to zero.
 }
 
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

 \item{\code{na_rm_lines}}{
   logical. If \code{TRUE} then a line of the data that contains a
   \code{NA} value is deleted. Otherwise it is tried to deal with the
   \code{NA} value at higher costs of computing time.

   Default: \code{FALSE}.
 }

\item{\code{vdim_close_together}}{logical. Used in functions only that
  create covariance matrices. If the model is multivariate, then two
  ways of ordering the matrix exists. To consider first all variables at
  a certain location (\code{vdim_close_together=TRUE}) or to consider first
  all locations keeping the variable fixed
  (\code{vdim_close_together=FALSE}). 

   Default: \code{FALSE}.
 }

 \item{\code{expected_number_simu}}{positive integer which is usally set
   internally as the value of the parameter \code{n} in
   \command{\link{RFsimulate}}. The parameter \code{expected_number_simu}
   should be set only by an advanced users and only if 
   \command{\link{RFsimulate}} will be called with parameter \code{n} allone.
 }

 \item{\code{xyz_notation}}{logical, NA

   \code{NA} : automatic choice (if possible)
   
   \code{false} : notation (x, y) should not be understood as
   as kernel definition, not as xyz notation
   
   \code{true}: xyz notation used

   % \code{2}:this value is only used by calls of RFcov and should not
   % be used by a user
 }

 \item{\code{coordinate_system}}{
   One of the values \code{"auto"}, \code{"cartesian"},
%   \code{"polar"}, \code{"cylindric"}, \code{"spherical"},
   \code{"earth"}

   Default: \code{"auto"}
 }

 \item{\code{coord_units}}{any string.
     
    If \code{coordinate_system = "earth"} and longitude and latitude
    are transformed to 3d cartesian coordinates, \code{coord_units}
    determines 
    whether the radius is given in kilometers (\code{"km"}) or miles
    (\code{"miles"}).
      If empty, then \code{"km"} is chosen.
 
    Default: \code{""}
  }
  

 \item{\code{new_coord_units}}{internal and should not be set
  internal.

  Default: \code{""}
}
  
   
  \item{\code{variab_units}}{vector of characters.
    The default units of the variables
    
    Default: \code{""}
  }

  \item{seed}{integer. If not \code{NULL} and not \code{NA}
    the \code{\link[base]{set.seed}(seed)} is set
    before simulations are performed, e.g. by
    \command{\link{RFsimulate}} or \command{\link{RFdistr}}.
    Otherwise  \command{\link{set.seed}} is \bold{not} called.

    If the parameter is set locally, i.e. within a function,
    it has the usual local effect. If it is set globally, i.e. by
    \command{RFoptions} the \code{seed} is fixed
    for all subsequent calls.

    Note that \command{\link{RFratiotest}} has its own parameter
    \code{seed} with a slightly different meaning.
  }
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{2. Gaussian random field simulation}

 \describe{
 
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

 \item{\code{approx_zero}}{
   Value below which a correlation is considered to be essentially zero.
   This parameter is used to determine the practical range of covariance
   function with non-compact support.
 
   Default: \code{0.05}
 }

 \item{\code{direct_bestvar}}{integer.
 When searching for an appropriate simuation method
 the matrix decomposition method (\code{method="direct"})
 is preferred if the number of variables is less than or equal to 
 \code{direct_bestvariables}.

 Default is \code{800}.
 }

 \item{\code{loggauss}}{
   See \command{\link{RPgauss}}
 }
}

 
 \bold{3. Options for Kriging}\cr
 \describe{
 \item{\code{method}}{
   \code{'S'} : simple kriging\cr
   \code{'O'} : ordinary kriging\cr
   \code{'M'} : kriging the mean\cr
   \code{'U'} : universal kriging\cr
   \code{'I'} : intrinsic kriging\cr
   \code{'A'} : automatic choice determined on the given
   model. In some situation an automatic choice is not possible and
   the \code{method} must be given explicitely.
 
   Default: \code{'A'}.
 }
 

 \item{\code{return.variance}}{logical.
 If \code{FALSE} the kriged field is
 returned. If \code{TRUE} a list of two elements, \code{estim} and
 \code{var}, i.e. the kriged field and the kriging variances,
 is returned.

 Default: \code{FALSE}. 
 }
 
 \item{\code{locmaxn}}{ % to do

   Default: \code{8000}. 
 }
 
 \item{\code{locsplitn}}{ % to do

   Default: \code{c(5000, 200, 1000)}. 
 }

 \item{\code{locsplitfactor}}{ % to do

   Default: \code{2}. 
 }
 
 \item{\code{fillall}}{ % to do

  Default: \code{TRUE}.
 }
 

 
 \item{\code{cholesky_R}}{
 if \code{TRUE} cholesky decomposition is used instead of LU.

 Default: \code{FALSE}. 
 }
 }


 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{4. Options for simulating with the standard circulant embedding
 method, including cutoff, intrinsic}\cr
 \describe{
 \item{\code{force}}{See \command{\link{RPcirculant}}}
 \item{\code{mmin}}{See \command{\link{RPcirculant}}}
 \item{\code{strategy}}{See \command{\link{RPcirculant}}}
 \item{\code{maxmem}}{See \command{\link{RPcirculant}}}
 \item{\code{tolIm}}{See \command{\link{RPcirculant}}}
 \item{\code{tolRe}}{See \command{\link{RPcirculant}}}
 \item{\code{trials}}{See \command{\link{RPcirculant}}}
 \item{\code{useprimes}}{See \command{\link{RPcirculant}}}
 \item{\code{dependent}}{See \command{\link{RPcirculant}}}
 \item{\code{approx_step}}{See \command{\link{RPcirculant}}}
 \item{\code{approx_maxgrid}}{See \command{\link{RPcirculant}}}
 }
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{5. Options for simulating by simple matrix decomposition}\cr
 \describe{
 \item{\code{root_method}}{See \command{\link{RPdirect}}}
 \item{\code{svdtolerance}}{See \command{\link{RPdirect}}}
 \item{\code{max_variab}}{See \command{\link{RPdirect}}}
 }
 


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{6. Options for simulating nugget effects}\cr
 Simulating a nugget effect is per se trivial.
 However, it gets complicated
 and best methods (including \code{direct} and \code{circulant
 embedding}!) fail if zonal anisotropies are considered,
 where sets of points have to be identified that belong to the
 same subspace of eigenvalue 0 of the anisotropy matrix.
 \describe{
 \item{\code{tol}}{
   see \command{\link{RPnugget}}
 } 
 }
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{7. Options for using the sequential method}\cr
 \describe{
 \item{\code{max_variables}}{See \command{\link{RPsequential}}}
 \item{\code{back_steps}}{See \command{\link{RPsequential}}}
 \item{\code{initial}}{See \command{\link{RPsequential}}}
 }


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{8. Options for simulating with spectral (turning bands) method}\cr
 \describe{
 \item{\code{sp_lines}}{
   see \command{\link{RPspectral}}
 }
 \item{\code{sp_grid}}{   see \command{\link{RPspectral}}
 }
 \item{\code{ergodic}}{
 In case of an additive model and \code{ergodic=FALSE},
 the additive component are chosen proportional to their
 variance. In total \code{lines} are simulated. If
 \code{ergodic=TRUE}, the components are simulated
 separately and then added.
 
 Default: \code{FALSE} .
 }

 
 \item{\code{prop_factor}}{see \command{\link{RPspectral}}

 }
 
 \item{\code{sigma}}{see \command{\link{RPspectral}}
 }

 }

 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{9. Options for simulating with a turning bands method}\cr
 
 \describe{
 \item{\code{tbmdim}}{
 if positiv integer, then the value itself. If negativ, then
 the value is substracted from fulldim.
 }
 
 \item{\code{fulldim}}{
 positiv integer.
 }
 
 \item{\code{center}}{Scalar or vector.
 If not \code{NA}, the \code{center} is used as the center of
 the turning bands for \code{TBM2} and \code{TBM3}.
 Otherwise the center is determined
 automatically such that the line length is minimal.
 See also \code{points} and the examples below.
 
 Default: \code{NA} .
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
 
 Default: \code{0} .
 }
 
 \item{\code{lines}}{
 Number of lines used.
 
 Default: \code{60} .
 }
 
 \item{\code{linesimufactor}}{ \code{linesimufactor} or
 \code{linesimustep} must be non-negative; if
 \code{linesimustep}
 is positive then \code{linesimufactor} is ignored.
 If both
 parameters are naught then \code{points} is used (and must be
 positive).
 The grid on the line is \code{linesimufactor}-times
 finer than the smallest distance. 
 See also \code{linesimustep}.
 
 Default: \code{2.0} .
 }
 \item{\code{linesimustep}}{
 If \code{linesimustep} is positive the grid on the line has lag
 \code{linesimustep}. 
 See also \code{linesimufactor}.
 
 Default: \code{0.0} .
 }
 \item{\code{layers}}{
 Logical or integer. If \code{TRUE} then the turning layers are used whenever
 a time component is given.
 If \code{NA} the turning layers are used only when the
 traditional TBM is not applicable.
 If \code{FALSE} then turning layers may never be used.
 
 Default: \code{TRUE} .
 }
 \item{\code{grid}}{Logical.
 The angle of the lines is random if
 \code{grid=FALSE}, 
 and \eqn{k\pi/}{k*pi/}\code{lines}
 for \eqn{k}{k} in \code{1:lines},
 otherwise. This parameter is only considered
 if the spectral measure, not the density is used.
 
 Default: \code{TRUE} .
 }
 }
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{10. Options for simulating by random coins (shot noise)}\cr
 \describe{
 \item{\code{n_estim_E}}{integer. Number of draws from the
 distribution of the scale to estimate the mean of the distribution.
 This is used only if the mean of the scale distribution
 is not explicitely given.
 Default: \code{50000} .
 }
 \item{\code{intensity}}{real.
 Number of superposed
 realisations (to approximate the normal distribution; total number
 for all (additive) components with same anisotropy);
 if \code{intensity<=0.0} then only a single
 value is simulated (for checking).
 
 Default: \code{100} .
 }
 
 \item{\code{about_zero}}{Functions that
 do not have 
 compact support are set to zero outside the ball outside for which the
 function has absolute values less than \code{about_zero}.
 
 Default: \code{0.001} .
 }
 
 }
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{11. Options for simulating hyperplane tessellations}\cr
 \describe{
 \item{\code{superpos}}{integer.
 number of superposed hyperplane tessellations.
 
 Default: \code{300} .
 }
 \item{\code{maxlines}}{integer.
 Maximum number of allowed lines.
 
 Default: \code{1000} .
 }
 \item{\code{mar_distr}}{integer.
   This parameter should not be changed yet.
 
   It  codes the marginal distribution used in the
 simulation:
 
 0 : uniform distribution\cr
 1 : Frechet distribution with form parameter \code{mar.param}\cr
 2 : Bernoulli distribution (Binomial with \eqn{n=1}) with
     parameter \code{mar.param}

 
 Default: \code{0} .
 }
 \item{\code{mar_param}}{Parameter used for the marginal
 distribution. The parameter should not be changed yet.
 
 Default: \code{NA} .
 }
 }

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{12. Options specific to simulating max-stable random fields}
 \describe{
 \item{\code{maxGauss}}{
 The simulation of the max-stable process based on random fields uses
 a stopping rule that necessarily needs a finite upper endpoint
 of the marginal distribution of the random field.
 In the case of Brown-Resnick processes and extremal Gaussian fields
 (\command{\link{RPbrownresnick}} and
 \command{\link{RPschlather}}, respectively), the upper endpoint is
 approximated by \code{standardmax}. 
 
 Default: \code{3.0} .
 }
 \item{\code{maxpoints}}{
  positive integer; the maximal number of Poisson points to be simulated
  for one realization of the max-stable random field. This option will
  not be considered for most of the users.
  
  Default: \code{2e9} .
 }

 \item{\code{xi}}{
   Extreme value index.
   
 Default: \code{1.0} .
 }

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
   value in \eqn{[0,1]}. This parameter is considered only
   if \code{flat=-1} and the simulation is performed on a grid.
   Then, the ratio between the highest and the lowest value is
   calculated within the convex hull of the grid. If the
   value is less than \code{density_ratio} then the grid points
   are considered separately. Else the density is considered to be
   constant in the convex hull of the grid.
   
   
   Default: .
 }
 
 \item{\code{flat}}{\code{-1, FALSE, TRUE}.
   The parameter is considered onlyu if the simulation is performed on a
   grid. 
   If \code{flat} is logical, then the density
   is considered to \code{flat} in the convex hull of the grid.
   If \code{flat=-1} the choice is done automatically.
     
 Default: -1 .
 }
 


 
 }
 
 
 
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  \bold{13. Specific method options for Brown-Resnick Fields}
  \describe{
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
    \item{\code{lowerbound_optimarea}}{
      the lower bound \eqn{x_0} of the domain \eqn{[x_0, \infty] \times E}
      {[x_0, Inf] x E} of the Poisson point process which is used for the
      simulation of the shape functions in the M3 representation of BR
      processes (cf. Oesting et al., 2012); only used for simulation of
      BR processes via \code{RPbrmixed}.
      
      Default: \code{-2.0} .
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
    \item{\code{optim_mixed}}{\code{0, 1, 2}; only used for simulation of BR
      processes via \code{RPbrmixed}.\cr
      If \code{optim_mixed=0}, the arguments \code{lambda} and
      \code{areamat} of \command{\link{RPbrshifted}} are used for the
      simulation.\cr
      If \code{optim_mixed=1}, \code{lambda} is estimated for
      \code{areamat=1}.\cr
      If If \code{optim_mixed=2}, \code{areamat} is optimized and
      \code{lambda} is estimated.
      
      Default: \code{0} .
    }
    \item{\code{optim_mixed_tol}}{
      value in \eqn{[0,1]}; only used for simulation of BR processes via
      \code{RPbrmixed} with \code{optim_mixed=2}. In this case,
      \code{areamat} is optimized under the constraint that the
      probability of drawing the shape function incorrectly is bounded by
      \code{optim_mixed_tol} (cf. Oesting et al., 2012).
      
      Default: \code{0.01} .
    }
    \item{\code{optim_mixed_maxpoints}}{
      positive integer; only used for simulation of BR processes via
      \code{RPbrmixed} with \code{optim_mixed>0}. Maximal number of Poisson
      points used for the optimization of \code{areamat} and the estimation
      of \code{lambda}.
      
      Default: \code{1000000} .
    }
    
    \item{\code{variobound}}{
      positive; the shape functions in the mixed moving maxima
      representation are cut off where the variogram belonging
      to \code{phi} exceeds \code{variobound}.
      
      Default: \code{3.0} .
    }
    
 }
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{14. Technical options for \command{\link{RRrectangular}}}
 \describe{
  \item{\code{safety}}{ % to do
 
 Default: \code{0.08} .
 }

 \item{\code{minsteplen}}{
   Technical parameters to run an algorithm to simulate from the
   \command{\link{RRrectangular}} distribution.
   
 Default: \code{0} (i.e. not used as a criterion.)
 }
 
 \item{\code{innermin}}{
   Technical parameters to run an algorithm to simulate from the
   \command{\link{RRrectangular}} distribution.
 
 Default: \code{1e} .
 }


 \item{\code{outermax}}{
   Technical parameters to run an algorithm to simulate from the
   \command{\link{RRrectangular}} distribution.
 
 Default: \code{1e5} .
 }
 
 \item{\code{maxsteps}}{
   Technical parameters to run an algorithm to simulate from the
   \command{\link{RRrectangular}} distribution.
 
 Default: \code{1000} .
 }
 
 \item{\code{parts}}{
   Technical parameters to run an algorithm to simulate from the
   \command{\link{RRrectangular}} distribution.
 
 Default: \code{8} .
 }
 
 \item{\code{maxit}}{
   Technical parameters to run an algorithm to simulate from the
   \command{\link{RRrectangular}} distribution.
 
 Default: \code{20} .
 }
 
 \item{\code{mcmc_n}}{
   In case of the use of MCMC it leaves out \eqn{n-1}
   member of the Markov chain bevor the \eqn{n} member
   is returned.
 
   Default: \code{15} .
 }
 


 }
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{15. Options for \command{\link{RFfit}} and \command{\link{RFratiotest}}}
 The following comments on \command{\link{RFfit}} apply
 also to \command{\link{RFratiotest}}.

 \describe{
 \item{\code{bin_dist_factor}}{
 
 Default: \code{0.5}.
 }
 \item{\code{upperbound_scale_factor}}{
 The upper bound for the scale is determined
 as \code{upperbound_scale_factor}\eqn{ * (maximum distance
	between all pairs of points)}.
 
 Default: \code{3}.
 }
 \item{\code{lowerbound_scale_factor}}{
 The lower bound for the scale is determined as
 \eqn{(minimum distance
	between different pairs of points) / }
 \code{lowerbound_scale_factor}.
 
 Default: \code{3}.
 }
 \item{\code{lowerbound_scale_ls_factor}}{ For the LSQ target
 function a different lower bound 
 for the scale is used. It is determined as
 \eqn{(minimum distance between different pairs of points) / }
 \code{lowerbound_scale_ls_factor}.
 
 Default: \code{5}.
 }
 \item{\code{upperbound_var_factor}}{ The upper bound for the
 variance and the nugget is determined 
 as \code{upperbound_var_factor} * var(\code{data})
 
 Default: \code{10}.
 }
 \item{\code{lowerbound_var_factor}}{
 The lower bound for the nugget and the variance is determined
 as var(\code{data}) / \code{lowerbound_var_factor}.
 If a standard model definition is given and
 either the nugget or the variance is fixed,
 the parameter to be estimated
 must also be greater than \code{lowerbound_sill}.
 
 Default: \code{10000}.
 }
 \item{\code{lowerbound_sill}}{absolute lower bound for variance
 and nugget. See \code{lowerbound_var_factor}.
 
 Default: \code{1E-10}.
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
 \item{\code{approximate_functioncalls}}{
 In case the parameter vector is too close to the given
 bounds, the ML target function is evaluated on a grid
 to get a new initial value for the ML estimation. 
 The number of points of the grid is approximately
 \code{approximate_functioncalls}.
 
 Default: \code{50}
 }
 \item{\code{refine_onborder}}{
 }
 \item{\code{minmixedvar}}{
 lower bound for variance in a mixed model;
 so, the covariance model for mixed model part might
 be calibrated appropriately
 }
 \item{\code{maxmixedvar}}{
 upper bound for variance in a mixed model;
 so, the covariance model for mixed model part might
 be calibrated appropriately
 }
 \item{\code{solvesigma}}{Boolean -- experimental stage!
 If a mixed effect part is present where the variance
 has to be estimated, then this variance parameter is solved
 iteratively within the profile likelihood function, if
 \code{solvesigma=TRUE}.This makes sense
 if the number of independent variables is very small.
 If \code{solvesigma=FALSE} then the variance parameter is
 treated as any other parameter to be estimated.}
% \item{see also \code{allowdistanceZero}}{ }
% \item{\code{na_rm}}{boolean -- experimental stage. 
% Only the data may have missing values.
% If \code{na_rm=TRUE} then lines of (repeated) data are deleted if at
% least one missing value appears. If \code{na_rm=FALSE} then the
% repetitions are treated separately.
% }
 \item{\code{bc_lambda_lb}}{
 lower bound for the Box-Cox transformation
 
 Default: \code{-10}.
 }
 \item{\code{bc_lambda_ub}}{
 upper bound for the Box-Cox transformation
 
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
 
 Note that this parameter does not influence
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
 \item{\code{bins}}{vector of explicit boundaries for the bins or the
 number of bins for the empirical variogram (used in the
 LSQ target function, which is described at the beginning
 of the Details).
 Note that for anisotropic models, the value of \code{bins} might
 be enlarged.
 
 Default: \code{20}.
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
 \item{\code{sill}}{
 Additionally to estimating \code{nugget} and \code{variance}
 separately, they may also be estimated together under the
 condition that \code{nugget} + \code{variance} = \code{sill}.
 For the latter a finite value for \code{sill} has to be supplied,
 and \code{nugget} and \code{variance} are set to \code{NA}.
 
 \code{sill} is only used for the standard model. 
 
 Default: \code{NA}. 
 } 
 \item{\code{only_users}}{boolean.
 If true then only \code{users_guess} is used as a
 starting point for the fitting algorithms
 
 Default: \code{FALSE}. 
 }
 \item{\code{optim_var_estimation}}{This parameter takes the values
 \code{'never'}, \code{'respect bound'}, \code{'try'},
 \code{'yes'}, and should only be 
 set by the advanced user.

 The meaning of the values is as follows.
 \itemize{
	\item{\code{'never'}}{
	 A global variance is never tried to be eliminated
	}
	\item{\code{'respect bound'}}{
	 A global variance is eliminated if such a variance
	 is detected and the user did not indicate bounds
	 for the parameters.
	}
	\item{\code{'try'}}{
	 A global variance is eliminated if such a variance
	 is detected.
	}
	\item{\code{'yes'}}{
	 A global variance is tried to be eliminated altough
	 the algorithm dis not find an indication. Here, the full
	 responsibility is left to the user. (This option might
	 make sense if \code{transform} is given.)
	 This option is only overwritten when it does not make sense,
	 e.g. no variance is estimated.
	}
      }
      Default: \code{'respect bound'}.
 }
 \item{\code{shortnamelength}}{
 The names of the variables in the returned table are
 abbreviated by taking the first \code{shortnamelength}
 letters.

 Default: \code{4}. 
 }
 \item{\code{use_spam}}{
 Should the package \code{spam} (sparse matrices)
 be used for matrix calculations?
 
 Default: \code{NA}. 
 }
 \item{\code{split}}{logical.
 If \code{TRUE} then \command{\link{RFfit}} checks whether a
 space-time covariance model or a multivariate covariance model
 can be split into components, so that certain parameters
 can be estimated separately.
 
 Default: \code{TRUE}. 
 }
 \item{\code{scale_ratio}}{
 \command{\link{RFfit}} uses \code{parscale} and \code{fnscale}
 in the calls of \command{\link{optim}}. As these parameters should
 have the magnitude of the estimated values, \command{\link{RFfit}}
 checks this by calculating the absolute log ratios.
 If they are larger than \code{scale_ratio},
 \code{parscale} and \code{fnscale} are reset and the optimisation
 is redone. 
 
 Default: \code{0.1}. 
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
 of the parameters.
 
 If \code{critical>=3} then a kind of profile likelihood
 optimization is always done (which takes about \code{n_crit}
 times extra time) for all the parameters.
 
 If \code{critical<0} then none of the refined methods
 are performed.
 
 Default: \code{TRUE}.
 }

 \item{\code{n_crit}}{integer.
 The approximate profiles that are considered.
 
 Default: \code{10}.
 }
 
 \item{\code{max_neighbours}}{integer.
 Maximum number of locations (with depending values)
 that are allowed.
 
 Default: \code{5000}.
 }
 
 \item{\code{splitn_neighbours}}{integer.
 In case the maximum number of locations \code{maxn}
 is exceeded, then \command{\link{RFfit}} tries to split the data set
 into parts of size \code{split} or less, but never more than 
 \code{maxn}. %to do
 
 Default: \code{c(3000, 200, 1000)}.
 }
 
 \item{\code{splitfactor_neighbours}}{% to do 3x
Default: \code{2}.
  }
 
 \item{\code{smalldataset}}{
 Default: \code{2000}.
 }
 
 \item{\code{min_diag}}{
 Default: \code{1e-7}.
 }
 

 \item{\code{reoptimise}}{logical.
 If \code{TRUE && !only_users} then at a very last step,
 the optimisation is redone with currently best parameters
 and likelihood as scale parameter for \command{\link{optim}}.

 Default: \code{TRUE}. 
 }
 
 \item{\code{optimiser}}{% to do 2x
Default: \code{0}.
  }
 
 \item{\code{algorithm}}{
Default: \code{-1}.
  }
 
 }



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{16. Options for calculating the empirical variogram}\cr
 \describe{ % to do
 \item{\code{phi0}}{
 
 Default: \code{0}.
 }
 \item{\code{theta0}}{
 
 Default: \code{0}.
 }
 \item{\code{tol0}}{
 
 Default: \code{1e-13}.
 }
 \item{\code{pseudovariogram}}{
 
 Default: \code{FALSE}.
 }
 \item{\code{fft}}{
 
 Default: \code{TRUE}.
 }
 }
 

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{17. Options for calculating the \command{\link{RFgui}}}\cr
 \describe{%to do
 \item{\code{alwaysSimulate}}{
 
 Default: \code{TRUE}.
 }
 \item{\code{simu_method}}{
 
 Default: \code{"RPcirculant"}.
 }
 \item{\code{size}}{
 
 Default: \code{c(512, 64)}.
 }
 }

 \bold{18. Options for graphical output}\cr
 \describe{%to do
   \item{\code{always_close_screen}}{logical.
     If \code{FALSE} the current device is kept as it is.
     Otherwise the action depends on the value of
     \code{height}: if \code{height} is not positive
     then \command{\link[graphics]{close.screen}}
     is performed on the current device.
     Else the current device is closed.

     Default: \code{TRUE}.
 }
 \item{\code{grPrintlevel}}{
   integer values 0, 1, 2. The higher the more text is shown in the plot.
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

   
   Default: \code{0}.
 }

 \item{\code{increase_upto}}{
   See \code{height}.
   
   Default: \code{c(3,4)}.
 }

 
 
 }

 \bold{19. Options for warnings and messages}\cr
 \describe{
 \item{\code{oldstyle}}{internal
 }
 \item{\code{newstyle}}{internal
 }
 \item{\code{newAniso}}{internal
 }
 \item{\code{ambiguous}}{internal
 }
\item{\code{normal_mode}}{internal
 }
\item{\code{warn_mode }}{internal
 }

 }
 
 \bold{General comments}\cr

 
 
 Most of the global parameters determine the basic settings of a
 simulation, e.g. \code{root_method} (which chooses the method to
 calculate a square
 root of a positive definite matrix). The values of
 such parameters are read by
 \command{\link{RFsimulate}} and stored in an internal register
 and used in preliminary, deterministic calculations (e.g. calculation
 of the square of a covariance matrix)
% These parameters have the flag \dQuote{}.
 Some few other parameters like \code{lines} (which determines the number of
 i.i.d. processes to be simulated on the line)
 are only relevant when generating
 random numbers.
 These parameters are read only when actually simulations are performed, and
  are marked by \dQuote{[do]}.

 The differnce of these two classes of parameters might become
 important only if \code{general.storing = TRUE} and
 \command{\link{RFsimulate}} is called only if with the parameter
 \code{n}.

}

 
 
\value{
 \code{NULL} if any parameter is given, and the full list of
 parameters, otherwise.

 if \code{no.readonly=FALSE} then additionally, the last element
 of the list is itself a list containing
 \cr
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
 Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University. 
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}}

\seealso{\command{\link{RFsimulate}},
 \code{\link[=RandomFields-package]{RandomFields}},
 and \command{\link{RFgetMethodNames}}.}

\examples{
RFoptions(seed=0)
RFoptions(storing=TRUE)
str(RFoptions())


############################################################
## ##
## use of exactness ##
## ## 
############################################################
% library(RandomFields, lib="~/TMP")
x <- seq(0, 1, if (interactive()) 1/30 else 0.5)
model <- RMgauss()

for (exactness in c(NA, FALSE, TRUE)) { 
  readline(paste("\n\nexactness: `", exactness, "'; press return"))
  z <- RFsimulate(model, x, x, grid=TRUE, exactness=exactness,
                  stationary_only=NA, storing=TRUE)
  Print(RFgetModelInfo(which="keys")$key$name)
}


\dontshow{\dontrun{
 #############################################################
 ## The following gives a tiny example on the advantage of  ##
 ## local.dependent=TRUE (and dependent=TRUE) if in a       ##
 ## study most of the time is spent with simulating the     ##
 ## Gaussian random fields. Here, the covariance at a pair  ##
 ## of points is estimated for n independentent repetitions ##
 ## and 2*n locally dependent dependent repetitions .       ##
 ## To get the precision, the procedure is repeated m times.##
 #############################################################

# In the example below, local.dependent speeds up the simulation
# by about factor 16 at the price of an increased variance of
# factor 1.5

x <- seq(0, 1, len=10)
y <- seq(0, 1, len=10)
grid.size <- c(length(x), length(y))
meth <- RPcirculant
model <- RMexp(var=1.1, Aniso=matrix(nc=2, c(2,0.1,1.5,1)))
#(RMmodel(matrix(c(1, -1), ncol=2), model=model)) 

m <- if (interactive()) 100 else 2
n <- if (interactive()) 100 else 10

# using local.dependent=FALSE (which is the default)
c1 <- numeric(m)
time <- unix.time(
  for (i in 1:m) {
    cat("", i)
    z <- RFsimulate(meth(model), x, y, grid=TRUE, n=n, pch="", 
                    dependent=FALSE, spConform=FALSE, trials=5)
    c1[i] <- cov(z[1,length(y), ], z[length(x), 1, ])
}) # many times slower than with local.dependent=TRUE below
Print(time, mean(c1), sd(c1))

# using local.dependent=TRUE...
c2 <- numeric(m)
time <- unix.time(
  for (i in 1:m) {
    cat("", i)
    z <- RFsimulate(meth(model), x, y, grid=TRUE, n=2 * n, pch="", 
                    dependent=TRUE, spConform=FALSE, trials=5)
    c2[i] <- cov(z[1,length(y),], z[length(x), 1 , ])
})
Print(time, mean(c2), sd(c2))
# the sd is samller (using more locally dependent realisations)
## but it is (much) faster! Note that for n=n2 instead of n=2 * n, 
## the value of sd(c2) would be larger due to the local dependencies 
## in the realisations.
}}



\dontshow{\dontrun{
############################################################
## ##
## use of points and center ##
## ##
############################################################
## The following example shows that the same realisation 
## can be obtained on different grid geometries (or point 
## configurations, i.e. grid, non-grid) using TBM3 (or TBM2)



x1 <- seq(-150,150,1)
y1 <- seq(-15, 15, 1)
x2 <- seq(-50, 50, 1)
model <- RPtbm(RMexp(scale=10))

RFoptions(storing=TRUE)
mar <- c(2.2, 2.2, 0.1, 0.1)
points <- 700

###### simulation of a random field on long thin stripe
z1 <- RFsimulate(model, x1, y1, grid=TRUE, center=0, seed=0,
                 points=points, storing=TRUE, spConform=FALSE)
do.call(getOption("device"), list(height=1.55, width=12))
par(mar=mar)
image(x1, y1, z1, col=rainbow(100))
polygon(range(x2)[c(1,2,2,1)], range(y1)[c(1,1,2,2)],
        border="red", lwd=3)


###### definition of a random field on a square of shorter diagonal
z2 <- RFsimulate(model, x2, x2, grid=TRUE, register=1, seed=0,
                 center=0, points=points, spConform=FALSE)
do.call(getOption("device"), list(height=4.3, width=4.3))
par(mar=mar)
image(x2, x2, z2, zlim=range(z1), col=rainbow(100))
polygon(range(x2)[c(1,2,2,1)], range(y1)[c(1,1,2,2)],
        border="red", lwd=3)
tbm.points <- RFgetModelInfo(RFsimulate)$loc$totpts
Print(tbm.points, empty.lines=0) # number of points on the line
}}

\dontshow{RFoptions(seed=NA)}

}

\keyword{spatial}

