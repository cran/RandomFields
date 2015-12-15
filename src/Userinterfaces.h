


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#ifndef RFuser_H
#define RFuser_H 1


#define R_PRINTLEVEL 1
#define C_PRINTLEVEL 1
#define NAT_SCALE 0
#define MAX_CE_MEM 16777216
#define generalN 23
// IMPORTANT: all names of general must be at least 3 letters long !!!
extern const char *general[generalN];
#define GENERAL_MODUS 0 
#define GENERAL_STORING 2
#define GENERAL_CPRINT 9
#define GENERAL_EXACTNESS 10
#define GENERAL_CLOSE 13
typedef struct general_param {
  char pch; /*  character shown after each simulation
    just for entertainment of the user
    except for "!", then the numbers are shown
	   */
  bool 
    allowdist0, na_rm_lines, skipchecks, vdim_close_together, storing,
  /* true: intermediate results are stored: might be rather memory consuming,
         but the simulation can (depending on the method chosen) be much faster
	 when done the second time with exactly the same parameters
	 do not forget to call DeleteAllKeys when all the simulation are done.
   false: intermediate results are not stored if SimulateRF is called or
         stored only until DoGauss is called.
	 minimum memory consumption, but might be rather slow if many
	 simulations are performed with exactly the same parameters
	 to be safe call DeleteAllKeys when all the simulations are done
	 Warning! init_foo may depend on GENERAL_STORING! This may cause problems
	 if GENERAL_STORING is changed between init_foo and do_foo from false to
	 true and if do_foo is called more than once! (SimulateRF is safe in
	 this case)

	 indifferent behaviour for the simulation methods if parameters are 
	 changed after each simulation.
	 In case, only MAXFIELDS different parameter sets are used, but after each 
	 step the parameter set is changed, use different keynr for each 
	 parametere set, and STORING==true, to get fast simulation 
  */
    sp_conform, /* should the simulation result be return in 
		   as an sp class or in the old form ?
		   --- getting obsolete in future --- todo
		*/
    detailed_output,
    asList,
    returncall;
 
  int  mode, /* hightailing, fast, normal, save, pedantic */
    output, /* output mode, #alternative to mode that changes 
	       several other parameters;
	       vice versa, spConform sets 'output'
	     */
    reportcoord,
    Rprintlevel,
    Cprintlevel; /* 
		   see convert.T PL_* and RF.h PL_*
		  */
  int naturalscaling;
  /* is called PracticalRange in R (see Chiles&Delfiner?!) 
   Note: RFparameters() allows only 0 and 1 as values!

   has an effect only if cov (and not only cov.local, e.g. Brownian Motion) is 
   defined
   0 : using the covariance function as defined
   1 : rescaling of cov fctn such that cov(1)~=0.05, if rescaling function 
       does not exist then failure, e.g. for Bessel model
   2 : exact or approximate value (e.g. expPLUScirc)
   3 : MLE (special needs taken into account, long memory covariance functions
       have too long tails, which are shortened (so threshold 0.05 is changed to
       higher values in this case))
  +10: if any of the above fails : numerical evaluation!
   else : using rescaling if rescaling function exists, otherwise without 
          rescaling
  */  

  int expected_number_simu,
    every,  // every 'every' line is announced if every>0

    // bool aniso;  // currerntly cannot be changed by the user !!

  // siehe InternalCov.cc line 350, simu.cc line 394
  /* expand scale immediately to aniso;
		 this allows higher flexibility in the choice of 
		 simulation methods (tbm2 does not like dimension reduction),
		 but it is slower
	      */
    seed,
    Ttriple, 
    set;
 
  double gridtolerance, exactness;

} general_param;
#define general_START \
  {pch[NM],								\
      allowdistance0[NM], false, skipchecks[NM], false, false, true, /* 6 */ \
      false, true, false,						\
      startmode/* mode */ , output_sp, reportcoord_warnings,		\
      R_PRINTLEVEL, C_PRINTLEVEL, NAT_SCALE, 1, 0,			\
      NA_INTEGER, NA_INTEGER, 0,					\
      1e-6, exactness[NM],						\
     }


#define gaussN 6
extern const char *gauss[gaussN];
#define GAUSS_BEST_DIRECT 3
#define GAUSS_BOXCOX_OPTION 5
typedef struct gauss_param{
  double stationary_only, // logical + NA
    approx_zero, 
    boxcox[2 * MAXGAUSSVDIM];
  bool paired, loggauss;
  int direct_bestvariables;
} gauss_param;
#define gauss_START  {RF_NA, 0.05,					\
      {RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0,		\
	  RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0},	\
      false, false, 1200}

#define MINSPLITN(locmaxn) (locmaxn / 40)
#define MEDSPLITN(locmaxn) (locmaxn / 8)
#define MAXSPLITN(locmaxn) ((locmaxn * 5) / 8)

#define krigeN 5
#define KRIGE_SPLITN 2
extern const char *krige[krigeN];
typedef struct krige_param {
  bool ret_variance,  fillall;
  int locmaxn,
    locsplitn[3], // 0:when splitting is done; 1:min pts in a neighbourhood ; 2:max pts when still neighbourhoods are included
    locsplitfactor;
} krige_param;
#define krige_START {false, true, locmaxn[NM], \
      {MINSPLITN(locmaxn[NM]), MEDSPLITN(locmaxn[NM]), MAXSPLITN(locmaxn[NM])},\
      2}



#define CEN 12
extern const char *CE[CEN];
typedef struct ce_param {
  bool force, useprimes, dependent;
  char strategy;
  int trials, maxgridsize, maxmem;
  double maxGB,
     tol_re, tol_im, mmin[MAXCEDIM],
    approx_grid_step;
} ce_param;
#define ce_START				\
  {ce_force[NM], true, ce_dependent[NM],	\
      TRIVIALSTRATEGY,				\
      3, MAX_CE_MEM, MAXINT,						\
      1, ce_tolRe[NM], ce_tolIm[NM], {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, \
      ce_approx_step[NM]}

#define spectralN 4
#define SPECTRAL_PROPFACTOR 2
extern const char * spectral[spectralN];
typedef struct spectral_param {
  bool grid;
  double prop_factor, sigma;
  int lines[MAXTBMSPDIM];
} spectral_param;
#define spectral_START {sp_grid[NM], 50, 0.0, {2500, 2500, 2500, 2500}}



#define pTBMN 9
extern const char * pTBM[pTBMN];
typedef struct tbm_param {
  bool grid;
  int tbmdim, fulldim, points,
    lines[MAXTBMSPDIM];          // number of lines simulated
  double layers; // 0, 1, NA
  double linesimufactor, /*factor by which simulation on the line is finer 
			     than the grid */
    linesimustep,      /* grid lag defined absolutely */
  //  bool  tbm2num; 
    center[MAXTBMSPDIM];
} tbm_param;
#define tbm_START \
  {sp_grid[NM], -2, 3, 0, {1, 60, 500}, RF_NA, tbm_linefactor[NM], 0.0, \
   {RF_NA, RF_NA, RF_NA, RF_NA}}



#define directN 1
#define DIRECT_MAXVAR_PARAM 0
#define DIRECT_ORIG_MAXVAR 8192
extern const char *direct[directN];
typedef struct direct_param {
  //  InversionMethod inversionmethod;
  //  double svdtolerance;
  int maxvariables;
} direct_param;
#define direct_START { DIRECT_ORIG_MAXVAR } // {Cholesky, direct_tol[NM], 8192 }

#define sequN 3
extern const char * sequ[sequN];
typedef struct sequ_param{
  int max, back, initial;
} sequ_param;
#define sequ_START {5000, 10, -10}



typedef struct ave_param {
} ave_param;
#define ave_START {}


#define pnuggetN 1
extern const char * pnugget[pnuggetN];
typedef struct nugget_param {
  double tol;
} nugget_param;
#define nugget_START {nugget_tol[NM]}


#define mppN 6
extern const char * mpp[mppN];
typedef struct mpp_param {
  int  n_estim_E,
    scatter_max[MAXMPPDIM];
  double intensity[MAXMPPDIM], // intensity factor for e.g. unif_initu
    about_zero,
    shape_power,
    scatter_step[MAXMPPDIM];
} mpp_param;
#define mpp_START					\
  {50000,						\
   {NA_INTEGER, NA_INTEGER, NA_INTEGER, NA_INTEGER},			\
   {mpp_intensity[NM], mpp_intensity[NM], mpp_intensity[NM],mpp_intensity[NM]},\
     mpp_zero[NM], /* about zero */					\
       2.0,								\
	 {RF_NA, RF_NA, RF_NA, RF_NA}					\
  }


#define hyperN 4
extern const char * hyper[hyperN];
typedef struct hyper_param {
  int superpos, maxlines, mar_distr;
  double mar_param;
} hyper_param;
#define hyper_START {700, 1000, HYPER_UNIFORM, RF_NA}


#define extremeN 10
extern const char * extreme[extremeN];
#define EXTREME_FLAT 5
typedef struct extremes_param {
  int maxpoints, check_every, flat, min_n_zhou, max_n_zhou, mcmc_zhou;
  double standardmax, GEV_xi, density_ratio, eps_zhou;
} extremes_param;
#define extreme_START				\
  {MAXINT, 30, FALSE, /* FLAT_UNDETERMINED */	\
      1000, 10000000, 20,			\
      max_max_gauss[NM],			\
      1.0, 0.0, 0.01}


#define brN 9
extern const char * br[brN];
typedef struct br_param {
  int BRmaxmem, BRvertnumber, BRoptimmaxpoints, BRoptim, deltaAM;
  double BRmeshsize, BRoptimtol, variobound, corr_factor;
} br_param;
#define br_START	       \
  {10000000, 7, 10000, 2, 300, \
      0.1, 0.01, 8.0, 0.1}     \


#define distrN 9
extern const char * distr[distrN];
typedef struct distr_param{
  double safety, minsteplen, innermin, outermax;
  int maxsteps, parts, maxit, mcmc_n, repetitions;
} distr_param;
#define distr_START  {0.08, 0, 1e-20, 1e5, 1000, 8, 20, 15, 1000}  // distr (rectangular) // todo should be 500 and better algorithm for approximation!



#define MINCLIQUE(locmaxn) (locmaxn / 25)
#define MEDCLIQUE(locmaxn) (locmaxn / 5)
#define MAXCLIQUE(locmaxn) ((locmaxn * 3) / 5)

#define FLAT_UNDETERMINED -1
#define GRIDEXPAND_AVOID -1
#define fitN 42
#define FIT_MAXNEIGHBOUR 25
#define FIT_BC_LB 12
#define FIT_BC_UB 13
#define FIT_CLIQUE 26
extern const char * fit[fitN];
typedef struct fit_param{
  double bin_dist_factor, upperbound_scale_factor, lowerbound_scale_factor, 
    lowerbound_scale_LS_factor, upperbound_var_factor, lowerbound_var_factor, 
  //lowerbound_sill, 
    scale_max_relative_factor, minbounddistance, 
    minboundreldist,   minmixedvar, maxmixedvar, 
    BC_lambdaLB[2 * MAXGAUSSVDIM], BC_lambdaUB[2 * MAXGAUSSVDIM], scale_ratio,
    min_diag, pgtol, pgtol_recall,
    factr, factr_recall;// logical + NA // 19
  int approximate_functioncalls, bins, nphi, 
    ntheta, ntime, critical, 
    n_crit, locmaxn, 
    locsplitn[3], // 0:when splitting is done; 1:min pts in a neighbourhood ; 2:max pts when still neighbourhoods are included
    locsplitfactor, smalldataset, optimiser,
    algorithm, likelihood, estimate_variance,
    split; // 16
  bool  use_naturalscaling, onlyuser, 
   reoptimise, ratiotest_approx,
    split_refined, cross_refit;      // 8
  char lengthshortname; // 1..255
} fit_param;
#define fit_START\
  {0.4, 3.0, 5.0, 3.0, 10.0, 1e4, /* 6 */				\
      1000.0, 0.001, 0.02, 1/1000, 1000,  /* 11 */			\
      {-10, -10.0, -10, -10.0, -10, -10.0, -10, -10.0, -10, -10.0,	\
	  -10, -10.0, -10, -10.0, -10, -10.0, -10, -10.0, -10, -10.0},	\
	{10, 10.0, 10, 10.0, 10, 10.0, 10, 10.0, 10, 10.0,		\
	    10, 10.0, 10, 10.0, 10, 10.0, 10, 10.0, 10, 10.0},		\
       0.1, 1e-7, fit_pgtol[NM], fit_pgtol_recall[NM], /* 17 */	\
      fit_factr[NM], fit_factr_recall[NM],				\
      /* int: */							\
      50, 20, 1, 1, 20, 0 /* critical */,            /* 6 */		\
	 5 /* ncrit */ , maxclique[NM],					\
	 {MINCLIQUE(maxclique[NM]), MEDCLIQUE(maxclique[NM]),\
	     MAXCLIQUE(maxclique[NM]) }, 2, 2000,  /* 11 */		\
      0, -1 /* algorithm */, 0, NA_INTEGER, fit_split[NM],	/* 16 */ \
      !true, false, fit_reoptimise[NM],  /* 5 */\
      fit_ratiotest_approx[NM], true, fit_cross_refit[NM], \
   12} // fit


#define empvarioN 5
extern const char * empvario[empvarioN];
typedef struct empvario_param{
  double phi0, theta0, tol;
  bool pseudovariogram, fft;
} empvario_param;
#define empvario_START {0.0, 0.0, 1e-13, false, true}

#define guiN 3
extern const char * gui[guiN];
#define GUI_SIZE 2
typedef struct gui_param{
  bool alwaysSimulate;
  int method, size[2];
} gui_param;
#define gui_START {true, CircEmbed, {1024, 64}}

#define graphicsN 9
extern const char *graphics[graphicsN];
#define GRAPHICS_UPTO 3
typedef struct graphics_param{
  bool onefile;
  double height, resolution;
  int PL, increase_upto[2], always_open, always_close, number;
  char filename[100];
} graphics_param;
#define graphics_START  {false,  6.0, 72, 1, {3, 4}, true, false, 0, ""}


#define registersN 3
extern const char *registers[registersN];
typedef struct registers_param {
  int keynr, predict, likelihood;
} registers_param;
#define register_START {0, MODEL_PREDICT, MODEL_USER}


#define internalN 20
extern const char * internals[internalN];
#define INTERNALS_NEWANISO 2
#define INTERNALS_ONGRID 9
#define INTERNALS_COORD_CHANGE 12
#define INTERNALS_ZENIT 14
typedef struct internal_param{ 
  // if changed, CHANGE ALSO RestWarnings in 'userinterfaces.cc';
  bool warn_oldstyle, warn_newstyle, warn_Aniso, 
    warn_ambiguous, warn_normal_mode, warn_mode, 
    stored_init, warn_scale, warn_coordinates,
    warn_on_grid, warn_new_definitions, warn_aspect_ratio, 
    warn_coord_change, warn_color_palette, warn_zenit,
    do_tests, warn_constant, warn_negvar, warn_onlyvar
    ;//   
  int examples_reduced;
} internal_param;
#define internal_START\
  {true, true, true, false, true,    \
      true, false, true, true, true, \
      true, true, true, true, true,  \
      DO_TESTS, true, true, true, 0}

#define coordsN 10
#define COORDS_XYZNOTATION 0
#define ZENIT 8
extern const char *coords[coordsN];
typedef struct coords_param{
  double xyz_notation, zenit[2];
  coord_sys_enum coord_system, new_coord_system;
  char newunits[MAXUNITS][MAXUNITSCHAR], // only to read
    curunits[MAXUNITS][MAXUNITSCHAR],
    varunits[MAXUNITS][MAXUNITSCHAR],
  // 2 user variables for data.frames
    data_names[MAXUNITS][MAXCHAR],///
    x_names[MAXUNITS][MAXCHAR];
  // auxiliary variables
  int data_nr_names, x_nr_names, // needed by data_names and x_names
    data_idx[2], x_idx[2]; // alternatives for data_names and x_namwe
  bool polar_coord;
} coords_param;
#define coords_START\
  {RF_NA, {1, RF_NA}, coord_auto, coord_keep, \
  {""}, {""}, {""}, {""}, {""}, 0, 0, {NA_INTEGER, NA_INTEGER}, \
   {NA_INTEGER, NA_INTEGER}, false}


#define specialN 1
extern const char * special[specialN];
typedef struct special_param {
   int multcopies;
} special_param;
#define special_START {10}



#define obsoleteN 10
extern const char * obsolete[obsoleteN];


typedef struct globalparam{
  general_param general;
  gauss_param gauss;
  krige_param krige;
  ce_param ce;
  spectral_param spectral;
  tbm_param tbm;
  direct_param direct;
  sequ_param sequ;
  ave_param ave;
  nugget_param nugget;
  mpp_param mpp;
  hyper_param hyper;
  extremes_param extreme;
  br_param br;
  distr_param distr;
  fit_param fit;
  empvario_param empvario;
  gui_param gui;
  graphics_param graphics;
  registers_param registers;
  internal_param internal;
  coords_param coords;
  special_param special;
  solve_param solve;
} globalparam;
extern globalparam GLOBAL;


#endif
