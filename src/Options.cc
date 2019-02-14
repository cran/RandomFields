
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 20017 -- 2017 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

//#include <R.h>
//#include <Rdefines.h>
//#include <R_ext/Linpack.h>
//#include <Rmath.h>  
//#include <stdio.h>  
//#include <unistd.h>
//#include <string.h>
//#include "RF.h"
//#include "primitive.h"
//#include "kleinkram.h"
// #include "Operator.h"
#include "RF.h"
#include "xport_import.h"
#include "init.h"

#define ADD(ELT) SET_VECTOR_ELT(sublist, k++, ELT);
#define ADDCHAR(ELT) x[0] = ELT; ADD(ScalarString(mkChar(x)));


#define nOptimVar 4
const char * OPTIM_VAR_NAMES[nOptimVar] = 
  {"never", "respect bounds", "try", "yes"}; // keep yes last !!


void ResetWarnings(int * allwarnings) {
  internal_param *w = &(GLOBAL.internal);
  w->warn_oldstyle = w->warn_newstyle = w->warn_Aniso = 
    w->warn_normal_mode = w->warn_mode = w->warn_coordinates =
    w->warn_on_grid = w->warn_new_definitions = w->warn_aspect_ratio = 
    w->warn_coord_change = w->warn_color_palette = w->warn_zenit =
    w->warn_constant = w->warn_negvar = w->warn_onlyvar =
    w->warn_parallel = true;
  w->allow_duplicated_loc = false;   
  w->warn_mathdef =  Nan;
  if (*allwarnings) w->warn_ambiguous = true;
}


#define startmode normal
#define NM ((int) startmode)   /* normal mode */
//       {careless, sloppy, easygoing, normal, precise, pedantic, neurotic}
bool 
  allowdistance0[nr_modes] = {true, false, false, false, false, false, false},
  skipchecks[nr_modes] =  {true, false, false, false, false, false, false},
  ce_force[nr_modes] =       {true, true, true, false, false, false, false},
  ce_dependent[nr_modes] =   {true, true, true, false, false, false, false},
  sp_grid[nr_modes] =           {true, true, true, true, true, false, false},
  fit_reoptimise[nr_modes] = {false, false, false, false, true, true, true},
  fit_ratiotest_approx[nr_modes]={true, true, true, true, false, false, false},
  fit_cross_refit[nr_modes] = {false, false, false, false, true, true, true}
  ;
char pch[nr_modes] =         {'\0',  '\0',  '\0',  '.',   '.',   '.',   '.'}
  ;
int locmaxn[nr_modes] =      {3000, 4000, 6000, 8000, 9000, 10000000, 10000000},
  ce_trials[nr_modes] =      {1,    2,    3,    3,    4,    5,   6},
  spectral_lines[nr_modes] = {300, 500, 1500,  2500, 3500, 5000, 10000},
  tbm_lines[nr_modes] =      {30,   40,   50,   60,   80,   100,   200},
  mpp_n_estim_E[nr_modes] =  {10000, 10000, 20000,50000,80000,100000,1000000},
  hyper_superpos[nr_modes] = {200, 300, 450, 700, 1000, 2000, 5000},
  maxclique[nr_modes] =      {1000, 2000, 3500, 5000, 6000, 7000, 10000000},
  fit_critical[nr_modes] =   {-1,  -1,   -1,   0,   1,   2,   3}, 
  fit_ncrit[nr_modes] =      { 2,   4,    5,   5,   5,  10,  20}, 
  //  fit_optim_var[nr_modes] =  { 2,   2,    2,   2,   1,   1,   1} ,
  fit_split[nr_modes] =      {10000,  5,  4, 4, 4, 3, 3},
  sequ_back_initial[nr_modes]={3,   3,    5,  10,  20,  30,  40}
  ;

usr_bool exactness[nr_modes] = {False, False, False, Nan, True, True, True};

double 
  ce_tolRe[nr_modes] =       {-1e2,  -1e1,  -1e-3, -1e-7,  -1e-14, 0, 0},
  ce_tolIm[nr_modes] =       {1e2,  1e1,  1e-1,  1e-3,   1e-7, 0, 0},
  ce_approx_step[nr_modes] = {1.0,  1.0,   RF_NA, RF_NA, RF_NA, 0.0, 0.0},
#ifdef SCHLATHERS_MACHINE
  svd_tol[nr_modes] =     {1e-5, 1e-6, 1e-7, 1e-8, 1e-10, 1e-12, 0},
#else
  svd_tol[nr_modes] =     {0, 0, 0, 0, 0, 0, 0},
#endif
  nugget_tol[nr_modes] =     {1e-8, 1e-8, 1e-12, 0,     0,     0,     0},
  tbm_linefactor[nr_modes] = {1.5,  1.5,  1.7,   2.0,  3.0,  5.0,  6.0},
  mpp_intensity[nr_modes] =  {50, 50, 80, 100, 200, 500, 1000},
  mpp_zero[nr_modes] =       {1e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7},
  max_max_gauss[nr_modes] =  {2,    2,    3,    4,    5,   6,   6},
  fit_pgtol[nr_modes]       ={10, 1e-2, 1e-3, 1e-4, 1e-4, 1e-6, 0},
  fit_pgtol_recall[nr_modes]={10, 1e-1, 1e-2, 1e-3, 1e-3, 1e-4, 0},
  fit_factr[nr_modes]       ={1e14, 1e13, 1e12, 1e11, 1e11, 1e9, 1e7},
  fit_factr_recall[nr_modes]={1e15, 1e14, 1e13, 1e12, 1e12, 1e10, 1e8}
   ;

const char *f_opt[nr_modes] = {"optim", "optim", "optim", "optim", "optim", "optim", "optim"}; // to do optimx


int 
PL = 1,
  CORES = INITCORES, // FINITE 
  NS = NAT_SCALE;
globalparam GLOBAL = {
  general_START,
  gauss_START,
  krige_START, 
  ce_START,
  spectral_START,
   tbm_START,
  direct_START,
  sequ_START,
  ave_START,
  nugget_START,
  mpp_START,

  hyper_START,
  extreme_START,
  br_START,
  distr_START,
  fit_START, 
  
  empvario_START, 
  gui_START,  
  graphics_START,
  register_START,
  internal_START,
  coords_START, 
  special_START
  
};
utilsparam *GLOBAL_UTILS;
  

void SetDefaultOutputModeValues(output_modes mode){
  general_param *gp = &(GLOBAL.general);
  gp->output = mode;
  gp->sp_conform = mode == output_sp;
  gp->returncall = mode == output_geor;
  gp->reportcoord = 
    mode == output_geor ? reportcoord_always : reportcoord_warnings;
  /*
  switch(mode) {
  case output_sp : 
     break;
  case output_rf :
  break;
  case output_geor :
     break;
  default: BUG;
  }
  */
}


void SetDefaultModeValues(int old, int m){
  //                      high  fast, normal, save, pedantic
  GLOBAL_UTILS->basic.skipchecks = skipchecks[m];
  GLOBAL.general.pch = pch[m];
  GLOBAL.general.exactness = exactness[m];
  GLOBAL.general.allowdist0 = allowdistance0[m];
  GLOBAL.krige.locmaxn = locmaxn[m];
  GLOBAL.krige.locsplitn[0]  = MINSPLITN(locmaxn[m]);
  GLOBAL.krige.locsplitn[1]  = MEDSPLITN(locmaxn[m]);
  GLOBAL.krige.locsplitn[2]  = MAXSPLITN(locmaxn[m]);
  GLOBAL.ce.force = ce_force[m];
  GLOBAL.ce.tol_im = ce_tolIm[m];
  GLOBAL.ce.tol_re = ce_tolRe[m];
  GLOBAL.ce.trials = ce_trials[m];
  GLOBAL.ce.dependent = ce_dependent[m];
  GLOBAL.ce.approx_grid_step = ce_approx_step[m];
  GLOBAL_UTILS->solve.svd_tol = svd_tol[m];
  GLOBAL.nugget.tol = nugget_tol[m];
  for(int i=0; i<4; GLOBAL.spectral.lines[i++] = spectral_lines[m]);
  GLOBAL.spectral.grid = sp_grid[m];
  GLOBAL.tbm.lines[1] = tbm_lines[m];
  GLOBAL.tbm.lines[2] = (tbm_lines[m] * 25) / 6;
  GLOBAL.tbm.linesimufactor = tbm_linefactor[m];
  GLOBAL.tbm.grid = sp_grid[m];
  GLOBAL.mpp.n_estim_E = mpp_n_estim_E[m];
  for(int i=0; i<4; GLOBAL.mpp.intensity[i++] = mpp_intensity[m]);
  GLOBAL.mpp.about_zero = mpp_zero[m];
  GLOBAL.hyper.superpos = hyper_superpos[m];
  GLOBAL.extreme.standardmax = max_max_gauss[m];
  GLOBAL.fit.locmaxn = maxclique[m];
  GLOBAL.fit.locsplitn[0] = MINCLIQUE(maxclique[m]);
  GLOBAL.fit.locsplitn[1] = MEDCLIQUE(maxclique[m]);
  GLOBAL.fit.locsplitn[2] = MAXCLIQUE(maxclique[m]);
  fit_param *f = &(GLOBAL.fit);
  f->split = fit_split[m];
  f->reoptimise = fit_reoptimise[m];
  f->ratiotest_approx = fit_ratiotest_approx[m];
  f->cross_refit = fit_cross_refit[m];
  f->critical = fit_critical[m];
  f->n_crit = fit_ncrit[m];
  f->pgtol = fit_pgtol[m];
  f->pgtol_recall = fit_pgtol_recall[m];
  f->factr = fit_factr[m];
  f->factr_recall = fit_factr_recall[m];
  //  f->optim_var_estim = fit_optim_var[m];
   char dummy[10];
  STRCPY(dummy, f_opt[m]);
  f->optimiser = Match(dummy, OPTIMISER_NAMES, nOptimiser);
  GLOBAL.sequ.initial = -(GLOBAL.sequ.back = sequ_back_initial[m]);

  //  printf("optimiser %d %.50s\n", f->optimiser,  OPTIMISER_NAMES[f->optimiser]);

  internal_param *w = &(GLOBAL.internal);
  w->stored_init = false;
  if (m < normal) {
    w->warn_Aniso = w->warn_ambiguous = 
      w->warn_on_grid = w->warn_aspect_ratio = w->warn_color_palette = false;
    w->warn_mathdef = False;
  } else if (m>old) { 
    w->warn_oldstyle = w->warn_Aniso = w->warn_new_definitions = 
      w->warn_aspect_ratio = true;
    w->warn_mathdef = Nan;
    if (m > normal) w->warn_ambiguous = w->warn_negvar = true;
  }
  if (m != normal && w->warn_mode) {
    // hier:  Syscall param sched_setaffinity(mask) points to unaddressable byte
     PRINTF("Note that the option '%.50s' is still in an experimental stage, so that the behaviour may change (slightly) in future.\n", general[GENERAL_MODUS]);
    w->warn_mode = false;
  }
}
 

void getUnits(SEXP el, char VARIABLE_IS_NOT_USED *name, 
	      char units[MAXCOOORDNAMES][MAXUNITSCHAR], 
	      char units2[MAXCOOORDNAMES][MAXUNITSCHAR]) {  
  int i, j, 
    l = length(el);
  if (TYPEOF(el) != NILSXP && TYPEOF(el) == STRSXP && l >= 1) {
    for (i=j=0; i<MAXCOOORDNAMES; i++, j=(j + 1) % l) {
      strcopyN(units[i], CHAR(STRING_ELT(el, j)), MAXUNITSCHAR);
      if (units2!=NULL) {
	strcopyN(units2[i], CHAR(STRING_ELT(el, j)), MAXUNITSCHAR);
      }
    }
  } else RFERROR("invalid units");
}


SEXP UNITS(char units[MAXCOOORDNAMES][MAXUNITSCHAR]) {
  SEXP unitnames;
  int nn;
  PROTECT(unitnames = allocVector(STRSXP, MAXCOOORDNAMES)); 
  for (nn=0; nn<MAXCOOORDNAMES; nn++) {
    SET_STRING_ELT(unitnames, nn, mkChar(units[nn]));
  }
  UNPROTECT(1);
  return unitnames;
}


void CE_set(SEXP el, int j, char *name, ce_param *cp, bool isList) {
  switch(j) {
  case 0: cp->force = LOGI; break;
  case 1:
    Real(el, name, cp->mmin, MAXCEDIM) ;
    int d;
    for (d=0; d<MAXCEDIM; d++) {
      if (cp->mmin[d]<0.0) {
	if ((cp->mmin[d]>-1.0)) {
	  cp->mmin[d] = - 1.0;
	  WARN1("'%.50s' set to -1.0.\n", name);
	}
      }
    }
    break;
  case 2: int strat;
    strat = INT;
    if (strat>LASTSTRATEGY) {  
      WARN2("%.50s <= %d not satisfied\n", name, LASTSTRATEGY);
    } else cp->strategy= (char) strat;          
    break;
  case 3: {
    cp->maxGB = POSNUM;
    if (!isList) cp->maxmem = MAXINT;	
  } break;
  case 4: {
    cp->maxmem = POSINT;  
    if (!isList) cp->maxGB = RF_INF;	
  }  break;
  case 5: cp->tol_im = POS0NUM; break;
  case 6: cp->tol_re = NEG0NUM; break;
  case 7: cp->trials = NUM;
    if (cp->trials<1) {
      cp->trials=1;
      WARN1("%.50s is set to 1\n", name);
    }
    break;
  case 8: cp->useprimes = LOGI; break;
  case 9: cp->dependent = LOGI; break;
  case 10: cp->approx_grid_step = POS0NUM; break;
  case 11: cp->maxgridsize = POS0INT; break;
  default: RFERROR("unknown parameter for GLOBAL.general");
  }
}
 

 

#define prefixN 23
const char * prefixlist[prefixN] = 
  {"general", "gauss", "krige", 
   "circulant", "direct",  "nugget",
   "sequ", "spectral", "tbm",
   "mpp", "hyper", "maxstable",  // 11
   "br", "distr", "fit",         // 14
   "empvario", "gui", "graphics",// 17 
   "registers", "internal", "coords", // 20
   "special", OBSOLETENAME  // 22
};


// IMPORTANT: all names of general must be at least 3 letters long !!!
const char *general[generalN] =
  { "modus_operandi", "storing", "every", "gridtolerance", // 3
    "pch", "practicalrange", "spConform",  // 6
    "exactness", "allowdistanceZero", "na_rm_lines",
    "vdim_close_together", "expected_number_simu", // 11
    "detailed_output", "Ttriple",
    "returncall", "output", "reportcoord", "set", // 17
    "seed_incr"
};


const char *gauss[gaussN]= {"paired", "stationary_only", "approx_zero", 
			    "direct_bestvar", "loggauss", "boxcox"};

const char *krige[krigeN] = { "return_variance", "locmaxn",
			     "locsplitn", "locsplitfactor", "fillall"};

const char *CE[CEN] = {"force", "mmin", "strategy", "maxGB",
		       "maxmem", "tolIm","tolRe", 
		       "trials", "useprimes", "dependent", 
		       "approx_step", "approx_maxgrid"};

const char *direct[directN] = {//"root_method", "svdtolerance", 
			       "max_variab"};

const char * pnugget[pnuggetN] ={"tol"};

const char * sequ[sequN] ={"back_steps", "initial"};

const char * spectral[spectralN] = {"sp_lines", "sp_grid",  
				    "prop_factor", "sigma"};

const char * pTBM[pTBMN] = {"reduceddim", "fulldim", "center", 
			    "points", "lines", "linessimufactor", 
			    "linesimustep", "layers", "grid"};

const char * mpp[mppN] = {"n_estim_E", // n to determine E by simulation
			  "intensity", 
			  // "refradius_factor", 
			  "about_zero", "shape_power",
			  "scatter_max", "scatter_step",
			  // "plus", 
			  // "samplingdist", "samplingr",// MPP_cc
			  //"p", // Gneiting_cc
                         };

const char * hyper[hyperN] = {"superpos", "maxlines", "mar_distr",
			      "mar_param"};

const char * extreme[extremeN] = 
  {"max_gauss", "maxpoints", "xi",
   "density_ratio", "check_every", "flathull",
   "min_n_zhou", "max_n_zhou", "eps_zhou",
   "mcmc_zhou", "min_shape_gumbel", "scatter_method"};

const char * br[brN] = 
  {"maxtrendmem", "meshsize", 
   "vertnumber", "optim_mixed", "optim_mixed_tol", 
   "variobound", "deltaAM"};

const char * distr[distrN] = 
  {"safety", "minsteplen", "maxsteps", 
   "parts", "maxit",  "innermin", 
   "outermax", "mcmc_n", "repetitions"};

const char * fit[fitN] = 
  {"bin_dist_factor", "upperbound_scale_factor", "lowerbound_scale_factor", //2
   "lowerbound_scale_ls_factor","upperbound_var_factor","lowerbound_var_factor",
   //"lowerbound_sill",
   "scale_max_relative_factor", "minbounddistance", "minboundreldist",  //8
   "approximate_functioncalls",  "boxcox_lb",  "boxcox_ub", //11
   "use_naturalscaling",  "bins",  "nphi", //14
   "ntheta", "ntime", "only_users", //17
   "shortnamelength",  "split",   "scale_ratio", //20
   "critical",  "n_crit",  "max_neighbours",  //23
   "cliquesize",  "splitfactor_neighbours",  "smalldataset", 
   "min_diag",  "reoptimise",   "optimiser", 
   "algorithm",  "likelihood",  "ratiotest_approx", 
   "split_refined",  "cross_refit",  "estimate_variance_globally", 
   "pgtol", "pgtol_recall",  "factr",
   "factr_recall" ,"use_pseudovariogram"
  };


const char * empvario[empvarioN] = 
  {"phi0", "theta0", "tol0",
   "pseudovariogram", "fft", "halvenangles"};

const char * gui[guiN] = 
  {"alwaysSimulate", "simu_method", "size"};



const char *graphics[graphicsN]= 
  {"close_screen" ,"grPrintlevel", "height", 
   "increase_upto", "always_open_device", "file", 
   "onefile", "filenumber", "resolution", "split_screen", "width",
   "always_close_device", "grDefault"};

const char *registers[registersN] = 
  {"register", "predict_register", "likelihood_register"};

const char * internals[internalN] =  {
  // Achtung ! warn parameter werden nicht pauschal zurueckgesetzt
  "warn_oldstyle", "warn_newstyle", "warn_newAniso", 
  "warn_ambiguous", "warn_normal_mode",  "warn_mode",
  "stored.init",  "warn_scale", "warn_coordinates",
  "warn_on_grid", "warn_new_definitions", "warn_aspect_ratio",
  "warn_coord_change", "warn_colour_palette", "warn_missing_zenit",
  "do_tests", "warn_constant", "warn_negvar",
  "warn_onlyvar", "warn_definite",  "examples_reduced",
  "warn_parallel", "allow_duplicated_locations", "warn_seed",
  "declare_PL", // unused ?!!
  "warn_modus_operandi"
};

const char *coords[coordsN] =
  { "xyz_notation", "coord_system", "new_coordunits", "coordunits", "varunits",
    "varidx", "varnames", "coordidx", "coordnames", "new_coord_system",
    "zenit", // "radius",
    "polar_coord", "angleunits", "allow_earth2cart"};

const char * special[specialN] = {"multicopies"};

const char * obsolete[obsoleteN] = 
  { "oldstyle", "newstyle",  "newAniso", 
    "ambiguous", "normal_mode", "solvesigma",
    "optim_var_elimination", "sill",
    "matrix_inversion", "matrix_tolerance"
};



const char **all[prefixN] = {general, gauss, krige, CE, direct, 
			     pnugget, sequ, spectral, pTBM, mpp,
			     hyper, extreme, br, distr, fit, 
			     empvario, gui, graphics, registers, internals, 
			     coords, special, obsolete};

int allN[prefixN] = {generalN, gaussN, krigeN, CEN, directN,
		 pnuggetN,  sequN, spectralN, pTBMN, mppN,
		 hyperN, extremeN, brN, distrN, fitN, 
		 empvarioN, guiN, graphicsN, registersN, internalN, 
		 coordsN, specialN, obsoleteN};


globalparam *WhichOptionList(int local) {  
  globalparam *ans = &GLOBAL;
  if (local != isGLOBAL) {
    KEY_type *KT = KEYT();
    ans = &(KT->global);
  }
  return ans;
}

void setparameter(int i, int j, SEXP el, char name[200], bool isList,
		  int local) {
  //  printf("setparameter %.50s %d %d\n", name, i, j);
  // printf("sp %d %d ", i, j);

  globalparam *options = WhichOptionList(local);
  //bool isList = options == &GLOBAL;
  
  switch(i) {
  case 0: {// general
    general_param *gp;
    gp = &(options->general);
    switch(j) {
    case GENERAL_MODUS: {
      int old_mode = gp->mode;
      gp->mode = GetName(el, name, MODE_NAMES, nr_modes, normal);
#ifndef SCHLATHERS_MACHINE
      if (old_mode != gp->mode && GLOBAL.internal.warn_modus_operandi) {
	  PRINTF("\n"); //verstehe wer will -- es gab hier eine valgrind warning
	 PRINTF("Note that behaviour of 'modus_operandi' has changed within 'RFfit' in version 3.1.0 of RandomFields. Roughly:\nwhat was called 'careless' is now called 'sloppy';\n");
PRINTF("what was called 'sloppy' is now called 'easygoing';\nwhat was called 'easygoing' is now called 'normal';\nwhat was called 'normal' is now called 'precise';\netc.\n");
	GLOBAL.internal.warn_modus_operandi = false;
   }
#endif  
     SetDefaultModeValues(old_mode, gp->mode);
  }
       break;
   case GENERAL_STORING: {
      gp->storing = LOGI;
      if (!gp->storing) {
	if  (length(el) > 1) {
	  int len = length(el);
	  for (int elnr=1; elnr < len; elnr++) {
	    int nr = Integer(el, (char*) "storing (register)", elnr);
	    if (nr != NA_INTEGER && nr>=0 && nr<=MODEL_MAX) {
	      model **key = KEY() + nr;
	      if (*key != NULL) COV_DELETE(key);
	    }
	  }
	} else {	  
	  // delete all keys
	  for (int kn=0; kn<PIDMODULUS; kn++) {
	    KEY_type *KT = PIDKEY[kn];
	    while (KT != NULL) {
	      KEY_type *q = KT;
	      KT = KT->next;
	      KEY_type_DELETE(&q);
	    }
	    PIDKEY[kn] = NULL;
	  }
	}
      }
   }
      break;
   case 2: gp->every = POS0INT;      break;
    case 3: gp->gridtolerance = NUM; break;
    case 4: gp->pch = CHR;           break;
    case 5: 
      int n;
      n =INT;  
      if (n!=0 && n!=1 && n!=2 && n!=3)
	RFERROR("PracticalRange out of range. It should be TRUE or FALSE.");
      NS = gp->naturalscaling = n;
      break;
    case 6: 
      SetDefaultOutputModeValues(LOGI ? output_sp : output_rf);
      break;
    case GENERAL_EXACTNESS: gp->exactness = USRLOG; break; 
    case 8: gp->allowdist0 = LOGI; break;
    case 9: gp->na_rm_lines = LOGI; break;
    case GENERAL_CLOSE: gp->vdim_close_together = LOGI;    
      if (gp->vdim_close_together) {
	gp->vdim_close_together = false;
 	RFERROR1("'%.50s' will likely not be programmed", general[GENERAL_CLOSE]);
	// Achtung! gausslikeli.cc waere stark betroffen!!
      }
      break;
    case 11: gp->expected_number_simu = POS0INT; break;
    case 12: gp->detailed_output = LOGI; break;
    case 13: gp->Ttriple = INT; break;
    case 14: gp->returncall = LOGI; break;
    case 15:
      SetDefaultOutputModeValues((output_modes) 
				 GetName(el, name, OUTPUTMODE_NAMES, 
					 nr_output_modes, gp->output));
      break;
    case GENERAL_REPORTCOORD:
      gp->reportcoord = GetName(el, name, REPORTCOORD_NAMES, 
				nr_reportcoord_modes, gp->reportcoord);
      break;      
    case 17: gp->set = POSINT - 1; break;
    case 18: gp->seed_incr = INT; break;
    default: BUG;
    }}
    break;
  case 1: { // gauss
    gauss_param *gp;
    gp = &(options->gauss);
    switch(j) {
    case 0: gp->paired = LOGI; break;
    case 1: gp->stationary_only = USRLOG; break;
    case 2: gp->approx_zero = POS0NUM; break;
    case GAUSS_BEST_DIRECT: gp->direct_bestvariables = POS0INT;    break;
    case 4: gp->loggauss = LOGI;
      if (gp->loggauss) gp->boxcox[0] = RF_INF;
      break;
    case GAUSS_BOXCOX_OPTION: {
      int len = length(el);
      double boxcox[2 * MAXBOXCOXVDIM];      
      if (len < 1 || (len > 2 && len % 2 != 0)) 
	RFERROR1("'%.50s' not a scalar or vector/matrix of multiple of length 2", 
	     gauss[GAUSS_BOXCOX_OPTION]);
      boxcox[1] = 0.0;
      for (int L=0; L<len; L++) boxcox[L] = Real(el, name, L);
      if (len == 1) len = 2;
      for (int k=0; k<len; k++) {
	if (R_FINITE(boxcox[k])) {
	  if (boxcox[k] < options->fit.BC_lambdaLB[k]) 
	    RFERROR6("%10g=%.50s[%d] < %.50s[%d]=%10g.", 
		     boxcox[k], gauss[GAUSS_BOXCOX_OPTION], k, fit[FIT_BC_LB],k,
		     options->fit.BC_lambdaLB[k]);
	  if (boxcox[k] > options->fit.BC_lambdaUB[k]) 
	    RFERROR6("%10g=%.50s[%d] > %.50s[%d]=%10g.", 
		     boxcox[k],gauss[GAUSS_BOXCOX_OPTION], k, fit[FIT_BC_UB], k,
		     options->fit.BC_lambdaUB[k]); 
	}
      }
      for (int k=0; k<len; k++) gp->boxcox[k] = boxcox[k];
      if (ISNA(gp->boxcox[0]) || gp->boxcox[0] != RF_INF) gp->loggauss = false;	
    }
      break;
    default: BUG;
    }}
    break;
  case 2: { // krige
    krige_param *kp;
    kp = &(options->krige);
    switch(j) { 
    case 0: kp->ret_variance = LOGI; break;
    case 1: kp->locmaxn = POS0INT; break;
    case KRIGE_SPLITN: {
      int  newval[3],
	len = length(el);
      if (len < 3) RFERROR1("'%.50s' must have 3 components", krige[KRIGE_SPLITN]);
      for (int ii=0; ii<len; ii++) {
	newval[ii] = Integer(el, name, ii); 
      }
      if (newval[0] > newval[1] || newval[1] > newval[2]) {
	RFERROR6("%.50s[1] <= %.50s[2] <= %.50s[3] not satisfied [ %d %d %d ]",
	     krige[KRIGE_SPLITN], krige[KRIGE_SPLITN], krige[KRIGE_SPLITN],
	     newval[0], newval[1], newval[2]
	     ); 
      } 
      for (int ii=0; ii<len; ii++) kp->locsplitn[ii] = newval[ii];
    }
      break;
    case 3: kp->locsplitfactor = POS0INT; break;
    case 4: kp->fillall = LOGI; break;
    default: BUG;
    }}
    break;
  case 3: // CE
    CE_set(el, j, name, &(options->ce), isList); 
    break;
  case 4: { //direct
    direct_param *dp;
    dp = &(options->direct);
    switch(j) {
    case DIRECT_MAXVAR_PARAM : {
      int mv = POS0INT;   
 
      if (mv <= options->gauss.direct_bestvariables) {
	RFERROR3("'%.50s' must be greater than '%.50s' = %d.\n", 
	     direct[DIRECT_MAXVAR_PARAM], gauss[GAUSS_BEST_DIRECT], 
	     options->gauss.direct_bestvariables);
      }
#define MAX_DIRECT_MAXVAR 30000
      if (mv > MAX_DIRECT_MAXVAR && dp->maxvariables <= MAX_DIRECT_MAXVAR) {
	WARN2("'%.50s' should better not exceed %d.\n", 
	     direct[DIRECT_MAXVAR_PARAM], MAX_DIRECT_MAXVAR);
      }
      dp->maxvariables = mv;
    }
     break;
    default: BUG;
    }}
    break;
  case 5:  {// pnugget, 
    nugget_param *np;
    np = &(options->nugget) ;
    switch(j) {
    case 0: np->tol = POS0NUM; break;
    default: BUG;
    }}
    break;
  case 6: {//sequ,
    sequ_param *sp;
    sp = &(options->sequ) ;
    switch(j) {
    case 0: sp->back = INT; if (sp->back < 1) sp->back=1;   break;
    case 1: sp->initial = INT; break;
    default: BUG;
    }}
    break;
  case 7: { // spectral,
    spectral_param *sp;
    sp = &(options->spectral) ;
    switch(j) {
    case 0: Integer(el, name, sp->lines, MAXTBMSPDIM); break;
    case 1: sp->grid = LOGI; break;
    case SPECTRAL_PROPFACTOR: sp->prop_factor = POS0NUM;
      if (sp->prop_factor <= 0.1) {
	sp->prop_factor=0.1;
	WARN1("'%.50s' less than 0.1. Set to 0.1.", spectral[SPECTRAL_PROPFACTOR]);
      }
      break;
    case 3: sp->sigma = NUM; break;
    default:  BUG;
    }}
    break;
  case 8: {// TBM
    tbm_param *tp;
    tp = &(options->tbm) ;
    switch(j) {
    case 0: tp->tbmdim = INT; break;
    case 1: tp->fulldim = POS0INT; break;
    case 2: Real(el, name, tp->center, MAXTBMSPDIM);  break;
    case 3: tp->points = POS0INT; break;
    case 4: Integer(el, name, tp->lines, MAXTBMDIM); break;
    case 5: 
      tp->linesimufactor = Real(el, name, 0);
      
      if (tp->linesimufactor < 0.0) {
	warning("Both 'linesimufactor' and 'linesimustep' must be non-negative.");
	tp->linesimufactor = 0.0;
      }
      if (tp->linesimufactor>0.0 && tp->linesimustep>0.0) tp->linesimustep=0.0;
      break;
    case 6: 
      tp->linesimustep = Real(el, name, 0);  
      if (tp->linesimustep < 0.0) {
	warning("Both 'linesimufactor' and 'linesimustep' must be non-negative.");
	tp->linesimustep = 0.0;
      }
      if (tp->linesimufactor>0.0 && tp->linesimustep>0.0) 
	tp->linesimufactor=0.0;
      break;
    case 7:
      tp->layers = UsrBool(el, name, 0);
      break;      
    case 8: 
      tp->grid = LOGI; break;
    default:  BUG;
    }}
    break;
  case 9: {//  mpp, 
    mpp_param *mp;
    mp = &(options->mpp);
    switch(j) {
    case 0: mp->n_estim_E = POS0INT; break;
    case 1: Real(el, name, mp->intensity, MAXMPPDIM); break;
      // case 2: mp->refradius_factor = POS0NUM; break;
    case 2: mp->about_zero = POS0NUM; break;
    case 3: mp->shape_power = NUM; break;
    case 4: {
      Integer(el, name, mp->scatter_max, MAXMPPDIM) ;
    }
      break;
    case 5: {
      Real(el, name, mp->scatter_step, MAXMPPDIM) ;
    }
      break;
    default: BUG;
   }}
    break;
  case 10: {//hyper, 
    hyper_param *hp;
    hp = &(options->hyper);
    switch(j) {
    case 0: hp->superpos = POS0INT;  break;
    case 1: hp->maxlines = POS0INT;  break;
    case 2: hp->mar_distr = INT; break;
    case 3: hp->mar_param = NUM; break;
    default:  BUG;
    }}
    break;
  case 11: {//		 extreme
    extremes_param *ep;
    ep = &(options->extreme);
    switch(j) {
    case 0: ep->standardmax = POS0NUM; break;
    case 1: ep->maxpoints = POS0INT; break;
    case 2: ep->GEV_xi = NUM; break;
    case 3: ep->density_ratio = POS0NUM; break;
    case 4: ep->check_every = POS0INT; break;
    case EXTREME_FLAT: ep->flathull = USRLOG; 
      //     if (ep->flathull != False) 
      // RFERROR1("currently only '%.50s=FALSE' allowed", extreme[EXTREME_FLAT]);
      // Programmierfehler in Huetchen.c
      break;
    case 6: ep->min_n_zhou = POSINT; break;
    case 7: ep->max_n_zhou = POSINT; break;
    case 8: ep->eps_zhou = POSNUM; break;
    case 9: ep->mcmc_zhou = POSINT; break;
    case 10: ep->min_shape_gumbel = NEG0NUM; break;
    case 11: ep->scatter_method =   
	(coord_sys_enum) GetName(el, name, POISSON_SCATTER_NAMES,
				 nPOISSON_SCATTER, POISSON_SCATTER_ANY);
      break; 
   default:  BUG;
    }}
    break;
  case 12 : { // br
    br_param *ep;
    ep = &(options->br);
    switch(j) {
    case 0: ep->BRmaxmem = POSINT; break;
    case 1: ep->BRmeshsize = POSNUM; break;
    case 2: ep->BRvertnumber = POSINT; break;
    case 3: ep->BRoptim = POS0INT; break;
    case 4: ep->BRoptimtol = POS0NUM; break;
    case 5: ep->variobound = NUM; break;
    case 6: ep->deltaAM = POSINT; break;
    default:  BUG;
    }}
    break;
  case 13 : {// distr
    distr_param *ep;
    ep = &(options->distr);
    switch(j) { 
    case 0: ep->safety=POSNUM; break;
    case 1: ep->minsteplen=POS0NUM; break;
    case 2: ep->maxsteps=POSINT; break;
    case 3: ep->parts=POSINT; break;
    case 4: ep->maxit=POSINT; break;
    case 5: ep->innermin=POSNUM; break;
    case 6: ep->outermax=POSNUM; break;
    case 7: ep->mcmc_n=POSNUM; break;
    case 8: ep->repetitions=POSNUM; break;
    default:  BUG;
    }}
    break;
  case 14: { // fit
    fit_param *ef;
    ef = &(options->fit);
    switch(j) {
    case 0: ef->bin_dist_factor = POS0NUM; break; //
    case 1: ef->upperbound_scale_factor = POS0NUM; break; //
    case 2: ef->lowerbound_scale_factor = POS0NUM; break; //
    case 3: ef->lowerbound_scale_LS_factor = POS0NUM; break; //
    case 4: ef->upperbound_var_factor = POS0NUM; break; //
    case 5: ef->lowerbound_var_factor = POS0NUM; break;//
      //    case 6: ef->lowerbound_sill = POS0NUM; break; 
    case 6: ef->scale_max_relative_factor = POS0NUM; break; //
    case 7: ef->minbounddistance = POS0NUM; break;//
    case 8: ef->minboundreldist = POS0NUM; break;// 
    case 9: ef->approximate_functioncalls = POS0INT; break; //
    case FIT_BC_LB: Real(el, name, ef->BC_lambdaLB, 2 * MAXBOXCOXVDIM); break;//
    case FIT_BC_UB: Real(el, name, ef->BC_lambdaUB, 2 * MAXBOXCOXVDIM); break;//
    case 12: ef->use_naturalscaling = LOGI; break;
    case 13: ef->bins = POS0INT; break;//
    case 14: ef->nphi = POS0INT; break;//
    case 15: ef->ntheta = POS0INT; break; //
    case 16: ef->ntime = POS0INT; break;//
    case 17: ef->onlyuser = LOGI; break;
    case 18: { // mle
      int ii=POS0INT; 
      if (ii==0) { ii = 1; warning("shortnamelength set to 1"); }
      if (ii >= LENMSG) {
	ii = LENMSG - 1;  
	WARN1("shortnamelength set to %d", ii); 
      }
      ef->lengthshortname=ii;
    }
    break;
    case 19: ef->split = INT; break;
    case 20: ef->scale_ratio = NUM; break;//
    case 21: ef->critical = INT; break;//
    case 22: ef->n_crit = POS0INT; break;///
    case FIT_MAXNEIGHBOUR: ef->locmaxn = POS0INT; break;//
    case FIT_CLIQUE: {
      int  newval[3],
	len = length(el);      
      if (len < 1 || len > 3) 
	RFERROR1("'%.50s' must have between 1 and 3 components", fit[FIT_CLIQUE]);
      for (int ii=0; ii<len; ii++) newval[ii] = Integer(el, name, ii); 
      if (len == 1) newval[1] = newval[2] = newval[0];
      else if (len == 2) {
	newval[2] = newval[1];
	newval[1] = (int) SQRT((double) newval[1] * newval[0]);
      }

      if (newval[0] > newval[1] || newval[1] > newval[2]) 
   	RFERROR3("%.50s[1] <= %.50s[2] <= %.50s[3] not satisfied",
	     fit[FIT_CLIQUE], fit[FIT_CLIQUE], fit[FIT_CLIQUE]
	     );      
      for (int ii=0; ii<3; ii++) ef->locsplitn[ii] = newval[ii];
      break;
    }
    case 25: ef->locsplitfactor = POS0INT; break;//
    case 26: ef->smalldataset = POS0INT; break;//
    case 27: ef->min_diag = NUM; break;//
    case 28: ef->reoptimise = LOGI; break;
    case 29: ef->optimiser = GetName(el, name, OPTIMISER_NAMES, nOptimiser, 0); //
      break;
    case 30: 
      switch(ef->optimiser) {
      case 3 : ef->algorithm = GetName(el, name, NLOPTR_NAMES, nNLOPTR); //
	break;
      default: ef->algorithm = UNSET;
      }
      break;
    case 31: ef->likelihood =
	GetName(el, name, LIKELIHOOD_NAMES, nLikelihood); ; break; 
    case 32: ef->ratiotest_approx = LOGI; break;
    case 33: ef->split_refined = LOGI; break;
    case 34: ef->cross_refit = LOGI; break;
    case 35: ef->estimate_variance = USRLOG; break;
    case 36: ef->pgtol = POSNUM; break;//
    case 37: ef->pgtol_recall = POSNUM; break;//
    case 38: ef->factr = POSNUM; break;//
    case 39: ef->factr_recall = POSNUM; break;//
    case 40: ef->use_pseudovariogram = LOGI; break;
    default:  BUG;
    }}
    break;
  case 15: { // empvario
    empvario_param *ep;
    ep = &(options->empvario);
    switch(j) {
    case 0: ep->phi0=NUM; break;
    case 1: ep->theta0=NUM; break;    
    case 2: ep->tol=NUM; break;    
    case 3: ep->pseudovariogram = LOGI; break;
    case 4: ep->fft = LOGI; break;
    case 5: ep->halveangles = LOGI; break;
   default: BUG;
    }}
    break;
  case 16: { // gui
    gui_param *gp;
    gp = &(options->gui);
    switch(j) {
    case 0: gp->alwaysSimulate = LOGI; break;
    case 1: 
      gp->method = GetName(el, name, METHOD_NAMES, Forbidden + 1, Nothing);
      break;
    case GUI_SIZE: {
      int sizedummy[2];
      if (length(el) != 2) RFERROR1("length of '%.50s' must be 2", gui[GUI_SIZE]);
      Integer(el, name, sizedummy, 2);
      for (int ii=0; ii<2; ii++) {	
	if (sizedummy[ii] <= 1) 
	  RFERROR("grid in RFgui must contain at least 2 points");
      }
      for (int ii=0; ii<2; ii++) {	 gp->size[ii] = sizedummy[ii]; }
    }
      break;
    default: BUG;
    }}
    break;
 case 17: { // graphics
   //   printf("not ok\n"); break;
    graphics_param *gp = &(options->graphics);
    switch(j) {
    case 0 : gp->close_screen = LOGI; break;
    case 1 : gp->PL = INT; break;
    case 2 : gp->height = NUM; break;
    case GRAPHICS_UPTO : {
      int uptodummy[2];
      if (length(el) != 2)
	RFERROR1("length of '%.50s' must be 2", graphics[GRAPHICS_UPTO]);
      Integer(el, name, uptodummy, 2);
      for (int ii=0; ii<2; ii++) {
	if (uptodummy[ii] <= 0)  RFERROR("increase_upto must be positive");
      } 
      for (int ii=0; ii<2; ii++) { gp->increase_upto[ii] = uptodummy[ii]; }
    }
      break;
    case 4 : gp->always_open = USRLOG;
      break;
    case 5 :
      if (!isList) {
	char old[100];
	strcopyN(old, gp->filename, 100);
	STR(gp->filename, 100);
	if (!gp->onefile && STRCMP(gp->filename, old) !=0) gp->number = 0;
      }
      break;
    case 6 :  if (!isList) {
      bool onefile = gp->onefile;
      gp->onefile = LOGI;
      if (!gp->onefile && onefile) gp->number = 0;
    } break;
    case 7 :  if (!isList) { 
	gp->number = INT; 
      }
      break;
    case 8 : gp->resolution = POSNUM; break;
    case 9 : gp->split_screen = LOGI; break;
    case 10 : gp->width = NUM; break;
    case 11 : gp->always_close = USRLOG; break;
    case 12 : gp->grDefault = LOGI; break;
   default: BUG;
    }}
   break;
       /*
 case 20: { // empvario
    empvario_param *ep;
    ep = &(options->empvario);
    switch(j) {
    case 0: 
      
     default: RFERROR("unknown option (empvario)");
   }}
    break;
    */ 
  case 18 : {
    registers_param *rp = &(options->registers);
    switch(j) {
    case 0: { // simu
      int keynr = INT;
      if ((keynr<0) || (keynr>MODEL_MAX)) RFERROR("register number out of range");
      rp->keynr=keynr; }
      break;
   case 1: { //  predict
      int predict = INT;
      if ((predict<0) || (predict>MODEL_MAX))
	RFERROR("register number out of range");
      rp->predict=predict; }
      break;
    case 2: { //  likelihood
      int likelihood = INT;
      if ((likelihood<0) || (likelihood>MODEL_MAX))
	RFERROR("register number out of range");
      rp->likelihood=likelihood; }
      break;
   default: BUG;
    }}
   break;

  case 19: {
    internal_param *wp = &(options->internal);
    // ACHTUNG internal wird nicht pauschal zurueckgesetzt !
    if (!isList) {
      switch(j) {
      case 0: wp->warn_oldstyle = LOGI;       break;
      case 1: wp->warn_newstyle = LOGI;       break;
      case INTERNALS_NEWANISO: wp->warn_Aniso = LOGI;       break;
      case 3: wp->warn_ambiguous = LOGI;       break;
      case 4: wp->warn_normal_mode = LOGI;       break;
      case 5: wp->warn_mode = LOGI;       break;
      case INTERNALS_STORED_INIT: wp->stored_init = LOGI;       break;
      case 7: wp->warn_scale = LOGI;       break;
      case 8: wp->warn_coordinates = LOGI;	break;
      case INTERNALS_ONGRID: wp->warn_on_grid = LOGI;       break;
      case 10: wp->warn_new_definitions = LOGI;       break;
      case 11: wp->warn_aspect_ratio = LOGI;       break;
      case INTERNALS_COORD_CHANGE: wp->warn_coord_change = LOGI;       break;
      case 13: wp->warn_color_palette = LOGI;       break;
      case INTERNALS_ZENIT: wp->warn_zenit = LOGI;       break;
      case INTERNALS_DO_TESTS: wp->do_tests = LOGI;       break;
      case 16: wp->warn_constant = LOGI;       break;
      case 17: wp->warn_negvar = LOGI;       break;
      case 18: wp->warn_onlyvar = LOGI;       break;
      case 19: wp->warn_mathdef = USRLOG; break;
      case INTERNALS_EX_RED: wp->examples_reduced= POS0INT;  
	if (wp->examples_reduced > 0 && 
	    wp->examples_reduced < MAX_LEN_EXAMPLES)
	  wp->examples_reduced = MAX_LEN_EXAMPLES;
	break;
      case INTERNALS_PARALLEL : wp->warn_parallel = LOGI;       break;
      case 22: wp->allow_duplicated_loc = LOGI; break;
      case 23: wp->warn_seed = LOGI; break;
      case 24: wp->declare_PL = POS0INT; break; 
      case 25: wp->warn_modus_operandi = LOGI; break;
     default: BUG; 
      }
    } else {
      switch(j) {
      case INTERNALS_STORED_INIT: wp->stored_init = LOGI;       break;
      case INTERNALS_ONGRID :  wp->warn_on_grid = LOGI; break;
      case INTERNALS_COORD_CHANGE: wp->warn_coord_change = LOGI;       break;
      case INTERNALS_ZENIT: wp->warn_zenit = LOGI;       break;
      case INTERNALS_DO_TESTS: wp->do_tests = LOGI;       break;
      case INTERNALS_EX_RED: wp->examples_reduced= POS0INT;  
	if (wp->examples_reduced > 0 && 
	    wp->examples_reduced < MAX_LEN_EXAMPLES)
	  wp->examples_reduced = MAX_LEN_EXAMPLES;
	break;
      case INTERNALS_PARALLEL : wp->warn_parallel = LOGI;       break;
  default : {} // none
      }      
    }
  }
    break;

  case 20: {
    coords_param *cp = &(options->coords);
    switch(j) {
    case COORDS_XYZNOTATION: cp->xyz_notation = USRLOG; break;
    case 1: {
      coord_sys_enum coord =
	(coord_sys_enum) GetName(el, name, COORD_SYS_NAMES, nr_coord_sys,
				 coord_auto);
      if (coord == coord_auto || coord == earth || coord == sphere || 
	  coord == cartesian)
	cp->coord_system = coord;
      else RFERROR2("Cannot take '%.50s' as a value for the initial '%.50s'", 
		CHAR(STRING_ELT(el, 0)), coords[j]);
      if ((coord == earth && cp->anglemode != degree) ||
	  (coord != earth && cp->anglemode == degree)
	  ) {
	cp->anglemode = coord == earth ? degree : radians;
	if (PL > 0 && !isList) {
	  PRINTF("Angle mode switches to '%.50s'.\n", ANGLE_NAMES[cp->anglemode]);
	}
      }
    } 
     break; 
    case 2: getUnits(el, name, cp->newunits, NULL);
      break;
    case 3: getUnits(el, name, cp->curunits, isList ? NULL : cp->newunits);
      break;
    case 4: getUnits(el, name, cp->varunits, NULL);      
      break;
    case COORDS_DATAIDX:
      // either exct column number must be given
      //                      or a set of potential variable names
      Integer2(el, name, cp->data_idx);
      if (cp->data_idx[1] != NA_INTEGER) {
	int n = cp->data_idx[1] - cp->data_idx[0] + 1;
	if (n > MAXDATANAMES)
	  RFERROR1("maximum number named variables is %d", MAXDATANAMES);
	cp->data_nr_names = 0;
      }
      break;
    case COORDS_DATANAMES:
      cp->data_nr_names = length(el);
      if (cp->data_nr_names > 0)
	String(el, name, cp->data_names, cp->data_nr_names);
      break;
    case COORDS_XIDX:
      Integer2(el, name, cp->x_idx);
      if (cp->x_idx[1] != NA_INTEGER) {
	int n =  cp->x_idx[1] - cp->x_idx[0] + 1;
	if (n > MAXDATANAMES)
	  RFERROR1("maximum number named variables is %d", MAXCOOORDNAMES);
	cp->x_nr_names = 0;
      }
      break;
    case COORDS_XNAMES:
      cp->x_nr_names = length(el);
      if (cp->x_nr_names > 0) String(el, name, cp->x_names, cp->x_nr_names);
      break;
    case 9: 
      cp->new_coord_system =
	(coord_sys_enum) GetName(el, name, COORD_SYS_NAMES, nr_coord_sys,
				 coord_keep);
      if (cp->new_coord_system == coord_auto) cp->new_coord_system = coord_keep;
      break; 

    case ZENIT: {
      Real(el, name, cp->zenit, 2);
    }
    break; 
 
    case 11: cp->polar_coord = LOGI; 
      break;

    case 12: {
      cp->anglemode = (angle_modes) GetName(el, name, ANGLE_NAMES,
					    last_angle_mode + 1, -1);
    }
     break;
     
    case 13: cp->allow_earth2cart = LOGI; 
      break;
  
   default: BUG;
    }}
    break;

  case 21: {
    special_param *sp = &(options->special);
    switch(j) {
    case 0: sp->multcopies = POSINT;
	break;      
    default: BUG;
    }}
    break;


  case 22: { // obsolete
    internal_param *wp = &(options->internal);
    switch(j) {
    case 0: wp->warn_oldstyle = LOGI;       break;
    case 1: wp->warn_newstyle = LOGI;       break;
    case 2: wp->warn_Aniso = LOGI;       break;
    case 3: wp->warn_ambiguous = LOGI;       break;
    case 4: wp->warn_normal_mode = LOGI;       break;
    case 8: RFERROR("'matrix_inversion' has been changed to 'matrix_methods'" );
      break;
    case 9: RFERROR("'matrix_tolerance' has been changed to 'svd_tol'" );
      break;

    default:
      RFERROR1("RFoption '%.50s' has been removed.", obsolete[j]);
    }}
    break;
    

  default: BUG;
  }
  // printf("end setparameter\n");
}

  
void getRFoptions(SEXP sublist, int i, int local) {
  int k=0;
  char x[2]=" ";
  //#define ADD(ELT) {printf(#ELT"\n");SET_VECTOR_ELT(sublist, k++, ELT);}
  globalparam *options = WhichOptionList(local);
   switch(i) {
  case 0 : { 
    // printf("OK %d\n", i);
     general_param *p = &(options->general);
    ADD(ScalarString(mkChar(MODE_NAMES[p->mode])));
    ADD(ScalarLogical(p->storing));
    ADD(ScalarInteger(p->every));    
    ADD(ScalarReal(p->gridtolerance));
    ADDCHAR(p->pch);
    if (p->naturalscaling==0 || p->naturalscaling==1) {
      ADD(ScalarLogical(p->naturalscaling));
    } else {
      ADD(ScalarInteger(p->naturalscaling));
    }
    ADD(ScalarLogical(p->sp_conform));
    ADD(ExtendedBooleanUsr(p->exactness));    
    ADD(ScalarLogical(p->allowdist0));   
    ADD(ScalarLogical(p->na_rm_lines));   
    ADD(ScalarLogical(p->vdim_close_together));
    ADD(ScalarInteger(p->expected_number_simu));    
    ADD(ScalarLogical(p->detailed_output));   
    ADD(ScalarLogical(p->Ttriple == NA_INTEGER ? NA_LOGICAL
		      : p->Ttriple != 0));
    ADD(ScalarLogical(p->returncall));   
    ADD(ScalarString(mkChar(OUTPUTMODE_NAMES[p->output])));
    ADD(ScalarString(mkChar(REPORTCOORD_NAMES[p->reportcoord])));
    ADD(ScalarInteger(p->set + 1));   
    ADD(ScalarInteger(p->seed_incr));   
  }
  break;
  
  //  printf("OK %d\n", i);
  
  case 1 : {
    gauss_param *p = &(options->gauss);
    // nachfolgend sollte immer >= 0 sein
    ADD(ScalarLogical(p->paired));   
    ADD(ExtendedBooleanUsr(p->stationary_only));
    ADD(ScalarReal(p->approx_zero));    
    ADD(ScalarInteger(p->direct_bestvariables));    
    ADD(ScalarLogical(p->loggauss));   
    SET_VECTOR_ELT(sublist, k++, Num(p->boxcox, 2 * MAXBOXCOXVDIM, 
					2 * MAXBOXCOXVDIM)); 
  }
   break;

  case 2 : {
    krige_param *p = &(options->krige);
    ADD(ScalarLogical(p->ret_variance));   
    ADD(ScalarInteger(p->locmaxn));
    SET_VECTOR_ELT(sublist, k++, Int(p->locsplitn, 3, 3));
    ADD(ScalarInteger(p->locsplitfactor));
    ADD(ScalarLogical(p->fillall));   
  }
   break;
   
  case 3 : {
    ce_param *p = &(options->ce);
    ADD(ScalarLogical(p->force)); 
    SET_VECTOR_ELT(sublist, k++, Num(p->mmin, MAXCEDIM, MAXCEDIM)); 
    ADD(ScalarInteger(p->strategy)); 
    ADD(ScalarReal(p->maxGB)); 
    ADD(ScalarReal(p->maxmem)); 
    ADD(ScalarReal(p->tol_im)); 
    ADD(ScalarReal(p->tol_re));	      
    ADD(ScalarInteger(p->trials));    
    ADD(ScalarLogical(p->useprimes)); 
    ADD(ScalarLogical(p->dependent)); 
    ADD(ScalarReal(p->approx_grid_step)); 
    ADD(ScalarInteger(p->maxgridsize));    
  } 
   break;
   
  case 4: {
    direct_param *p = &(options->direct);
    //ADD(ScalarInteger(p->inversionmethod));    
    //ADD(ScalarReal(p->svdtolerance));     
    ADD(ScalarInteger(p->maxvariables));    
  }
    break;

  case 5:  {
    nugget_param *p = &(options->nugget);
    ADD(ScalarReal(p->tol));
    // ADD(ScalarLogical(p->meth));
  }
    break;
    
  case 6 : {
    sequ_param *p = &(options->sequ);
    ADD(ScalarInteger(p->back));    
    ADD(ScalarInteger(p->initial));    
  }
   break;
   
  case 7 :; {
    spectral_param *p = &(options->spectral);
    SET_VECTOR_ELT(sublist, k++, Int(p->lines, MAXTBMSPDIM,MAXTBMSPDIM));
    ADD(ScalarLogical(p->grid));
    ADD(ScalarReal(p->prop_factor));
    ADD(ScalarReal(p->sigma));
  }
   break;

  case 8: {
    tbm_param *p = &(options->tbm);
    // nachfolgend sollte immer >= 0 sein
    ADD(ScalarInteger(p->tbmdim));
    ADD(ScalarInteger(p->fulldim));
    SET_VECTOR_ELT(sublist, k++, Num(p->center, MAXTBMSPDIM, MAXTBMSPDIM));
    ADD(ScalarInteger(p->points));     
    // ADD(p->method>=0 ? ScalarString(mkChar(METHOD_NAMES[p->method])) : 
    //	R_NilValue);
     SET_VECTOR_ELT(sublist, k++, Int(p->lines, MAXTBMDIM,MAXTBMDIM));
    ADD(ScalarReal(p->linesimufactor));
    ADD(ScalarReal(p->linesimustep));
    ADD(ExtendedBooleanUsr(p->layers));
    ADD(ScalarLogical(p->grid));
  }
   break;


  case 9: {
    mpp_param *p = &(options->mpp);
     ADD(ScalarInteger(p->n_estim_E));
    SET_VECTOR_ELT(sublist, k++,
		   Num(p->intensity, MAXMPPDIM ,MAXMPPDIM));
    //ADD(ScalarReal(p->refradius_factor));
    ADD(ScalarReal(p->about_zero));
    ADD(ScalarReal(p->shape_power));
    SET_VECTOR_ELT(sublist, k++, Num(p->scatter_step, MAXMPPDIM,MAXMPPDIM)); 
    SET_VECTOR_ELT(sublist, k++, Int(p->scatter_max, MAXMPPDIM,MAXMPPDIM)); 
     //    SET_VECTOR_ELT(sublist, k++, Num(p->plus, MAXMPPDIM ,MAXMPPDIM));
    //    ADD(ScalarReal(p->approxzero));
    //   ADD(ScalarReal(p->samplingdist));
    //  ADD(ScalarReal(p->samplingr));
     //  ADD(ScalarReal(p->p));
  }
   break;

  case 10: {
    hyper_param *p = &(options->hyper);
    ADD(ScalarInteger(p->superpos));    
    ADD(ScalarInteger(p->maxlines));    
    ADD(ScalarInteger(p->mar_distr));
    ADD(ScalarReal(p->mar_param));    
  }
   break;

  case 11 : {
    extremes_param *p = &(options->extreme);
    ADD(ScalarReal(p->standardmax));
    ADD(ScalarInteger(p->maxpoints));           
    ADD(ScalarReal(p->GEV_xi));
    ADD(ScalarReal(p->density_ratio));
    ADD(ScalarInteger(p->check_every));
    ADD(ExtendedBooleanUsr(p->flathull));
    ADD(ScalarInteger(p->min_n_zhou));
    ADD(ScalarInteger(p->max_n_zhou));
    ADD(ScalarReal(p->eps_zhou));
    ADD(ScalarInteger(p->mcmc_zhou));
    ADD(ScalarReal(p->min_shape_gumbel));
    ADD(ScalarString(mkChar(POISSON_SCATTER_NAMES[p->scatter_method])));
 
  }
   break;

  case 12: {
    br_param *p = &(options->br);
    ADD(ScalarInteger(p->BRmaxmem)); 
    ADD(ScalarReal(p->BRmeshsize));
    ADD(ScalarInteger(p->BRvertnumber));
    ADD(ScalarInteger(p->BRoptim));
    ADD(ScalarReal(p->BRoptimtol));
    ADD(ScalarReal(p->variobound));
    ADD(ScalarInteger(p->deltaAM));
  }
   break;

  case 13 : {
    distr_param *p = &(options->distr);
    ADD(ScalarReal(p->safety));
    ADD(ScalarReal(p->minsteplen));
    ADD(ScalarInteger(p->maxsteps));
    ADD(ScalarInteger(p->parts));
    ADD(ScalarInteger(p->maxit));
    ADD(ScalarReal(p->innermin));
    ADD(ScalarReal(p->outermax));
    ADD(ScalarInteger(p->mcmc_n));
    ADD(ScalarInteger(p->repetitions));
  }
   break;
   
  case 14 : {
    fit_param *p = &(options->fit);
    ADD(ScalarReal(p->bin_dist_factor));
    ADD(ScalarReal(p->upperbound_scale_factor)); 
    ADD(ScalarReal(p->lowerbound_scale_factor)); 
    ADD(ScalarReal(p->lowerbound_scale_LS_factor)); 
    ADD(ScalarReal(p->upperbound_var_factor)); 
    ADD(ScalarReal(p->lowerbound_var_factor));
    ADD(ScalarReal(p->scale_max_relative_factor)); 
    ADD(ScalarReal(p->minbounddistance));
    ADD(ScalarReal(p->minboundreldist));  
    ADD(ScalarInteger(p->approximate_functioncalls)); 
    SET_VECTOR_ELT(sublist, k++, Num(p->BC_lambdaLB, 2 * MAXBOXCOXVDIM, 
					2 * MAXBOXCOXVDIM));
    SET_VECTOR_ELT(sublist, k++, Num(p->BC_lambdaUB, 2 * MAXBOXCOXVDIM, 
					2 * MAXBOXCOXVDIM));
    ADD(ScalarLogical(p->use_naturalscaling));
    ADD(ScalarInteger(p->bins));
    ADD(ScalarInteger(p->nphi));
    ADD(ScalarInteger(p->ntheta)); 
    ADD(ScalarInteger(p->ntime));
    //   ADD(ScalarReal(p->sill));
    ADD(ScalarLogical(p->onlyuser));
    //    ADD(ScalarString(mkChar(OPTIM_VAR_NAMES[p->optim_var_estim])));
    ADD(ScalarInteger(p->lengthshortname));
    ADD(ScalarInteger(p->split));
    ADD(ScalarReal(p->scale_ratio));
    ADD(ScalarInteger(p->critical));
    ADD(ScalarInteger(p->n_crit));
    ADD(ScalarInteger(p->locmaxn));
    SET_VECTOR_ELT(sublist, k++, Int(p->locsplitn, 3, 3));
    ADD(ScalarInteger(p->locsplitfactor));
    ADD(ScalarInteger(p->smalldataset));
    ADD(ScalarReal(p->min_diag));
    ADD(ScalarLogical(p->reoptimise));
    ADD(p->optimiser>=0 ? ScalarString(mkChar(OPTIMISER_NAMES[p->optimiser])) 
	: R_NilValue);
    ADD(p->algorithm == UNSET ? R_NilValue :
	ScalarString(mkChar(p->optimiser == 3 ? NLOPTR_NAMES[p->algorithm]
			    : "")
		     ));
    ADD(ScalarString(mkChar(LIKELIHOOD_NAMES[p->likelihood])));
    ADD(ScalarLogical(p->ratiotest_approx));
    ADD(ScalarLogical(p->split_refined));
    ADD(ScalarLogical(p->cross_refit));
    ADD(ExtendedBooleanUsr(p->estimate_variance));
    ADD(ScalarReal(p->pgtol)); 
    ADD(ScalarReal(p->pgtol_recall)); 
    ADD(ScalarReal(p->factr)); 
    ADD(ScalarReal(p->factr_recall)); 
    ADD(ScalarLogical(p->use_pseudovariogram));
  } 
   break;
   
  case 15: {
    empvario_param *p = &(options->empvario);   
    ADD(ScalarReal(p->phi0));
    ADD(ScalarReal(p->theta0));
    ADD(ScalarReal(p->tol));
    ADD(ScalarLogical(p->pseudovariogram));
    ADD(ScalarLogical(p->fft));
    ADD(ScalarLogical(p->halveangles));
  }
   break;

  case 16: {
    gui_param *p = &(options->gui);   
    ADD(ScalarLogical(p->alwaysSimulate));
    ADD(p->method>=0 ? ScalarString(mkChar(METHOD_NAMES[p->method])) 
	: R_NilValue);
    SET_VECTOR_ELT(sublist, k++, Int(p->size, 2, 2));
  }
   break;

  case 17 : {
    graphics_param *p = &(options->graphics);   
    ADD(ScalarLogical(p->close_screen));
    ADD(ScalarInteger(p->PL));
    ADD(ScalarReal(p->height));
    SET_VECTOR_ELT(sublist, k++, Int(p->increase_upto, 2, 2));
    ADD(ExtendedBooleanUsr(p->always_open));
    ADD(ScalarString(mkChar(p->filename)));
    ADD(ScalarLogical(p->onefile));
    ADD(ScalarInteger(p->number));
    ADD(ScalarReal(p->resolution));
    ADD(ScalarLogical(p->split_screen));
    ADD(ScalarReal(p->width));
    ADD(ExtendedBooleanUsr(p->always_close));
    ADD(ScalarLogical(p->grDefault));
   }
   break;

  case 18: {
    registers_param *p = &(options->registers);
    ADD(ScalarInteger(p->keynr));    
    ADD(ScalarInteger(p->predict));    
    ADD(ScalarInteger(p->likelihood));    
  }
   break;
  
  case 19 : {
    internal_param *p = &(options->internal);
    ADD(ScalarLogical(p->warn_oldstyle));
    ADD(ScalarLogical(p->warn_newstyle));
    ADD(ScalarLogical(p->warn_Aniso));
    ADD(ScalarLogical(p->warn_ambiguous));
    ADD(ScalarLogical(p->warn_normal_mode));
    ADD(ScalarLogical(p->warn_mode));
    ADD(ScalarLogical(p->stored_init));
    ADD(ScalarLogical(p->warn_scale));
    ADD(ScalarLogical(p->warn_coordinates));
    ADD(ScalarLogical(p->warn_on_grid));
    ADD(ScalarLogical(p->warn_new_definitions));
    ADD(ScalarLogical(p->warn_aspect_ratio));
    ADD(ScalarLogical(p->warn_coord_change));
    ADD(ScalarLogical(p->warn_color_palette));
    ADD(ScalarLogical(p->warn_zenit));
    ADD(ScalarLogical(p->do_tests));
    ADD(ScalarLogical(p->warn_constant));
    ADD(ScalarLogical(p->warn_negvar));
    ADD(ScalarLogical(p->warn_onlyvar));
    ADD(ExtendedBooleanUsr(p->warn_mathdef));
    ADD(ScalarInteger(p->examples_reduced));
    ADD(ScalarLogical(p->warn_parallel));
    ADD(ScalarLogical(p->allow_duplicated_loc));
    ADD(ScalarLogical(p->warn_seed));
    ADD(ScalarInteger(p->declare_PL));
    ADD(ScalarLogical(p->warn_modus_operandi));
   }
   break;


  case 20 : {
    coords_param *p = &(options->coords);
    ADD(ExtendedBooleanUsr(p->xyz_notation));    
    ADD(ScalarString(mkChar(COORD_SYS_NAMES[p->coord_system])));
    ADD(UNITS(p->newunits));
    ADD(UNITS(p->curunits));
    ADD(UNITS(p->varunits)); // 4
    int idx_default[2] = {1, NA_INTEGER};
    if (p->data_nr_names == 0) {
      SET_VECTOR_ELT(sublist, k++, Int(p->data_idx, 2, 2));
      SET_VECTOR_ELT(sublist, k++, String(NULL, 0, 0));
    } else {      
      SET_VECTOR_ELT(sublist, k++, Int(idx_default, 2, 2));
      SET_VECTOR_ELT(sublist, k++,
		     String(p->data_names, p->data_nr_names, p->data_nr_names));
    }
    if (p->x_nr_names == 0) { // 7
      SET_VECTOR_ELT(sublist, k++, Int(p->x_idx, 2, 2));
      SET_VECTOR_ELT(sublist, k++, String(NULL, 0, 0));
    } else {
      SET_VECTOR_ELT(sublist, k++, Int(idx_default, 2, 2));
      SET_VECTOR_ELT(sublist, k++,
		     String(p->x_names, p->x_nr_names, p->x_nr_names));
    }
    ADD(ScalarString(mkChar(COORD_SYS_NAMES[p->new_coord_system]))); //9
    SET_VECTOR_ELT(sublist, k++, Num(p->zenit, 2, 2));     
    ADD(ScalarLogical(p->polar_coord));
    ADD(ScalarString(mkChar(ANGLE_NAMES[p->anglemode])));
    ADD(ScalarLogical(p->allow_earth2cart));
 }
   break;

  case 21 : {
    special_param *p = &(options->special);
    ADD(ScalarInteger(p->multcopies));    
  }
   break;
   
  default : BUG;
  }
 

}


void finalparameter(int VARIABLE_IS_NOT_USED local) {
  PL = GLOBAL_UTILS->basic.Cprintlevel - PLoffset;
  CORES = GLOBAL_UTILS->basic.cores;
}

void attachRFoptionsRandomFields(int *show) { // no print commands!!!
   includeXport();
   // printf("GLO RU start %l/d\n", (long  int) GLOBAL_UTILS);
  Ext_getUtilsParam(&GLOBAL_UTILS);
  // printf("GLO RU %ld\n", (long int) GLOBAL_UTILS);
  //  GLOBAL_UTILS->basic.cores = 2;
  //  PRINTF("Some default options in pkg 'RandomFieldsUtils' for matrix decomposition have been changed by 'RandomFields'. Only if you use 'RandomFieldsUtils' directly, you might notice these changes."); 
  GLOBAL_UTILS->solve.max_chol = 16384;
  GLOBAL_UTILS->solve.max_svd = 6555;
  GLOBAL_UTILS->solve.pivot = PIVOT_AUTO;
  GLOBAL_UTILS->solve.pivot_check = Nan;
  Ext_attachRFoptions(prefixlist, prefixN, all, allN,
		      setparameter, finalparameter, getRFoptions, NULL,
		      PLoffset, true);
   finalparameter(isGLOBAL);
   InitModelList();
   if (*show) {
     PRINTF("You get the old former graphic style by 'RFoption(grDefault=FALSE)'"); // OK
#ifdef DO_PARALLEL
   PRINTF("'RandomFields' will use OMP\n");
#else
   PRINTF("'RandomFields' will NOT use OMP\n");
#endif
   }
#ifdef SCHLATHERS_MACHINE   
   //   PRINTF("Martin todo: on.exit(RFoption) streichen; LOCAL=isGLOBAL, COPY_GLOBAL+LOCAL; nur LOCAL\n");
#endif   
  //  printf("end  ssss\n");
}


void detachRFoptionsRandomFields() {
  Ext_detachRFoptions(prefixlist, prefixN);
}

void RelaxUnknownRFoptions(int *RELAX) {
  Ext_relaxUnknownRFoption((bool) *RELAX); 
}

