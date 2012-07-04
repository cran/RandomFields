#ifndef RFsimu_H
#define RFsimu_H 1

#include "auxiliary.h"
#include "RandomFields.h"
#include "includeMarkov.h"
#include "error.h"
#include <string.h>
#include <R_ext/Complex.h>
#include <errno.h>


///////////////////////////////////////////////////////////////////////
// max. dimension of the field. Standard is 4
// in some application higher dimensions are needed.
///////////////////////////////////////////////////////////////////////
/// !!!! *** Reihenfolge unbedingt beibehalten *** !!!
#define MAXCOVDIM 11000
#define MAXMLEDIM MAXCOVDIM
#define MAXSIMUDIM MAXCOVDIM

#define MAXCEDIM 13
#define MAXTBMSPDIM 4
#define MAXMPPDIM 4
#define MAXHYPERDIM 4
#define MAXNUGGDIM 20
#define MAXVARIODIM 20
// !!! if list get longer change GetrfParameters in userinterface.cc
// !!! und RFparameters

#define MAXMAKEEXPLICITE 20

///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////
// type of covariance functions that need distinguished treatment
// Reihenfolge beachten! Hoehere Nummern sind in kleinere umwandelbar

// IMPORTANT ISO and STAT must start with 0 and may not have gaps
// see GetModelNames in getNset.R !

// NOTE!!!
// if any definition is changed check GetModelNames in R !!!
#define ISOTROPIC (isotropy_type) 0
#define SPACEISOTROPIC (isotropy_type) 1 /* fully symmertric */
#define ZEROSPACEISO (isotropy_type) 2   /* isotropic if time diff is zero */
#define ANISOTROPIC (isotropy_type) 3 
#define PREVMODELI (isotropy_type) 4  /* type taken from submodels */
#define ISO_MISMATCH  (isotropy_type) 9 /* always last ! */
#define LAST_ISO 9

// *****************
// NOTE!!!
// if any definition is changed check setbackward() AND GetModelNames in R !!!
// IN - variables !
#define STATIONARY (stationary_type) 0
#define VARIOGRAM (stationary_type) 1
#define IRF (stationary_type) 2  /* but not intrinsically stationary */
#define AUXMATRIX (stationary_type) 3 /* not a genuine covariance function;
			takes only x vector, but is not statonary */
#define AUXVECTOR (stationary_type) 4 /* not a genuine covariance function;
			takes only x vector, but is not statonary */
#define GENERALISEDCOVARIANCE (stationary_type) 5 /* "stationary",
    but "cov fct" does not have finite value at the origin, e.g. 1/x */
#define COVARIANCE (stationary_type) 6 /* x,y; FIRST non-stationary;
					this fact is used allover !! */
#define GEN_VARIOGRAM (stationary_type) 7 /* x,y but not stationary */
#define PREVMODELS (stationary_type) 8 /* type taken from calling models */
#define STAT_MISMATCH (stationary_type) 9 /* always last ! */
#define LAST_STAT 9


// vdim is "strictly" passed from bottom to top
#define UNIVARIATE 1
/* 0 sicherheitshalber nirgendswo nehmen */
#define PARAMETERM -1 /* parameter dependent; for CovList only */
#define SUBMODELM -2
//#define PARAMSUBMODELM -2 /* for CovList only */
#define M_MISMATCH -3

// way of implementing simulation methods:
#define NOT_IMPLEMENTED 0 /* do not change this value except to false */
#define IMPLEMENTED 1     /* do not change this value except to true */
#define NUM_APPROX 2 
#define GIVEN_METH_IGNORED 3
#define HYPERIMPLEMENTED 4



///////////////////////////////////////////////////////////////////////
// BASIC DIMENSIONS AND VARIABLES
///////////////////////////////////////////////////////////////////////
#define DISTRMAXCHAR 10
#define DISTRNR 3
#define MAXCHAR 16 /* max number of characters for covariance name  */
#define PARAMMAXCHAR MAXCHAR
#define METHODMAXCHAR MAXCHAR /* max number of character to describe a Method (including \0)*/

#define MAXSUB 10     /* max number of submodels -- NEVER MORE THAN 255 */
#define MAXTBMDIM 3
#/* if explicite simulation method is not given by ????????
   user, the programme will first try a methods 
   according to "first"; 
   then it tries the other methods according to
   a default list: (CircEmbed only if the simulation is on a grid)
   1 dim: CircEmbed, Direct
   2 dim: CircEmbed, TBM2, spectralTBM, TBM3, 
   AdditiveMpp, Direct
   3 dim: CircEmbed, TBM3, Direct
   the method actually used is stored in actually_used. 
*/ 
#define MAXKEYS 10
#define MAXNRCOVFCTS 150
#define MAXPARAM 6
#define MAXELEMENTS 10
#define MAX_NA 30


#define LISTOF 100  /* may not interfere with define of *SXP in Rinternal.h*/


// #define MAXINTERNALPARAM 2
extern double ZERO[MAXSIMUDIM];
extern int True;
extern int False;

typedef char name_type[][MAXCHAR];
typedef char NAname_type[MAX_NA][255];

///////////////////////////////////////////////////////////////////////
// PARAMETER AND THEIR LOCATIONS
///////////////////////////////////////////////////////////////////////


// parameters used in local cov fcts
#define LOCAL_R 0 
#define LOCAL_MSG LOCAL_R + 1
#define CUTOFF_B LOCAL_MSG + 1
#define CUTOFF_ASQRTR CUTOFF_B + 1
#define CUTOFF_THEOR CUTOFF_ASQRTR + 1 
#define CUTOFF_MAX CUTOFF_THEOR + 1   /* size of vector q */
#define INTRINSIC_A0 CUTOFF_MAX + 1
#define INTRINSIC_A2 INTRINSIC_A0 + 1
#define INTRINSIC_B INTRINSIC_A2 + 1
#define INTRINSIC_MAX INTRINSIC_B + 1 /* size of vector q */
#define LOCAL_MAX INTRINSIC_MAX

#define pLOC_DIAM 0 // parameter p
#define pLOC_A 1

#define DIAG1 0

#define DVAR 0
#define DSCALE 1
#define DANISO 2
#define DALEFT 3
#define DPROJ 4
#define DMAX DPROJ

#define BETAMIXED 1

///////////////////////////////////////////////////////////////////////
// BASICS
///////////////////////////////////////////////////////////////////////

typedef enum SimulationType {
  /* never change the ordering -- at least check with "compatible" and 
     InternalGetKeyInfo, 
     variables: METHODNAMES, Besselupperbound, selectlocal, do_comp_meth,
                initmethod, do_incomp_meth, Standard, LastDefault,
		DirectFst, DirectLst, raw, 
		Allowed, allowed
 */
  CircEmbed,     /* Circulant embedding -- must always be the first one! */
  CircEmbedCutoff,
  CircEmbedIntrinsic,
  TBM2, TBM3,    /* Turning Bands, performed in 2 or 3 dimensions */
  SpectralTBM,   /* Spectral turning bands, only 2 dimensional */
  Direct,        /* directly, by matrix inversion */
  Sequential,    /* sequential simulation */
  Markov,        /* Gaussian Markov Field */
  Average,       /* Random spatial averages */
  Nugget,        /* just simulate independent variables */
  RandomCoin,   /* "additive MPP (random coins)", will be only in the plane*/
  Hyperplane,    /* by superposition of Poisson hyperplane tessellations;
		    only in the plane! Not implemented yet*/
  Nothing,       /* must always be the last of the list of methods for
		    simulation Gaussian random fields;
		    used in preference list whether variogram/covariance fct can be
		    calculated
		 */
  MaxMpp,
  ExtremalGauss,
  Forbidden      /* must always be the last one */
} SimulationType;
/*
  NOTE: change: Besserlupperbound in CovFct.cc,
  empty_idx and header in getNset.cc 
  METHODNAMES, InitModelList and ErrorMessage in initNerror.cc
  addXXX and createmodel in getNset
  init_method_list in simu.cc
*/
typedef SimulationType LocalType;
#define SimulationTypeLength ((int) Forbidden) - 1


///////////////////////////////////////////////////////////////////////
// USER CHANGABLE PARAMETERS
///////////////////////////////////////////////////////////////////////
typedef struct ce_param {
  bool force, useprimes, dependent;
  char strategy;
  int trials, method;
  double maxmem, tol_re, tol_im, mmin[MAXCEDIM];
} ce_param;

typedef struct tbm_param {
  int lines,           // number of lines simulated
      layers;
  double linesimufactor, /*factor by which simulation on the line is finer 
			     than the grid */
      linesimustep;      /* grid lag defined absolutely */
//  bool  tbm2num; 
} tbm_param;

typedef struct tbmgen_param {
  SimulationType method;
  double center[MAXTBMSPDIM];
  int points;
} tbmgen_param;


typedef struct spectral_param {
  bool grid, metro, ergodic;
  int nmetro;
  double sigma;
  int lines[MAXTBMSPDIM];
} spectral_param;

typedef enum InversionMethod { Cholesky, SVD, TCholesky, TSVD,
			       NoFurtherInversionMethod } 
InversionMethod;
typedef struct direct_param {
  InversionMethod inversionmethod;
  int bestvariables, maxvariables;
  double svdtolerance;
} direct_param;

typedef struct sequ_param{
  int max, back, initial;
} sequ_param;

typedef struct markov_param{
  int neighbours, cyclic, maxmem;
  double precision;
} markov_param;

typedef struct ave_param {
} ave_param;

typedef struct nugget_param {
  double tol;
  bool meth;
} nugget_param;

typedef struct mpp_param{
  int locations;
  double intens[MAXMPPDIM], //  add_realisations;
    plus[MAXMPPDIM],  // mpp.radius
    relplus[MAXMPPDIM],
//    scale[MAXMPPDIM],
    approxzero,
    samplingdist,
    samplingr,
    p, beta;
} mpp_param;	

#define HYPER_UNIFORM 0   
typedef struct hyper_param {
  int superpos, maxlines, max_distr;
  double mar_param;
} hyper_param;

typedef struct special_param {
//  bool none;
} special_param;

typedef struct extremes_param {
  double standardmax;
} extremes_param;

// others
typedef struct general_param {
  int skipchecks;
  char pch;/*  character printed after each simulation
    just for entertainment of the user
    except for "!", then the numbers are shown
	   */
  bool storing;
  /* true: intermediate results are stored: might be rather memory consuming,
         but the simulation can (depending on the method chosen) be much faster
	 when done the second time with exactly the same parameters
	 do not forget to call DeleteAllKeys when all the simulation are done.
   false: intermediate results are not stored if SimulateRF is called or
         stored only until DoSimulateRF is called.
	 minimum memory consumption, but might be rather slow if many
	 simulations are performed with exactly the same parameters
	 to be safe call DeleteAllKeys when all the simulations are done
	 Warning! init_foo may depend on GENERAL_STORING! This may cause problems
	 if GENERAL_STORING is changed between init_foo and do_foo from false to
	 true and if do_foo is called more than once! (SimulateRF is safe in
	 this case)

	 indifferent behaviour for the simulation methods if parameters are 
	 changed after each simulation.
	 In case, only MAXKEYS different parameter sets are used, but after each 
	 step the parameter set is changed, use different keynr for each 
	 parametere set, and STORING==true, to get fast simulation 
  */
  int printlevel; /* 0:no messages; 
   1:error messages; 
   2: hints when algorithm is at a forcation; 
   >=3: more and more debugging information
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

  int every;          // every `every' line is announced if every>0

    // bool aniso;  // currently cannot be changed by the user !!

  // siehe InternalCov.cc line 350, simu.cc line 394
  /* expand scale immediately to aniso;
		 this allows higher flexibility in the choice of 
		 simulation methods (tbm2 does not like dimension reduction),
		 but it is slower
	      */
  
} general_param;

typedef struct decision_param{
  int stationary_only, exactness;
} decision_param;

typedef struct globalparam{
  general_param general;
  decision_param decision;
  ce_param ce, localce, tbmce;
  tbm_param TBM2, TBM3;
  tbmgen_param TBMgen;
  spectral_param spectral;
  direct_param direct;
  sequ_param sequ;
  markov_param markov;
  ave_param ave;
  nugget_param nugget;
  mpp_param mpp;
  hyper_param hyper;
  special_param special;
  extremes_param extremes;
} globalparam;
extern globalparam GLOBAL;
extern bool NAOK_RANGE;



///////////////////////////////////////////////////////////////////////
// GENERAL PARAMETERS FOR THE SIMULATIONMETHODS
///////////////////////////////////////////////////////////////////////
#define SET_DESTRUCT(A)\
  assert(meth->S==NULL && meth->destruct==NULL);\
  meth->destruct = A;

typedef struct method_type method_type;
typedef int (*init_meth)(method_type *meth);
typedef void (*do_meth)(method_type *meth, res_type *res);
typedef void (*destructor)(void **);
typedef struct cov_model cov_model;
typedef double *coord_type[MAXSIMUDIM];
typedef double static_grid[MAXSIMUDIM][3];
typedef enum matrix_type {TypeIso, TypeDiag, TypeTimesep, TypeAny} matrix_type;

typedef struct trend_type {
  int lTrendFct, // only used for storing -- used in R
    lLinTrend, // only used for storing -- used in R
    TrendModus;                 /* spatial dimension of RF */
  char *TrendFunction; // only used for storing -- used in R
  double mean,         // only used if storagemodus=0 (see rf.R)
    *LinearTrend; // only used for storing -- used in R
} trend_type;

typedef struct location_type {
  int 
    timespacedim,      /* total dimension, including time and space */
    length[MAXSIMUDIM], /* grid==true: number of pixels of each direction 
		       grid==false: length[0]=number of spatial points in total
		     */
      lx,
      spatialdim;
  long spatialtotalpoints,
    totalpoints;      /* number of pixels in total;
			    == lx if grid==false
			    == lengthx * lengthy * lengthz if grid==true
			 */
  bool grid,  /* simulation on a grid required by user? */
    Time;             /* is a time component given ? */

  coord_type xgr;      /* the coordinates */
  double *x,          /* sortiert 1. Koord. 1. Punkt, 2. Koord 1. Punkt, ... */
    T[3];             /* gridtriple definition for the time components --
			 doubled in last coordinate of *x
		       */
} location_type;

typedef struct simu_type {
  bool active,   /* has the init_procedure been called successfully? */
    stop; /* should the simulations be continued? -- appearing in
	     circulant embedding when cepar->severalrealisation=false,
	     may not be merged with active -- otherwise trends cannot be added
	     afterwards
	  */
  int distribution;
  int expected_number_simu;
} simu_type;

typedef struct method_type {
  /* key information */
  globalparam *gp, *gpdo;
  simu_type *simu;
  location_type *loc;


  /* method information */
  int nr;  // negative means primitive (-code-1); positive refers to 
  //          respective covariance function
  bool compatible;  // determines whether simulation result must be stored
    //                 in an intermediate storage since partial results cannot
    //                 be added directly
    int nsub;
  method_type *sub[MAXSUB];
  do_meth domethod;
  destructor destruct;  /* function called to delete intermediate results */
  void *S;       /* pointer to intermediate results;
			    delicate : either ==NULL or should be a genuine 
			               pointer !
			    no hanging pointers allowed!
			 */ 

  /* cov information */
  cov_model *cov; // *redcov;
  double *caniso, cscale, cvar; 
  int *cproj;
  int xdimout;

  matrix_type type;   // von oben
  cov_model *hanging; // at some instance # must be jumped;
  //                     this pointer refers to the the first jump in the row
    // should be obsolete somewhen !


  /* coordinate information */
  double 
    *space,  // trafo of spatial components
    *sptime; // trafo of all components
  static_grid grani;  // grid transformed by (diagonal) aniso
} method_type;





extern char METHODNAMES[][METHODMAXCHAR];
extern char ISONAMES[LAST_ISO + 1][15];
extern char STATNAMES[LAST_STAT + 1][15];
extern init_meth init_method[SimulationTypeLength];
extern bool compatible_primitive[SimulationTypeLength];
extern do_meth do_method[SimulationTypeLength];
SEXP GetMethodInfo(method_type *meth, 
		   bool ignore_active, int depth, int max, long*mem);




///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////

typedef char stationary_type;
typedef char isotropy_type;

/* 
   different definitions for the "dimension" of the field is used
      dim : the formal dimension (given by the user)
      effectivedim: the dimension used in the evaluation of the covariance fct
                    E.G. for isotropic models it is 1
      reduceddim : *  dim of lin. independent subspace for ISOTROPIC models;
                   *  1 + dim of lin. indep. spatial subspace for SPACEISOTROPIC
                   *  dim for ANISOTROPIC covariance models 
*/


#define MAX_RANGES 5
typedef int pref_type[Forbidden + 1];
typedef int pref_shorttype[Nothing + 1];
typedef struct range_type {
  double min[MAXPARAM];
  double max[MAXPARAM];
  bool openmin[MAXPARAM];
  bool openmax[MAXPARAM];
  double pmin[MAXPARAM];
  double pmax[MAXPARAM];
  int maxdim; // default ?! Inf!
  bool finiterange;
} range_type;

typedef struct range_arraytype {
  int n;
  range_type ranges[MAX_RANGES];
} range_arraytype;


typedef double *param_type[MAXPARAM];
typedef struct listoftype {
    double *p[MAXELEMENTS];
    int ncol[MAXELEMENTS], nrow[MAXELEMENTS];
} listoftype;
typedef double *naptr_type[MAX_NA];
// typedef double *internal_type[MAXINTERNALPARAM];

typedef double (*spectral_density)(double *, cov_model *cov); 
typedef struct spectral_storage {
  double phistep2d, phi2d;
  bool grid, ergodic; 
} spectral_storage;


typedef struct spec_covstorage {
  spectral_density density;
  double sigma;
  int nmetro;
  double sub_sd_cum[MAXSUB];
} spec_covstorage;

typedef enum tri {TriFalse, TriTrue, TriMaxFalse, TriBeyond} tri;
/* Tri hat nichts mehr zu bedeuten;
   Wird nur im Zusammenhang mit struct-Eintrag anyNAscaleup und 
   anyNAdown verwendet
   
   
   False:
   True:
   MaxFalse:
   Beyond:
*/


#define DEL_COV -100
typedef struct cov_model {
  // user given information
  int nr; /* number of the covariance function in CovList; 
	     obtained by function GetCovFunction 
	  */
  param_type p;       // 24 b
  int nrow[MAXPARAM], // 24 bytes
    ncol[MAXPARAM];   // 24 bytes 
  double *q;
  int qlen;
  int nsub; /* current number of submodels */
  cov_model *sub[MAXSUB], *calling; 

// Information used and created within simuations
//  SimulationType usermethod; // method given by the user
//    method /* the current method (out of SimulationsType) which 
//			   is tried at the moment or which has been 
//			   successfully initialized */
  spec_covstorage spec;


  /////////////////////////////////////////////////////////
  // VARIABLES PASSED DOWNWARDS
  /////////////////////////////////////////////////////////
  int 
    tsdim, /* current timespacedim  */
    xdim,  /* current xdim (when entering),
	      which may differ from tsdim for isotropic or
	      space-isotropic models   */
    vdim;  /* multivariate dimension passed upwards */


  /////////////////////////////////////////////////////////
  // VARIABLES PASSED UPWARDS
  /////////////////////////////////////////////////////////
  // forward analysis of user's information
  stationary_type statIn; /* formal stationary, s. method->type f. Inhalt
			     IN parameter is usually set fix; 
			  */
                         /*  OUT must be managed by the check function,
			     since OUT may vary among the submodels, by 
			     calling GATTER 
			  */
    isotropy_type isoIn;     /* formal isotropic parameter */

  // backward analysis of user's information
  int maxdim, /* maxdim of model combined with information of the submodels
	        */
  derivatives;  /* number of derivatives of model combined with information 
  		of the submodels */

  bool normalmix, /* for simple model: normal mix model iff maxdim = INFDIM
		     for hypermodel we need nonetheless this parameter
		     parameter set by getinfo(), not within info of the cov 
		     fcts comining information of submodels  
		   */
    finiterange, /* also information obtained by model and submodels 
		   */
    diag,       /*  */
    semiseparatelast, /*  */
    separatelast,     /*  */
    tbm2num,  /* "time" component separately?  */
      hess;  /* can a hessian matrix be provided? */
    

  pref_type pref, /* 0 : not possible; 
		     5 : best possible
		     (including subs)!
		   */
    user; /* 0 : not possible; 
	     5 : best possible
	     (including subs)!
		  */
    double *MLE;
    tri anyNAdown, anyNAscaleup;
    short int manipulating_x; // passed=-1, false, true
} cov_model;

#define X_PASSED -1 /* plus, mal */

#define PREF_BEST 5
#define PREF_NONE 0
// larger than Pref * Nothing
#define LOC_PREF_BEST 9999
#define LOC_PREF_NONE -9999

#define MODEL_USER 0  /* for user call of Covariance etc. */
#define MODEL_SIMU 1  /* for GaussRF etc */ 
#define MODEL_INTERN 2 /* for kriging, etc; internal call of covariance */
#define MODEL_MLE  3  /* do not change !! */
#define MODEL_BOUNDS 4 /* MLE, lower, upper */
#define MODEL_MAX MODEL_BOUNDS + 1
extern cov_model *STORED_MODEL[MODEL_MAX];

#define MPP_MISMATCH -2
#define MPP_AUTO -1
#define MPP_GRID 0
#define MPP_POISS 1
#define MPP_UNIF 2
#define MPP_GAUSS 3

#define MAXLOCALINSTANCES 3
typedef struct localinfotype{
  int instances;
  int msg[MAXLOCALINSTANCES];
  double value[MAXLOCALINSTANCES];
} localinfotype;

typedef void (*rangefct)(cov_model *cov, range_arraytype* ra);
typedef int (*checkfct)(cov_model *cov); 
typedef double (*natscalefct)(cov_model *cov); /* parameters,; natural 
						 scaling */
typedef void (*covfct)(double *, cov_model*, double*); /* h, cov, result */ 
typedef void (*nonstat_covfct)(double *, double*, 
			      cov_model*,  double*); /* x,y, cov, result */
typedef void (*aux_covfct)(double *, double*, double, 
			      cov_model*,  double*); /* x,y, Aux, cov, result */
typedef void (*getlocalparam)(cov_model *, localinfotype *);
typedef bool (*altlocalparam)(cov_model *);
typedef double (*realfct)(cov_model *);

typedef void (*ave_logfct)(double *u, cov_model *, int dim, 
			   double *logg, // log(|g|)
			   double *sign); // sign(g) (not checked if in -1,0,1)
typedef double (*ave_fct)(double *, cov_model *, int dim); /* h,parameters */
typedef void (*avesampling)(double *, double*, int, double*); 
typedef double (*ave_h)(double *, double*, double*, int);
typedef int (*spectral_init)(cov_model *);
typedef void (*spectral_do)(cov_model *, spectral_storage *, double *);

typedef double (*ave_cnst)(double *, double*); /* h,parameters */

typedef struct mpp_storage mpp_storage;
typedef double (*draw_logmix) (cov_model *cov, mpp_storage *s);
typedef double (*log_mixweight)(double *x, cov_model *cov, mpp_storage *s);
typedef void (*mpp_init)(mpp_storage *s, cov_model *cov); /* mpp_init
							   */
typedef void (*sd_fct)(mpp_storage *s, cov_model *cov);
		       				   
typedef void (*mpp_model)(mpp_storage *, cov_model *cov);  /* parameters for 
							  average, RANDOM */
typedef res_type (*mpp_get)(double *, double*, cov_model *cov,
			  mpp_storage *); /* h,parameters */
typedef res_type (*mpp_getstat)(double *, cov_model *cov,
			  mpp_storage *); /* h,parameters */
typedef void (*pointsampler)(mpp_storage *, cov_model *cov);

typedef int (*hyper_pp_fct)(double, double*, double*, cov_model *, bool, 
			    double**, double**, double**);
typedef void (*size_fct)(int i, cov_model *cov, int *nr, int *nc);


// here all the different method for simulating a RF of a specified covariance
// function are stored
typedef struct cov_fct {  
  char name[MAXCHAR];
  char kappas, /* number of parameters  */
    minsub, maxsub; // ==0
  stationary_type stationary;
  isotropy_type isotropy;
    int vdim;
  char kappanames[MAXPARAM][PARAMMAXCHAR],
    subnames[MAXSUB][PARAMMAXCHAR];
  SEXPTYPE kappatype[MAXPARAM];
  size_fct kappasize;

  rangefct range;
  checkfct check;
  int implemented[SimulationTypeLength];
   bool  internal;
  
  pref_shorttype pref;

  natscalefct naturalscale;
    covfct cov, D, D2, D3, D4, tbm2, invSq;
  int derivs;
    nonstat_covfct nonstat_cov;
    covfct nabla, hess;  /* CircEmbed, direct */
  aux_covfct aux_cov; // complicated cov-model that can be used only
  //                     as submodels and that needs an auxiliary argument
  //                     for the evaluation
  getlocalparam coinit, ieinit; // set within primitives
  altlocalparam alternative; //getparam: guess for good local 
    // param (cutoff, intrinsic); alternative gives alternative in a 
    // second or third try (used by co and Stein)

  spectral_init initspectral;
  spectral_do spectral;


  draw_logmix drawlogmix;
  log_mixweight logmixweight;

  mpp_init mppinit;
  mpp_model randomcoin;
  mpp_getstat mppgetstat;
  mpp_get mppget;
  int mpplocations;
  bool timesep;
  ave_fct avef;
  ave_logfct avelogg;
  realfct mineigenvalue;
  sd_fct sd;

  hyper_pp_fct hyperplane;      /* hyperplane tessellations */        
  init_meth initmethod;    /* any other exotic way */
  do_meth domethod;      /* which may be given by a programmer */

  bool primitive;

} cov_fct;

extern int currentNrCov;
typedef struct cov_fct cov_fct;
extern cov_fct *CovList;
extern int GENERALISEDCAUCHY, STABLE, WHITTLE, BROWNIAN, CAUCHY,
    EXPONENTIAL, MATERN, GAUSS, NUGGET, PLUS, MLEPLUS, MIXEDEFFECT, 
    MLEMIXEDEFFECT, MIXX, TBM2NR, VECTOR,
    ELEMENTNR_PLUS, LIST_ELEMENT;
extern char OP_SIGN[][3];
extern int COVLISTALL[MAXSUB];
extern char CovNames[MAXNRCOVFCTS][MAXCHAR];
extern char InternalName[];

// extern int etc, if other are enabled to use these functions?
int IncludePrim(const char *name, int kappas,  
		stationary_type stationary, isotropy_type isotropy,
		checkfct check, rangefct range);
int IncludePrim(const char *name, int kappas,  size_fct kappasize,
		stationary_type stationary, isotropy_type isotropy,
		checkfct check, rangefct range);
int IncludePrim(const char *name, int kappas,  
		stationary_type stationary, isotropy_type isotropy,
		checkfct check, rangefct range, int vdim);
int IncludePrim(const char *name, int kappas,  size_fct kappasize,
		stationary_type stationary, isotropy_type isotropy,
		checkfct check, rangefct range, int vdim);

int IncludePrim(const char *name, int kappas, 
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range, pref_type pref,
		 int vdim);
int IncludePrim(const char *name, int kappas, size_fct kappasize,
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range,pref_type pref, int vdim);

void make_internal();


int IncludeModel(const char *name, char minsub, char maxsub, int kappas,
		  size_fct kappasize,
		  stationary_type stationary, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref, 
		  init_meth initmethod, do_meth domethod,
		  bool internal, int vdim);
int IncludeModel(const char *name, char minsub, char maxsub, int kappas, 
		  size_fct kappasize,
		  stationary_type stationary, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref);
int IncludeModel(const char *name, char minsub, char maxsub, int kappas, 
		  stationary_type stationary, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref);

void addCov(covfct cf, covfct D, natscalefct naturalscale);
void addCov(covfct cf, covfct D, covfct D2, natscalefct naturalscale);
void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    natscalefct naturalscale);
void addCov(nonstat_covfct cf);
void addCov(aux_covfct auxcf);
void nablahess(covfct nabla, covfct hess);
void addSpaceD(covfct spaceD);
void addLocal(getlocalparam coinit, getlocalparam ieinit);
void addCallLocal(altlocalparam alt);
//void addTBM(covfct tbm2, covfct spaceD);
int addTBM(covfct tbm2);
void addTBM(spectral_init initspectral, spectral_do spectral);
void addTBM(covfct tbm2, spectral_init initspectral, spectral_do spectral);
void addInv(covfct invSq);
void addHyper(hyper_pp_fct hyper_pp);
void addMarkov(int*variable);
void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_get mppget, 
	    sd_fct sd, int loc, bool timesep);
void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_get mppget, 
	    sd_fct sd, int loc);
void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_get mppget, 
	    sd_fct sd, int loc,
	    ave_fct avef, ave_logfct avelogg);
void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_getstat mppget, 
	    sd_fct sd, int loc, bool timesep);
void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_getstat mppget, 
	    sd_fct sd, int loc);
void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_getstat mppget, 
	    sd_fct sd, int loc,
	    ave_fct avef, ave_logfct avelogg);
void addSpecial(realfct mineigen);
void addGaussMixture(draw_logmix drawlogmix, log_mixweight logmixweight);

int addFurtherCov(covfct cf, covfct D);
int addFurtherCov(covfct cf, covfct D, covfct D2);
int addFurtherCov(nonstat_covfct cf);

void addusersfunctions();

int check_within_range(cov_model *cov);
int check_within_range(cov_model *cov, bool NAOK);

/* function necessary to set CircEmbed correctly if the function is not 
   even in some coordinates! 
*/
void InitModelList();   
/* absolutely necessary to call this function at the very beginning !
   but done automatically if SimulateRF, etc. is called for the first time
   (direct call necessary if further covariance functions ate to be added by 
   user!)
   initiating CovList with a dozen standard covariance models
*/



///////////////////////////////////////////////////////////////////////
// GLOBAL PARAMETERS AND STRUCTURE
///////////////////////////////////////////////////////////////////////
#define DECISION_TRUE 1
#define DECISION_FALSE 0
#define DECISION_CASESPEC -1
#define MAXERRORSTRING 100 + MAXNRCOVFCTS * (MAXCHAR+1)
#define nErrorLoc 1000
typedef struct key_type {
  globalparam gp, gpdo;
  simu_type simu;
  location_type loc;
  trend_type trend;
  cov_model *cov;
  method_type *meth;
} key_type;
extern char DISTRNAMES[][DISTRMAXCHAR];
extern int PL, NS;
extern double GENERAL_PRECISION;
extern key_type KEY[MAXKEYS]; 
extern char ERRORSTRING[MAXERRORSTRING], ERRORSTRING_OK[MAXERRORSTRING],
  ERRORSTRING_WRONG[MAXERRORSTRING], ERROR_LOC[nErrorLoc];
extern int ERRORMODELNUMBER;
void errorMSG(int error, char* EM);
void ErrorMessage(SimulationType m, int error);
void KEY_DELETE(key_type *key);
void KEY_NULL(key_type *key);
void KEY_DELETE_NOTTREND(key_type *key);
void KEY_NULLNOTTREND(key_type *key);
void TREND_DELETE(trend_type *trend);
void TREND_NULL(trend_type *trend);
void LOC_DELETE(location_type *loc);
void LOC_NULL(location_type *loc);
void METHOD_DELETE(method_type **meth);
void METHOD_NULL( method_type *meth);
void COV_DELETE(cov_model **cov);
void COV_NULL(cov_model *cov);
void COV_DELETE_WITHOUTSUB(cov_model **Cov);
int InternalGetGridSize(double *x[MAXSIMUDIM], int dim, int *lx);
int internal_InitSimulateRF(double *x, double *T, 
			    int spatialdim, /* spatial dim only ! */
			    int lx, bool grid, bool Time, 
			    int distr, /* still unused */
			    key_type *key, globalparam *gp,
			    int expected_number_simu,
			    cov_model **cov);
SEXP InternalGetKeyInfo(key_type *key, bool ignore_active, int depth, int max);
void matrixrotat(double *paramaniso, int col, int row,
		 int *truedim, double* aniso);
SEXP TooLarge(int *n, int l); /* if printout is larger than given size */




///////////////////////////////////////////////////////////////////////
// PARAMETERS AND FUNCTIONS FOR THE SPECIFIC SIMULATION METHODS
///////////////////////////////////////////////////////////////////////
#define TRIVIALSTRATEGY 0
#define STRATEGY_PARTIAL 1
#define LASTSTRATEGY STRATEGY_PARTIAL


// see circembed.cc
typedef struct FFT_storage { double* work; int *iwork, nseg; } FFT_storage;
typedef struct CE_storage {
  int m[MAXCEDIM], trials,
    halfm[MAXCEDIM], nn[MAXCEDIM], cumm[MAXCEDIM+1], 
    cur_square[MAXCEDIM], max_squares[MAXCEDIM], /* !!!! **** */ 
    vdim; //  added by PM 12/08

  long mtot, square_seg[MAXCEDIM];
    double **c, **d, smallestRe, largestAbsIm, *aniso;
  bool positivedefinite,
    new_simulation_next,
    cur_call_odd,
     dependent; // eigentlich braucht es nicht waehrend der initialisierung
    // festgelegt zu werden. Ist aber wesentlich einfacher zu handhaben,
    // da sonst bei internal_dosimulate die parameter alle RFparameter alle
    // nochmals gesetzt werden muessen
  FFT_storage FFT;
} CE_storage;
typedef struct localCE_storage {
    key_type key;
    double *aniso;
    void *correction;
} localCE_storage;
unsigned long NiceFFTNumber(unsigned long nn);
int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *S);
int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *S);
void FFT_destruct(FFT_storage *S);
void FFT_NULL(FFT_storage *S);
void SetParamCE( int *action, int *force, double *tolRe, double *tolIm,
		 int *trials, 
		 double *mmin, /* vector of length MAXCEDIM */
		 int *useprimes, int *strategy, double *maxmem, int* dependent,
		 ce_param *CE, const char *name);
int init_circ_embed(method_type *meth);
void do_circ_embed(method_type *meth, res_type *res);
int init_circ_embed_intr(method_type *meth);
void do_circ_embed_intr(method_type *meth, res_type *res );
int init_circ_embed_cutoff(method_type *meth);
void do_circ_embed_cutoff(method_type *meth, res_type *res );

// see tbm.cc			  
typedef struct TBM_storage {
  // bool genuine_dim[MAX362DIM];
  int ce_dim, simuspatialdim, spatialtotalpts;
  double center[MAXTBMSPDIM],  linesimuscale, *aniso;
  res_type *simuline; 
  key_type key;
} TBM_storage;
int init_turningbands2(method_type *meth);
int init_turningbands3(method_type *meth);
void do_tbm2(method_type *meth, res_type *res);
void do_tbm3(method_type *meth, res_type *res);
extern SimulationType TBM_METHOD;


// see spectral.cc
int init_spectral(method_type *meth);
void do_spectral(method_type *meth, res_type *res );
void metropolis(cov_model *cov, double *x);
int search_metropolis(cov_model *cov);
// see direct.cc
typedef struct direct_storage {
  InversionMethod method;
  double *U, *G;
} direct_storage;
int init_directGauss(method_type *meth);
void do_directGauss(method_type *meth, res_type *res);

// see sequential.cc
typedef struct sequential_storage {
  int back, totpnts, spatialpnts, ntime, initial;
  double *U22, *G, *MuT, *U11,  *Cov21, *Inv22;
    res_type *res0;
} sequential_storage;
int init_sequential(method_type *meth);
void do_sequential(method_type *meth, res_type *res);


// see Rue.cc
// for typedef struct markov_storage see rue.cc	   
int init_markov(method_type *meth);
void do_markov(method_type *meth, res_type *res );

// see averages.cc
#define MPP_VV 1
#define MPP_SPECTRAL 2
#define MPP_EIGEN 3

#define AVE_SPECTRAL 4
#define AVE_VV 5

#define DRAWMIX_EXPONENT 6
#define DRAWMIX_ALPHA 7
#define DRAWMIX_FACTOR 8
#define DRAWMIX_P 9
#define DRAWMIX_BETA 10

#define LASTC 10


// nugget
typedef struct nugget_storage {
    //double sqrtnugget;
  bool simple, simugrid;
  int *pos, reduced_dim[MAXNUGGDIM], prod_dim[MAXNUGGDIM + 1];
  res_type *red_field;
  
} nugget_storage;
int init_nugget(method_type *meth);
void do_nugget(method_type *meth, res_type *res );

// ave.cc
int init_ave(method_type *meth);
void do_ave(method_type *meth, res_type *res);

// see MPP.cc	
typedef struct mpp_storage {
  // general
  // avoid the same names as in addBool_storage, RFbool.h
  double integral, integralsq, // add specific
    effectiveRadius, effectivearea,
    plus, relplus, lensimu[MAXMPPDIM], 
    min[MAXMPPDIM], max[MAXMPPDIM],
    minsimu[MAXMPPDIM], maxsimu[MAXMPPDIM],
    mean[MAXMPPDIM], sdgauss[MAXMPPDIM],
    u[MAXMPPDIM],
    c[LASTC + 1], /* local memory for diverse constants,
	     depending on the specific model */
    logapproxzero, samplingdist, samplingr,
    average, factor, trueintensity, intensity, griddist,
    *aniso;
 
  // extremes specific AND for correction

  int dim, nmax[MAXMPPDIM], n[MAXMPPDIM],
    ntot;

  // variable parst and  model specific parts
  // ----------------------------------------

  // standard models
    double r, rsq, spherical_adddist;
 
  // non-standard models
  double invscale, loginvscale, logInvSqrtDens;

  // extremes
  double maxheight, integralpos;

  pointsampler location;
  sd_fct getsd;

} mpp_storage;
int init_mpp(method_type *meth);
void do_addmpp(method_type *meth, res_type *res );
int init_mppave(method_type *meth, int dim, bool timesep);
void approx_integral(ave_logfct logg, cov_model *cov, mpp_storage *s);


// dummy version, hyperplane
typedef struct hyper_storage{
    double rx[MAXHYPERDIM], center[MAXHYPERDIM], radius, *aniso;
  hyper_pp_fct hyperplane;
} hyper_storage;
int init_hyperplane(method_type *meth);
void do_hyperplane(method_type *meth, res_type *res);


//extremes
typedef struct extremes_storage{
    double inv_mean_pos, assumedmax;
    res_type *rf;
    key_type key;
    SimulationType meth;
} extremes_storage;


///////////////////////////////////////////////////////////////////////
// POSITIONS OF X-VALUES
///////////////////////////////////////////////////////////////////////
#define XSTART 0
#define XEND   1
#define XSTEP  2

void Transform2NoGrid(method_type *meth, bool timesep);
void Transform2NoGridExt(method_type *m, bool timesep, bool gridexpand, 
			 double **Space, double **Sptime);


///////////////////////////////////////////////////////////////////////
// TREND SPECIFICATION
///////////////////////////////////////////////////////////////////////
#define TREND_MEAN 0
#define TREND_LINEAR 1
#define TREND_FCT 2
#define TREND_PARAM_FCT 3



///////////////////////////////////////////////////////////////////////
// Fortran subroutines
///////////////////////////////////////////////////////////////////////

//EXTERN void F77_NAME(ave1f)(double *x, double *u,
//			    double *h, double *Mz,
//			    double *r, int *d, double *res);

extern int DOLLAR, LASTDOLLAR, GATTER,  LASTGATTER, MLE_ENDCOV, OUT, SPLIT,
    ISO2ISO,SP2SP, SP2ISO, S2ISO, S2SP, S2S, NATSC, TREND, CONSTANT;

extern char MSG[1000];
extern char NEWMSG[1000];
#define UERR if (PL>4) { PRINTF("\n\n================================\n");\
  if (STORED_MODEL[MODEL_USER]!=NULL) \
    PrintModelInfo(STORED_MODEL[MODEL_USER]); \
  if (STORED_MODEL[MODEL_INTERN]!=NULL) \
    PrintModelInfo(STORED_MODEL[MODEL_INTERN]); \
  if (STORED_MODEL[MODEL_SIMU]!=NULL) \
    PrintModelInfo(STORED_MODEL[MODEL_SIMU]); }\
 PRINTF("---------------------------------\n");\

#define ERR(X) {UERR; sprintf(MSG, "%s\n%s", ERROR_LOC, X); error(MSG);}
#define XERR(X) {errorMSG(X, MSG); \
 sprintf(NEWMSG, "in `%s' error %d: %s", ERROR_LOC, X, MSG); error(NEWMSG);}
#define PERR(X) {UERR;sprintf(MSG, "%s\n%s: %s", ERROR_LOC, param_name, X); error(MSG);}

void ErrCov(double *x, cov_model *cov, double *v);
void ErrCovNonstat(double *x, double *y, cov_model *cov, double *v);
void GetNaturalScaling(cov_model *cov, double *natscale);
void analyse_matrix(double *aniso, int row, int col,
		    bool *diag, bool *quasidiag, int *idx,
		    bool *semiseparatelast, bool *separatelast);

//void covcpy(cov_model **localcov, cov_model *cov);
void covcpy(cov_model **localcov, cov_model *cov, bool insertgatter,
	    bool keepuser);
int getmodelnr(char *name);

typedef int (*set_local_q_type) (cov_model *next, double a, double d, 
				  double *q);
//void setdefault(cov_model *cov);
void setbackward(cov_model *cov, cov_model *sub);

int getmodelnr(char *name);

void kappanames(const char* n1, SEXPTYPE t1);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6);
void subnames(const char* n1);
void subnames(const char* n1, const char* n2);
void subnames(const char* n1, const char* n2, const char* n3);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5, const char* n6);

void setAniso(method_type *meth);
SEXP CheckModelInternal(SEXP model, SEXP tsdim, SEXP xdim, cov_model **cov,
    int maxdim);
void orderingInt(int *d, int len, int dim, int *pos);
int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int lx,  bool Time, bool grid,
	    location_type *loc, globalparam *gp);
int Match(char *name, name_type List, int n);

extern name_type STANDARDPARAM, STANDARDSUB;

void range_default(range_arraytype *ra);
void CMbuild(SEXP model, int tsdim, int xdim, int level, cov_model **Cov,
	     cov_model previous);
void CheckModelInternal(SEXP model, int tsdim, int xdim,
			bool stationary, 
			cov_model **Cov,
			int maxdim);
void CheckModel(SEXP model, SEXP tsdim, SEXP xdim, SEXP stationary, 
		cov_model **Cov,
		int maxdim);
void PrintModelInfo(cov_model *cov);

int initstandard(method_type *meth);
void dostandard(method_type *meth, double *res);

void cpyMethod(method_type *meth, method_type *sub, bool cp_aniso);
//double * Aniso(double *x, method_type *meth);
//void Aniso(double *x, method_type *meth, double *res);
//void Aniso(double *x, double *caniso, int orgidim, int dim,  double *res);
matrix_type Type(double *m, int nrow, int ncol);
double GetDiameter(double *origcenter, double *min, double *max, 
		   double *aniso, int tsdim, int dim);
double GetDiameter(double *origcenter, double *oirgmin, double *origmax, 
		   double *aniso, int origdim, int spatialdim,
		   double *min, double *max, double *center);
double * matrixmult(double *m1, double *m2, int dim1, int dim2, int dim3);
void GetMinMax(method_type *meth, double *min, double *max, double *center,
    int MaxDim);
void addModel(cov_model **pcov, int covnr);


#define MULTIPLEMATCHING -2
#define NOMATCHING -1
#define MATCHESINTERNAL -3

#define UNSET -1

int init_nothing(method_type *meth);
void do_nothing(method_type *meth, res_type *res );

void PrintMethodInfo(method_type *meth);
void PrintModelList();


void E1(spectral_storage *s, double A, double *e);
void E2(spectral_storage *s, double A, double *e);
void E3(spectral_storage *s, double A, double *e);
void E(int dim, spectral_storage *s, double A, double *e);


void gauss_initu(mpp_storage *s); 
void gauss_u(mpp_storage *s, cov_model *cov); 
void unif_initu(mpp_storage *s); 
void unif_u(mpp_storage *s, cov_model *cov);
void poissgauss_initu(mpp_storage *s);
void poiss_initu(mpp_storage *s); 
void grid_initu(mpp_storage *s);
void grid_u(mpp_storage *s, cov_model *cov);

int checkkappas(cov_model *cov);
int checkkappas(cov_model *cov, bool errornull);
double gaussInt(int d, int xi, double sigma, double R);
cov_model * covOhneGatter(cov_model *cov);

void updatepref(cov_model *cov, cov_model *sub);

extern char PREF_FAILURE[100 * Nothing];
void leer(int level);
typedef enum sortsofparam { // never change ordering; just add new ones !!!!!!
  VARPARAM, SCALEPARAM, DIAGPARAM, ANISOPARAM, // 0..3
  INTERNALPARAM, ANYPARAM, TRENDPARAM, NUGGETVAR, MIXEDVAR, // 4..8
  REGRESSION // 9 
} sortsofparam;
typedef enum sortsofeffect{
  deteffect, fixedeffect, randomeffect, rvareffect,
  largeeffect, lvareffect, spaceeffect, remaining, eff_error
} sortsofeffect;

void get_ranges(cov_model *cov, cov_model **min, cov_model **max, 
	       cov_model **openmin, cov_model **openmax,
		bool practical);
int check_recursive_range(cov_model *cov, bool NAOK);
int internal_DoSimulateRF(key_type *key, int nn, res_type *orig_res);

void strcopyN(char *dest, const char *src, int n);
void expandgrid(coord_type xgr, int *len, double **xx, int nrow);
void xtime2x(double *x, int nx, double *T, int len, double **newx, int nrow);
void removeOnly(cov_model **Cov);

SEXP Logi(bool* V, int n, int max, long *mem) ;
SEXP Num(double* V, int n, int max, long *mem) ;
SEXP Int(int *V, int n, int max, long *mem) ;
SEXP Char(char *V, int n, int max, long *mem) ;
SEXP Mat(double* V, int row, int col, int max, long *mem);
SEXP MatInt(int* V, int row, int col, int max, long *mem) ;
SEXP Array3D(int** V, int depth, int row, int col, int max, long *mem) ;
SEXP TooLarge(int *n, int l);

//int Checking(cov_model **Cov);
int check2X(cov_model *cov, int tsdim, int xdim,
	    stationary_type statprev, isotropy_type isoprev, int vdim);

void kdefault(cov_model *cov, int i, double v);
void DeleteGatter(cov_model **Cov);

extern pref_type PREF_ALL, PREF_NOTHING;

#define MAXDEFMATRIX 3
//extern double *userdefinedCovMatrix[MAXDEFMATRIX][MAXMAKEEXPLICITE];
// extern int userdefinedCM_RC[MAXDEFMATRIX][MAXMAKEEXPLICITE];
extern int CovMatrixRow, CovMatrixCol, CovMatrixTotal;
// extern bool CovMatrixIndexActive;
void zeronugget(int* zeronugget);
				
void EinheitsMatrix(double **mem, int dim);
void cum_dollar(method_type *meth, int origdim, cov_model *cov,
		method_type *sub);
				
extern "C" double *OutX, *OutY;

void CovarianceMatrix(double *x, bool dist, int lx, cov_model *cov, 
		      double *result);

// extern int CumIdxMakeExplicite[MAXMAKEEXPLICITE];
double *getAnisoMatrix(method_type *meth);
double *getAnisoMatrix(cov_model *cov, int dim);



#define AveMaxDim 10 /* nur technisch */
#define CoxMaxDim 3  /* nur technisch ?! */
#define StpMaxDim 10  /* nur technisch ?! */
#define EaxxaMaxDim 10  /* nur technisch */
#define UmatrixMaxDim MAXMLEDIM /* nur technisch */
#define ShiftMaxDim 10  /* nur technisch ?! */


void vstop();
void vassert();
void vassert(bool b);


double XkCXtl(double *X, double *C, int nrow, int dim, int k, int l);
void XCXt(double *X, double *C, double *V, int nrow, int dim);
void xA(double *x, double*A, int nrow, int ncol, double *y);
void Ax(double *A, double*x, int nrow, int ncol, double *y);

SEXP Param(void* p, int nrow, int ncol, SEXPTYPE type, bool drop, long* mem);
void CovIntern(double *x, double *y, bool y_not_given,  int lx, double *result);
int is_positive_definite(double *C, int dim);

// typedef long long POINTER;
typedef long int POINTER;
	
extern int MAX_PMI;				


//extern void
//F77_NAME(cpotf2)(int* uplo, int* n, Rcomplex* a, int* lda, int* info);

//extern void
//F77_NAME(zpotf)(char* uplo, int* n, Rcomplex* a, int* lda, int* info);

//extern void F77_NAME(zpotf)(int* info);
extern "C" void F77_CALL(zpotf2)(char *name, int *row, Rcomplex *U, int *xxx,
				 int *Err);
				

				
#endif


 
