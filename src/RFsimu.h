#ifndef RFsimu_H
#define RFsimu_H 1

#include "auxiliary.h"
#include "RandomFields.h"
#include <string.h>
#include "error.h"

///////////////////////////////////////////////////////////////////////
// BASIC DIMENSIONS AND VARIABLES
///////////////////////////////////////////////////////////////////////
#define DISTRMAXCHAR 10
#define DISTRNR 3
#define METHODMAXCHAR 31 /* max number of character to describe a Method (including \0)*/
#define COVMAXCHAR 14 /* max number of characters for covariance name*/
#define MAXCOV 10     /* max number of covariance functions for a random 
			 field -- NEVER MORE THAN 255 */
#define MAXDIM 4          
#define MAXDIMSQ MAXDIM * MAXDIM
#define ZWEIHOCHMAXDIM 16 /* 2^MAXDIM */
/* if explicite simulation method is not given by
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
#define MAXNRCOVFCTS 40
extern double ZERO[MAXDIM];
extern double UNIT[MAXDIM];


typedef struct key_type;
typedef int covlist_type[MAXCOV];



///////////////////////////////////////////////////////////////////////
// PARAMETER AND THEIR LOCATIONS
///////////////////////////////////////////////////////////////////////
#define VARIANCE   0   /* without NUGGET */
#define KAPPA      1   /* Kappas should be always the last ones of transferred 
			  parameters */ 
#define KAPPAI     KAPPA
#define KAPPAII    2  
#define KAPPAIII   3
#define KAPPAIV    4
#define KAPPAV     5
#define KAPPAVI    6
#define KAPPAVII   7
#define LASTKAPPA  KAPPAVII
#define TBM2NUM          LASTKAPPA + 1
#define HYPERINTERNALI   LASTKAPPA + 2
#define HYPERINTERNALII  LASTKAPPA + 3
#define HYPERINTERNALIII LASTKAPPA + 4
#define HYPERINTERNALIV  LASTKAPPA + 5
#define LASTHYPERINTERNAL HYPERINTERNALIV
#define SCALE LASTHYPERINTERNAL + 1
#define ANISO SCALE
#define TOTAL_PARAM ANISO + MAXDIM * MAXDIM
/*  LASTKAPPA + MAXDIM^2 + 1
    if major changings check CheckAndBuildCov, InitSimulateRF (both C & R), 
    hyperbolic, TBM3hyperbolic, twodimfractalbrownian 
    R if changed! */

#define HYPERNR KAPPAI
#define HYPERKAPPAI KAPPAII
#define HYPERKAPPAII KAPPAIII
#define HYPERKAPPAIIII KAPPAIV
#define HYPERKAPPAIV KAPPAV
#define LASTHYPERKAPPA LASTKAPPA

// parameters used in local cov fcts
#define LOCAL_R HYPERINTERNALI
#define CUTOFF_B HYPERINTERNALII
#define CUTOFF_ASQRTR HYPERINTERNALIII
#define CUTOFF_THEOR HYPERINTERNALIV
#define INTRINSIC_A0 HYPERINTERNALII
#define INTRINSIC_A2 HYPERINTERNALIII
#define INTRINSIC_B HYPERINTERNALIV

#define DIAMETER HYPERKAPPAI /* used in local circulant embedding models */
#define CUTOFF_A HYPERKAPPAII
#define INTRINSIC_RAWR HYPERKAPPAII

typedef double param_type[TOTAL_PARAM];
typedef double param_type_array[MAXCOV][TOTAL_PARAM];
typedef double aniso_type[MAXDIMSQ];



///////////////////////////////////////////////////////////////////////
// GENERAL PARAMETERS FOR THE SIMULATIONMETHODS
///////////////////////////////////////////////////////////////////////
#define SET_DESTRUCT(A, m)\
  assert(key->active==false);\
  assert((key->meth[m].S==NULL) && (key->meth[m].destruct==NULL));\
  key->meth[m].destruct = A;
typedef enum SimulationType {
  /* never change the ordering -- at least check with "compatible" and 
     InternalGetKeyInfo*/
  CircEmbed,     /* Circulant embedding -- must always be the first one! */
  CircEmbedCutoff,
  CircEmbedIntrinsic,
  TBM2, TBM3,    /* Turning Bands, performed in 2 or 3 dimensions */
  SpectralTBM,   /* Spectral turning bands, only 2 dimensional */
  Direct,        /* directly, by matrix inversion */
  Nugget,        /* just simulate independent variables */
  AdditiveMpp,   /* "additive MPP (random coins)", will be only in the plane*/
  Hyperplane,    /* by superposition of Poisson hyperplane tessellations;
		    only in the plane! Not implemented yet*/
  Special,       /* left open */
  Nothing,       /* must always be the last of the list of methods for
		    simulation Gaussian random fields
		 */
  MaxMpp,
  ExtremalGauss,
  Forbidden      /* must always be the last one */
} SimulationType; 
typedef SimulationType LocalType;
#define SimulationTypeLength ((int) Forbidden) - 1

typedef void (*destructor)(void **);
typedef struct methodvalue_type {
  void *S;       /* pointer to intermediate results;
			    delicate : either ==NULL or should be a genuine 
			               pointer !
			    no hanging pointers allowed!
			 */ 
  destructor destruct;  
                        /* function called to delete intermediate results */
  SimulationType unimeth;
  bool unimeth_alreadyInUse, incompatible;
  int actcov;
  covlist_type covlist;  
} methodvalue_type;
typedef methodvalue_type methodvalue_arraytype[MAXCOV];
typedef int (*initmethod)(key_type*, int); /* h,parameters */
typedef void (*do_comp_meth)(key_type*, int, double*); /* h,parameters */
typedef void (*do_incomp_meth)(key_type*, int, double*); /* h,parameters */
extern char METHODNAMES[][METHODMAXCHAR];
extern initmethod init_method[SimulationTypeLength];
extern do_comp_meth do_compatiblemethod[SimulationTypeLength];
extern do_incomp_meth  do_incompatiblemethod[SimulationTypeLength];
SEXP GetMethodInfo(key_type *key, methodvalue_arraytype keymethod);



///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////
// type of covariance functions that need distinguished treatment
#define FULLISOTROPIC 0
#define SPACEISOTROPIC 1
#define ANISOTROPIC 2
#define ISOHYPERMODEL 3
#define ANISOHYPERMODEL 4

// way of implementing simulation methods:
#define NOT_IMPLEMENTED 0 /* do not change this value except to false */
#define IMPLEMENTED 1     /* do not change this value except to true */
#define NUM_APPROX 2 
#define GIVEN_METH_IGNORED 3
#define HYPERIMPLEMENTED 4

typedef struct covinfo_type {
  SimulationType method;
  /* the current method (out of SimulationsType) which 
     is tried at the moment or which has been 
     successfully initialized */
  int dim, truetimespacedim,
    length[MAXDIM], /* if simugrid: what is the genuine extension of the grid?
		       That is, length is one where the diag element is 0
		    */
    idx[MAXDIM], /* if ANISO=diag, what is the new column in aniso,
		    where the columns with 0 diag elements are put
		    to the end */
    nr, /* number of the covariance function in CovList; 
	   obtained by function GetCovFunction */
    op; /* operator after the specified covariance function */
  
  bool genuine_dim[MAXDIM], /* if ANISO==diag, is d-th diag element != 0? */
    genuine_last_dimension, /* 
			       i.e. time component indicated and at least
			       one component of the last column of 
			       aniso is different from 0;
			       extra treating and not within aniso (where
			       time component is not reduced then) necessary
			       for SPACEISOTROPIC covriance functions,
				and the MaStein hyper model*/
    simugrid, /* can the simulation technically be performed on a grid ?*/
    left;        /* left to be considered in the simulation
		    initialisations ! */ 
  param_type param;
  aniso_type aniso;
  double *x; // this pointer may only be used for the coordinates
  // transformed by Transform2NoGrid; since this operation is time and 
  // memory consuming but not needed by all methods, it is performed only
  // when needed, but then kept for the next method if the calling 
  // method fails
} covinfo_type;
typedef covinfo_type covinfo_arraytype[MAXCOV];
typedef double (*CovFctType)(double *, int, covinfo_arraytype,
			       covlist_type covlist, int, bool);

typedef void (*infofct)(double *, int *, int *); /* param, maxdim, 
      CEbadlybehaved? --  to create preference list */
typedef void (*rangefct)(int, int *, double*);
/* DO NOT change double to anything else since memcpy  fcts are used
   dim, class parameter like in nsst, 
   return 2 x length(param) x {theor, pract }*/
typedef int (*checkfct)(double*, int, SimulationType); /* h,parameters; 
							covariance fct, nr. */
typedef int (*checkhyper)(covinfo_arraytype, covlist_type, int, 
			  SimulationType);
typedef double (*natscalefct)(double *, int); /* parameters, ; natural 
						 scaling */
typedef double (*covfct)(double *, double* , int); /* h,parameters */
typedef int (*getparamfct)(covinfo_type *, param_type, int);
typedef bool (*alternativeparamfct)(covinfo_type *, int);
typedef double (*isofct)(double*, double*); /* h,parameters */
typedef double (*randommeasure)(double *);     /* parameters, RANDOM */
typedef struct mpp_storage;
typedef double (*mppmodel)(double *); /* point of random field */
typedef void (*initmppfct)(mpp_storage*, int dim, 
			   param_type param); /* h,parameters */
typedef void (*MPPRandom)(mpp_storage *s,
			  double *min, double *max,
			  mppmodel*  /* the random function/model, returned */
                          );
typedef int (*hyper_pp_fct)(double, double*, double*, int, bool, 
			    double**, double**, double**);
typedef int (*generalSimuInit)(key_type*);   
typedef void (*generalSimuMethod)(key_type*, double*); 
typedef int (*kappas_type)(int);

// here all the different method for simulating a RF of a specified covariance
// function are stored
typedef struct cov_fct {  
  char name[COVMAXCHAR];
  unsigned int exception;
  kappas_type kappas; /* number of parameters additional to the standard ones;
		 standard ones are variance and scale */
  char type;
  infofct info;
  bool variogram, even, odd[MAXDIM]; // even and odd not used yet, but set
  rangefct range;
  checkfct check;
  checkhyper checkNinit;             
  int implemented[SimulationTypeLength];
 
  natscalefct naturalscale;
  covfct cov;                             /* CircEmbed, direct */
  LocalType localtype;  
    getparamfct getparam;
  alternativeparamfct alternative; //getparam: guess for good local 
    // param (cutoff, intrinsic); alternative gives alternative in a 
    // second or third try.

  isofct derivative, secondderivt, cov_tbm2, cov_tbm3; // local ce and TBM 
  randommeasure spectral;       /* spectral tbm */

  initmppfct add_mpp_scl;
  MPPRandom add_mpp_rnd; 
  hyper_pp_fct hyperplane;      /* hyperplane tessellations */        
  generalSimuInit initother;    /* any other exotic way */
  generalSimuMethod other;      /* which may be given by a programmer */
} cov_fct;

extern int currentNrCov;
typedef struct cov_fct;
extern cov_fct *CovList;
extern int GENERALISEDCAUCHY, STABLE, WHITTLEMATERN, BROWNIAN, CAUCHY,
  EXPONENTIAL;
extern char OP_SIGN[][3];
extern int COVLISTALL[MAXCOV];

// extern int etc, if other are enabled to use these functions?
int IncludeModel(char *name, int kappas, checkfct check,
			int isotropic, bool variogram, infofct info,
			rangefct range);	
extern int IncludeModel(char *name, kappas_type kappas, checkfct check,
			int type, bool variogram, infofct info,
			rangefct range);
int IncludeHyperModel(char *name, int kappas, checkhyper check,
		      int type, bool variogram, infofct info, rangefct range);
int IncludeHyperModel(char *name, kappas_type kappas, checkhyper check,
		      int type, bool variogram, infofct info, rangefct range);
void modelexception(int nr, bool scale, bool aniso);
void addSimu(int nr, SimulationType r1,SimulationType r2,SimulationType r3);
void addCov(int nr, covfct cov, isofct derivative, natscalefct naturalscale);
void addLocal(int nr, bool cutoff, isofct secondderiv, int *variable);
void addInitLocal(int nr, getparamfct getparam,
			 alternativeparamfct alternative, LocalType type);
void addTBM(int nr, isofct cov_tbm2, isofct cov_tbm3, randommeasure spectral);
void addOther(int nr, initmppfct add_mpp_scl, MPPRandom add_mpp_rnd, 
	      hyper_pp_fct hyper_pp, generalSimuInit initother, 
	      generalSimuMethod other);
void addodd(int nr, int dim); 
void addusersfunctions();

/* function necessary to set CircEmbed correctly if the function is not 
   even in some coordinates! 
*/
int FirstCheck_Cov(key_type *key, int m, bool MultiplyandHyper);
double CovFct(double *x, int dim, covinfo_arraytype keycov,
		 int covlist[MAXCOV], int ncov, bool anisotropy);
double DerivCovFct(double *x, int dim, covinfo_arraytype keycov,
		   covlist_type covlist, int ncov);
double SndDerivCovFct(double *x, int dim, covinfo_arraytype keycov,
		      covlist_type covlist, int ncov);
double CovFctTBM2(double *x, int dim, covinfo_arraytype keycov,
		  covlist_type covlist, int ncov, bool anisotropy);
void InitModelList();   
/* absolutely necessary to call this function at the very beginning !
   but done automatically if SimulateRF, etc. is called for the first time
   (direct call necessary if further covariance functions ate to be added by 
   user!)
   initiating CovList with a dozen standard covariance models
*/
SEXP GetModelInfo(covinfo_arraytype keycov, int nc, int totalparam, 
		  long totalpoints);


int kappaZero(int dim);
int kappaOne(int dim); 
int kappaTwo(int dim); 
int kappaThree(int dim); 
int kappaFour(int dim); 
int kappaFive(int dim);
int kappaSix(int dim); 
int kappaSeven(int dim); 
int kappaFalse(int dim); 


///////////////////////////////////////////////////////////////////////
// GLOBAL PARAMETERS AND STRUCTURE
///////////////////////////////////////////////////////////////////////
#define DECISION_TRUE 1
#define DECISION_FALSE 0
#define DECISION_CASESPEC -1
#define MAXERRORSTRING 100 + MAXNRCOVFCTS * (COVMAXCHAR+1)
#define nErrorLoc 100
typedef struct decision_param{
  int stationary_only, exactness;
} decision_param;
typedef struct key_type {
  bool active,   /* has the init_procedure been called successfully? */
    stop, /* should the simulations be continued? -- appearing in
	     circulant embedding when cepar->severalrealisation=false,
	     may not be merged with active -- otherwise trends cannot be added
	     afterwards
	  */
    anisotropy, /* is it an isotropic covariance function ? */
    compatible,     // determines whether simulation result must be stored
    //                 in an intermediate storage since partial results cannot
    //                 be added directly
    grid,  /* simulation on a grid required ? */
    storing, /* value of GENERAL_STORING at initialisation time */
    Time;             /* is a time component given ? */
  
  int distribution,
    ncov,             /* number of covariance function actually used */
    timespacedim,      /* total dimension, including time and space */
    totalparam, 
    length[MAXDIM], /* grid==true: number of pixels of each direction 
		       grid==false: length[0]=number of spatial points in total
		     */
    spatialdim,
    NML[MAXCOV], IML[MAXCOV], /* NML: number of considered methods: 
				 IML: array index for current method */ 
    lTrendFct, // only used for storing -- used in R
    lLinTrend, // only used for storing -- used in R
    TrendModus;                 /* spatial dimension of RF */
  unsigned char n_unimeth;
  char *TrendFunction; // only used for storing -- used in R
  long spatialtotalpoints,
    totalpoints;      /* number of pixels in total;
			    == lx if grid==false
			    == lengthx * lengthy * lengthz if grid==true
			 */

  double *x[MAXDIM],       /* the coordinates */
    mean,         // only used if storagemodus=0 (see rf.R)
    *LinearTrend, // only used for storing -- used in R
    T[3];             /* gridtriple definition for the time components --
			 doubled in last coordinate of *x
		       */
  // int timelength;        /* number of points in the time direction */

  CovFctType covFct; // tbm covariance function or the usual one ?
  covinfo_arraytype cov;
  methodvalue_arraytype meth;
  SimulationType ML[MAXCOV][SimulationTypeLength];/* for storing the sequence 
						     of considered  methods */ 
} key_type;
extern char DISTRNAMES[][DISTRMAXCHAR];
extern int GENERAL_PRINTLEVEL;
extern int GENERAL_NATURALSCALING;
extern int GENERAL_STORING;
extern double GENERAL_PRECISION;
extern decision_param DECISION_PARAM;
extern char GENERAL_PCH[2]; 
extern key_type KEY[MAXKEYS]; 
extern char ERRORSTRING_OK[MAXERRORSTRING],
  ERRORSTRING_WRONG[MAXERRORSTRING], ERROR_LOC[nErrorLoc];
extern int ERRORMODELNUMBER;
void ErrorMessage(SimulationType m, int error);
void DeleteKeyNotTrend(key_type *key);
void KEY_NULL(key_type *key);
int internal_InitSimulateRF(double *x, double *T, 
			    int spatialdim, /* spatial dim only ! */
			    int lx, bool grid, bool Time,
			    int *covnr, double *ParamList, int nParam,
			    int ncov, bool anisotropy, int *op,
			    int *method, 
			    int distr, /* still unused */
			    key_type *key,
			    /* bool storing, */ 
			    int natural_scaling,
			    CovFctType covFct);
int internal_DoSimulateRF(key_type *key, int nn, double *res);
SEXP InternalGetKeyInfo(key_type *key, bool ignore_active, int depth, int max);
void printkey(key_type *key);
void GetTrueDim(bool anisotropy, int timespacedim, param_type param,
		char type, bool *Time, int *truetimespacedim, 
		aniso_type aniso);
SEXP TooLarge(int *n, int l); /* if printout is larger than given size */




///////////////////////////////////////////////////////////////////////
// PARAMETERS AND FUNCTIONS FOR THE SPECIFIC SIMULATION METHODS
///////////////////////////////////////////////////////////////////////
#define TRIVIALSTRATEGY 0
#define STRATEGY_PARTIAL 1
#define LASTSTRATEGY STRATEGY_PARTIAL
typedef struct ce_param {
  bool force, useprimes, dependent, severalrealisations;
  char strategy;
  int trials;
  double maxmem, tol_re, tol_im, mmin[MAXDIM];
} ce_param;
extern ce_param CIRCEMBED, LOCAL_CE;
typedef struct local_user_param{
  double cutoff_a, intrinsic_r;
} local_user_param;


// see RFcircembed.cc
typedef struct FFT_storage { double* work; int *iwork, nseg; } FFT_storage;
typedef struct CE_storage {
  int m[MAXDIM], halfm[MAXDIM], nn[MAXDIM], cumm[MAXDIM+1], 
    mtot, cur_square[MAXDIM], max_squares[MAXDIM]; /* !!!! **** */
  long square_seg[MAXDIM];
  double *c,*d;
  bool new_simulation_next,
     dependent; // eigentlich braucht es nicht waehrend der initialisierung
    // festgelegt zu werden. Ist aber wesentlich einfacher zu handhaben,
    // da sonst bei internal_dosimulate die parameter alle RFparameter alle
    // nochmals gesetzt werden muessen
  FFT_storage FFT;
} CE_storage;
typedef struct localCE_storage {
    key_type key;
    void *correction[MAXCOV];
} localCE_storage;
unsigned long NiceFFTNumber(unsigned long nn);
int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *S);
int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *S);
void FFT_destruct(FFT_storage *S);
void FFT_NULL(FFT_storage *S);
void SetParamCE( int *action, int *force, double *tolRe, double *tolIm,
		 int *trials, int *severalrealisations,
		 double *mmin, /* vector of length MAXDIM */
		 int *useprimes, int *strategy, double *maxmem, int* dependent,
		 ce_param *CE, char *name);
int init_circ_embed(key_type * key, int m);
void do_circ_embed(key_type *key, int m, double *res);
int init_circ_embed_local(key_type * key, int m);
void do_circ_embed_local(key_type *key,  int m, double *res );
#define nLocalCovList 10
extern int LocalCovList[nLocalCovList];
extern int nlocal;


// see RFtbm.cc			  
typedef struct TBM_storage {
  aniso_type aniso;
  bool simugrid, genuine_dim[MAXDIM];
  int simuspatialdim, ce_dim, truetimespacedim;
  double center[MAXDIM], *simuline, *x;
  key_type key;
} TBM_storage;
int init_turningbands2(key_type *key, int m);
int init_turningbands3(key_type *key, int m);
void do_turningbands(key_type *key, int m, double *res );
extern SimulationType TBM_METHOD;
extern bool TBM2NUMERIC;


// see RFspectral.cc			   
typedef struct spectral_storage {
  randommeasure randomAmplitude[MAXCOV];
} spectral_storage;
int init_simulatespectral(key_type *key, int m);
void do_simulatespectral(key_type *key, int m, double *res );


// see RFdirect.cc
extern int DIRECTGAUSS_BESTVARIABLES;
extern int DIRECTGAUSS_MAXVARIABLES;
typedef enum InversionMethod { Cholesky, SVD, NoFurtherInversionMethod } 
InversionMethod;
typedef struct direct_storage {
  InversionMethod method;
  double *U,*G;
} direct_storage;
int init_directGauss(key_type *key, int m);
void do_directGauss(key_type *key, int m, double *res);


// nugget
typedef struct nugget_storage {
  double sqrtnugget;
  bool simple, simugrid;
  int *pos, reduced_dim[MAXDIM], prod_dim[MAXDIM + 1];
  double diag[MAXDIM], *red_field;
  
} nugget_storage;
extern double NUGGET_TOL;
int init_nugget(key_type * key, int m);
void do_nugget(key_type *key, int m, double *res );


// see RFbool.cc		
typedef struct mpp_storage {
  // avoid the same names as in addBool_storage, RFbool.h
  double integral, integralsq, // add specific
    integralpos, factor, maxheight,  // extremes specific
    effectiveRadius, effectivearea,
    addradius, min[MAXDIM], length[MAXDIM],
    c[6]; /* local memory for diverse constants,
	     depending on the specific model */
    MPPRandom MppFct;
    int dim;
} mpp_storage;
int init_mpp(key_type * key, int m);
void do_addmpp(key_type *key, int m, double *res );
extern double MPP_APPROXZERO;


// dummy version, hyperplane
typedef struct hyper_storage{
  double rx[MAXDIM], center[MAXDIM], radius;
  hyper_pp_fct hyperplane;
} hyper_storage;
int init_hyperplane(key_type *key, int m);
void do_hyperplane(key_type *key, int m, double *res);


// dummy version, Special Fcts
int init_special(key_type *key, int m); 
void do_special(key_type *key, int m, double *res );


//extremes
typedef struct extremes_storage{
  double inv_mean_pos, assumedmax, *rf;
  key_type key;
} extremes_storage;


// others





///////////////////////////////////////////////////////////////////////
// POSITIONS OF X-VALUES
///////////////////////////////////////////////////////////////////////
#define XSTART 0
#define XEND   1
#define XSTEP  2
#define XSTARTDIM1 XSTART 
#define XENDDIM1 XEND 
#define XSTEPDIM1 XSTEP
#define XSTARTDIM2 XSTART+3
#define XENDDIM2 XEND+3
#define XSTEPDIM2 XSTEP+3
#define XSTARTDIM3 XSTART+6
#define XENDDIM3 XEND+6
#define XSTEPDIM3 XSTEP+6
#define XSTARTDIM4 XSTART+9
#define XENDDIM4 XEND+9
#define XSTEPDIM4 XSTEP+9
extern int XSTARTD[MAXDIM], XENDD[MAXDIM], XSTEPD[MAXDIM];
int Transform2NoGrid(key_type *key, int v);
int Transform2NoGrid(key_type *key, aniso_type aniso, int truetimespacedim,
		     bool simugrid, double **x);
void GetCornersOfElement(double *x[MAXDIM], int timespacedim,
			 covinfo_type *cov, double *sxx);
void GetCornersOfGrid(key_type *key, int Stimespacedim, double *aniso, 
		      double *sxx);
void GetRangeCornerDistances(key_type *key, double *sxx, int Stimespacedim,
			  int Ssimuspatialdim, double *min, double *max);
void GetCenterAndDiameter(key_type *key, bool simugrid, int simuspatialdim, 
			  int truetimespacedim, double *x, aniso_type aniso,
			  double *center, double *lx, double *diameter);


///////////////////////////////////////////////////////////////////////
// TREND SPECIFICATION
///////////////////////////////////////////////////////////////////////
#define TREND_MEAN 0
#define TREND_LINEAR 1
#define TREND_FCT 2
#define TREND_PARAM_FCT 3
void DeleteKeyTrend(key_type *key);



#endif


 
