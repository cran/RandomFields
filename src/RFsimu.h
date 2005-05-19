
#ifndef RFsimu_H
#define RFsimu_H 1

#include "auxiliary.h"
#include "RandomFields.h"
#include <string.h>

// global definitions used:
// RF_FLOAT, 


#define EPSILON     0.00000000001
#define EPSILON1000 0.000000001

#define DISTRMAXCHAR 10
#define DISTRNR 3
#define METHODMAXCHAR 31 /* max number of character to describe a Method (including \0)*/
#define COVMAXCHAR 14 /* max number of characters for covariance name*/
#define MAXCOV 10     /* max number of covariance functions for a random 
			 field -- NEVER MORE THAN 255 */
#define MAXDIM 4          
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
#define MAXERRORSTRING 100
#define MAXKEYS 10
#define MAXNRCOVFCTS 30
#define GetNotPrint 1 /* if parameters are questionned, then also printed if =0 */
// Error codes & messages
// positive values are error codes
// negative values are messages
#define USEOLDINIT -1            /* not an error, just a message !*/
#define NOERROR 0                 
#define ERRORNOTDEFINED 1        /* the specification for the  covariance and 
				    method is not given/known, e.g. TBM2 for 
				    many covariance functions  */
#define ERRORMETHODNOTALLOWED 2  /* wrong method for the specified covariance 
				    function or grid */
#define ERRORNOTPROGRAMMED 3     
#define ERRORCOVNOTALLOWED 4     /* wrong dimension for the specified covariance 
				    function */
#define ERRORFAILED 5            /* method didn't work for the specified 
				    parameters */
#define ERRORMEMORYALLOCATION 6  /* malloc returned NULL pointer */
#define ERRORNOTINITIALIZED 7    /* key.active==false in DoSimulateRF; is only 
				    checked there !!! */
#define ERRORKEYNROUTOFRANGE 8   /* wrong keynr */
#define ERRORDECOMPOSITION 9     /* direct simulation; matrix inversion failed */
#define ERRORPRECISION 10        /* direct simulation; matrix inversion didn't 
				    reach the desired presicion */
#define ERRORRESCALING 11        /* attempt to rescale where not possible */
#define ERRORVARIOGRAMONLY 12    /* attempt to get the covariance whereas only 
				    the variogram is defined */
#define ERRORFOURIER 13          /* specific error in Fourier transformation */
#define ERRORCOVFAILED 14        /* covariance model not defined for specified 
				    parameters */
#define ERRORREGISTER 15         /* wrong register number */
#define ERRORCOORDINATES 16      /* coordinates are not given or invalid grid 
				    specification */
#define ERRORNEGATIVEVAR 17
#define ERRORPARAMNUMBER 18
#define ERRORDIM 19              /* dim<1 or dim>MAXDIM */
#define ERRORNEGATIVESCALE 20
#define ERRORWAVING 21           /* numerical determination of practical range
				    fails */
#define ERRORWRONGINIT 22        /* dosimulate and initsimulate do not match */
#define ERRORNN 23               /* nn in TBM too large */
#define ERRORANISOTROPIC 24      /* trying to call anisotropic fct with 
				    isotropic parameters */
#define ERRORMETHODMIX 25        /* diff. methods for multiplication 
				    covariance */
#define ERRORTIMENOTANISO 26      /* if time is used then, anisotropic must
				    be true */
#define ERRORNOMULTIPLICATION 27 /* method not allowed for multiplicative 
				    covariance functions*/
#define ERRORANISOMIX 28         /* diff. anisotropies for multiplication 
				    covariance and tbm */
#define ERRORISOTROPICMETHOD 29  /* TBM might be called only by isotropic
                                    functions */
#define ERRORTIMESEPARATE 30     /* last line of of matrix 
				  */
#define ERRORUNKNOWNMETHOD 31  

#define ERRORWITHOUTTIME 32     /* spatially isotropic covariance functions
				   must always be called with a time 
				   component */
#define ERRORCOVNROUTOFRANGE 33 /* not allowed covariance number */
#define ERRORWRONGDIM 34        /* (tbm) dimension couldn't be reduced to be 
                                   appropriate for the respective method */
#define ERRORNOTANISO 35        /* make sure that SPACEISOTROPIC and ANISOTROPIC
                                   covariance function are called by anisotropy
				   definition */
#define ERRORDIMMISMATCH 36     /* logical dim does not match physical dim in the
                                   covariance function call */
#define ERRORTOOMANYPOINTS 37   /* two many points to determine smallest 
				   distance between within TBM! */
#define ERRORTIMENOTALLOWED 38  /* time component not allowed for the specified 
                                   method */
#define ERRORCOVNUMERICAL 39    /* covariance function of derivative cannot
				   be calculated exactly, e.g. in TBM2 */
#define ERRORMAXMEMORY 40       /* in circulant embedding: circulant vector too
				   large */
#define ERRORTOOMANYLINES 41    /* Hyperplane tesselation: estimated simulated
				   lines passes given threshold*/
#define ERRORANYLEFT 42         /* some covariances have not been simulated yet
				 */
#define ERRORMETHODEXCLUDED 43  /* method explicitely excluded by other 
				   parameter; here: stationary_only=TRUE */
#define ERROROUTOFMETHODLIST 44 /* no further method we can try is available */
				 
#define NOERROR_REPEAT 50
#define NOERROR_ENDOFLIST 51

#define ERRORDUMMY 99            /* no error, only for tracing purposes */

#define ERRORSILLNULL 101        /* extremes : sill must be positive */
#define ERRORPAIRS 102          /* extremes : paired not allowed */
#define ERRORTREND 103           /* extremes: no trend allowed */

#define ERRORVARMEAN 201         /* Poisson distr.: mean<>variance */
#define ERRORPOKETBMdim 300     /* mixing of two initialisations for TBM not 
				    allowed: dimensions differ */
#define ERRORPOKETBMmeth 301     /* methods differ */
#define ERRORPOKETBMparam 302     /* model parameters differ -- 
				     not extensively tested */

/* do not use numbers 800 -- 900 : reserved to MPP package */

#define ERRORTRIVIAL 998        /* tbm does not allow that the matrix `a' of 
                                   a product matrix are all zero
				 */ 
#define ERRORUNSPECIFIED 999   
/* parameters: */
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
#define SCALE      8
#define INVSCALE   9 /* always SCALE + 1 */
#define ANISO      8
#define TIMEINVSCALE 8 /* used in tbm3, where invscale and timeinvscale 
			  are used */
#define DIAMETER ANISO + MAXDIM * MAXDIM /* usually the diameter of the field
					 in scale units -- used in local
					 circulant embedding */
#define TBM2NUM DIAMETER /* numerical evaluation for TBM2 ? */
#define TOTAL_PARAM 1 + LASTKAPPA + MAXDIM * MAXDIM + 1 /* 25 */
/*  LASTKAPPA + MAXDIM^2 + 1
    check CheckCovariance, InitSimulateRF (both C & R), 
    and, hyperbolic, TBM3hyperbolic, twodimfractalbrownian 
    R if changed! */


// parameters used in local cov fcts
#define LOCAL_R 0

#define CUTOFF_R LOCAL_R
#define CUTOFF_A 1
#define CUTOFF_B 2

#define INTRINSIC_R LOCAL_R
#define INTRINSIC_A0 1
#define INTRINSIC_A2 2
#define INTRINSIC_B 3

#define TOTAL_LOCALPARAM 3+1 //see INTRINSIC_B

typedef enum local_strategy_type
{
  // Do NOT change the ordering as they have relative priorities used in overall_strategy
  // in RFcircembed.cc

  // TheoGuaranteed, JustTry and UserSpecified are used in cutoff
  // TheoGuaranteed, NumeGuaranteed, SearchR, and UserSpecified are used in intrinsic

  TheoGuaranteed,       /*  parameters are 
                            within the range of analytical results, so the
                            covariance matrix (with the specified R(intrinsic) or a(cutoff)) 
                            is ensured to be nonnegative definite */
  JustTry,              /*  used in cutoff, when alpha>1 */
  NumeGuaranteed,        /*  parameters are 
                            within the range of numerical results, so the
                            covariance matrix (with the specified R) 
                            is believed to be nonnegative definite */
  SearchR,               /* parameters are 
                            outside the range of both analytical and numerical results,
                            so we'll try several possible values of R */
  UserSpecified,

  // below are only for communications between functions
  // like in the case of alpha<= 0.5 in cutoff

  TellMeTheStrategy,
  CallMeAgain,
  IncreaseCutoffA
} local_strategy_type;


typedef enum EmbedType
{
  Standard, 
  Cutoff, 
  Intrinsic
} EmbedType;

typedef local_strategy_type (*IntrinsicStrategyFct)(double *, double, int, 
						    double *);
typedef local_strategy_type (*CutoffStrategyFct)(double *, double*, 
						 local_strategy_type);

#define XSTART 0
#define XEND   1
#define XSTEP  2

#define ANISOTROPIC 0
#define SPACEISOTROPIC 1
#define FULLISOTROPIC 2

// way of implementing simulation methods:
#define NOT_IMPLEMENTED 0 /* do not change this value except to false */
#define IMPLEMENTED 1     /* do not change this value except to true */
#define NUM_APPROX 2 
#define GIVEN_METH_IGNORED -1

extern int currentNrCov;
extern char ERRORSTRING_OK[MAXERRORSTRING],ERRORSTRING_WRONG[MAXERRORSTRING];
typedef struct cov_fct;
extern cov_fct *CovList;

typedef enum InversionMethod { Cholesky, SVD, NoFurtherInversionMethod } 
InversionMethod;

typedef double param_type[MAXCOV][TOTAL_PARAM];
typedef double local_param_type[MAXCOV][TOTAL_LOCALPARAM]; 

extern char METHODNAMES[][METHODMAXCHAR];
extern char DISTRNAMES[][DISTRMAXCHAR];
extern char OP_SIGN[][2];
typedef enum SimulationType {
  /* never change the ordering -- at least check with "compatible" */
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
  Forbidden      /* must always be the last one */
} SimulationType; 
// #define SimulationTypeLength 13
#define SimulationTypeLength ((int) Forbidden) - 1

typedef enum action_type{Store,Storing,UseOldOne,NoAction} action_type;

#define TREND_MEAN 0
#define TREND_LINEAR 1
#define TREND_FCT 2
#define TREND_PARAM_FCT 3

// in the following structure, user specifications and the intermediate results 
// by the specific init_procedures is stored

typedef void (*destructor)(void **);
typedef struct key_type {
  bool 
  storing, /* value of GENERAL_STORING at initialisation time */
    grid,  /* simulation on a grid ? */
    active,   /* has the init_procedure been called successfully? */
    anisotropy; /* is it an isotropic covariance function ? */
  //   nuggetincluded; /* is the method (see below) able to simulate an additive
  //		       nugget effect? */
  SimulationType method[MAXCOV], 
                         /* the current method (out of SimulationsType) which 
			    is tried at the moment or which has been 
			    successfully initialized */
    tbm_method;  /* if method==tbm2/3, then the method for the line is 
			    specified, this method can only be chosen by the 
			    programmer of the covariance functions, not by the 
			    user (except by global settings)!
			    maybe, a search algorithm among all 1-dim 
			    simulation methods should be implemented, in case 
			    the programmers specification fails
			 */
  double *x[MAXDIM];       /* the coordinates */
  param_type param;      /* parameters */
  int totalparam;
  void *S[MAXCOV];       /* pointer to intermediate results;
			    delicate : either ==NULL or should be a genuine 
			               pointer !
			    no hanging pointers allowed!
			 */ 
  destructor destruct[MAXCOV];  
                        /* function called to delete intermediate results */
  void *SExtremes;      /* pointer to intermediate results, second type,
			   like extremes -- not used yet */
  destructor destructExtremes;  
  int op[MAXCOV],
    covnr[MAXCOV],     /* number of the covariance function in CovList; 
			    obtained by function GetCovFunction */
    lx,                  /* number of elements in *x[0][] */
    length[MAXDIM], /* only if grid==true: number of pixels of each 
				direction */
    spatialdim;                 /* spatial dimension of RF */
  long spatialtotalpoints,
    totalpoints;      /* number of pixels in total;
			    == lx if grid==false
			    == lengthx * lengthy * lengthz if grid==true
			 */
  //  int naturalscaling;
  int distribution,
    ncov,             /* number of covariance function actually used */
    timespacedim;      /* total dimension, including time and space */
  bool Time,             /* is a time component given ? */
    left[MAXCOV];        /* left to be considered in the simulation
			    initialisations ! */
  double T[3];             /* gridtriple definition for the time components */
  // int timelength;        /* number of points in the time direction */
  bool compatible;     // determines whether simulation result must be stored
  //                   in an intermediate storage since partial results cannot
  //                   be added directly
  SimulationType unimeth[MAXCOV];
  bool unimeth_alreadyInUse[MAXCOV];
  unsigned char n_unimeth;
  double mean;         // only used if storagemodus=0 (see rf.R)
  double *LinearTrend; // only used for storing -- used in R
  int lLinTrend;              
  char *TrendFunction; // only used for storing -- used in R
  int lTrendFct;
  int TrendModus;      // dito
} key_type;

// how a covariance function, a spectral measure, a function specifying the
// additive MPP model  should look like
typedef void (*infofct)(double *, int *, int *); /* param, maxdim, 
      CEbadlybehaved? --  to create preference list */
typedef void (*rangefct)(int, int *, double*);
/* DO NOT change double to anything else since memcpy  fcts are used
   dim, class parameter like in nsst, 
   return 2 x length(param) x {theor, pract }*/
typedef double (*covfct)(double *, double* , int); /* h,parameters */
typedef double (*isofct)(double*, double*, int); /* h,parameters */
typedef double (*scalefct)(double *); /* parameters, Brownian motion  */
typedef double (*natscalefct)(double *, int); /* parameters, ; natural scaling */
typedef int (*checkfct)(double*, int, SimulationType); /* h,parameters; covariance fct, nr. */
typedef int (*initfct)(key_type*, int); /* h,parameters */
typedef double (*CovFctType)(double*, int, int*, int*, param_type, int, bool);
typedef void (*do_comp_fct)(key_type*, bool, int, double*); /* h,parameters */
typedef void (*do_incomp_fct)(key_type*, int, double*); /* h,parameters */
typedef double (*randommeasure)(double *);     /* parameters, RANDOM */

typedef double (*mppmodel)(double *); /* point of random field */
/* integral, integralOFsq, effectiveRadius */
typedef struct mpp_storage;
typedef void (*MPPRandom)(mpp_storage *s,
			  int v,
			  double *min, double *max,
			  mppmodel*  /* the random function/model,
					returned */
                          );
typedef void (*initmppfct)(mpp_storage*, int v); /* h,parameters */
// how any simulation method should look like, they return the error value
typedef int (*hyper_pp_fct)(double*, double*, int, bool, double**, double**, 
			    double**);

typedef struct mpp_storage{
  // avoid the same names as in addBool_storage, RFbool.h
  double integral[MAXCOV], integralsq[MAXCOV], integralpos[MAXCOV],
    effectiveRadius[MAXCOV], effectivearea[MAXCOV], 
    maxheight[MAXCOV],
    addradius[MAXCOV];
  double min[MAXDIM], length[MAXDIM],
    c[MAXCOV][6]; /* local memory for diverse constants,
	      depending on the specific model */
  param_type param;
  MPPRandom MppFct[MAXCOV];
  int actcov,simuspatialdim, timespacedim;
  bool grid;
  double *x;
} mpp_storage;

#define TRIVIALSTARTEGY 0
#define STRATEGY_PARTIAL 1
#define LASTSTRATEGY STRATEGY_PARTIAL
typedef struct ce_param{
  bool force, userfft;
  char strategy;
  double tol_re, tol_im;
  int trials;
  double maxmem;
  int mmin[MAXDIM];
} ce_param;

#define DECISION_TRUE 1
#define DECISION_FALSE 0
#define DECISION_CASESPEC -1
typedef struct decision_param{
  int stationary_only, exactness;
} decision_param;

typedef struct local_user_param{
  double cutoff_a, intrinsic_r;
} local_user_param;

typedef int (*generalSimuInit)(key_type*);   
typedef void (*generalSimuMethod)(key_type*,double* ); /* key, RF, 
								    RANDOM */
// specification for simulation method on line for TBM2/3:
typedef void (*simu1dim)(key_type*, int m, double* ); /* h,parameters */

// here all the different method for simulating a RF of a specified covariance
// function are stored
typedef struct cov_fct {  
  char name[COVMAXCHAR];
  covfct cov;                             /* CircEmbed, direct */
  natscalefct naturalscale;
  int implemented[SimulationTypeLength];
  // since other methods make use of other information to determine whether they 
  // are implemented for the cov model, this bool vector is currently only used by 
  // local circ embeds; but it also can be used by future methods
  IntrinsicStrategyFct intrinsic_strategy;
  CutoffStrategyFct cutoff_strategy;

  isofct cov_tbm2, cov_tbm3, derivative, secondderivt; 
  SimulationType tbm_method;
  randommeasure spectral;       /* spectral tbm */
  initmppfct add_mpp_scl;
  MPPRandom add_mpp_rnd; 
  hyper_pp_fct hyperplane;      /* hyperplane tessellations */        
  generalSimuInit initother;    /* any other exotic way */
  generalSimuMethod other;      /* which may be given by a programmer */
  infofct info;
  /// SimulationType first[MAXDIM]; 
  /* OBSOLETE -- replaced by method
     the preferred method for a  specific dimension 
     no distinction between grid/no grid -- lack! 
  */
  rangefct range;
  int kappas; /* number of parameters additional to the standard ones, 
		 MEAN..SCALE */
  checkfct check;
  char isotropic;
  bool variogram, even, odd[MAXDIM];
  //
} cov_fct;


// the following function enables the programmer to get a new covariance 
// function stored
extern int IncludeModel(char *name, int kappas, checkfct check,
			int isotropic, bool variogram, infofct info,
			rangefct range);	
extern void addSimu(int nr, 
		    SimulationType r1,SimulationType r2,SimulationType r3);
extern void addCov(int nr, covfct cov, isofct derivative,
		   natscalefct naturalscale);
extern void addLocal(int nr, IntrinsicStrategyFct intrinsic_strategy, 
		     CutoffStrategyFct cutoff_strategy, isofct secondderivt);
extern void addTBM(int nr, isofct cov_tbm2, isofct cov_tbm3,
		   SimulationType tbm_method, randommeasure spectral);
extern void addOther(int nr, initmppfct add_mpp_scl, MPPRandom add_mpp_rnd, 
		     hyper_pp_fct hyper_pp,
		     generalSimuInit initother, generalSimuMethod other);
extern void addodd(int nr, int dim); 
/* function necessary to set CircEmbed correctly if the function is not 
   even in some coordinates! 
*/


int FirstCheck_Cov(key_type *key, SimulationType Method, param_type param, 
		   bool No_Multiply, int *covnr, int *multiply, 
		   unsigned short int *actcov);

void ErrorMessage(SimulationType m, int error);
void printkey(key_type *key);
extern void InitModelList();   
/* absolutely necessary to call this function at the very beginning !
   but done automatically if SimulateRF, etc. is called for the first time
   (direct call necessary if further covariance functions ate to be added by user!)
   initiating CovList with a dozen standard covariance models
*/


double CovFct(double *x, int dim, int *covnr, int *op, param_type param, 
	    int ncov, bool anisotropy);
void GetTrueDim(bool anisotropy, int timespacedim, double* param, 
	        int *TrueDim, bool *no_last_comp,
		int *start_param, int *index_dim);
int Transform2NoGrid(key_type *key, double* param,  int TrueDim, 
		     int *start_param, double **xx);
void GetCornersOfGrid(key_type  *key, int Stimespacedim, int* start_param, 
		      double* param, double *sxx);
void GetCornersOfElement(key_type  *key, int Stimespacedim, int* start_param, 
		      double* param, double *sxx);
void GetRangeCornerDistances(key_type *key, double *sxx, int Stimespacedim,
			  int Ssimuspatialdim, double *min, double *max);
extern int eigenvalues(double *C, int dim, double *ev);
// all simulation methods have the following structure
// int init_foo(key_type *key): 
//      returns errorcode
//      This function does the most initialization settings,
//      often time consuming and errors may occur
//      intermediate results are stored by means the pointer key->S
//      the function expects that key->S NULL when called
//      the function is expected to return  key->S on error
//      
// void do_foo(key_type *key, double *res ) 
//      returns the simulation result in res ;
//      only basics checks like assert(RANDOM!=NULL);
//      any other error may not occur, as this should be put into init_foo
//      key_active is set to GENERAL_STORING; key->S is freed and set to NULL 
//      if !key_active


// see RFspectral.cc			   
int init_simulatespectral(key_type *key, int m);
void do_simulatespectral(key_type *key, int m, double *res );


// see RFdirect.cc
typedef struct direct_storage {
  InversionMethod method;
  double *U,*G;
} direct_storage;
extern int DIRECTGAUSS_BESTVARIABLES;
extern int DIRECTGAUSS_MAXVARIABLES;
int init_directGauss(key_type *key, int m);
int internal_init_directGauss(direct_storage **S, bool grid,
			      int spatialdim, bool Time, double** x, int* length,
			      int totalpoints, double *T, CovFctType CovFct,
			      int *covnr, int *op, param_type param,
			      int actcov, bool anisotropy);
void do_directGauss(key_type *key,  bool add, int m, double *res);
void internal_do_directGauss(direct_storage *S, bool add, long totpnts,
			     double *res);
void direct_destruct(void ** S);


// see RFcircembed.cc
typedef struct FFT_storage { double* work; int *iwork, nseg; } FFT_storage;
unsigned long NiceFFTNumber(unsigned long nn);
int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *S);
int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *S);
void FFT_destruct(FFT_storage *S);
void FFT_NULL(FFT_storage *S);
void SetParamCE( int *action, int *force, double *tolRe, double *tolIm,
		 int *trials, 
		 int *mmin, /* vector of length MAXDIM */
		 int *userfft, int *strategy, double *maxmem,
		 ce_param *CE, char *name);
int internal_init_circ_embed(double *steps, bool anisotropy, 
			     int *covnr, int *op,
			     param_type param, 
			     int *nn, int *m, int *cumm, int *halfm, 
			     int dim, int actcov,
			     CovFctType CovFct,
			     ce_param *cepar,
			     FFT_storage *FFT,
			     long *twoRealmtot,
			     double **cc); //!
void internal_do_circ_embed(int *nn, int *m, int *cumm, int *halfm,
			    double *c, double *d, int Ntot, int dim,
			    FFT_storage *FFT_heap, bool add,  double *res );
int init_circ_embed(key_type * key, int m);
void do_circ_embed(key_type *key, bool add, int m, double *res );
int init_circ_embed_intrinsic(key_type * key, int m);
void do_circ_embed_intrinsic(key_type *key, bool add, int m, double *res );
int init_circ_embed_cutoff(key_type * key, int m);


// see RFtbm.cc			  
int init_turningbands2(key_type *key, int m);
int init_turningbands3(key_type *key, int m);
void do_turningbands(key_type *key, int m, double *res );
extern SimulationType TBM_METHOD;


// see RFbool.cc		
int init_mpp(key_type * key, int m);
void do_addmpp(key_type *key, int m, double *res );

// dummy version, hyperplane
int init_hyperplane(key_type *key, int m);
void do_hyperplane(key_type *key, int m, double *res);

// dummy version, Special Fcts
int init_special(key_type *key, int m); 
void do_special(key_type *key, int m, double *res );

// nugget
int init_nugget(key_type * key, int m);
void do_nugget(key_type *key, bool add, int m, double *res );

extern int GENERAL_PRINTLEVEL;
extern int GENERAL_NATURALSCALING;
extern int GENERAL_STORING;
extern double GENERAL_PRECISION;
extern decision_param DECISION_PARAM;
extern char ERRORSTRING_OK[MAXERRORSTRING];
extern char ERRORSTRING_WRONG[MAXERRORSTRING];
extern char GENERAL_PCH[2]; 


extern double MPP_APPROXZERO;
extern key_type KEY[MAXKEYS]; 
extern double ZERO[MAXDIM];
extern double UNIT[MAXDIM];
extern initfct init_method[SimulationTypeLength];
extern do_comp_fct do_compatiblemethod[SimulationTypeLength];
extern do_incomp_fct  do_incompatiblemethod[SimulationTypeLength];


#define SET_DESTRUCT(A)\
  assert(key->active==false);\
  assert((key->S[m]==NULL) && (key->destruct[m]==NULL));\
  key->destruct[m] = A;
void DeleteKeyTrend(key_type *key);


#endif


 
