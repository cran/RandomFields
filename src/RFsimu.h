#ifndef RFsimu_H
#define RFsimu_H 1
// #define RFBOOL 1

#include "RFsimu.public.h"
#include "GSLvsR.h"
#include <string.h>

// global definitions used:
// RF_FLOAT, RF_GSL


#ifdef RF_FLOAT
// we use Real instead of double or float, to get
// flexibility (double->precision, float:reduced memory allocation)
// R needs compilation with Real==double
// tolerances are defined correpondingly
#define EPSILON    0.0000001
#define EPSILON1000    0.0001

#else
#define EPSILON    0.00000000001
#define EPSILON1000    0.00000001
#endif


#define METHODMAXCHAR 31 /* max number of character to describe a Method (including \0)*/
#define COVMAXCHAR 14 /* max number of characters for covariance name*/
#define MAXDIM 3    /* if explicite simulation method is not given by
		       user, the programme will first try a methods 
		       according to "first"; 
		1       then it tries the other methods according to
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

#define TWOPI 6.2831853071795864769252867666
#define INVSQRTTWO 0.70710678118654752440084436210
#define INVSQRTTWOPI 0.39894228040143270286
#define SQRTTWOPI 2.5066282746310002416
#define SQRTPI 1.7724538509055158819
#define SQRTTWO 1.414213562373095048801688724209698
#define INVLOG2 1.442695040888963

//#define TWOSQRTPI 3.5449077018110317638
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
#define ERRORNN 23 /* nn in TBM too large */

#define ERRORDUMMY 99            /* no error, only for tracing purposes */

#define ERRORSILLNULL 101        /* extremes : sill must be positive */
#define ERRORVARMEAN 201          /* Poisson distr.: mean<>variance */

/* parameters: */
#define MEAN       0
#define VARIANCE   1   /* without NUGGET */
#define NUGGET     2
#define SCALE      3 
#define KAPPA      4   /* Kappas should be always the last ones of transferred 
			  parameters */ 
#define KAPPAX     5  
#define KAPPAXX    6
#define LASTKAPPA  7
#define INVSCALE   8   /* 1.0 / SCALE, calculated by InitSimulateRF */
#define SILL       9   /* VARIANCE + NUGGET */
#define TOTAL_PARAM 10 /* check CheckCovariance, InitSimulateRF (both C & R), 
			  and, hyperbolic, TBM3hyperbolic, twodimfractalbrownian 
			  R if changed! */

/* Natural Scaling -- in contrast to the RFparameters(PracticalRange)
   GENERAL_NATURALSCALING may take the following values (furthermore,
   0: no rescaling, +10 numerical determination allowed
*/
#define NATSCALE_EXACT 1
#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check mleRF when changing !! */

/* 
   check with InitSimulateRF in case of any changes !!
   both InitSimulateRF and DoSimulateRF
*/
#define DISTR_GAUSS 0   
#define DISTR_POISSON 1
#define DISTR_MAXSTABLE 2

#define XSTART 0
#define XEND   1
#define XSTEP  2

extern int currentNrCov;
extern char ERRORSTRING_OK[MAXERRORSTRING],ERRORSTRING_WRONG[MAXERRORSTRING];

typedef enum InversionMethod { Cholesky, SVD, NoFurtherInversionMethod } 
InversionMethod;


extern char METHODNAMES[][METHODMAXCHAR];
typedef enum SimulationType {
  CircEmbed,     /* Circulant embedding */
  CircEmbedLocal,/* Circulant embedding on an enlarged grid; (fract. Brown.), 
		    not implemented yet */
  TBM2, TBM3,    /* Turning Bands, performed in 2 or 3 dimensions */
  SpectralTBM,   /* Spectral turning bands, only 2 dimensional */
  Direct,        /* directly, by matrix inversion */
  Nugget,        /* just simulate independent variables */
  AdditiveMpp, /* "additive MPP (random coins)", will be only in the plane*/
  Hyperplane,    /* by superposition of Poisson hyperplane tessellations;
		    only in the plane! Not implemented yet*/
  Special,       /* Not implemented yet, any other (exoctic) method not 
		    included in the above ones */
  Nothing,       /* must always be the last of the list of methods for
		    simulation Gaussian random fields
		 */
  MaxMpp,
  Forbidden      /* must always be the last one */
} SimulationType; 

typedef enum action_type{Store,Storing,UseOldOne,NoAction} action_type;


typedef struct cov_fct;
// in the following structure, user specifications and the intermediate results 
// by the specific init_procedures is stored

typedef void (*destructor)(void **);
typedef struct key_type {
  bool grid,  /* simulation on a grid ? */
    active,   /* has the init_procedure been called successfully? */
    nuggetincluded; /* is the method (see below) able to simulate an additive 
		       nugget effect? */
  SimulationType method, /* the current method (out of SimulationsType) which 
			    is tried at the moment or which has been 
			    successfully initialized */
    tbm_method;          /* if method==tbm2/3, then the method for the line is 
			    specified, this method can only be chosen by the 
			    programmer of the covariance functions, not by the 
			    user (except by global settings)!
			    maybe, a search algorithm among all 1-dim 
			    simulation methods should be implemented, in case 
			    the programmers specification fails
			 */
  Real *x[MAXDIM],       /* the coordinates */
    param[TOTAL_PARAM];  /* parameters */
  void *S;               /* pointer to intermediate results;
			    delicate : either ==NULL or should be a genuine 
			               pointer !
			    no hanging pointers allowed!
			 */ 
  destructor destruct;  /* function called to delete intermediate results */
  void *SX;             /* pointer to intermediate results, second type,
			   like extremes -- not used yet */
  destructor destructX;  
  cov_fct *cov;          /* pointer to CovList[covnr], where details about the
			    covariance functions can be found */
  int covnr,             /* number of the covariance function in CovList; 
			    obtained by function GetCovFunction */
    lx,                  /* number of elements in *x[0][] */
    length[MAXDIM], /* only if grid==true: number of pixels of each 
				direction */
    dim;                 /* dimension of RF */
  long totallength;      /* number of pixels in total;
			    == lx if grid==false
			    == lengthx * lengthy * lengthz if grid==true
			 */
  int naturalscaling;
  int distribution;
} key_type;
extern key_type KEY[MAXKEYS]; 

typedef struct mpp_storage{
  // avoid the same names as in addBool_storage, RFbool.h
  Real integral, integralsq, integralpos,
    effectiveRadius, effectivearea, 
    maxheight, addradius;
  Real min[MAXDIM], length[MAXDIM];
  Real c[10]; /* local memory for diverse constants,
		 depending on the specific model */
  int dim;
} mpp_storage;

// how a covariance function, a spectral measure, a function specifying the
// additive MPP model  should look like
typedef Real (*covfct)(Real,Real *); /* h,parameters */
typedef Real (*scalefct)(Real *, int); /* h,parameters */
typedef int (*checkfct)(key_type*); /* h,parameters */
typedef Real (*parameterfct)(Real *,Real,Real *); /* ?? h,parameters */
typedef Real (*randommeasure)(Real * END_WITH_GSL);     /* parameters, RANDOM */

typedef Real (*mppmodel)(Real *); /* point of random field */
typedef void (*MPPScales)(key_type *); 
/* integral, integralOFsq, effectiveRadius */
typedef void (*MPPRandom)(mpp_storage *bs,
				    Real *min, Real *max,
				    mppmodel*  /* the random function/model,
						    returned */
				    END_WITH_RANDOM);
// how any simulation method should look like, they return the error value
typedef int (*generalSimuInit)(key_type*);   
typedef void (*generalSimuMethod)(key_type*,Real* END_WITH_GSL); /* key, RF, 
								    RANDOM */
// specification for simulation method on line for TBM2/3:
typedef void (*simu1dim)(key_type*,Real* END_WITH_GSL); /* h,parameters */

// here all the different method for simulating a RF of a specified covariance
// function are stored
typedef struct cov_fct {  
  char name[COVMAXCHAR];
  covfct cov;                             /* CircEmbed, direct */
  scalefct naturalscale;
  /* CircEmbedLocal, Brownian motion!: cov_loc : local covariance function
     Real extension_factor(Real *p, Real x, Real *newp): 
             returns the range of the local covariance 
             function (if true one has range one)
             furthermore, it returns a new parameter set, 
	     that is adapted to the local simulation
	     x is usually the largest distance to be simulated
	     square_factor    : see Stein, p.3 (1), parameter $c_2$.
   */
  covfct cov_loc; parameterfct  extension_factor;scalefct square_factor;  
  covfct cov_tbm2,cov_tbm3; SimulationType tbm_method;
  randommeasure spectral;       /* spectral tbm */
  MPPScales add_mpp_scl;
  MPPRandom add_mpp_rnd; 
  void *hyperplane;             /* not programmed yet */
  generalSimuInit initother;    /* any other exotic way */
  generalSimuMethod other;      /* which may be given by a programmer */
  SimulationType first[MAXDIM]; /* the preferred method for a specific dimension 
				   no distinction between grid/no grid -- lack! 
				*/
  int kappas; /* number of parameters additional to the standard ones, 
		 MEAN..SCALE */
  checkfct check;
  //
} cov_fct;


void ErrorMessage(SimulationType m, int error);
extern void InitModelList();   
/* absolutely necessary to call this function at the very beginning !
   but done automatically if SimulateRF, etc. is called for the first time
   (direct call necessary if further covariance functions ate to be added by user!)
   initiating CovList with a dozen standard covariance models
*/

// the following function enables the programmer to get a new covariance 
// function stored
extern void IncludeModel(char *name, 
			  covfct cov, 
			  scalefct naturalscale, 
			  covfct cov_loc, 
			  parameterfct extension_factor,
			  scalefct square_factor,
			  covfct cov_tbm2, 
			  covfct cov_tbm3, 
			  SimulationType tbm_method,
			  randommeasure spectral, 
			  MPPScales add_mpp_scl, 
			  MPPRandom add_mpp_rnd, 
			  generalSimuInit initother,
			  generalSimuMethod other,
			  int kappas,
			  SimulationType r1,SimulationType r2,SimulationType r3,
			  checkfct check
			  ); 
extern void IncludeModel(char *name, covfct cov,scalefct naturalscale,
			  covfct cov_loc, Real extension_factor,
			  covfct cov_tbm2, covfct cov_tbm3, 
			  SimulationType tbm_method, 
			  randommeasure spectral, 
			  int kappas,
			  SimulationType r1,SimulationType r2,SimulationType r3,
			  checkfct check
			  );
extern void IncludeModel(char *name, covfct cov, scalefct naturalscale,
			  covfct cov_tbm2, covfct cov_tbm3, 
			  SimulationType tbm_method,
			  randommeasure spectral, 
			  int kappas,
			  SimulationType r1,SimulationType r2,SimulationType r3,
			  checkfct check
			  );

			 
// all simulation methods have the following structure
// int init_foo(key_type *key): 
//      returns errorcode
//      This function does the most initialization settings,
//      often time consuming and errors may occur
//      intermediate results are stored by means the pointer key->S
//      the function expects that key->S NULL when called
//      the function is expected to return  key->S on error
//      
// void do_foo(key_type *key, Real *res END_WITH_GSLRNG) 
//      returns the simulation result in res ;
//      only basics checks like assert(RANDOM!=NULL);
//      any other error may not occur, as this should be put into init_foo
//      key_active is set to GENERAL_STORING; key->S is freed and set to NULL 
//      if !key_active

// see RFspectral.cc			   
int init_simulatespectral(key_type *key);
void do_simulatespectral(key_type *key, Real *res END_WITH_GSLRNG);

// see RFdirect.cc
int init_directGauss(key_type *key);
void do_directGauss(key_type *key, Real *res END_WITH_GSLRNG);

// see RFcircembed.cc			  
typedef struct FFT_storage { double* work; int *iwork, nseg; } FFT_storage;
int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *S);
void FFT_destruct(FFT_storage *S);
void FFT_NULL(FFT_storage *S);
int internal_init_circulantembedding(key_type * key);
int init_circulantembedding(key_type * key);
void do_circulantembedding(key_type *key, Real *res END_WITH_GSLRNG);
int init_circ_embed_local(key_type * key);
void do_circ_embed_local(key_type *key, Real *res END_WITH_GSLRNG);


// see RFtbm.cc			  
int init_circulantembedding1dim(key_type *key, Real *p);
void do_circulantembedding1dim(key_type *key,Real *res END_WITH_GSLRNG);
int init_turningbands(key_type *key);
void do_turningbands(key_type *key,Real *res END_WITH_GSLRNG);
extern SimulationType TBM_METHOD;


// see RFbool.cc		
int init_mpp(key_type * key);
void do_addmpp(key_type *key, Real *res END_WITH_GSLRNG);


void printKEY(int keyNr);
extern int GENERAL_PRINTLEVEL;
extern int GENERAL_NATURALSCALING;
extern int GENERAL_STORING;
extern char ERRORSTRING_OK[MAXERRORSTRING];
extern char ERRORSTRING_WRONG[MAXERRORSTRING];

extern Real MPP_APPROXZERO;

#endif











