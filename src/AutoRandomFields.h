#ifndef AutoRandomFields_H
#define AutoRandomFields_H 1
#include <R.h>

#define MAXCHAR_RF MAXCHAR // max number of characters for (covariance) names  
#define METHODMAXCHAR MAXCHAR // max number of character to describe a Method (including \0)
#define MAXPARAM 20
#define MAXSUB 10

#define MAXCEDIM 13
#define MAXTBMSPDIM 4
#define MAXMPPDIM 4
#define MAXMPPVDIM 10
#define MAXVDIM MAXMPPVDIM
#define MAXHYPERDIM 4
#define MAXVARIODIM 20
#define MAXTBMVDIM 5
#define MAXGAUSSVDIM 10
#define MAXACCEPTED 1e300
#define MAXVARIANTS 6
#define MAXSYSTEMS 1
#define MAXDATANAMES 5

#define MAXBOXCOXVDIM 10

#define MAXCOOORDNAMES 4
typedef enum units_enum {
  units_none, units_km, units_miles, units_time, units_user} units_enum;
#define nr_units (units_user + 1)

// see also which_system in rf_interface.cc
typedef enum coord_sys_enum {
  coord_auto, coord_keep, cartesian, earth, sphere, gnomonic, orthographic,
  coord_mix
} coord_sys_enum;
#define nr_coord_sys (coord_mix + 1)

typedef enum reportcoord_modes {reportcoord_always, 
				reportcoord_warnings_orally, 
				reportcoord_warnings, 
				reportcoord_none} reportcoord_modes;
#define nr_reportcoord_modes (reportcoord_none + 1)

typedef enum modes {careless, sloppy, easygoing, normal, precise, 
		    pedantic, neurotic} modes;
#define nr_modes (neurotic + 1)

typedef enum output_modes {output_sp, output_rf, output_geor} output_modes;
#define nr_output_modes (output_geor + 1)

typedef enum angle_modes {radians, degree} angle_modes;
#define last_angle_mode degree

// 0 sicherheitshalber nirgendswo nehmen
// JEGLICHE AENDERUNG AUCH BEI monotone_type DURCHZIEHEN
#define PARAM_DEP -1 // parameter dependent fuer vdim, finiterange, etc; for DefList only, falls von submodel und param abhaengig, so submodel_dep 
#define PREVMODEL_DEP -2 // NIE AENDERN!
#define SUBMODEL_DEP -3  // immer, wenn zumindest submodel depedence
#define MISMATCH -4
#define UNSET -5 // NIE AENDERN! -- vgl. Keyinfo, montone, mpp.moment

// NOTE!!!
// if any definition is changed check setbackward() 
typedef enum domain_type { // IF CHANGED CHANGE ALSO RFgetModelNames
  XONLY, KERNEL, PREVMODEL_D, SUBMODEL_D, PARAMDEP_D, KEEPCOPY_DOM, 
  DOMAIN_MISMATCH // always last ! 
} domain_type;
#define FIRST_DOMAIN XONLY
#define LAST_DOMAINUSER KERNEL
#define LAST_DOMAIN DOMAIN_MISMATCH


// NOTE!!
// !!! CHANGE ALSO RF_GLOBALS.R !!!
typedef enum isotropy_type { 
  ISOTROPIC,  // RC genauer : rotation invariant!!
  DOUBLEISOTROPIC , // RC fully symmetric 
  VECTORISOTROPIC ,
  SYMMETRIC , //"stationary" only if XONLY; in the covariance sense if multivariate!
  CARTESIAN_COORD,  // RC
  GNOMONIC_PROJ ,   // RC 5, projection on the plane
  ORTHOGRAPHIC_PROJ ,  // RC dito; to do: implement further projections
  //
  SPHERICAL_ISOTROPIC , 
  SPHERICAL_SYMMETRIC , 
  SPHERICAL_COORD  , // RC
  //
  EARTH_ISOTROPIC  , // 10
  EARTH_SYMMETRIC  ,
  EARTH_COORD  , // RC
  //
  //  LOGCART_SYMMETRIC,
  //  LOGCART_COORDS,
  CYLINDER_COORD, // unused
  UNREDUCED ,  // full cartesian or earth coord 
  PREVMODEL_I ,  // 15, type taken from submodels
  SUBMODEL_I, // if model is both, submodel_dep and param_dep set submodel_i;
              // setDI wird erst vor domain-setzen aufgerufen
  PARAMDEP_I, // purely paramdep, not submodel_dep
  
  KEEPCOPY_ISO,
  ISO_MISMATCH // always last ! 
} isotropy_type;


#define FIRST_CARTESIAN ISOTROPIC
#define FIRST_ISOUSER FIRST_CARTESIAN 
#define LAST_REDUCEDXDIM_CART DOUBLEISOTROPIC
#define FIRST_PROJECTION GNOMONIC_PROJ
#define LAST_PROJECTION ORTHOGRAPHIC_PROJ // danach MUSS sphaerisch kommen!!
#define LAST_CARTESIAN LAST_PROJECTION // danach MUSS sphaerisch kommen!!
#define FIRST_SPHERICAL SPHERICAL_ISOTROPIC
#define FIRST_ANYSPHERICAL FIRST_SPHERICAL
#define LAST_TRUESPHERICAL SPHERICAL_COORD
#define FIRST_EARTH EARTH_ISOTROPIC
#define LAST_EARTH EARTH_COORD
#define LAST_ANYSPHERICAL LAST_EARTH
//#define FIRST_LOGCART LOGCART_SYMMETRIC
//#define LAST_LOGCART LOGCART_COORDS
#define LAST_ISOUSER UNREDUCED
#define LAST_ISO ISO_MISMATCH

  

typedef enum monotone_type {
  MON_UNSET = UNSET,
  MON_MISMATCH = MISMATCH,
  MON_SUB_DEP = SUBMODEL_DEP,
  MON_PREV_DEP = PREVMODEL_DEP,
  MON_PARAMETER = PARAM_DEP,
  NOT_MONOTONE = 0,
  MONOTONE = 1,
  GNEITING_MON = 2, // Euclid's hat, Gneiting, J.Mult.Anal. 69, 1999
  NORMAL_MIXTURE = 3,
  COMPLETELY_MON = 4,
  BERNSTEIN = 5
} monotone_type; // is different from everything else before
#define MONOTONE_TOTAL (BERNSTEIN - MON_UNSET + 1)

  

#define MAXFIELDS 10
#define MODEL_USER (MAXFIELDS + 0)  // for user call of Covariance etc.
#define MODEL_COV (MAXFIELDS + 1)
#define MODEL_COVMATRIX (MAXFIELDS + 2)
#define MODEL_VARIOGRAM (MAXFIELDS + 3)
#define MODEL_PSEUDO (MAXFIELDS + 4)
#define MODEL_FCTN (MAXFIELDS + 5)
#define MODEL_DISTR (MAXFIELDS + 6)
#define MODEL_CALC (MAXFIELDS + 7)


#define LAST_MODEL_USER (MAXFIELDS + 9)
#define FIRST_INTERNAL (LAST_MODEL_USER + 1)
#define MODEL_AUX  (FIRST_INTERNAL + 0)  // auxiliary in fitgauss.R & plot
#define MODEL_INTERN (FIRST_INTERNAL+1)//for kriging, etc; internal call of cov 
#define MODEL_SPLIT (FIRST_INTERNAL + 2)// split cov model & other aux methods 
#define MODEL_GUI (FIRST_INTERNAL + 3)   // RFgui 
#define MODEL_MLE (FIRST_INTERNAL + 4) // mle covariance model 
#define MODEL_MLESPLIT (FIRST_INTERNAL + 5)  // ="= 
#define MODEL_LSQ (FIRST_INTERNAL + 6)  // used in fitgauss for variogram calc
#define MODEL_BOUNDS (FIRST_INTERNAL + 7)  // MLE, lower, upper 
#define MODEL_KRIGE  (FIRST_INTERNAL + 8)
#define MODEL_PREDICT  (FIRST_INTERNAL + 9) 
#define MODEL_ERR  (FIRST_INTERNAL + 10)
#define MODEL_MAX MODEL_ERR  


typedef enum Types { // IF CHANGED, CHANGE RFgetModelNames!
  // types within model definitions (without processes)
  TcfType, PosDefType, VariogramType, NegDefType,
  PointShapeType, // i.e. shape function for max-stable process + process
  ShapeType, // 5, i.e., no special propeties with structure of covariance fctn
  TrendType, // 
  RandomOrShapeType, // only for parameters, e.g. nu in Whittle
  ManifoldType, // Bedeutung: C->Type:TypeFct existiert; cov->typus:ungesetzt
                 // frame : unknown
  ProcessType, GaussMethodType, // # 10
  NormedProcessType, // currently only used of binary processes
  BrMethodType, // change also rf_globals.R if deleted

  SmithType, // 13
  SchlatherType, 
  PoissonType,
  PoissonGaussType, // 16

  RandomType, // 17
  InterfaceType, 
  MathDefType, // 19
  OtherType, // 20 always last for usual use
  //
  // special types, only for internal use
  BadType,
  SameAsPrevType, // currently unused
  LikelihoodType,
  EvaluationType, //
  //
  // special inputs
  MixedInputType,//
  CharInputType, 
  NN1, NN2, // 
  NN3, NN4 // 30, for use in generateMmodels only
} Types;

/* FRAMES: calls from within the respective processes
  TrendType, // also any other rather arbitrary function
  ProcessType, GaussMethodType, // 
  NormedProcessType, //
  BrMethodType, // 
  RandomType, // with a random variable construction
  InterfaceType, // 
  EvaluationType, // likelihood, RFcov, etc
  SmithType,  // 
  SchlatherType,  // 
  PoissonType,
  PoissonGaussType, /
*/


#define LASTTYPE NN4

#define nOptimiser 8
#define nNLOPTR 15
#define nLikelihood 4
#define nNamesOfNames 12
#define nVAR2COV_METHODS 4
#define nProjections 2

typedef enum sortsofeffect { // ! always compare with convert.R!
  DetTrendEffect, // trend, nichts wird geschaetzt
  FixedTrendEffect, // linearer Parameter wird geschaetzt  
  FixedEffect, // trend is also converted to FixedEffect ?
  DataEffect, // box cox
  RandomEffect, // random efffect : ALWAYS first covariance model part
  ErrorEffect, // error term: orig Remaining...
  effect_error //
} sortsofeffect;

//#define FirstMixedEffect FixedEffect
//#define LastMixedEffect SpVarEffect

//////////////////////////////////////////////////////////////////////
// the different types of parameters
typedef enum sortsofparam { // never change ordering; just add new ones !!!!!!
  VARPARAM,  // 0
  SIGNEDVARPARAM, //1
  SDPARAM, //2
  SIGNEDSDPARAM, // 0..3
  SCALEPARAM, //4
  DIAGPARAM, // 5
  ANISOPARAM, // 4..6
  INTEGERPARAM,//7
  ANYPARAM, //8
  TRENDPARAM, // unused , // 7..9
  NUGGETVAR,
  CRITICALPARAM, 
  DONOTVERIFYPARAM,  // lists arguments and worse; not used in MLE!
  ONLYRETURN, // not recognised in mle, but returned
  FORBIDDENPARAM, // still a possible user input!! just estimation is forbidden

  UNKNOWNPARAM,
  
  VARONLYMLE, // same as VARPARAM, but not returned
  CRITONLYMLE, // same as CRITICAL, but not returned
  ONLYMLE, // not returned  -- never delete this value from the list
  IGNOREPARAM, // neither recognised in MLE nor returned by key

  STANDARD, // for use in R: all up to  LASTUSERSORTOF C
  INCLUDENOTRETURN, // STANDARD + those marked with IGNOREPARAM
  INTERNALPARAMETERS, //only those that have the INTERNAL_PARAM name
  ALLPARAMETERS, //really all parameters
  NOPARAMETERS // no parameter at all
} sortsofparam;
#define LASTRETURNED FORBIDDENPARAM
#define FIRSTONLYMLE VARONLYMLE
#define LASTONLYMLE ONLYMLE
#define LASTUSERSORTOF FORBIDDENPARAM
#define LASTSORTOF NOPARAMETERS


// never change the ordering -- at least chech nich Standard in location_rules
//     in gauss.cc
//
typedef enum Methods {
  CircEmbed,     // Circulant embedding - must always be the first one!  //0
  CircEmbedCutoff,
  CircEmbedIntrinsic,
  TBM,           // Turning Bands, performed in 2 or 3 dimensions 
  SpectralTBM,   // Spectral turning bands, only 2 dimensional
  
  Direct,        // directly, by matrix inversion  // 5
  Sequential,    // sequential simulation 
  Trendproc,        // Trend evaluation
  Average,       // Random spatial averages 
  Nugget,        // just simulate independent variables
  
  RandomCoin,    // "additive MPP (random coins)"  // 10
  Hyperplane,    // by superposition of Poisson hyperplane tessellations;  only in the plane! Not implemented yet
  Specific,      // Specific Methods 
  Nothing,       // (*)
  Forbidden      // must always be the last one 
} Methods;

// (*) [continued]  must always be the last of the list of methods for 
//		    simulation of Gaussian random fields;
//
//		    NOTE: Nothing has three meanings:
//		    (1) used in preference list whether variogram/covariance fct
//		        can be calculated
//		    (2) method equals Nothing, if no method is preferred
//                    (3) method is not shown, if ErrorMessage is called with
//		        Nothing

typedef enum sort_origin {
  original, mle_conform, all_origins
} sort_origin;

#define INTERNAL_PARAM "internal"



#define GETMODEL_AS_SAVED 0    
#define GETMODEL_DEL_NATSC 1   
#define GETMODEL_SOLVE_NATSC 2
#define GETMODEL_DEL_MLE 3
#define GETMODEL_SOLVE_MLE 4

#define MIXED_X_NAME "X"
#define MIXED_BETA_NAME "beta"

#define COVARIATE_C_NAME "data"
#define COVARIATE_X_NAME "x"
#define COVARIATE_RAW_NAME "raw"
#define COVARIATE_ADDNA_NAME "addNA"

#define CONST_A_NAME "x"

#define COVARIATE_DATA_NAME "data"

#define MINMAX_PMIN 1
#define MINMAX_PMAX 2
#define MINMAX_TYPE 3
#define MINMAX_NAN 4
#define MINMAX_MIN 5
#define MINMAX_MAX 6
#define MINMAX_OMIN 7
#define MINMAX_OMAX 8
#define MINMAX_ROWS 9
#define MINMAX_COLS 10
#define MINMAX_BAYES 11
#define MINMAX_ENTRIES MINMAX_BAYES



#define XLIST_X 0
#define XLIST_Y 1
#define XLIST_T 2
#define XLIST_GRID 3
#define XLIST_SPATIALDIM 4
#define XLIST_TIME 5
#define XLIST_DIST 6
#define XLIST_RESTOT 7
#define XLIST_L 8
#define XLIST_UNITS 9
#define XLIST_NEWUNITS 10
#define XLIST_ENTRIES  (XLIST_NEWUNITS + 1)


#define PROJ_SPACE -1
#define PROJ_TIME -2
#define PROJECTIONS 2

#define VAR2COV_EXTREMAL -1
#define VAR2COV_ORIGIN -2
#define VAR2COV_CENTER -3
#define VAR2COV_METHODS 3


#define INTERN_SHOW 2


#define METHOD_VARIOGRAM 0
#define METHOD_PSEUDO 1
#define METHOD_COVARIANCE 2
#define METHOD_PSEUDOMADOGRAM 3
#define METHOD_MADOGRAM 4

#define EV_FFT_EV 0
#define EV_FFT_N 1
#define EV_EV 0
#define EV_SDSQ 1
#define EV_N 2

#define BRP_UNIF 0
#define BRP_COV 1
#define BRP_KRIG 2

//#define RMCOV_X 0
//#define RMCOV_A 1
#define VAR2COV_X 0
#define VAR2COV_C 1
#define VAR2COV_SPECIAL 1


#define POISSON_SCATTER_OPTIM 0
#define POISSON_SCATTER_STANDARD 1
#define POISSON_SCATTER_ANY 2
#define nPOISSON_SCATTER (POISSON_SCATTER_ANY + 1)

extern const char *ISO_NAMES[LAST_ISO + 1],
  *EQ_NAMES[6],
  *OPTIMISER_NAMES[nOptimiser], *NLOPTR_NAMES[nNLOPTR],
  *LIKELIHOOD_NAMES[nLikelihood],
  *DOMAIN_NAMES[LAST_DOMAIN + 1],
  *TYPE_NAMES[LASTTYPE + 1],
  *MONOTONE_NAMES[MONOTONE_TOTAL] ,
  *MODE_NAMES[nr_modes], *OUTPUTMODE_NAMES[nr_output_modes], 
  *ANGLE_NAMES[last_angle_mode + 1],
  *REPORTCOORD_NAMES[nr_reportcoord_modes],
  *UNITS_NAMES[nr_units], 
  *COORD_SYS_NAMES[nr_coord_sys],
  *CARTESIAN_SYS_NAMES[3],
  *PROJECTION_NAMES[nProjections],
  *TYPEOF_PARAM_NAMES[LASTSORTOF + 1],
  *NAMES_OF_NAMES[nNamesOfNames],
  *RMCOV_X[nVAR2COV_METHODS],
  *POISSON_SCATTER_NAMES[nPOISSON_SCATTER];
 

#endif
