#ifndef AutoRandomFields_H
#define AutoRandomFields_H 1
#include <R.h>

#define MAXCHAR 18 // max number of characters for (covariance) names  
#define METHODMAXCHAR MAXCHAR // max number of character to describe a Method (including \0)


#define MAXCOVDIM 11000
#define MAXMLEDIM MAXCOVDIM
#define MAXSIMUDIM MAXCOVDIM
#define MAXSUB 10

#define MAXCEDIM 13
#define MAXTBMSPDIM 4
#define MAXMPPDIM 4
#define MAXMPPVDIM 10
#define MAXHYPERDIM 4
#define MAXNUGGETDIM 20
#define MAXVARIODIM 20
#define MAXTBMVDIM 5
#define MAXGETNATSCALE 5
#define MAXGAUSSVDIM 10

#define MAXPARAM 20


typedef int domain_type;
typedef int isotropy_type;

#define MAXUNITS 4
typedef enum units_enum {
  units_none, units_km, units_miles, units_time, units_user} units_enum;
#define nr_units (units_user + 1)

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


// 0 sicherheitshalber nirgendswo nehmen 
#define PARAM_DEP -1 // parameter dependent fuer vdim, finiterange, etc; for CovList only, falls von submodel und param abhaengig, so submodel_dep 
#define PREVMODEL_DEP -2 
#define SUBMODEL_DEP -3  // immer, wenn zumindest submodel depedence
//#define PARAMSUBMODELM -2 // for CovList only 
#define MISMATCH -4
//#define NON_QUADRATIC -5

// NOTE!!!
// if any definition is changed check setbackward() 
#define XONLY (domain_type) 0   // TRANS_INV
#define KERNEL (domain_type) 1 // not a genuine covariance function;
                              // takes only x vector, but is not statonary 
#define PREVMODELD (domain_type) 2 // type taken from calling models 
#define DOMAIN_MISMATCH (domain_type) 3// always last ! 
#define LAST_DOMAIN DOMAIN_MISMATCH

// NOTE!!
// !!! CHANGE ALSO RF_GLOBALS.R !!!
#define ISOTROPIC (isotropy_type) 0 // RC genauer : rotation invariant!!
#define SPACEISOTROPIC (isotropy_type) 1 // RC fully symmetric 
#define ZEROSPACEISO (isotropy_type) 2   // isotropic if time diff is zero 
#define VECTORISOTROPIC (isotropy_type) 3
#define SYMMETRIC (isotropy_type) 4 //"stationary" only if XONLY; in the covariance sense if multivariate!
#define CARTESIAN_COORD (isotropy_type) 5 // RC
#define GNOMONIC_PROJ (isotropy_type) 6  // RC projection on the plane
#define ORTHOGRAPHIC_PROJ (isotropy_type) 7 // RC dito; to do: implement further projections
#define LAST_CARTESIAN ORTHOGRAPHIC_PROJ // danach MUSS sphaerisch kommen!!
#define SPHERICAL_ISOTROPIC (isotropy_type) 8
#define SPHERICAL_SYMMETRIC (isotropy_type) 9
#define SPHERICAL_COORDS (isotropy_type) 10 // RC
#define EARTH_ISOTROPIC (isotropy_type) 11
#define EARTH_SYMMETRIC (isotropy_type) 12
#define EARTH_COORDS (isotropy_type) 13 // RC 
#define LAST_ANYSPHERICAL EARTH_COORDS
#define CYLINDER_COORD (isotropy_type) 14 // unused
#define UNREDUCED (isotropy_type) 15  //CARTESIAN_COORD or EARTH_COORD  
#define PREVMODELI (isotropy_type) 16  // type taken from submodels 
#define ISO_MISMATCH  (isotropy_type) 17 // always last ! 
#define LAST_ISO ISO_MISMATCH

#define MON_PARAMETER -1
#define NOT_MONOTONE 0
#define MONOTONE 1
#define GNEITING_MON 2 // Euclid's hat, Gneiting, J.Mult.Anal. 69, 1999
#define NORMAL_MIXTURE 3
#define COMPLETELY_MON 4
#define BERNSTEIN 5 // is different from everything else before


#define MAXFIELDS 10
#define MODEL_USER (MAXFIELDS + 0)  // for user call of Covariance etc. 
#define MODEL_AUX  (MAXFIELDS + 1)  // auxiliary in fitgauss.R
#define MODEL_INTERN (MAXFIELDS + 2) // for kriging, etc; internal call of cov 
#define MODEL_SPLIT (MAXFIELDS + 3) // split covariance model and other auxiliary methods 
#define MODEL_GUI (MAXFIELDS + 4)   // RFgui 
#define MODEL_MLE (MAXFIELDS + 5) // mle covariance model 
#define MODEL_MLESPLIT (MAXFIELDS + 6)  // ="= 
#define MODEL_LSQ (MAXFIELDS + 7)  // unused
#define MODEL_BOUNDS (MAXFIELDS + 8)  // MLE, lower, upper 
#define MODEL_KRIGE  (MAXFIELDS + 9)
#define MODEL_PREDICT  (MAXFIELDS + 10) 
#define MODEL_ERR  (MAXFIELDS + 11)
#define MODEL_MAX MODEL_ERR


//////////////////////////////////////////////////////////////////////
// the different levels of printing

#define PL_IMPORTANT 1 
#define PL_SUBIMPORTANT 2
#define PL_RECURSIVE 3
#define PL_REC_DETAILS 4 
#define PL_STRUCTURE 5 // see also initNerror.ERROROUTOFMETHOD
#define PL_ERRORS  6 // only those that are caught internally

#define PL_FCTN_DETAILS 7  // R
#define PL_FCTN_SUBDETAILS 8

#define PL_COV_STRUCTURE 7 // C
#define PL_DIRECT_SEQU 8
#define PL_DETAILS 9
#define PL_SUBDETAILS 10

typedef enum Types {
  TcfType, PosDefType, VariogramType, NegDefType, 
  ProcessType, GaussMethodType, 
  BrMethodType, // change also rf_globals.R if deleted
  PointShapeType, RandomType, ShapeType, TrendType, InterfaceType, 
  RandomOrShapeType, // only for parameters, e.g. nu in Whittle
  UndefinedType, // Bedeutung: C->Type:TypeFct existiert; cov->typus:ungesetzt
  MathDefinition,
  OtherType, // always last for usual use
  NN1, NN2, NN3, NN4 // for use in generateMmodels only
} Types;

#define nOptimiser 8
#define nNLOPTR 15
#define nLikelihood 4

extern const char *ISONAMES[LAST_ISO + 1], 
  *OPTIMISER_NAMES[nOptimiser], *NLOPTR_NAMES[nNLOPTR],
  *LIKELIHOOD_NAMES[nLikelihood],
  *DOMAIN_NAMES[LAST_DOMAIN + 1], *TYPENAMES[OtherType + 1],
  *MONOTONE_NAMES[BERNSTEIN + 1 - MISMATCH] ,
  *MODENAMES[nr_modes], *OUTPUTMODENAMES[nr_output_modes], 
  *REPORTCOORDNAMES[nr_reportcoord_modes],
  *UNITS_NAMES[nr_units], 
  *COORD_SYS_NAMES[nr_coord_sys];

typedef enum sortsofeffect { // ! always compare with convert.R!
  DetTrendEffect, // trend, nichts wird geschaetzt
  FixedTrendEffect, // linearer Parameter wird geschaetzt  
  FixedEffect, // trend is also converted to FixedEffect ?
  RandomEffect, // b is random, no variance is estimated; cov. matrix
  RvarEffect, // b is random, variance is estimated; covariance matrix
  LargeEffect, //  wie RandomEffect, aber gross, somit keine Optimierung
  LVarEffect, //   wie RVarEffekt, aber gross, somit keine Optimierung
  SpaceEffect, // spatial covariance model for random effect
  SpVarEffect, //spatial covariance model for random effect
  DataEffect, // box cox
  RemainingError, // error term:
  effect_error //
} sortsofeffect;

#define FirstMixedEffect FixedEffect
#define LastMixedEffect SpVarEffect

//////////////////////////////////////////////////////////////////////
// the different types of parameters
typedef enum sortsofparam { // never change ordering; just add new ones !!!!!!
  VARPARAM, SIGNEDVARPARAM, SDPARAM, SIGNEDSDPARAM, // 0..3
  SCALEPARAM, DIAGPARAM, ANISOPARAM, // 4..6
  INTEGERPARAM, ANYPARAM, TRENDPARAM, // unused , // 7..9
  NUGGETVAR, MIXEDVAR,  // unused   // 10..12
  CRITICALPARAM, 
  IGNOREPARAM, //  not recognised in mle, e.g. when some parameter is doubled by another, e.g. biWM
  DONOTVERIFYPARAM, 
  DONOTRETURNPARAM, // some as IGNOREPARAM, but also not returned by key
  FORBIDDENPARAM
} sortsofparam;


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
  TrendEval,        // Trend evaluation
  Average,       // Random spatial averages 
  Nugget,        // just simulate independent variables 
  RandomCoin,    // "additive MPP (random coins)"  // 10
  Hyperplane,    // by superposition of Poisson hyperplane tessellations;  only in the plane! Not implemented yet
  Specific,      // Specific Methods 
  Nothing,       // (*)
  Forbidden      // must always be the last one 
} SimulationMeth;

// (*) [continued]  must always be the last of the list of methods for 
//		    simulation of Gaussian random fields;
//
//		    NOTE: Nothing has three meanings:
//		    (1) used in preference list whether variogram/covariance fct
//		        can be calculated
//		    (2) method equals Nothing, if no method is preferred
//                    (3) method is not shown, if ErrorMessage is called with
//		        Nothing
		 
#define INTERNAL_PARAM "internal"

#define MAXVARIANTS 6


#define GETMODEL_AS_SAVED 0    
#define GETMODEL_DEL_NATSC 1   
#define GETMODEL_SOLVE_NATSC 2
#define GETMODEL_DEL_MLE 3
#define GETMODEL_SOLVE_MLE 4

#define MIXED_X_NAME "X"
#define MIXED_BETA_NAME "beta"

#define COVARIATE_C_NAME "c"
#define COVARIATE_X_NAME "x"
#define COVARIATE_ADDNA_NAME "addNA"

#define CONST_A_NAME "a"


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


#define PROJ_SPACE -2 
#define PROJ_TIME -1
#define PROJECTIONS 2

#define INTERN_SHOW 2


#endif
