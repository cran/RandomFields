#ifndef RFsimu_H
#define RFsimu_H 1

#include "auxiliary.h"
#include "RandomFields.h"
#include "error.h"
#include <string.h>

//   intptr_t and uintptr_t fuer: umwandlung pointer in int und umgekehrt

#define ASSERT_GATTERONLY(Cov) assert(TrafoOK(Cov, false))
#define ASSERT_GATTER(Cov) assert(TrafoOK(Cov, true))
#define ASSERT_CHECKED(Cov) assert(Cov->checked)

#define DOPRINT !true

// #define RF_DEBUGGING 1



#ifndef RF_DEBUGGING

//#define ERRLINE 
#define NICK(COV) (isDollar(COV) ? CovList[(COV)->sub[0]->nr].nick : CovList[(COV)->nr].nick)

#define ERRLINE assert({PRINTF("(ERROR %s, line %d)\n", __FILE__, __LINE__); true;})

#define MEMCOPY(A,B,C) memcpy(A,B,C)
//#define MEMCOPY(A,B,C) memory_copy(A, B, C)

#define MALLOC malloc
#define CALLOC calloc

#define CHECK(C,T,X,type,D,I,V,R) check2X(C,T,X,type,D,I,V,R)
#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,R) check2X(C,T,X,type,D,I,V0,V1,R)
#define CHECKPD2ND(N,D1,D2,I,V,R) CheckPD2ND(N,D1,D2,I,V,R)
#define INIT(M, Moments, S) INIT_intern(M, Moments, S)
#define REINIT(M, Moments, S) REINIT_intern(M, Moments, S)
#define INIT_RANDOM(M, Moments, S, P) INIT_RANDOM_intern(M, Moments, S, P)
#define STRUCT(Cov, NM) CovList[(Cov)->gatternr].Struct(Cov, NM)

#define PARAM(P, IDX) (P)->px[IDX]
#define PARAMINT(P, IDX) ((int *) (P)->px[IDX])
#define PARAMENV(P, IDX) ((sexp_type *) (P)->px[IDX])
#define PARAMLIST(P, IDX) ((listoftype *) (P)->px[IDX])

#define PARAM0(P, IDX) PARAM(P, IDX)[0]
#define PARAM0INT(P, IDX) PARAMINT(P, IDX)[0]

#define PCOPY(TO, FROM, IDX) 						\
  MEMCOPY((TO)->px[IDX], (FROM)->px[IDX],				\
	  ((FROM)->nrow[IDX]) * ((FROM)->ncol[IDX]) *			\
	  (CovList[FROM->nr].kappatype[IDX]==REALSXP ? sizeof(double) : \
	    CovList[FROM->nr].kappatype[IDX]==INTSXP ? sizeof(int) :	\
	    -1))
 
#define DEBUGINFOERR
#define DEBUGINFO 

#else
#define NICK(COV) (CovList[(COV)->nr].nick)
#define ERRLINE assert({PRINTF("(ERROR %s, line %d)\n", __FILE__, __LINE__); true;})


#define MEMCOPY(A,B,C) ({ assert((A)!=NULL && (B)!=NULL); memcpy(A,B,C); })
//#define MEMCOPY(A,B,C) memory_copy(A, B, C)

//
#define MALLOC(X) ({DOPRINTF("(MALLOC %s, line %d)\n", __FILE__, __LINE__);assert(X>0 && X<=536870912);malloc(X);})
//
#define CALLOC(X, Y) ({DOPRINTF("(CALLOC %s, line %d)\n",__FILE__, __LINE__);assert(X>0 && X<1e8 && Y>0 && Y<=64); calloc(X,Y);})
//#define MALLOC malloc
//#define CALLOC calloc

#define XX(C) assert((C)->simu.expected_number_simu >= 0|| ({DOPRINTF("Start.\n"); false;}))
#define YY(C) assert((C)->simu.expected_number_simu >= 0 || ({DOPRINTF("End.\n"); false;}))

#define LLPRINT(SIGN, cov, Z) {			\
    cov_model *lprint_z=cov;				\
    int lprint_i=0;					\
    while (lprint_z->calling != NULL && lprint_i<10) {	\
      lprint_z=lprint_z->calling;			\
      DOPRINTF(SIGN); DOPRINTF(" ");				\
      lprint_i++;}					\
    if (lprint_i==100) {				\
      PRINTF("LPRINT i=%d\n", lprint_i);		\
      PMI(cov);						\
      assert(false);					\
    }							\
  }DOPRINTF("(");DOPRINTF("%s", Z);DOPRINTF(" %s, line %d : %s\n", __FILE__, __LINE__, NAME(cov))

#define CHECKSIGN "_"
#define INITSIGN "_"
#define STRUCTSIGN "_"
#define CHECK(C,T,X,type,D,I,V,R) ({LLPRINT(CHECKSIGN, C, "CHECK"); assert(type!=RandomType); XX(C); int _x = check2X(C,T,X,type,D,I,V,R); YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED");} _x;})
#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,R) ({LLPRINT(CHECKSIGN, C, "CHECKVDIM"); XX(C); int _x = check2X(C,T,X,type,D,I,V0,V1,R);YY(C); _x;})
#define CHECKPD2ND(N,D1,D2,I,V,R) ({LLPRINT(CHECKSIGN, N, "CHECKPD2ND"); XX(N); int _x = CheckPD2ND(N,D1,D2,I,V,R);YY(N); _x;})
#define INIT(M, Moments, S) ({LLPRINT(INITSIGN, M, "INIT");  XX(M); int _x = INIT_intern(M, Moments, S);YY(M); LLPRINT(INITSIGN, M, "INIT DONE"); _x;})
#define REINIT(M, Moments, S) ({LLPRINT(INITSIGN, M, "INIT");  XX(M); int _x = REINIT_intern(M, Moments, S); YY(M); _x;})

#define INIT_RANDOM(M, Moments, S, P) ({LLPRINT(INITSIGN, M, "INITRANDOM");  XX(M); int _x = INIT_RANDOM_intern(M, Moments, S, P);YY(M); _x;})

#define STRUCT(Cov, NM)  ({LLPRINT(STRUCTSIGN, Cov, "STRUCT"); ASSERT_GATTER(Cov);  XX(Cov); int _x = CovList[(Cov)->gatternr].Struct(Cov, NM);YY(Cov); if (_x==NOERROR){LLPRINT(STRUCTSIGN, Cov, "STRUCT DONE\n");}else{LLPRINT(STRUCTSIGN, Cov, "STRUCT FAILED");}_x;})


#define PARAM(P, IDX) ({assert(CovList[(P)->nr].kappatype[IDX] == REALSXP); (P)->px[IDX];})
#define PARAMINT(P, IDX) ({assert(CovList[(P)->nr].kappatype[IDX] == INTSXP); ((int *) (P)->px[IDX]);})
#define PARAMENV(P, IDX) ({assert(CovList[(P)->nr].kappatype[IDX] == LANGSXP); ((sexp_type *) (P)->px[IDX]);})
#define PARAMLIST(P, IDX) ({assert(CovList[(P)->nr].kappatype[IDX] == LANGSXP); ((listoftype *) (P)->px[IDX]);})

#define PARAM0(P, IDX) ({assert(PARAM(P, IDX) != NULL); PARAM(P, IDX)[0];})
#define PARAM0INT(P, IDX)  ({assert(PARAMINT(P, IDX) != NULL); PARAMINT(P, IDX)[0];})


#define PCOPY(TO, FROM, IDX) {						\
    if (!((TO) != NULL && (FROM) != NULL &&  (TO)->px[IDX] != NULL && (FROM)->px[IDX] &&IDX > 0 && (FROM)->ncol[IDX] == (TO)->ncol[IDX] &&   (FROM)->nrow[IDX] == (TO)->nrow[IDX] && CovList[(FROM)->nr].kappatype[IDX]==CovList[(TO)->nr].kappatype[IDX] &&  (CovList[(FROM)->nr].kappatype[IDX]==REALSXP ||	     CovList[(FROM)->nr].kappatype[IDX]==INTSXP)  )) BUG; \
    MEMCOPY((TO)->px[IDX], (FROM)->px[IDX],				\
	    ((FROM)->nrow[IDX]) * ((FROM)->ncol[IDX]) *			\
	   (CovList[(FROM)->nr].kappatype[IDX]==REALSXP ? sizeof(double) : \
	    CovList[(FROM)->nr].kappatype[IDX]==INTSXP ? sizeof(int) :	\
	    -1));							\
}

   /*  printf("//%ld %ld %d %d %d idx=%d\n", (TO)->px[IDX], (FROM)->px[IDX], \
       (FROM)->nrow[IDX], (FROM)->ncol[IDX],				\
       CovList[(FROM)->nr].kappatype[IDX]==REALSXP ? sizeof(double) :	\
       CovList[(FROM)->nr].kappatype[IDX]==INTSXP ? sizeof(int) :	\
       -1, IDX);							\
				 */

#define DEBUGINFOERR {							\
    char dummy[MAXERRORSTRING]; strcpy(dummy, ERRORSTRING);		\
    sprintf(ERRORSTRING, "%s (%s, line %d)\n", dummy, __FILE__, __LINE__); \
  }
#define DEBUGINFO DOPRINTF("(info:  %s, line %d)\n", __FILE__, __LINE__)

#endif



#define P(IDX) PARAM(cov, IDX) 
#define PINT(IDX) PARAMINT(cov, IDX)
#define P0(IDX) PARAM0(cov, IDX) 
#define P0INT(IDX) PARAM0INT(cov, IDX)
#define PENV(IDX) PARAMENV(cov, IDX) 
#define PSEXP PENV
#define PLIST(IDX) PARAMLIST(cov, IDX)

#define PARAMFREE(P, IDX) if ((P)->px[IDX] != NULL) { free((P)->px[IDX]); (P)->px[IDX] = NULL; (P)->ncol[IDX] = (P)->nrow[IDX] = 0;}


#define PARAMALLOC(P, IDX, ROW, COL) {					\
    int _PARAMsize;							\
    switch(CovList[(P)->nr].kappatype[IDX]) {				\
    case REALSXP : _PARAMsize = sizeof(double); break;			\
    case INTSXP : _PARAMsize = sizeof(int); break;			\
    default : BUG;							\
    }									\
    assert((P)->px[IDX]==NULL);						\
    (P)->nrow[IDX] = ROW; (P)->ncol[IDX] = COL;				\
    if (((P)->px[IDX] =							\
	 (double*) CALLOC((ROW) * (COL), _PARAMsize)) == NULL) {	\
      XERR(ERRORMEMORYALLOCATION)					\
    }									\
  }
#define PALLOC(IDX, ROW, COL) PARAMALLOC(cov, IDX, ROW, COL)
#define PINRALLOC(IDX, ROW, COL) PARAMINTALLOC(cov, IDX, ROW, COL)

#define PFREE(IDX) PARAMFREE(cov, IDX)
#define PARAMisNULL(P, IDX) ((P)->px[IDX] == NULL)
#define PisNULL(IDX) PARAMisNULL(cov, IDX)

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
#define MAXMPPVDIM 10
#define MAXHYPERDIM 4
#define MAXNUGGDIM 20
#define MAXVARIODIM 20
#define MAXTBMVDIM 5
#define MAXGETNATSCALE 5
/// !!! if list get longer change GetrfParameters in userinterface.cc
// !!! und RFparameters

#define MAXMAKEEXPLICITE 20

///////////////////////////////////////////////////////////////////////
// basic types of (dynamic) model membership
///////////////////////////////////////////////////////////////////////

#define ROLE_BASE 0 // to check the most basic things
#define ROLE_COV 1
#define ROLE_GAUSS 2 // including  alpha-stable?!
#define ROLE_MAXSTABLE 3  // reserved for indicating a shape function
//                           in the max-stable (Smith) context
#define ROLE_BROWNRESNICK 4
#define ROLE_SMITH 5  // process itself? and tcf
#define ROLE_SCHLATHER 6
#define ROLE_POISSON 7
#define ROLE_POISSON_GAUSS 8
#define ROLE_BERNOULLI 9
#define ROLE_DISTR 10
#define ROLE_FAILED 11
#define ROLE_UNDEFINED 12
#define ROLE_LAST ROLE_UNDEFINED



///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////
// type of covariance functions that need distinguished treatment
// Reihenfolge beachten! Hoehere Nummern sind meist in kleinere umwandelbar


// *****************
// NOTE!!!
// if any definition is changed check setbackward() 
// NOTE!!
// !!! CHANGE ALSO RF_GLOBALS.R !!!
#define XONLY (domain_type) 0   // TRANS_INV
#define KERNEL (domain_type) 1 /* not a genuine covariance function;
			takes only x vector, but is not statonary */
#define PREVMODELD (domain_type) 2 /* type taken from calling models */
#define STAT_MISMATCH (domain_type) 3 /* always last ! */
#define LAST_STAT STAT_MISMATCH

// NOTE!!
// !!! CHANGE ALSO RF_GLOBALS.R !!!
#define ISOTROPIC (isotropy_type) 0 // genauer : rotation invariant!!
#define SPACEISOTROPIC (isotropy_type) 1 /* fully symmertric */
#define ZEROSPACEISO (isotropy_type) 2   /* isotropic if time diff is zero */
#define VECTORISOTROPIC (isotropy_type) 3
#define SYMMETRIC (isotropy_type) 4 // in the covariance sense if multivariate!
#define CARTESIAN_COORD (isotropy_type) 5
#define EARTH_COORD (isotropy_type) 6
#define SPHERICAL_COORD (isotropy_type) 7
#define CYLINDER_COORD (isotropy_type) 8
#define UNREDUCED (isotropy_type) 9  /*CARTESIAN_COORD  or EARTH_COORD  */
#define PREVMODELI (isotropy_type) 10  /* type taken from submodels */
#define ISO_MISMATCH  (isotropy_type) 11 /* always last ! */
#define LAST_ISO ISO_MISMATCH

 //#define ISOTROPIC (isotropy_type) 0
 //#define SPACEISOTROPIC (isotropy_type) 1 /* fully symmertric */
 //#define ZEROSPACEISO (isotropy_type) 2   /* isotropic if time diff is zero */
 //#define VECTORISOTROPIC (isotropy_type) 3
 //#define ANISOTROPIC (isotropy_type) 4
 //#define PREVMODELI (isotropy_type) 5  /* type taken from submodels */
 //#define ISO_MISMATCH  (isotropy_type) 9 /* always last ! */
 //#define LAST_ISO 9


// vdim is "strictly" passed from bottom to top
#define SCALAR 1 
/* 0 sicherheitshalber nirgendswo nehmen */
#define PARAM_DEP -1 /* parameter dependent fuer vdim, finiterange, etc; for CovList only, falls von submodel und param abhaengig, so submodel_dep */
#define PREVMODEL_DEP -2 
#define SUBMODEL_DEP -3  // immer, wenn zumindest submodel depedence
//#define PARAMSUBMODELM -2 /* for CovList only */
#define MISMATCH -4
//#define NON_QUADRATIC -5

// way of implementing simulation methods:
#define NOT_IMPLEMENTED 0 /* do not change this value except to false */
#define IMPLEMENTED 1     /* do not change this value except to true */
#define NUM_APPROX 2 
#define GIVEN_METH_IGNORED 3
#define HYPERIMPLEMENTED 4

#define NOT_MONOTONE 0
#define MONOTONE 1
#define GNEITING_MON 2 // Euclid's hat, Gneiting, J.Mult.Anal. 69, 1999
#define NORMAL_MIXTURE 3
#define COMPLETELY_MON 4
#define BERNSTEIN 5 // is different from everything else before


///////////////////////////////////////////////////////////////////////
// BASIC DIMENSIONS AND VARIABLES
///////////////////////////////////////////////////////////////////////
#define MAXCHAR 17 /* max number of characters for covariance name  */
#define PARAMMAXCHAR MAXCHAR
#define METHODMAXCHAR MAXCHAR /* max number of character to describe a Method (including \0)*/
#define MAXLONGCHAR 40 // for RFoptions that are string
#define MAXUNITSCHAR 10
#define MAXUNITS 4

#define MAXNRCOVFCTS 200
#define MAXPARAM 20
#define MAXELEMENTS 100
#define MAX_NA 30
#define MAX_MLE_ELMNTS 10
#define MAXSUB 10
/* max number of submodels -- NEVER MORE THAN 255 */
#define MAXTBMDIM 3
#define MAXTAYLOR 3
#/* if explicite simulation method is not given by ????????
   user, the programme will first try a methods 
   according to "first"; 
   then it tries the other methods according to
   a default list: (CircEmbed only if the simulation is on a grid)
   1 dim: CircEmbed, Direct
   2 dim: CircEmbed, TBM, spectralTBM, 
   AdditiveMpp, Direct
   3 dim: CircEmbed, TBM, Direct
   the method actually used is stored in actually_used. 
*/ 
#define MAXFIELDS 10
#define MODEL_USER (MAXFIELDS+0)  /* for user call of Covariance etc. */
#define MODEL_UNUSED  (MAXFIELDS+1)  /* unused */ 
#define MODEL_INTERN (MAXFIELDS+2) /* for kriging, etc; internal call of cov */
#define MODEL_SPLIT (MAXFIELDS+3) /* split covariance model and other auxiliary
				     methods */
#define MODEL_GUI (MAXFIELDS+4)   /* RFgui */
#define MODEL_MLE (MAXFIELDS+5) /* mle covariance model */
#define MODEL_MLESPLIT (MAXFIELDS+6)  /* ="= */
#define MODEL_MLETREND (MAXFIELDS+7)  /* mle trend model !! */
#define MODEL_BOUNDS (MAXFIELDS+8)  /* MLE, lower, upper */
#define MODEL_KRIGE  (MAXFIELDS+9)
#define MODEL_COND  (MAXFIELDS+10)
#define MODEL_ERR  (MAXFIELDS+11)
#define MODEL_MAX MODEL_ERR



#define LISTOF 100  /* may not interfere with define of *SXP in Rinternal.h*/


// #define MAXINTERNALPARAM 2
#define UNKNOWN RF_NAN;
extern double ZERO[MAXSIMUDIM], ONE;
extern int True;
extern int False;
extern char *FREEVARIABLE;

typedef char name_type[][MAXCHAR];
typedef char longname_type[][MAXLONGCHAR];

typedef char NAname_type[MAX_NA][255];

#define GETMODEL_AS_SAVED 0
#define GETMODEL_DEL_NATSC 1
#define GETMODEL_SOLVE_NATSC 2
#define GETMODEL_DEL_MLE 3
#define GETMODEL_SOLVE_MLE 4


///////////////////////////////////////////////////////////////////////
// PARAMETER AND THEIR LOCATIONS FOR VARIOUS MODELS
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// Dollar
#define DOLLAR_SUB 0
#define DVAR 0
#define DSCALE 1
#define DANISO 2 /* internal */
#define DALEFT 3 /* user defined */
#define DPROJ 4
#define DMAX DPROJ

// POWER_DOLLAR
#define POWVAR 0
#define POWSCALE 1
#define POWPOWER 2
#define POW_SUB 0

///////////////////////////////////////////////////////////////////////
//  null_model
#define NULL_TYPE 0

///////////////////////////////////////////////////////////////////////
//  fractal Brown
#define BROWN_ALPHA 0
#define BROWN_GEN_BETA 1


///////////////////////////////////////////////////////////////////////
// whittle + matern
#define WM_NU 0
#define WM_NOTINV 1


///////////////////////////////////////////////////////////////////////
// truncsupport
#define TRUNC_RADIUS 0


///////////////////////////////////////////////////////////////////////
// Angle
#define ANGLE_ANGLE 0
#define ANGLE_RATIO 1
#define ANGLE_DIAG 2


///////////////////////////////////////////////////////////////////////
// RRspheric
#define SPHERIC_SPACEDIM 0
#define SPHERIC_BALLDIM 1



///////////////////////////////////////////////////////////////////////
// models for method "average" (mixture of spectral and Poisson)
#define AVERAGE_YPHASE 0
#define AVERAGE_YFREQ 1


///////////////////////////////////////////////////////////////////////
// point shape function
#define PGS_FCT 0 // never change
#define PGS_LOC 1

#define PGS_RATIO 0
#define PGS_FLAT 1
#define PGS_INFTY_SMALL 2
#define PGS_NORMED 3
#define PGS_ISOTROPIC 4

///////////////////////////////////////////////////////////////////////
// trend
#define TREND_MEAN 0
#define TREND_LINEAR 1
#define TREND_POLY 2
#define TREND_PARAM_POLY 3
#define TREND_FCT 4
#define TREND_PARAM_FCT 5

#define MIXED_ELMNT 0
#define MIXED_X 1
#define MIXED_BETA 2
#define MIXED_COORD 3
#define MIXED_DIST 4
#define MIXED_DIM 5

///////////////////////////////////////////////////////////////////////
// gaussproc
#define COMMON_GAUSS -1
#define GAUSSPROC_STATONLY  (COMMON_GAUSS + 1)
#define GAUSSPROC_LAST GAUSSPROC_STATONLY



///////////////////////////////////////////////////////////////////////
// TBM
#define TBM_FULLDIM (COMMON_GAUSS + 1)
#define TBM_TBMDIM (COMMON_GAUSS + 2)
#define TBM_LAYERS (COMMON_GAUSS + 3)

#define TBMOP_FULLDIM 0
#define TBMOP_TBMDIM 1
#define TBMOP_LAYERS 2

#define TBMOP_COV 0


///////////////////////////////////////////////////////////////////////
// constant
#define CONSTANT_ELMNT 0
#define CONSTANT_M 1
#define CONSTANT_VDIM 2


///////////////////////////////////////////////////////////////////////
// coins
#define COIN_COV 0
#define COIN_SHAPE 1

///////////////////////////////////////////////////////////////////////
// local ce covariance operator und q-values for ce method / ce operator
#define pLOC_DIAM 0 // parameter p
#define pLOC_A 1 // always last of pLOC_*
#define pLOC_R pLOC_A // always last of pLOC_*

#define LOCAL_R 0 
#define LOCAL_MSG (LOCAL_R + 1)

#define CUTOFF_B (LOCAL_MSG + 1)
#define CUTOFF_ASQRTR (CUTOFF_B + 1)
#define CUTOFF_THEOR (CUTOFF_ASQRTR + 1) // muss immer == 4 sein
#define CUTOFF_MAX (CUTOFF_THEOR + 1)   /* size of vector q */

#define INTRINSIC_A0 (LOCAL_MSG + 1)
#define INTRINSIC_A2 (INTRINSIC_A0 + 1)
#define INTRINSIC_B (INTRINSIC_A2 + 1)
#define INTRINSIC_MAX (INTRINSIC_B + 1) /* size of vector q */
#define LOCAL_MAX INTRINSIC_MAX


///////////////////////////////////////////////////////////////////////
// random poisson hyperplane polygon
#define POLYGON_BETA 0
#define POLYGON_SAFETY 1


///////////////////////////////////////////////////////////////////////
// Families

#define GAUSS_DISTR_MEAN 0
#define GAUSS_DISTR_SD 1
#define GAUSS_DISTR_LOG 2

#define UNIF_MIN 0
#define UNIF_MAX 1
#define UNIF_NORMED 2

#define SETPARAM_LOCAL 0
#define SET_PERFORMDO 0
//#define SETPARAM_VARIANT 1
//#define SETPARAM_FROM 1
//#define SETPARAM_SET 2

#define LOC_LOC 0
#define LOC_SCALE 1
#define LOC_POWER 2


///////////////////////////////////////////////////////////////////////
// max-stable
#define GEV_XI 0
#define GEV_MU 1
#define GEV_S 2
#define LAST_MAXSTABLE GEV_S

#define MPP_SHAPE 0
#define MPP_TCF 1

///////////////////////////////////////////////////////////////////////
// ARCSQRT
#define ARCSQRT_SCALE 0


///////////////////////////////////////////////////////////////////////
// rectangular

#define RECT_SAFETY 0
#define RECT_MINSTEPLENGTH 1
#define RECT_MAXSTEPS 2
#define RECT_PARTS 3
#define RECT_MAXIT 4
#define RECT_INNERMIN 5
#define RECT_OUTERMAX 6
#define RECT_MCMC_N 7
#define RECT_NORMED 8
#define RECT_APPROX 9
#define RECT_ONESIDED 10

   
///////////////////////////////////////////////////////////////////////
// (mle) select;
#define SELECT_SUBNR 0

///////////////////////////////////////////////////////////////////////
// Preferences for the choice of the gaussian simulation method
#define PREF_BEST 5
#define PREF_NONE 0
#define PREF_NUGGET 1
#define PREF_PENALTY Nothing
#define PREF_FACTOR Nothing
// larger than Pref * Nothing
#define LOC_PREFB_BEST 9999
#define LOC_PREF_NONE -10000


//////////////////////////////////////////////////////////////////////
// the different levels of printing

#define PL_IMPORTANT 1 
#define PL_SUBIMPORTANT 2
#define PL_RECURSIVE 3
#define PL_REC_DETAILS 4 
#define PL_STRUCTURE 5 // see also initNerror.ERROROUTOFMETHOD
#define PL_ERRORS  6 // only those that are caught internally
#define PL_COV_STRUCTURE 7
#define PL_DIRECT_SEQU 8
#define PL_DETAILS 9
#define PL_SUBDETAILS 10

///////////////////////////////////////////////////////////////////////
// Kappas
#define SIZE_NOT_DETERMINED 0 // do not change the value !

//////////////////////////////////////////////////////////////////////
// the different types of parameters
typedef enum sortsofparam { // never change ordering; just add new ones !!!!!!
  VARPARAM, SIGNEDVARPARAM, SDPARAM, SIGNEDSDPARAM, // 0..3
  SCALEPARAM, DIAGPARAM, ANISOPARAM, // 4..6
  INTEGERPARAM, ANYPARAM, TRENDPARAM /* unused */, // 7..9
  NUGGETVAR, MIXEDVAR,  MIXEDPARAM, /* unused */  // 10..12
  CRITICALPARAM, IGNOREPARAM, LISTPARAM, 
  DONOTRETURNPARAM
} sortsofparam;
typedef enum sortsofeffect{ // ! always compare with convert.R!
  dettrend, deteffect, fixedtrend, fixedeffect, randomeffect, rvareffect,
  largeeffect, lvareffect,  spaceeffect, spvareffect,
  primitive, simple, remaining, eff_error
} sortsofeffect;

///////////////////////////////////////////////////////////////////////
// BASICS
///////////////////////////////////////////////////////////////////////

typedef enum Types {
  // !!! CHANGE ALSO RF_GLOBALS.R !!!
  TcfType, PosDefType, NegDefType, ProcessType, GaussMethodType, 
  BrMethodType, // change also rf_globals.R if deleted
  PointShapeType, RandomType, ShapeType, TrendType, InterfaceType, 
  UndefinedType, // Bedeutung: C->Type:TypeFct existiert; cov->typus:ungesetzt
  OtherType // always last
} Types;


// posdef \subset negdef  ok 
// methodtype \subset processtype
// posdef, negdef \subset shapetype


typedef enum Methods {
  /* never change the ordering -- at least chech nich Standard in location_rules
     in gauss.cc

     variables: METHODNAMES, Besselupperbound, selectlocal, do_comp_meth,
                initmethod, do_incomp_meth, Standard, LastDefault,
		DirectFst, DirectLst, raw, 
		Allowed, allowed
 */
  CircEmbed,     /* Circulant embedding - must always be the first one!*/  //0
  CircEmbedCutoff,
  CircEmbedIntrinsic,
  TBM,           /* Turning Bands, performed in 2 or 3 dimensions */
  SpectralTBM,   /* Spectral turning bands, only 2 dimensional */
  Direct,        /* directly, by matrix inversion */ // 5
  Sequential,    /* sequential simulation */
  Markov,        /* Gaussian Markov Field -- currently not programmed yet */
  Average,       /* Random spatial averages */
  Nugget,        /* just simulate independent variables */
  RandomCoin,    /* "additive MPP (random coins)" */ // 10
  Hyperplane,    /* by superposition of Poisson hyperplane tessellations; 
		    only in the plane! Not implemented yet*/
  Specific,      /* Specific Methods */
  Nothing,       /* must always be the last of the list of methods for
		    simulation Gaussian random fields;

		    NOTE: Nothing has three meanings:
		    (1) used in preference list whether variogram/covariance fct
		        can be calculated
		    (2) method equals Nothing, if no method is preferred
                    (3) method is not shown, if ErrorMessage is called with
		        Nothing
		 */
  Forbidden      /* must always be the last one */
  
} SimulationMeth;
/*
  NOTE: change: Besserlupperbound in CovFct.cc,
  empty_idx and header in getNset.cc 
  METHODNAMES, InitModelList and ErrorMessage in initNerror.cc
  addXXX and createmodel in getNset
*/

#define BrOriginal (Forbidden + 1)
#define BrShifted (Forbidden + 2)
#define BrMixed (Forbidden + 3) // M3
#define BrHat (Forbidden + 4)
#define BrCheck (Forbidden + 5)
#define TotalNumberOfMethods (BrCheck + 1)

typedef enum BRmethod {
  BRorig,
  BRshift,
  BRmixed
} BRmethod;


#define HUETCHEN_OWN_GRIDSIZE 2
typedef signed char ext_bool; // PARAM_DEP, SUBMODEL_DEP, false, true
//                               HUETCHEN_OWN_GRIDSIZE 
typedef char domain_type;
typedef char isotropy_type;

typedef enum TaylorCoeff {TaylorConst, TaylorPow, TaylorExpConst, 
			  TaylorExpPow} TaylorCoeff;
// Const * x^Pow * e^{ExpConst * x^ExpPow} //
// Falls finite range, so wird in TaylorConst der range abgelegt


///////////////////////////////////////////////////////////////////////
// USER CHANGABLE PARAMETERS
///////////////////////////////////////////////////////////////////////
// others

typedef enum units_enum {
  units_none, units_km, units_miles, units_time, units_user
} units_enum;
#define nr_units (units_user + 1)

typedef enum coord_sys_enum {
  coord_auto, cartesian, earth // earth always last
} coord_sys_enum;
#define nr_coord_sys (earth + 1)

#define MAXINVERSIONS 2
typedef enum modes {careless, sloppy, easygoing, normal, precise, 
		    pedantic, neurotic} modes;
#define nr_modes (neurotic + 1)
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
		*/
    detailed_output;
 
  int  mode, /* hightailing, fast, normal, save, pedantic */
    Rprintlevel,
    Cprintlevel; /* 
		   see convert.T PL.* and RF.h PL_*
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
    every,  // every `every' line is announced if every>0

    // bool aniso;  // currently cannot be changed by the user !!

  // siehe InternalCov.cc line 350, simu.cc line 394
  /* expand scale immediately to aniso;
		 this allows higher flexibility in the choice of 
		 simulation methods (tbm2 does not like dimension reduction),
		 but it is slower
	      */
    matrix_inversion[MAXINVERSIONS], // 0:cholesky, 1:QR, 2:SVD; negativ
    seed;
 
  double gridtolerance, matrixtolerance, exactness;

} general_param;

typedef struct gauss_param{
  double stationary_only, // logical + NA
    approx_zero;
  bool paired, loggauss;
  int direct_bestvariables;
} gauss_param;

typedef struct krige_param{
  char method;
  bool cholesky, ret_variance,  fillall;
  int locmaxn,
    locsplitn[3], // 0:when splitting is done; 1:min pts in a neighbourhood ; 2:max pts when still neighbourhoods are included
    locsplitfactor;
} krige_param;

typedef struct ce_param {
  bool force, useprimes, dependent;
  char strategy;
  int trials, maxgridsize;
  double maxmem, tol_re, tol_im, mmin[MAXCEDIM],
    approx_grid_step;
} ce_param;


typedef enum InversionMethod { Cholesky, SVD, NoFurtherInversionMethod } 
InversionMethod;
typedef struct direct_param {
  InversionMethod inversionmethod;
  double svdtolerance;
  int maxvariables;
} direct_param;

typedef struct nugget_param {
  double tol;
} nugget_param;

typedef struct sequ_param{
  int max, back, initial;
} sequ_param;


typedef struct special_param {
   int multcopies;
} special_param;

typedef struct spectral_param {
  bool grid, ergodic;
  double prop_factor, sigma;
  int lines[MAXTBMSPDIM];
} spectral_param;

 
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



typedef struct markov_param{
  int neighbours, cyclic, maxmem;
  double precision;
} markov_param;

typedef struct ave_param {
} ave_param;


typedef struct mpp_param {
  int  n_estim_E;
  double intensity[MAXMPPDIM], // intensity factor for e.g. unif_initu
    about_zero;
} mpp_param;	

#define HYPER_UNIFORM 0   
#define HYPER_FRECHET 1
#define HYPER_BERNOULLI 2
typedef struct hyper_param {
  int superpos, maxlines, mar_distr;
  double mar_param;
} hyper_param;


typedef struct extremes_param {
  int maxpoints, check_every, flat, min_n_zhou, max_n_zhou, mcmc_zhou;
  double standardmax, GEV_xi, density_ratio, eps_zhou;
} extremes_param;

typedef struct br_param {
  int BRmaxmem, BRvertnumber, BRoptimmaxpoints, BRoptim, deltaAM;
  double BRmeshsize, BRoptimtol, variobound, corr_factor;
} br_param;

typedef struct distr_param{
  double safety, minsteplen, innermin, outermax;
  int maxsteps, parts, maxit, mcmc_n, repetitions;
} distr_param;



#define FLAT_UNDETERMINED -1
#define GRIDEXPAND_AVOID -1

typedef struct fit_param{
  double bin_dist_factor, upperbound_scale_factor, lowerbound_scale_factor, 
    lowerbound_scale_LS_factor, upperbound_var_factor, lowerbound_var_factor, 
    lowerbound_sill, scale_max_relative_factor, minbounddistance, 
    minboundreldist, minmixedvar, maxmixedvar, 
    solvesigma, BC_lambdaLB, BC_lambdaUB, 
    sill, usespam, scale_ratio,  // logical + NA
    min_diag; //
  int approximate_functioncalls, bins, nphi, ntheta, ntime, critical, 
    n_crit, locmaxn,
    locsplitn[3], // 0:when splitting is done; 1:min pts in a neighbourhood ; 2:max pts when still neighbourhoods are included
    locsplitfactor, smalldataset,
    optimiser, algorithm, optim_var_estim, likelihood;
  bool refine_onborder, use_naturalscaling, onlyuser, 
    split, reoptimise, ratiotest_approx, split_refined,
    cross_refit;
  char lengthshortname; // 1..255
} fit_param;


typedef struct empvario_param{
  double phi0, theta0, tol;
  bool pseudovariogram, fft;
} empvario_param;

typedef struct gui_param{
  bool alwaysSimulate;
  int method, size[2];
} gui_param;

typedef struct graphics_param{
  bool always_close;
  double height;
  int PL, increase_upto[2];
} graphics_param;


typedef struct registers_param {
   int keynr, interpolregister, condregister, errregister, guiregister;
} registers_param;


typedef struct internal_param{
  // if changed, CHANGE ALSO RestWarnings in 'userinterfaces.cc';
  bool warn_oldstyle, warn_newstyle, warn_Aniso, warn_ambiguous, 
    warn_normal_mode, warn_mode, warn_colour, stored_init, warn_scale,
    warn_coordinates, warn_on_grid
    ;// 
} internal_param;


typedef struct coords_param{
  double xyz_notation;
  coord_sys_enum coord_system;
  char newunits[MAXUNITS][MAXUNITSCHAR], curunits[MAXUNITS][MAXUNITSCHAR],
    varunits[MAXUNITS][MAXUNITSCHAR];
} coords_param;

typedef double *coord_type[MAXSIMUDIM];
typedef struct location_type {
  int 
    timespacedim,      /* total dimension, including time and space */
    length[MAXSIMUDIM], /* grid==true: number of pixels of each direction 
		       grid==false: length[0]=number of spatial points in total
		     */
    lx, ly, /*
	      the total number of locations of the user's x vector, i.e.
	      * except for grids where lx=3
              * physical length = xdimOZ * lx * (lx - 1) /2 for distances
	      * physical length = xdimOZ * lx else
	     */
    spatialdim, 
    xdimOZ; // physical xdim, without time
  long spatialtotalpoints,
    totalpoints;      /* number of pixels in total;
			    == lx if grid==false
			    == lengthx * lengthy * lengthz if grid==true
			 */
  bool grid,  /* simulation on a grid required by user? */
    delete_x, /* is loc->x only a reference? -- used in MLE */
    distances, /* are distances given, and not locations ? -- only for matrices*/
    Time;             /* is a time component given ? */

  coord_type xgr, ygr;      /* the coordinates */  

  double *x,   /* sortiert 1. Koord. 1. Punkt, 2. Koord 1. Punkt, ... */
    *y, /* only used by Covariance function evualation */
    T[3],      /* gridtriple definition for the time components --
		  doubled in last coordinate of *x
	       */
    *caniso;   /* only used for grid. Cummulates the 
		  anisotropy matrices when initS is called
		*/

  // domain_type domain;  
 
  // int CovMatrixCol, CovMatrixRow, CovMatrixTotal;
  int i_row, i_col, cani_ncol, cani_nrow;

} location_type;

typedef struct globalparam{
  general_param general;
  gauss_param gauss;
  krige_param krige;
  ce_param ce;
  spectral_param spectral;
  tbm_param tbm;
  direct_param direct;
  sequ_param sequ;
  markov_param markov;
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
} globalparam;
extern globalparam GLOBAL;

//typedef struct globalorig{
//  bool set;
//  globalparam param;
//} globalorig;
//extern globalorig GLOBALORIG;



extern bool NAOK_RANGE;



///////////////////////////////////////////////////////////////////////
// GENERAL PARAMETERS FOR THE SIMULAIONMETHODS
///////////////////////////////////////////////////////////////////////

#define Nick(COV) (CovList[(COV)->nr].nick)
#define NAME(COV) CovList[(COV)->nr].name
#define KNAME(NAME) CovList[cov->nr].kappanames[NAME]




#define SET_DESTRUCT(A)\
  assert(meth->S==NULL && meth->destruct==NULL);\
  meth->destruct = A;

typedef struct cov_model cov_model;
typedef struct storage storage;
typedef int (*structfct)(cov_model *cov, cov_model **newcov);
typedef int (*initfct)(cov_model *cov, storage *s);
typedef void (*dofct)(cov_model *cov, storage *s);
typedef void (*do_random_fct)(cov_model *cov, double *v);
typedef void (*param_set_fct)(cov_model *to, cov_model *from, int variant);
typedef bool (*type_fct)(Types required, cov_model *cov);

//typedef double static_grid[MAXSIMUDIM][3];
typedef enum matrix_type {TypeMiso, TypeMdiag, TypeMtimesepproj, TypeMtimesep,
			  TypeMproj, 
			  TypeMany} matrix_type;


typedef struct simu_type {
  bool active, pair;   /* has the init_procedure been called successfully? */
  int expected_number_simu;
} simu_type;


extern const char *METHODNAMES[Forbidden+1], *ISONAMES[LAST_ISO + 1], 
  *STATNAMES[LAST_STAT + 1], *ROLENAMES[ROLE_LAST+1], *TYPENAMES[OtherType + 1],
  *CAT_TYPENAMES[OtherType + 1], *REGNAMES[MODEL_MAX+1], 
  *MONOTONE_NAMES[BERNSTEIN + 1 - MISMATCH] ,
  *MODENAMES[nr_modes], *UNITS_NAMES[nr_units], 
  *COORD_SYS_NAMES[nr_coord_sys];



//extern char BRMETHODNAMES[][METHODMAXCHAR];
//extern char PROCESSESNAMES[ROLE_LAST + 1];

SEXP GetMethodInfo(cov_model *cov, 
		   bool ignore_active, int depth, int max, long*mem);




///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////


/* 
   different definitions for the "dimension" of the field is used
      dim : the formal dimension (given by the user)
      effectivedim: the dimension used in the evaluation of the covariance fct
                    E.G. for isotropic models it is 1
      reduceddim : *  dim of lin. independent subspace for ISOTROPIC models;
                   *  1 + dim of lin. indep. spatial subspace for SPACEISOTROPIC
                   *  dim for ANISOTROPIC covariance models 
*/


typedef int pref_type[Forbidden+1];
typedef int pref_shorttype[Nothing + 1];
typedef struct range_type {
  double min[MAXPARAM];
  double max[MAXPARAM];
  bool openmin[MAXPARAM];
  bool openmax[MAXPARAM];
  double pmin[MAXPARAM];
  double pmax[MAXPARAM];
} range_type;

typedef double *param_type[MAXPARAM];
typedef struct sexp_type {
  bool Delete;
  SEXP sexp;
} sexp_type;
typedef struct listoftype {
  bool deletelist;
  double *p[MAXELEMENTS];
  int ncol[MAXELEMENTS], // in case of random parameter, they might be given
    nrow[MAXELEMENTS];   // although p[.] is NULL; latter set in CHECK
} listoftype;
typedef double *naptr_type[MAX_NA];
typedef cov_model *covptr_type[MAX_NA];
typedef int *elptr_type[MAX_MLE_ELMNTS];
// typedef double *internal_type[MAXINTERNALPARAM];

typedef double (*spectral_density)(double *, cov_model *cov); 
typedef struct spectral_storage {
  double phistep2d, phi2d, prop_factor;
  bool grid, ergodic; 
} spectral_storage;


typedef enum tri {TriFalse, TriTrue, TriMaxFalse, TriBeyond} tri;
/* Tri hat nichts mehr zu bedeuten;
   Wird nur im Zusammenhang mit struct-Eintrag anyNAscaleup und 
   anyNAdown verwendet
   
   
   False:
   True:
   MaxFalse:
   Beyond:
*/

typedef enum user_given_type {ug_explicit, ug_implicit, ug_internal}
  user_given_type;

extern cov_model *KEY[MODEL_MAX + 1];
extern int MEM_NAS[MODEL_MAX+1];


#define MAXLOCALINSTANCES 3
typedef struct localinfotype{
  int instances;
  int msg[MAXLOCALINSTANCES];
  double value[MAXLOCALINSTANCES];
} localinfotype;



typedef void (*rangefct)(cov_model *cov, range_type* ra);
typedef int (*checkfct)(cov_model *cov); 
// typedef double (*inverse_fct)(double *v, cov_model *cov, double *x); /* parameters,; natural scaling */
typedef void (*covfct)(double *, cov_model*, double*); /* h, cov, result */ 
typedef void (*nonstat_covfct)(double *, double*, 
			      cov_model*,  double*); /* x,y, cov, result */
typedef void (*nonstat_inv)(double *, cov_model*, double*, 
			       double*); /* fctvalue cov, inv.x inv.y */
typedef void (*logfct)(double *, cov_model*, double* v, double* sign);
typedef void (*nonstat_logfct)(double *, double*, 
			      cov_model*,  double* v, double* sign); 
typedef void (*aux_covfct)(double *, double*, double, 
			      cov_model*,  double*); /* x,y, Aux, cov, result */
typedef void (*return_fct)(cov_model*, double*); /* cov, result */ 
typedef char (*ext_bool_ret_fct)(cov_model*); /* covt */ 
typedef void (*returnX_fct)(cov_model*, int *, int, double*);/* cov,..,result */
typedef void (*getlocalparam)(cov_model *, localinfotype *);
typedef bool (*altlocalparam)(cov_model *);
typedef void (*minmaxfct)(cov_model *, double *);


typedef void (*spectral_do)(cov_model *, storage *, double *);

typedef void (*draw_random) (cov_model *cov, double *random);
typedef double (*log_mixdens)(double *x, double logV, cov_model *cov);

typedef void (*sd_fct)(storage *s, cov_model *cov);
		       				   
typedef int (*hyper_pp_fct)(double, double*, double*, cov_model *, bool, 
			    double**, double**, double**);
typedef void (*size_fct)(int i, cov_model *cov, int *nr, int *nc);
typedef sortsofparam  (*paramtype_fct)(int k, int row, int col);

// here all the different method for simulating a RF of a specified covariance
// function are stored


typedef struct cov_fct {  
  char name[MAXCHAR], nick[MAXCHAR];
  char kappas, /* number of parameters  */
    minsub, maxsub; // ==0
  domain_type domain;
  isotropy_type isotropy;
  int vdim, maxdim, maxmoments, Monotone;  //
  char kappanames[MAXPARAM][PARAMMAXCHAR],
    subnames[MAXSUB][PARAMMAXCHAR];
 
  SEXPTYPE kappatype[MAXPARAM];
  Types kappaParamType[MAXPARAM];
  size_fct kappasize;
  paramtype_fct paramtype;

  rangefct range;
  checkfct check;
  int implemented[Forbidden], Specific;
  ext_bool finiterange;
  bool internal, 
    subintern[MAXSUB];   // do any subnames match exactly a parameter name?
  
  pref_shorttype pref;

  covfct cov, D, D2, D3, D4, tbm2, inverse, nabla, hess, random, logD /* Vtlgen */;
  logfct log; 
  nonstat_covfct nonstat_cov, nonstat_D, nonstat_random;
  nonstat_inv nonstat_inverse, nonstat_loginverse, nonstat_inverse_D;
  nonstat_logfct nonstatlog;
  int F_derivs, RS_derivs;
  param_set_fct param_set;
  aux_covfct aux_cov; // complicated cov-model that can be used only
  //                     as submodels and that needs an auxiliary argument
  //                     for the evaluation
  getlocalparam coinit, ieinit; // set within primitives
  altlocalparam alternative; //getparam: guess for good local 
    // param (cutoff, intrinsic); alternative gives alternative in a 
    // second or third try (used by co and Stein)

  spectral_do spectral;

  draw_random drawmix;
  log_mixdens logmixdens;
 
  structfct Struct;
  initfct Init;
  dofct Do;
  do_random_fct DoRandom;
  minmaxfct minmaxeigenvalue;
 
  hyper_pp_fct hyperplane;      /* hyperplane tessellations */        

  bool primitive;
  Types Type;
  type_fct TypeFct;

  int density;
  
  double Taylor[MAXTAYLOR][TaylorPow + 1], 
    Tail[MAXTAYLOR][TaylorExpPow + 1]; 
  int TaylorN, TailN;

  return_fct covariance, covmatrix, inversecovmatrix, variogram,
    pseudovariogram;
  ext_bool_ret_fct is_covariance, is_covmatrix, is_inversecovmatrix, 
    is_variogram, is_pseudovariogram, is_selectedcovmatrix;
  returnX_fct selectedcovmatrix;
} cov_fct;

extern int currentNrCov, currentRegister;
typedef struct cov_fct cov_fct;
extern cov_fct *CovList;
extern int GENERALISEDCAUCHY, STABLE, BROWNIAN, CAUCHY,
  GAUSS, NUGGET, PLUS, MLEPLUS, MIXEDEFFECT, BALL, ECF, MULT,
  DISTRIBUTION,  DETERM_DISTR, GAUSS_DISTR, SETPARAM, SET_DISTR,
  COVMATRIX, RFGET, COVFCTN,
  DOLLAR, POWER_DOLLAR,  MLE_ENDCOV,  SPLIT, 
  ISO2ISO,SP2SP, SP2ISO, S2ISO, S2SP, S2S, SId, EARTHKM2CART, EARTHMILES2CART,
  NATSC_USER, NATSC_INTERN, TREND, CONSTANT,
  LOC, SET_DISTR, SCALESPHERICAL, TRUNCSUPPORT, BROWNRESNICK, 
  STROKORB_MONO, STROKORB_BALL_INNER, POLYGON, RECTANGULAR, MULT_INVERSE,
  SHAPESTP, SHAPEAVE, SPHERICAL, UNIF, MPPPLUS, PTS_GIVEN_SHAPE,
  STATIONARY_SHAPE, STANDARD_SHAPE,
  VARIOGRAM_CALL, 
  BRORIGINAL_USER, BRMIXED_USER, BRSHIFTED_USER,
  ARCSQRT_DISTR,
  
  BRORIGINAL_INTERN, BRMIXED_INTERN, BRSHIFTED_INTERN,
  MISSING_COV, NULL_MODEL, TBM_OP, USER, 
  DOLLAR_PROC, PLUS_PROC,
  BINARYPROC, BROWNRESNICKPROC, GAUSSPROC, POISSONPROC,
  SCHLATHERPROC, SMITHPROC, CHI2PROC,
  NUGGET_USER,  NUGGET_INTERN,
  CIRCEMBED, CUTOFF, STEIN, SPECTRAL_PROC_USER, SPECTRAL_PROC_INTERN, 
  DIRECT, SEQUENTIAL, 
  AVERAGE_USER, AVERAGE_INTERN, HYPERPLANE_USER, HYPERPLANE_INTERN,
  RANDOMCOIN_USER, CE_CUTOFFPROC_USER, CE_CUTOFFPROC_INTERN, 
  CE_INTRINPROC_USER, CE_INTRINPROC_INTERN, TBM_PROC_USER, TBM_PROC_INTERN, 
  SPECIFIC, SELECT,
  BRSHIFTED, BRMIXED, BRORIGINAL, EXTREMALGAUSSIAN, RANDOMSIGN,
  TBM2NR, VECTOR;
extern char OP_SIGN[][3];
extern char CovNames[MAXNRCOVFCTS][MAXCHAR],
  CovNickNames[MAXNRCOVFCTS][MAXCHAR];
extern char InternalName[];

// extern int etc, if other are enabled to use these functions?
int IncludePrim(const char *name, Types type, int kappas,  
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range, int maxdim, 
		ext_bool finiterange, int monotonicity);
int IncludePrim(const char *name, Types type, int kappas,  size_fct kappasize,
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range, int maxdim, 
		ext_bool finiterange, int monotonicity);
int IncludePrim(const char *name, Types type, int kappas,  
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range, int vdim, 
		int maxdim, ext_bool finiterange, int monotonicity);
int IncludePrim(const char *name,  Types type, int kappas,  size_fct kappasize,
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range, int vdim,
		int maxdim, ext_bool finiterange, int monotonicity);

int IncludePrim(const char *name,Types type,  int kappas, 
		 domain_type domain, isotropy_type isotropy,	
		 checkfct check, rangefct range, pref_type pref,
		 int vdim, int maxdim, ext_bool finiterange, int monotonicity);
int IncludePrim(const char *name, Types type, int kappas, size_fct kappasize,
		 domain_type domain, isotropy_type isotropy,	
		 checkfct check, rangefct range, pref_type pref,
		int vdim, int maxdim, ext_bool finiterange, int monotonicity);

void make_internal();


int IncludeModel(const char *name, Types type, char minsub, char maxsub,
		 int kappas,
		  size_fct kappasize,
		  domain_type domain, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref, 
		  bool internal, int vdim, int maxdim, ext_bool finiterange,
		 int monotonicity);
int IncludeModel(const char *name, Types type, char minsub, char maxsub, 
		 int kappas, 
		  domain_type domain, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref, 
		 int maxdim, ext_bool finiterange, int monotonicity);
int CopyModel(const char *name, int which);
int CopyModel(const char *name, int which, Types type);
int CopyModel(const char *name, int which, checkfct check);
void nickname(const char *nick);
void addCov(covfct cf, covfct D, covfct inverse);
void addCov(covfct cf, covfct D, covfct D2, covfct inverse);
void addCov(covfct cf, covfct D, covfct D2, nonstat_inv inverse);
void addCov(covfct cf, covfct D, covfct D2, 
	    covfct inverse, nonstat_inv nonstat_inverse);
void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse);
void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse,  nonstat_inv nonstat_inverse);
void addCov( int F_derivs, covfct cf, covfct D, covfct inverse);
void addCov(int F_derivs, covfct cf, covfct D, covfct D2, covfct inverse,
	    nonstat_inv nonstatinverse);
void addCov( int F_derivs, covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse);
void addCov(nonstat_covfct cf);
void addCov( int F_derivs , nonstat_covfct cf);
void addCov(aux_covfct auxcf);
void addCov(covfct distrD, covfct logdistrD, nonstat_inv Dinverse,
	    covfct distrP, nonstat_covfct distrP2sided,
	    covfct distrQ, covfct distrR, nonstat_covfct distrR2sided);

void addlogCov(logfct log, nonstat_logfct nonstatlog, 
	       nonstat_inv lognonstat_inverse);
void addlogCov(logfct log);
void addlogCov(nonstat_logfct nonstatlog);


void addlogD(covfct logD);

void nablahess(covfct nabla, covfct hess);
void addSpaceD(covfct spaceD);
void addLocal(getlocalparam coinit, getlocalparam ieinit);
void addCallLocal(altlocalparam alt);
//void addTBM(covfct tbm2, covfct spaceD);
int addTBM(covfct tbm2);
void addTBM(covfct tbm2, initfct Init, spectral_do spectral);
void addTBM(initfct Init, spectral_do spectral);
void addHyper(hyper_pp_fct hyper_pp);
void addMarkov(int*variable);
void addSpecific(int cov);
void RandomShape(int, structfct Struct, initfct init, dofct Do,
		 do_random_fct DoRandom, 
		 bool average, bool randomcoin, bool specific);
void RandomShape(int, structfct Struct, initfct Init, dofct Do, 
		 bool average, bool randomcoin, bool specific);
void RandomShape(int, structfct Struct, initfct Init, do_random_fct DoRandom);
void RandomShape(int, structfct Struct, initfct init, dofct Do);
void RandomShape(structfct Struct, bool average);
void RandomShape(structfct Struct);
void RandomShape(int, initfct init, dofct Do, bool average);
void RandomShape(int, initfct init, dofct Do);
void addSpecial(minmaxfct minmaxeigen);
void addGaussMixture(draw_random drawmix, log_mixdens logmixdens);

int addFurtherCov(covfct cf, covfct D);
int addFurtherCov(covfct cf, covfct D, covfct D2);
int addFurtherCov(nonstat_covfct cf);
void add_paramtype(paramtype_fct paramtype);
void addReturns(return_fct Covariance, ext_bool_ret_fct isCovariance, 
		return_fct CovMatrix, ext_bool_ret_fct isCovMatrix,
		return_fct InverseCovMatrix,ext_bool_ret_fct isInverseCovMatrix,
		return_fct Variogram, ext_bool_ret_fct isVariogram,
		return_fct PseudoVariogram, ext_bool_ret_fct isPseudoVariogram,
		returnX_fct SelectedCovMatrix, 
		ext_bool_ret_fct isSelectedCovMatrix);
void addTypeFct(type_fct TypeFct);


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
// GLBBAL PARAMETERS AND STRUCTURE
///////////////////////////////////////////////////////////////////////
#define DECISION_TRUE 1
#define DECISION_FALSE 0
// #define DECISION_CASESPEC -1
#define MAXERRORSTRING (100 + (MAXNRCOVFCTS) * (MAXCHAR+1))
#define nErrorLoc 1000

extern int PL, NS;
extern double GENERAL_PRECISION;
extern char ERRORSTRING[MAXERRORSTRING], ERRORSTRING_OK[MAXERRORSTRING],
  ERRORSTRING_WRONG[MAXERRORSTRING], ERROR_LOC[nErrorLoc];
extern int ERRORMODELNUMBER;
void errorMSG(int error, char* EM);
void errorMSG(int err, char* m, int len);

//void ExitInit(int err, bool *);
//void EnterInit(char *name);
void matrixrotat(double *paramaniso, int col, int row,
		 int *truedim, double* aniso);
SEXP TooLarge(int *n, int l); /* if vector is larger than given size */


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
  complex *gauss1, *gauss2;
  bool positivedefinite,
    stop,
    new_simulation_next,
    cur_call_odd,
     dependent; // eigentlich braucht es nicht waehrend der initialisierung
    // festgelegt zu werden. Ist aber wesentlich einfacher zu handhaben,
    // da sonst bei internal_dosimulate die parameter alle RFparameter alle
    // nochmals gesetzt werden muessen
  FFT_storage FFT;
} CE_storage;

typedef struct localCE_storage {
    void *correction;
} localCE_storage;

typedef struct ce_approx_storage {
  int *idx;
} ce_approx_storage;
unsigned long NiceFFTNumber(unsigned long nn);
int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *S);
int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *S);
void FFT_destruct(FFT_storage *S);
void FFT_NULL(FFT_storage *S);


typedef struct mixed_storage {
  // bool genuine_dim[MAX362DIM];
  double *Xb, *mixedcov;
  bool initialized;
  // hiermit kann unterschieden werden, ob random effects nur
  // einmal simuliert werden, wenn RFsimulate(n > 1)
} mixed_storage;

typedef struct trend_storage {
  double *x;
  int *xi;
  double *evalplane;
  int *powmatrix;
} trend_storage;

// see tbm.cc			  
typedef struct TBM_storage {
  // bool genuine_dim[MAX362DIM];
  int ce_dim, simuspatialdim, method, spatialtotalpts, err;
  double center[MAXTBMSPDIM],  
    linesimuscale, linesimufactor, xline_length; 
} TBM_storage;

// see spectral.cc

void metropolis(cov_model *cov, storage *S, double *x);
int search_metropolis(cov_model *cov, storage *S);
// see direct.cc
typedef struct direct_storage {
  InversionMethod method;
  double *U, *G;
} direct_storage;


// see sequential.cc
typedef struct sequential_storage {
  int back, totpnts, spatialpnts, ntime, initial;
  double *U22, *G, *MuT, *U11,  *Cov21, *Inv22;
    res_type *res0;
} sequential_storage;


// nugget
typedef struct nugget_storage {
    //double sqrtnugget;
  bool simple, simugrid;
  int *pos, reduced_dim[MAXNUGGDIM], prod_dim[MAXNUGGDIM + 1];
  res_type *red_field;
} nugget_storage;


// dummy version, hyperplane
typedef struct hyper_storage{
  double rx[MAXHYPERDIM], center[MAXHYPERDIM], radius;
  hyper_pp_fct hyperplane;
} hyper_storage;


typedef struct plus_storage{
  cov_model *keys[MAXSUB];
  int struct_err[MAXSUB];
} plus_storage;


typedef struct integral_type {
  double E, M2, positive;
} integral_type;


typedef struct BR_storage {
  int *locindex,
    trendlen,
    *loc2mem,
    *mem2loc,
    memcounter,
    vertnumber,
    idx,
    maxidx,
    next_am_check,
    zeropos[MAXSUB],
    **countvector
   ;
  double **trend, *newx,
    **areamatrix,
    *shiftedloc,
    minradius,
    *logvertnumber,
    *lowerbounds[MAXSUB],
    radii[MAXSUB], thresholds[MAXSUB],
    *locmin, *locmax, *loccentre, // only dummy variable!  
    *suppmin, *suppmax; // only dummy variable!  
  cov_model *vario, *submodel, 
    *sub[MAXSUB]
    ;
} BR_storage;


typedef struct get_storage {
  cov_model *orig, *get_cov;
  int param_nr, size, vdim2[2], *idx;
  bool all;
} get_storage;


typedef struct pgs_storage {
  // urpsprunglich nur fuer pts_given_shape; jedoch allgemein
  // fuer shape functionen
  bool flat, estimated_zhou_c, logmean; 
  coord_type xgr;
  double old_zhou, // for mcmc only
  totalmass, // (inverser) Normierungsfaktor, um Raeumliche Funktion
  //            zu einer Dichte zu machen
    currentthreshold, log_density, globalmin,
    *single, *total, *halfstepvector,  // global
    *supportmin, *supportmax, *supportcentre, // global inkl. dompp
    *own_grid_start, *own_grid_step,
    *own_grid_length, // only for HUETCHEN_OWN_GRIDSIZE
    *v, *x, *y,  // local
    *localmin, *localmax, // local
    *minmean, *maxmean, // standard_shape
    *xstart, *inc,// local dompp
    intensity;
  int 
    own_grid_cumsum[MAXMPPDIM],
    size, *len, // global
     *pos, // local
    *min, *max,  *gridlen, *start, *end, *delta, *nx;// local dompp
  // not used in pgs, but in variogramAndCo.cc
  int *endy, *startny, *ptrcol, *ptrrow;
  double *C0x, *C0y, *z, *cross, **Val;
  param_set_fct param_set;
  cov_model *cov;

  long double sum_zhou_c, sq_zhou_c;
  long int n_zhou_c;
  double zhou_c; // c in oesting, s, zhoy
} pgs_storage;

typedef struct set_storage {
  cov_model *remote;
  param_set_fct set;
  //  void **valueRemote,
  //    **valueLocal;
  int variant;
  //     *bytes; // to be transfered
} set_storage;

#define MAXDIM_POLY 2
typedef struct vertex {
  double x[2];				// vertex coordinates
} vertex;
typedef struct edge {
  double u[2];				// normal vector
  double p;	       			// distance from the origin
} edge;
typedef struct polygon {
	int n;			// number of vertices = number of edges
        vertex *v;		// array of vertices
        edge *e;		// array of edges
	double box0[2];		// coordinates of the left lower vertex of smallest box containing the polygon
	double box1[2];		// coordinates of the right upper vertex of smallest box containing the polygon
} polygon;
typedef struct polygon_storage {
  polygon *P;
  double **vdual;
  vertex *vprim;
  int n_vdual, n_vertex, n_v;
} polygon_storage;




typedef struct rect_storage {
  double inner, inner_const, inner_pow,
    outer, outer_const, outer_pow, outer_pow_const,
    step, *value, *weight, 
    *tmp_weight, *right_endpoint, *ysort, *z;
    int nstep, *squeezed_dim, *assign, *i, tmp_n;
} rect_storage;



typedef struct dollar_storage {
  double *z, *z2, *y, *save_aniso, *inv_aniso;
  int *cumsum, *nx, *total, *len;
} dollar_storage;

typedef struct gatter_storage {
  double *z, *zsys;
} gatter_storage;

typedef struct biwm_storage {
  bool nudiag_given, cdiag_given;
  double a[3], lg[3], aorig[3], nunew[3],
    scale[4], gamma[4], c[4];
} biwm_storage;


typedef struct inv_storage {
  double *v, *wert;
} inv__storage;


typedef struct spec_properties {
  spectral_density density;
  double sigma, E[MAXTBMSPDIM];
  int nmetro;
  double sub_sd_cum[MAXSUB];
} spec_properties;


typedef struct mpp_properties {
  double 
    unnormedmass, // RRrectangle: mass of function that is at least
    //                as large as the given one
   // or of f / g (Oesting, Sch    lather, Zhou), depending on the 
  // function (SHAPE_FCT)
    maxheights[MAXMPPVDIM], // maximum of f resp. of \d F   
    *mM, *mMplus // = int f^k \D \lambda, falls keine Verteilungsfamilie
    //                      und falls f kein stochastischer Prozess
    //       = \EE X^k = \int x^k \F(x) falls Verteilungsfamilie oder falls
    //             \EE X(0)^k, falls stochastischer Prozess
    // refradius , refsd,
    // totalmass 
    ;
  int // methnr,
    moments;
  //bool loc_done;
} mpp_properties;


typedef struct storage { // cov_storage, do_storage, init_storage
  // wird in Struct initialisiert, sofern INIT aufgerufen wird;
  // abgespeichert wird es immer im aktuell aufrufenden (fuer next oder key)
  bool check; // used in biWM, BiGneiting
  spec_properties spec;       // used in init
  spectral_storage Sspectral; // used in do
} storage;

///////////////////////////////////////////////////////////////////////
// POSITIONS OF X-VALUES
///////////////////////////////////////////////////////////////////////
#define XSTART 0
#define XSTEP 1
#define XLENGTH 2

void Transform2NoGrid(cov_model *cov, bool timesep, int gridexpand);
void Transform2NoGrid(cov_model *cov, double **xx);
void Transform2NoGrid(cov_model *cov, double **xx, double **yy);
void Transform2NoGridExt(cov_model *cov, bool timesep, int gridexpand, 
			 double ** grani, double **SpaceTime,
			 double **caniso, int *nrow, int *ncol,
			 bool *Time, bool *grid, int *newdim, bool takeX);


///////////////////////////////////////////////////////////////////////
// UNAUFGERAEUMT:
///////////////////////////////////////////////////////////////////////

//EXTERN void F77_NAME(ave1f)(double *x, double *u,
//			    double *h, double *Mz,
//			    double *r, int *d, double *res);

#define LENERRMSG 2000
extern char MSG[LENERRMSG];
extern char NEWMSG[LENERRMSG];



#define STOP // assert(false)
#define ERR(X) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); STOP; error(MSG);}
#define XXERR(X) {/* UERR; */ ERRLINE;errorMSG(X, MSG); 	\
 sprintf(NEWMSG, "in `%s' error %d: %s", ERROR_LOC, X, MSG); STOP; error(NEWMSG);}
#define AERR(X) {ERRLINE; PRINTF("AERR: "); errorMSG(X, MSG); 	\
    if (PL<PL_ERRORS) PRINTF("%s%s\n", ERROR_LOC, MSG); assert(false);}
#define MERR(X) { PRINTF("MERR: "); errorMSG(X, MSG, 50);		\
    if (PL<PL_ERRORS) PRINTF("%s%s\n", ERROR_LOC, MSG);}
#define XERR(X) {/* UERR; */ ERRLINE;errorMSG(X, MSG); 	\
    sprintf(NEWMSG, "%s%s", ERROR_LOC, MSG); STOP; error(NEWMSG);}
#define PERR(X) {sprintf(MSG, "%s '%s': %s", ERROR_LOC, param_name, X); STOP;error(MSG);}
#define QERR(X) {sprintf(ERRORSTRING, "'%s' %s", param_name, X); DEBUGINFOERR; return ERRORM;}
#define QERRX(ERR, X) {errorMSG(ERR, MSG); sprintf(ERRORSTRING, "'%s' %s (%s)", param_name, X, MSG); DEBUGINFOERR; return ERRORM;}
#define QERRC(NR,X) {sprintf(ERRORSTRING, "%s '%s': %s", ERROR_LOC, CovList[cov->nr].kappanames[NR], X); DEBUGINFOERR; return ERRORM;}
#define SERR(X) { strcpy(ERRORSTRING, X); /* PRINTF("serr:%s\n", X); */ DEBUGINFOERR; return ERRORM;}
#define SERR1(X,Y) { sprintf(ERRORSTRING, X, Y); /* PRINTF("serr:%s\n", X); */DEBUGINFOERR;  return ERRORM;}
#define SERR2(X, Y, Z) { sprintf(ERRORSTRING, X, Y, Z); /* PRINTF("serr:%s\n", X); */ DEBUGINFOERR; return ERRORM;}
#define SERR3(X, Y, Z, A) { sprintf(ERRORSTRING, X, Y, Z, A); /* PRINTF("serr:%s\n", X); */ DEBUGINFOERR; return ERRORM;}
#define SERR4(X, Y, Z, A, B) { sprintf(ERRORSTRING, X, Y, Z, A, B); /* PRINTF("serr:%s\n", X); */DEBUGINFOERR;  return ERRORM;}
#define SERR5(X, Y, Z, A, B, C) { sprintf(ERRORSTRING, X, Y, Z, A, B, C); /* PRINTF("serr:%s\n", X); */ DEBUGINFOERR; return ERRORM;}
#define SERR6(X, Y, Z, A, B, C, D) { sprintf(ERRORSTRING, X, Y, Z, A, B, C, D); /* PRINTF("serr:%s\n", X); */ DEBUGINFOERR; return ERRORM;}
#define SERR7(X, Y, Z, A, B, C, D, E) { sprintf(ERRORSTRING, X, Y, Z, A, B, C, D, E); /* PRINTF("serr:%s\n", X); */ DEBUGINFOERR; return ERRORM;}
#define GERR(X) { strcpy(ERRORSTRING, X); /* PRINTF("gerr:%s\n", X); */ err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR1(X,Y) { sprintf(ERRORSTRING, X, Y); /* PRINTF("gerr:%s\n", X); */ err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR2(X,Y,Z) { sprintf(ERRORSTRING, X, Y, Z); /* PRINTF("gerr:%s\n", X); */ err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR3(X,Y,Z,A) { sprintf(ERRORSTRING, X, Y, Z, A); /* PRINTF("gerr:%s\n", X); */ err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR4(X,Y,Z,A,B) { sprintf(ERRORSTRING, X, Y, Z, A, B); /* PRINTF("gerr:%s\n", X); */ err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR5(X,Y,Z,A,B,C) { sprintf(ERRORSTRING, X, Y, Z, A, B, C); /* PRINTF("gerr:%s\n", X); */ err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}

extern char BUG_MSG[250];
#define BUG {								\
    sprintf(BUG_MSG, "Severe error occured in function '%s' (file '%s', line %d). Please contact maintainer martin.schlather@math.uni-mannheim.de .", \
	    __FUNCTION__, __FILE__, __LINE__);				\
    error(BUG_MSG);							\
  }									

#define NotProgrammedYet(X) {						\
    { if (strcmp(X, "") == 0)						\
      sprintf(BUG_MSG,							\
	      "function '%s' (file '%s', line %d) not programmed yet.", \
	      __FUNCTION__, __FILE__, __LINE__);			\
    else								\
      sprintf(BUG_MSG, "'%s' in '%s' (file '%s', line %d) not programmed yet.",\
	      X, __FUNCTION__, __FILE__, __LINE__);			\
      error(BUG_MSG);							\
    }\
    }									

#define nNotProgrammedYet(X) {						\
    { if (strcmp(X, "") == 0)						\
      sprintf(BUG_MSG,							\
	      "function '%s' (file '%s', line %d) not programmed yet.", \
	      __FUNCTION__, __FILE__, __LINE__);			\
    else								\
      sprintf(BUG_MSG, "'%s' in '%s' (file '%s', line %d) not programmed yet.",\
	      X, __FUNCTION__, __FILE__, __LINE__);			\
      error(BUG_MSG);							\
    }\
    }									


void contact(const char *msg);

void ErrCov(double *x, cov_model *cov, double *v);
void ErrCovNonstat(double *x, double *y, cov_model *cov, double *v);
void GetNaturalScaling(cov_model *cov, double *natscale);
void analyse_matrix(double *aniso, int row, int col,
		    bool *diag, bool *quasidiag, int *idx,
		    bool *semiseparatelast, bool *separatelast);

void paramcpy(cov_model *current, cov_model *cov, bool freeing,
	      bool allocating, bool copy_lists, bool recursive, 
	      bool copy_mpp);
int covcpy(cov_model **localcov, cov_model *cov);
int covcpy(cov_model **localcov, cov_model *cov, bool copy_lists);
int covcpy(cov_model **localcov, bool sub, cov_model *cov,
	   location_type *prevloc);
int covcpy(cov_model **localcov, bool sub, cov_model *cov, 
	   location_type *prevloc, location_type *ownloc, bool copy_lists);
int covcpy(cov_model **localcov, cov_model *cov,
	   double *x, double *T, int spatialdim, int xdim, int lx, bool Time, 
	   bool grid, bool distances);
int covcpyWithoutRandomParam(cov_model **localcov, cov_model *cov);
int covcpy(cov_model **localcov, bool sub, cov_model *cov, // err
	   location_type *prevloc, location_type *ownloc,
	   bool copy_lists, bool copy_randomparam, 
	   bool allowCopyingInterface);
void Ssetcpy(cov_model *localcov, cov_model *remotecov, cov_model *cov,
	     cov_model *rmt);

int newmodel_covcpy(cov_model **localcov, int model, cov_model *cov);
int newmodel_covcpy(cov_model **localcov, int model, cov_model *cov,
		    double *x, double *y, double *T, int spatialdim,
	            int xdim, int lx, int ly, bool Time, bool grid,
	            bool distances);

void TaylorCopy(cov_model *to, cov_model *from);

int getmodelnr(char *name);

typedef int (*set_local_q_type) (cov_model *next, double a, double d, 
				  double *q);
//void setdefault(cov_model *cov);
void setbackward(cov_model *cov, cov_model *sub);

int getmodelnr(char *name);

void addkappa(int i, const char *n, SEXPTYPE t, Types ParamType);
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
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15, const char* n16, SEXPTYPE t16);

void subnames(const char* n1);
void subnames(const char* n1, const char* n2);
void subnames(const char* n1, const char* n2, const char* n3);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5, const char* n6);

void setAniso(cov_model *cov);

void orderingInt(int *d, int len, int dim, int *pos);
int partial_loc_set(location_type *loc, double *x, double *y,
		    int lx, int ly, bool dist, int xdim, double *T, 
		    bool grid, bool cpy);
//void partial_loc_set(cov_model *cov, double *x, double *y,
//		     int lx, bool dist, int xdim, double *T);
int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdim, int lx, bool Time, bool grid,
	    bool distances,
	    location_type **loc);
int loc_set(double *x, double *y, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdim,
	    int lx, int ly, bool Time, bool grid,
	    bool distances,
	    location_type **Loc);
int loc_set(cov_model *cov, int totalpoints);
int add_y_zero(location_type *loc);


int Match(char *name, const char * List[], int n);
int Match(char *name, name_type List, int n);

extern name_type STANDARDPARAM, STANDARDSUB;

void CMbuild(SEXP model, int tsdim, int xdim, int level,
	     cov_model **Cov, cov_model previous);
void CheckModelInternal(SEXP model, double *x, double *y, double *T, 
			int spatialdim, /* spatial dim only ! */
			int xdim,
			int lx,  int ly,
			bool grid,
			bool distances,
			bool Time, 
			cov_model **Cov);
void CheckModel(SEXP model, SEXP tsdim, SEXP xdim, SEXP domain,
		cov_model **Cov,
		int maxdim);
void pmi(cov_model *cov);
void pmi(cov_model *cov, const char *msg);
void pmi(cov_model *cov, char all);

matrix_type Type(double *m, int nrow, int ncol);
double GetDiameter(location_type *loc);
double GetDiameter(location_type *loc, 
		   double *min, double *max, double *center);
double * matrixmult(double *m1, double *m2, int dim1, int dim2, int dim3);
void addModel(cov_model **pcov, int covnr);
int addUnifModel(cov_model *cov, double radius, cov_model **newmodel);

void addVariable(char *name, double *x, int nrow, int ncol, SEXP env);

//void addModel(cov_model **pcov, int covnr, bool takeoverloc);



#define MULTIPLEMATCHING -2
#define NOMATCHING -1
#define MATCHESINTERNAL -3

#define UNSET -1


void PrintModelList();


void E1(spectral_storage *s, double A, double *e);
void E2(spectral_storage *s, double A, double *e);
void E12(spectral_storage *s, int dim, double A, double *e);
void E3(spectral_storage *s, double A, double *e);
void E(int dim, spectral_storage *s, double A, double *e);

int checkkappas(cov_model *cov);
int checkkappas(cov_model *cov, bool errornull);
//double gaussInt(int d, int xi, double sigma, double R);

void updatepref(cov_model *cov, cov_model *sub);

extern char PREF_FAILURE[100 * Nothing];
void leer(int level);

int get_ranges(cov_model *cov, cov_model **min, cov_model **max,
		cov_model **pmin, cov_model **pmax, 
		cov_model **openmin, cov_model **openmax);
int check_recursive_range(cov_model *cov, bool NAOK);
int internal_DoGauss(cov_model  *key, int nn, res_type *orig_res);

void strcopyN(char *dest, const char *src, int n);
void expandgrid(coord_type xgr, int *len, double **xx, int nrow);
void xtime2x(double *x, int nx, double *T, int len, double **newx, int nrow);
void removeOnly(cov_model **Cov);

SEXP Logi(bool* V, int n, int max) ;
SEXP Num(double* V, int n, int max) ;
SEXP Int(int *V, int n, int max) ;
SEXP Char(const char **V, int n, int max) ;
SEXP Mat(double* V, int row, int col, int max);
SEXP MatInt(int* V, int row, int col, int max) ;
SEXP Array3D(int** V, int depth, int row, int col, int max) ;

SEXP Logi(bool* V, int n) ;
SEXP Num(double* V, int n) ;
SEXP Int(int *V, int n) ;
SEXP Char(const char **V, int n) ;
SEXP Mat(double* V, int row, int col);
SEXP MatInt(int* V, int row, int col) ;
SEXP Array3D(int** V, int depth, int row, int col) ;

SEXP TooLarge(int *n, int l);

//int Checking(cov_model **Cov);
int check2X(cov_model *cov, int tsdim, int tsxdim, Types type, 
	    domain_type domprev, isotropy_type isoprev, 
	    int vdim, int role);
int check2X(cov_model *cov, int tsdim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
	    int *vdim, int role);
int check2X(cov_model *cov, int tsdim, int tsxdim, Types type, 
	    domain_type domprev, isotropy_type isoprev,
	    int vdim0, int vdim1, int role);


void kdefault(cov_model *cov, int i, double v);
void DeleteGatter(cov_model **Cov);

extern pref_type PREF_ALL, PREF_NOTHING;

#define MAXDEFMATRIX 3
//extern double *userdefinedCovMatrix[MAXDEFMATRIX][MAXMAKEEXPLICITE];
// extern int userdefinedCM_RC[MAXDEFMATRIX][MAXMAKEEXPLICITE];
//extern int CovMatrixRow, CovMatrixCol, CovMatrixTotal;
// extern bool CovMatrixIndexActive;
void zeronugget(int* zeronugget);
				
double *EinheitsMatrix(int dim);
				
extern "C" double *OutX, *OutY;

void CovVario(cov_model *cov, bool is_cov, bool pseudo, 
	      double *value);
void CovIntern(int reg, double *x, double *value);
void CovIntern(int reg, double *x, double *y, int lx, int ly, double *value);
void PseudovariogramIntern(int reg, double *x, double *y,
			   int lx, int ly, double *value);
void PseudovariogramIntern(int reg, double *x, double *value);
void CovarianceMatrix(cov_model *cov,double *COV);
void InverseCovMatrix(cov_model *cov, double *v);
void SelectedCovMatrix(cov_model *cov, int *selected /* nr! */, int nsel, 
		       double *v);
void StandardCovariance(cov_model *cov, double *v);
void StandardCovMatrix(cov_model *cov, double *v);
void StandardInverseCovMatrix(cov_model *cov, double *v);
void StandardVariogram(cov_model *cov, double *v);
void StandardPseudoVariogram(cov_model *cov, double *v);
void StandardSelectedCovMatrix(cov_model *cov, int *sel, int nsel, double *v);

// extern int CumIdxMakeExplicite[MAXMAKEEXPLICITE];
double *getAnisoMatrix(cov_model *cov, int *nrow, int *ncol);
double *getAnisoMatrix(cov_model *cov, bool null_if_id, int *nrow, int *ncol);



#define AveMaxDim 10 /* nur technisch ! MAXDIM !*/
#define CoxMaxDim 3  /* nur technisch ?! MAXDIM  ! */
#define StpMaxDim 10  /* nur technisch ?! MAXDIM !*/
#define EaxxaMaxDim 10  /* nur technisch */
//#define UmatrixMaxDim MAXMLEDIM /* nur technisch */
#define ShiftMaxDim 10  /* nur technisch ?! */
#define ParsWMMaxVDim 10 


void vstop();
void vassert();
void vassert(bool b);


double XkCXtl(double *X, double *C, int nrow, int dim, int k, int l);
void XCXt(double *X, double *C, double *V, int nrow, int dim);
void AtA(double *a, int nrow, int ncol, double *A);
void xA(double *x, double*A, int nrow, int ncol, double *y);
void Ax(double *A, double*x, int nrow, int ncol, double *y);
void AxResType(double *A, res_type *x, int nrow, int ncol, double *y);

SEXP Param(void* p, int nrow, int ncol, SEXPTYPE type, bool drop, long* mem);
int is_positive_definite(double *C, int dim);

// typedef long long POINTER;
typedef long int POINTER;

extern int PrInL;				

void SetTemporarilyGlobal(globalparam *gp);
void ReSetGlobal();
				
int change_coordinate_system(isotropy_type callingprev, isotropy_type isoprev,
			     int *nr, isotropy_type  *newisoprev,
			     int *newtsdim, int *newxdim);
int SetGatter(domain_type domprev, domain_type statnext, 
	      isotropy_type isoprev, isotropy_type isonext, 
	      int *nr, int *delflag);
				
int cpy_add_gatter(cov_model *cov, cov_model **sm);

double Real(SEXP p, char *name, int idx);
int Integer(SEXP p, char *name, int idx);
bool Logical(SEXP p, char *name, int idx);

int binomialcoeff(int n, int k);
void poly_basis(cov_model *cov, storage *s);
double evalpoly(double *x, int *powmatrix, double *polycoeff, int basislen,
	      int dim);

void MultiDimRange(cov_model *cov, double *natscale);

#define DEL_COV -100 // #->tsdim wird auf del_cov gesetzt, wenn zu loeschen  
typedef struct cov_model {
  // user given information
  int nr, /* number of the covariance function in CovList; 
	     obtained by function GetCovFunction 
	  */
    gatternr,/* trafo function for the coordinates */
    secondarygatternr, /* used if first one is a cooord system trafo */
    zaehler; /* for debugging only */
  param_type px;       // 24 b
  int nrow[MAXPARAM], // 24 bytes
    ncol[MAXPARAM];   // 24 bytes
  double *q;
  int qlen;
  int nsub; /* current number of submodels */
  cov_model *sub[MAXSUB], *kappasub[MAXPARAM], *calling; 
  char **ownkappanames;
  user_given_type user_given;
 
  /////////////////////////////////////////////////////////
  // VARIABLES PASSED DOWNWARDS
  /////////////////////////////////////////////////////////
  Types typus;  /* CovList->Type may change in case of operators */
  int role, /* current role of the model */
     tsdim,  /* current timespacedim  */
    xdimprev,  /* current xdim *including time* (when entering # or >),
		  which may differ from tsdim for isotropic or
		  space-isotropic models   */
    xdimgatter, /* only differs from xdimprev if > is called and 
		 only used by #*/
    xdimown,  /* current xdim *including time* (when leaving #),
		 which may differ from tsdim for isotropic or
		 space-isotropic models   */
  // vdim, obsolete because of next line /* multivariate dimension passed upwards */
    vdim2[2]; /* for non-quadratic returns */


  /////////////////////////////////////////////////////////
  // VARIABLES PASSED UPWARDS
  /////////////////////////////////////////////////////////
  // forward analysis of user's information
  domain_type domprev, domown;
                    /* formal domain, s. method->type f. Inhalt
			     IN parameter is usually set fix; 
		    */
                         /*  OUT must be managed by the check function,
			     since OUT may vary among the submodels, by 
			     calling "#" through gatternr
			  */
  isotropy_type isoprev, isoown;     /* formal isotropic parameter */
  double logspeed; /* 
		      logspeed = lim_{h->infty} \gamma(h)/log(h) 
		      in case of isotropic model and RF_NAN otherwise
		   */

  // backward analysis of user's information
  //  double  total_sum; // only for RandomType
  int
  //    total_n, // only for RandomType
    delflag,//currently unused: if != 0 dann sollte gatternr = nr gesetzt werden
    maxdim, /* maxdim of model combined with information of the submodels
	        */
                 
    full_derivs, rese_derivs; /* rechtsseitige bzw. vollstaendige Anzahl 
				 Ableitungen number of derivatives of model 
				 combined with information of the submodels */

  ext_bool monotone, /* for simple model: normal mix model iff maxdim = INFDIM
		     for hypermodel we need nonetheless this parameter
		     parameter set by getinfo(), not within info of the cov 
		     fcts comining information of submodels  
		     */
    finiterange; /* also information obtained by model and submodels 
		  */
  bool deterministic,
    matrix_indep_of_x,  
    diag,       /*  */
    semiseparatelast, /*  */
    separatelast,     /*  */
    tbm2num,  /* "time" component separately?  */
    hess;  /* can a hessian matrix be provided? */
    

  pref_type pref; /* 0 : not possible; 
		     5 : best possible
		     (including subs)!
		   */
  //   user; /* 0 : not possible; 
  //	     5 : best possible
  //	     (including subs)!
  //		  */
  double *MLE;

  // Information used and created within simuations
  Methods method; /* the current method (out of SimulationsType) which 
			    is tried at the moment or which has been 
			    successfully initialized */
  mpp_properties mpp;

  simu_type simu;
  location_type *prevloc, *ownloc;
  cov_model *key; // this one should be deleted
  // double *caniso, cscale, cvar; 
  //  int *cproj;
  //  int xdimout;
  // matrix_type type;   // von oben
  
  bool initialised, origrf, checked;
  ext_bool loggiven, fieldreturn;
  res_type
    *rf; // for storing random shapes; this is a pointer pointing
  //               to *rf in storage (logically rf only belongs to 
  //               a single atom. But programmed as such, this would take
  //               an enormous amount of memory in case of M3 processes)

  double taylor[MAXTAYLOR][TaylorPow + 1], 
    tail[MAXTAYLOR][TaylorExpPow + 1]; 
  int taylorN, // number of summands in the taylor expansion -- whatever the
  // exponents of the terms are
    tailN; // dito


  CE_storage *SCE;
  localCE_storage *SlocalCE;
  ce_approx_storage *SapproxCE;
  direct_storage *Sdirect;
  hyper_storage *Shyper;
  mixed_storage *Smixed;
  nugget_storage *Snugget;
  plus_storage *Splus;
  sequential_storage *Sseq;
  TBM_storage *Stbm;
  trend_storage *Strend;
  BR_storage *SBR;
  get_storage *Sget;
  pgs_storage *Spgs;
  set_storage *Sset;
  polygon_storage *Spolygon;
  rect_storage *Srect;
  dollar_storage *Sdollar;
  gatter_storage *S2;
  biwm_storage *Sbiwm;
  inv_storage *Sinv;
  //select_storage *Sselect;
  storage *stor; // only once for the whole model tree
} cov_model;


void LOC_DELETE(location_type **Loc);
void LOC_NULL(location_type *loc);
void COV_DELETE(cov_model **cov);
void COV_NULL(cov_model *cov);
void COV_DELETE_WITHOUTSUB(cov_model **Cov);
void COV_DELETE_WITHOUT_LOC(cov_model **Cov);
void CE_NULL(CE_storage* x);
void CE_DELETE(CE_storage **S);
void LOCAL_NULL(localCE_storage* x);
void LOCAL_DELETE(localCE_storage**S);
void CE_APPROX_NULL(ce_approx_storage* x);
void CE_APPROX_DELETE(ce_approx_storage **S);
void DIRECT_NULL(direct_storage  *x);
void DIRECT_DELETE(direct_storage  ** S);
void HYPER_NULL(hyper_storage* x);
void HYPER_DELETE(hyper_storage  **S);
void NUGGET_NULL(nugget_storage *x);
void NUGGET_DELETE(nugget_storage ** S);
void PLUS_NULL(plus_storage *x);
void PLUS_DELETE(plus_storage ** S);
void SEQU_NULL(sequential_storage *x);
void SEQU_DELETE(sequential_storage **S);
void SPECTRAL_NULL(sequential_storage *x);
void SPECTRAL_DELETE(sequential_storage **S);
void TBM_DELETE(TBM_storage **S); 
void TBM_NULL(TBM_storage* x);
void BR_DELETE(BR_storage **S); 
void BR_NULL(BR_storage* x);
void GET_STORAGE_NULL(get_storage *S);
void GET_STORAGE_DELETE(get_storage **S);

void PGS_DELETE(pgs_storage **S); 
void PGS_NULL(pgs_storage* x);

void SET_DELETE(set_storage **S); 
void SET_NULL(set_storage* x);

void POLYGON_DELETE(polygon_storage **S); 
void POLYGON_NULL(polygon_storage* x);

void RECT_DELETE(rect_storage **S); 
void RECT_NULL(rect_storage* x);

void DOLLAR_DELETE(dollar_storage **S); 
void DOLLAR_NULL(dollar_storage* x);

void GATTER_DELETE(gatter_storage **S); 
void GATTER_NULL(gatter_storage* x);

void BIWM_DELETE(biwm_storage **S); 
void BIWM_NULL(biwm_storage* x);

void INV_DELETE(inv_storage **S);
void INV_NULL(inv_storage* x);


void TREND_DELETE(trend_storage ** S);
void TREND_NULL(trend_storage* x);
void MIXED_DELETE(mixed_storage **S);
void MIXED_NULL(mixed_storage* x);
void STORAGE_NULL(storage *x);
void STORAGE_DELETE(storage **S);

//void SELECT_NULL(select_storage *S);
//void SELECT_DELETE(select_storage **S);


///////////////////////////////////////////
#define Gettotalpoints(cov)(cov->ownloc != NULL ? cov->ownloc->totalpoints : \
			    cov->prevloc != NULL ? cov->prevloc->totalpoints :0)
#define Gettimespacedim(cov)(cov->ownloc != NULL ? cov->ownloc->timespacedim: \
			     cov->prevloc!=NULL ? cov->prevloc->timespacedim :0)
#define GetTime(cov)(cov->ownloc != NULL ? cov->ownloc->Time :	\
                     cov->prevloc != NULL ? cov->prevloc->Time : false)
#define Getgrid(cov) (cov->ownloc != NULL ? cov->ownloc->grid:		\
		      cov->prevloc != NULL ? cov->prevloc->grid : false)
#define Loc(cov) (cov->ownloc != NULL ? cov->ownloc : cov->prevloc)

void crash(cov_model *cov);
void crash();

sortsofparam ignoreall_paramtype(int k, int row, int col);
int checktbm_basics(cov_model *cov, bool tbmop);
void SetLoc2NewLoc(cov_model *cov, location_type *loc);

extern int gaussmethod[Forbidden + 1];
void FFT_destruct(FFT_storage *FFT);

int FieldReturn(cov_model *cov);
void DoPPfields(int *n);
int initmpp(double *x, double *T, int dim, int lx, bool grid, bool distances, 
	    bool Time, cov_model  *key, globalparam *gp,
	    int expected_number_simu, cov_model **COV);

void ErrLogCov(double *x, cov_model *cov, double *v, double *sign);
void ErrLogCovNonstat(double *x, double *y, cov_model *cov, double *v, 
		      double *sign);
int struct_failed(cov_model *cov, cov_model **atom);
int init_failed(cov_model *cov, storage *s);
int init_statiso(cov_model *cov, storage *s);
void do_failed(cov_model *cov, storage *s);
void do_statiso(cov_model *cov, storage VARIABLE_IS_NOT_USED *s);
int checkOK(cov_model *cov);
int checkplusmal(cov_model *cov);
void ErrInverse(double *v, cov_model *cov, double *x);
void InverseIsotropic(double *U, cov_model *cov, double *inverse);


void AtA(double *a, int nrow, int ncol, double *A) ;
void xA(double *x, double*A, int nrow, int ncol, double *y) ;
void AxResType(double *A, res_type *x, int nrow, int ncol, double *y) ;
void Ax(double *A, double*x, int nrow, int ncol, double *y) ;
double XkCXtl(double *X, double *C, int nrow, int dim, int k, int l) ;
double xUy(double *x, double *U, double *y, int dim) ;
double xUxz(double *x, double *U, int dim, double *z) ;
double x_UxPz(double *x, double *U, double *z, int dim) ;
double xUx(double *x, double *U, int dim) ;
double detU(double *C, int dim); // strictly pos def matrices
void det_UpperInv(double *C, double *det, int dim) ;
void matmult(double *A, double *B, double *C, int l, int m, int n) ;
void matmulttransposed(double *A, double *B, double *C, int m, int l, int n) ;
int invertMatrix(double *M, int size) ;
//void solveMatrix(int method, double *M, double *v, double *res, int size) ;
double getMinimalAbsEigenValue(double *Aniso, int dim);
double getDet(double *Aniso, int dim); // arbitrary matrix

void iexplDollar(cov_model *cov, bool MLEnatsc_only);


void BRtrend_destruct(double **BRtrend, int trendlen);
int StructBR(cov_model *cov, storage *s, cov_model **atom);

sortsofparam uncritical_paramtype(int k, int row, int col);

void indextrafo(int onedimindex, int *length, int dim, int *multidimindex);

int check_within_range(cov_model *cov, bool NAOK);
extern bool RELAX_UNKNOWN_RFOPTION;


int INIT_RANDOM_intern(cov_model *M, int moments, storage *s, double *p);
int INIT_intern(cov_model *M, int moments, storage *s);
int REINIT_intern(cov_model *M, int moments, storage *s);
void set_initialised_false(cov_model *cov, bool init_deterministic);
  
int alloc_mpp_M(cov_model *cov, int moments);
void free_mpp_M(cov_model *cov);

int alloc_pgs(cov_model *cov, int dim);
int alloc_pgs(cov_model *cov);
int alloc_cov(cov_model *cov, int dim, int rows, int cols);


#define DO(Cov, S) {assert(!TypeConsistency(RandomType, Cov)); \
    assert(Cov->initialised);				       \
    ASSERT_GATTER(Cov);					       \
    PL--;						       \
    CovList[(Cov)->gatternr].Do(Cov, S);		       \
    PL++;						       \
  }
   
#define DORANDOM(Cov, S) {\
    assert(TypeConsistency(RandomType, Cov));			      \
    ASSERT_GATTER(Cov);						      \
    PL--;							      \
    CovList[(Cov)->gatternr].DoRandom(Cov, S);			      \
    PL++;							      \
  }
#define COV(X, Cov, V) {  ASSERT_GATTER(Cov); CovList[(Cov)->gatternr].cov(X, Cov, V);}
#define LOGCOV(X, Cov, V, S) {ASSERT_GATTER(Cov);CovList[(Cov)->gatternr].log(X, Cov, V, S);}
#define SHAPE COV
#define FCTN COV
#define LOGSHAPE LOGCOV
#define VTLG_D(X, Cov, V) { ASSERT_CHECKED(Cov); CovList[(Cov)->nr].D(X, Cov, V);} // kein gatter notw.
#define VTLG_DLOG(X, Cov, V) { ASSERT_CHECKED(Cov); CovList[(Cov)->nr].logD(X, Cov, V);}
#define VTLG_P(X, Cov, V) { ASSERT_CHECKED(Cov); CovList[(Cov)->nr].cov(X, Cov, V);}
#define VTLG_P2SIDED(X, Y, Cov, V) { ASSERT_CHECKED(Cov); CovList[(Cov)->nr].nonstat_cov(X, Y, Cov, V);} /* nicht gatter, da X=NULL sein kann !*/ 
#define VTLG_Q(V, Cov, X) { ASSERT_CHECKED(Cov); CovList[(Cov)->nr].inverse(V, Cov, X);}
#define VTLG_R(X, Cov, V) { ASSERT_CHECKED(Cov); CovList[(Cov)->nr].random(X, Cov, V);} /* dito */
#define VTLG_R2SIDED(X, Y, Cov, V) { ASSERT_CHECKED(Cov); CovList[(Cov)->nr].nonstat_random(X, Y, Cov, V);}
#define NONSTATINVERSE_D(V, Cov, X, Y)		       \
  { ASSERT_CHECKED(Cov);  CovList[(Cov)->nr].nonstat_inverse_D(V, Cov, X, Y);}

#define CHECK_R(cov, vdim)				     \
  CHECK_VDIM(cov, vdim, vdim, RandomType, KERNEL, CARTESIAN_COORD, \
	  vdim, 1, ROLE_DISTR)			       

//#define DENSITY(X, Dens, V) CovList[(Cov)->gatternr].approx_dens(X, Dens, V)
#define NONSTATCOV(X, Y, Cov, V) {ASSERT_GATTER(Cov);CovList[(Cov)->gatternr].nonstat_cov(X, Y, Cov,V);}
#define LOGNONSTATCOV(X, Y, Cov, V, S) {ASSERT_GATTER(Cov);CovList[(Cov)->gatternr].nonstatlog(X, Y, Cov,V,S);}
#define Abl1(X, Cov, V) {ASSERT_GATTER(Cov);CovList[(Cov)->gatternr].D(X, Cov, V);}
#define Abl2(X, Cov, V) {ASSERT_GATTER(Cov);CovList[(Cov)->gatternr].D2(X, Cov, V);}
#define Abl3(X, Cov, V) {ASSERT_GATTER(Cov);CovList[(Cov)->nr].D3(X, Cov, V);} /* OK ? */
#define Abl4(X, Cov, V) {ASSERT_GATTER(Cov);CovList[(Cov)->nr].D4(X, Cov, V);} /* OK ? */
#define SPECTRAL(Cov, S, E) {ASSERT_GATTER(Cov); CovList[(Cov)->nr].spectral(Cov, S, E);} /* not gatter */
#define TBM2CALL(X, Cov, V) {ASSERT_GATTER(Cov); CovList[(Cov)->nr].tbm2(X, Cov, V);}
#define INVERSE(V, Cov, X) {ASSERT_GATTER(Cov);CovList[(Cov)->gatternr].inverse(V, Cov, X);}
#define NONSTATINVERSE(V, Cov, X, Y) {ASSERT_GATTER(Cov);\
				      CovList[(Cov)->gatternr].nonstat_inverse(V, Cov, X, Y);}
#define NONSTATLOGINVERSE(V, Cov, X, Y) {ASSERT_GATTER(Cov);\
				      CovList[(Cov)->gatternr].nonstat_loginverse(V, Cov, X, Y);}
#define HESSE(X, Cov, V) {ASSERT_GATTER(Cov);CovList[(Cov)->nr].hess(X, Cov, V);}
#define NABLA(X, Cov, V) {ASSERT_GATTER(Cov);CovList[(Cov)->nr].nabla(X, Cov, V);}


#define ROLE_ASSERT(Role) \
  if (!(cov->role == ROLE_BASE || cov->role == Role)){/* NICHT! : '(Role)' */ \
    assert(({PMI(cov) ; true;}));					\
    SERR2("Role '%s' not recognised by '%s'.",				\
	  ROLENAMES[cov->role], NICK(cov));				\
}

#define ILLEGAL_ROLE							\
  SERR4("cannot initiate '%s' by role '%s' [debug info: '%s' at line %d]", \
	NICK(cov), ROLENAMES[cov->role], __FILE__, __LINE__)

#define ILLEGAL_ROLE_STRUCT \
  SERR2("cannot restructure '%s' by role '%s'", NICK(cov), ROLENAMES[cov->role])

#define ROLE_ASSERT_GAUSS if (cov->role != ROLE_GAUSS) ILLEGAL_ROLE 
   
#define ASSERT_GAUSS_METHOD(METHOD) DEBUGINFO;				\
  if(cov->role != GAUSS || cov->method != Average)			\
    SERR2("Gaussian field for '%s' only possible with '%s' as specific method",\
	 NICK(cov),							\
	 CovList[METHOD == RandomCoin ? RANDOMCOIN_USER	:		\
		 gaussmethod[METHOD] -					\
		 CovList[gaussmethod[METHOD]].internal].nick)

#define ASSERT_NEWMODEL_NOT_NULL\
  if (newmodel == NULL) {\
    SERR1("unexpected call of struct_%s", NICK(cov));\
  }

#define ASSERT_NEWMODEL_NULL\
  if (newmodel != NULL) {\
    SERR1("Unexpected call of struct_%s", NICK(cov));\
  }

#define ASSERT_ONE_SUBMODEL(COV)					\
  DEBUGINFO;								\
  if (!(((COV)->sub[0] == NULL) xor ((COV)->sub[1] == NULL)))		\
    SERR2("either '%s' or '%s' must be given", CovList[(COV)->nr].subnames[0],\
      CovList[(COV)->nr].subnames[1])

#define ASSERT_ROLE_DEFINED(SUB) DEBUGINFO;			\
  if (role == ROLE_UNDEFINED)					\
    SERR1("'%s' not allowed as shape function.", NICK(SUB));	\

void PSTOR(cov_model *cov, storage *x);
int structOK(cov_model *cov, cov_model **newmodel);
int initOK(cov_model *cov, storage *s);
void doOK(cov_model *cov, storage *s);
void do_random_ok(cov_model *cov, double *v);
void do_random_failed(cov_model *cov, double *v);


#define APMI(cov) {\
    PRINTF("\n'%s', line %d", __FILE__, __LINE__); \
    pmi(cov); assert(false);}
#define AAPMI(cov, txt) { pmi(cov, txt); assert(false); }

#define PMI \
  PRINTF("\n(PMI '%s', line %d)", __FILE__, __LINE__);	\
  pmi


#define PMI \
  PRINTF("\n(PMI '%s', line %d)", __FILE__, __LINE__);	\
  pmi


#define PLE \
  PRINTF("\n(PLE '%s', line %d)", __FILE__, __LINE__);	\
  ple
void ple(cov_model *cov);
void ple(char *name);



void memory_copy(void *dest, void *src, int bytes);
int invertMatrix(double *M, int vdimnf);

int role_of_process(int nr);

int addShapeFct(cov_model **Cov);


#define BALL_RADIUS 1.0 // nicht mehr aendern!!
double SurfaceSphere(int d, double r);
double VolumeBall(int d, double r);

double intpow(double x, int p);

bool isDummyInit(initfct Init);

bool isMtimesep(matrix_type type);
bool isMproj(matrix_type type);
bool isMdiag(matrix_type type);
bool isMiso(matrix_type type);

void ScaleDollarToLoc(cov_model *to, cov_model *from, int depth);
bool ScaleOnly(cov_model *cov);
bool isDollar(cov_model *cov);
bool isDollarProc(cov_model *cov);
bool isAnyDollar(cov_model *cov);
bool isNatsc(cov_model *cov);

bool isTcf(cov_model *cov);
bool isPosDef(cov_model *cov);
bool isNegDef(cov_model *cov);
bool isProcess(cov_model *cov);
bool isGaussMethod(cov_model *cov);
bool isBRuserProcess(cov_model *cov);
bool isPointShape(cov_model *cov);
bool isRandom(cov_model *cov);
bool isShape(cov_model *cov);
bool isTrend(cov_model *cov);
bool isInterface(cov_model *cov);
bool isOtherType(cov_model *cov);

bool isTcf(Types type);
bool isPosDef(Types type);
bool isNegDef(Types type);
bool isDefinite(Types type);
bool isProcess(Types type);
bool isGaussMethod(Types type);
bool isBRuserProcess(Types type);
bool isPointShape(Types type);
bool isRandom(Types type);
bool isShape(Types type);
bool isTrend(Types type);
bool isInterface(Types type);
bool isUndefinedType(Types type);
bool isOtherType(Types type);

bool isMonotone(cov_model *cov);
bool isMonotone(int montone);
bool isNormalMixture(int montone);
bool isNormalMixture(cov_model *cov);
bool isBernstein(cov_model *cov);
bool isBernstein(int monotone);
bool isGneiting(cov_model *cov);
bool isGneiting(int monotone);

bool isGaussProcess(cov_model *cov);
bool isGaussMethod(cov_model *cov);
bool isBernoulliProcess(cov_model *cov);
bool isGaussBasedProcess(cov_model *cov);
bool isBrownResnickProcess(cov_model *cov);
bool isMaxStable(cov_model *cov);

bool isCov(cov_model *cov);

bool isIsotropic(isotropy_type iso);
bool isSpaceIsotropic(isotropy_type iso);
bool isZeroSpaceIsotropic(isotropy_type iso);
bool isVectorIsotropic(isotropy_type iso);
bool isSymmetric(isotropy_type iso);
bool isNotRotatInv(isotropy_type iso);
bool isCartesian(isotropy_type iso);
bool isSpherical(isotropy_type iso);
bool isCylinder(isotropy_type iso);
bool isEarth(isotropy_type iso);
bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2);
bool TrafoOK(cov_model *cov, bool all);

bool hasNoRole(cov_model *cov);
bool hasMaxStableRole(cov_model *cov);
bool hasExactMaxStableRole(cov_model *cov);
bool hasPoissonRole(cov_model *cov);
bool hasAnyShapeRole(cov_model *cov);
bool hasDistrRole(cov_model *cov);

bool TypeConsistency(Types requiredtype, Types deliveredtype);
bool TypeConsistency(Types requiredtype, cov_model *cov);
int CheckPD2ND(cov_model *cov, int tsdim, int tsxdim, isotropy_type isoprev,
	       int vdim, int role);

int role_of_process(int nr);
void Taylor(double c, double pow);
void TailTaylor( double t, double tpow, double texp, double texppow);
void Taylor(double c, double pow, double c1, double pow1);
void Taylor(double c, double pow, double c1, double pow1, 
	    double c2, double pow2);
int searchFirstGreater(double *v, int len, double z);
int CeilIndex(double x, double *cum, int size);
double searchInverse(covfct fct, cov_model *cov, 
		     double start, double value, double releps);
double searchInverse(covfct fct, cov_model *cov, 
		     double start, double min, double value, double releps);
// double gamma(double x);
double incomplete_gamma(double start, double end, double s);
void ErrInverseNonstat(double *v, cov_model *cov, double *x, double *y);
void StandardInverseNonstat(double *v, cov_model *cov,
			    double *left, double *right);
double GetTotal(cov_model* cov);
int addressbits(void *addr);


#define MIN(A,B) ((A) < (B) ? (A) : (B));
#define MAX(A,B) ((A) > (B) ? (A) : (B));


#define EXTRA_STORAGE						\
  if (cov->S2 != NULL && cov->S2->z != NULL)			\
    GATTER_DELETE(&(cov->S2));					\
  if (cov->S2 == NULL) {					\
    cov->S2 = (gatter_storage*) MALLOC(sizeof(gatter_storage));	\
    GATTER_NULL(cov->S2);					\
  }								\
  assert(cov->S2->z == NULL);


#define ALLOC_EXTRA1(Z, SIZE)						\
  double *Z = cov->S2->z;						\
  if (Z == NULL) Z = cov->S2->z = (double*) MALLOC(sizeof(double)* SIZE)


// hier keine klammern drum rum
#define ALLOC_EXTRA2(Z, SIZE)						\
  double *Z = cov->S2->zsys;						\
  if (Z == NULL) Z = cov->S2->zsys = (double*) MALLOC(sizeof(double)* (SIZE))


// Polygon Functions
double scProd(const double *x, const double *y);
int compareAngles(const void * a, const void * b);
void rTriangle(double *phi);
int rPoissonPolygon(struct polygon *P, double lambda, bool do_centering);
int rPoissonPolygon2(polygon_storage *S, double lambda, bool do_centering);
void freePolygon(struct polygon *P);
bool isInside(struct polygon *P, double *x);
double getArea(struct polygon *P);


int addPointShape(cov_model **Key, cov_model *shape, cov_model *pts, 
		  cov_model *local_pts,int dim, int vdim);

int addPointShape(cov_model **Key, cov_model *shape, cov_model *pts, 
		  int dim, int vdim);

int PointShapeLocations(cov_model *key, cov_model *shape);




//extern void F77_NAME(zpotf)(int* info);
extern "C" void F77_CALL(zpotf2)(char *name, int *row, complex *U, int *xxx,
				 int *Err);
				

typedef int (*inituser)(double *x, double *T, int dim, int lx, bool grid, 
			bool distances, bool Time, globalparam *gp,
			int expected_number_simu, cov_model **COV);
typedef void (*douser)(int *n, double *res);


#endif
