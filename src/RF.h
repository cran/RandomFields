
#ifndef RFsimu_H
#define RFsimu_H 1

#define showfree !true 
#define DOPRINT true
//// 1


#include "error.h"
#include <string.h>
#include "auxiliary.h"
#include "RandomFields.h"
#include "AutoRandomFields.h"

//   intptr_t and uintptr_t fuer: umwandlung pointer in int und umgekehrt

#define ASSERT_GATTERONLY(Cov) assert(TrafoOK(Cov, false))
#define ASSERT_GATTER(Cov) assert(TrafoOK(Cov, true))
#define ASSERT_CHECKED(Cov) assert(Cov->checked)

#define DOPRINTF  if (DOPRINT) Rprintf
#define KPRINT leer(PrInL);Rprintf
#define LPRINT {cov_model *lprint_z=cov; int lprint_i=0; while (lprint_z->calling != NULL && lprint_i<10) {lprint_z=lprint_z->calling; if (DOPRINT) {Rprintf(DOT); Rprintf(" ");} lprint_i++;} if (lprint_i==100) {Rprintf("LPRINT i=%d\n", lprint_i);PMI(cov); assert(false);}} if (DOPRINT) Rprintf


#ifndef RANDOMFIELDS_DEBUGGING

#define ERRLINE 

#define NICK(COV) (isDollar(COV) ? CovList[(COV)->sub[0]->nr].nick : CovList[(COV)->nr].nick)
 


#define MEMCOPY(A,B,C) memcpy(A,B,C)
//#define MEMCOPY(A,B,C) memory_copy(A, B, C)

#define MALLOC malloc
#define CALLOC calloc
#define FREE(X) if ((X) != NULL) {free(X); (X)=NULL;}
#define UNCONDFREE(X) {free(X); (X)=NULL;}

#define CHECK(C,T,X,type,D,I,V,R) check2X(C,T,X,type,D,I,V,R)
#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,R) check2X(C,T,X,type,D,I,V0,V1,R)
#define CHECKPD2ND(N,D1,D2,I,V,R) CheckPD2ND(N,D1,D2,I,V,R)
#define INIT(M, Moments, S) INIT_intern(M, Moments, S)
#define REINIT(M, Moments, S) REINIT_intern(M, Moments, S)
#define INIT_RANDOM(M, Moments, S, P) INIT_RANDOM_intern(M, Moments, S, P)
#define STRUCT(Cov, NM) CovList[(Cov)->gatternr].Struct(Cov, NM)


#define PARAM(P, IDX) ((double *) (P)->px[IDX])
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
 
#define QcovALLOC(cov, nr) {\
  cov->qlen = (int) (nr);						\
  if ((cov->q = (double*) CALLOC((int) (nr), sizeof(double))) == NULL)	\
    ERR("memory allocation error for local memory");		\
}

#define DEBUGINFOERR
#define DEBUGINFO 

#else
#define ERRLINE ERRLINES

#define NICK(COV) (CovList[(COV)->nr].nick)


#define MEMCOPY(A,B,C) ({ assert((A)!=NULL && (B)!=NULL); memcpy(A,B,C); })
//#define MEMCOPY(A,B,C) memory_copy(A, B, C)

//
#define MALLOC(X) ({DOPRINTF("(MALLOC %s, line %d)\n", __FILE__, __LINE__);assert(X>0 && X<=668467200);malloc(X);})
//
#define CALLOC(X, Y) ({DOPRINTF("(CALLOC %s, line %d)\n",__FILE__, __LINE__);assert((X)>0 && X<1e8 && (Y)>0 && (Y)<=64); calloc(X,Y);})
//#define MALLOC malloc
//#define CALLOC calloc
#define FREE(X) { if (showfree) DOPRINTF("(free in %s, line %d)\n", __FILE__, __LINE__); if ((X) != NULL) free(X); (X)=NULL; }
#define UNCONDFREE(X) { if (showfree) DOPRINTF("(free in %s, line %d)\n", __FILE__, __LINE__); free(X); (X)=NULL;}
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
  }DOPRINTF("(");DOPRINTF("%s", Z);DOPRINTF(" %s, line %d : %s)\n", __FILE__, __LINE__, NAME(cov))

#define CHECKSIGN "_"
#define INITSIGN "_"
#define STRUCTSIGN "_"
#define CHECK(C,T,X,type,D,I,V,R) ({LLPRINT(CHECKSIGN, C, "CHECK"); assert(type!=RandomType); XX(C); int _x = check2X(C,T,X,type,D,I,V,R); YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED"); errorMSG(_x, MSG, LENERRMSG); PRINTF("%s%s\n", ERROR_LOC, MSG);} _x;})
#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,R) ({LLPRINT(CHECKSIGN, C, "CHECKVDIM"); XX(C); int _x = check2X(C,T,X,type,D,I,V0,V1,R);YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED");} _x;})
#define CHECKPD2ND(C,D1,D2,I,V,R) ({LLPRINT(CHECKSIGN, C, "CHECKPD2ND"); XX(C); int _x = CheckPD2ND(C,D1,D2,I,V,R);YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED");} _x;})
#define INIT(C, Moments, S) ({LLPRINT(INITSIGN, C, "INIT");  XX(C); int _x = INIT_intern(C, Moments, S);YY(C); if (_x==NOERROR){LLPRINT(STRUCTSIGN, C, "INIT DONE");}else{LLPRINT(STRUCTSIGN, C, "INIT FAILED");}_x;})
#define REINIT(C, Moments, S) ({LLPRINT(INITSIGN, C, "INIT");  XX(C); int _x = REINIT_intern(C, Moments, S); YY(C); _x;})

#define INIT_RANDOM(C, Moments, S, P) ({LLPRINT(INITSIGN, C, "INITRANDOM");  XX(C); int _x = INIT_RANDOM_intern(C, Moments, S, P);YY(C); if (_x==NOERROR){LLPRINT(STRUCTSIGN, C, "INIT DONE");}else{LLPRINT(STRUCTSIGN, C, "INIT FAILED");}_x;})
#define STRUCT(C, NM)  ({LLPRINT(STRUCTSIGN, C, "STRUCT"); ASSERT_GATTER(C);  XX(C); int _x = CovList[(C)->gatternr].Struct(C, NM);YY(C); if (_x==NOERROR){LLPRINT(STRUCTSIGN, C, "STRUCT DONE");}else{LLPRINT(STRUCTSIGN, C, "STRUCT FAILED");}_x;})


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

#define QcovALLOC(cov, nr) {\
  assert(cov->q == NULL);   \
  cov->qlen = (int) (nr);						\
  if ((cov->q = (double*) CALLOC((int) (nr), sizeof(double))) == NULL)	\
    ERR("memory allocation error for local memory");		\
}

#define DEBUGINFOERR {							\
    char dummy[MAXERRORSTRING]; strcpy(dummy, ERRORSTRING);		\
    sprintf(ERRORSTRING, "%s (%s, line %d)\n", dummy, __FILE__, __LINE__); \
  }
#define DEBUGINFO DOPRINTF("(currently at  %s, line %d)\n", __FILE__, __LINE__)

#endif



#define P(IDX) PARAM(cov, IDX) 
#define PINT(IDX) PARAMINT(cov, IDX)
#define P0(IDX) PARAM0(cov, IDX) 
#define P0INT(IDX) PARAM0INT(cov, IDX)
#define PENV(IDX) PARAMENV(cov, IDX) 
#define PSEXP PENV
#define PLIST(IDX) PARAMLIST(cov, IDX)

#define PARAMFREE(P, IDX) if ((P)->px[IDX] != NULL) { UNCONDFREE((P)->px[IDX]); (P)->ncol[IDX] = (P)->nrow[IDX] = 0;}


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
#define PINTALLOC(IDX, ROW, COL) PARAMINTALLOC(cov, IDX, ROW, COL)

#define PFREE(IDX) PARAMFREE(cov, IDX)
#define PARAMisNULL(P, IDX) ((P)->px[IDX] == NULL)
#define PisNULL(IDX) PARAMisNULL(cov, IDX)
#define PARAMtoNULL(P, IDX) (P)->px[IDX] = NULL
#define PtoNULL(IDX) PARAMtoNULL(cov, IDX)

#define QALLOC(nr) QcovALLOC(cov, nr)
#define QFREE					\
  if (cov->q != NULL) {				\
    UNCONDFREE(cov->q);				\
    cov->qlen = 0;				\
  }


#define FIRST_INIT(init)				\
  int xerr;						\
  gen_storage s;					\
  gen_NULL(&s);						\
  s.check = true;					\
  if ((xerr=init(cov, &s)) != NOERROR) return xerr;


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
#define ROLE_POISSON_GAUSS 8 // umsetzung von gauss nach random coin model
#define ROLE_HYPER 9
#define ROLE_BERNOULLI 10
#define ROLE_DISTR 11
#define ROLE_FAILED 12
#define ROLE_UNDEFINED 13
#define ROLE_LAST ROLE_UNDEFINED



///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////
// type of covariance functions that need distinguished treatment
// Reihenfolge beachten! Hoehere Nummern sind meist in kleinere umwandelbar




// vdim is "strictly" passed from bottom to top
#define SCALAR 1 

// way of implementing simulation methods:
#define NOT_IMPLEMENTED 0 /* do not change this value except to false */
#define IMPLEMENTED 1     /* do not change this value except to true */
#define NUM_APPROX 2 
#define GIVEN_METH_IGNORED 3
#define HYPERIMPLEMENTED 4



///////////////////////////////////////////////////////////////////////
// BASIC DIMENSIONS AND VARIABLES
///////////////////////////////////////////////////////////////////////
#define PARAMMAXCHAR MAXCHAR
#define MAXLONGCHAR 40 // for RFoptions that are string
#define MAXUNITSCHAR 10

#define MAXNRCOVFCTS 300
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


#define NATSCALE_EXACT 1 /* or approx or mle */
#define NATSCALE_ORNUMERIC 2
//#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check fitvario when changing; best not to change !!! */


#define LISTOF 100  /* may not interfere with define of *SXP in Rinternal.h*/


// #define MAXINTERNALPARAM 2
extern double ZERO[MAXSIMUDIM], ONE;
extern int True;
extern int False;
extern char *FREEVARIABLE;

typedef char name_type[][MAXCHAR];
typedef char longname_type[][MAXLONGCHAR];

typedef char NAname_type[MAX_NA][255];



///////////////////////////////////////////////////////////////////////
// PARAMETER AND THEIR LOCATIONS FOR VARIOUS MODELS
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// Dollar
#define DOLLAR_SUB 0
#define DVAR 0
#define DSCALE 1
#define DANISO 2 /* internal */
#define DAUSER 3 /* user defined */
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
#define ANGLE_LATANGLE 1
#define ANGLE_RATIO 2
#define ANGLE_DIAG 3


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

#define ELEMENT "element"
#define RFOPTIONS "RFoptions"

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
// Power, shapepow
#define POW_ALPHA 0

///////////////////////////////////////////////////////////////////////
// trafo
#define TRAFO_ISO 0


///////////////////////////////////////////////////////////////////////
// Scatter
#define SCATTER_STEP 0
#define SCATTER_MAX 1

///////////////////////////////////////////////////////////////////////
// fix
#define FIX_ELMNT 0
#define FIX_M 1
#define FIX_VDIM 2
#define FIX_VAR 3



///////////////////////////////////////////////////////////////////////
// coins
#define COIN_COV 0
#define COIN_SHAPE 1

///////////////////////////////////////////////////////////////////////
// hyper
#define HYPER_UNIFORM 0   
#define HYPER_FRECHET 1
#define HYPER_BERNOULLI 2

///////////////////////////////////////////////////////////////////////
// CE
#define CE_FORCE (COMMON_GAUSS + 1)
#define CE_MMIN  (COMMON_GAUSS + 2)
#define CE_STRATEGY  (COMMON_GAUSS + 3)
#define CE_MAXGB  (COMMON_GAUSS + 4)
#define CE_MAXMEM  (COMMON_GAUSS + 5)
#define CE_TOLIM  (COMMON_GAUSS + 6)
#define CE_TOLRE  (COMMON_GAUSS + 7)
#define CE_TRIALS (COMMON_GAUSS +  8)
#define CE_USEPRIMES  (COMMON_GAUSS + 9)
#define CE_DEPENDENT  (COMMON_GAUSS + 10)
#define CE_APPROXSTEP  (COMMON_GAUSS + 11)
#define CE_APPROXMAXGRID  (COMMON_GAUSS + 12) // alway last of CE_*
#define CE_LAST CE_APPROXMAXGRID 


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
// direct
#define DIRECT_METHOD (COMMON_GAUSS + 1)
#define DIRECT_SVDTOL (COMMON_GAUSS + 2)
#define DIRECT_MAXVAR (COMMON_GAUSS + 3)



///////////////////////////////////////////////////////////////////////
// random poisson hyperplane polygon
#define POLYGON_BETA 0
#define POLYGON_SAFETY 1


///////////////////////////////////////////////////////////////////////
// Families
#define DISTR_DX 0
#define DISTR_PX 1
#define DISTR_QX 2
#define DISTR_RX 3
#define DISTR_NCOL 4
#define DISTR_NROW 5
#define DISTR_ENV 6
#define DISTR_NAME 7
#define DISTR_LAST DISTR_NAME

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


///////////////////////////////////////////////////////////////////////
// Kappas
#define SIZE_NOT_DETERMINED 0 // do not change the value !


///////////////////////////////////////////////////////////////////////
// BASICS
///////////////////////////////////////////////////////////////////////


// posdef \subset negdef  ok 
// methodtype \subset processtype
// posdef, negdef \subset shapetype



/*
  NOTE: change: Besserlupperbound in CovFct.cc,
  empty_idx and header in getNset.cc 
  METHODNAMES, InitModelList and ErrorMessage in initNerror.cc
  addXXX and createmodel in getNset
*/
extern double EIGENVALUE_EPS; // used in GetTrueDim

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

typedef enum TaylorCoeff {TaylorConst, TaylorPow, TaylorExpConst, 
			  TaylorExpPow} TaylorCoeff;
// Const * x^Pow * e^{ExpConst * x^ExpPow} //
// Falls finite range, so wird in TaylorConst der range abgelegt


///////////////////////////////////////////////////////////////////////
// USER CHANGABLE PARAMETERS
///////////////////////////////////////////////////////////////////////
// others


#define MAXINVERSIONS 2


#define generalN 24
// IMPORTANT: all names of general must be at least 3 letters long !!!
extern const char *general[generalN];
#define GENERAL_MODUS 0 
#define GENERAL_STORING 2
#define GENERAL_CPRINT 9
#define GENERAL_EXACTNESS 10
#define GENERAL_CLOSE 15
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
    matrix_inversion[MAXINVERSIONS], // 0:cholesky, 1:QR, 2:SVD; negativ
    seed,
    Ttriple;
 
  double gridtolerance, matrixtolerance, exactness;

} general_param;

#define gaussN 5
extern const char *gauss[gaussN];
typedef struct gauss_param{
  double stationary_only, // logical + NA
    approx_zero;
  bool paired, loggauss;
  int direct_bestvariables;
} gauss_param;

#define krigeN 7
extern const char *krige[krigeN];
typedef struct krige_param{
  char method;
  bool cholesky, ret_variance,  fillall;
  int locmaxn,
    locsplitn[3], // 0:when splitting is done; 1:min pts in a neighbourhood ; 2:max pts when still neighbourhoods are included
    locsplitfactor;
} krige_param;

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


typedef enum InversionMethod { Cholesky, SVD, NoFurtherInversionMethod } 
InversionMethod;
#define directN 3
extern const char *direct[directN];
typedef struct direct_param {
  InversionMethod inversionmethod;
  double svdtolerance;
  int maxvariables;
} direct_param;


#define pnuggetN 1
extern const char * pnugget[pnuggetN];
typedef struct nugget_param {
  double tol;
} nugget_param;


#define sequN 3
extern const char * sequ[sequN];
typedef struct sequ_param{
  int max, back, initial;
} sequ_param;

#define spectralN 4
#define SPECTRAL_PROPFACTOR 2
extern const char * spectral[spectralN];
typedef struct spectral_param {
  bool grid;
  double prop_factor, sigma;
  int lines[MAXTBMSPDIM];
} spectral_param;

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



typedef struct ave_param {
} ave_param;

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

#define hyperN 4
extern const char * hyper[hyperN];
typedef struct hyper_param {
  int superpos, maxlines, mar_distr;
  double mar_param;
} hyper_param;

#define extremeN 10
extern const char * extreme[extremeN];
#define EXTREME_FLAT 5
typedef struct extremes_param {
  int maxpoints, check_every, flat, min_n_zhou, max_n_zhou, mcmc_zhou;
  double standardmax, GEV_xi, density_ratio, eps_zhou;
} extremes_param;

#define brN 9
extern const char * br[brN];
typedef struct br_param {
  int BRmaxmem, BRvertnumber, BRoptimmaxpoints, BRoptim, deltaAM;
  double BRmeshsize, BRoptimtol, variobound, corr_factor;
} br_param;

#define distrN 9
extern const char * distr[distrN];
typedef struct distr_param{
  double safety, minsteplen, innermin, outermax;
  int maxsteps, parts, maxit, mcmc_n, repetitions;
} distr_param;



#define FLAT_UNDETERMINED -1
#define GRIDEXPAND_AVOID -1
#define fitN 43
extern const char * fit[fitN];
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


#define empvarioN 5
extern const char * empvario[empvarioN];
typedef struct empvario_param{
  double phi0, theta0, tol;
  bool pseudovariogram, fft;
} empvario_param;

#define guiN 3
extern const char * gui[guiN];
#define GUI_SIZE 2
typedef struct gui_param{
  bool alwaysSimulate;
  int method, size[2];
} gui_param;

#define graphicsN 9
extern const char *graphics[graphicsN];
#define GRAPHICS_UPTO 3
typedef struct graphics_param{
  bool onefile;
  double height, resolution;
  int PL, increase_upto[2], always_open, always_close, number;
  char filename[100];
} graphics_param;

#define registersN 1
extern const char *registers[registersN];
typedef struct registers_param {
  int keynr;
} registers_param;

#define internalN 18
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
    do_tests, warn_constant, warn_var
    ;// 
} internal_param;

#define coordsN 10
#define COORDS_XYZNOTATION 0
#define ZENIT 8
extern const char *coords[coordsN];
typedef struct coords_param{
  double xyz_notation, zenit[2];
  coord_sys_enum coord_system, new_coord_system;
  char newunits[MAXUNITS][MAXUNITSCHAR], 
    curunits[MAXUNITS][MAXUNITSCHAR],
    varunits[MAXUNITS][MAXUNITSCHAR];
  // 2 user variables for data.frames
  char data_names[MAXUNITS][MAXCHAR], x_names[MAXUNITS][MAXCHAR];
  int data_nr_names, x_nr_names, data_idx[2], x_idx[2];
  bool polar_coord;

} coords_param;

#define specialN 1
extern const char * special[specialN];
typedef struct special_param {
   int multcopies;
} special_param;

#define obsoleteN 5
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
} globalparam;
extern globalparam GLOBAL;





typedef double *coord_type[MAXSIMUDIM];
typedef struct location_type {
  int 
    timespacedim,      /* total dimension, including time and space */
  //length[MAXSIMUDIM], /* grid==true: number of pixels of each direction 
  //		       grid==false: length[0]=number of spatial points in total
  //		       	     */
    spatialdim, 
    xdimOZ; // physical xdim, without time
  long lx, ly, /*
	      the total number of locations of the user's x vector, i.e.
	      * except for grids where lx=3
              * physical length = xdimOZ * lx * (lx - 1) /2 for distances
	      * physical length = xdimOZ * lx else
	     */
    spatialtotalpoints,
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



//typedef struct globalorig{
//  bool set;
//  globalparam param;
//} globalorig;
//extern globalorig GLOBALRIG;



extern bool NAOK_RANGE;



///////////////////////////////////////////////////////////////////////
// GENERAL PARAMETERS FOR THE SIMULAIONMETHODS
///////////////////////////////////////////////////////////////////////

#define Nick(COV) (CovList[(COV)->nr].nick)
#define NAME(COV) CovList[(COV)->nr].name
#define KNAME(NAME) CovList[cov->nr].kappanames[NAME]
#define SNAME(NAME) CovList[cov->nr].subnames[NAME]




#define SET_DESTRUCT(A)\
  assert(meth->S==NULL && meth->destruct==NULL);\
  meth->destruct = A;

typedef struct cov_model cov_model;
typedef struct gen_storage gen_storage;
typedef int (*structfct)(cov_model *cov, cov_model **newcov);
typedef int (*initfct)(cov_model *cov, gen_storage *s);
typedef void (*dofct)(cov_model *cov, gen_storage *s);
typedef void (*do_random_fct)(cov_model *cov, double *v);
typedef void (*param_set_fct)(cov_model *to, cov_model *from, int variant);
typedef bool (*type_fct)(Types required, cov_model *cov, int depth);

//typedef double static_grid[MAXSIMUDIM][3];
typedef enum matrix_type {//TypeMid,
                          TypeMiso, TypeMdiag, 
			  TypeMtimesepproj, // TypeMtimesep and TypeMproj
			  TypeMtimesep, // last column is zero, but last entry
			  TypeMproj, 
			  TypeMany} matrix_type;


typedef enum ptwise_type {pt_posdef, pt_indef, pt_negdef, 
			  pt_zero, pt_paramdep, pt_submodeldep,
			  pt_undefined, pt_unknown, pt_mismatch} 
  ptwise_type;
#define last_pt_definite  pt_zero

typedef struct simu_type {
  bool active, pair;   /* has the init_procedure been called successfully? */
  int expected_number_simu;
} simu_type;


extern const char *METHODNAMES[Forbidden+1],  *ROLENAMES[ROLE_LAST+1], 
  *CAT_TYPENAMES[OtherType + 1],
  *REGNAMES[MODEL_MAX+1], *POSITIVITY_NAMES[pt_mismatch + 1] ;


// extern int SYS_TO_ISO[nr_coord_sys_proj];

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
  bool grid; 
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
typedef void (*logfct)(double *, cov_model*, double* v, double* Sign);
typedef void (*nonstat_logfct)(double *, double*, 
			      cov_model*,  double* v, double* Sign); 
typedef void (*aux_covfct)(double *, double*, double, 
			      cov_model*,  double*); /* x,y, Aux, cov, result */
typedef void (*return_fct)(cov_model*, double*); /* cov, result */ 
typedef void (*return_covmat)(cov_model*, double*, int*); /* cov, result */ 
typedef char (*ext_bool_ret_fct)(cov_model*); /* covt */ 
typedef void (*returnX_covmat)(cov_model*, int *, int, double*, int*);/* cov,..,result */
typedef void (*getlocalparam)(cov_model *, localinfotype *);
typedef bool (*altlocalparam)(cov_model *);
typedef void (*minmaxfct)(cov_model *, double *);


typedef void (*spectral_do)(cov_model *, gen_storage *, double *);

typedef void (*draw_random) (cov_model *cov, double *random);
typedef double (*log_mixdens)(double *x, double logV, cov_model *cov);

typedef void (*sd_fct)(gen_storage *s, cov_model *cov);
		       				   
typedef int (*hyper_pp_fct)(double, double*, double*, cov_model *, bool, 
			    double**, double**, double**);
typedef void (*size_fct)(int i, cov_model *cov, int *nr, int *nc);
typedef sortsofparam  (*paramtype_fct)(int k, int row, int col);

// here all the different method for simulating a RF of a specified covariance
// function are stored

#define ONEARGUMENT_NAME 'k'
#define DEFAULT_SUBNAME 'C'
typedef struct cov_fct {  
  char name[MAXCHAR], nick[MAXCHAR];
  int kappas, /* number of parameters  */
    minsub, maxsub, variants; // ==0
  domain_type domain;
  isotropy_type Isotropy[MAXVARIANTS];
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
  ptwise_type ptwise_definite; 

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
  Types Typi[MAXVARIANTS];
  type_fct TypeFct;

  int density;
  
  double Taylor[MAXTAYLOR][TaylorPow + 1], 
    Tail[MAXTAYLOR][TaylorExpPow + 1]; 
  int TaylorN, TailN;

  return_covmat covmatrix;
  return_fct covariance, inversecovmatrix, variogram,
    pseudovariogram;
  ext_bool_ret_fct is_covariance, is_covmatrix, is_inversecovmatrix, 
    is_variogram, is_pseudovariogram, is_selectedcovmatrix;
  returnX_covmat selectedcovmatrix;
} cov_fct;

extern int currentNrCov, currentRegister;
typedef struct cov_fct cov_fct;
extern cov_fct *CovList;
extern int GENERALISEDCAUCHY, STABLE, BROWNIAN, CAUCHY,  GENNSST_INTERN,
  GAUSS, NUGGET, PLUS, MLEPLUS, MIXEDEFFECT, BALL, ECF, MULT, C_FACTOR,
  DISTRIBUTION,  DETERM_DISTR, GAUSS_DISTR, SETPARAM, SET_DISTR, PROD,
  COVMATRIX, RFGET, COVFCTN,
  DOLLAR, POWER_DOLLAR,  MLE_ENDCOV,  SPLIT, SCATTER, MCMC_PGS, MCMC,
  ISO2ISO,SP2SP, SP2ISO, S2ISO, S2SP, S2S, SId, E2E, E2EIso, 
  E2Sph, E2SphIso, Sph2Sph, Sph2SphIso, LASTGATTER,
  FIRST_PLANE, LAST_PLANE, EARTHKM2CART, EARTHMILES2CART,
  EARTHKM2GNOMONIC, EARTHMILES2GNOMONIC,
  EARTHKM2ORTHOGRAPHIC, EARTHMILES2ORTHOGRAPHIC, 
  NATSC_USER, NATSC_INTERN, TREND, FIX, 
  LOC, SET_DISTR, SCALESPHERICAL, TRUNCSUPPORT, BROWNRESNICK, 
  STROKORB_MONO, STROKORB_BALL_INNER, POLYGON, RECTANGULAR,
  IDCOORD, MULT_INVERSE,
  SHAPESTP, SHAPEAVE, SPHERICAL, UNIF, MPPPLUS, PTS_GIVEN_SHAPE,
  STATIONARY_SHAPE, STANDARD_SHAPE, TRAFO,
  DENSITY, VARIOGRAM_CALL, SIMULATE,
  BRORIGINAL_USER, BRMIXED_USER, BRSHIFTED_USER,
  ARCSQRT_DISTR,SHAPEPOW, POW, SCATTER, GNEITING_INTERN,
  
  BRORIGINAL_INTERN, BRMIXED_INTERN, BRSHIFTED_INTERN,
  MISSING_COV, NULL_MODEL, TBM_OP, USER, 
  DOLLAR_PROC, PLUS_PROC,
  MPPPLUS_PROC, MULT_PROC, TREND_PROC,
  BINARYPROC, BROWNRESNICKPROC, GAUSSPROC, POISSONPROC,
  SCHLATHERPROC, SMITHPROC, CHI2PROC, EXTREMALTPROC, TRENDEVAL,
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


int IncludeModel(const char *name, Types type, int minsub, int maxsub,
		 int kappas,
		  size_fct kappasize,
		  domain_type domain, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref, 
		  bool internal, int vdim, int maxdim, ext_bool finiterange,
		 int monotonicity);
int IncludeModel(const char *name, Types type, int minsub, int maxsub, 
		 int kappas, 
		  domain_type domain, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref, 
		 int maxdim, ext_bool finiterange, int monotonicity);
int CopyModel(const char *name, int which);
int CopyModel(const char *name, int which, Types type);
int CopyModel(const char *name, int which, checkfct check);
void nickname(const char *nick);
void AddVariant(Types type, isotropy_type iso);
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
		return_covmat CovMatrix, ext_bool_ret_fct isCovMatrix,
		return_fct InverseCovMatrix,ext_bool_ret_fct isInverseCovMatrix,
		return_fct Variogram, ext_bool_ret_fct isVariogram,
		return_fct PseudoVariogram, ext_bool_ret_fct isPseudoVariogram,
		returnX_covmat SelectedCovMatrix, 
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
typedef struct ce_storage {
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
} ce_storage;

typedef struct localCE_storage {
    void *correction;
} localCE_storage;

typedef struct approxCE_storage {
  int *idx;
} approxCE_storage;
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
typedef struct tbm_storage {
  // bool genuine_dim[MAX362DIM];
  int ce_dim, simuspatialdim, method, spatialtotalpts, err;
  double center[MAXTBMSPDIM],  
    linesimuscale, linesimufactor, xline_length; 
} tbm_storage;

// see spectral.cc

void metropolis(cov_model *cov, gen_storage *S, double *x);
int search_metropolis(cov_model *cov, gen_storage *S);
// see direct.cc


typedef struct direct_storage {
  InversionMethod method;
  double *U, *G;
} direct_storage;


// see sequential.cc
typedef struct sequ_storage {
  int back, totpnts, spatialpnts, ntime, initial;
  double *U11, *U22, *MuT,  *G,   *Cov21, *Inv22;
  res_type *res0;
} sequ_storage;


// nugget
typedef struct nugget_storage {
    //double sqrtnugget;
  bool simple, simugrid;
  int *pos, reduced_dim[MAXNUGGETDIM], prod_dim[MAXNUGGETDIM + 1];
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



typedef struct br_storage {
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
} br_storage;


typedef struct get_storage {
  cov_model *orig, *get_cov;
  int param_nr, size, vdim[2], *idx;
  bool all;
} get_storage;


typedef struct pgs_storage {
  // urpsprunglich nur fuer pts_given_shape; jedoch allgemein
  // fuer shape functionen
  bool flat, estimated_zhou_c, logmean; 
  double old_zhou, // for mcmc only
  totalmass, // (inverser) Normierungsfaktor, um Raeumliche Funktion
  //            zu einer Dichte zu machen
    currentthreshold, log_density, globalmin, intensity, alpha;
  int 
    own_grid_cumsum[MAXMPPDIM],
    size; //  *len, // global
 
  long double sum_zhou_c, sq_zhou_c;
  long int n_zhou_c;
  double zhou_c; // c in oesting, s, zhoy

  // Huetchen
  double *v, *y;  // local
  coord_type xgr;
  int *pos, // local
     *min, *max;// local dompp
  double *single, *total, *halfstepvector,  // global
    *localmin, *localmax, // local
    *minmean, *maxmean; // standard_shape

 
  // rf_interface.cc
   double
     *supportmin, *supportmax, *supportcentre, // global inkl. dompp
      *own_grid_start, *own_grid_step,
     *own_grid_len; // only for HUETCHEN_OWN_GRIDSIZE
  int *gridlen, *end, *start, *delta, *nx;
  double *xstart, *x, *inc;// local dompp


  // not used in pgs, but in variogramAndCo.cc
  int *endy, *startny, *ptrcol, *ptrrow;
  double *C0x, *C0y, *cross, *z,  **Val;
  // param_set_fct param_set;
  cov_model *cov;

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
    int nstep, *squeezed_dim, *asSign, *i, tmp_n;
} rect_storage;



typedef struct dollar_storage {
  matrix_type type;
  double *z, *y, *z2, *y2, *save_aniso, *inv_aniso;
  int *cumsum, *nx, *total, *len;
  bool simplevar;
} dollar_storage;

typedef struct gatter_storage {
  double *z;
} gatter_storage;

typedef struct earth_storage {
  double 
  *X, *Y, *U, *V, 
    P[9], cart_zenit[3]; // earth2cart u.ae.
} earth_storage;

typedef struct extra_storage {
  double *a, *b, *c;
} extra_storage;

typedef struct biwm_storage {
  bool nudiag_given, cdiag_given;
  double a[3], lg[3], aorig[3], nunew[3],
    scale[4], gamma[4], c[4];
} biwm_storage;


typedef struct inv_storage {
  double *v, *wert;
} inv_storage;


typedef struct scatter_storage{
  int vdim, dim, *min, *max, *nx;
  double *step, *x, *xmin, *value;
} scatter_storage;


typedef struct mcmc_storage{
  //bool done;
  //  int ;
  double integral, posvalue, *pos, *deltapos, *proposed, *propdelta;
} mcmc_storage;


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
   // or of f / g ( Oesting, Sch    lather, Zhou), depending on the 
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


typedef struct gen_storage { // cov_storage, do_storage, init_storage
  // wird in Struct initialisiert, sofern INIT aufgerufen wird;
  // abgespeichert wird es immer im aktuell aufrufenden (fuer next oder key)
  bool check; // used in biWM, BiGneiting
  spec_properties spec;       // used in init
  spectral_storage Sspectral; // used in do
} gen_storage;

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
extern char MSG2[LENERRMSG];
extern char NEWMSG[LENERRMSG];

#define INDENT if (PL >= PL_RECURSIVE) LPRINT("")

#define WARNING1(X, Y) {sprintf(MSG, X, Y); warning(MSG); }
#define ERR(X) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); error(MSG);}
#define ERR1(X, Y) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X); \
    sprintf(MSG2, MSG, Y);					 \
    error(MSG2);}
#define ERR2(X, Y, Z) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X);\
    sprintf(MSG2, MSG, Y, Z);					\
    error(MSG2);}
#define ERR3(X, Y, Z, ZZ) {ERRLINE;sprintf(MSG, "%s %s", ERROR_LOC, X);	\
    sprintf(MSG2, MSG, Y, Z, ZZ);					\
    error(MSG2);}
#define XXERR(X) {/* UERR; */ ERRLINE;errorMSG(X, MSG); 	\
 sprintf(NEWMSG, "in '%s' error %d: %s", ERROR_LOC, X, MSG); error(NEWMSG);}
#define AERR(X) {ERRLINE; PRINTF("AERR: "); errorMSG(X, MSG); 	\
    if (PL<PL_ERRORS) PRINTF("%s%s\n", ERROR_LOC, MSG); assert(false);}
#define MERR(X) {INDENT; PRINTF("error: ");				\
    errorMSG(X, MSG, LENERRMSG);					\
    if (PL<PL_ERRORS) PRINTF("%s%s\n", ERROR_LOC, MSG);}
#define XERR(X) {/* UERR; */ ERRLINE;errorMSG(X, MSG); 	\
    sprintf(NEWMSG, "%s%s", ERROR_LOC, MSG); error(NEWMSG);}
#define PERR(X) {sprintf(MSG, "%s '%s': %s", ERROR_LOC, param_name, X); error(MSG);}
#define PERR1(X,Y) {sprintf(MSG, "%s '%s': %s", ERROR_LOC, param_name, X); sprintf(MSG2, MSG, Y); error(MSG2);}
#define QERR(X) {sprintf(ERRORSTRING, "'%s' %s", param_name, X); DEBUGINFOERR; return ERRORM;}
#define QERRX(ERR, X) {errorMSG(ERR, MSG); sprintf(ERRORSTRING, "'%s' %s (%s)", param_name, X, MSG); DEBUGINFOERR; return ERRORM;}
#define QERRC(NR,X) {sprintf(ERRORSTRING, "%s '%s': %s", ERROR_LOC, CovList[cov->nr].kappanames[NR], X); DEBUGINFOERR; return ERRORM;}
#define QERRC1(NR,X,Y) {sprintf(MSG, "%s '%s': %s", ERROR_LOC, CovList[cov->nr].kappanames[NR], X); sprintf(ERRORSTRING, MSG, KNAME(Y)); DEBUGINFOERR; return ERRORM;}
#define QERRC2(NR,X,Y,Z) {sprintf(MSG, "%s '%s': %s", ERROR_LOC, CovList[cov->nr].kappanames[NR], X); sprintf(ERRORSTRING, MSG, KNAME(Y), KNAME(Z)); DEBUGINFOERR; return ERRORM;}
#define FERR(X) strcpy(ERRORSTRING, X); DEBUGINFOERR;
#define SERR(X) { FERR(X); return ERRORM;}
#define CERR(X) { FERR(X); err=ERRORM; continue;}
#define FERR1(X,Y) sprintf(ERRORSTRING, X, Y); DEBUGINFOERR;  
#define SERR1(X,Y) { FERR1(X, Y); return ERRORM;}
#define CERR1(X,Y) { FERR1(X, Y); err=ERRORM; continue; }
#define SERR2(X, Y, Z) { sprintf(ERRORSTRING, X, Y, Z);  DEBUGINFOERR; return ERRORM;}
#define CERR2(X, Y, Z) { sprintf(ERRORSTRING, X, Y, Z);  err=ERRORM; continue;}
#define SERR3(X, Y, Z, A) { sprintf(ERRORSTRING, X, Y, Z, A); DEBUGINFOERR; return ERRORM;}
#define CERR3(X, Y, Z, A) { sprintf(ERRORSTRING, X, Y, Z, A); err=ERRORM; continue;}
#define SERR4(X, Y, Z, A, B) { sprintf(ERRORSTRING, X, Y, Z, A, B); DEBUGINFOERR;  return ERRORM;}
#define SERR5(X, Y, Z, A, B, C) { sprintf(ERRORSTRING, X, Y, Z, A, B, C); DEBUGINFOERR; return ERRORM;}
#define SERR6(X, Y, Z, A, B, C, D) { sprintf(ERRORSTRING, X, Y, Z, A, B, C, D);  DEBUGINFOERR; return ERRORM;}
#define SERR7(X, Y, Z, A, B, C, D, E) { sprintf(ERRORSTRING, X, Y, Z, A, B, C, D, E); DEBUGINFOERR; return ERRORM;}
#define GERR(X) { strcpy(ERRORSTRING, X); err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR1(X,Y) { sprintf(ERRORSTRING, X, Y); err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR2(X,Y,Z) { sprintf(ERRORSTRING, X, Y, Z); err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR3(X,Y,Z,A) { sprintf(ERRORSTRING, X, Y, Z, A);  err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR4(X,Y,Z,A,B) { sprintf(ERRORSTRING, X, Y, Z, A, B);  err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}
#define GERR5(X,Y,Z,A,B,C) { sprintf(ERRORSTRING, X, Y, Z, A, B, C);  err = ERRORM; DEBUGINFOERR; goto ErrorHandling;}

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
	   double *x, double *T, int spatialdim, int xdim, long lx, bool Time, 
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
	            int xdim, long lx, long ly, bool Time, bool grid,
	            bool distances);

void TaylorCopy(cov_model *to, cov_model *from);

int getmodelnr(char *name);

typedef int (*set_local_q_type) (cov_model *next, double a, double d, 
				  double *q);
void setdefault(cov_model *cov);
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
		    long lx, long ly, bool dist, int xdim, double *T, 
		    bool grid, bool cpy);
//void partial_loc_set(cov_model *cov, double *x, double *y,
//		     long lx, bool dist, int xdim, double *T);
int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdim, long lx, bool Time, bool grid,
	    bool distances,
	    location_type **loc);
int loc_set(double *x, double *y, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdim,
	    long lx, long ly, bool Time, bool grid,
	    bool distances,
	    location_type **Loc);
int loc_set(cov_model *cov, long totalpoints);
//int add_y_zero(location_type *loc);


int Match(char *name, const char * List[], int n);
int Match(char *name, name_type List, int n);

extern name_type STANDARDPARAM, STANDARDSUB;

void CMbuild(SEXP model, int level, cov_model **Cov);
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
void addModel(cov_model *pcov, int subnr, int covnr);
void addModelKappa(cov_model *pcov, int subnr, int covnr);
void addModel(cov_model **pcov, int covnr, cov_model *calling);
void addModel(cov_model **pcov, int covnr, cov_model *calling, bool nullOK); 

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
void expandgrid(coord_type xgr, double **xx, int nrow);
void xtime2x(double *x, int nx, double *T, int len, double **newx, int nrow);
void removeOnly(cov_model **Cov);

SEXP Logi(bool* V, int n, int max) ;
SEXP Num(double* V, int n, int max) ;
SEXP Int(int *V, int n, int max) ;
SEXP Char(const char **V, int n, int max) ;
SEXP Mat(double* V, int row, int col, int max);
SEXP MatInt(int* V, int row, int col, int max) ;
SEXP Array3D(int** V, int depth, int row, int col, int max) ;
SEXP String(char [MAXUNITS][MAXCHAR], int n, int max);


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

extern pref_type PREF_ALL, PREF_NOTHING, PREF_TREND;

#define MAXDEFMATRIX 3
//extern double *userdefinedCovMatrix[MAXDEFMATRIX][MAXMAKEEXPLICITE];
// extern int userdefinedCM_RC[MAXDEFMATRIX][MAXMAKEEXPLICITE];
//extern int CovMatrixRow, CovMatrixCol, CovMatrixTotal;
// extern bool CovMatrixIndexActive;
void zeronugget(int* zeronugget);
				
double *EinheitsMatrix(int dim);
				
extern "C" double *OutX, *OutY;


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
void xA(double *x1, double *x2,  double*A, int nrow, int ncol, double *y1,
	double *y2);
void Ax(double *A, double*x, int nrow, int ncol, double *y);
void Ax(double *A, double*x1, double*x2, int nrow, int ncol, double *y1,
	double *y2);
void AxResType(double *A, res_type *x, int nrow, int ncol, double *y);

SEXP Param(void* p, int nrow, int ncol, SEXPTYPE type, bool drop, long* mem);
int is_positive_definite(double *C, int dim);

// typedef long long POINTER;
typedef long int POINTER;

extern int PrInL;				

void SetTemporarilyGlobal(globalparam *gp);
void ReSetGlobal();
				
int change_coordinate_system(isotropy_type callingprev, isotropy_type isoprev,
			     int tsdim, int xdimprev,
			     int *nr, isotropy_type  *newisoprev,
			     int *newtsdim, int *newxdim);
int SetGatter(domain_type domprev, domain_type statnext, 
	      isotropy_type isoprev, isotropy_type isonext, // bool first, 
	      int *nr, int *delflag);
				
int cpy_add_gatter(cov_model *cov, cov_model **sm);

int binomialcoeff(int n, int k);
void poly_basis(cov_model *cov, gen_storage *s);
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
    vdim[2]; /* for non-quadratic returns */


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

  ptwise_type ptwise_definite;


  ce_storage *Sce;
  localCE_storage *SlocalCE;
  approxCE_storage *SapproxCE;
  direct_storage *Sdirect;
  hyper_storage *Shyper;
  mixed_storage *Smixed;
  nugget_storage *Snugget;
  plus_storage *Splus;
  sequ_storage *Ssequ;
  tbm_storage *Stbm;
  trend_storage *Strend;
  br_storage *Sbr;
  get_storage *Sget;
  pgs_storage *Spgs;
  set_storage *Sset;
  polygon_storage *Spolygon;
  rect_storage *Srect;
  dollar_storage *Sdollar;
  gatter_storage *Sgatter;
  earth_storage *Searth;
  extra_storage *Sextra;
  biwm_storage *Sbiwm;
  inv_storage *Sinv;
  scatter_storage *Sscatter;
  mcmc_storage *Smcmc;
  //check in getNset:  COV_DELETE_WITHOUTSUB, COV_ALWAYS_NULL; add *_DEL, *_NULL
  
  //select_storage *Sselect;
  gen_storage *Sgen; // only once for the whole model tree
} cov_model;


void LOC_DELETE(location_type **Loc);
void LOC_NULL(location_type *loc);
void COV_DELETE(cov_model **cov);
void COV_NULL(cov_model *cov);
void COV_DELETE_WITHOUTSUB(cov_model **Cov);
void COV_DELETE_WITHOUT_LOC(cov_model **Cov);
void ce_NULL(ce_storage* x);
void ce_DELETE(ce_storage **S);
void localCE_NULL(localCE_storage* x);
void localCE_DELETE(localCE_storage**S);
void approxCE_NULL(approxCE_storage* x);
void approxCE_DELETE(approxCE_storage **S);
void direct_NULL(direct_storage  *x);
void direct_DELETE(direct_storage  ** S);
void hyper_NULL(hyper_storage* x);
void hyper_DELETE(hyper_storage  **S);
void nugget_NULL(nugget_storage *x);
void nugget_DELETE(nugget_storage ** S);
void plus_NULL(plus_storage *x);
void plus_DELETE(plus_storage ** S);
void sequ_NULL(sequ_storage *x);
void sequ_DELETE(sequ_storage **S);
void spectral_NULL(sequ_storage *x);
void spectral_DELETE(sequ_storage **S);
void tbm_DELETE(tbm_storage **S); 
void tbm_NULL(tbm_storage* x);
void br_DELETE(br_storage **S); 
void br_NULL(br_storage* x);
void get_NULL(get_storage *S);
void get_DELETE(get_storage **S);
void pgs_DELETE(pgs_storage **S); 
void pgs_NULL(pgs_storage* x);
void set_DELETE(set_storage **S); 
void set_NULL(set_storage* x);
void polygon_DELETE(polygon_storage **S); 
void polygon_NULL(polygon_storage* x);
void rect_DELETE(rect_storage **S); 
void rect_NULL(rect_storage* x);
void dollar_DELETE(dollar_storage **S); 
void dollar_NULL(dollar_storage* x);
void gatter_DELETE(gatter_storage **S); 
void gatter_NULL(gatter_storage* x);
void earth_DELETE(earth_storage **S); 
void earth_NULL(earth_storage* x);
void extra_DELETE(extra_storage **S); 
void extra_NULL(extra_storage* x);
void biwm_DELETE(biwm_storage **S); 
void biwm_NULL(biwm_storage* x);
void inv_DELETE(inv_storage **S);
void inv_NULL(inv_storage* x);
void scatter_DELETE(scatter_storage **S);
void scatter_NULL(scatter_storage* x);
void mcmc_DELETE(mcmc_storage **S);
void mcmc_NULL(mcmc_storage* x);
void trend_DELETE(trend_storage ** S);
void trend_NULL(trend_storage* x);
void mixed_DELETE(mixed_storage **S);
void mixed_NULL(mixed_storage* x);
void gen_NULL(gen_storage *x);
void gen_DELETE(gen_storage **S);

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

void ErrLogCov(double *x, cov_model *cov, double *v, double *Sign);
void ErrLogCovNonstat(double *x, double *y, cov_model *cov, double *v, 
		      double *Sign);
int struct_failed(cov_model *cov, cov_model **atom);
int init_failed(cov_model *cov, gen_storage *s);
int init_statiso(cov_model *cov, gen_storage *s);
void do_failed(cov_model *cov, gen_storage *s);
void do_statiso(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
int checkOK(cov_model *cov);
int checkplusmal(cov_model *cov);
void ErrInverse(double *v, cov_model *cov, double *x);
void InverseIsotropic(double *U, cov_model *cov, double *inverse);
bool verysimple(cov_model *cov);


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
int StructBR(cov_model *cov, gen_storage *s, cov_model **atom);

sortsofparam uncritical_paramtype(int k, int row, int col);

int check_within_range(cov_model *cov, bool NAOK);
extern bool RELAX_UNKNOWN_RFOPTION;


int INIT_RANDOM_intern(cov_model *M, int moments, gen_storage *s, double *p);
int INIT_intern(cov_model *M, int moments, gen_storage *s);
int REINIT_intern(cov_model *M, int moments, gen_storage *s);
void set_initialised_false(cov_model *cov, bool init_deterministic);
  
int alloc_mpp_M(cov_model *cov, int moments);
void free_mpp_M(cov_model *cov);

int alloc_pgs(cov_model *cov, int dim);
int alloc_pgs(cov_model *cov);
int alloc_cov(cov_model *cov, int dim, int rows, int cols);


#define DO(Cov, S) {assert(!TypeConsistency(RandomType, Cov, 0)); \
    assert(Cov->initialised);				       \
    ASSERT_GATTER(Cov);					       \
    PL--;						       \
    CovList[(Cov)->gatternr].Do(Cov, S);		       \
    PL++;						       \
  }
   
#define DORANDOM(Cov, S) {\
    assert(TypeConsistency(RandomType, Cov, 0));			      \
    ASSERT_GATTER(Cov);						      \
    PL--;							      \
    CovList[(Cov)->gatternr].DoRandom(Cov, S);			      \
    PL++;							      \
  }
#define COV(X, Cov, V) { ASSERT_GATTER(Cov); CovList[(Cov)->gatternr].cov(X, Cov, V);}
#define LOGCOV(X, Cov, V, S) {ASSERT_GATTER(Cov);CovList[(Cov)->gatternr].log(X, Cov, V, S);}
#define SHAPE COV
#define FCTN COV
#define ABSFCTN(X, Cov, V) { COV(X, Cov, V); *(V) = fabs(*(V)); }
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

//#define DENSITYFCT(X, Dens, V) CovList[(Cov)->gatternr].approx_dens(X, Dens, V)
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
#define ROLE_ASSERT_HYPER if (cov->role != ROLE_HYPER) ILLEGAL_ROLE 
  
#define ASSERT_GAUSS_METHOD(METHOD) DEBUGINFO;				\
  if(cov->role != ROLE_GAUSS || cov->method != METHOD) 			\
    SERR4("Gaussian field for '%s' only possible with '%s' as method. Got role '%s' and method '%s'.", \
	 NICK(cov),							\
	 CovList[METHOD == RandomCoin ? RANDOMCOIN_USER	:		\
		 gaussmethod[METHOD] -					\
		 CovList[gaussmethod[METHOD]].internal].nick,		\
	  ROLENAMES[cov->role],						\
	  CovList[cov->method == RandomCoin ? RANDOMCOIN_USER	:	\
		  gaussmethod[cov->method] -				\
		  CovList[gaussmethod[cov->method]].internal].nick	\
	  )

#define ASSERT_NEWMODEL_NOT_NULL\
  if (newmodel == NULL) {\
    SERR1("unexpected call of struct_%s", NAME(cov));\
  }

#define ASSERT_NEWMODEL_NULL\
  if (newmodel != NULL) {\
    SERR1("Unexpected call of struct_%s", NAME(cov));\
  }

#define ASSERT_ONE_SUBMODEL(COV)					\
  DEBUGINFO;								\
  if (!(((COV)->sub[0] == NULL) xor ((COV)->sub[1] == NULL)))		\
    SERR2("either '%s' or '%s' must be given", CovList[(COV)->nr].subnames[0],\
      CovList[(COV)->nr].subnames[1])

#define ASSERT_ROLE_DEFINED(SUB) DEBUGINFO;			\
  if (role == ROLE_UNDEFINED)					\
    SERR1("'%s' not allowed as shape function.", NICK(SUB));	\

void PSTOR(cov_model *cov, gen_storage *x);
int structOK(cov_model *cov, cov_model **newmodel);
int initOK(cov_model *cov, gen_storage *s);
void doOK(cov_model *cov, gen_storage *s);
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


#define PLE PRINTF("\n(PLE '%s', line %d)", __FILE__, __LINE__); ple_
void ple_(cov_model *cov);
void ple_(char *name);



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

bool isPosDef(cov_model *cov);
bool isNegDef(cov_model *cov);
bool isVariogram(cov_model *cov);
bool isProcess(cov_model *cov);
bool isGaussMethod(cov_model *cov);
bool isBRuserProcess(cov_model *cov);
bool isPointShape(cov_model *cov);
bool isRandom(cov_model *cov);
bool isShape(cov_model *cov);
bool isTrend(cov_model *cov);
bool isInterface(cov_model *cov);
bool isUndefined(cov_model *cov);
bool isOther(cov_model *cov);

bool isTcf(Types type);
bool isPosDef(Types type);
bool isVariogram(Types type);
bool isNegDef(Types type);
//bool isDefinite(Types type);
bool isProcess(Types type);
bool isGaussMethod(Types type);
bool isBRuserProcess(Types type);
bool isPointShape(Types type);
bool isRandom(Types type);
bool isShape(Types type);
bool isTrend(Types type);
bool isInterface(Types type);
bool isUndefined(Types type);
bool isOther(Types type);
bool isUndefined(cov_fct *C);

bool isMonotone(cov_model *cov);
bool isMonotone(int montone);
bool isNormalMixture(int montone);
bool isNormalMixture(cov_model *cov);
bool isBernstein(cov_model *cov);
bool isBernstein(int monotone);
bool isCompletelyMonotone(int monotone);
bool isGneiting(cov_model *cov);
bool isGneiting(int monotone);

bool isGaussProcess(cov_model *cov);
bool isGaussMethod(cov_model *cov);
bool isBernoulliProcess(cov_model *cov);
bool isGaussBasedProcess(cov_model *cov);
bool isBrownResnickProcess(cov_model *cov);
bool isMaxStable(cov_model *cov);

bool isCov(cov_model *cov);
typedef bool (*typusfct)(Types type); /* h, cov, result */ 
bool is_any(typusfct t, cov_fct *C);
bool is_all(typusfct t, cov_fct *C);


bool hasNoRole(cov_model *cov);
bool hasMaxStableRole(cov_model *cov);
bool hasExactMaxStableRole(cov_model *cov);
bool hasPoissonRole(cov_model *cov);
bool hasAnyShapeRole(cov_model *cov);
bool hasDistrRole(cov_model *cov);

bool TypeConsistency(Types requiredtype, Types deliveredtype);
int TypeConsistency(Types requiredtype, cov_model *cov, int depth);
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



#define NEW_COV_STORAGE(cov, new) {		\
  if ((cov)->S##new != NULL) {					\
    new##_DELETE(&((cov)->S##new));					\
    assert((cov)->S##new == NULL);					\
  }								\
  if ((cov)->S##new == NULL) {					\
    (cov)->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    new##_NULL((cov)->S##new);					\
    if ((cov)->S##new == NULL) BUG;				\
  }}								

#define NEW_STORAGE(new)	\
  NEW_COV_STORAGE(cov, new)


#define COND_NEW_STORAGE(new, WHAT) {	\
  if (cov->S##new != NULL && cov->S##new->WHAT != NULL) {		\
    new##_DELETE(&(cov->S##new));					\
    assert(cov->S##new == NULL);					\
  }								\
  if (cov->S##new == NULL) {					\
   cov->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    new##_NULL(cov->S##new);					\
    if (cov->S##new == NULL) BUG;					\
  }								\
  assert(cov->S##new->WHAT == NULL);				\
  }


#define ALLOC_NEW(Snew, Z, SIZE, WHAT)					\
  assert(cov->Snew != NULL);						\
  double *Z = cov->Snew->WHAT;						\
  if (Z == NULL)							\
    Z = cov->Snew->WHAT = (double*) MALLOC(sizeof(double) * (SIZE)) 
  

#define ALLOC_NEWINT(Snew, Z, SIZE, WHAT)				\
  assert(cov->Snew != NULL);						\
  int *Z = cov->Snew->WHAT;						\
  if (Z == NULL)							\
    Z = cov->Snew->WHAT = (int*) MALLOC(sizeof(int) * (SIZE)) 


#define EXTRA_STORAGE COND_NEW_STORAGE(extra, a)
#define ALLOC_EXTRA(Z, SIZE) ALLOC_NEW(Sextra, Z, SIZE, a) 
#define ALLOC_EXTRA2(Z, SIZE) ALLOC_NEW(Sextra, Z, SIZE, b)     
#define ALLOC_EXTRA3(Z, SIZE) ALLOC_NEW(Sextra, Z, SIZE, c)     


#define DOLLAR_STORAGE COND_NEW_STORAGE(dollar, z)
#define ALLOC_DOLLAR(Z, SIZE) ALLOC_NEW(Sdollar, Z, SIZE, z) 
#define ALLOC_DOLLARY(Z1, Z2, SIZE) \
     ALLOC_DOLLAR(Z1, SIZE); ALLOC_NEW(Sdollar, Z2, SIZE, y)

#define ALLOC_DOLLAR2(Z, SIZE) ALLOC_NEW(Sdollar, Z, SIZE, z2)     
#define ALLOC_DOLLAR3(Z, SIZE) ALLOC_NEW(Sdollar, Z, SIZE, y)     
#define ALLOC_DOLLAR4(Z, SIZE) ALLOC_NEW(Sdollar, Z, SIZE, y2)     


#define WARN_NEWDEFINITIONS						\
  if (GLOBAL.internal.warn_new_definitions) {				\
  warning("Note that in Version 3.0.33 some definitions have changed (and some typos corrected), see 'RMbernoulli', 'RMbrownresnick', 'RMbr2bg' and 'RMbr2eg'.\nNote that in Version 3.0.43 some typos have been corrected in 'RMS' influencing the result."); \
  GLOBAL.internal.warn_new_definitions = false;				\
  }


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


#define AVL_FUNC_TYPES 1
typedef struct cell_type {
    unsigned int *code;
    res_type colour;
} cell_type;
typedef int (*avl_comparison_func) (cell_type *, cell_type *, int *);
typedef void (*avl_node_func) (cell_type*, int *);
typedef cell_type *(*avl_copy_func) (cell_type*, int *);

bool CallingSet(cov_model *cov);
void includeStandardMath();

#define I_COL_NA -1

void setptwise(ptwise_type pt);


#ifdef LOCAL_MACHINE
#define STOPAFTER(COND, DO) 				\
  static bool __stopafter__ = false;			\
  if (__stopafter__ && !(COND)) { DO; BUG;}		\
  __stopafter__  = COND;				
#else 
#define STOPAFTER(COND, DO) 	
#endif
  							
    
bool TrafoOK(cov_model *cov, bool all);

#endif
