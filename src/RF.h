

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

#ifndef RFsimu_H
#define RFsimu_H 1

// ACHTUNG : REIHENFOLGE WICHTIG
#include <Basic_utils.h>
#include "basic.h"

//
// 1
//// 1



#define showfree !true 
#define DOPRINT true


// in /home/schlather/TMP/RandomFieldsUtils/include/utils.h:
// #define MEMCOPY(A,B,C) memory_copy(A, B, C)

#include <string.h>
#include  <utils.h>
#include "error.h"
#include "auxiliary.h"
#include "RandomFields.h"
#include "AutoRandomFields.h"
#include "xport.h"
#include "Userinterfaces.h"

//   intptr_t and uintptr_t fuer: umwandlung pointer in int und umgekehrt

#ifndef RFERROR
#define RFERROR error
#endif



//////////////////////////////////////////////////////////////////////
// CHECKING 
//////////////////////////////////////////////////////////////////////
//#define ASSERT_GATTERONLY(Cov) assert(TrafoOK(Cov, false))
#define ASSERT_GATTER(Cov) assert(TrafoOK(Cov))
#define ASSERT_CHECKED(Cov) assert((Cov)->checked)
#ifdef SCHLATHERS_MACHINE
#define _no_loc_ SERR2("locations not initialised (%s line %d).", __FILE__, __LINE__)
//#undef NULL
//#define NULL nullptr
#else
#define _no_loc_ SERR("locations not initialised.")
#endif

#define ASSERT_LOC_GIVEN if (loc == NULL) {PMI(cov); _no_loc_;}

#ifdef __GNUC__
#define HIDE_UNUSED_VARIABLE 1
#endif


#ifndef SCHLATHERS_MACHINE
#define SET_IDX(P, IDX) (GLOBAL.general.set % (P)->nrow[IDX])

#define CHECK(C,T,X,type,D,I,V,R) check2X(C,T,X,type,D,I,V,R)
#define CHECK_NO_TRAFO(C,T,X,type,D,I,V,R) check2Xnotrafo(C,T,X,type,D,I,V,R)
#define STRUCT(Cov, NM) CovList[(Cov)->gatternr].Struct(Cov, NM)

#define PARAM(P, IDX) ((double *) (P)->px[IDX])
#define PARAMINT(P, IDX) ((int *) (P)->px[IDX])
#define PARAMVEC(P, IDX) ((sexp_type *) (P)->px[IDX])->sexp
#define PARAMENV(P, IDX) ((sexp_type *) (P)->px[IDX])
#define PARAMLANG(P, IDX) ((sexp_type *) (P)->px[IDX])
#define PARAMLIST(P, IDX) ((listoftype *) (P)->px[IDX])
#define LPARAM(P, IDX) ((double *) (PARAMLIST(P, IDX)->lpx[SET_IDX(P, IDX)]))
#define LPARAMINT(P, IDX) ((int *) (PARAMLIST(P, IDX)->lpx[SET_IDX(P, IDX)]))
//#define LISTLIST(P, IDX) ((listoftype *) (P)->px[IDX])

#define PARAM0(P, IDX) PARAM(P, IDX)[0]
#define PARAM0INT(P, IDX) PARAMINT(P, IDX)[0]
#define LPARAM0(P, IDX) LPARAM(P, IDX)[0]
#define LPARAM0INT(P, IDX) LPARAMINT(P, IDX)[0]

#define PCOPY(TO, FROM, IDX) 						\
  MEMCOPY((TO)->px[IDX], (FROM)->px[IDX],				\
	  ((FROM)->nrow[IDX]) * ((FROM)->ncol[IDX]) *			\
	  (CovList[(FROM)->nr].kappatype[IDX]==REALSXP ? sizeof(double) : \
	   CovList[(FROM)->nr].kappatype[IDX]==INTSXP ? sizeof(int) :	\
	    -1))

#define QcovALLOC(cov, nr) {\
    (cov)->qlen = (int) (nr);						\
  assert((cov)->q == NULL);						\
  if (((cov)->q = (double*) CALLOC((int) (nr), sizeof(double))) == NULL) \
    ERR("memory allocation error for local memory");		\
}

#define B2D(X) (double) (X)
#define B2I(X) (int) (X)
#define E2D(X) (X)
#define E2I(X) (int) (X)
#define E2B(X) (bool) (X)
#define I2D(X) (double) (X)

////////////////////////////////////////////////////////////////////////////
#else // SCHLATHERS_MACHINE
#define SET_IDX(P, IDX) __extension__({assert((P) != NULL && (P)->nrow[IDX] >=0 && (P)->nrow[IDX] <= 1000000); (GLOBAL.general.set % (P)->nrow[IDX]);})

#define CHECK(C,T,X,type,D,I,V,R) __extension__({assert((type)!=RandomType); check2X(C,T,X,type,D,I,V,R);})
#define CHECK_NO_TRAFO(C,T,X,type,D,I,V,R) __extension__({assert((type)!=RandomType); check2Xnotrafo(C,T,X,type,D,I,V,R);})
#define STRUCT(C, NM)  __extension__({ASSERT_GATTER(C);  CovList[(C)->gatternr].Struct(C, NM);})

// #define NN
#define PARAM(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == REALSXP);  (P)->px[IDX];})
#define PARAMINT(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == INTSXP);  ((int *) (P)->px[IDX]);})
//#define PARAMVEC(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == VECSXP);  ((SEXP) (P)->px[IDX]);})
#define PARAMVEC(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == VECSXP);  ((sexp_type *) (P)->px[IDX])->sexp;})
#define PARAMENV(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == ENVSXP);  ((sexp_type *) (P)->px[IDX]);})
#define PARAMLANG(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == LANGSXP);  ((sexp_type *) (P)->px[IDX]);})
#define PARAMLIST(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] >= LISTOF);  ((listoftype *) (P)->px[IDX]);})
#define LPARAM(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == LISTOF + REALSXP); assert(PARAMLIST(P, IDX) != NULL); ((double *) (PARAMLIST(P, IDX)->lpx[SET_IDX(P, IDX)]));})
#define LPARAMINT(P, IDX) __extension__({assert(CovList[(P)->nr].kappatype[IDX] == LISTOF + INTSXP); assert(PARAMLIST(P, IDX) != NULL); ((int *) (PARAMKLIST(P, IDX)->lpx[SET_IDX(P, IDX)]));})
//#define LISTLIST(P, IDX) __extension__({assert((P)->px != NULL); (freelistoftype *) (P)->px[IDX];})

#define PARAM0(P, IDX) __extension__({assert(PARAM(P, IDX) != NULL); PARAM(P, IDX)[0];})
#define PARAM0INT(P, IDX)  __extension__({assert(PARAMINT(P, IDX) != NULL); PARAMINT(P, IDX)[0];})
#define LPARAM0(P, IDX) __extension__({assert(LPARAM(P, IDX) != NULL); LPARAM(P, IDX)[0];})
#define LPARAM0INT(P, IDX)  __extension__({assert(LPARAMINT(P, IDX) != NULL); LPARAMINT(P, IDX)[0];})


#define PCOPY(TO, FROM, IDX) 					\
  {if (!((TO) != NULL && (FROM) != NULL &&			\
	 (TO)->px[IDX] != NULL && (FROM)->px[IDX] &&IDX > 0 &&	\
	 (FROM)->ncol[IDX] == (TO)->ncol[IDX] &&		\
	 (FROM)->nrow[IDX] == (TO)->nrow[IDX] &&			\
	 CovList[(FROM)->nr].kappatype[IDX]==CovList[(TO)->nr].kappatype[IDX]&&\
	 (CovList[(FROM)->nr].kappatype[IDX]==REALSXP ||		\
	  CovList[(FROM)->nr].kappatype[IDX]==INTSXP)  )) BUG;		\
    MEMCOPY((TO)->px[IDX], (FROM)->px[IDX],				\
	    ((FROM)->nrow[IDX]) * ((FROM)->ncol[IDX]) *			\
	    (CovList[(FROM)->nr].kappatype[IDX]==REALSXP ? sizeof(double) : \
	     CovList[(FROM)->nr].kappatype[IDX]==INTSXP ? sizeof(int) :	\
	     -1));							\
  }

#define QcovALLOC(cov, nr) {\
    assert((cov)->q == NULL);						\
    (cov)->qlen = (int) (nr);						\
    if (((cov)->q = (double*) CALLOC((int) (nr), sizeof(double))) == NULL) \
    ERR("memory allocation error for local memory");		\
}

#define B2D(X) ((X) != true) && ((X) != false) ? ({BUG; 0.0;}) : (double) (X)

#endif // SCHLATHERS_MACHINE

#define SET_OUT_OF(D) (D)->lpx[GLOBAL.general.set]
#define NROW_OUT_OF(D) (D)->nrow[GLOBAL.general.set]
#define NCOL_OUT_OF(D) (D)->ncol[GLOBAL.general.set]



//////////////////////////////////////////////////////////////////////
// DEDUGGING INFORMATION
//////////////////////////////////////////////////////////////////////

#define KPRINT leer(PrInL);Rprintf
#define LPRINT {cov_model *lprint_z=cov; int lprint_i=0; while (lprint_z->calling != NULL && lprint_i<10) {lprint_z=lprint_z->calling; if (DOPRINT) {Rprintf(DOT); Rprintf(" ");} lprint_i++;} if (lprint_i==100) {Rprintf("LPRINT i=%d\n", lprint_i);PMI(cov); assert(false);}} if (DOPRINT) Rprintf


#ifndef RANDOMFIELDS_DEBUGGING

#define COV_DELETE COV_DELETE_

#define NICK(COV) (isDollar(COV) ? CovList[(COV)->sub[0]->nr].nick : CovList[(COV)->nr].nick)
 

#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,R) check2X(C,T,X,type,D,I,V0,V1,R,true)
#define CHECKPD2ND(N,D1,D2,I,V,R) CheckPD2ND(N,D1,D2,I,V,R)
#define INIT(M, Moments, S) INIT_intern(M, Moments, S)
#define REINIT(M, Moments, S) REINIT_intern(M, Moments, S)
#define INIT_RANDOM(M, Moments, S, P) INIT_RANDOM_intern(M, Moments, S, P)
 



#else // of RANDOMFIELDS_DEBUGGING

#define COV_DELETE(COV) {					\
    PRINTF("\nCOV_DELETE: '%s', line %d", __FILE__, __LINE__);	\
    COV_DELETE_(COV);}

#define NICK(COV) (CovList[(COV)->nr].nick)


#define XX(C) assert((C)->simu.expected_number_simu >= 0|| __extension__({DOPRINTF("Start.\n"); false;}))
#define YY(C) assert((C)->simu.expected_number_simu >= 0 || __extension__({DOPRINTF("End.\n"); false;}))

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
#undef CHECK
#undef CHECK_NO_TRAFO
#define CHECK(C,T,X,type,D,I,V,R) ({LLPRINT(CHECKSIGN, C, "CHECK"); assert((type)!=RandomType); XX(C); int _x = check2X(C,T,X,type,D,I,V,R); YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED"); errorMSG(_x, MSG); PRINTF("%s%s\n", ERROR_LOC, MSG);} _x;})
#define CHECK_NO_TRAFO(C,T,X,type,D,I,V,R) ({LLPRINT(CHECKSIGN, C, "CHECK"); assert((type)!=RandomType); XX(C); int _x = check2Xnotrafo(C,T,X,type,D,I,V,R); YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED"); errorMSG(_x, MSG); PRINTF("%s%s\n", ERROR_LOC, MSG);} _x;})
#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,R) ({LLPRINT(CHECKSIGN, C, "CHECKVDIM"); XX(C); int _x = check2X(C,T,X,type,D,I,V0,V1,R,true);YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED");} _x;})
#define CHECKPD2ND(C,D1,D2,I,V,R) ({LLPRINT(CHECKSIGN, C, "CHECKPD2ND"); XX(C); int _x = CheckPD2ND(C,D1,D2,I,V,R);YY(C); if (_x==NOERROR){LLPRINT(CHECKSIGN, C, "CHECK DONE");}else{LLPRINT(CHECKSIGN, C, "CHECK FAILED");} _x;})
#define INIT(C, Moments, S) ({LLPRINT(INITSIGN, C, "INIT");  XX(C); int _x = INIT_intern(C, Moments, S);YY(C); if (_x==NOERROR){LLPRINT(STRUCTSIGN, C, "INIT DONE");}else{LLPRINT(STRUCTSIGN, C, "INIT FAILED");}_x;})
#define REINIT(C, Moments, S) ({LLPRINT(INITSIGN, C, "INIT");  XX(C); int _x = REINIT_intern(C, Moments, S); YY(C); _x;})
#define INIT_RANDOM(C, Moments, S, P) ({LLPRINT(INITSIGN, C, "INITRANDOM");  XX(C); int _x = INIT_RANDOM_intern(C, Moments, S, P);YY(C); if (_x==NOERROR){LLPRINT(STRUCTSIGN, C, "INIT DONE");}else{LLPRINT(STRUCTSIGN, C, "INIT FAILED");}_x;})
#undef STRUCT
#define STRUCT(C, NM)  ({LLPRINT(STRUCTSIGN, C, "STRUCT"); ASSERT_GATTER(C);  XX(C); int _x = CovList[(C)->gatternr].Struct(C, NM);YY(C); if (_x==NOERROR){LLPRINT(STRUCTSIGN, C, "STRUCT DONE");}else{LLPRINT(STRUCTSIGN, C, "STRUCT FAILED");}_x;})

   /*  printf("//%ld %ld %d %d %d idx=%d\n", (TO)->px[IDX], (FROM)->px[IDX], \
       (FROM)->nrow[IDX], (FROM)->ncol[IDX],				\
       CovList[(FROM)->nr].kappatype[IDX]==REALSXP ? sizeof(double) :	\
       CovList[(FROM)->nr].kappatype[IDX]==INTSXP ? sizeof(int) :	\
       -1, IDX);							\
				 */


#endif //  RANDOMFIELDS_DEBUGGING



#define P(IDX) PARAM(cov, IDX) 
#define PINT(IDX) PARAMINT(cov, IDX)
#define P0(IDX) PARAM0(cov, IDX) 
#define P0INT(IDX) PARAM0INT(cov, IDX)
#define PVEC(IDX) PARAMVEC(cov, IDX) 
#define PENV(IDX) PARAMENV(cov, IDX) 
#define PLANG(IDX) PARAMLANG(cov, IDX) 
#define PLIST(IDX) PARAMLIST(cov, IDX)
#define PARAMSEXP(P, IDX) ((sexp_type *) (P)->px[IDX]) /* kein assert! */
#define PSEXP(IDX) PARAMSEXP(cov, IDX)
#define NROW(IDX) cov->nrow[IDX]
#define NCOL(IDX) cov->ncol[IDX]
#define LP(IDX) LPARAM(cov, IDX) 
#define LPINT(IDX) LPARAMINT(cov, IDX)
#define LP0(IDX) LPARAM0(cov, IDX) 
#define LP0INT(IDX) LPARAM0INT(cov, IDX)
#define LNRC_(I, rc) PLIST(I)->rc[SET_IDX(cov, I)]
#ifndef SCHLATHERS_MACHINE
#define LNRC(I, rc) LNRC_(I, rc)
#else
#define LNRC(I,rc) __extension__ ({assert(PLIST(I) != NULL); PLIST(I)->rc[SET_IDX(cov, I)];})
#endif
#define LNROW(IDX) LNRC(IDX, nrow)
#define LNCOL(IDX) LNRC(IDX, ncol)
#define LPARAMNROW(P, IDX) PARAMLIST(P, IDX)->nrow[SET_IDX(P, IDX)]
#define LPARAMNCOL(P, IDX) PARAMLIST(P, IDX)->ncol[SET_IDX(P, IDX)]

#define PARAMFREE(P, IDX) if ((P)->px[IDX] != NULL) {	\
  if (CovList[(P)->nr].kappatype[IDX] < LISTOF) {		\
    UNCONDFREE((P)->px[IDX]);				\
  } else {						\
    LIST_DELETE( (listoftype **) ((P)->px + IDX));	\
  }							\
  (P)->ncol[IDX] = (P)->nrow[IDX] = 0;			\
  }


#define PARAMALLOC(P, IDX, ROW, COL) {					\
    int _PARAMsize;							\
    switch(CovList[(P)->nr].kappatype[IDX]) {				\
    case REALSXP : _PARAMsize = sizeof(double); break;			\
    case INTSXP : _PARAMsize = sizeof(int); break;			\
    default : 						\
      if ((P)->kappasub[IDX]!=NULL && (P)->kappasub[IDX]->nr==DISTRIBUTION) { \
	  ERR("argument value recognized as distribution family although it should not. Maybe the error is caused by a non-existing variable."); \
      } else BUG;							\
    }									\
    assert((P)->px[IDX]==NULL);						\
    (P)->nrow[IDX] = ROW; (P)->ncol[IDX] = COL;				\
    if (((P)->px[IDX] =							\
	 (double*) CALLOC((ROW) * (COL), _PARAMsize)) == NULL) {	\
      XERR(ERRORMEMORYALLOCATION)					\
    }									\
  }
#define PALLOC(IDX, ROW, COL) PARAMALLOC(cov, IDX, ROW, COL)

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
#define ROLE_LIKELIHOOD 12
#define ROLE_FAILED 13
#define ROLE_UNDEFINED 14
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

#define MAX_NA 30
//#define MAX_MLE_ELMNTS 10
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
#define MAXLILIGRIDDIM 10


#define NATSCALE_EXACT 1 /* or approx or mle */
#define NATSCALE_ORNUMERIC 2
//#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check fitvario when changing; best not to change !!! */


#define LISTOF 1000  /* may not interfere with define of *SXP in Rinternal.h*/


extern bool LOCAL_DEBUG;


// #define MAXINTERNALPARAM 2
extern double ZERO[MAXSIMUDIM], ONE;
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

#define PARSnudiag 0



///////////////////////////////////////////////////////////////////////
// stp
#define STP_S 0
#define STP_Z 1
#define STP_M 2


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
#define TREND_FCT 4 // CLOSXP
#define TREND_PARAM_FCT 5

#define MIXED_X 0
#define MIXED_BETA 1
#define MIXED_COORD 2
#define MIXED_DIST 3
#define MIXED_DIM 4

///////////////////////////////////////////////////////////////////////
// Covariate models
#define COVARIATE_C 0
#define COVARIATE_X 1
#define COVARIATE_RAW 2
#define COVARIATE_ADDNA 3
#define COVARIATE_FACTOR 4

///////////////////////////////////////////////////////////////////////
// fixed covariance
#define FIXCOV_M COVARIATE_C
#define FIXCOV_X COVARIATE_X 
#define FIXCOV_RAW COVARIATE_RAW




///////////////////////////////////////////////////////////////////////
// QAM
#define QAM_THETA 0

///////////////////////////////////////////////////////////////////////
// NSST
#define NSST_DELTA 0

///////////////////////////////////////////////////////////////////////
// nugget
#define NUGGET_TOL 0
#define NUGGET_VDIM 1
#define NUGGET_PROC_TOL (COMMON_GAUSS + 1)
#define NUGGET_PROC_VDIM (COMMON_GAUSS + 2)


///////////////////////////////////////////////////////////////////////
//
#define RFOPTIONS "RFoptions"


///////////////////////////////////////////////////////////////////////
// gengneiting
#define GENGNEITING_K 0
#define GENGNEITING_MU 1


///////////////////////////////////////////////////////////////////////
// bigneiting
#define GNEITING_K GENGNEITING_K    // important to keep !
#define GNEITING_MU 1
#define GNEITING_S 2
#define GNEITING_SRED 3
#define GNEITING_GAMMA 4
#define GNEITING_CDIAG 5
#define GNEITING_RHORED 6
#define GNEITING_C 7

///////////////////////////////////////////////////////////////////////
// biwm
#define BInudiag 0
#define BInured 1
#define BInu 2
#define BIs 3
#define BIcdiag 4
#define BIrhored 5
#define BIc 6
#define BInotinvnu 7


///////////////////////////////////////////////////////////////////////
// Power, shapepow
#define POW_ALPHA 0

///////////////////////////////////////////////////////////////////////
// trafo
#define TRAFO_ISO 0


///////////////////////////////////////////////////////////////////////
// general for many R.fctns
#define MATH_FACTOR 2

/////////////////////////////////////
// for model RMderivative
#define PARTIAL 0
#define ALGEBRA 1

/////////////////////////////////////
// for model RMhelmholtz
#define COMPONENT  0
#define ANISO_VAR 1


///////////////////////////////////////////////////////////////////////
// constant
#define CONST_C 0


///////////////////////////////////////////////////////////////////////
// R.is
#define IS_IS 1


///////////////////////////////////////////////////////////////////////
// R.p
#define PROJ_PROJ 0
#define PROJ_ISO 1
#define PROJ_FACTOR 2


///////////////////////////////////////////////////////////////////////
// Scatter
#define SCATTER_STEP 0
#define SCATTER_MAX 1


///////////////////////////////////////////////////////////////////////
// likelihood
#define LIKELIHOOD_DATA 0
#define LIKELIHOOD_NA_VAR 1
#define LIKELIHOOD_BETASSEPARATE 2
#define LIKELIHOOD_IGNORETREND 3
#define LIKELIHOOD_LAST LIKELIHOOD_IGNORETREND
#define LIKELI_NA_INTEGER false
#define LIKELI_EXCLUDE_TREND true


///////////////////////////////////////////////////////////////////////
// gaussproc -- GENERAL
#define GAUSS_BOXCOX 0
#define COMMON_GAUSS 0


///////////////////////////////////////////////////////////////////////
// gaussproc
#define GAUSSPROC_STATONLY  (COMMON_GAUSS + 1)
#define GAUSSPROC_LAST GAUSSPROC_STATONLY


///////////////////////////////////////////////////////////////////////
// coins
#define COIN_COV 0
#define COIN_SHAPE 1

///////////////////////////////////////////////////////////////////////
// hyper
#define HYPER_UNIFORM  0
#define HYPER_FRECHET  1 
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
// direct
#define DIRECT_MAXVAR (COMMON_GAUSS + 1)
#define MAX_DIRECT_MAXVAR 30000



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
// local ce covariance operator und q-values for ce method / ce operator
#define pLOC_DIAM 0 // parameter p
#define pLOC_A 1 // always last of pLOC_*
#define pLOC_R pLOC_A // always last of pLOC_*

#define LOCAL_R 0 
#define LOCAL_MSG (LOCAL_R + 1)

#define CUTOFF_B (LOCAL_MSG + 1)
#define CUTOFF_ASQRTR (CUTOFF_B + 1)
#define CUTOFF_THEOR (CUTOFF_ASQRTR + 1) // muss immer == 4 sein

#define CUTOFF_CUBE_A (CUTOFF_THEOR + 1)
#define CUTOFF_CUBE_B (CUTOFF_CUBE_A + 1)
#define CUTOFF_CUBE_C (CUTOFF_CUBE_B + 1)
#define CUTOFF_CONSTANT (CUTOFF_CUBE_C + 1)

#define CUTOFF_MAX (CUTOFF_CONSTANT + 1)   /* size of vector q */


#define CUTOFF_THIRD_CONDITION 3



#define INTRINSIC_A0 (LOCAL_MSG + 1)
#define INTRINSIC_A2 (INTRINSIC_A0 + 1)
#define INTRINSIC_B (INTRINSIC_A2 + 1)
#define INTRINSIC_MAX (INTRINSIC_B + 1) /* size of vector q */
#define LOCAL_MAX (INTRINSIC_MAX + CUTOFF_MAX) /* sicher ist sicher */



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
#define DISTR_NROW 4
#define DISTR_NCOL 5
#define DISTR_ENV 6
#define DISTR_LAST DISTR_ENV

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


//////////////////////////////////////////////////////////////////////
// the different levels of printing

#ifndef PL_IMPORTANT
#define PL_IMPORTANT 1 // crucial messages
#define PL_BRANCHING 2 // user relevant RPgauss branching etc
#define PL_DETAILSUSER 3
#define PL_RECURSIVE 4 
#define PL_STRUCTURE 5 // see also initNerror.ERROROUTOFMETHOD
#define PL_ERRORS  6 // only those that are caught internally

#define PL_FCTN_DETAILS 7  // R
#define PL_FCTN_SUBDETAILS 8

#define PL_COV_STRUCTURE 7 // C
#define PL_DIRECT_SEQU 8
#define PL_DETAILS 9
#define PL_SUBDETAILS 10
#endif


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





typedef double *coord_type[MAXSIMUDIM];
typedef struct location_type {
  int 
    timespacedim,      /* total dimension, including time and space */
  //length[MAXSIMUDIM], /* grid==true: number of pixels of each direction 
  //		       grid==false: length[0]=number of spatial points in total
  //		       	     */
    spatialdim, 
    xdimOZ, // physical xdim, without time
    len; // nur beim ersten Element sinnvoll. Aber der 
  // Einfachheit nicht nochmals einen strukt drumrumgemacht, da sowieso
  // meist == 1
 
  int  lx, ly, /*
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
    delete_x,delete_y,  /* is loc->x only a reference? -- used in MLE */
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
			  TypeMproj, // including scale: values need not be 1 
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
  double **lpx;
  int type, // LISTOF + REALSXP
    len, // identical to nrow of calling structure
    *ncol, // in case of random parameter, they might be given
    *nrow;   // although p[.] is NULL; latter set in CHECK
} listoftype;

typedef double *naptr_type[MAX_NA];
typedef cov_model *covptr_type[MAX_NA];
//typedef int *elptr_type[MAX_MLE_ELMNTS];
// typedef double *internal_type[MAXINTERNALPARAM];
extern naptr_type MEMORY[MODEL_MAX+1];
extern int MEM_NAS[MODEL_MAX+1];

typedef double (*spectral_density)(double *, cov_model *cov); 
typedef struct spectral_storage {
  double phistep2d, phi2d, prop_factor;
  bool grid; 
} spectral_storage;

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
typedef void (*tworeturns_fct)(cov_model*, double*, double*); /* cov, result */ 
typedef void (*return_covmat)(cov_model*, double*); /* cov, result */ 
typedef char (*ext_bool_ret_fct)(cov_model*); /* covt */ 
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
typedef sortsofparam  (*sortof_fct)(int k, int row, int col);


// here all the different method for simulating a RF of a specified covariance
// function are stored

#define ONEARGUMENT_NAME 'k'
#define DEFAULT_SUBNAME 'C'
typedef struct cov_fct {  
  char name[MAXCHAR], nick[MAXCHAR];
  int kappas, // number of parameters  
    minsub, maxsub, variants; // ==0
  domain_type domain;
  isotropy_type Isotropy[MAXVARIANTS];
  int vdim, maxdim, maxmoments, Monotone;  //
  char kappanames[MAXPARAM][PARAMMAXCHAR],
    subnames[MAXSUB][PARAMMAXCHAR];
 
  SEXPTYPE kappatype[MAXPARAM]; // REALSXP, VECSXP, etc, what is expected
  //                                  from user
  Types kappaParamType[MAXPARAM];//
  //                      RandomType : parameter might be random
  //                      ShapeType : parameter must be fixed
  //                      NN2 : special purpose type
  size_fct kappasize;  // function that gives the size of the argument;
  //                     matrix or vector or scale 
  //                     (all represented as respecitve matrices)
  //
  sortsofparam sortof_tab[MAXPARAM]; // the constants used by the
  //                                     standard paramtype
  sortof_fct sortof; // particularly for MLE, function returning
  //                                   VARPARAM : ANYPARAM,
  //                          but also FORBIDDENPARAM, 

  rangefct range;
  checkfct check;
  int implemented[Forbidden], Specific;
  ext_bool finiterange;
  ptwise_type ptwise_definite; 

  int internal;
  bool 
    subintern[MAXSUB];   // do any subnames match exactly a parameter name?
  
  pref_shorttype pref;

  covfct cov, D, D2, D3, D4, tbm2, inverse, nabla, hess, random, logD; // Vtlgen 
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
 
  hyper_pp_fct hyperplane;      // hyperplane tessellations         

  bool primitive;
  Types Typi[MAXVARIANTS];
  type_fct TypeFct;

  int density;
  
  double Taylor[MAXTAYLOR][TaylorPow + 1], 
    Tail[MAXTAYLOR][TaylorExpPow + 1]; 
  int TaylorN, TailN;

  return_covmat covmatrix;
  tworeturns_fct inversecovmatrix;
  return_fct covariance, variogram, pseudovariogram;
  ext_bool_ret_fct is_covariance, is_covmatrix, is_inversecovmatrix, 
    is_variogram, is_pseudovariogram;
} cov_fct;

extern int currentNrCov, currentRegister;
typedef struct cov_fct cov_fct;
extern cov_fct *CovList;
extern int GENERALISEDCAUCHY, STABLE, BROWNIAN, CAUCHY,  GENNSST_INTERN,
  GAUSS, NUGGET, PLUS, MLEPLUS, MIXEDEFFECT, BALL, ECF, MULT, CONST, BIND,
  DISTRIBUTION,  DETERM_DISTR, GAUSS_DISTR, SETPARAM, SET_DISTR, PROD,
  COVMATRIX, RFGET, COVFCTN,LIKELIHOOD_CALL, LINEARPART_CALL,
  PREDICT_CALL,
  DOLLAR, POWER_DOLLAR,  MLE_ENDCOV,  SPLIT, SCATTER, MCMC_PGS, MCMC,
  ISO2ISO,SP2SP, SP2ISO, S2ISO, S2SP, S2S, SId, E2E, E2EIso, 
  E2Sph, E2SphIso, Sph2Sph, Sph2SphIso, LASTGATTER,
  FIRST_PLANE, LAST_PLANE, EARTHKM2CART, EARTHMILES2CART,
  EARTHKM2GNOMONIC, EARTHMILES2GNOMONIC,
  EARTHKM2ORTHOGRAPHIC, EARTHMILES2ORTHOGRAPHIC, 
  NATSC_USER, NATSC_INTERN, TREND,
  LOC, SET_DISTR, SCALESPHERICAL, TRUNCSUPPORT, BROWNRESNICK, 
  STROKORB_MONO, STROKORB_BALL_INNER, POLYGON, RECTANGULAR,
  IDCOORD, MULT_INVERSE,
  SHAPESTP, SHAPEAVE, SPHERICAL, UNIF, MPPPLUS, PTS_GIVEN_SHAPE,
  STATIONARY_SHAPE, STANDARD_SHAPE, TRAFO, COVARIATE,
// DENSITY, 
  VARIOGRAM_CALL, SIMULATE,
  BRORIGINAL_USER, BRMIXED_USER, BRSHIFTED_USER,
  ARCSQRT_DISTR,SHAPEPOW, POW, SCATTER, GNEITING_INTERN,
  
  BRORIGINAL_INTERN, BRMIXED_INTERN, BRSHIFTED_INTERN,
  MISSING_COV, NULL_MODEL, TBM_OP, USER, 
  DOLLAR_PROC, PLUS_PROC,
  MPPPLUS_PROC, MULT_PROC, // TREND_PROC,
  BINARYPROC, BROWNRESNICKPROC, GAUSSPROC, POISSONPROC,
  SCHLATHERPROC, SMITHPROC, CHI2PROC, EXTREMALTPROC, TRENDEVAL,
  NUGGET_USER,  NUGGET_INTERN,
  CIRCEMBED, CUTOFF, STEIN, SPECTRAL_PROC_USER, SPECTRAL_PROC_INTERN, 
  DIRECT, SEQUENTIAL, 
  AVERAGE_USER, AVERAGE_INTERN, HYPERPLANE_USER, HYPERPLANE_INTERN,
  RANDOMCOIN_USER, CE_CUTOFFPROC_USER, CE_CUTOFFPROC_INTERN, 
  CE_INTRINPROC_USER, CE_INTRINPROC_INTERN, TBM_PROC_USER, TBM_PROC_INTERN, 
  SPECIFIC, SELECTNR,
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
		  int internal, int vdim, int maxdim, ext_bool finiterange,
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
void add_sortof(sortof_fct sortof);
void change_sortof(int i, sortsofparam sort);
void change_typeof(int i, Types type);
bool ParamAllowsRandomness(cov_model *cov, int i);
bool ParamAllowsShape(cov_model *cov, int i);
bool ParamIsTrend(cov_model *cov, int i);
sortsofparam SortOf(cov_model *cov,int k, int row, int col);


void addReturns(return_fct Covariance, ext_bool_ret_fct isCovariance, 
		return_covmat CovMatrix, ext_bool_ret_fct isCovMatrix,
		tworeturns_fct InverseCovMatrix,
		ext_bool_ret_fct isInverseCovMatrix,
		return_fct Variogram, ext_bool_ret_fct isVariogram,
		return_fct PseudoVariogram, ext_bool_ret_fct isPseudoVariogram);
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

extern int PL, NS;
extern double GENERAL_PRECISION;

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
  double *G;
} direct_storage;


// see sequential.cc
typedef struct sequ_storage {
  int back, totpnts, spatialpnts, ntime, initial;
  double *U11, *U22, *MuT,  *G,   *Cov21, *Inv22;
  double *res0;
} sequ_storage;


// nugget
typedef struct nugget_storage {
    //double sqrtnugget;
  bool simple, simugrid;
  int *pos, reduced_dim[MAXNUGGETDIM], prod_dim[MAXNUGGETDIM + 1];
  double *red_field;
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
  // fuer shape functionen und zur Berechnung der Covariance/Variogram
  bool flat, estimated_zhou_c, logmean; 
  double old_zhou, // for mcmc only
  totalmass, // (inverser) Normierungsfaktor, um Raeumliche Funktion
  //            zu einer Dichte zu machen
    currentthreshold, log_density, globalmin, intensity, alpha;
  int 
  rowscols, /// just for control
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


#define MAX_LIN_COMP (MAXSUB * MAXSUB)
#define model_undefined -1
#define model_morethan1 -2
typedef struct likelihood_info {
  int varmodel, NAs, nas[MAX_LIN_COMP],
    effect[MAX_LIN_COMP];
  cov_model *Var; // ja nicht free, da nur pointer
  double *Matrix, 
    *pt_variance; // ja nicht free, da nur pointer
  bool trans_inv, isotropic, globalvariance;
  int newxdim, neffect;
  NAname_type names;
} likelihood_info;

typedef struct likelihood_storage {
  listoftype *datasets;  
  double **X, **YhatWithoutNA, *XCY, *XtX, *XitXi, *C, *CinvXY, *matrix,
    *betavec, **where_fixed, *sumY, *work, *Cwork,  *Xwork;
  int sets, fixedtrends, dettrends, random, max_total_data, *data_nas, maxbeta,
    betas[MAX_LIN_COMP + 1], nas[MAX_LIN_COMP], nas_det[MAX_LIN_COMP], 
    nas_fixed[MAX_LIN_COMP], nas_random[MAX_LIN_COMP], nas_boxcox,
    nas_boxcox_mu; 
  bool dettrend_has_nas, fixedtrend_has_nas, random_has_nas, 
    data_has_nas,
    betas_separate, ignore_trend;
  char *(betanames[MAX_LIN_COMP]);
  cov_model *cov_fixed[MAX_LIN_COMP], *cov_det[MAX_LIN_COMP], 
    *cov_random[MAX_LIN_COMP];
  likelihood_info info;
} likelihood_storage;



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
  double *z, *y, *z2, *y2, *save_aniso, *inv_aniso, *var;
  int *cumsum, *nx, *total, *len, proj;
  isotropy_type isoown;
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
  double *a, *b, *c, *d;
  location_type **loc; // in trafoproc, for instance
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
  bool check, dosimulate,
    prodproc_random; // used in biWM, BiGneiting
  spec_properties spec;       // used in init
  spectral_storage Sspectral; // used in do
} gen_storage;



typedef struct covariate_storage { // cov_storage, do_storage, init_storage
  // wird in Struct initialisiert, sofern INIT aufgerufen wird;
  // abgespeichert wird es immer im aktuell aufrufenden (fuer next oder key)
  location_type **loc;
  int matrix_err;
} covariate_storage;




///////////////////////////////////////////////////////////////////////
// POSITIONS OF X-VALUES
///////////////////////////////////////////////////////////////////////
#define XSTART 0
#define XSTEP 1
#define XLENGTH 2

void TransformLoc(cov_model *cov, bool timesep, int gridexpand, 
		      bool involvedollar);
void TransformLocReduce(cov_model *cov, bool timesep, int gridexpand, 
			bool involvedollar);
int TransformLoc(cov_model *cov, double **xx, bool involvedollar);
int TransformLoc(cov_model *cov, double **xx, double **yy,
		     bool involvedollar);

///////////////////////////////////////////////////////////////////////
// UNAUFGERAEUMT:
///////////////////////////////////////////////////////////////////////

//EXTERN void F77_NAME(ave1f)(double *x, double *u,
//			    double *h, double *Mz,
//			    double *r, int *d, double *res);

extern char NEWMSG[LENERRMSG];

#define WARNING1(X, Y) {sprintf(MSG, X, Y); warning(MSG); }
#define AERR(X) {ERRLINE; PRINTF("AERR: "); errorMSG(X, MSG); 	\
    if (PL<PL_ERRORS) PRINTF("%s%s\n", ERROR_LOC, MSG); assert(false);}
#define MERR(X) {LPRINT("error: ");				\
    errorMSG(X, MSG);					\
    if (PL<PL_ERRORS) PRINTF("%s%s\n", ERROR_LOC, MSG);}
#define XERR(X) {/* UERR; */ errorMSG(X, MSG); ERR(MSG);}
#define PERR(X) {sprintf(MSG, "'%s': %s", param_name, X); ERR(MSG);}
#define PERR1(X,Y) {sprintf(MSG, "'%s': %s", param_name, X); sprintf(MSG2, MSG, Y); ERR(MSG2);}
#define QERR(X) {sprintf(ERRORSTRING, "'%s' : %s", param_name, X); DEBUGINFOERR; return ERRORM;}
#define QERRX(ERR, X) {errorMSG(ERR, MSG); sprintf(ERRORSTRING, "'%s' : %s (%s)", param_name, X, MSG); DEBUGINFOERR; return ERRORM;}
#define QERRC(NR,X) {sprintf(ERRORSTRING, "%s '%s': %s", ERROR_LOC, CovList[cov->nr].kappanames[NR], X); DEBUGINFOERR; return ERRORM;}
#define QERRC1(NR,X,Y) {sprintf(MSG, "%s '%s': %s", ERROR_LOC, CovList[cov->nr].kappanames[NR], X); sprintf(ERRORSTRING, MSG, KNAME(Y)); DEBUGINFOERR; return ERRORM;}
#define QERRC2(NR,X,Y,Z) {sprintf(MSG, "%s '%s': %s", ERROR_LOC, CovList[cov->nr].kappanames[NR], X); sprintf(ERRORSTRING, MSG, KNAME(Y), KNAME(Z)); DEBUGINFOERR; return ERRORM;}

#define NotProgrammedYet(X) {						\
    { if (strcmp(X, "") == 0)						\
      sprintf(BUG_MSG,							\
	      "function '%s' (file '%s', line %d) not programmed yet.", \
	      __FUNCTION__, __FILE__, __LINE__);			\
    else								\
      sprintf(BUG_MSG, "'%s' in '%s' (file '%s', line %d) not programmed yet.",\
	      X, __FUNCTION__, __FILE__, __LINE__);			\
      RFERROR(BUG_MSG);							\
    }\
    }									


void contact(const char *msg);

void ErrCov(double *x, cov_model *cov, double *v);
void ErrCovNonstat(double *x, double *y, cov_model *cov, double *v);
void GetNaturalScaling(cov_model *cov, double *natscale);
void analyse_matrix(double *aniso, int row, int col,
		    bool *diag, bool *quasidiag, int *idx,
		    bool *semiseparatelast, bool *separatelast);


listoftype *LIST_CREATE(int len, int type);
void LIST_DELETE(listoftype **x);
void listcpy(listoftype **To, listoftype *p, bool force_allocating);
void paramcpy(cov_model *current, cov_model *cov, bool freeing,
	      bool allocating, bool copy_lists, bool recursive, bool copy_mpp);
int covCpy(cov_model **localcov, cov_model *cov);
int covCpy(cov_model **localcov, cov_model *cov, bool copy_lists);
int covCpy(cov_model **localcov, bool sub, cov_model *cov,
	   location_type **prevloc);
int covCpy(cov_model **localcov, cov_model *cov,
	   double *x, double *T, int spatialdim, int xdim, long lx, bool Time, 
	   bool grid, bool distances);
int covCpyWithoutRandomParam(cov_model **localcov, cov_model *cov);
int covCpy(cov_model **localcov, bool sub, cov_model *cov, // err
	   location_type **prevloc, location_type **ownloc,
	   bool copy_lists,  bool copy_randomparam, bool allowCopyingInterface);
void Ssetcpy(cov_model *localcov, cov_model *remotecov, cov_model *cov,
	     cov_model *rmt);

int newmodel_covCpy(cov_model **localcov, int model, cov_model *cov);
int newmodel_covCpy(cov_model **localcov, int model, cov_model *cov,
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
int setgrid(coord_type xgr, double *x, long lx, int spatialdim);
int partial_loc_set(location_type *loc, double *x, double *y,
		    long lx, long ly, bool dist, int xdim, double *T, 
		    bool grid, bool cpy);
//void partial_loc_set(cov_model *cov, double *x, double *y,
//		     long lx, bool dist, int xdim, double *T);
int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdim, long lx, bool Time, bool grid,
	    bool distances, location_type **loc);
int loc_set(double *x, double *y, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdim,
	    long lx, long ly, bool Time, bool grid,
	    bool distances,
	    location_type **Loc);

int loc_set(double *x, double *y, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ,
	    long lx, long ly, bool Time, bool grid,
	    bool distances,
	    cov_model *cov);

int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ, /* original ! */
	    long lx, bool Time, bool grid,
	    bool distances,
	    cov_model *cov);
location_type ** loc_set(SEXP xlist, bool dist_ok);
int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ, /* original ! */
	    long lx, bool Time, bool grid,
	    bool distances, int n,
	    location_type ***Loc);
void SetLoc2NewLoc(cov_model *cov, location_type **loc);

//int loc_set(cov_model *cov, long totalpoints);
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
			SEXP xlist,
			cov_model **Cov);
void CheckModel(SEXP model, SEXP tsdim, SEXP xdim, SEXP domain,
		cov_model **Cov,
		int maxdim);


matrix_type Type(double *m, int nrow, int ncol);
double GetDiameter(location_type *loc);
double GetDiameter(location_type *loc, double *min, double *max,double *center);
double GetDiameter(location_type *loc, double *min, double *max,
		   double *center, bool docaniso);
 
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
int internal_DoGauss(cov_model  *key, int nn, double *orig_res);

void strcopyN(char *dest, const char *src, int n);
//void expandgrid(coord_type xgr, double **xx, int nrow);
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
int check2Xnotrafo(cov_model *cov, int tsdim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
		   int vdim, int role);
  int check2X(cov_model *cov, int tsdim, int tsxdim, Types type, 
	    domain_type domprev, isotropy_type isoprev, 
	    int vdim, int role);
int check2X(cov_model *cov, int tsdim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
	    int *vdim, int role);
int check2X(cov_model *cov, int tsdim, int tsxdim, Types type, 
	    domain_type domprev, isotropy_type isoprev,
	    int vdim0, int vdim1, int role, bool coordinate_trafo);


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


void AxResType(double *A, double *x, int nrow, int ncol, double *y);

SEXP Param(int nr, void* p, int nrow, int ncol, SEXPTYPE type, bool drop, long* mem);
int is_positive_definite(double *C, int dim);

// typedef long long POINTER;
typedef long int POINTER;

extern int PrInL;				

void SetTemporarilyGlobal(globalparam *gp);
void ReSetGlobal();
				
int change_coord_system(isotropy_type callingprev, isotropy_type isoprev,
			int tsdim, int xdimprev,
			int *nr, isotropy_type  *newisoprev,
			int *newtsdim, int *newxdim, bool time);
int SetGatter(domain_type domprev, domain_type statnext, 
	      isotropy_type isoprev, isotropy_type isonext, // bool first, 
	      int *nr, int *delflag);
				
int cpy_add_gatter(cov_model *cov, cov_model **sm);

int binomialcoeff(int n, int k);
void poly_basis(cov_model *cov, gen_storage *s);
double evalpoly(double *x, int *powmatrix, double *polycoeff, int basislen,
	      int dim);

void MultiDimRange(int set, cov_model *cov, double *natscale);

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
  location_type **prevloc, **ownloc;
  cov_model *key; // this one should be deleted
  // double *caniso, cscale, cvar; 
  //  int *cproj;
  //  int xdimout;
  // matrix_type type;   // von oben
  
  bool initialised, origrf, checked;
  ext_bool loggiven, fieldreturn;
  double
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
  solve_storage *Ssolve;
  biwm_storage *Sbiwm;
  inv_storage *Sinv;
  scatter_storage *Sscatter;
  mcmc_storage *Smcmc;
  gen_storage *Sgen; // only once for the whole model tree
  likelihood_storage *Slikelihood;
  covariate_storage *Scovariate; // only once for the whole model tree
 
  //check in getNset:  COV_DELETE_WITHOUTSUB, COV_ALWAYS_NULL; add *_DEL, *_NULL
  
  //select_storage *Sselect;
 } cov_model;


void LOC_DELETE(location_type ***Loc);
void LOC_NULL(location_type **loc);
void LOCLIST_CREATE(cov_model *cov, int n);
location_type **LOCLIST_CREATE(int n);
void COV_DELETE_(cov_model **cov);
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
void likelihood_NULL(likelihood_storage *x);
void likelihood_DELETE(likelihood_storage **S);
//likelihood_storage * likelihood_CREATE(int n);
void likelihood_info_NULL(likelihood_info *x);
void likelihood_info_DELETE(likelihood_info *x);
void covariate_NULL(covariate_storage *x);
void covariate_DELETE(covariate_storage **S);


void trend_DELETE(trend_storage ** S);
void trend_NULL(trend_storage* x);
void mixed_DELETE(mixed_storage **S);
void mixed_NULL(mixed_storage* x);
void gen_NULL(gen_storage *x);
void gen_DELETE(gen_storage **S);


//void SELECT_NULL(select_storage *S);
//void SELECT_DELETE(select_storage **S);


///////////////////////////////////////////
#define LocLoc(loc) (loc != NULL ? loc[GLOBAL.general.set % loc[0]->len] : NULL)
#define Loc(cov) (cov->ownloc!=NULL ?LocLoc(cov->ownloc) :LocLoc(cov->prevloc))
#define PLoc(cov) (cov->ownloc != NULL ? cov->ownloc : cov->prevloc)
#define PrevLoc(cov) (cov->prevloc == NULL ? NULL : cov->prevloc[GLOBAL.general.set % cov->prevloc[0]->len])
#define OwnLoc(cov) (cov->ownloc == NULL ? NULL : cov->ownloc[GLOBAL.general.set % cov->ownloc[0]->len])
#define Getgrid(cov) (PLoc(cov) != NULL ? Loc(cov)->grid : false)
#define GetTime(cov) (PLoc(cov) != NULL ? Loc(cov)->Time : false)
#define Gettimespacedim(cov) (PLoc(cov) != NULL ? Loc(cov)->timespacedim : 0)
#define Gettotalpoints(cov) (PLoc(cov) != NULL ? Loc(cov)->totalpoints : -1)
#define GET_LOC_SETS(cov) (PLoc(cov) != NULL ? PLoc(cov)[0]->len : 0)


void crash();

int checktbm_basics(cov_model *cov, bool tbmop);

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


void AxResType(double *A, double *x, int nrow, int ncol, double *y) ;

double detU(double *C, int dim); // strictly pos def matrices
void det_UpperInv(double *C, double *det, int dim) ;
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
    assert((Cov)->initialised);					  \
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


void PrintLoc(int level, location_type *loc, bool own);

void pmi(cov_model *cov, char all, int level, int maxlevel); 
void pmi(cov_model *cov, int maxlevel);
void pmi(cov_model *cov);
void pci(int nr);
void pci();
#define APMI(cov) {\
    PRINTF("\n'%s', line %d", __FILE__, __LINE__); \
    pmi(cov); assert(false);}
#define AAPMI(cov, txt) { pmi(cov, txt); assert(false); }

#define PMI(COV) {						\
    PRINTF("\n(PMI '%s', line %d)", __FILE__, __LINE__);	\
    pmi(COV);							\
  }

#define PMIL(COV, ML) {						\
  PRINTF("\n(PMI '%s', line %d)", __FILE__, __LINE__);		\
  pmi(COV, ML);							\
}



#define PCI							\
    PRINTF("\n(PMI '%s', line %d)", __FILE__, __LINE__);	\
    pci
      




#define PLE PRINTF("\n(PLE '%s', line %d)", __FILE__, __LINE__); ple_
void ple_(cov_model *cov);
void ple_(char *name);



void memory_copy(void *dest, void *src, int bytes);
//int solvePosDef(double *M, int vdimnf);

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
bool isRObject(int type);


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


#define EXT_NEW_COV_STORAGE(cov, new) {				\
  if ((cov)->S##new != NULL) {					\
    Ext_##new##_DELETE(&((cov)->S##new));		    	\
      assert((cov)->S##new == NULL);			       	\
  }								\
  if ((cov)->S##new == NULL) {					\
    (cov)->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    Ext_##new##_NULL((cov)->S##new);					\
    if ((cov)->S##new == NULL) BUG;				\
  }} 							

#define EXT_NEW_STORAGE(new)	\
  EXT_NEW_COV_STORAGE(cov, new)

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


#define CONDCOV_NEW_STORAGE(cov, new, WHAT) {				\
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

#define COND_NEW_STORAGE(new, WHAT) CONDCOV_NEW_STORAGE(cov, new, WHAT)


#define ALLOCCOV_NEW(cov, Snew, Z, SIZE, WHAT)				\
  assert(cov->Snew != NULL);						\
  double *Z = cov->Snew->WHAT;						\
  if (Z == NULL)							\
    Z = cov->Snew->WHAT = (double*) MALLOC(sizeof(double) * (SIZE)) 
  

#define ALLOC_NEW(Snew, Z, SIZE, WHAT) ALLOCCOV_NEW(cov, Snew, Z, SIZE, WHAT)
  

#define ALLOC_NEWINT(Snew, Z, SIZE, WHAT)				\
  assert(cov->Snew != NULL);						\
  int *Z = cov->Snew->WHAT;						\
  if (Z == NULL)							\
    Z = cov->Snew->WHAT = (int*) MALLOC(sizeof(int) * (SIZE)) 


#define SOLVE_STORAGE EXT_NEW_STORAGE(solve)


#define EXTRA_STORAGE COND_NEW_STORAGE(extra, a)
#define ALLOC_EXTRA(Z, SIZE) ALLOC_NEW(Sextra, Z, SIZE, a) 
#define ALLOC_EXTRA2(Z, SIZE) ALLOC_NEW(Sextra, Z, SIZE, b)     
#define ALLOC_EXTRA3(Z, SIZE) ALLOC_NEW(Sextra, Z, SIZE, c)     
#define ALLOC_EXTRA4(Z, SIZE) ALLOC_NEW(Sextra, Z, SIZE, d)     


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
    double colour;
} cell_type;
typedef int (*avl_comparison_func) (cell_type *, cell_type *, int *);
typedef void (*avl_node_func) (cell_type*, int *);
typedef cell_type *(*avl_copy_func) (cell_type*, int *);

bool CallingSet(cov_model *cov);
void includeStandardMath();

#define I_COL_NA -1

void setptwise(ptwise_type pt);


#ifdef SCHLATHERS_MACHINE
#define STOPAFTER(COND, DO) 				\
  static bool __stopafter__ = false;			\
  if (__stopafter__ && !(COND)) { DO; BUG;}		\
  __stopafter__  = COND;				
#else 
#define STOPAFTER(COND, DO)
#endif
  							
    
bool TrafoOK(cov_model *cov);

int SetAndGetModelInfo(cov_model *cov, int shortlen, 
		       int allowforintegerNA, bool excludetrend, // IN
		       int newxdim,  //IN
		       int globvar,
		       likelihood_info *info); // OUT		

void boxcox_inverse(double boxcox[], int vdim, double *res, int pts, 
		    int repet);	

#define SAVE_GAUSS_TRAFO    

#define BOXCOX_TRAFO boxcox_trafo(boxcox, res, Gettotalpoints(cov) * cov->vdim[0]); 

#define BOXCOX_INVERSE_PARAM(GAUSS_BOXCOX)				\
  boxcox_inverse(P(GAUSS_BOXCOX), cov->vdim[0], res, Gettotalpoints(cov), 1)

#define BOXCOX_INVERSE BOXCOX_INVERSE_PARAM(GAUSS_BOXCOX)


int kappaBoxCoxParam(cov_model *cov, int BC);


#define  GAUSS_COMMON_RANGE_PARAM(GAUSS_BOXCOX)	\
  range->min[GAUSS_BOXCOX] = RF_NEGINF;		\
  range->max[GAUSS_BOXCOX] = RF_INF;		\
  range->pmin[GAUSS_BOXCOX] = 0;		\
  range->pmax[GAUSS_BOXCOX] = 2;		\
  range->openmin[GAUSS_BOXCOX] = false;		\
  range->openmax[GAUSS_BOXCOX] = false

#define  GAUSS_COMMON_RANGE GAUSS_COMMON_RANGE_PARAM(GAUSS_BOXCOX)

#define MAX_LEN_EXAMPLES 4




#endif

