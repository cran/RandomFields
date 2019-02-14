
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part

 Copyright (C) 2001 -- 2017 Martin Schlather

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



// to do: V3.1+: Randverteilungen der RPs in cov, D etc implementieren;
//        do -> random part; init wird innerhalb von do aufgerufen,
//        falls nicht initialisiert oder random submodels??
// to do: MLE: random parameters einsammeln


#include <stdio.h>  
#include <string.h>
#include "RF.h"
#include "operator.h"
#include "startGetNset.h"
#include "primitive.h"
#include "init.h"
#include "questions.h"
#include "shape.h"
#include "primitive.others.h"


KEY_type *PIDKEY[PIDMODULUS];
int parentpid;


defn *DefList=NULL;
int currentNrCov=UNSET,
  gaussmethod[Forbidden+1];
double ONE = 1;
char CovNames[MAXNRCOVFCTS][MAXCHAR], CovNickNames[MAXNRCOVFCTS][MAXCHAR],
  *FREEVARIABLE= (char*) "...";


//                        CE         CO         CI        TBM       Sp
//                        di         sq       trend        av       n
//                        mpp         Hy       spf        any     forbidden
pref_type PREF_ALL = {PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, 
		      PREF_BEST, PREF_BEST, PREF_NONE, PREF_BEST, PREF_BEST,
		      PREF_BEST, PREF_BEST, PREF_NONE, // specific
		                                       PREF_BEST, PREF_BEST},
  PREF_NOTHING = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST, // nothing
		                                              PREF_NONE},
  PREF_TREND =  {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_BEST, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST, // nothing
		                                              PREF_NONE},
  PREF_AUX = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE},
  PREF_MATHDEF = {PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST,
		  PREF_BEST, PREF_BEST, PREF_BEST, PREF_NONE, PREF_NONE,
		  PREF_BEST, PREF_BEST, PREF_NONE, PREF_BEST, PREF_NONE};
 



const char *METHOD_NAMES[Forbidden+1]={"circulant", //0
				      "cutoff",
				      "intrinsic",
				      "tbm", 
				      "spectral", //4
				      "direct",
				      "sequential",
				      "trend",
				      "average",
				      "nugget", //9
				      "coins",
				      "hyperplane",
				      "specific",
				      "any method", // nothing
				      "forbidden"},
  *CAT_TYPE_NAMES[OtherType + 1] = {
  // TcfType, PosDefType, VariogramType, NegDefType,
  // PointShapeType, ShapeType, TrendType, RandomOrShape, Manifold,
  // ProcessType, GaussMethodType, NormedProcessType, BrMethodType,
  // SmithType, SchlatherType, PoissonType, PoissonGaussType,
  // RandomType, InterfaceType, MathDefType, OtherType
  // not including BadType, ...
    "RM", "RM", "RM", "RM",
    "RM", "RM", "RM", "RM", "RM", 
    "RP", "RP", "RP", "RP",
    "RP", "RP", "RP", "RP",
    "RR",  "RF", "R.", "RO"},
  
  *REG_NAMES[MODEL_MAX+1] = {"reg0", "reg1", "reg2", "reg3", "reg4", 
			    "reg5", "reg6", "reg7", "reg8", "reg9",
			     "user", "cov", "covmatrix", "variogram",
			     "pseudovariogram",
			     "function", "distribution", "calculation",
			     "unused", "unused",
			     "auxiliary", "intern", "split", "gui", "mle",
			     "mlesplit", "mletrend", "mlebounds",
		 	    "kriging", "conditional", "error model"},
  *POSITIVITY_NAMES[(int) pt_mismatch + 1] = 
			      {"pt-wise positive definite", 
			       "pt-wise indefinite", 
			       "pt-wise negative definite", "pt-wise zero", 
			       "pt-wise param dependent", 
			       "pt-wise submodel dependent", 
			       "pt-wise undefined", "pt-wise unknown", 
			       "pt-wise mismatch"},
  **LIST_OF_NAMES[nNamesOfNames] = // see also NAME_OF_NAMES in AutoRandomF*.cc
				 {EQ_NAMES, ISO_NAMES, DOMAIN_NAMES,
				  TYPE_NAMES, MONOTONE_NAMES, MODE_NAMES,
				  OUTPUTMODE_NAMES, REPORTCOORD_NAMES,
				  UNITS_NAMES, COORD_SYS_NAMES,
				  CARTESIAN_SYS_NAMES, TYPEOF_PARAM_NAMES};

//  int SYS_TO_ISO[nr_coord_sys_proj] =  {GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ};



char
  STANDARDPARAM[MAXPARAM][MAXCHAR],
  STANDARDSUB[MAXSUB][MAXCHAR];

model **KEY() { return KEYT()->KEY; }
KEY_type *KEYT() {  
  int mypid;
  Ext_pid(&mypid);
  bool parallel = mypid != parentpid;

  if (parallel && GLOBAL.internal.warn_parallel) {
    GLOBAL.internal.warn_parallel = false;
    PRINTF("Do not forget to run 'RFoptions(storing=FALSE)' after each call of a parallel command (e.g. from packages 'parallel') that calls a function in 'RandomFields'. (OMP within RandomFields is not affected.) This message appears only once per session."); //
  }

  KEY_type *p = PIDKEY[mypid % PIDMODULUS];
  if (p == NULL) {
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    assert(neu != NULL);
    assert(neu->zerox == NULL);
    PIDKEY[mypid % PIDMODULUS] = neu;
    neu->visitingpid = mypid;    
    if (PIDKEY[mypid % PIDMODULUS] != neu) { // another process had the
      //                                        same idea
      FREE(neu);
      return KEYT(); // ... and try again
    }
    neu->pid = mypid;
    neu->visitingpid = 0;
    neu->ok = true;
    if (PIDKEY[mypid % PIDMODULUS] != neu) BUG;
    KEY_type_NULL(neu);
     return neu;
  }
  while (p->pid != mypid && p->next != NULL) p = p->next;
  if (p->pid != mypid) {
    if (!p->ok || p->visitingpid != 0) {
      if (PL >= PL_ERRORS) {
	PRINTF("pid collision %d %d\n",  p->ok, p->visitingpid);
      }
      //    BUG;
     return KEYT();
    }
    p->visitingpid = mypid;
    p->ok = false;
    if (p->visitingpid != mypid || p->ok) return KEYT();
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    neu->currentRegister = UNSET;
    neu->pid = mypid;
    if (!p->ok && p->visitingpid == mypid) {
      p->next = neu;
      p->visitingpid = 0;
      p->ok = true;
      
       return neu;
    }
    FREE(neu);
    p->visitingpid = 0;
    p->ok = true;
    KEY_type_NULL(neu);
    return KEYT();
  }
  MEMCOPY(&(p->global), &GLOBAL, sizeof(globalparam));
  p->error_causing_cov = NULL;
  return p;
}

int currentRegister() {
  KEY_type *p = KEYT();
  if (p == NULL) BUG;
  return p->currentRegister;
}


void set_currentRegister(int cR) {
  KEY_type *p = KEYT();
  if (p == NULL) BUG;
  p->currentRegister = cR;
}

int zaehler = 0;
double *ZERO(int dim, KEY_type *KT) {
  //   printf("zero enter %lu %lu\n", KT, KT->zerox);
  assert(KT != NULL);
  if (dim > KT->nzero) {
    FREE(KT->zerox); 
    KT->nzero = dim;
    KT->zerox = (double*) CALLOC(KT->nzero, sizeof(double));
  }
  //  printf("%lu %lu %d %10g\n", KT, KT->zerox, KT->nzero,  KT->zerox[0]);
  // if (zaehler++) crash();
  return KT->zerox;
}

double *ZERO(model *cov) {
  assert(cov != NULL && cov->base != NULL);
  return ZERO(PREVTOTALXDIM + 1, cov->base); // i_ncol
}


bool TrafoOK(model *cov, const char *file, int line) {// check other things, too, besides gatter ?
  //  PMI0(cov);
  //   if (!PREV_INITIALISED) crash();
  //  assert(PREV_INITIALISED);
    
   bool ok = ((FIRSTGATTER <= GATTERNR && GATTERNR <= LASTGATTER)
	     &&
	     ((TRAFONR >= FIRST_TRAFO && TRAFONR <= LAST_TRAFO) ||
	      TRAFONR == UNSET))
    && cov->checked;
  if (!false && !ok) {
    PMI0(cov->calling); //
    PMI0(cov);//
    PRINTF("%.50s: not: %d < %d <= %d UND (%d <= %d <= %d oder %d == %d) UND checked=%d in %.50s, line %d\n", 
	   NAME(cov), FIRSTGATTER, GATTERNR, LASTGATTER, 
	   //
	    FIRST_TRAFO, TRAFONR, LAST_TRAFO, TRAFONR, UNSET, 
	   //
	   cov->checked, file, line);
    // crash();
    // BUG;
  }
  return ok;
}


int
  PLUS,
  USER,    
  DOLLAR, LASTDOLLAR, 
  ISO2ISO0, SP2SP0, SP2ISO0, S2ISO0, S2SP0, S2S0, SId0, E2E0, E2EIso0, 
  E2Sph0, E2SphIso0, Sph2Sph0, Sph2SphIso0,  FIRSTGATTER0, LASTGATTER0,
  FIRST_PLANE, LAST_PLANE, EARTHKM2CART, EARTHMILES2CART,
  EARTHKM2GNOMONIC, EARTHMILES2GNOMONIC,
  EARTHKM2ORTHOGRAPHIC, EARTHMILES2ORTHOGRAPHIC,  
  FIRST_TRAFO, LAST_TRAFO;


   

bool CheckListmodel(){
  assert((bool) 5 == 1);
  assert(MODEL_MAX == 30); // otherwise change REG_NAMES
  assert(OtherType == 20 && LASTTYPE == 30); // otherwise change TYPE_NAMES, 
  //                                           CAT_TYPE_NAMES[OtherType + 1]
  //                                           Type-Consistency in question.cc
  assert(LAST_ISO == 19); // otherwise change ISO_NAMES
  //  assert(MAXMPPDIM <= MAXSIMUDIM); // ZERO
  // assert(CUTOFF_THEOR == 4);/* do not change this value as used in RFmethods.Rd */
  

  int nr, err=NOERROR;
  defn *C;
  for (nr=0; nr<currentNrCov; nr++) {     
    C = DefList + nr; // nicht gatternr    
    //  printf("%.50s\n", C->name);
    if (C->vdim != SCALAR && equalsSymmetric(C->systems[0][0].iso)
	&& C->check != checkmqam
	 ) BUG;
    
    err = 1;
    if (isManifold(SYSTYPE(C->systems[0], 0)) && C->TypeFct == NULL) 
      goto ErrorHandling;

    err = 2;
    if (isMathDef(DefList + nr)) {
      if (DefList[nr].variants > 2 && nr != CONST) goto ErrorHandling;
    }

    err = 3;
    if ((equalsParamDepI(C->systems[0][0].iso) ||
	 equalsParamDepD(C->systems[0][0].dom))
	&& (C->TypeFct == NULL || C->setDI == NULL)) goto ErrorHandling;

    err = 4;
    if (C->setDI == NULL) {
      if (isParamDepI(C) || isParamDepD(C)) goto ErrorHandling;
    } else {
      if (!(isParamDepI(C) || isParamDepD(C) || isSubModelI(C)))
	goto ErrorHandling;
    }

    err = 5;
    if (C->TypeFct != NULL) {
      Types type = SYSTYPE(C->systems[0], 0);
      if(isProcess(type) || isPointShape(type) || equalsRandom(type) ||
	 equalsInterface(type)) goto ErrorHandling;
    }

    err = 6;
    for (int k=0; k<C->kappas; k++) {
      if (C->kappanames[k][0] == ONEARGUMENT_NAME
	  && C->kappanames[k][1] >= '0'
	  && C->kappanames[k][1] <= '9') goto ErrorHandling;
    }

    err = 7;
    if ( !(!equalsUnreduced(ISO(C->systems[0], 0)) || C->variants == 1 ||
	   (C->variants == 2 && equalsUnreduced(ISO(C->systems[1], 0))) ||
	   ((C->variants >= 3) &&
	    C->Iallowed != NULL && C->Iallowed != allowedIstandard &&
	    C->Iallowed != allowedPrevModelI &&
	    equalsIsotropic(ISO(C->systems[1], 0)) &&
	    equalsEarthIsotropic(ISO(C->systems[2], 0)) &&
	    (C->variants == 3 || equalsUnreduced(ISO(C->systems[3], 0))))
	   )) goto ErrorHandling;

  }
  return true;

 ErrorHandling:
  PRINTF("\n\nFehler in CheckListcov: name = %.50s %d (err=%d)\n\n", C->nick, nr, err);
  // printf("%d\n", (C->setDI != NULL) ^ (isParamDepI(C) || isParamDepD(C)));
  return false;
}

void InitModelList() {
  assert(currentNrCov == UNSET); // otherwise something went wrong with the call
 
  for (int i=0; i<MAXPARAM; i++) SPRINTF(STANDARDPARAM[i], "k%d", i+1);
  for (int i=0; i<MAXSUB; i++) SPRINTF(STANDARDSUB[i], "u%d", i+1);
  /* ja nicht setzen !! macht riesen aerger, da RF opt ions InitModel
    nicht aufruft:
    for (i=0; i<MAX UNITS; i++) {
    STRCPY(GLOBAL.general.newunits[i], "");
    STRCPY(GLOBAL.general.curunits[i], "");
    STRCPY(GLOBAL.general.varunits[i], "");
  }
  */


   // init models
  Ext_pid(&parentpid);
  
  for (int i=0; i<PIDMODULUS; i++) PIDKEY[i] = NULL; 

  if (DefList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  DefList = (defn*) MALLOC(sizeof(defn) * (MAXNRCOVFCTS+1));
  // + 1 is necessary because of COVINFO_NULL that uses the last + 
  currentNrCov = 0;

 
  // *******************
  // **** RO-models ****
  // *******************


  FIRSTGATTER0 =  // 0 -- always first
    IncludeModel("#",  OtherType, 1, 1, 0, NULL, PREVMODEL_D, PREVMODEL_I,
  		 checkNotOK, NULL, PREF_NOTHING, true, SUBMODEL_DEP,
  		 SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  assert(FIRSTGATTER0 == FIRSTGATTER);
  addCov(stat2, D_2, DD_2, inverse2, nonstatinverse2);
  addCov(nonstat2);// 
  addlogCov(logstat2, lognonstat2, nonstat_loginverse2);
  // addCov(iso2iso, D_2, DD_2, inverse2, nonstatinverse2);
  RandomShape(INFTY, struct2, init2, do2, dorandom2, true, true, false); 

  
  ISO2ISO0 = addFurtherCov(ErrCov, ErrD, ErrD); // 1
  SP2SP0 = addFurtherCov(ErrCov, ErrD, ErrD); // 2
  SP2ISO0 = addFurtherCov(ErrCov, ErrD, ErrD); // 3
  S2ISO0 = addFurtherCov(ErrCov, ErrD, ErrD); // 4
  S2S0 = addFurtherCov(ErrCov, ErrD, ErrD);// 5
  SId0 = addFurtherCov(ErrCov, ErrD, ErrD);// 6
  S2SP0 = addFurtherCov(ErrCov, ErrD, ErrD);// 7
  E2EIso0 = addFurtherCov(ErrCov, ErrD);// 8
  E2E0 = addFurtherCov(ErrCov, ErrD);// 9
  E2SphIso0 = addFurtherCov(ErrCov, ErrD);// 10
  E2Sph0 = addFurtherCov(ErrCov, ErrD);// 11
  Sph2SphIso0 = addFurtherCov(ErrCov, ErrD);// 12
  Sph2Sph0 = addFurtherCov(ErrCov, ErrD);// 13

  LASTGATTER0 = Sph2Sph0;
  assert(LASTGATTER0 == LASTGATTER);
  
  FIRST_TRAFO = EARTHKM2CART= // 14
      IncludeModel(">",  OtherType, 1, 1, 0, NULL,
		   PREVMODEL_D, PREVMODEL_I, // dummy values
		   checkEarth, NULL, PREF_NOTHING, true, SUBMODEL_DEP,
		   4, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  
  addCov(EarthKM2CartStat, NULL, NULL);
  addlogCov(EarthKM2Cart);//

  
  
  EARTHMILES2CART = addFurtherCov(EarthMiles2CartStat, ErrD);// 15
  addlogCov(EarthMiles2Cart);// 

  FIRST_PLANE = EARTHKM2GNOMONIC =
    addFurtherCov(Earth2GnomonicStat, ErrD);// 16
  addlogCov(Earth2Gnomonic);// 

  EARTHMILES2GNOMONIC =  CopyModel(">", EARTHKM2GNOMONIC); // 17

  EARTHKM2ORTHOGRAPHIC = addFurtherCov(EarthKM2OrthogStat, ErrD);// 18
  addlogCov(EarthKM2Orthog);// 

  EARTHMILES2ORTHOGRAPHIC = addFurtherCov(EarthMiles2OrthogStat, ErrD);// 19
  addlogCov(EarthMiles2Orthog);// 
 
  LAST_PLANE = EARTHMILES2ORTHOGRAPHIC;
  LAST_TRAFO =  EARTHMILES2ORTHOGRAPHIC;
  assert(LAST_TRAFO == 19);

  pref_type pplus =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Tr av n mpp Hy spf any
  PLUS = 
    IncludeModel("+", ManifoldType, 1, MAXSUB, 0, NULL, SUBMODEL_D, SUBMODEL_I,
		 checkplus, NULL, pplus, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP,
		 (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  // Achtung in covmatrix_plus wird SELECT_SUBNR verwendet!
  nickname("plus");
  addCov(plusStat, Dplus, DDplus, NULL, NULL);
  addCov(plusNonStat);
  addTBM(NULL, spectralplus);
  RandomShape(0, structplus, initplus, doplus);
  addReturns(covmatrix_plus, iscovmatrix_plus);
  setptwise(pt_submodeldep);
  addTypeFct(Typeplus);
  setDI(allowedDplus, allowedIplus, NULL);

  pref_type pmal =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 4, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  MULT = IncludeModel("*", ManifoldType,  1, MAXSUB, 0, NULL,
		      SUBMODEL_D, SUBMODEL_I,
		      checkmal, NULL, pmal, false, SUBMODEL_DEP,
		      SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  nickname("mult");
  addCov(malStat, Dmal, NULL);
  addCov(malNonStat);
  addlogCov(logmalStat, logmalNonStat, NULL);
  setptwise(pt_submodeldep);
  addTypeFct(Typemal);
  setDI(allowedDplus, allowedIplus, NULL);


  pref_type pS=  {5, 0, 0,  5, 5, 5, 5, 0, 0, 5, 0, 0, 1, 5};
  //              CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  DOLLAR = IncludeModel("$",  ManifoldType, // to do: tcftype durch einen allgemeinen Type ersetzen, da auch Trend dem "$" folgen kann. Z.Z. nicht moeglich.
			1, 1, 5, kappaS, // kappadollar,
			SUBMODEL_D, SUBMODEL_I, checkS, rangeS, pS,
			false, SUBMODEL_DEP, SUBMODEL_DEP,
			(ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  // do not change Order!!
  nickname("S");
  kappanames("var", REALSXP, "scale", REALSXP, "anisoT", REALSXP,
	     "Aniso", REALSXP, "proj", INTSXP);
  change_typeof(DVAR, RandomOrShapeType);
  change_typeof(DSCALE, RandomOrShapeType);
  change_typeof(DAUSER, ShapeType);
  subnames("phi");
  addCov(Siso, DS, DDS, D3S, D4S, inverseS, nonstatinverseS); // unterscheidung nur wegen der 
  //  geschwindigkeit, also Siso ist ein sehr haeufiger Spezialfall von Sstat
  addCov(Snonstat);
  addlogCov(logSiso, NULL, nonstat_loginverseS);
  addLocal(coinitS, ieinitS);  
  addTBM(tbm2S, NULL, spectralS);
  nablahess(nablaS, hessS);
  RandomShape(INFTY, structS, initS, doS, true, true, false);
  addReturns(covmatrixS, iscovmatrixS);
  Taylor(RF_NA, RF_NA, RF_NA, RF_NA);
  TailTaylor(RF_NA, RF_NA, RF_NA, RF_NA);
  setptwise(pt_submodeldep); 
  addTypeFct(TypeS);
  setDI(allowedDS, allowedIS, NULL); //setS
  
   
  LASTDOLLAR = addFurtherCov(Sstat, DS, DDS); // 20.8.14 aus ErrCov (wieder)
  //                                        D2 gemacht
  addCov(Snonstat);
  addlogCov(logSstat, logSnonstat, NULL);

  // printf("%d\n",  currentNrCov); BUG;

  pref_type pPowS=  {5, 0, 0,  5, 5, 5, 5, 0, 0, 5, 0, 0, 1, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  POWER_DOLLAR = 
    IncludeModel("$power", ManifoldType, // to do: tcftype durch einen allgemeinen Type ersetzen, da auch Trend dem "$" folgen kann. Z.Z. nicht moeglich.
		 1, 1, 3, NULL, // kappadollar,
		 SUBMODEL_D, SUBMODEL_I, checkPowS, rangePowS, pPowS,
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_SUB_DEP);
  // do not change Order!!
  nickname("Spower");
  kappanames("var", REALSXP, "scale", REALSXP, "pow", REALSXP);
  subnames("phi");
  addCov(PowSstat, NULL, inversePowS); // unterscheidung nur wegen der 
  //  geschwindigkeit, also Siso ist ein sehr haeufiger Spezialfall von Sstat
  addCov(PowSnonstat);
  addlogCov(logSstat, logSnonstat, NULL);
  // addLocal(coinitS, ieinitS);  
  RandomShape(INFTY, structPowS, initPowS, doPowS, true, true, true);
  Taylor(RF_NA, RF_NA, RF_NA, RF_NA);
  TailTaylor(RF_NA, RF_NA, RF_NA, RF_NA);
  addTypeFct(TypePowS);
 

  ////////////////////////////////////////////////////////////
  includeCovModels();  
  includeOtherModels();
  ////////////////////////////////////////////////////////////

   IncludeModel("minus", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMath, PREF_TREND,
	       false,SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP, "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathminus, NULL, NULL);
  AddVariant(TrendType, PREVMODEL_I);
  set_type(DefList[currentNrCov-1].systems[0], 0, ShapeType);
 
  IncludeModel("plus", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMath, PREF_MATHDEF, 
	      false,SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP, "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathplus, NULL, NULL);
  AddVariant(TrendType, PREVMODEL_I);

  IncludeModel("div", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMath, PREF_MATHDEF, 
	      false,SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP,  "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathdiv, NULL, NULL);
  AddVariant(TrendType, PREVMODEL_I);

  IncludeModel("mult", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMath, PREF_MATHDEF, 
	      false,SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP,  "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathmult, NULL, NULL);
  AddVariant(TrendType, PREVMODEL_I);

  CONST =
  IncludePrim("const", MathDefType, 1, NULL, XONLY, PREVMODEL_I,
	      check_c, rangec, PREF_MATHDEF, 
	      SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE);
  kappanames(CONST_A_NAME, REALSXP);
  change_sortof(CONST_C, TRENDPARAM);
  addCov(Mathc, NULL, NULL);
  AddVariant(TrendType, PREVMODEL_I);
  AddVariant(NegDefType, PREVMODEL_I);
  AddVariant(TcfType, PREVMODEL_I);
  setDI(NULL, allowedItrue, NULL);

  IncludeModel("p", MathDefType, 0, 0, 3, NULL, XONLY, PARAMDEP_I,
	       checkproj, rangeproj, PREF_MATHDEF, 
	      false, SCALAR, INFDIM-1, falsch, NOT_MONOTONE);
  kappanames("proj",  INTSXP, "new", INTSXP, "factor",  REALSXP);
  change_typeof(PROJ_ISO, NN2);// "new"
  change_sortof(PROJ_FACTOR , TRENDPARAM);
  addCov(proj, NULL, NULL);
  AddVariant(TrendType, PREVMODEL_I);
  setDI(NULL, allowedIp, setproj);
  addTypeFct(Typeproj);
  
  BIND = 
  IncludeModel("c", MathDefType, 0, 0, 18, NULL, SUBMODEL_D, SUBMODEL_I,
	       check_bind, rangeMath, PREF_TREND, 
	      false, PARAM_DEP, 1, falsch, NOT_MONOTONE);
  kappanames("a", REALSXP, "b", REALSXP, "c", REALSXP, 
	     "d", REALSXP, "e", REALSXP, "f", REALSXP,
	     "g", REALSXP, "h", REALSXP, "i", REALSXP,
	     "j", REALSXP, "l", REALSXP, "m", REALSXP,
 	     "n", REALSXP, "o", REALSXP, "p", REALSXP, 
 	     "q", REALSXP, // R: 16, C: 15
	     "ncol", INTSXP, "factor", REALSXP);
  change_sortof(DefList[BIND].kappas - 1, TRENDPARAM);
  addCov(Mathbind, NULL, NULL);
  AddVariant(TrendType, SUBMODEL_I);
  assert(BIND_VARIABLES == 16);
  assert(DefList[BIND].kappas == BIND_VARIABLES + 2);
  set_type(DefList[currentNrCov-1].systems[0], 0, ShapeType);
  setDI(allowedDbind, allowedIbind, NULL);
  addTypeFct(Typebind);

 
  IncludeModel("is", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMathIs, PREF_TREND, 
	      false, SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "is", INTSXP, "y", REALSXP);
  change_typeof(IS_IS, NN1);
  addCov(Mathis, NULL, NULL);
  AddVariant(TrendType, PREVMODEL_I);
  set_type(DefList[currentNrCov-1].systems[0], 0, ShapeType);

  includeStandardMath();
  
  assert(CheckListmodel());

}
