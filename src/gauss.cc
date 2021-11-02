/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

 Copyright (C) 2005 -- 2017 Martin Schlather

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

#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "extern.h"
#include "questions.h"
#include "Processes.h"
#include "Coordinate_systems.h"
#include "variogramAndCo.h"
#include "shape.h"

void location_rules(model *cov, pref_type pref) {

  // 1. rules that depend on the the locations and the user's preferences only,
  // but not on the covariance model#
  if (COVNR != GAUSSPROC &&  COVNR != BINARYPROC) BUG;

  location_type *loc = Loc(cov);
  usr_bool exactness = GLOBAL.general.exactness; //

  unsigned int maxmem=500000000;
  int i;

  Methods Standard[Nothing] = {
     CircEmbed, CircEmbedIntrinsic, CircEmbedCutoff, SpectralTBM, TBM,
     Direct, Specific, Sequential, Trendproc, Average, Nugget, RandomCoin,
     Hyperplane
  };
  for (i=0; i<Nothing; i++) {
    pref[Standard[i]] = Nothing - i; 
  }
  
  if (P0INT(GAUSSPROC_STATONLY) == (int) True) {
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 1;
  }

  if (exactness == True) {
   pref[TBM] = pref[SpectralTBM] = pref[Average] = pref[RandomCoin] = 
      pref[Sequential] = pref[Hyperplane] 
      = LOC_PREF_NONE - 2;
  }

  if (loc->timespacedim == 1) pref[TBM] -= 2 * PREF_PENALTY;

  if (loc->distances) {
    if (loc->grid) BUG;
     for (i=CircEmbed; i<Nothing; i++) pref[i] = (i==Direct) * LOC_PREF_NONE;
  } else if (loc->grid) {
    if (exactness != True && 
	loc->totalpoints * (1 << loc->timespacedim) * sizeof(double) > maxmem){
      pref[CircEmbed] -= PREF_PENALTY;
      pref[CircEmbedIntrinsic] -= PREF_PENALTY;
      pref[CircEmbedCutoff] -= PREF_PENALTY;
    }
  } else {
   if (exactness == True){
       pref[CircEmbed] = pref[CircEmbedIntrinsic] = pref[CircEmbedCutoff] = -3;
    } else {
      pref[CircEmbed] -= PREF_PENALTY;
      pref[CircEmbedIntrinsic] = -3; 
      pref[CircEmbedCutoff] -= PREF_PENALTY;
    } 
    if (!loc->Time) pref[Sequential] = LOC_PREF_NONE;
  }
}

void mixed_rules(model *cov, pref_type locpref, 
		 pref_type pref, int *order) {
  
  assert(COVNR == GAUSSPROC);
 
 
  model *sub = cov->sub[0];
  location_type *loc = Loc(cov);
  
  int i,
    *covpref = sub->pref,
    totalpref[Nothing],
    best_dirct=GLOBAL.gauss.direct_bestvariables,
    max_variab= MAX(GLOBAL.direct.maxvariables, GLOBAL_UTILS->solve.max_chol),
    vdim = VDIM0;
  
  for (i=0; i<Nothing; i++) {
    totalpref[i] =  i == Nugget ? PREF_NUGGET * (PREF_FACTOR + 1): PREF_BEST;
    
    //     printf("%d %d", totalpref[i], covpref[i]);
    
    if (totalpref[i] > covpref[i]) totalpref[i] = covpref[i];
    
    if (totalpref[i] > PREF_NONE && locpref[i] > LOC_PREF_NONE)
      pref[i] = totalpref[i] * PREF_FACTOR + locpref[i];
    else 
      pref[i] = locpref[i] <= LOC_PREF_NONE ? locpref[i] : LOC_PREF_NONE - 4;
    
    //     printf("i=%d %d [%d, %d]\n", i, pref[i], locpref[i], totalpref[i]);
  }
  
  int vdimtot = loc->totalpoints * vdim;
  if (vdimtot > max_variab &&
      (sub->finiterange == falsch || GLOBAL_UTILS->solve.sparse == False))
      pref[Direct] = LOC_PREF_NONE - 0;

  //  printf("spam %10g %d %10g %d %d\n", GLOBAL_UTILS->solve.sparse, GLOBAL_UTILS->solve.sparse == false, sub->finiterange,sub->finiterange==false,  pref[Direct]);

  if (vdimtot <= best_dirct && totalpref[Direct] == PREF_BEST) {
    pref[Direct] = (PREF_BEST + 1 + (int) (best_dirct>=256 && vdimtot<256)) * 
      PREF_FACTOR;
  }
  else if (pref[Direct] >= PREF_NONE && GLOBAL_UTILS->solve.sparse != True) {
    double ratio = -0.1;
    if (max_variab <= DIRECT_ORIG_MAXVAR)
      ratio = (double) (vdimtot - best_dirct) / (double) max_variab;
    ratio *= FABS(ratio);
    //printf("pref[direct] = %d %d %d\n", pref[Direct],max_variab, DIRECT_ORIG_MAXVAR );
    pref[Direct] = (int) POW((double) pref[Direct], 1.0 - ratio);
    if (pref[Direct] < PREF_BEST) 
      pref[Direct] = sub->finiterange == wahr ? PREF_BEST : PREF_BEST / 2;
    //printf("ratio =%f\n", ratio);
  }
  
  if (P0INT(GAUSSPROC_STATONLY) == NA_LOGICAL && isnowPosDef(cov)) {
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 6;
  } 

  if (!isCartesian(OWNISO(0))) {
    pref[CircEmbedIntrinsic] = pref[CircEmbed] = pref[CircEmbedCutoff] =
				    LOC_PREF_NONE - 7;
    if (isAnySpherical(OWNISO(0)) && OWNTOTALXDIM < 3)
     pref[Sequential] = LOC_PREF_NONE - 8;
  }

  Ext_orderingInt(pref, Nothing, 1, order);
}

bool NAequal(double X, double Y) {
  return (ISNAN(X) && ISNAN(Y)) || X == Y;
}

int kappaBoxCoxParam(model *cov, int BC) {
  int vdim_2 = VDIM0 * 2;					
  if (PisNULL(BC)) {						
    PALLOC(BC, 2, VDIM0);				
    if (GLOBAL.gauss.loggauss) {					
      for (int i=0; i<vdim_2; i++) P(BC)[i] = 0.0;	
      GLOBAL.gauss.loggauss = false;					
    } else {								
      for (int i=0; i<vdim_2; i++)					
	P(BC)[i] = GLOBAL.gauss.boxcox[i];			
    }									
  } else {								
    int total = cov->nrow[BC] * cov->ncol[BC];	
    if (total == 1) {							
      double _boxcox = P0(BC);				
      PFREE(BC);						
      PALLOC(BC, 2, 1);					
      P(BC)[0] = _boxcox;					
      P(BC)[1] =  GLOBAL.gauss.boxcox[1];			
      total = 2;							
    }									
    if (total < vdim_2) SERR("too few parameters for boxcox given");	
    if (total % 2 != 0) SERR("number of parameters of boxcox must even"); 
    if (cov->nrow[BC]  > 2 && cov->ncol[BC]  > 1) SERR("parameters of boxcox must be given as a matrix with 2 rows"); 
    cov->ncol[BC] = total / 2;				
    cov->nrow[BC] = 2;					
    bool notok = false;						
    if (GLOBAL.gauss.loggauss) {					
      for (int i=0; i<VDIM0; i++) {				
	if ((notok = P(BC)[2*i] != RF_INF &&		
	     (P(BC)[2*i] != 0.0 || P(BC)[2*i+1] != 0.0)))
	  break;							
      }									
    } else {								
      for (int i=0; i<VDIM0; i++) {
	//printf("%d  %10g %10g   %10g %10g\n", i, GLOBAL.gauss.boxcox[2 * i], GLOBAL.gauss.boxcox[2 * i + 1], P(BC)[2 * i], P(BC)[2 * i +1]);:
 	if ((notok = (GLOBAL.gauss.boxcox[2 * i] != RF_INF &&
		      !NAequal(P(BC)[2 * i], GLOBAL.gauss.boxcox[2 * i])) ||   
	     !NAequal(P(BC)[2 * i +1], GLOBAL.gauss.boxcox[2 * i + 1])))
	  break;							
      }									
    }									
    if (notok)								
      SERR("Box Cox transformation is given twice, locally and through 'RFoptions'"); 
  }									
  for (int i=0; i<VDIM0; i++) {				
    GLOBAL.gauss.boxcox[2 * i] = RF_INF;				
    GLOBAL.gauss.boxcox[2 * i + 1] = 0.0;				
  }									
  RETURN_NOERROR;
}


void kappaGProc(int i, model *cov, int *nr, int *nc){
  *nc = i == GAUSS_BOXCOX ? SIZE_NOT_DETERMINED : 1;
  *nr = i == GAUSS_BOXCOX ? SIZE_NOT_DETERMINED 
    : i < DefList[COVNR].kappas ? 1 : -1;
}

int checkgaussprocess(model *cov) {
  //  printf("\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n");
  
  ASSERT_ONESYSTEM;
  model 
    *next=cov->sub[cov->sub[0] == NULL],
    *key = cov->key;
  //location_type *loc=Loc(cov);
  int err,
    xdim = OWNXDIM(0), // could differ from logicaldim in case of distances!
    dim = OWNLOGDIM(0);
  gauss_param *gp  = &(GLOBAL.gauss);  

  assert((Loc(cov)->distances && xdim==1) || xdim == dim);

  if (!hasAnyProcessFrame(cov) &&
      !(hasInterfaceFrame(cov) && cov->calling != NULL &&
	cov->calling->calling == NULL) &&
      !hasAnyEvaluationFrame(cov)) {
    ILLEGAL_FRAME;
  }
 						      
  kdefault(cov, GAUSSPROC_STATONLY, (int) gp->stationary_only);

  if (MAX(GLOBAL.direct.maxvariables, GLOBAL_UTILS->solve.max_chol) <
      gp->direct_bestvariables)
    SERR("maximum variables less than bestvariables for direct method");
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);

  set_maxdim(OWN, 0, INFDIM);

   if (key == NULL) {
    if (isGaussMethod(next)) {// waere GaussMethodType 
      // wird verboten
      SERR1("%.50s may not call a method", NICK(cov));
    } else {
      Types frame = hasAnyEvaluationFrame(cov) ? cov->frame :
	EvaluationType; // sehr schwach, sonst kommst aber trend nicht durch
	//GaussMethodType;

      //      printf("@@@@@@@@@@@@ %d\n", hasEvaluationFrame(cov) );
      
      if ((err = CHECKPOS2VAR(next, SUBMODEL_DEP, frame, KERNEL)) != NOERROR){

	//	if (NEXTTYPE(0)!=TrendType) {APMI(cov);BUG} else printf("jetzt C HECK_GEN\n");
	set_type(PREVSYSOF(next), 0, TrendType);
	set_iso(PREVSYSOF(next), 0, OWNISO(0));
	set_dom(PREVSYSOF(next), 0, XONLY);

	if ((err = CHECK_GEN(next,  SUBMODEL_DEP, // err mocht verwemdet
			     SUBMODEL_DEP, TrendType, false)) != NOERROR) {
	  RETURN_ERR(err); // previous error
	}
	//   	printf("jff etzt C HECK_GEN\n");
      }

      /*
      if ((err = CHECKPD2ND(next, dim, xdim, 
			    Sy mmetricOf(OWNISO(0)),// Jan 2015 S YMMETRIC,
			    SUBMODEL_DEP, frame))
	  != NOERROR) {
	if (CHECK(next, dim, dim, TrendType, XONLY, OWNISO(0),  // err mocjt verwemdet
		  SUBMODEL_DEP, frame)) RETURN_ERR(err); // previous error
      */
    }
  } else {
    if (PL >= PL_COV_STRUCTURE) {
      PRINTF("checking key in gauss process  ...\n");
    }
    Types frame = hasAnyEvaluationFrame(cov) ? cov->frame : GaussMethodType;
    if ((err = CHECK(key, dim, xdim, ProcessType, XONLY, OWNISO(0), 
		     SUBMODEL_DEP, frame)) != NOERROR) RETURN_ERR(err);
  }
  
  model *sub = cov->key==NULL ? next : key;  
  setbackward(cov, sub);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov, true)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}




void rangegaussprocess(model  VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;
  booleanRange(GAUSSPROC_STATONLY);
 }

int gauss_init_settings(model *cov) {
  model 
    *next = cov->sub[cov->sub[0] == NULL],
    *sub = cov->key == NULL ? next : cov->key;
  double sigma, var,  meanDsigma, 
    *mean = NULL,
    *variance = NULL;
  int w, v, nmP1,
    err = NOERROR,
    vdim = next->vdim[0],
    vdimSq = vdim * vdim,
    vdimP1 = vdim + 1;

  if ((variance = (double*) MALLOC(sizeof(double) * vdimSq)) == NULL ||
      (mean = (double*) CALLOC(vdim, sizeof(double))) == NULL
      ) goto ErrorHandling;

  GetInternalMean(next, vdim, mean);
  if (ISNAN(mean[0])) // GetInternalMean currently only allows ...
    GERR("Mean equals NA. Did you try a multivariate model or an incomplete (mixed) model?");
  
  //  if (false && cov->key != NULL) { // sind diese Zeilen notwendig?
  //    err = CHECK(next, OWNLOGDIM(0), OWNXDIM(0),  PosDefType, KERNEL,
  //		 , cov->vdim, EvaluationType);
  // if (err != NOERROR) goto ErrorHandling;
  // }

  // if (!next->checked) crash();
  
  assert(next->checked);
 
  if (isXonly(PREVSYSOF(next))) {
    COV(ZERO(next), next, variance);
  } else {
    assert(isKernel(PREVSYSOF(next)));
    //double *zero = ZERO(next);
    //  NONSTATCOV(zero, zero, next, variance);
    for (v=0; v<vdimSq; variance[v++] = 0.0); // just a dummy
  }

  if (cov->q == NULL) QALLOC(vdim);
  if ((err = alloc_mpp_M(cov, 2)) != NOERROR) goto ErrorHandling;
 
  nmP1 = cov->mpp.moments + 1;
  for (v=w=0; w<vdimSq; w+=vdimP1, v++) { 
    sigma = SQRT(variance[w]);
    meanDsigma = sigma == 0 ? RF_INF : mean[v] / sigma;
    
    //the following line is SQRT(2 pi s^2) int_0^infty x EXP(-0.5 (x-mu)^2/s^2)
    var = sigma * INVSQRTTWOPI * EXP(-0.5 * meanDsigma * meanDsigma) +
      mean[v] * /* correct next lines the 2nd factor is given */
      pnorm(0.0, mean[v], sigma, false, false);
    var = 1.0 / (var * var);
    cov->q[v] = var;
    
    // trend um factor var korrigieren in rfmodel
    // eigentliches model um factor var korrigieren
    // so dass ein Gausssches Zufallsfeld rauskommt, dessen Positivteil
    // zumindest im Ursprung den Erwartungswert 1 hat.

    cov->mpp.maxheights[v] = 
      GLOBAL.extreme.standardmax * sigma + ((mean[v]>0) ? mean[v] : 0);//approx!
  
    int idx = v * nmP1;
    cov->mpp.mM[idx + 0] = cov->mpp.mMplus[idx + 0] = 1.0;       
    cov->mpp.mMplus[idx + 1] = 
      sigma * INVSQRTTWOPI * EXP(-0.5 * mean[v] * mean[v]) 
      + mean[v] * pnorm(-mean[v], 0.0, 1.0, false, false);
    //(2pi)^-0.5 int x EXP(-(x-m)/2s^2)
    cov->mpp.mM[idx + 1] = 0.0;
    cov->mpp.mM[idx + 2] = variance[w];
    // todo: cov->mpp.mMplus[2] berechnen
  }
 
   ReturnOtherField(cov, sub);

 ErrorHandling:
  FREE(variance);
  FREE(mean);

  RETURN_ERR(err);
}


int struct_gaussmethod(model *cov, model **newmodel) {
  // uebernimmt struct- und init-Aufgaben !!
  model *next = cov->sub[0]; // nicht cov->sub[cov->sub != NULL]
  SPRINTF(cov->base->error_loc, "simulation procedure for %.50s", NAME(next));
  ASSERT_ONESYSTEM;
  location_type *loc=Loc(cov);
  int
    nr = COVNR,
    err = ERROROUTOFMETHODLIST,
    xdim = OWNXDIM(0), // could differ from logicaldim in case of distances!
    logicaldim = OWNLOGDIM(0);
  assert((loc->distances && xdim==1) || xdim == logicaldim);

  cov->fieldreturn = wahr;
  ASSERT_NEWMODEL_NULL;
  assert(cov!=NULL);

  if ((logicaldim != PREVXDIM(0) || logicaldim != xdim) &&
      (!loc->distances || PREVXDIM(0)!=1)) {
    BUG; // RETURN_ERR(ERRORDIM);
  }

  if (next!= NULL && !isnowVariogram(next)) {
    SERR("submodel not a covariance function");
  }
  
  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
  if ((err = covcpy(&(cov->key), cov)) != NOERROR) RETURN_ERR(err); //!!cov


  defn *C = DefList + nr; 
  set_nr(SYSOF(cov->key), C->Specific);


  // cov!! nicht cov->key!!  ??? Nein, sonst ist Coin u.ae. nicht gesetzt.
  model 
    *firstmethod = cov->key; // z.B. coins_intern
  int idx = firstmethod->sub[0] == NULL;
  model
    *sub = firstmethod->sub[idx]; // z.B. oberste fkt shape von coins

  if ((err = CHECK_PASSTF(cov->key, // 11.1.19 cov WARUM ??
			  GaussMethodType, VDIM0, 
			// frame // 19.5.2013 auskommentiert und ersetzt durch
			  // WARUM ??
			  GaussMethodType)) != NOERROR && !isAnyDollar(sub)) {
     RETURN_ERR(err);
  }
  //  if ((err = CHECK( c o v, dim, xdim, GaussMethodType, OWNDOM(0),
  //		   OWNISO(0), cov->vdim, 
  //		   // frame // 19.5.2013 auskommentiert und ersetzt durch
  //		   AnyType  // wegen spectral.cc
  //		   )) != NOERROR) {
  //     RETURN_ERR(err);
  //  }


  if (err == NOERROR) err = STRUCT(cov->key, NULL);
  if (err != NOERROR && !isAnyDollar(sub)) RETURN_ERR(err);
  
  model *dollar,
    *lastmethod= equalsnowGaussMethod(sub) ? sub : firstmethod;
  
  if (err != NOERROR) {
    assert(isAnyDollar(sub));
    dollar = cov->key = lastmethod->sub[idx]; // $ jetzt im cov->key
    lastmethod->sub[idx] = dollar->sub[0]; // eigentl. Modell jetzt anhaengen
    SET_CALLING(lastmethod->sub[idx], lastmethod);   
    dollar->sub[0] = firstmethod;       // auf $ setzen
    SET_CALLING(firstmethod, dollar); 
    SET_CALLING(dollar, cov);              // $ mit cov verbinden
    dollar->prevloc = cov->prevloc;
    set_nr(SYSOF(dollar), DOLLAR_PROC);

    Types frame =
      nr == AVERAGE_USER || nr == RANDOMCOIN_USER
      ? PoissonGaussType
      : GaussMethodType;
  /// cov nicht cov->key !!! ???? OK. Da sonst cov u.U. nicht gesetzt
    if ((err = CHECK_PASSTF(cov->key, // 11.1.19 cov WARUM ??
			    GaussMethodType, //?! 8.12.17
			    VDIM0, frame)) != NOERROR) {
      RETURN_ERR(err);
    }
    if ((err = STRUCT(cov->key, NULL)) != NOERROR) { RETURN_ERR(err); }
  }

  Methods m;
  for (m=CircEmbed; m <= Hyperplane; m = (Methods) ((int) (m + 1))) {
    if (gaussmethod[(Methods) m] == COVNR) break;
  }
  cov->key->method = m;
  
  RETURN_NOERROR;
}


  
int struct_gaussprocess(model *cov, model **newmodel) {
  // uebernimmt struct- und init-Aufgaben !!
  KEY_type *KT = cov->base;  
  model *next = cov->sub[cov->sub[0] == NULL];
  SPRINTF(KT->error_loc, "simulation method for %.50s", NAME(next));
  ASSERT_NEWMODEL_NULL;
  assert(cov!=NULL);

  if (hasAnyEvaluationFrame(cov)) {
    if (VDIM0 > MAXGAUSSVDIM)
      SERR2("The number of variables, %d, exceeds the maximum number %d",
	    VDIM0,  MAXGAUSSVDIM);
    return struct_gauss_logli(cov); //loglihood
  } else if (!hasAnyProcessFrame(cov) &&
	     !(hasInterfaceFrame(cov) && cov->calling != NULL &&
	       cov->calling->calling == NULL)
	     ) RETURN_ERR(ERRORFAILED);

 cov->fieldreturn = wahr;
  location_type *loc = Loc(cov);
  pref_type locpref, pref;
  int order[Nothing], i;
  char dummy[100 * (int) Nothing];
  int 
    xdim = OWNXDIM(0), // could differ from tsdgetim in case of distances!
    logicaldim = OWNLOGDIM(0);
  assert((loc->distances && xdim==1) || xdim == logicaldim);
#define  MAXFAILMSG 9
  char FailureMsg[MAXFAILMSG][80] = {"unknown reason",
     "total number of points > maximum value",
     "non-domain model not allowed",
     "not an exact method as required by user",
     "locations not on a grid",
     "denied by cov. model or by user's pref.",
     "no analytic solution of Abelian integral",
     "intrinsic fields denied by user",
     "Unknown reason"
    };
  bool all_PREF_NONE;
  Methods unimeth = Forbidden;
  int 
    meth, 
    zaehler = 0,
    nr = NEXTNR;
#define nsel 4

  ASSERT_ONESYSTEM;
  if( (logicaldim != PREVXDIM(0) || logicaldim != xdim) &&
      (!loc->distances || PREVXDIM(0) !=1)) {
    BUG; //RETURN_ERR(ERRORDIM);
  }

  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
  if (!isnowVariogram(next) && !equalsnowTrend(next)) {
    SERR("submodel must be a covariance function");
  }
 

  //////////////////////////////////////////////////////////////////////
  // GaussMethodType; cov->key not given

  { // only for error reporting !
    model *sub = cov; 
    while(isAnyDollar(sub)) sub = sub->sub[0];
    SPRINTF(KT->error_loc, "Searching a simulation method for '%.50s'",
	    NICK(sub));
  }

  location_rules(cov, locpref); // hier cov
  mixed_rules(cov, locpref, pref, order);

  if (PL >= PL_RECURSIVE) {
    defn *N = DefList + nr;
    LPRINT("\n%.50s:%.50s\n", NAME(cov), NAME(next));
    for (i=0; i < Nothing; i++) {
      LPRINT("%-10s: base=%1d/%1d  locpref=%3d  totpref=%5d\n", 
	     METHOD_NAMES[i], N->pref[i], next->pref[i], 
	     (2 * locpref[i] < LOC_PREF_NONE) ? LOC_PREF_NONE + 1 : locpref[i],
	     pref[i]);
    }
    LPRINT("initstandard %.50s (%d) [pref value %d]\n", 
	   DefList[nr].name, nr, pref[order[Nothing -1]]);
  }
  all_PREF_NONE = true;
  for (i=0; i<Nothing; i++) all_PREF_NONE &= pref[i] == PREF_NONE;
  
  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
  int err;
  if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);

  if (all_PREF_NONE) {
    err = pref[Nothing] == PREF_NONE ? ERRORINVALIDMODEL : ERRORODDMODEL;
    goto ErrorHandling;
  } 

  err = CERROROUTOFMETHODLIST; // in case none of the methods are available 
  for (i = Nothing - 1; i>=0 && pref[order[i]] > PREF_NONE; i--) {

    //  assert(i ==Nothing-1);
 
    if (PL >= PL_BRANCHING && i < Nothing - 1) {
      LPRINT("'%.50s%.50s' failed for '%.50s'\n", 
	     CAT_TYPE_NAMES[GaussMethodType], METHOD_NAMES[unimeth],NICK(next));
    }

    zaehler++;
    unimeth = (Methods) order[i];

    if (PL >= PL_BRANCHING) {
      LPRINT("trying '%.50s%.50s' for '%.50s'\n", 
	     CAT_TYPE_NAMES[GaussMethodType], METHOD_NAMES[unimeth],NICK(next));
    }
    
    meth = gaussmethod[unimeth];      
  
    addModel(&(cov->key), meth);
    SET_CALLING(cov->key, cov);
    
    model *key = cov->key;
    key->prevloc = PLoc(cov);
    Types frame =
      MODELNR(key) == AVERAGE_INTERN ? PoissonGaussType :      
      hasAnyEvaluationFrame(cov) ? cov->frame : GaussMethodType;

    if ((err = CHECK(key, logicaldim, xdim, GaussMethodType, OWNDOM(0),
		     OWNISO(0), cov->vdim, frame)) == NOERROR) {

      //PMI(key);
      
      if ((err = STRUCT(key, NULL)) <= NOERROR) {

        //PMI(cov);
       //  assert(SYSMODEL(SYSOF(key)) != DIRECT); printf("xx\n");
      
	if (!key->initialised) {
	  if ((err = CHECK(key, logicaldim, xdim, GaussMethodType, OWNDOM(0),
			   OWNISO(0), cov->vdim, frame)) == NOERROR){
	    CHECKED;
	    key->method = unimeth;
	    
	    NEW_STORAGE(gen);
	    err = INIT(key, isPoissonGauss(frame) ? 2 : 0, cov->Sgen);

	    if (err == NOERROR) {
	      if (PL >= PL_RECURSIVE) {
		LPRINT("returning to higher level...\n");
	      }
	      goto Initialised;
	    }
	  }
	}     
      }
   
      if (PL >= PL_RECURSIVE) {  
	LPRINT("back from %.50s: err =%d [%d, %.50s]\n",
	       METHOD_NAMES[unimeth], err, unimeth, METHOD_NAMES[unimeth]);     
	
	if (PL > PL_ERRORS) {
	  char info[LENERRMSG]; errorMSG(err, info);
	  PRINTF("error in gauss (init): %.50s\n",info);
	  // BUG;//
	}
      } // else if (PL >= PL_BRANCHING) M ERR(err);
    } else {
      if (PL > PL_ERRORS) {
	char info[LENERRMSG]; errorMSG(err,info);
	PRINTF("error in check: %.50s\n", info);
      }
    }

    key = key->sub[0];
    COV_DELETE_WITHOUTSUB(&(cov->key), cov);  
    
    cov->key=key;
    SET_CALLING(key, cov);
  }
  
  if (PL >= PL_RECURSIVE) {
    char *PREF_FAILURE = KEYtypeOf(cov)->PREF_FAILURE;
    STRCPY(PREF_FAILURE, "");
    int p, lp;
#define NMAX 14
    char lpd[255], pd[255], names[NMAX];
    names[NMAX-1] = '\0';
    for (i=0; i<Nothing; i++) {
      lp = locpref[i] - LOC_PREF_NONE;
      p = pref[i] - LOC_PREF_NONE;
      if (lp>0) {
	STRCPY(lpd, "");
      } else {
	SPRINTF(lpd, "%.50s (locp) ", FailureMsg[MAX(0, MIN(MAXFAILMSG-1, 1-lp))]);
      }
      if (p>0 || p == lp) {
	if (lp>0) STRCPY(pd, "(method specific failure)");
	else STRCPY(pd, ""); // STRCPY(pd, "(method specific failure)");
      } else {
	SPRINTF(pd, "%.50s (pref)", FailureMsg[MAX(0, MIN(MAXFAILMSG-1, 1-p))]);
      }
      strcopyN(names, METHOD_NAMES[i], NMAX-1);
      SPRINTF(dummy, "%.*s %-13.13s: %.5s%.70s\n", 100 * (int) Nothing -100,
	      PREF_FAILURE, names, lpd, pd);
      STRCPY(PREF_FAILURE, dummy);
    }
  }
   
  if (err != NOERROR) goto ErrorHandling;
  
 Initialised:
  // nachfolgende Zeile muss stehen, da Init in struct aufgerufen!
  err = gauss_init_settings(cov);

  cov->simu.active = cov->initialised = err == NOERROR;	  

 ErrorHandling:
  if (err <= NOERROR || zaehler==0){
    if (err <= NOERROR) {
      if (cov->key == NULL) next->simu.active = true;
      else cov->key->simu.active = true;
    }
    RETURN_ERR(err);
  }

  if (zaehler==1) {
    SPRINTF(KT->error_loc, "Only 1 method found for '%.50s', namely '%.50s', that comes into question.\nHowever, this method failed (err=%d)", NICK(next), METHOD_NAMES[unimeth], err);
    RETURN_ERR(err);
  }
  
  { // only for error reporting !
    model *sub = next; 
    while(isAnyDollar(sub)) sub = sub->sub[0];
    SPRINTF(KT->error_loc, "searching a simulation method for '%.50s'",
	    NICK(sub));
  }

  RETURN_ERR(ERROROUTOFMETHODLIST);
}


int init_gaussprocess(model *cov, gen_storage *s) {
  /// ACHTUNG struct_gaussprocess hat bereits init aufgerufen!

  // direkter Aufruf, u.a. durch RPtbm, average, cutoff, intrinsic, hyperplane,
  //                    nugget, spectral !!!
  if (hasAnyEvaluationFrame(cov)) {
    assert(cov->key == NULL);
    if (isnowVariogram(cov->sub[0])) return NOERROR;
    return INIT(cov->sub[0], 0, s);
  }

  model 
    *key = cov->key;
  int err;

  assert(key != NULL);

  if ((err = INIT(key, 0, s)) != NOERROR) RETURN_ERR(err);
 	    
 if ((err = gauss_init_settings(cov)) != NOERROR) {
   RETURN_ERR(err);
  }
  
  cov->simu.active = key->simu.active = true;

  RETURN_NOERROR;
}


   

void do_gaussprocess(model *cov, gen_storage *s) {
  assert(s != NULL);
  // reopened by internal_dogauss
  errorloc_type errorloc_save;
  double *res = cov->rf;
  int i,
    vdimtot = Gettotalpoints(cov) * VDIM0 ;
  model *key = cov->key ;
  KEY_type *KT = cov->base;
  SAVE_GAUSS_TRAFO;
 
  STRCPY( errorloc_save,  KT->error_loc);

  if (cov->simu.pair) {
    for (i=0; i<vdimtot; i++) res[i] = -res[i];
    cov->simu.pair = false;
    return;  
  } else {
    cov->simu.pair = GLOBAL.gauss.paired;
  }

  assert(key != NULL);

  // falls gaussprocess als solcher aufgerufen, dann cov->Sgen != NULL
  // falls ueber *_USER, dann i.a.(?!) nicht gesetzt, sondern s:
  DO(key, cov->Sgen == NULL ? s : cov->Sgen); 
 
  // (x^\lambda_1-1)/\lambda_1+\lambda_2

  BOXCOX_INVERSE;
  STRCPY( KT->error_loc, errorloc_save);
}




//////////////////////////////////////////////////////////////////////
// Binary Gauss
#define BINARY_THRESHOLD (GAUSSPROC_LAST + 1)

void kappa_binaryprocess(int i, model VARIABLE_IS_NOT_USED *cov, 
			 int *nr, int *nc){
  kappaGProc(i, cov, nr, nc);
  if (i == BINARY_THRESHOLD) *nr = SIZE_NOT_DETERMINED;
}

int checkbinaryprocess(model *cov) {
  // Achtung! binary hat Doppelbedeutung: binary-Gaussian als default
  // und binary angewandt auf x-bel. process.
  model
    *key = cov->key,
    *next = cov->sub[0],
    *sub = key != NULL ? key : next;
  double v;
  int err = NOERROR;
  if (PisNULL(BINARY_THRESHOLD)) kdefault(cov, BINARY_THRESHOLD, 0);

  if (key == NULL && !isProcess(next)) {
    if ((err = checkgaussprocess(cov)) != NOERROR) {
      RETURN_ERR(err);
    }
    
    COV(ZERO(sub), sub, &v);
    if (v != 1.0) 
      SERR("binaryian requires a correlation function as submodel.");
  } else {
    if ((err = CHECK_PASSTF(sub, ProcessType, SUBMODEL_DEP,
			    hasAnyEvaluationFrame(cov) ? cov->frame :
			    NormedProcessType
			    )) != NOERROR) 
      RETURN_ERR(err);
    //    if ((err = CHECK(sub, cov->tsdim, cov->xdimprev, ProcessType, 
    //		     OWNDOM(0), OWNISO(0),
    //		     SUBMODEL_DEP, 
    //		     frame)) != NOERROR) RETURN_ERR(err);
    setbackward(cov, sub);
  }

  VDIM0 = sub->vdim[0];
  VDIM1 = sub->vdim[1];
  RETURN_NOERROR;
}

int struct_binaryprocess(model *cov, model VARIABLE_IS_NOT_USED **newmodel) {
  model
    *next = cov->sub[0];
  int err;
  
  if (isnowVariogram(next)) {
    assert(cov->key == NULL);
    err = covcpy(&(cov->key), cov);

    // darauffolgende Zeile absichern:
    if (DefList[COVNR].kappas != 3 || DefList[GAUSSPROC].kappas != 2) BUG;
    if (cov->key != NULL && !PARAMisNULL(cov->key, BINARY_THRESHOLD)) {
      PARAMFREE(cov->key, BINARY_THRESHOLD);
      PARAMtoNULL(cov->key, BINARY_THRESHOLD);
    }

    if (err != NOERROR)  RETURN_ERR(err); //!!cov

    SET_NR(cov->key, GAUSSPROC);
    err = CHECK_PASSTF(cov->key, ProcessType, SUBMODEL_DEP, GaussMethodType);
    //    err = CHECK(cov->key, cov->tsdim, cov->xdimprev, ProcessType, 
    //		OWNDOM(0), OWNISO(0), SUBMODEL_DEP, GaussMethodType);
    if (err != NOERROR)  RETURN_ERR(err);
    err = STRUCT(cov->key, NULL);
    RETURN_ERR(err);
  }
  else return STRUCT(next, NULL);
}

int init_binaryprocess( model *cov, gen_storage *s) {
  double sigma, 
    *mean = NULL, 
    *variance = NULL,
    *p = P(BINARY_THRESHOLD);
  model 
    *key = cov->key,
    *next=cov->sub[0],
    *sub = key != NULL ? key : next;
  int v, w, nmP1, pi,
    npi = cov->nrow[BINARY_THRESHOLD],
    err = NOERROR,
     vdim = next->vdim[0],
    vdimSq = vdim * vdim,
    vdimP1 = vdim + 1;
  assert(next->vdim[0] == next->vdim[1]);
 
  if ((variance = (double*) MALLOC(sizeof(double) * vdimSq)) == NULL ||
      (mean = (double*) CALLOC(vdim, sizeof(double))) == NULL
      ) goto ErrorHandling;

  if ((err = INIT(sub, 0, s)) != NOERROR) goto ErrorHandling;
  cov->rf = sub->rf;
  cov->origrf = false;
 
  if (isnowVariogram(next) || NEXTNR == GAUSSPROC) {
    GetInternalMean(next, vdim, mean);
    if (ISNAN(mean[0])) // GetInternalMean currently only allows ...
      GERR1("'%.50s' currently only allows scalar fields - NA returned",
	    NICK(cov));
     if (cov->mpp.moments >= 1)  {
      model *sub0 = NEXTNR==GAUSSPROC ? next->sub[0]: next;
      COV(ZERO(sub0), sub0, variance);
    }
    nmP1 = cov->mpp.moments + 1;
    for (pi=v=w=0; w<vdimSq; w+=vdimP1, v++, pi = (pi + 1) % npi ) { 
      int idx = v * nmP1;
      cov->mpp.maxheights[v] = 1.0;
      if (cov->mpp.moments >= 0) { 
	assert(cov->mpp.mM != NULL);
	cov->mpp.mM[idx + 0] = cov->mpp.mMplus[idx + 0] = 1.0; 
	if (cov->mpp.moments >= 1) {
	  if (variance[w]==0.0)
	    GERR1("Vanishing sill not allowed in '%.50s'", NICK(next));
	  sigma = SQRT(variance[w]);
	  cov->mpp.mM[idx + 1] = cov->mpp.mMplus[idx + 1] = 
	    pnorm(p[pi], mean[v], sigma, false, false);
	  int i;
	  for (i=2; i<= cov->mpp.moments; i++)
	    cov->mpp.mM[idx + i] = cov->mpp.mMplus[idx + i] = cov->mpp.mM[idx + 1];
	}
      }
    }
  }

  cov->fieldreturn = wahr;
  cov->simu.active = err == NOERROR;

 ErrorHandling:
  FREE(variance);
  FREE(mean);
 
  RETURN_ERR(err);
}



void do_binaryprocess(model *cov, gen_storage *s){
  // reopened by internal_dogauss
  long j,
    tot = Loc(cov)->totalpoints,
    endfor = tot;
  assert(cov->ownloc == NULL);
  int  i,
    pi, 
    npi = cov->nrow[BINARY_THRESHOLD],
    vdim = VDIM0;
  double 
    threshold,
    *p = P(BINARY_THRESHOLD),
    *rf = cov->rf;
  model 
    *next=cov->sub[0];
  assert(VDIM0 == VDIM1);
  if (isnowVariogram(next)) {
    do_gaussprocess(cov, s);
  } else DO(next, s);

  for (j=pi=i=0; i<vdim; i++, pi = (pi + 1) % npi, endfor += tot) {
    threshold = p[pi];
    // printf("endfor %d %10g\n", endfor, threshold);
    if (threshold > RF_NEGINF && threshold < RF_INF)
      for ( ; j<endfor; j++) rf[j] = (double) (rf[j] >= threshold);
  }
}



void rangebinaryprocess(model *cov, range_type *range) { 
  rangegaussprocess(cov, range);

  range->min[BINARY_THRESHOLD] = RF_NEGINF;
  range->max[BINARY_THRESHOLD] = RF_INF;
  range->pmin[BINARY_THRESHOLD] = -3;
  range->pmax[BINARY_THRESHOLD] = 3;
  range->openmin[BINARY_THRESHOLD] = false;
  range->openmax[BINARY_THRESHOLD] = false;
}





//////////////////////////////////////////////////////////////////////
// Chisq Gauss
#define CHISQ_DEGREE 1
#define T_NU CHISQ_DEGREE


int checkchisqprocess(model *cov) {
  ASSERT_ONESYSTEM;
  model
    *key = cov->key,
    *next = cov->sub[0];
  double *v = NULL;
  int err = NOERROR,
   xdim = OWNXDIM(0), // could differ from logicaldim in case of distances!
    logicaldim = OWNLOGDIM(0);
  assert((Loc(cov)->distances && xdim==1) || xdim == logicaldim);

  if (PisNULL(CHISQ_DEGREE)) SERR("degree of freedom must be given");
  if (key == NULL) {
    //    if (!isGaussMethod(next) && !isVariogram(next))
    //      SERR1("Gaussian process required, but '%.50s' obtained", NICK(cov));
    if ((err = CHECK(next, logicaldim, xdim, ProcessType, XONLY, OWNISO(0), 
		     SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
      if ((err = CHECK(next, logicaldim, xdim, 
		       isCartesian(OWNISO(0)) ? VariogramType : PosDefType, 
		       KERNEL, 
		       CoordinateSystemOf(OWNISO(0)),
		       SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
      RETURN_ERR(err);
      }
    }
    
    int 
      vdim = next->vdim[0],
      vdimSq = vdim * vdim;
    assert(vdim > 0);

    if ((v = (double*) MALLOC(sizeof(double) * vdimSq)) == NULL)
      RETURN_ERR(ERRORMEMORYALLOCATION);

    model *sub = next;
    while (sub != NULL && isnowProcess(sub)) sub=sub->sub[0];
    if (sub == NULL) BUG;
    if (isXonly(PREVSYSOF(sub))) {
      COV(ZERO(next), next, v);
    } else {
      assert(isKernel(PREVSYSOF(sub)));
      double *zero = ZERO(next);
      NONSTATCOV(zero, zero, next, v);
    }
    
    int w,
      vdimP1 = vdim + 1;
    for (w=0; w<vdimSq; w+=vdimP1) {
      if (v[w] != 1.0) {
	FREE(v);
	SERR("chisq requires a correlation function as submodel.");
      }
    }
    FREE(v);
    
    VDIM0 = sub->vdim[0];
    VDIM1 = sub->vdim[1];
     //no need to setbackward
  } else {
    if ((err = CHECK_PASSTF(key, ProcessType, SUBMODEL_DEP, GaussMethodType))
	!= NOERROR) RETURN_ERR(err);
    //   if ((err = CHECK(key, cov->tsdim,  cov->xdimprev, ProcessType,
    //		       OWNDOM(0), OWNISO(0),
    //		       SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
    setbackward(cov, key);
  }

  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
 
  RETURN_NOERROR;
}
 
int struct_chisqprocess(model *cov, 
			model VARIABLE_IS_NOT_USED **newmodel) {
  model
    *next = cov->sub[0];
  int err;

  if (isnowVariogram(next)) {
    assert(cov->key == NULL);
    if ((err = covcpy(&(cov->key), next)) > NOERROR) RETURN_ERR(err);
    addModel(&(cov->key), GAUSSPROC);
    SET_CALLING(cov->key, cov);
    if ((err = CHECK_PASSFRAME(cov->key, GaussMethodType)) != NOERROR) RETURN_ERR(err);
    //    if ((err = CHECK(cov->key, cov->tsdim,  cov->tsdim, ProcessType,
    //		     cov->domprev, cov->isoprev,
    //		     SUBMODEL_DEP, GaussMethodType)) != NOERROR) RETURN_ERR(err);
   return STRUCT(cov->key, NULL);
  } else {
    return STRUCT(next, NULL);
  }
}

int init_chisqprocess(model *cov, gen_storage *s) {
  double // sigma,
    mean, m2, 
    variance;
  model 
    *key = cov->key,
    *next = cov->sub[0],
    *sub = key == NULL ? next : key;
  int i, err, nmP1,
    nmP1sub = sub->mpp.moments + 1,
    vdim = VDIM0;

  cov->simu.active = false;
  if ((err = INIT(sub, 
		  DefList[COVNR].range == rangechisqprocess ? 2 :
		  DefList[COVNR].range == rangetprocess ? 1 : 9999
		  , s)) != NOERROR) RETURN_ERR(err);		  
 
  nmP1 = cov->mpp.moments + 1;
  for (i=0; i<vdim; i++) { 
    int 
      idx = i * nmP1,
      idxsub = i * nmP1sub;
    mean = sub->mpp.mM[idxsub + 1];
    m2 = sub->mpp.mM[idxsub + 2];
    variance =  m2 - mean * mean; 
    if (variance==0.0) 
      SERR1("Vanishing sill not allowed in '%.50s'", NICK(sub));
    if (ISNAN(mean)) // GetInternalMean currently only allows ...
      SERR1("'%.50s' currently only allows scalar fields -- NA returned", 
	   NICK(cov));
    cov->mpp.maxheights[i] = GLOBAL.extreme.standardmax * 
      GLOBAL.extreme.standardmax * m2;

    if (cov->mpp.moments >= 0) { 
      assert(cov->mpp.mM != NULL);
      cov->mpp.mM[idx + 0] = cov->mpp.mMplus[idx + 0] = 1.0; 
      if (cov->mpp.moments >= 1) {
	cov->mpp.mMplus[idx + 1] = 		  
	  DefList[COVNR].range == rangechisqprocess ? m2 
	  : DefList[COVNR].range == rangetprocess ? RF_NAN
	  : RF_NAN;
	cov->mpp.mM[idx + 1] = RF_NA; // to do: richtiger Wert
	if (cov->mpp.moments >= 2)
	  cov->mpp.mM[idx + 2] = 3 * variance * RF_NA;//check the 4th moment of a Gaussian, to do
      }
    }
  }


  if (DefList[COVNR].range == rangechisqprocess) ReturnOwnField(cov);
  else if (DefList[COVNR].range == rangetprocess) ReturnOtherField(cov, sub);
  else BUG;
  cov->simu.active = true;
 
  RETURN_NOERROR;
}



void do_chisqprocess(model *cov, gen_storage *s){
  // reopened by internal_dogauss
  int  i, f,
    degree = P0INT(CHISQ_DEGREE),
    totvdim = Loc(cov)->totalpoints * VDIM0;
 assert(cov->ownloc == NULL);
  model 
    *next=cov->sub[0],
    *key = cov->key,
    *sub = key == NULL ? next : key;
  double 
    *subrf = sub->rf,
    *res = cov->rf;

  for (i=0; i<totvdim; i++) res[i] = 0.0;
  for (f=0; f<degree; f++) {
    DO(sub, s);   
    for (i=0; i<totvdim; i++) {
      double dummy = subrf[i];
      res[i] += dummy * dummy;
    }
  }
  BOXCOX_INVERSE;
}

void rangechisqprocess(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;
  
  range->min[CHISQ_DEGREE] = 1;
  range->max[CHISQ_DEGREE] = RF_INF;
  range->pmin[CHISQ_DEGREE] = 1;
  range->pmax[CHISQ_DEGREE] = 1000;
  range->openmin[CHISQ_DEGREE] = false;
  range->openmax[CHISQ_DEGREE] = true;
}



//////////////////////////////////////////////////////////////////////
//  multivariate t distribution





void do_tprocess(model *cov, gen_storage *s){
  // reopened by internal_dogauss
  int  i,
    n = Loc(cov)->totalpoints * VDIM0;
  model 
    *next=cov->sub[0],
    *key = cov->key,
    *sub = key == NULL ? next : key;
  double 
    nu = P0(T_NU),
    factor = SQRT(nu / rgamma(0.5 * nu, 0.5)), // note inverse!! see def rgamma
    *res = cov->rf;

  DO(sub, s);  
  for (i=0; i<n; res[i++] *= factor);

  BOXCOX_INVERSE;
}

void rangetprocess(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;
 
  range->min[T_NU] = 0;
  range->max[T_NU] = RF_INF;
  range->pmin[T_NU] = 1;
  range->pmax[T_NU] = 1000;
  range->openmin[T_NU] = true;
  range->openmax[T_NU] = true;
}

