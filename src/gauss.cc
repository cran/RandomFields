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

#include <Rmath.h>
#include "RF.h"
#include "Operator.h"
#include "shape_processes.h"
#include "Coordinate_systems.h"
#include "variogramAndCo.h"



void location_rules(cov_model *cov, pref_type pref) {

  // 1. rules that depend on the the locations and the user's preferences only,
  // but not on the covariance model#
  //PMI(cov);
  if (cov->nr != GAUSSPROC &&  cov->nr != BINARYPROC) BUG;

  location_type *loc = Loc(cov);
  double exactness = GLOBAL.general.exactness; //

  unsigned int maxmem=500000000;
  int i;

  Methods Standard[Nothing] = {
     CircEmbed, CircEmbedIntrinsic, CircEmbedCutoff, SpectralTBM, TBM,
     Direct, Specific, Sequential, TrendEval, Average, Nugget, RandomCoin,
     Hyperplane
  };
  for (i=0; i<Nothing; i++) {
    pref[Standard[i]] = Nothing - i; 
  }
  
  if (P0INT(GAUSSPROC_STATONLY) > 0) {
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 1;
  }

  if (!ISNA(exactness) && exactness) {
    pref[TBM] = pref[SpectralTBM] = pref[Average] = pref[RandomCoin] = 
      pref[Sequential] = pref[Hyperplane] 
      = LOC_PREF_NONE - 2;
  }

  if (loc->timespacedim == 1) pref[TBM] -= 2 * PREF_PENALTY;

  if (loc->distances) {
    if (loc->grid) BUG;
     for (i=CircEmbed; i<Nothing; i++) 
       if (i!=Direct // && i!=Nugget
	  ) pref[i] = LOC_PREF_NONE;
  } else if (loc->grid) {
    if ((!ISNA(exactness) || !exactness) && 
	loc->totalpoints * (1 << loc->timespacedim) * sizeof(double) > maxmem){
      pref[CircEmbed] -= PREF_PENALTY;
      pref[CircEmbedIntrinsic] -= PREF_PENALTY;
      pref[CircEmbedCutoff] -= PREF_PENALTY;
    }
  } else {
    if (ISNA(exactness) || !exactness){
      pref[CircEmbed] -= PREF_PENALTY;
      pref[CircEmbedIntrinsic] = -3; 
      pref[CircEmbedCutoff] -= PREF_PENALTY;
    } else {
      pref[CircEmbed] = pref[CircEmbedIntrinsic] = pref[CircEmbedCutoff] = -3;
    }
    if (!loc->Time) pref[Sequential] = LOC_PREF_NONE;
  }
}

void mixed_rules(cov_model *cov, pref_type locpref, 
		 pref_type pref, int *order) {
  
  assert(cov->nr == GAUSSPROC);
 
 
  cov_model *sub = cov->sub[0];
  location_type *loc = Loc(cov);
  
  int i,
    *covpref = sub->pref,
    totalpref[Nothing],
    best_dirct=GLOBAL.gauss.direct_bestvariables,
    max_variab= GLOBAL.direct.maxvariables,
    vdim = cov->vdim[0];
  
  for (i=0; i<Nothing; i++) {
     totalpref[i] =  i == Nugget ? PREF_NUGGET * (PREF_FACTOR + 1): PREF_BEST;

    if (totalpref[i] > covpref[i]) {
       totalpref[i] = covpref[i];
     }

     pref[i] = (totalpref[i] > PREF_NONE && locpref[i] > LOC_PREF_NONE)
       ? totalpref[i] * PREF_FACTOR + locpref[i]
       : locpref[i] <= LOC_PREF_NONE ? locpref[i] :  LOC_PREF_NONE - 4;
   }


  if (loc->totalpoints * vdim > max_variab &&
      (sub->finiterange == false || GLOBAL_UTILS->solve.sparse == false))
      pref[Direct] = LOC_PREF_NONE - 0;

  //  printf("spam %f %d %f %d %d\n", GLOBAL_UTILS->solve.sparse, GLOBAL_UTILS->solve.sparse == false, sub->finiterange,sub->finiterange==false,  pref[Direct]);
  //  printf("spam %f %d %d %d %d\n", GLOBAL_UTILS->solve.sparse, GLOBAL_UTILS->solve.sparse == false, sub->finiterange == true, pref[Direct], PREF_FACTOR);
  //PMI(sub);

  int vdimtot = loc->totalpoints * vdim;
  if (vdimtot <= best_dirct && totalpref[Direct] == PREF_BEST) {
    pref[Direct] = (PREF_BEST + 1) * PREF_FACTOR;
  }
  else if (pref[Direct] >= PREF_NONE && GLOBAL_UTILS->solve.sparse != true) {
#define MONTRAFO(x) ((x) * (x))
    //#define MONTRAFO log   
    //   bool orig_max = max_variab == DIRECT_ORIG_MAXVAR;
    double ratio = max_variab <= DIRECT_ORIG_MAXVAR 
      ? (double) (vdimtot - best_dirct) / (double) max_variab
      : -0.1;
      ratio *= FABS(ratio);
    pref[Direct] = (int) POW((double) pref[Direct], 1.0 - ratio);
    if (pref[Direct] < PREF_BEST) 
      pref[Direct] = sub->finiterange ? PREF_BEST : PREF_BEST / 2;
  } 

  if (P0INT(GAUSSPROC_STATONLY) < 0 && isPosDef(cov)) {
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 6;
  } 

  if (!isCartesian(cov->isoown)) {
    pref[CircEmbedIntrinsic] = pref[CircEmbed] = pref[CircEmbedCutoff] =
				    LOC_PREF_NONE - 7;
    if (isAnySpherical(cov->isoown) && cov->xdimown < 3)
      pref[Sequential] = LOC_PREF_NONE - 8;
  }

   
  RU_orderingInt(pref, Nothing, 1, order);
  

}

bool NAequal(double X, double Y) {
  return ((ISNA(X) || ISNAN(X)) && (ISNA(Y) || ISNAN(Y))) || X == Y;
}

int kappaBoxCoxParam(cov_model *cov, int BC) {
  //  PMI(cov);
  //printf("kapap bc %s %d\n", NAME(cov), BC);
  int vdim_2 = cov->vdim[0] * 2;					
  if (PisNULL(BC)) {						
    PALLOC(BC, 2, cov->vdim[0]);				
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
      for (int i=0; i<cov->vdim[0]; i++) {				
	if ((notok = P(BC)[2*i] != RF_INF &&		
	     (P(BC)[2*i] != 0.0 || P(BC)[2*i+1] != 0.0)))
	  break;							
      }									
    } else {								
      for (int i=0; i<cov->vdim[0]; i++) {
	//printf("%d  %f %f   %f %f\n", i, GLOBAL.gauss.boxcox[2 * i], GLOBAL.gauss.boxcox[2 * i + 1], P(BC)[2 * i], P(BC)[2 * i +1]);
 	if ((notok = (GLOBAL.gauss.boxcox[2 * i] != RF_INF &&
		      !NAequal(P(BC)[2 * i], GLOBAL.gauss.boxcox[2 * i])) ||   
	     !NAequal(P(BC)[2 * i +1], GLOBAL.gauss.boxcox[2 * i + 1])))
	  break;							
      }									
    }									
    if (notok)								
      SERR("Box Cox transformation is given twice, locally and through 'RFoptions'"); 
  }									
  for (int i=0; i<cov->vdim[0]; i++) {				
    GLOBAL.gauss.boxcox[2 * i] = RF_INF;				
    GLOBAL.gauss.boxcox[2 * i + 1] = 0.0;				
  }									
  return NOERROR;
}


void kappaGProc(int i, cov_model *cov, int *nr, int *nc){
  *nc = i == GAUSS_BOXCOX ? SIZE_NOT_DETERMINED : 1;
  *nr = i == GAUSS_BOXCOX ? SIZE_NOT_DETERMINED 
    : i < CovList[cov->nr].kappas ? 1 : -1;
}

int checkgaussprocess(cov_model *cov) {
  cov_model 
    *next=cov->sub[cov->sub[0] == NULL],
    *key = cov->key;
  //location_type *loc=Loc(cov);
  int err, role,
    xdim = cov->xdimown, // could differ from tsdim in case of distances!
    dim = cov->tsdim;
  gauss_param *gp  = &(GLOBAL.gauss);  
  
  assert((Loc(cov)->distances && xdim==1) || xdim == dim);
    
  ROLE_ASSERT(ROLE_GAUSS || cov->role == ROLE_BERNOULLI || 
	      cov->role == ROLE_MAXSTABLE || cov->role == ROLE_LIKELIHOOD);
 						      
   kdefault(cov, GAUSSPROC_STATONLY,
	   gp->stationary_only >= 0 ? gp->stationary_only : -1);

  if (GLOBAL.direct.maxvariables < gp->direct_bestvariables) 
    SERR("maximum variables less than bestvariables for direct method");
  if ((err = checkkappas(cov, false)) != NOERROR) return err;

  cov->maxdim=INFDIM;
  role = isVariogram(next) ?  ROLE_COV 
    : isTrend(next) ? ROLE_GAUSS
    : isGaussMethod(next) ? ROLE_GAUSS 
    : ROLE_UNDEFINED;
  ASSERT_ROLE_DEFINED(next);
  
  if (key == NULL) {
    if (isGaussMethod(next)) {
      // wird verboten
      SERR1("%s may not call a method", NICK(cov));
    } else {
          
      if ((err = CHECKPD2ND(next, dim, xdim, 
			    SymmetricOf(cov->isoown),// Jan 2015 SYMMETRIC,
			    SUBMODEL_DEP, role))
	  != NOERROR) {

	if (CHECK(next, dim, dim, TrendType, XONLY, cov->isoown,  // err mocjt verwemdet
		  SUBMODEL_DEP, role)) return err; // previous error
      }
      
    }
    // next->delflag = DEL_COV - 21;
  } else {
    if (PL >= PL_COV_STRUCTURE) {
      PRINTF("checking key in gauss process  ...\n");
    }

    if ((err = CHECK(key, dim, xdim, ProcessType, XONLY, cov->isoown, 
		     SUBMODEL_DEP, 
		     cov->role == ROLE_BASE ? ROLE_BASE : ROLE_GAUSS))
	!= NOERROR) return err;
  }
  
  cov_model *sub = cov->key==NULL ? next : key;
  setbackward(cov, sub);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) return err;
  if ((err = checkkappas(cov, true)) != NOERROR) return err;
  return NOERROR;
}




void rangegaussprocess(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;
  range->min[GAUSSPROC_STATONLY] = -1;
  range->max[GAUSSPROC_STATONLY] = 1;
  range->pmin[GAUSSPROC_STATONLY] = -1;
  range->pmax[GAUSSPROC_STATONLY] = 1;
  range->openmin[GAUSSPROC_STATONLY] = false;
  range->openmax[GAUSSPROC_STATONLY] = false; 

}

int gauss_init_settings(cov_model *cov) {
  cov_model 
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
  
  if (false && cov->key != NULL) { // sind diese Zeilen notwendig?
     err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, KERNEL, SYMMETRIC,
		cov->vdim, ROLE_COV);
    if (err != NOERROR) goto ErrorHandling;
  }
  
  if (next->domprev == XONLY) {
    COV(ZERO, next, variance);
  } else {
    assert(next->domprev == KERNEL);
    NONSTATCOV(ZERO, ZERO, next, variance);
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
 
  cov->fieldreturn = true;
  cov->origrf = false;
  cov->rf = sub->rf;

 ErrorHandling:
  FREE(variance);
  FREE(mean);


  return err;
}


int struct_extractdollar(cov_model *cov, cov_model **newmodel) {
  // uebernimmt struct- und init-Aufgaben !!
  cov_model *next = cov->sub[0]; // nicht cov->sub[cov->sub != NULL]
  location_type *loc=Loc(cov);
  int role,
    nr = cov->nr,
    err = ERROROUTOFMETHODLIST,
    xdim = cov->xdimprev, // could differ from tsdim in case of distances!
    dim = cov->tsdim;
  assert((loc->distances && xdim==1) || xdim == dim);

  cov->fieldreturn = true;
  ASSERT_NEWMODEL_NULL;
  assert(cov!=NULL);
  ROLE_ASSERT_GAUSS;

  if ((cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) &&
      (!loc->distances || cov->xdimprev!=1)) {
    return ERRORDIM;
  }

  if (next!= NULL && !isVariogram(next)) {
    SERR("submodel not a covariance function");
  }
  
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if ((err = covCpy(&(cov->key), cov)) != NOERROR) return err; //!!cov
  
  cov->key->nr =
    nr == AVERAGE_USER       ? AVERAGE_INTERN :
    nr == CE_CUTOFFPROC_USER ? CE_CUTOFFPROC_INTERN : 
    nr == CE_INTRINPROC_USER ? CE_INTRINPROC_INTERN :
    nr == HYPERPLANE_USER    ? HYPERPLANE_INTERN :
    nr == NUGGET_USER        ? NUGGET_INTERN :
    nr == RANDOMCOIN_USER    ? AVERAGE_INTERN :
    nr == SPECTRAL_PROC_USER ? SPECTRAL_PROC_INTERN :
    nr == TBM_PROC_USER      ? TBM_PROC_INTERN : MISSING_COV;

  role = nr == AVERAGE_USER || nr == RANDOMCOIN_USER ?  ROLE_POISSON_GAUSS
    : nr == HYPERPLANE_USER ? ROLE_GAUSS : ROLE_GAUSS;
 
  // cov!! nicht cov->key!!  ??? Nein, sonst ist Coin u.ae. nicht gesetzt.
  if ((err = CHECK(cov, dim, xdim, GaussMethodType, cov->domown,
		   cov->isoown, cov->vdim, 
		   // role // 19.5.2013 auskommentiert und ersetzt durch
		   ROLE_BASE  // wegen spectral.cc
		   )) != NOERROR) {
     return err;
  }
   
  err = STRUCT(cov->key, NULL);
  cov->role = ROLE_GAUSS; // wichtig !! -> check_nugget_proc

  cov_model *dollar,
    *firstmethod = cov->key; // z.B. coins_intern
  int idx = firstmethod->sub[0] == NULL;
  cov_model
    *sub = firstmethod->sub[idx]; // z.B. oberste fkt shape von coins
  cov_model *lastmethod= isGaussMethod(sub->typus) ? sub : firstmethod;
  if (err != NOERROR && (err != ERRORPREFNONE || !isAnyDollar(sub))){
    
    return err;
  }
  
  if (err != NOERROR) {   
    dollar = cov->key = lastmethod->sub[idx]; // $ jetzt im cov->key
    lastmethod->sub[idx] = dollar->sub[0]; // eigentl. Modell jetzt anhaengen
    lastmethod->sub[idx]->calling = lastmethod;     
    dollar->sub[0] = firstmethod;       // auf $ setzen
    firstmethod->calling = dollar; 
    dollar->calling = cov;              // $ mit cov verbinden
    dollar->prevloc = cov->prevloc;


    /// cov nicht cov->key !!! ???? OK. Da sonst cov u.U. nicht gesetzt
    if ((err = CHECK(cov, dim, xdim, ProcessType, cov->domown, cov->isoown, 
		     cov->vdim, role)) != NOERROR) {
      return err;
    }

    if ((err = STRUCT(cov->key, NULL)) != NOERROR) { return err; }
  }

  Methods m;
  for (m=CircEmbed; m <= Hyperplane; m = (Methods) ((int) (m + 1))) {
    if (gaussmethod[(Methods) m] == cov->nr) break;
  }
  cov->key->method = m;
  
  return NOERROR;
}


  
int struct_gaussprocess(cov_model *cov, cov_model **newmodel) {
  // uebernimmt struct- und init-Aufgaben !!
  ASSERT_NEWMODEL_NULL;
  assert(cov!=NULL);
  if (cov->role == ROLE_LIKELIHOOD) {
    if (cov->vdim[0] > MAXGAUSSVDIM)
      SERR2("The number of variables, %d, exceeds the maximum number %d",
	    cov->vdim[0],  MAXGAUSSVDIM);
    return struct_gauss_logli(cov);
  } else if (cov->role!= ROLE_GAUSS) return ERRORFAILED;

  cov->fieldreturn = true;
  location_type *loc = Loc(cov);
  cov_model *next = cov->sub[cov->sub[0] == NULL];
  pref_type locpref, pref;
  int order[Nothing], i,
    err = ERROROUTOFMETHODLIST;
  char dummy[100 * Nothing];
  int 
    xdim = cov->xdimprev, // could differ from tsdgetim in case of distances!
    dim = cov->tsdim;
  assert((loc->distances && xdim==1) || xdim == dim);
#define  MAXFAILMSG 9
  static char FailureMsg[MAXFAILMSG][80] = 
    {"unknown reason",
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
    nr = next->nr;
#define nsel 4

  if( (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) &&
    (!loc->distances || cov->xdimprev!=1)) {
    return ERRORDIM;
  }

  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (!isVariogram(next) && !isTrend(next)) {
    SERR("submodel must be a covariance function");
  }
 

  //////////////////////////////////////////////////////////////////////
  // ROLE_GAUSS; cov->key not given

  { // only for error reporting !
    cov_model *sub = cov; 
    while(isAnyDollar(sub)) sub = sub->sub[0];
    SPRINTF(ERROR_LOC, "Searching a simulation method for '%s': ",
	    NICK(sub));
  }

  location_rules(cov, locpref); // hier cov 
  mixed_rules(cov, locpref, pref, order);

  if (PL >= PL_RECURSIVE) {
    cov_fct *N = CovList + next->nr;
    LPRINT("\n%s:%s\n", NAME(cov), NAME(next));
    for (i=0; i < Nothing; i++) {
      LPRINT("%-15s: base=%1d  covprev=%1d   locpref=%4d   totpref=%6d\n", 
	     METHODNAMES[i], N->pref[i], next->pref[i], 
	     (locpref[i] < -9000) ? -999 : locpref[i], pref[i]);
     }  
    LPRINT("initstandard %s (%d) [pref value %d]\n", 
	   CovList[nr].name, nr, pref[order[Nothing -1]]);
  }


  all_PREF_NONE = true;
  for (i=0; i<Nothing; i++) all_PREF_NONE &= pref[i] == PREF_NONE;
  
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if ((err = covCpy(&(cov->key), next)) != NOERROR) {
    return err;
  }

  if (all_PREF_NONE) {
    err = pref[Nothing] == PREF_NONE ? ERRORINVALIDMODEL : ERRORODDMODEL;
    goto ErrorHandling;
  } 

  err = ERROROUTOFMETHODLIST; // in case none of the methods are available 
  for (i = Nothing - 1; i>=0 && pref[order[i]] > PREF_NONE; i--) {  

    if (PL >= PL_BRANCHING && i < Nothing - 1) {
      LPRINT("'%s%s' failed for '%s'\n", 
	     CAT_TYPENAMES[GaussMethodType], METHODNAMES[unimeth], NICK(next)); 
    }

    zaehler++;
    unimeth = (Methods) order[i];

    if (PL >= PL_BRANCHING) {
      LPRINT("trying '%s%s' for '%s'\n", 
	     CAT_TYPENAMES[GaussMethodType], METHODNAMES[unimeth], NICK(next));
    }
    
    meth = gaussmethod[unimeth];      
  
    addModel(&(cov->key), meth);
 
    cov->key->calling = cov;
    
    cov_model *key = cov->key;
    key->prevloc = PLoc(cov);
    int role = cov->role != ROLE_GAUSS || key->nr != AVERAGE_INTERN 
      ? cov->role
      : ROLE_POISSON_GAUSS;

    if ((err = CHECK(key, dim, xdim, GaussMethodType, cov->domown,
		     cov->isoown, cov->vdim, role)) == NOERROR) {
    
      if ((err = STRUCT(key, NULL)) <= NOERROR) {
	if (!key->initialised) {
	  if ((err = CHECK(key, dim, xdim, GaussMethodType, cov->domown,
			   cov->isoown, cov->vdim, role)) == NOERROR){
	    key->method = unimeth;
	    
	    NEW_STORAGE(gen);
	    
	    err = INIT(key, role == ROLE_POISSON_GAUSS ? 2 : 0, cov->Sgen);
	    if (err == NOERROR) {
	      if (PL >= PL_RECURSIVE) {
		LPRINT("returning to higher level...\n");
	      }
	      goto Initialized;
	    }
	  }
	}
      }
     
      if (PL >= PL_RECURSIVE) {  
	LPRINT("back from %s: err =%d [%d, %s]\n",
	       METHODNAMES[unimeth], err, unimeth, METHODNAMES[unimeth]);     
	
	if (PL > PL_ERRORS) {
	  char msg[LENERRMSG]; FinalErrorMSG(err, msg); PRINTF("error in init: %s\n",msg);
	}
      } else if (PL >= PL_BRANCHING) MERR(err) //
    } else {
      if (PL > PL_ERRORS) {
	char msg[LENERRMSG]; FinalErrorMSG(err,msg); PRINTF("error in check: %s\n",msg);
      }
    }
    key = key->sub[0];
  
    COV_DELETE_WITHOUTSUB(&(cov->key));
  
    cov->key=key;
    key->calling = cov;      
  }
  
  if (PL >= PL_RECURSIVE) {
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
	SPRINTF(lpd, "%s (locp) ", FailureMsg[MAX(0, MIN(MAXFAILMSG-1, 1-lp))]);
      }
      if (p>0 || p == lp) {
	if (lp>0) STRCPY(pd, "(method specific failure)");
	else STRCPY(pd, ""); // STRCPY(pd, "(method specific failure)");
      } else {
	SPRINTF(pd, "%s (pref)", FailureMsg[MAX(0, MIN(MAXFAILMSG-1, 1-p))]);
      }
      strcopyN(names, METHODNAMES[i], NMAX-1);
      SPRINTF(dummy, "%s %-13s: %s%s\n", PREF_FAILURE, names, lpd, pd);
      STRCPY(PREF_FAILURE, dummy);
    }
  }
   
  if (err != NOERROR) goto ErrorHandling;
  
 Initialized:
  // nachfolgende Zeile muss stehen, da Init in struct aufgerufen!
  err = gauss_init_settings(cov);

  cov->initialised = err == NOERROR;	  

 ErrorHandling:
  /*
  if (false)
  if (err > 0 && ISNA(GLOBAL.general.exactness)) {
    if (*newmodel != NULL) COV_DELETE(newmodel);
    GLOBAL.general.exactness = FALSE;
    int newerr = struct_gaussprocess(cov, newmodel);
    GLOBAL.general.exactness = RF_NA;
    if (newerr == NOERROR) {
      char msg[LENERRMSG];
      SPRINTF(msg, "A somehow rougher approximation has been used. Set '%s=FALSE' to avoid this warning. See 'RFoptions' for more details.", general[GENERAL_EXACTNESS]);
      warning(msg);
      return NOERROR;
    } else return err;
  }
  */

  if (err <= NOERROR || zaehler==0){
    if (err <= NOERROR) {
      if (cov->key == NULL) next->simu.active = true;
      else cov->key->simu.active = true;
    }
    return err;
  }

  if (zaehler==1) {
    SPRINTF(ERROR_LOC, "Only 1 method found for '%s', namely '%s', that comes into question.\nHowever", NICK(next), METHODNAMES[unimeth]);
    return err;
  }
  
  { // only for error reporting !
    cov_model *sub = next; 
    while(isAnyDollar(sub)) sub = sub->sub[0];
    SPRINTF(ERROR_LOC, "searching a simulation method for '%s': ",
	    NICK(sub));
  }

  return ERROROUTOFMETHODLIST;
}


int init_gaussprocess(cov_model *cov, gen_storage *s) {
  /// ACHTUNG struct_gaussprocess hat bereits init aufgerufen!

  // direkter Aufruf, u.a. durch RPtbm, average, cutoff, intrinsic, hyperplane,
  //                    nugget, spectral !!!

  cov_model 
    *key = cov->key;
  int err;

  ROLE_ASSERT_GAUSS;
  assert(key != NULL);


  if ((err = INIT(key, 0, s)) != NOERROR) return err;
 	    
 if ((err = gauss_init_settings(cov)) != NOERROR) {
   return err;
  }
  
  key->simu.active = true;

  return NOERROR;
}


   

void do_gaussprocess(cov_model *cov, gen_storage *s) {
  assert(s != NULL);
  // reopened by internal_dogauss
  errorloc_type errorloc_save;
  double *res = cov->rf;
  int i,
    vdimtot = Gettotalpoints(cov) * cov->vdim[0] ;
  cov_model *key = cov->key ;
  SAVE_GAUSS_TRAFO;
 
  STRCPY( errorloc_save,  ERROR_LOC);

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
  STRCPY( ERROR_LOC, errorloc_save);
}




//////////////////////////////////////////////////////////////////////
// Binary Gauss
#define BINARY_THRESHOLD (GAUSSPROC_LAST + 1)

void kappa_binaryprocess(int i, cov_model VARIABLE_IS_NOT_USED *cov, 
			 int *nr, int *nc){
  kappaGProc(i, cov, nr, nc);
  if (i == BINARY_THRESHOLD) *nr = SIZE_NOT_DETERMINED;
}

int checkbinaryprocess(cov_model *cov) {
  // Achtung! binary hat Doppelbedeutung: binary-Gaussian als default
  // und binary angewandt auf x-bel. process.
  cov_model
    *key = cov->key,
    *next = cov->sub[0],
    *sub = key != NULL ? key : next;
  double v;
  int role,
    err = NOERROR;
  if (PisNULL(BINARY_THRESHOLD)) kdefault(cov, BINARY_THRESHOLD, 0);

  // PMI(cov); printf("%d %d\n", key==NULL, isVariogram(next));
    
  if (key == NULL && isVariogram(next)) {
    if ((err = checkgaussprocess(cov)) != NOERROR) {
      return err;
    }
    
    COV(ZERO, sub, &v);
    if (v != 1.0) 
      SERR("binaryian requires a correlation function as submodel.");

  } else if (isProcess(sub)) {
     role = cov->role == ROLE_BASE ? cov->role : role_of_process(sub->nr);
    
    if ((err = CHECK(sub, cov->tsdim, cov->xdimprev, ProcessType, 
		     cov->domown, cov->isoown,
		     SUBMODEL_DEP, 
		     role)) != NOERROR) return err;
    setbackward(cov, sub);
  } else SERR1("process type model required, but '%s' obtained", NICK(sub));

  cov->vdim[0] = sub->vdim[0];
  cov->vdim[1] = sub->vdim[1];
  return NOERROR;
}

int struct_binaryprocess(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  cov_model
    *next = cov->sub[0];
  int err;

  ROLE_ASSERT(ROLE_BERNOULLI);
  if (isVariogram(next)) {
    assert(cov->key == NULL);
    err = covCpy(&(cov->key), cov);

    // darauffolgende Zeile absichern:
    if (CovList[cov->nr].kappas != 3 || CovList[GAUSSPROC].kappas != 2) BUG;
    if (cov->key != NULL && !PARAMisNULL(cov->key, BINARY_THRESHOLD)) {
      PARAMFREE(cov->key, BINARY_THRESHOLD);
      PARAMtoNULL(cov->key, BINARY_THRESHOLD);
    }

    if (err != NOERROR)  return err; //!!cov

    cov->key->nr = GAUSSPROC;
    err = CHECK(cov->key, cov->tsdim, cov->xdimprev, ProcessType, 
		cov->domown, cov->isoown, SUBMODEL_DEP, ROLE_GAUSS);
    if (err != NOERROR)  return err;
    err = STRUCT(cov->key, NULL);
    return err;
  }
  else return STRUCT(next, NULL);
}

int init_binaryprocess( cov_model *cov, gen_storage *s) {
  double sigma, 
    *mean = NULL, 
    *variance = NULL,
    *p = P(BINARY_THRESHOLD);
  cov_model 
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
 
  if (isVariogram(next) || next->nr == GAUSSPROC) {
    GetInternalMean(next, vdim, mean);
    if (ISNAN(mean[0])) // GetInternalMean currently only allows ...
      GERR1("'%s' currently only allows scalar fields - NA returned",
	    NICK(cov));
     if (cov->mpp.moments >= 1) 
      COV(ZERO, next->nr==GAUSSPROC ? next->sub[0]: next, variance);
    nmP1 = cov->mpp.moments + 1;
    for (pi=v=w=0; w<vdimSq; w+=vdimP1, v++, pi = (pi + 1) % npi ) { 
      int idx = v * nmP1;
      cov->mpp.maxheights[v] = 1.0;
      if (cov->mpp.moments >= 0) { 
	assert(cov->mpp.mM != NULL);
	cov->mpp.mM[idx + 0] = cov->mpp.mMplus[idx + 0] = 1.0; 
	if (cov->mpp.moments >= 1) {
	  if (variance[w]==0.0)
	    GERR1("Vanishing sill not allowed in '%s'", NICK(next));
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

  cov->fieldreturn = true;
  cov->simu.active = err == NOERROR;

 ErrorHandling:
  FREE(variance);
  FREE(mean);
 
  return err;
}



void do_binaryprocess(cov_model *cov, gen_storage *s){
  // reopened by internal_dogauss
  long j,
    tot = Loc(cov)->totalpoints,
    endfor = tot;
  assert(cov->ownloc == NULL);
  int  i,
    pi, 
    npi = cov->nrow[BINARY_THRESHOLD],
    vdim = cov->vdim[0];
  double 
    threshold,
    *p = P(BINARY_THRESHOLD),
    *rf = cov->rf;
  cov_model 
    *next=cov->sub[0];
  assert(cov->vdim[0] == cov->vdim[1]);
  if (isVariogram(next)) {
    do_gaussprocess(cov, s);
  } else DO(next, s);

  for (j=pi=i=0; i<vdim; i++, pi = (pi + 1) % npi, endfor += tot) {
    threshold = p[pi];
    // printf("endfor %d %f\n", endfor, threshold);
    if (threshold > RF_NEGINF && threshold < RF_INF)
      for ( ; j<endfor; j++) rf[j] = (double) (rf[j] >= threshold);
  }
}



void rangebinaryprocess(cov_model *cov, range_type *range) { 
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


int checkchisqprocess(cov_model *cov) {
  cov_model
    *key = cov->key,
    *next = cov->sub[0];
  double *v = NULL;
  int err = NOERROR,
    xdim = cov->xdimprev, // could differ from tsdim in case of distances!
    dim = cov->tsdim;
  assert((Loc(cov)->distances && xdim==1) || xdim == dim);

  if (PisNULL(CHISQ_DEGREE)) SERR("degree of freedom must be given");
  if (key == NULL) {
    if (!isGaussProcess(next) && !isVariogram(next))
      SERR1("Gaussian process required, but '%s' obtained", NICK(cov));
    if ((err = CHECK(next, dim, xdim, ProcessType, XONLY, cov->isoown, 
		     SUBMODEL_DEP, cov->role)) != NOERROR) {
      if ((err = CHECK(next, dim, xdim, 
		       isCartesian(cov->isoown) ? VariogramType : PosDefType, 
		       KERNEL, 
		       SymmetricOf(cov->isoown),
		       SUBMODEL_DEP, ROLE_COV)) != NOERROR) {
      return err;
      }
    }
    
    int 
      vdim = next->vdim[0],
      vdimSq = vdim * vdim;
    assert(vdim > 0);

    if ((v = (double*) MALLOC(sizeof(double) * vdimSq)) == NULL)
      return ERRORMEMORYALLOCATION;

    cov_model *sub = next;
    while (sub != NULL && isProcess(sub)) sub=sub->sub[0];
    if (sub == NULL) BUG;
    if (sub->domprev == XONLY) {
      COV(ZERO, next, v);
    } else {
      assert(sub->domprev == KERNEL);
      NONSTATCOV(ZERO, ZERO, next, v);
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
    
    cov->vdim[0] = sub->vdim[0];
    cov->vdim[1] = sub->vdim[1];
     //no need to setbackward
  } else {
    if ((err = CHECK(key, cov->tsdim,  cov->xdimprev, ProcessType,
		       cov->domown, cov->isoown,
		       SUBMODEL_DEP, cov->role)) != NOERROR) return err;
    setbackward(cov, key);
  }

  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) return err;
 
  return NOERROR;
}
 
int struct_chisqprocess(cov_model *cov, 
			cov_model VARIABLE_IS_NOT_USED **newmodel) {
  cov_model
    *next = cov->sub[0];
  int err;

  ROLE_ASSERT_GAUSS;

  if (isVariogram(next)) {
    assert(cov->key == NULL);
    if ((err = covCpy(&(cov->key), next)) > NOERROR) return err;
    addModel(&(cov->key), GAUSSPROC);
    cov->key->calling = cov;
    if ((err = CHECK(cov->key, cov->tsdim,  cov->tsdim, ProcessType,
		     cov->domprev, cov->isoprev,
		     SUBMODEL_DEP, ROLE_GAUSS)) != NOERROR) return err;
   return STRUCT(cov->key, NULL);
  } else {
    return STRUCT(next, NULL);
  }
}

int init_chisqprocess(cov_model *cov, gen_storage *s) {
  double // sigma,
    mean, m2, 
    variance;
  cov_model 
    *key = cov->key,
    *next = cov->sub[0],
    *sub = key == NULL ? next : key;
  int i, err, nmP1,
    nmP1sub = sub->mpp.moments + 1,
    vdim = cov->vdim[0];

  cov->simu.active = false;
  if ((err = INIT(sub, 
		  CovList[cov->nr].range == rangechisqprocess ? 2 :
		  CovList[cov->nr].range == rangetprocess ? 1 : 9999
		  , s)) != NOERROR) return err;		  
 
  nmP1 = cov->mpp.moments + 1;
  for (i=0; i<vdim; i++) { 
    int 
      idx = i * nmP1,
      idxsub = i * nmP1sub;
    mean = sub->mpp.mM[idxsub + 1];
    m2 = sub->mpp.mM[idxsub + 2];
    variance =  m2 - mean * mean; 
    if (variance==0.0) 
      SERR1("Vanishing sill not allowed in '%s'", NICK(sub));
    if (ISNAN(mean)) // GetInternalMean currently only allows ...
      SERR1("'%s' currently only allows scalar fields -- NA returned", 
	   NICK(cov));
    cov->mpp.maxheights[i] = GLOBAL.extreme.standardmax * 
      GLOBAL.extreme.standardmax * m2;

    if (cov->mpp.moments >= 0) { 
      assert(cov->mpp.mM != NULL);
      cov->mpp.mM[idx + 0] = cov->mpp.mMplus[idx + 0] = 1.0; 
      if (cov->mpp.moments >= 1) {
	cov->mpp.mMplus[idx + 1] = 		  
	  CovList[cov->nr].range == rangechisqprocess ? m2 
	  : CovList[cov->nr].range == rangetprocess ? RF_NAN
	  : RF_NAN;
	cov->mpp.mM[idx + 1] = RF_NA; // to do: richtiger Wert
	if (cov->mpp.moments >= 2)
	  cov->mpp.mM[idx + 2] = 3 * variance * RF_NA;//check the 4th moment of a Gaussian, to do
      }
    }
  }


  if (CovList[cov->nr].range == rangechisqprocess) FieldReturn(cov);
  else if (CovList[cov->nr].range == rangetprocess){
    cov->fieldreturn = true;
    cov->origrf = false;
    cov->rf = sub->rf;
  } else BUG;
  cov->simu.active = true;
 
  return NOERROR;
}



void do_chisqprocess(cov_model *cov, gen_storage *s){
  // reopened by internal_dogauss
  int  i, f,
    degree = P0INT(CHISQ_DEGREE),
    totvdim = Loc(cov)->totalpoints * cov->vdim[0];
 assert(cov->ownloc == NULL);
  cov_model 
    *next=cov->sub[0],
    *key = cov->key,
    *sub = key == NULL ? next : key;
  double 
    *subrf = sub->rf,
    *res = cov->rf;

  for (i=0; i<totvdim; i++) res[i] = 0.0;

  //  APMI(cov);

  for (f=0; f<degree; f++) {
    //  printf("f=%d %d\n", f, degree);
    DO(sub, s);   
    for (i=0; i<totvdim; i++) {
      double dummy = subrf[i];
      res[i] += dummy * dummy;
    }
  }
  BOXCOX_INVERSE;
}

void rangechisqprocess(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
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





void do_tprocess(cov_model *cov, gen_storage *s){
  // reopened by internal_dogauss
  int  i,
    n = Loc(cov)->totalpoints * cov->vdim[0];
  cov_model 
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

void rangetprocess(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;
 
  range->min[T_NU] = 0;
  range->max[T_NU] = RF_INF;
  range->pmin[T_NU] = 1;
  range->pmax[T_NU] = 1000;
  range->openmin[T_NU] = true;
  range->openmax[T_NU] = true;
}

