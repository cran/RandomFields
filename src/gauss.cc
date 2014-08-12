/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

 Copyright (C) 2005 -- 2014 Martin Schlather

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

#include <math.h>
#include "RF.h"
#include "Covariance.h"



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
     Direct, Specific, Sequential, Markov, Average, Nugget, RandomCoin,
     Hyperplane
  };
  for (i=0; i<Nothing; i++) {
    pref[Standard[i]] = Nothing - i; 
  }

//  if (loc->totalpoints > max_dirct) pref[Direct] = LOC_PREF_NONE - 0;
  
  // APMI(cov);
  

  if (P0INT(GAUSSPROC_STATONLY) > 0) {
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 1;
  }

  if (!ISNA(exactness) && exactness) {
    pref[TBM] = pref[SpectralTBM] = pref[Average] = pref[RandomCoin] = 
      pref[Markov] = pref[Sequential] = pref[Hyperplane] 
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
      pref[CircEmbedIntrinsic] = -3; // not programmed yet
      pref[CircEmbedCutoff] -= PREF_PENALTY;
    } else {
      pref[CircEmbed] = pref[CircEmbedIntrinsic] = pref[CircEmbedCutoff] = -3;
    }
    pref[Markov] = LOC_PREF_NONE - 3;
    if (!loc->Time) pref[Sequential] = LOC_PREF_NONE;
  }
  //  print("%d\n", pref[TBM3]);
}


void mixed_rules(cov_model *cov, pref_type locpref, 
		 pref_type pref, int *order) {
  
  assert(cov->nr == GAUSSPROC);
  
  // printf("mixedrules %ld\n\n", (long int) cov); PMI(cov);
  //printf("back to mixed\n");
  
  cov_model *sub = cov->sub[0];
  location_type *loc = Loc(cov);
  //  double exactness = GLOBAL.general.exactness; //
  
  int *covpref = sub->pref;
  //  cov_fct *C = CovList + cov->nr;// nicht gatternr
  // printf("herab\n");
  int i,
    totalpref[Nothing],
    best_dirct=GLOBAL.gauss.direct_bestvariables,
    max_dirct= GLOBAL.direct.maxvariables,
    vdim = cov->vdim2[0];
  
  //printf("hera %s\n", NICK(sub));
  //for(i=0; i<Nothing; i++) printf("tot i=%d %d\n", i, CovList[TREND].pref[i]);

  for (i=0; i<Nothing; i++) {
     totalpref[i] =  i == Nugget ? PREF_NUGGET * (PREF_FACTOR + 1): PREF_BEST;
     //printf("totalpref %d ----- ", totalpref[i]);
    if (totalpref[i] > covpref[i]) {
       totalpref[i] = covpref[i];
     }
     // printf("i=%s %d\n", METHODNAMES[i], covpref[i]);
   
     // print("ok\n"); printf("loc i=%d %d\n", i, locpref[i]);
     // 
    //printf("tot i=%d %d (%d, %d)\n", i, totalpref[i],
    //CovList[sub->nr].pref[i], CovList[TREND].pref[i]);

     pref[i] = (totalpref[i] > PREF_NONE && locpref[i] > LOC_PREF_NONE)
       ? totalpref[i] * PREF_FACTOR + locpref[i]
       : locpref[i] <= LOC_PREF_NONE ? locpref[i] :  LOC_PREF_NONE - 4;
   }
  //  PLE(sub);

   //   printf("\n\n>>> her %d\n\n\n", loc->totalpoints ); 
 
  if (loc->totalpoints * vdim > max_dirct) pref[Direct] = LOC_PREF_NONE-0;

  if (loc->totalpoints * vdim <= best_dirct && 
      totalpref[Direct] == PREF_BEST)
    pref[Direct] = (PREF_BEST + 1) * PREF_FACTOR;
  

  //  if ((ISNA(exactness) || !exactness) && C->tbm2 == NULL)
  //    pref[TBM2] = LOC_PREF_NONE - 5;

  if (P0INT(GAUSSPROC_STATONLY) < 0 && isPosDef(cov)) {
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 6;
  }
    
  orderingInt(pref, Nothing, 1, order);


}





//////////////////////////////////////////////////////////////////////
// Gauss process


void gaussprocessDlog(double VARIABLE_IS_NOT_USED  *x, cov_model *cov, double VARIABLE_IS_NOT_USED *v) {
  BUG;
  // search true covariance function
  location_type *loc = Loc(cov);
  int 
    vdim = cov->vdim2[0];
  long
    totalpoints = loc->totalpoints,
    total = vdim * totalpoints,
    totalsq = total * total;

  ALLOC_EXTRA(M, totalsq);
  invertMatrix(M, total);  
}


int checkgaussprocess(cov_model *cov) {
  cov_model 
    *next=cov->sub[cov->sub[0] == NULL],
    *key = cov->key;
  location_type *loc=Loc(cov);
  int err, role,
    //    domown = cov->domown,
    // taken[MAX DIM],
    xdim = cov->xdimprev, // could differ from tsdim in case of distances!
    dim = cov->tsdim;
  gauss_param *gp  = &(GLOBAL.gauss);  
  assert((loc->distances && xdim==1) || xdim == dim);
     
  ROLE_ASSERT(ROLE_GAUSS || cov->role == ROLE_BERNOULLI || 
	      cov->role == ROLE_MAXSTABLE);
 
 						      
  kdefault(cov, GAUSSPROC_STATONLY,
	   gp->stationary_only >= 0 ? gp->stationary_only : -1);

  //  printf("direct %d %d\n", GLOBAL.direct.maxvariables,gp->direct_bestvariables);
  if (GLOBAL.direct.maxvariables < gp->direct_bestvariables) 
    SERR("maximum variables less than bestvariables for direct method");
  if ((err = checkkappas(cov, true)) != NOERROR) return err;

  //if (loc->distances) {
  //    //ssert(false);
  //    if (domown != STATIONARY || cov->isoown != ISOTROPIC)
  //      return ERRORINCORRECTSTATISO;
  //  } else {
  //    if ((domown != COVARIANCE && domown != STATIONARY && domown != VARIOGRAM) 
  //	|| cov->isoown != ANISOTROPIC) {
  //  //     printf(">>>>>> %d %d %d\n", cov->isoown, ANISOTROPIC, domown);
  //  //      crash(cov);
  //  return  ERRORINCORRECTSTATISO;
  //  }
  //}

  if ((cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) &&
      (!loc->distances || cov->xdimprev!=1)) {
    //    
    //PMI(cov->calling);
    return ERRORDIM;
  }

  cov->maxdim=INFDIM;
  role = isNegDef(next) ?  ROLE_COV 
    : isTrend(next) ? ROLE_GAUSS
    : isGaussMethod(next) ? ROLE_GAUSS 
    : ROLE_UNDEFINED;
  ASSERT_ROLE_DEFINED(next);
  
  if (key == NULL) {
    if (isGaussMethod(next)) {
      // wird verboten
      SERR1("%s may not call a method", NICK(cov));
    } else {
      
      // PMI(cov);
      //  assert(role == ROLE_COV); // oder trend
      
      if ((err = CHECKPD2ND(next, dim, xdim, SYMMETRIC, SUBMODEL_DEP, role))
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

    // PMI(cov);    printf("gauss %d %s\n", dim, NICK(key));
    if ((err = CHECK(key, dim, xdim, ProcessType, XONLY, cov->isoown, 
		     SUBMODEL_DEP, 
		     cov->role == ROLE_BASE ? ROLE_BASE : ROLE_GAUSS))
	!= NOERROR) return err;
  }
  
  cov_model *sub = cov->key==NULL ? next : key;
  setbackward(cov, sub);

  EXTRA_STORAGE;
 
  // print("\n\n\nEND Check GAUSS\n\n\n");

  return NOERROR;
}




void rangegaussprocess(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range) {
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
    vdim = next->vdim2[0],
    vdimSq = vdim * vdim,
    vdimP1 = vdim + 1;
   //  printf("ok3\n");

  if ((variance = (double*) MALLOC(sizeof(double) * vdimSq)) == NULL ||
      (mean = (double*) CALLOC(vdim, sizeof(double))) == NULL
      ) goto ErrorHandling;

  GetInternalMean(next, vdim, mean);
  if (ISNAN(mean[0])) // GetInternalMean currently only allows ...
    GERR("Mean equals NA. Did you try a multivariate model or an incomplete (mixed) model?");
  
  //  printf("ok4\n");  APMI(next);
 
  
  // PMI(next, "gaus_init_XXXX");
  if (false && cov->key != NULL) { // sind diese Zeilen notwendig?
    // APMI(next);
    err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, KERNEL, SYMMETRIC,
		cov->vdim2, ROLE_COV);
    //XERR(err);
    if (err != NOERROR) goto ErrorHandling;
  }
  
  //    PMI(next, "gaus_init"); printf("cov->key=%ld", (long int) cov->key);//assert(false);
  if (next->domprev == XONLY) {
    COV(ZERO, next, variance);
  } else {
    assert(next->domprev == KERNEL);
    NONSTATCOV(ZERO, ZERO, next, variance);
  }
  
  //    printf("end nonstatcov in gauss.cc\n");
  if (cov->q == NULL) cov->q = (double*) MALLOC(sizeof(double) * vdim);
  cov->qlen = vdim;
  if ((err = alloc_mpp_M(cov, 2)) != NOERROR) goto ErrorHandling;
 
  nmP1 = cov->mpp.moments + 1;
   for (v=w=0; w<vdimSq; w+=vdimP1, v++) { 
    sigma = sqrt(variance[w]);
    meanDsigma = sigma == 0 ? RF_INF : mean[v] / sigma;
    
    //the following line is sqrt(2 pi s^2) int_0^infty x exp(-0.5 (x-mu)^2/s^2)
    var = sigma * INVSQRTTWOPI * exp(-0.5 * meanDsigma * meanDsigma) +
      mean[v] * /* correct next lines the 2nd factor is given */
      // pnorm(meanDsigma, key->mean, sigma, 1, 0));
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

    //printf("%d %d %d %d\n", v, nmP1, idx, cov->mpp.moments);

    cov->mpp.mM[idx + 0] = cov->mpp.mMplus[idx + 0] = 1.0;       
    cov->mpp.mMplus[idx + 1] = 
      sigma * INVSQRTTWOPI * exp(-0.5 * mean[v] * mean[v]) 
      + mean[v] * pnorm(-mean[v], 0.0, 1.0, false, false);
    //(2pi)^-0.5 int x exp(-(x-m)/2s^2)
    cov->mpp.mM[idx + 1] = 0.0;
    cov->mpp.mM[idx + 2] = variance[w];
    // todo: cov->mpp.mMplus[2] berechnen
  }
 
  cov->fieldreturn = true;
  cov->origrf = false;
  cov->rf = sub->rf;

  // printf("ok5\n");
 ErrorHandling:
  if (variance != NULL) free(variance);
  if (mean != NULL) free(mean);


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
  // print("extract dollar A\n");
  ROLE_ASSERT_GAUSS;

  if ((cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) &&
      (!loc->distances || cov->xdimprev!=1)) {
    return ERRORDIM;
  }

  //  PMI(cov, "xxxx");
  //if (cov->role == ROLE_COV) return NOERROR;  else

  
  //printf("covnr = %d %d\n", cov->nr, RANDOMCOIN_USER);
  //printf("covnr = %ld\n",  cov->sub[COIN_SHAPE]); assert(false);

  if (next!= NULL && !isNegDef(next)) {
    SERR("submodel not a covariance function");
  }
  
  //  PMI(cov, "before");
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if ((err = covcpy(&(cov->key), cov)) != NOERROR) return err; //!!cov

  cov->key->nr =
    nr == AVERAGE_USER       ? AVERAGE_INTERN :
    nr == CE_CUTOFFPROC_USER ? CE_CUTOFFPROC_INTERN : 
    nr == CE_INTRINPROC_USER ? CE_INTRINPROC_INTERN :
    nr == HYPERPLANE_USER    ? HYPERPLANE_INTERN :
    nr == NUGGET_USER        ? NUGGET_INTERN :
    nr == RANDOMCOIN_USER    ? AVERAGE_INTERN :
    nr == SPECTRAL_PROC_USER ? SPECTRAL_PROC_INTERN :
    nr == TBM_PROC_USER      ? TBM_PROC_INTERN : MISSING_COV;
  role = nr == AVERAGE_USER || nr == RANDOMCOIN_USER ? 
    ROLE_POISSON_GAUSS : ROLE_GAUSS;

  // printf("\n extractdollar %d\n", role);

  // print("key %d nr=%d %s %d %d\n\n\n\n\n\n", role, nr, 
  //NICK(cov),AVERAGE_USER, RANDOMCOIN_USER);   
  //APMI(cov);
  // assert(false);
  
  // cov!! nicht cov->key!!  ??? Nein, sonst ist Coin u.ae. nicht gesetzt.
  if ((err = CHECK(cov, dim, xdim, GaussMethodType, cov->domown,
		   cov->isoown, cov->vdim2, 
		   // role // 19.5.2013 auskommentiert und ersetzt durch
		   ROLE_BASE  // wegen spectral.cc
		   )) != NOERROR) {
    // APMI(cov);//, "YUUUU");
    // crash(cov);
    // print("extract dollar A %d", err);
    //  XERR(err);
    return err;
  }
  //print("\n\n\n\n\n err YY  %d \n", err);
  //  PMI(cov, "YYY");

  //XERR(err);
  // print("extract dollar C\n");
  
  err = STRUCT(cov->key, NULL);
  cov->role = ROLE_GAUSS; // wichtig !! -> check_nugget_proc

  // print("extract dollar D %d\n", err);
  //   
  // APMI(cov);//, "YUUUU");

  cov_model *dollar,
    *firstmethod = cov->key; // z.B. coins_intern
  int idx = firstmethod->sub[0] == NULL;
  cov_model
    *sub = firstmethod->sub[idx]; // z.B. oberste fkt shape von coins
  cov_model *lastmethod= isGaussMethod(sub->typus) ? sub : firstmethod;
  if (err != NOERROR && (err != ERRORPREFNONE || !isAnyDollar(sub))){
    //print("extract dollar C%d\n", err);
    
    return err;
  }
  //   XERR(err);


  // print("extract dollar\n");


  // RO, Nugget, S, Nuggetintern, RMnugget
 
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
		     cov->vdim2, role)) != NOERROR) {
      return err;
    }
         
    //     print("\n\n\n\n\n\n");   PMI(cov, "extract"); //assert(false);

    if ((err = STRUCT(cov->key, NULL)) != NOERROR) { return err; }
    // print("end extract\n");    
  }

  // APMI(cov);
  // assert(false);

  Methods m;
  for (m=CircEmbed; m <= Hyperplane; m = (Methods) ((int) (m + 1))) {
    if (gaussmethod[(Methods) m] == cov->nr) break;
  }
  cov->key->method = m;
  
  //APMI(cov->calling->calling);
  
  return NOERROR;
}

  
int struct_gaussprocess(cov_model *cov, cov_model **newmodel) {
  // uebernimmt struct- und init-Aufgaben !!


  ASSERT_NEWMODEL_NULL;
  assert(cov!=NULL);
  
  //  cov_model *alt = cov;
  //print("struct_gaussprocess %s\n", NICK(cov));

  cov->fieldreturn = true;
  location_type *loc = Loc(cov);
  cov_model 
    *next = cov->sub[cov->sub[0] == NULL];
 
  pref_type locpref, pref;
  int order[Nothing], i,
    err = ERROROUTOFMETHODLIST;
  char dummy[100 * Nothing];
  int 
    xdim = cov->xdimprev, // could differ from tsdim in case of distances!
    dim = cov->tsdim;
  assert((loc->distances && xdim==1) || xdim == dim);
  static char FailureMsg[][80] = 
    {"total number of points > maximum value",
     "non-domain model not allowed",
     "not an exact method as required by user",
     "locations not on a grid",
     "denied by cov. model or by user's pref.",
     "no analytic solutn of Abelian integral",
     "intrinsic fields denied by user"
    };
  bool all_PREF_NONE;
  Methods unimeth = Forbidden;
  int 
    meth, 
    zaehler = 0,
    nr = next->nr;
#define nsel 4
  // int statselect[nsel]={STATIONARY, VARIOGRAM, COVARIANCE, GEN_VARIOGRAM};


  if (cov->role != ROLE_GAUSS){
    //     APMI(cov);
    return ERRORFAILED;
  }

  if( (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) &&
    (!loc->distances || cov->xdimprev!=1)) {
    return ERRORDIM;
  }

  if (cov->key != NULL) COV_DELETE(&(cov->key));

  //  printf("ok7\n");
  //  if (cov->role == ROLE_COV) return NOERROR; else 
  if (!isNegDef(next) && !isTrend(next)) {
    SERR("submodel must be a covariance function");
  }
 

  //////////////////////////////////////////////////////////////////////
  // ROLE_GAUSS; cov->key not given

  //  print("-------------gauspr\n");  PMI(cov);

  { // only for error reporting !
    cov_model *sub = cov; 
    while(isAnyDollar(sub)) sub = sub->sub[0];
    sprintf(ERROR_LOC, "Searching a simulation method for '%s': ",
	    NICK(sub));
  }

 
  location_rules(cov, locpref); // hier cov 
  mixed_rules(cov, locpref, pref, order);

  //    PMI(cov->calling->calling);
  //   print("PL: %d\n", PL); 

  //PMI(cov);
  if (PL >= PL_REC_DETAILS) {
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
  if ((err = covcpy(&(cov->key), next)) != NOERROR) return err;

  if (all_PREF_NONE) {
    //PMI(cov);
    err = pref[Nothing] == PREF_NONE ? ERRORINVALIDMODEL : ERRORODDMODEL;
    goto ErrorHandling;
  } 

   
  err = ERROROUTOFMETHODLIST; // in case none of the methods are available 
  for (i = Nothing - 1; i>=0 && pref[order[i]] > PREF_NONE; i--) {


    zaehler++;
    unimeth = (Methods) order[i];

    //    if (unimeth == Specific) {APMI(cov); XERR(err); APMI(cov);}
    //    assert(unimeth != Specific);

    if (PL >= PL_RECURSIVE) {
      LPRINT("%s : pref=%d\n", METHODNAMES[unimeth], pref[unimeth]);
    }
    
    meth = gaussmethod[unimeth];      
    //    printf("++**************** meth =%d %s\n", meth, METHODNAMES[unimeth]);
  
    addModel(&(cov->key), meth);
    cov->key->calling = cov;
    
    cov_model *key = cov->key;
    key->prevloc = loc;
    int role = cov->role != ROLE_GAUSS || key->nr != AVERAGE_INTERN 
      ? cov->role
      : ROLE_POISSON_GAUSS;


    //      PMI(cov);

    if ((err = CHECK(key, dim, xdim, GaussMethodType, cov->domown,
		     cov->isoown, cov->vdim2, role)) == NOERROR) {
      
  
      // PMI(key, "struct gauss");
    //     printf("checked %d\n", key->initialised);  
    //     assert(false);
      if ((err = STRUCT(key, NULL)) <= NOERROR) {

	//	PMI(cov);

        //PMI(key);	printf("initialised %d\n", key->initialised);
	  
	if (!key->initialised) {
	  
	  if ((err = CHECK(key, dim, xdim, GaussMethodType, cov->domown,
			   cov->isoown, cov->vdim2, role)) == NOERROR){
	    key->method = unimeth;
	    
	    NEW_STORAGE(stor, STORAGE, gen_storage);
	    
	    // 	    PMI(key);
	    err = INIT(key, role == ROLE_POISSON_GAUSS ? 2 : 0, cov->stor);

	    // printf("err ===== %d\n", err); if (err > 0) {XERR(err);}
	    
	    if (err == NOERROR) {
	      if (PL >= PL_REC_DETAILS) {
		LPRINT("returning to higher level...\n");
	      }
	      goto Initialized;
	    }
	  }
	}
      }

      //printf("back to gauss %d\n", err);
      
      if (PL >= PL_REC_DETAILS) {  
	LPRINT("back from %s: err =%d [%d, %s]\n",
	       METHODNAMES[unimeth], err, unimeth, METHODNAMES[unimeth]);
      }
      if (PL > PL_ERRORS) {
	char msg[2000]; errorMSG(err, msg); PRINTF("error in init: %s\n",msg);
      }
    } else {
      if (PL > PL_ERRORS) {
	char msg[2000]; errorMSG(err,msg); PRINTF("error in check: %s\n",msg);
      }
    }
    key = key->sub[0];
    //       printf("delete A back to gauss %d\n", err);
  
    COV_DELETE_WITHOUTSUB(&(cov->key));
    //   printf("donr A back to gauss %d\n", err);
  
    cov->key=key;
    key->calling = cov;      
  }
  //   printf("A back to gauss %d\n", err); assert(false);
  
  if (PL >= PL_REC_DETAILS) {
    // printf("B back to gauss %d\n", err);
    strcpy(PREF_FAILURE, "");
    int p, lp;
#define NMAX 14
    char lpd[255], pd[255], names[NMAX];
    names[NMAX-1] = '\0';
    for (i=0; i<Nothing; i++) {
      lp = locpref[i] - LOC_PREF_NONE;
      p = pref[i] - LOC_PREF_NONE;
      if (lp>0) {
	strcpy(lpd, "");
      } else {
	sprintf(lpd, "%s (locp) ", FailureMsg[-lp]);
      }
      if (p>0 || p == lp) {
	if (lp>0) strcpy(pd, "(method specific failure)");
	else strcpy(pd, ""); // strcpy(pd, "(method specific failure)");
      } else {
	sprintf(pd, "%s (pref)", FailureMsg[-p]);
      }
      strcopyN(names, METHODNAMES[i], NMAX-1);
      sprintf(dummy, "%s %-13s: %s%s\n", PREF_FAILURE, names, lpd, pd);
      strcpy(PREF_FAILURE, dummy);
    }
  }
  //   printf("C back to gauss %d\n", err);
  //PMI(cov);
   
  if (err != NOERROR) goto ErrorHandling;
  
 Initialized:
  // nachfolgende Zeile muss stehen, da Init in struct aufgerufen!
  err = gauss_init_settings(cov);
  cov->initialised = err == NOERROR;
  
 ErrorHandling:

  //printf("cov %ld %ld %ld\n", (long int) cov, (long int) alt, NULL);  assert(false);
  //   printf("E back to gauss %d\n", err);
 
  if (err <= NOERROR || zaehler==0){
    if (err <= NOERROR) {
      if (cov->key == NULL) next->simu.active = true;
      else cov->key->simu.active = true;
    }
    return err;
  }

  if (zaehler==1) {
    //   printf("F back to gauss %d\n", err);
    sprintf(ERROR_LOC, "Only 1 method found for '%s', namely '%s', that comes into question.\nHowever:", NICK(next), METHODNAMES[unimeth]);
    return err;
  }
  
  { // only for error reporting !
    cov_model *sub = next; 
    while(isAnyDollar(sub)) sub = sub->sub[0];
    sprintf(ERROR_LOC, "searching a simulation method for '%s': ",
	    NICK(sub));
  }

  // printf("provisorisch\n");  XERR(err);
  //    printf("G back to gauss %d\n", err);
 
  return ERROROUTOFMETHODLIST;
}


int init_gaussprocess(cov_model *cov, gen_storage *s) {
  /// ACHTUNG struct_gaussprocess hat bereits init aufgerufen!

  // direkter Aufruf, u.a. durch RPtbm, average, cutoff, intrinsic, hyperplane,
  //                    nugget, spectral !!!

  cov_model 
    //    *next = cov->sub[0],
    *key = cov->key;
  int err;

  ROLE_ASSERT_GAUSS;
  assert(key != NULL);


  if ((err = INIT(key, 0, s)) != NOERROR) return err;

  //  }
 	    
 if ((err = gauss_init_settings(cov)) != NOERROR) {
   return err;
  }
  //  assert(false);
  
  key->simu.active = true;

  return NOERROR;
}




void do_gaussprocess(cov_model *cov, gen_storage *s) {
  //  if (s == NULL) crash(cov);
  assert(s != NULL);
  // reopened by internal_dogauss
  location_type *loc = Loc(cov);		
  int i;
  long 
    total_pts = loc->totalpoints, 
    vdimtot = total_pts * cov->vdim2[0];
  char errorloc_save[nErrorLoc];
  double *res = cov->rf;
  cov_model 
    // *next = cov->sub[0],
    *key = cov->key ;
  bool loggauss = GLOBAL.gauss.loggauss;
 
  strcpy( errorloc_save,  ERROR_LOC);

  if (cov->simu.pair) {
    for (i=0; i<vdimtot; i++) res[i] = -res[i];
    cov->simu.pair = false;
    return;  
  } else {
    cov->simu.pair = GLOBAL.gauss.paired;
  }

  assert(key != NULL);
  //  assert(cov->stor != NULL);

  // falls gaussprocess als solcher aufgerufen, dann cov->stor != NULL
  // falls ueber *_USER, dann i.a.(?!) nicht gesetzt, sondern s:
  DO(key, cov->stor == NULL ? s : cov->stor);

  // (x^\lambda_1-1)/\lambda_1+\lambda_2

  if (loggauss) for (i=0; i<vdimtot; i++) res[i] = exp(res[i]);

  //   total = GetTotalPts(cov);
  
  strcpy( ERROR_LOC, errorloc_save);

}




//////////////////////////////////////////////////////////////////////
// Binary Gauss
#define BINARY_THRESHOLD GAUSSPROC_LAST + 1

void kappa_binaryprocess(int i, cov_model VARIABLE_IS_NOT_USED *cov, 
			 int *nr, int *nc){
  *nr = i == BINARY_THRESHOLD ? SIZE_NOT_DETERMINED : 
    i < BINARY_THRESHOLD ? 1 : -1;
  *nc = 1;
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

  // PMI(cov); printf("%d %d\n", key==NULL, isNegDef(next));
    
  if (key == NULL && isNegDef(next)) {
    if ((err = checkgaussprocess(cov)) != NOERROR) {
      return err;
    }
    
    COV(ZERO, sub, &v);
    if (v != 1.0) 
      SERR("binaryian requires a correlation function as submodel.");

  } else if (isProcess(sub)) {
    // if (sub->nr == GAUSSPROC) SERR("RTgauss is the default setting of RTbernoulli -- do not call it explicitely");
    role = cov->role == ROLE_BASE ? cov->role : role_of_process(sub->nr);
    
    if ((err = CHECK(sub, cov->tsdim, cov->xdimprev, ProcessType, 
		     cov->domown, cov->isoown,
		     SUBMODEL_DEP, 
		     role)) != NOERROR) return err;
    setbackward(cov, sub);
  } else SERR1("process type model required, but '%s' obtained", NICK(sub));

  cov->vdim2[0] = sub->vdim2[0];
  cov->vdim2[1] = sub->vdim2[1];
  return NOERROR;
}

int struct_binaryprocess(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  cov_model
    *next = cov->sub[0];
  int err;
  
  //PMI(cov);

  ROLE_ASSERT(ROLE_BERNOULLI);
  if (isNegDef(next)) {
    assert(cov->key == NULL);
    err = covcpy(&(cov->key), cov);

    // darauffolgende Zeile absichern:
    if (CovList[cov->nr].kappas != 2 || CovList[GAUSSPROC].kappas != 1) BUG;
    if (cov->key != NULL && PARAM(cov->key, BINARY_THRESHOLD) != NULL) {
      free(PARAM(cov->key, BINARY_THRESHOLD));
      cov->key->px[BINARY_THRESHOLD] = NULL;
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
     vdim = next->vdim2[0],
    vdimSq = vdim * vdim,
    vdimP1 = vdim + 1;
  assert(next->vdim2[0] == next->vdim2[1]);
 
  if ((variance = (double*) MALLOC(sizeof(double) * vdimSq)) == NULL ||
      (mean = (double*) CALLOC(vdim, sizeof(double))) == NULL
      ) goto ErrorHandling;

  if ((err = INIT(sub, 0, s)) != NOERROR) goto ErrorHandling;
  cov->rf = sub->rf;
  cov->origrf = false;
 
  if (isNegDef(next) || next->nr == GAUSSPROC) {
    GetInternalMean(next, vdim, mean);
    if (ISNAN(mean[0])) // GetInternalMean currently only allows ...
      GERR("'binaryprocess' currently only allows scalar fields - NA returned");
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
	    GERR("Vanishing sill not allowed in 'gaussprocess'");
	  sigma = sqrt(variance[w]);
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
  if (variance != NULL) free(variance);
  if (mean != NULL) free(mean);
 
  return err;
}



void do_binaryprocess(cov_model *cov, gen_storage *s){
  // reopened by internal_dogauss
  long j,
    tot = cov->prevloc->totalpoints,
    endfor = tot;
  int  i,
    pi, 
    npi = cov->nrow[BINARY_THRESHOLD],
    vdim = cov->vdim2[0];
  double 
    threshold,
    *p = P(BINARY_THRESHOLD),
    *rf = cov->rf;
  cov_model 
    *next=cov->sub[0];
  assert(cov->vdim2[0] == cov->vdim2[1]);
  if (isNegDef(next)) {
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
#define CHISQ_DEGREE 0


int checkchisqprocess(cov_model *cov) {
  cov_model
    *key = cov->key,
    *next = cov->sub[0];
  double *v = NULL;
  int err = NOERROR,
    xdim = cov->xdimprev, // could differ from tsdim in case of distances!
    dim = cov->tsdim;
  assert((Loc(cov)->distances && xdim==1) || xdim == dim);

  // PMI(cov);
  //  printf("name %d\n", cov->calling->sub[0]->nr);

  if (PisNULL(CHISQ_DEGREE)) SERR("degree of freedom must be given");
  if (key == NULL) {
    // if (next->nr >= FIRSTGAUSSPROC && next->nr <= LASTGAUSSPROC) {
    //   // wird verboten
    //   SERR("'%s' may not call a Gauss method", NICK(cov));
    // }
    if (!isGaussProcess(next) && !isNegDef(next))
      SERR1("Gaussian process required, but '%s' obtained", NICK(cov));
    if ((err = CHECK(next, dim, xdim, ProcessType, XONLY, cov->isoown, 
		     SUBMODEL_DEP, cov->role)) != NOERROR) {
      if ((err = CHECK(next, dim, xdim, NegDefType, KERNEL, SYMMETRIC,
		     SUBMODEL_DEP, ROLE_COV)) != NOERROR) {
      return err;
      }
    }
    
    int 
      vdim = next->vdim2[0],
      vdimSq = vdim * vdim;
    assert(vdim > 0);

    if ((v = (double*) MALLOC(sizeof(double) * vdimSq)) == NULL)
      return ERRORMEMORYALLOCATION;

    cov_model *sub = next;
    while (sub != NULL && isProcess(sub)) sub=sub->sub[0];
    if (sub == NULL) BUG;
    //printf("vdim = %d\n", vdim);
    if (sub->domprev == XONLY) {
      COV(ZERO, next, v);
    } else {
      assert(sub->domprev == KERNEL);
      NONSTATCOV(ZERO, ZERO, next, v);
    }
    
    int w,
      vdimP1 = vdim + 1;
    for (w=0; w<vdimSq; w+=vdimP1) {
      //printf("%d %f\n", w, v[w]);
      if (v[w] != 1.0) {
	free(v);
	SERR("chisq requires a correlation function as submodel.");
      }
    }
    free(v);  v=NULL;
    
    cov->vdim2[0] = sub->vdim2[0];
    cov->vdim2[1] = sub->vdim2[1];
     //no need to setbackward
  } else {
    if ((err = CHECK(key, cov->tsdim,  cov->xdimprev, ProcessType,
		       cov->domown, cov->isoown,
		       SUBMODEL_DEP, cov->role)) != NOERROR) return err;
    setbackward(cov, key);
  }
  //  printf("end name %d\n", cov->calling->sub[0]->nr);

  return NOERROR;
}
 
int struct_chisqprocess(cov_model *cov, 
			cov_model VARIABLE_IS_NOT_USED **newmodel) {
  cov_model
    *next = cov->sub[0];
  int err;

  ROLE_ASSERT_GAUSS;

  if (isNegDef(next)) {
    assert(cov->key == NULL);
    if ((err = covcpy(&(cov->key), next)) > NOERROR) return err;
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
    vdim = cov->vdim2[0];

  cov->simu.active = false;
  if ((err = INIT(sub, 2, s)) != NOERROR) return err;		  
 
  nmP1 = cov->mpp.moments + 1;
  for (i=0; i<vdim; i++) { 
    int 
      idx = i * nmP1,
      idxsub = i * nmP1sub;
    mean = sub->mpp.mM[idxsub + 1];
    m2 = sub->mpp.mM[idxsub + 2];
    variance =  m2 - mean * mean; 
    if (variance==0.0) SERR("Vanishing sill not allowed in 'gaussprocess'");
    if (ISNAN(mean)) // GetInternalMean currently only allows ...
      SERR("'chisqprocess' currently only allows scalar fields -- NA returned");
    cov->mpp.maxheights[i] = GLOBAL.extreme.standardmax * 
      GLOBAL.extreme.standardmax * m2;

    if (cov->mpp.moments >= 0) { 
      assert(cov->mpp.mM != NULL);
      cov->mpp.mM[idx + 0] = cov->mpp.mMplus[idx + 0] = 1.0; 
      if (cov->mpp.moments >= 1) {
	cov->mpp.mMplus[idx + 1] = m2;
	cov->mpp.mM[idx + 1] = RF_NA; // to do: richtiger Wert
	if (cov->mpp.moments >= 2)
	  cov->mpp.mM[idx + 2] = 3 * variance * RF_NA;//check the 4th moment of a Gaussian, to do
      }
    }
  }

  FieldReturn(cov);
  cov->simu.active = true;
 
  return NOERROR;
}



void do_chisqprocess(cov_model *cov, gen_storage *s){
  // reopened by internal_dogauss
  int  i, f,
    degree = P0INT(CHISQ_DEGREE),
    n = cov->prevloc->totalpoints * cov->vdim2[0];
  cov_model 
    *next=cov->sub[0],
    *key = cov->key,
    *sub = key == NULL ? next : key;
  double 
    *subrf = sub->rf,
    *rf = cov->rf;

  for (i=0; i<n; i++) rf[i] = 0.0;

  //  APMI(cov);

  for (f=0; f<degree; f++) {
    //  printf("f=%d %d\n", f, degree);
    DO(sub, s);   
    for (i=0; i<n; i++) {
      double dummy = subrf[i];
      rf[i] += dummy * dummy;
    }
  }
}

void rangechisqprocess(cov_model *cov, range_type *range) {
  rangegaussprocess(cov, range);
  range->min[CHISQ_DEGREE] = 1;
  range->max[CHISQ_DEGREE] = RF_INF;
  range->pmin[CHISQ_DEGREE] = 1;
  range->pmax[CHISQ_DEGREE] = 1000;
  range->openmin[CHISQ_DEGREE] = false;
  range->openmax[CHISQ_DEGREE] = true;
}


