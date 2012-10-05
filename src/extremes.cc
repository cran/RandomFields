/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of max-stable random fields

 Copyright (C) 2001 -- 2011 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include <math.h>
#include <stdio.h>
#include "RF.h"
#include "Covariance.h"
 
#include "MPPstandard.h"



void SetExtremes(int *action, double *standardGausMax)
{
  extremes_param *gp = &(GLOBAL.extremes);
  if (*action) {
    if (*standardGausMax<=0.0) {
       if (PL>0)
	 PRINTF("\nERROR! `standardGausMax' not a positive number. Ignored.");
    } else  gp->standardmax = *standardGausMax;
  } else {
    *standardGausMax= gp->standardmax;
  }
}


void extremes_destruct(void **SExtremes){
  if (*SExtremes!=NULL) {
    extremes_storage *x;
    x = *((extremes_storage**) SExtremes);
    if (x->rf != NULL) free(x->rf);
    KEY_DELETE(&(x->key));
    free(*SExtremes);
    *SExtremes = NULL;
    // do not free mpp since only pointed
  }
}

void doMaxMpp(method_type *Emeth, res_type * res) {
  
//   vassert(false);
    

  /* ***********************************************************
     MPP FUNCTIONS
     * *********************************************************** */
  
  // what to do with the mean?? -> ignore it
  // what to do with the variance?? -> ignore it
  GetRNGstate();
  extremes_storage* Es = (extremes_storage*) (Emeth->S);
  key_type *key = &(Es->key);
  method_type *meth = key->meth;

  // assert(false);

  // jump leading dollars
  // note that the variance parameter is and be ignored here, in contrast
  // the additive models !
  while (meth->nr != -RandomCoin-1 && meth->nsub==1 && meth->sub[0]!=NULL) {
      meth = meth->sub[0];
  }
  assert(meth->nr == -RandomCoin-1);

  STANDARDDEFINITIONEN(meth->cov);

  // mpp_model randomcoin = C->randomcoin;
//  printf("%ld", s);
  // PrintModelInfo(cov); // assert(false);

  double  poisson, invpoisson;
  res_type threshold;
  long control = 0;
//    ntot= s->ntot,
  int  dim = cov->tsdim,
    spatial = loc->totalpoints;
  mpp_get mppget = CovList[meth->cov->nr].mppget; 

//  print("name =%s %d %d %ld %ld %ld\n", 
//	CovList[meth->cov->nr].name, meth->cov->nr,  S2ISO, 
//	CovList[S2ISO].mppget, CovList[meth->cov->nr].mppget,
//	mppget);
//  assert(S2ISO == meth->cov->nr);
//  assert(mppget != NULL);

  STANDARDINIT;

// assert(false);
  poiss_initu(s);
  s->location = unif_u;
  poisson = rexp(1.0); 
  invpoisson = 1.0 / poisson;

  // assert(false);
 for(n=0; control < total; n++) { 
    // n nur zur Unterhaltung; loop allein ueber control !
//    printf("%d\n", n);

//     printf("%f\n", invpoisson);

    randomcoin(s, cov);
    Radius = s->effectiveRadius;
    
    if (loc->grid && meth->type <= TypeDiag) {
      STANDARDGRID;
      if (inside) { 
	while (true) {
	  res_type dummy;

//	  printf("mppget %ld %ld\n", mppget, mppget_Gauss);
	  dummy = (res_type) mppget(x, s->u, cov, s) * invpoisson;
//	  assert(false);
//	  printf("%f u=%f %f %f\n", x[0], s->u[0], dummy, invpoisson);

	  if (res[zaehler] < dummy) res[zaehler]= dummy;
	  STANDARDINKREMENT;
	}
      }
    } else {  // not a grid
      double *x=meth->space;
      // the following algorithm can greatly be improved !
	// but for ease, just the stupid algorithm
      for (zaehler=0; zaehler<spatial; zaehler++, x += dim) {
	  res_type dummy = (res_type) mppget(x, s->u, cov, s) * invpoisson;
	if (res[zaehler] < dummy) res[zaehler] = dummy;
      }
    }
    
    poisson += rexp(1.0);  
    invpoisson = 1.0 / poisson;
    threshold = (res_type) (s->maxheight * invpoisson);
    while ((control<total) && (res[control]>=threshold)) control++;

    if (n >= nthreshold) {
       nthreshold += deltathresh;
       PRINTF("%d %d-th position: value=%f threshold=%f \n",
	      n, control, (double) res[control], (double) threshold); 
    }
    R_CheckUserInterrupt();
       
  } // while control < total
 
 if ((factor = (res_type) s->factor) != 1.0)
    for (zaehler=0; zaehler<total; zaehler++) 
	res[zaehler] *= factor;
  PutRNGstate();
}

void doExtrGauss(method_type *Emeth, res_type * res) {
  /* ***********************************************************
     EXTREMAL GAUSSIAN
     * *********************************************************** */
  extremes_storage* Es = (extremes_storage*) (Emeth->S);
  key_type *key = &(Es->key);
  method_type *meth = key->meth;
  
  location_type *loc = meth->loc;
  globalparam *gp = meth->gp;
  int i, err,
    total = loc->totalpoints,
    every = gp->general.every,
    nthreshold = (every>0) ? every : MAXINT,
    deltathresh = nthreshold;     
  long n = 0, 
    control=0;
  double  poisson, invpoisson;
  res_type factor, threshold;

  if ((err = internal_DoSimulateRF(&(Es->key), 1, Es->rf)) != NOERROR)
    XERR(err);
  //  to get some consistency with GaussRF concerning Randomseed,
  //  DoSimulate must be called first; afterwards, it is called after 
  //  creation of Poisson variable 
    
  GetRNGstate();
  poisson = rexp(1.0); 
  PutRNGstate();
  invpoisson = 1.0 / poisson;
  
  while (true) {
    for (i=0; i<total; i++) {// 0, not control, since we approximate
      //                     and 0 as starting control increases precision
	res_type dummy = (res_type) (Es->rf[i] * invpoisson);
      if (res[i] < dummy) res[i]=dummy;
    }
    GetRNGstate();
    poisson += rexp(1.0); 
    PutRNGstate();
    invpoisson = 1.0 / poisson;
    threshold = (res_type) (Es->assumedmax * invpoisson);
    while ((control<total) && (res[control]>=threshold)) {control++;}
    if (control>=total) break;
    if ((err = internal_DoSimulateRF(&(Es->key), 1, Es->rf)) != NOERROR)
      XERR(err);

    if (n >= nthreshold) {
      PRINTF("%d %d-th position: value=%f threshold=%f \n",
	     n, control, res[control], (double) threshold); 
      nthreshold += deltathresh;
    }
    n++;

  } // while
  
  if ((factor = (res_type) Es->inv_mean_pos) != 1.0)
    for (i=0; i<total; i++) res[i] *= factor;
}


void InitMaxStableRF(double *x, double *T, int *dim, int *lx, int *grid, 
		     int *Time,
		     int *distr,
		     int *keyNr,
		     int expected_number_simu,
		     int *err)
// for all parameters  see InitSimulateRF in RFsimu.cc
{ 

  key_type *key = KEY + *keyNr;
  extremes_param *lp;
  cov_model *cov = STORED_MODEL[MODEL_SIMU];
  extremes_storage* Es;
  globalparam newparam;
  method_type *Emeth;
 
  strcpy(ERROR_LOC, "MaxStable");
  strcpy(PREF_FAILURE, "");
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *err=ERRORREGISTER; 
    goto ErrorHandling;
  }

//  PrintModelInfo(cov);
//      assert(false);    


  KEY_DELETE_NOTTREND(key);
  KEY_NULLNOTTREND(key);

  if (*Time) {*err=ERRORNOTPROGRAMMED; goto ErrorHandling;} 
  if (key->trend.TrendModus != TREND_MEAN) {
    *err=ERRORTREND; goto ErrorHandling;} 
  
  key->simu.distribution=*distr; 
  *err = loc_set(x, T, *dim, *lx, (bool) *Time, (bool) *grid, 
		 &(key->loc), &GLOBAL); 
  if (*err) goto ErrorHandling;
  memcpy(&(key->gp), &(GLOBAL), sizeof(globalparam));
  lp = &(key->gp.extremes);

  Emeth = key->meth = (method_type*) malloc(sizeof(method_type));
  METHOD_NULL(Emeth);
 
  Emeth->destruct = extremes_destruct; //  SET_DESTRUCT(extremes_destruct);

  if ((Emeth->S =(extremes_storage*) malloc(sizeof(extremes_storage)))==NULL){
    *err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  Es = (extremes_storage*) Emeth->S;    // key->meth->S
  KEY_NULL(&(Es->key));
  Es->rf = NULL;
  Emeth->nr = -1-(Es->meth = cov->user[MaxMpp] ? MaxMpp : ExtremalGauss);


//  printf("%d %d %d %d\n", MaxMpp, ExtremalGauss, Emeth->nr, cov->user[MaxMpp]); assert(false);

//  PrintModelInfo(cov);
//  printf("%d %d\n", cov->pref[MaxMpp], cov->pref[ExtremalGauss]);
//  assert(false);
    
  // I cannot see why a time component might be useful
  // hence cathegorically do not allow for the use!
  // if you allow for it check if it use is consistent with the
  // rest

  //PrintModelInfo(STORED_MODEL[MODEL_SIMU]);
  
  strcpy(ERROR_LOC, "");
  strcpy(PREF_FAILURE, "");
  memcpy(&newparam, &(key->gp), sizeof(globalparam));
  newparam.general.storing = true;    

 

  *err = internal_InitSimulateRF(x, T, *dim, *lx, (bool) *grid,
				 (bool) *Time, 
				 Es->meth == MaxMpp ? *distr : DISTR_GAUSS,
				 &(Es->key), &newparam, 1000, 
				 STORED_MODEL + MODEL_SIMU // hier cov
				 // wuerde zur Katastrophe fuehren
				 // da  STORED_MODEL + MODEL_SIMU  von
				 // internal_InitSimulateRF nicht auf 0 gesetzt
    );
  cov = Es->key.cov;
  covcpy(&(key->cov), cov, false, true);

  if (*err!=NOERROR) goto ErrorHandling;
//  genuinemeth = Es->key.meth;



 
//   printf("%ld %ld %d\n", genuinemeth->S, genuinemeth, genuinemeth->nr);
//   printf("XYOH!!\n");
// PrintMethodInfo(key->meth); assert(false);

//   printf("XOH!!\n");
//   PrintMethodInfo(genuinemeth); assert(false);

 
  if (Emeth->nr == -1-MaxMpp) {
    Emeth->domethod = doMaxMpp;
    mpp_storage *mpp;
    method_type *genuinemeth = Es->key.meth;

    
    //  PrintMethodInfo(genuinemeth);
    while (genuinemeth != NULL && genuinemeth->nr != -RandomCoin-1) {
	genuinemeth = genuinemeth->sub[0];
    }
    if (genuinemeth == NULL) {
	*err=ERRORNOGENUINEMETHOD;
	goto ErrorHandling;
    }
    

//  while (genuinemeth->nr != -RandomCoin-1 && genuinemeth->nsub==1 && genuinemeth->sub[0]!=NULL) {
//	if (genuinemeth->nr != DOLLAR)  {
//	    PrintMethodInfo(meth);
//	    ERR("unexpected method for shape function -- please contact author");
//	}
//	// printf("%d\n", genuinemeth->nr);
//	genuinemeth = genuinemeth->sub[0];
//    }

//   printf("%d %d %d %d\n", genuinemeth->nr, -RandomCoin, GATTER, LASTGATTER);
//    if (genuinemeth->nr != -RandomCoin-1) 
//	ERR("shape function could not be detected -- please contact author");
    mpp = (mpp_storage*) genuinemeth->S; 
    assert(mpp != NULL);

    
//     printf("XOH!! %d\n", mpp);
//   PrintMethodInfo(Es->key.meth);

    // MaxMpp may not be mixed with other methods ! 
    /* ***********************************************************
       MPP FUNCTIONS
     * *********************************************************** */   

//    for (m=0; m<=n; m++) key->meth[m].unimeth = MaxMpp;

    // the basic functions in MPPFcts.cc are kept quite simple (e.g. height 1.0)
    // hence height must be corrected to "unit volume"; effectiveara is 
    // involved since below only standard exponential variables are considered;
    // not that effectivearea is not always the simulation window !
    // xxx corrects for the max-superposition of the xxx maxstable rf

//    printf("%d\n", genuinemeth->S);
//    printf("%f\n", mpp->effectivearea);
//    printf("%f\n", mpp->integralpos);
 
    mpp->factor = mpp->effectivearea / (mpp->integralpos); 

//       printf("factor = %f %f %f\n",   
//	   mpp->factor, mpp->effectivearea, (mpp->integralpos)); 
       //  assert(false);

//    assert(false); // kkkk

    //    s = (mpp_storage*) key->meth[n].S;  
    //   s->factor = s->effectivearea / (s->integralpos * (double) (key->m));
    //    for (v=0; v<n; v++) {
    //     mpp_storage* sp1; 
    //     s = (mpp_storage*) key->meth[v].S;  
    //     sp1 = (mpp_storage*) key->meth[v+1].S;  
    //     s->factor = s->effectivearea / s->integralpos / 
    //	  (sp1->effectivearea / sp1->integralpos);
    //  }
    
  } else { 
    /* ***********************************************************
       EXTREMAL GAUSSIAN
     * *********************************************************** */

//    PrintModelInfo(cov);
//    printf("%d\n", cov->user[ExtremalGauss]);

//    assert(cov->user[ExtremalGauss]);
    assert(cov->user[MaxMpp] == PREF_NONE);
    assert(*distr == DISTR_MAXSTABLE);
 
    Emeth->domethod = doExtrGauss;

    double sigma, meanDsigma, mean;

    key->simu.distribution = *distr;
    CovList[cov->nr].cov(ZERO, cov, &sigma);
    if (sigma==0) {*err = ERRORSILLNULL; goto ErrorHandling;}
    sigma = sqrt(sigma);
    mean = Es->key.trend.mean;
    meanDsigma = mean / sigma;
    Es->inv_mean_pos = 
      1.0 / (sigma * INVSQRTTWOPI * exp(-0.5 * meanDsigma * meanDsigma) 
	     + mean * /* correct next lines the 2nd factor is given */
	     // pnorm(meanDsigma, key->mean, sigma,1,0));
	     pnorm(0.0, mean, sigma, 1, 0) 
	     );
  
    Es->assumedmax = lp->standardmax * sigma + ((mean>0) ? mean : 0);  
    *err=NOERROR;
    key->simu.active = true;

    if ((Es->rf =
      (res_type*) malloc(sizeof(res_type) * Es->key.loc.totalpoints))==NULL) {
	*err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
  }

  strcpy(ERROR_LOC, "");

  key->simu.active = true;  
  return;

 ErrorHandling:
  char ErrM[500], ErrM2[800];
  key->simu.active=false;
  // printf("!!!!!!! %s", PREF_FAILURE);

  errorMSG(*err, ErrM);
  sprintf(ErrM2, "%s%s", PREF_FAILURE, ErrM);
  ERR(ErrM2);
//  return;
}




void DoMaxStableRF(int *keyNr, int *n, int *pairs, res_type *res, int *err)
{
  key_type *key;
  long total;
  res_type *RES;
  int ni, i;
  method_type *meth;


  *err=0; 
  key = NULL;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *err=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  if (*pairs) {*err=ERRORPAIRS; goto ErrorHandling;}

  key = &(KEY[*keyNr]);
  if (!key->simu.active) {*err=ERRORNOTINITIALIZED; goto ErrorHandling;}
  total = ((extremes_storage*)key->meth->S)->key.loc.totalpoints;     

  meth = key->meth;
  for (RES = res, ni=0; ni<*n; ni++, RES += total) {
    for (i=0; i<total; i++) RES[i]=0.0;    
    meth->domethod(meth, RES);
  }


  return;
  
 ErrorHandling:
  if (PL>0) ErrorMessage(Nothing,*err);  
  if (key!=NULL) key->simu.active=false;  
}



