/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Collection of auxiliary functions

 Copyright (C) 2001 -- 2006  Martin Schlather, 

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
#include "RFsimu.h"
#include "MPPFcts.h"
#include <assert.h>


double EXTREMES_STANDARDMAX = 3.0; // should be about 3 or 4

void SetExtremes(int *action, double *standardGausMax)
{
  if (*action) {
    if (*standardGausMax<=0.0) {
       if (GENERAL_PRINTLEVEL>0) PRINTF("\nERROR! `standardGausMax' not a positive number. Ignored.");
    } else  EXTREMES_STANDARDMAX = *standardGausMax;
  } else {
    *standardGausMax= EXTREMES_STANDARDMAX;
  }
}


void extremes_destruct(void **SExtremes){
  if (*SExtremes!=NULL) {
    extremes_storage *x;
    x = *((extremes_storage**) SExtremes);
    if (x->rf != NULL) free(x->rf);
    DeleteKeyNotTrend(&(x->key));
    free(*SExtremes);
    *SExtremes = NULL;
  }
}

void InitMaxStableRF(double *x, double *T, int *dim, int *lx, int *grid, 
		     int *Time,
		     int *covnr, double *ParamList, int *nParam, 
		     int *ncov, int *anisotropy, int *op,
		     int *method, 
		     int *distr,
		     int *keyNr,
		     int *error)
// for all parameters  see InitSimulateRF in RFsimu.cc
{ 
  key_type *key;

  strcpy(ERROR_LOC, "MaxStable");
  key = NULL;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORREGISTER; 
    goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  // I cannot see why a time component might be useful
  // hence cathegorically do not allow for the use!
  // if you allow for it check if it use is consistent with the
  // rest
  if (*Time) {*error=ERRORNOTPROGRAMMED; goto ErrorHandling;} 
  if (key->TrendModus != TREND_MEAN) {
    *error=ERRORTREND; goto ErrorHandling;} 
  if (method[0] == (int) MaxMpp) {
    int init_method[MAXCOV], v, n, m;
    mpp_storage* s; 
    // MaxMpp may not be mixed with other methods ! 
    /* ***********************************************************
       MPP FUNCTIONS
     * *********************************************************** */    
    for (v=0; v<*ncov; v++) {
      init_method[v] = (int) AdditiveMpp;
      if (v<*ncov-1 && op[v]) {
	ERRORMODELNUMBER = v;	
	*error=ERRORNOMULTIPLICATION;
	goto InternalErrorHandling;
      }
    }
    *error = internal_InitSimulateRF(x, T, *dim, *lx, (bool) *grid,
				     (bool) *Time, covnr, ParamList, *nParam, 
				     *ncov, (bool) *anisotropy, op,
				     init_method, *distr, key,
				     GENERAL_NATURALSCALING, CovFct);
    strcpy(ERROR_LOC, "");
    if (*error!=NOERROR) goto ErrorHandling;
    for (v=0; v<*ncov; v++) key->cov[v].method = MaxMpp;
    for (m=0; m<key->n_unimeth; m++) key->meth[m].unimeth = MaxMpp;

    // the basic functions in MPPFcts.cc are kept quite simple (e.g. height 1.0)
    // hence height must be corrected to "unit volume"; effectiveara is 
    // involved since below only standard exponential variables are considered;
    // not that effectivearea is not always the simulation window !
    // actcov corrects for the max-superposition of the actcov maxstable rf
    n = key->n_unimeth - 1;
    s = (mpp_storage*) key->meth[n].S;  
    s->factor = s->effectivearea / (s->integralpos * (double) (key->n_unimeth));
    for (v=0; v<n; v++) {
      mpp_storage* sp1; 
      s = (mpp_storage*) key->meth[v].S;  
      sp1 = (mpp_storage*) key->meth[v+1].S;  
      s->factor = s->effectivearea / s->integralpos / 
	  (sp1->effectivearea / sp1->integralpos);
    }
  } else { 
    /* ***********************************************************
       EXTREMAL GAUSSIAN
     * *********************************************************** */
    double sigma, meanDsigma, mean;
    int d;
    ce_param ce_save, local_save;
    methodvalue_type *meth; 
    extremes_storage* s;
  
    assert(*distr == DISTR_MAXSTABLE);
    DeleteKeyNotTrend(key); // should be improved ...
    key->active=false;
    meth = &(key->meth[0]);      
    SET_DESTRUCT(extremes_destruct, 0);
    if ((meth->S = (extremes_storage*) malloc(sizeof(extremes_storage)))==NULL){
	*error=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    s = (extremes_storage*) meth->S;
    KEY_NULL(&(s->key));
    s->rf = NULL;
 
    memcpy(&ce_save, &CIRCEMBED, sizeof(ce_param));
    memcpy(&local_save, &LOCAL_CE, sizeof(ce_param));
    LOCAL_CE.severalrealisations = CIRCEMBED.severalrealisations = true;
    *error = internal_InitSimulateRF(x, T, *dim, *lx, (bool) *grid, 
				     (bool) *Time, covnr, ParamList, *nParam, 
				     *ncov, (bool) *anisotropy, op, method, 
				     DISTR_GAUSS, &(s->key), 
				     GENERAL_NATURALSCALING, CovFct);
    memcpy(&CIRCEMBED, &ce_save, sizeof(ce_param));
    memcpy(&LOCAL_CE, &local_save, sizeof(ce_param));
    if (*error!=NOERROR) goto ErrorHandling;

    key->distribution = *distr;
    sigma = sqrt(CovFct(ZERO, *dim, s->key.cov, COVLISTALL, *ncov, *anisotropy));
    if (sigma==0) {*error = ERRORSILLNULL; goto InternalErrorHandling;}
    mean = s->key.mean;
    meanDsigma = mean / sigma;
    s->inv_mean_pos = 
      1.0 / (sigma * INVSQRTTWOPI * exp(-0.5 * meanDsigma * meanDsigma) 
	     + mean * /* correct next lines the 2nd factor is given */
	     // pnorm(meanDsigma, key->mean, sigma,1,0));
	     pnorm(0.0, mean, sigma, 1, 0) 
	     );
    s->assumedmax = EXTREMES_STANDARDMAX * sigma + ((mean>0) ? mean : 0);  
    *error=NOERROR;
    key->totalpoints = s->key.totalpoints;// DoSimulateRF regards whether
    //                                       key->totalpoints is initialised
    if ((s->rf =
	 (double*) malloc(sizeof(double) * key->totalpoints))==NULL) {
	*error=ERRORMEMORYALLOCATION; goto InternalErrorHandling;
    }   
    key->totalpoints = s->key.totalpoints;
    key->spatialdim = s->key.spatialdim;
    key->timespacedim = s->key.timespacedim;
    for (d=0; d<key->timespacedim; d++) key->length[d] = s->key.length[d];
    key->grid = s->key.grid;
    key->meth[0].unimeth = key->cov[0].method = ExtremalGauss;
    key->n_unimeth = 1;
    key->meth[0].actcov = 0; 
    key->active = true;
  }
  strcpy(ERROR_LOC, "");
  return;

 InternalErrorHandling:
  ErrorMessage(Nothing,*error); 
 ErrorHandling:
  key->active=false;
  return;
}



void DoMaxStableRF(int *keyNr, int *n, int *pairs, double *res, int *error)
{
  key_type *key;
  long totalpoints, control, counter;
  double poisson, invpoisson, *RES, threshold;


  *error=0; 
  key = NULL;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  if (*pairs) {*error=ERRORPAIRS; goto ErrorHandling;}

  key = &(KEY[*keyNr]);
  if (!key->active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}
  totalpoints = key->totalpoints;     

  counter = 0;

  if (key->cov[0].method == (int) MaxMpp) {
    /* ***********************************************************
       MPP FUNCTIONS
       * *********************************************************** */
      
    // what to do with the mean?? -> ignore it
    // what to do with the variance?? -> ignore it
    int  dim, d, actcov, m, ni;
    double min[MAXDIM], max[MAXDIM], next=false;
    long segment[MAXDIM+1], i; // just for printing messages;

    mpp_storage *s;
    mppmodel model;
    MPPRandom MppFct;
    methodvalue_type *meth; 
    covinfo_type *kc;
    
    actcov = key->n_unimeth;
    assert(actcov == key->ncov);

    for (segment[0] = 1, d = 0; d < key->timespacedim; d++)// used if simugrid
      segment[d + 1] = segment[d] * key->length[d];

    GetRNGstate();
    for (RES = res, ni=0; ni<*n; ni++, RES += key->totalpoints) {
      for  (i=0; i<totalpoints; i++) RES[i]=0.0;      
      for (m=0; m<actcov; m++) {
	meth = &(key->meth[m]); 
	assert(meth->unimeth == MaxMpp);
	kc = &(key->cov[meth->covlist[m]]);
	s = (mpp_storage*) meth->S;  
	MppFct = s->MppFct;
	dim = s->dim;
	
	control = 0;
	poisson = rexp(1.0); 
	invpoisson = 1.0 / poisson;
	while (control < totalpoints) {   
	  // the following is essentially the same as the code in MPP.cc!!
	  // how to reorganise it to get rid of the duplication of the code??
	  //	  
	  MppFct(s, min, max, &model); 
	  if (kc->simugrid) {
	    int start[MAXDIM], end[MAXDIM], resindex, index[MAXDIM],
	      segmentdelta[MAXDIM], dimM1;
	    double coord[MAXDIM], startcoord[MAXDIM];
	    
	    resindex = 0;
	    for (d=0; d<key->timespacedim; d++) {	 
	      // determine rectangle of indices, where the mpp function
	      // is different from zero
	      if (kc->genuine_dim[d]) {
	        if ((next = ((min[kc->idx[d]] > kc->xsimugr[XENDD[d]]) ||
			     (max[kc->idx[d]] < kc->xsimugr[XSTARTD[d]])))) { 
		  break;
		}
		if (min[kc->idx[d]] < kc->xsimugr[XSTARTD[d]]) {start[d]=0;}
		else start[d] = (int)((min[kc->idx[d]] - kc->xsimugr[XSTARTD[d]])
				      / kc->xsimugr[XSTEPD[d]]);
		// "end[d] = 1 + *"  since inequalities are "* < end[d]" 
		// not "* <= end[d]" !
		end[d] = (int) ((max[kc->idx[d]] - kc->xsimugr[XSTARTD[d]]) / 
				kc->xsimugr[XSTEPD[d]]);
		if (end[d] > key->length[d]) end[d] = key->length[d];
	      } else {
		start[d] = 0;
		end[d] = key->length[d];
	      }
	      // prepare coord starting value and the segment vectors for res
	      index[d]=start[d];
	      segmentdelta[d] = segment[d] * (end[d] - start[d]);
	      resindex += segment[d] * start[d];
	      coord[kc->idx[d]] = startcoord[d] = kc->xsimugr[XSTARTD[d]] + 
		(double)start[d] * kc->xsimugr[XSTEPD[d]];
	    }	
	    if (next) continue;
	      
	    // "add" mpp function to res
	    dimM1 = key->timespacedim - 1;
	    d = 0; // only needed for assert(d<dim)
	    while (index[dimM1]<end[dimM1]) {
	      double dummy;
	      assert(d<dim); 
	      dummy =  model(coord) * invpoisson;
	      if (RES[resindex] < dummy) RES[resindex]=dummy;
	      d=0;
	      index[d]++;
	      coord[kc->idx[d]] += kc->xsimugr[XSTEPD[d]];
	      resindex++;
	      while (index[d] >= end[d] && d < dimM1) { 
		// loop never entered if dim=1
		index[d]=start[d];
		coord[kc->idx[d]] = startcoord[d];
		resindex -= segmentdelta[d];
		d++; // if (d>=dim) break;
		index[d]++;
		coord[kc->idx[d]] += kc->xsimugr[XSTEPD[d]];
		resindex += segment[d];
	      }
	    }
	  } else {  // not a grid
	    double y[MAXDIM];
	    long j;
	    // the following algorithm can greatly be improved !
	    // but for ease, just the stupid algorithm
	    for (d=dim; d<dim; d++) y[d]=0.0;
	    for (j=i=0; i<totalpoints; i++) {
	      for (d=0; d<dim; d++) y[d] = kc->x[j++];
	      double dummy;
	      dummy =  model(y) * invpoisson;
	      if (RES[i] < dummy) RES[i]=dummy;
	    }
	  }
	  
	  poisson += rexp(1.0);  
	  invpoisson = 1.0 / poisson;
	  threshold = s->maxheight * invpoisson;
	  while ((control<totalpoints) && (RES[control]>=threshold)) control++;
	  if (GENERAL_PRINTLEVEL>=3) {
	    counter++;
	    if ((counter % 100==0) && (control<totalpoints))
	      PRINTF("%d, %d-th position: value=%f threshold=%f \n",
		     counter, control, RES[control], threshold); 
	  }   
	} // while control < totalpoints

	for (i=0; i<totalpoints; i++) RES[i] *= s->factor;
      } //v
    } // for, n
    PutRNGstate();
  } else {
      /* ***********************************************************
	 EXTREMAL GAUSSIAN
	 * *********************************************************** */
    int ni, i;
    extremes_storage* s;

    if (key->distribution!=DISTR_MAXSTABLE) {
      *error=ERRORWRONGINIT; 
      goto ErrorHandling;
    }
    s = (extremes_storage*) key->meth[0].S;
    assert(s!=NULL);
     
    for (RES = res, ni=0; ni<*n; ni++, RES += key->totalpoints) {
      for  (i=0; i<totalpoints; i++) RES[i]=0.0;     
    
      control = 0;
      
      if ((*error = InternalSimulate(&(s->key), s->rf)) != NOERROR)
	  return;
      //  to get some consistency with GaussRF concerning Randomseed,
      //  DoSimulate must be called first; afterwards, it is called after 
      //  creation of Poisson variable 
      
      GetRNGstate();
      poisson = rexp(1.0); 
      PutRNGstate();
      invpoisson = 1.0 / poisson;
           
      while (true) {
	for (i=0; i<totalpoints; i++) {// 0, not control, since we approximate
	  //                     and 0 as starting control increases precision
	  s->rf[i] *= invpoisson;
	  if (RES[i] < s->rf[i]) RES[i]=s->rf[i];
	}
	GetRNGstate();
	poisson += rexp(1.0); 
	PutRNGstate();
	invpoisson = 1.0 / poisson;
	threshold = s->assumedmax * invpoisson;
	while ((control<totalpoints) && (RES[control]>=threshold)) {control++;}
	if (control>=totalpoints) break;
	if ((*error = InternalSimulate(&(s->key), s->rf)) != NOERROR)
	   return;
	if (GENERAL_PRINTLEVEL>=3) {
	  counter++;
	  if (counter % 10==0) 
	      PRINTF("%d, %d-th position: value=%f threshold=%f \n",
		     counter, control, RES[control], threshold); 
	}
      } // while
      
      for (i=0; i<totalpoints; i++) RES[i] *= s->inv_mean_pos;
    } // for, n
  }
  return;
  
 ErrorHandling:
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);  
  if (key!=NULL) key->active=false;  
}


