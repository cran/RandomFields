/* 
 Authors
 Martin Schlather, martin.schlather@cu.lu

 Collection of auxiliary functions

 Copyright (C) 2001 -- 2004  Martin Schlather, 

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


Real EXTREMES_STANDARDMAX = 3.0;


void SetExtremes(int *action, Real *standardGausMax)
{
  switch(*action) {
  case 0 :
    if (*standardGausMax<=0.0) {
       if (GENERAL_PRINTLEVEL>0) PRINTF("\nERROR! `standardGausMax' not a positive number. Ignored.");
    } else  EXTREMES_STANDARDMAX = *standardGausMax;
    break;
  case 1 :
    *standardGausMax= EXTREMES_STANDARDMAX;
    if (GetNotPrint) break;
  case 2 :
    PRINTF("\nMax-stable Random fields\n========================\nstandardGausMax=%e\n",
	   EXTREMES_STANDARDMAX);
    break;
  default : if (GENERAL_PRINTLEVEL>0) PRINTF("unknown action\n");
  }
}

typedef struct extremes_storage{
  Real *rf;
  Real inv_mean_pos, assumedmax;
} extremes_storage;

void extremes_destruct(void **SX){
  if (*SX!=NULL) {
    extremes_storage *x;
    x = *((extremes_storage**) SX);
    if (x->rf != NULL) free(x->rf);
    free(*SX);
    *SX = NULL;
  }
}

void InitMaxStableRF(Real *x, Real *T, int *dim, int *lx, int *grid, 
		     int *Time,
		     int *covnr, Real *ParamList, int *nParam, 
		     Real *mean,
		     int *ncov, int *anisotropy, int *op,
		     int *method, 
		     int *distr,
		     int *keyNr,
		     int *error)
/* for all parameters except precision see InitSimulateRF in the 
 * package rf!
 *
 * assumedmax: assumed maximum value of a Gaussian random field;
 *             should be about 3 or 4
 *     
 */
{ 
  Real sigma,meanDsigma;
  extremes_storage *es;
  int init_method[MAXCOV] ;
  int gauss_distr = DISTR_GAUSS;
  int i;
  key_type *key;
  bool storing;

  storing = GENERAL_STORING;
  GENERAL_STORING = true;    
  key = NULL;
  es = NULL;
  // I cannot see why a time component might be useful
  // hence cathegorically do not allow for the use!
  // if you allow for it check if it use is consistent with the
  // rest
  if (*Time) {*error=ERRORNOTPROGRAMMED; goto ErrorHandling;}  


  if (method[0] == (int) MaxMpp) {
    /* ***********************************************************
       MPP FUNCTIONS
     * *********************************************************** */    
    if (*ncov!=1) {*error=ERRORFAILED; goto InternalErrorHandling;}
    for (i=0; i<*covnr; i++) init_method[i] = (int) AdditiveMpp;

    InitSimulateRF(x, T, dim, lx, grid, Time, covnr, ParamList, nParam, 
		   mean, ncov, anisotropy, op,
		   init_method, distr, keyNr, error);
    key = &(KEY[*keyNr]);

    assert(*error>=0);
    if (*error!=NOERROR) {goto ErrorHandling;}
    key->method[0] = MaxMpp;
    assert( (key->SX==NULL) && (key->destructX==NULL));
  } else { 
    /* ***********************************************************
       EXTREMAL GAUSSIAN
     * *********************************************************** */
    InitSimulateRF(x, T, dim, lx, grid, Time, covnr, ParamList, nParam, 
		   mean, ncov, anisotropy, op,
		   method,  &gauss_distr, keyNr, error);
    key = &(KEY[*keyNr]);

    assert(*error>=0);
    if (*error!=NOERROR) {goto ErrorHandling;}
    key->distribution = DISTR_MAXSTABLE;
    
    if (key->SX==NULL) {
      key->destructX = extremes_destruct;
      if ((key->SX=
	   (extremes_storage*) malloc(sizeof(extremes_storage)))==NULL) {
	*error=ERRORMEMORYALLOCATION; goto InternalErrorHandling;
      }
      es = (extremes_storage*) key->SX;
      if ((es->rf =
	   (Real*) malloc(sizeof(Real) *  key->totalpoints))==NULL) {
	*error=ERRORMEMORYALLOCATION; goto InternalErrorHandling;
      }
    } else {es = (extremes_storage*) key->SX;}
    sigma = sqrt(CovFct(ZERO,*dim,covnr,op,key->param,*ncov,*anisotropy));
    if (sigma==0) {*error = ERRORSILLNULL; goto InternalErrorHandling;}
    meanDsigma = key->mean / sigma;
    es->inv_mean_pos = 
      1.0 / (sigma * INVSQRTTWOPI * exp(-0.5 * meanDsigma * meanDsigma) 
	     + key->mean * /* correct next lines the 2nd factor is given */
	     // pnorm(meanDsigma, key->mean, sigma,1,0));
	     pnorm(0.0, key->mean, sigma,1,0) 
	     );
  for (i=0;i<key->totalpoints;i++) es->rf[i]=0.0;
    
    es->assumedmax = EXTREMES_STANDARDMAX * sigma + 
      ((key->mean>0) ? key->mean : 0);  
    *error=NOERROR;
  } 
  GENERAL_STORING = key->storing = storing;  
  return;

 InternalErrorHandling:
  ErrorMessage(Nothing,*error); 
 ErrorHandling:
  GENERAL_STORING = storing;  
  key->active=false;
  return;
}



void DoMaxStableRF(int *keyNr, Real *res, int *error)
{
  key_type *key;
  extremes_storage *es;  
  long totalpoints,i,control;
  Real poisson, invpoisson, threshold, *zw, *RES;
  long counter=0; // just for printing messages
  bool next=false, storing;

  *error=0; 
  zw = NULL;
  key = NULL;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);

  if (!key->active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}
  storing = key->storing;
  key->storing = true;

  if (key->method[0] == (int) MaxMpp) {
     /* ***********************************************************
       MPP FUNCTIONS
     * *********************************************************** */
    
    // what to do with the mean?? -> ignore it
    // what to do with the variance?? -> ignore it
    
    int  dim, d, v;
    Real  min[MAXDIM], max[MAXDIM], factor;
    long segment[MAXDIM+1];
    mpp_storage *s;
    mppmodel model;
    MPPRandom MppFct;
    
    GetRNGstate();
   assert(key->ncov>0);
    if (key->ncov>1) {
      // not used
      if ((zw=(Real *)malloc(sizeof(Real)*key->totalpoints))==0){
	*error=ERRORMEMORYALLOCATION;goto ErrorHandling;}
      RES = zw;
      for  (i=0; i<key->totalpoints; i++) res[i]=zw[i]=0.0;
    } else {
      for  (i=0; i<key->totalpoints; i++) res[i]=0.0;
      RES = res;
    }

    s = (mpp_storage*) key->S[0];  
    dim = s->timespacedim;
    assert(dim>0);

    segment[0] = 1;
    for (d=0; d<dim; d++) 
      segment[d+1] = segment[d] * key->length[d];
    totalpoints = key->totalpoints;
     
     
    for (v=0; v<s->actcov; v++) {
      MppFct = s->MppFct[v];
      
      factor = s->effectivearea[v] / s->integralpos[v];
      
      control = 0;
      poisson = -log(1.0 - UNIFORM_RANDOM); /* it is important to write 1-U
					       as the random generator returns
					       0 inclusively; the latter would
					       lead to an error !
					      */
      invpoisson = 1.0 / poisson;
      while (control<totalpoints) {   
	// the following is essentially the same as the code in MPP.cc!!
	// how to reorganise it to get rid of the duplication of the code??
	//
	  
	MppFct(s, v, min, max, &model );

	if (s->grid) {
	  int start[MAXDIM], end[MAXDIM], resindex, index[MAXDIM],
	    segmentdelta[MAXDIM], dimM1;
	  Real coord[MAXDIM], startcoord[MAXDIM];
	  
	  // determine rectangle of indices, where the mpp function
	  // is different from zero
	  for (d=0; d<dim; d++) {	 
	    if (next = ((min[d] > key->x[d][XEND]) ||
			(max[d] < key->x[d][XSTART]))) { break;}
	    if (min[d]< key->x[d][XSTART]) {start[d]=0;}
	    else start[d] = (int) ((min[d] - key->x[d][XSTART]) / 
				   key->x[d][XSTEP]);
	    // "end[d] = 1 + *"  since inequalities are "* < end[d]" 
	    // not "* <= end[d]" !
	    end[d] = (int) ((max[d] - key->x[d][XSTART]) / 
			    key->x[d][XSTEP]);
	    if (end[d] > key->length[d]) end[d] = key->length[d];	  
	  }	
	  if (!next) {
	    // prepare coord starting value and the segment vectors for res
	    resindex = 0;
	    for (d=0; d<dim; d++) {
	      index[d]=start[d];
	      segmentdelta[d] = segment[d] * (end[d] - start[d]);
	      resindex += segment[d] * start[d];
	      coord[d] = startcoord[d] = 
		key->x[d][XSTART] + (Real)start[d] * key->x[d][XSTEP];
	    }
	    
	    // "add" mpp function to res
	    dimM1 = dim - 1;
	    d = 0; // only needed for assert(d<dim)
	    while (index[dimM1]<end[dimM1]) {
	      register Real dummy;
	      assert(d<dim); 
	      dummy =  model(coord) * invpoisson;
	      if (RES[resindex] < dummy) RES[resindex]=dummy;
	      d=0;
	      index[d]++;
	      coord[d]+=key->x[d][XSTEP];
	      resindex++;
	      while (index[d] >= end[d]) { 
		if (d>=dimM1) break;  // below not possible
		//                       as dim==1 will fail!
		index[d]=start[d];
		coord[d]=startcoord[d];
		resindex -= segmentdelta[d];
		d++; // if (d>=dim) break;
		index[d]++;
		coord[d]+=key->x[d][XSTEP];
		resindex += segment[d];
	      }
	    }
	  } /* next */ 
	} else {  // not agrid
	  Real y[MAXDIM];
	  long j;
	  // the following algorithm can greatly be improved !
	  // but for ease, just the stupid algorithm
	  for (j=i=0; i<totalpoints; i++) {
	    for (d=0; d<dim; d++,j++) { y[d] = s->x[j]; }
	    register Real dummy;
	    dummy =  model(y) * invpoisson;
	    if (RES[i] < dummy) RES[i]=dummy;
	  }
	}   
	
	poisson -= log(1.0 - UNIFORM_RANDOM);  
	invpoisson = 1.0 / poisson;
	threshold = s->maxheight[v] * invpoisson;
	while ((control<totalpoints) && (RES[control]>=threshold)) {control++;}
	if (GENERAL_PRINTLEVEL>=3) {
	  counter++;
	  if ((counter % 100==0) && (control<totalpoints))
	    PRINTF("%d: %d %f %f \n",counter,control,RES[control],
		   threshold); 
	}      
      } // while control < totalpoints
      if (key->ncov>1) {  // void -- maybe for future extensions
	assert(false);
	register Real dummy;
	for (i=0; i<totalpoints;i++) {
	  if (RES[i] < (dummy=RES[i]*factor)) RES[i]=dummy;
	  RES[i] = 0.0;
	}
      } else {
	for (i=0; i<totalpoints;i++) RES[i] *= factor;
      }
    } //v
  
    //if (key->param[NUGGET]>0) {
    //  Real nugget;
    //  nugget = key->param[NUGGET];
    //  for (i=0; i<totalpoints; i++) {
    //	register Real dummy;
    //	dummy = - nugget / log(1.0 - UNIFORM_RANDOM);
    //	if (res[i] < dummy) res[i]=dummy;
    //   }
    //}
    PutRNGstate();
  } else {
     /* ***********************************************************
       EXTREMAL GAUSSIAN
     * *********************************************************** */

    if (key->distribution!=DISTR_MAXSTABLE) {
      *error=ERRORWRONGINIT; 
      goto InternalErrorHandling;
    }
    es = (extremes_storage*) key->SX;
    assert(es!=NULL);
    totalpoints = key->totalpoints;
    
    control = 0;

    DoSimulateRF(keyNr, es->rf, error);
    //  to get some consistency with GaussRF concerning Randomseed,
    //  DoSimulate must be called first; afterwards, it is called after 
    //  creation of Poisson variable 
    if (*error)  {goto ErrorHandling;}
    if (!key->active){ *error=ERRORNOTINITIALIZED; goto InternalErrorHandling; }

    GetRNGstate();
    poisson = -log(1.0 - UNIFORM_RANDOM); 
    // it is important to write 1-U as the random generator returns
    //   0 inclusively; the latter would lead to an error ! 
    PutRNGstate();
    invpoisson = 1.0 / poisson;

    while (true) {
      for (i=0; i<totalpoints; i++) {// 0, not control, since we approximate
	//                     and 0 as starting control increases precision
	es->rf[i] *= invpoisson;
	if (res[i] < es->rf[i]) res[i]=es->rf[i];
      }
      GetRNGstate();
      poisson -= log(1.0 - UNIFORM_RANDOM);  // !! -= !! 
      PutRNGstate();
      invpoisson = 1.0 / poisson;
      threshold = es->assumedmax * invpoisson;
      while ((control<totalpoints) && (res[control]>=threshold)) {control++;}
      if (control>=totalpoints) break;
      DoSimulateRF(keyNr, es->rf, error);
      // first time no error; no reason to assume that an error occurs now;
      // so error is not checked 
     
      if (GENERAL_PRINTLEVEL>=3) {
	counter++;
	if (counter % 10==0) 
	  PRINTF("%d: %d %e %e \n",counter,control,invpoisson,threshold); 
      }
    } // while

    for (i=0; i<key->totalpoints; i++) 
      res[i] *= es->inv_mean_pos;
  }
  key->storing = storing;
  key->active=GENERAL_STORING && key->storing;
  if (zw!=NULL) free(zw);
  return;
  
 InternalErrorHandling:
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);  
 ErrorHandling:
  if (zw!=NULL) free(zw);
  if (key!=NULL) key->active=false;  
}


