#include <math.h>
#include <stdio.h>
#include "RFsimu.h"
#include "extremes.h"
#include "MPPFcts.h"
#include <Rmath.h>
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

void InitMaxStableRF(Real *x, Real *y, Real *z, int *dim, int *lx, 
		     int *grid, int *covnr, Real *ParamList, int *nParam, 
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
  int init_method = (int) AdditiveMpp;
  int gauss_distr = DISTR_GAUSS;
  int i;
 
  if (*method == (int) MaxMpp) {
    /* ***********************************************************
       MPP FUNCTIONS
     * *********************************************************** */
    
    InitSimulateRF(x,y,z, dim, lx, grid, covnr, ParamList, nParam, 
		   &init_method, distr, keyNr, error);
    if (*error>0) {goto ErrorHandling;}
    KEY[*keyNr].method = MaxMpp;
    assert( (KEY[*keyNr].SX==NULL) && (KEY[*keyNr].destructX==NULL));
  } else { 
    /* ***********************************************************
       EXTREMAL GAUSSIAN
     * *********************************************************** */
    InitSimulateRF(x,y,z, dim, lx, grid, covnr, ParamList, nParam, 
		   method, keyNr, &gauss_distr, error);
    if (*error>0) {goto ErrorHandling;}
    KEY[*keyNr].distribution = DISTR_MAXSTABLE;
    
    assert( ((*error<0) && ! ((KEY[*keyNr].SX==NULL) ^ 
			      (KEY[*keyNr].destructX==NULL)))
	    || 
	    ((*error==0) && (((KEY[*keyNr].SX==NULL) && 
			      (KEY[*keyNr].destructX==NULL)))) 
	    );
    if (KEY[*keyNr].SX==NULL) {
      KEY[*keyNr].destructX = extremes_destruct;
      if ((KEY[*keyNr].SX=
	   (extremes_storage*) malloc(sizeof(extremes_storage)))==NULL) {
	*error=ERRORMEMORYALLOCATION; goto InternalErrorHandling;
      }
      es = (extremes_storage*) KEY[*keyNr].SX;
      if ((es->rf =
	   (Real*) malloc(sizeof(Real) *  KEY[*keyNr].totallength))==NULL) {
	*error=ERRORMEMORYALLOCATION; goto InternalErrorHandling;
      }
    } else {es = (extremes_storage*) KEY[*keyNr].SX;}
    
    sigma = sqrt(KEY[*keyNr].param[SILL]);
    if (sigma==0) {*error = ERRORSILLNULL; goto InternalErrorHandling;}
    meanDsigma = KEY[*keyNr].param[MEAN] / sigma;
    es->inv_mean_pos = 
      1.0 / (sigma * INVSQRTTWOPI * exp(-0.5 * meanDsigma * meanDsigma) 
	     + KEY[*keyNr].param[MEAN] * 
#ifdef RF_GSL
	     ??
#else
	     // pnorm(meanDsigma, KEY[*keyNr].param[MEAN], sigma,1,0));
	     pnorm(0.0, KEY[*keyNr].param[MEAN], sigma,1,0));
#endif  
    for (i=0;i<KEY[*keyNr].totallength;i++) es->rf[i]=0.0;
    
    es->assumedmax = EXTREMES_STANDARDMAX * sqrt(KEY[*keyNr].param[SILL]) + 
      ((KEY[*keyNr].param[MEAN]>0) ? KEY[*keyNr].param[MEAN] : 0);  
    *error=NOERROR;
  } 
  return;

 InternalErrorHandling:
  ErrorMessage(Nothing,*error);  
  if (KEY[*keyNr].destructX!=NULL) {
    KEY[*keyNr].destructX(&(KEY[*keyNr].SX));
    KEY[*keyNr].destructX = NULL;
  }
 ErrorHandling:
  return;
}


void DoMaxStableRF(int *keyNr, Real *res, int *error)
{
  extremes_storage *es;  
  long totallength,i,control;
  Real poisson, invpoisson, threshold;
  long counter=0; // just for printing messages
  bool next;

  *error=0; 

  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  if (!KEY[*keyNr].active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}
#ifdef RF_GSL
  assert(RANDOM!=NULL);
#endif


  if (KEY[*keyNr].method == (int) MaxMpp) {
     /* ***********************************************************
       MPP FUNCTIONS
     * *********************************************************** */
    
    // what to do with the mean?? -> ignore it
    // what to do with the variance?? -> ignore it

    int  dim, d;
    Real  min[MAXDIM], max[MAXDIM], factor;
    long segment[MAXDIM+1];
    mpp_storage *bs;
    mppmodel model;
    MPPRandom MppFct;

#ifndef RF_GSL
      GetRNGstate();
#endif
    bs = (mpp_storage*) KEY[*keyNr].S;  
    MppFct = KEY[*keyNr].cov->add_mpp_rnd;
    dim = KEY[*keyNr].dim;
    
    //        printf("OK\n");
    factor = bs->effectivearea / bs->integralpos;
 
    segment[0] = 1;
    for (d=0; d<dim; d++) 
      segment[d+1] = segment[d] * KEY[*keyNr].length[d];
    for  (i=0; i<KEY[*keyNr].totallength; i++) res[i]=0.0;

    control = 0;
    totallength = KEY[*keyNr].totallength;
    poisson = -log(1.0 - UNIFORM_RANDOM); /* it is important to write 1-U
					     as the random generator returns
					     0 inclusively; the latter would
					     lead to an error !
					  */
    invpoisson = 1.0 / poisson;
    //    printf("OK1\n");
    while (control<totallength) {   
      // the following is essentially the same as the code in MPP.cc!!
      // how to reorganise it to get rid of the duplication of the code??
      //
      
      //         printf("OK1a\n");

      MppFct( (mpp_storage*) KEY[*keyNr].S, min, max, &model);

      //           printf("OK2\n");
     if (KEY[*keyNr].grid) {
	int start[MAXDIM], end[MAXDIM], resindex, index[MAXDIM],
	  segmentdelta[MAXDIM], dimM1;
	Real coord[MAXDIM], startcoord[MAXDIM];
	
	// determine rectangle of indices, where the mpp function
	// is different from zero
	for (d=0; d<dim; d++) {	 
	  if (next = ((min[d] > KEY[*keyNr].x[d][XEND]) ||
		       (max[d] < KEY[*keyNr].x[d][XSTART]))) { break;}
	  if (min[d]< KEY[*keyNr].x[d][XSTART]) {start[d]=0;}
	  else start[d] = (int) ((min[d] - KEY[*keyNr].x[d][XSTART]) / 
				 KEY[*keyNr].x[d][XSTEP]);
	  // "end[d] = 1 + *"  since inequalities are "* < end[d]" 
	  // not "* <= end[d]" !
	  end[d] = (int) ((max[d] - KEY[*keyNr].x[d][XSTART]) / 
			      KEY[*keyNr].x[d][XSTEP]);
	  if (end[d] > KEY[*keyNr].length[d]) end[d] = KEY[*keyNr].length[d];	  
	}	
	if (!next) {
	  // prepare coord starting value and the segment vectors for res
	  resindex = 0;
	  for (d=0; d<dim; d++) {
	    index[d]=start[d];
	    segmentdelta[d] = segment[d] * (end[d] - start[d]);
	    resindex += segment[d] * start[d];
	    coord[d] = startcoord[d] = 
	      KEY[*keyNr].x[d][XSTART] + (Real)start[d] * KEY[*keyNr].x[d][XSTEP];
	    //printf("dim=%d: %d %d %d %d %f %d\n",d,start[d],
	    //   segmentdelta[d],segment[d],end[d], coord[d],resindex);
	  }
	  //		   printf("OK3\n");
	  
	  // add mpp function to res
	  dimM1 = dim - 1;
	  d = 0; // only needed for assert(d<dim)
	  while (index[dimM1]<end[dimM1]) {
	    register Real dummy;
	    assert(d<dim); 
	    //	   if (GENERAL_PRINTLEVEL>5) printf("OK4\n");
	    dummy =  model(coord) * invpoisson;
	    //	   if (GENERAL_PRINTLEVEL>5) printf("OK5 %d %d \n",resindex,end[1]);
	    if (res[resindex] < dummy) res[resindex]=dummy;
	    d=0;
	    index[d]++;
	    coord[d]+=KEY[*keyNr].x[d][XSTEP];
	    resindex++;
	    while (index[d] >= end[d]) { 
	      if (d>=dimM1) break;  // below not possible
	      //                       as dim==1 will fail!
	      index[d]=start[d];
	      coord[d]=startcoord[d];
	      resindex -= segmentdelta[d];
	      d++; // if (d>=dim) break;
	      index[d]++;
	      coord[d]+=KEY[*keyNr].x[d][XSTEP];
	      resindex += segment[d];
	    }
	  }
	  //       	printf(" end %d %d\n",control,totallength);////////////////////
	} /* next */ //else {printf("next\n");}
     } else {
       Real y[MAXDIM];
	// the following algorithm can greatly be improved !
	// but for ease, just the stupid algorithm
	for (i=0; i<totallength; i++) {
	  for (d=0; d<dim; d++) { y[d] = KEY[*keyNr].x[d][i]; }
	  register Real dummy;
	  dummy =  model(y) * invpoisson;
	  if (res[i] < dummy) res[i]=dummy;
	}
      }   
     //   printf(" %d %d\n",control,totallength);/////////////////////

      poisson -= log(1.0 - UNIFORM_RANDOM);  
      invpoisson = 1.0 / poisson;
      threshold = bs->maxheight * invpoisson;
      while ((control<totallength) && (res[control]>=threshold)) {control++;}
      if (GENERAL_PRINTLEVEL>=3) {
	counter++;
	if ((counter % 50==0) && (control<totallength))
	  PRINTF("%d: %d %f %f \n",counter,control,res[control],threshold); 
      }      
    }
    for (i=0; i<totallength;i++) res[i] *= factor;
    if (KEY[*keyNr].param[NUGGET]>0) {
      Real nugget;
      nugget = KEY[*keyNr].param[NUGGET];
      for (i=0; i<totallength; i++) {
	register Real dummy;
	dummy = - nugget / log(1.0 - UNIFORM_RANDOM);
	if (res[i] < dummy) res[i]=dummy;
      }
    }
#ifndef RF_GSL
    PutRNGstate();
#endif  
  } else {
     /* ***********************************************************
       EXTREMAL GAUSSIAN
     * *********************************************************** */

    if (KEY[*keyNr].distribution!=DISTR_MAXSTABLE) {
      *error=ERRORWRONGINIT; 
      goto InternalErrorHandling;
    }
    es = (extremes_storage*) KEY[*keyNr].SX;
    totallength = KEY[*keyNr].totallength;
    
    control = 0;
#ifndef RF_GSL
    GetRNGstate();
#endif
    poisson = -log(1.0 - UNIFORM_RANDOM); /* it is important to write 1-U
					     as the random generator returns
					     0 inclusively; the latter would
					     lead to an error !
					  */
#ifndef RF_GSL
    PutRNGstate();
#endif
    invpoisson = 1.0 / poisson;
    while (control<totallength) {
      DoSimulateRF(keyNr,es->rf,error);
      if (*error)  {goto ErrorHandling;}
      if (!KEY[*keyNr].active){
	*error=ERRORNOTINITIALIZED;
	goto InternalErrorHandling;
      }
      for (i=0; i<totallength; i++) {// 0, not control, since we approximate
	//                              and 0 as starting control increases precision
	es->rf[i] *= invpoisson;
	if (res[i] < es->rf[i]) res[i]=es->rf[i];
      }
#ifndef RF_GSL
      GetRNGstate();
#endif
      poisson -= log(1.0 - UNIFORM_RANDOM);  
#ifndef RF_GSL
      PutRNGstate();
#endif
      invpoisson = 1.0 / poisson;
      threshold = es->assumedmax * invpoisson;
      while ((control<totallength) && (res[control]>=threshold)) {control++;}
      
      if (GENERAL_PRINTLEVEL>=3) {
      counter++;
      if (counter % 10==0) 
	PRINTF("%d: %d %f %f \n",counter,control,invpoisson,threshold); 
      }
    }
    for (i=0; i<KEY[*keyNr].totallength; i++) 
      res[i] *= es->inv_mean_pos;
  }
  return;
  
  InternalErrorHandling:
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);  
  ErrorHandling:
  KEY[*keyNr].active=false;  
  if (KEY[*keyNr].destructX!=NULL) {
    KEY[*keyNr].destructX(&(KEY[*keyNr].SX));
    KEY[*keyNr].destructX = NULL;
  }
}


