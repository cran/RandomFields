/*
 Authors 
 Martin Schlather, Martin.Schlather@uni-bayreuth.de 

 *********** this function is still under construction *********

 Copyright (C) 2001 Martin Schlather, 


This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
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
//#include "MathlibCC.h"
#include <R_ext/Mathlib.h>
#include <assert.h>
#include "RFsimu.h"
#include "MPPFcts.h"
//#include "RFCovFcts.h"
#include <assert.h>

Real MPP_RADIUS =0.0;
Real MPP_APPROXZERO  =0.001;
Real ADDMPP_REALISATIONS = 100.0; /* logically an integer, implementation allows for
				      positive real ! */

void SetMPP(int *action, Real *approxzero, Real *realisations, Real *radius)
{
  switch(*action) {
  case 0 :
    if (*approxzero<=0.0) {
       if (GENERAL_PRINTLEVEL>0) PRINTF("\nERROR! `approx.zero' not a positive number. Ignored.\n");
    } else  MPP_APPROXZERO = *approxzero;
    if (*realisations<=0.0) {
      if (GENERAL_PRINTLEVEL>0) PRINTF("\nERROR! `realisation' not a positive number. Ignored.\n");
    } else ADDMPP_REALISATIONS = *realisations;
    MPP_RADIUS = *radius;
    if ((MPP_RADIUS>0.0) && (GENERAL_PRINTLEVEL>1))
      PRINTF("\nWARNING! Positive `addradius' may have bad simulation results as consequence, difficult to detect.\n");
    break;
  case 1 :
    *approxzero = MPP_APPROXZERO;
    *realisations =  ADDMPP_REALISATIONS;
    *radius = MPP_RADIUS;
    if (GetNotPrint) break;
  case 2 :
    PRINTF("\nmarked point process methods\n============================\napproxzero=%e,\nrealisations=%e\n",
	   MPP_APPROXZERO,ADDMPP_REALISATIONS);
    break;
  default : if (GENERAL_PRINTLEVEL>0) PRINTF("unknown action\n");
  }
}

void printMPPinfo(int *keyNr)
{    
   mpp_storage *bs;
   if (((KEY[*keyNr].method!=AdditiveMpp) && (KEY[*keyNr].method!=MaxMpp))
       || (!KEY[*keyNr].active)) {
    PRINTF("no information available\n");
    return;
   }
   assert(KEY[*keyNr].S!=NULL);
   bs = (mpp_storage*) KEY[*keyNr].S;
   PRINTF("integral=%e\n",bs->integral);
   PRINTF("integralsq=%e\n",bs->integralsq);
   PRINTF("integralpos=%e\n",bs->integralpos);
   PRINTF("eff.Radius=%e\n",bs->effectiveRadius);
   PRINTF("eff.area=%e\n",bs->effectivearea);
   PRINTF("maxheight=%e\n",bs->maxheight);
   PRINTF("add.Radius=%e\n",bs->addradius);
}

void GetMPPInfo(int *keyNr, Real *integral, Real *integralsq,
		  Real *integralpos, Real *maxheight)
{
  mpp_storage *bs;
  if (((KEY[*keyNr].method!=AdditiveMpp) && (KEY[*keyNr].method!=MaxMpp))
      || (!KEY[*keyNr].active)) {
    *integral = *integralsq = *integralpos = *maxheight = RF_NAN;
    return;
  }
  assert(KEY[*keyNr].S!=NULL);
  bs = (mpp_storage*) KEY[*keyNr].S;
  *integral = bs->integral;
  *integralsq = bs->integralsq;
  *integralpos = bs->integralpos;
  *maxheight = bs->maxheight;
}

void mpp_destruct(void **S) {
  // mpp_storage in RFsimu.h !!
  if (*S!=NULL) {
    free(*S);
    *S = NULL;
  }
}

int init_mpp(key_type * key)
{
  int error, d, i;
  Real max[MAXDIM];  
  mpp_storage *bs;

  assert(key->active==false);
  assert((key->covnr>=0) &&(key->covnr<currentNrCov));
  //assert(key->cov->cov!=NULL); here is the first !! 
  //                              where cov is unknown, but should be simulated
  assert(key->param[VARIANCE]>=0.0);
  assert(key->param[SILL]>=key->param[VARIANCE]); 
  assert( fabs(key->param[SCALE]*key->param[INVSCALE]-1.0) < EPSILON);
  assert((key->S==NULL) && (key->destruct==NULL));
  assert((key->cov->add_mpp_scl!=NULL) && (key->cov->add_mpp_rnd!=NULL));

  if ((key->S=malloc(sizeof(mpp_storage)))==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  key->destruct = mpp_destruct;
  
  // determine the minimum area where random field values are to be generated
  bs = (mpp_storage*) key->S;
  if (key->grid) {
    for (d=0; d<key->dim; d++) {
      bs->min[d]   = key->x[d][XSTART];
      bs->length[d]= key->x[d][XSTEP] * (Real) (key->length[d] - 1);
    }
  } else {
    for (d=0;d<MAXDIM;d++) {bs->min[d]=RF_INF; max[d]=RF_NEGINF;}
    for (i=0;i<key->totallength;i++) {
      for (d=0;d<key->dim;d++) {
	if (key->x[d][i]<bs->min[d]) bs->min[d]=key->x[d][i];
	if (key->x[d][i]>max[d]) max[d]=key->x[d][i];
      }
    }     
    for (d=0;d<MAXDIM;d++) { bs->length[d] = max[d] - bs->min[d]; }
  }
  bs->dim = key->dim;
  bs -> addradius = MPP_RADIUS;
  if ((MPP_RADIUS>0.0) && (GENERAL_PRINTLEVEL>=2))
    PRINTF("Note: area has been enlarged by fixed value (%f)",bs->addradius);
  key->cov->add_mpp_scl(key);
  return NOERROR;
  
 ErrorHandling:
  if (key->destruct!=NULL) key->destruct(&key->S);
  key->destruct=NULL;
  key->active = false;
  return error; 
}


void do_addmpp(key_type *key, Real *res END_WITH_GSLRNG)
{ 
  int i, d, dim;
  Real lambda, min[MAXDIM], max[MAXDIM], factor, average, poisson;
  long segment[MAXDIM+1];
  mpp_storage *bs;
  mppmodel model;
  MPPRandom MppFct;

  assert(key->active);
#ifdef RF_GSL
  assert(RANDOM!=NULL);
#endif

  bs = (mpp_storage*) key->S;  
  MppFct = key->cov->add_mpp_rnd;
  dim = key->dim;
  switch (key->distribution) {
  case DISTR_GAUSS : 
    lambda = ADDMPP_REALISATIONS * bs->effectivearea;
    factor = 
      sqrt(key->param[VARIANCE] / (ADDMPP_REALISATIONS * bs->integralsq));
    average = ADDMPP_REALISATIONS * bs->integral;
    break;
  case DISTR_POISSON :
    lambda = key->param[VARIANCE] * bs->effectivearea; 
    factor = average = 0 /* replace by true value !! */;  assert(false);
    if(key->param[VARIANCE]!=key->param[MEAN]) {
      //      if (GENERAL_PRINTLEVEL>0)
      ErrorMessage(AdditiveMpp,ERRORVARMEAN); // print it always

      assert(false); // for debugging -- delete this later on

      // finally, the value of MEAN is ignored, and that of VARIANCE
      // is used for simulating the random field
    }
    break;
  default : assert(false);
  }
  
  segment[0] = 1;
  for (d=0; d<dim; d++) 
    segment[d+1] = segment[d] * key->length[d];
  for  (i=0; i<key->totallength; i++) res[i]=0.0;
  for(poisson = -log(1.0 - UNIFORM_RANDOM); 
      poisson < lambda;
      poisson -= log(1.0 - UNIFORM_RANDOM)) 
    {   
      MppFct( (mpp_storage*) key->S, min, max, &model);
      if (key->grid) {
	int start[MAXDIM], end[MAXDIM], resindex, index[MAXDIM],
	  segmentdelta[MAXDIM], dimM1;
	Real coord[MAXDIM], startcoord[MAXDIM];
	 
	// determine rectangle of indices, where the mpp function
	// is different from zero
	for (d=0; d<dim; d++) {	 
	  if (min[d]< key->x[d][XSTART]) {start[d]=0;}
	  else start[d] = (int) ( (min[d] - key->x[d][XSTART]) / key->x[d][XSTEP]);
	  // "end[d] = 1 + *"  since inequalities are "* < end[d]" 
	  // not "* <= end[d]" !
	  end[d] = 1 + (int) ((max[d] - key->x[d][XSTART]) / key->x[d][XSTEP]);
	  if (end[d] > key->length[d]) end[d] = key->length[d];	  
	}	

	// prepare coord starting value and the segment vectors for res
	resindex = 0;
	for (d=0; d<dim; d++) {
	  index[d]=start[d];
	  segmentdelta[d] = segment[d] * (end[d]-start[d]);
	  resindex += segment[d] * start[d];
	  coord[d] = startcoord[d] = 
	    key->x[d][XSTART] + (Real)start[d] * key->x[d][XSTEP];
	}

	// add mpp function to res
	dimM1 = dim -1;
	d = 0; // only needed for assert(d<dim)
	while (index[dimM1]<end[dimM1]) {
	  assert(d<dim); 
	  res[resindex] += model(coord);
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
      } else {
	Real y[MAXDIM];
	// the following algorithm can greatly be improved !
	// but for ease, just the stupid algorithm
	for (i=0; i<key->totallength; i++) {
	  for (d=0; d<dim; d++) { y[d] = key->x[d][i]; }
	  res[i] += model(y);
	}
      }      
    }
  if (key->distribution==DISTR_GAUSS) {
    for (i=0; i<key->totallength;i++) 
      res[i] = factor * (res[i] - average);
  }    
}












