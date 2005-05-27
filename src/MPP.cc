

/*
 Authors 
 Martin Schlather, schlath@hsu-hh.de 

 Copyright (C) 2001 -- 2005 Martin Schlather, 


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

/* TO DO
  * noch fast alles auf key->x programmiert statt s->x
  * time component: aenderungen siehe TBM !
*/

#include <math.h>
#include <assert.h>
#include "RFsimu.h"
#include "MPPFcts.h"

double MPP_RADIUS =0.0;
double MPP_APPROXZERO  =0.001;
double ADDMPP_REALISATIONS = 100.0; // logically an integer, implementation 
                                  //  allows for positive real ! 

void SetMPP(int *action, double *approxzero, double *realisations, double *radius)
{
  switch(*action) {
  case 0 :
    if (*approxzero<=0.0) {
       if (GENERAL_PRINTLEVEL>0) 
	 PRINTF("\nERROR! `approx.zero' not a positive number. Ignored.\n");
    } else  MPP_APPROXZERO = *approxzero;
   if (*realisations<=0.0) {
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nERROR! `realisation' not a positive number. Ignored.\n");
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

void mpp_destruct(void **S) {
  // mpp_storage in RFsimu.h !!
  if (*S!=NULL) {
    free(*S);
    *S = NULL;
  }
}


int init_mpp(key_type * key, int m)
{
  int error, d, i, v, start_param[MAXDIM], index_dim[MAXDIM];
  bool no_last_comp /* , any_time_exception*/;
  double max[MAXDIM];  
  mpp_storage *s;
  cov_fct *cov;
  unsigned short int actcov;
  int covnr[MAXCOV];

  if (key->anisotropy) {
    error=ERRORNOTPROGRAMMED; goto ErrorHandling; 
    // RFtest.new.features.1.01.R: MaxStableRF is excluded by if (FALSE)
    // excludes also any time component
    //
    // even in the case of isotropy the method is not checked yet!
  }
			 
  SET_DESTRUCT(mpp_destruct);
  assert(key->S[m]==NULL);
  if ((key->S[m]=malloc(sizeof(mpp_storage)))==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (mpp_storage*) key->S[m];

  actcov=0;
  { 
    int v;
    for (v=0; v<key->ncov; v++) {
      if ((key->method[v]==AdditiveMpp) && (key->left[v])
	  && (key->param[v][VARIANCE]>0)) {
	key->left[v]=false;
	assert(key->covnr[v]>=0 && key->covnr[v]<currentNrCov);
	covnr[actcov] = key->covnr[v];
	if (CovList[covnr[actcov]].add_mpp_scl==NULL)
	  {error=ERRORNOTDEFINED; goto ErrorHandling;}
	memcpy(s->param[actcov], key->param[v], sizeof(double) * key->totalparam);
	if ((v<key->ncov-1 && key->op[v]) || (v>0 && key->op[v-1])) { 
	  // no multiplicative models are allowed
	  error=ERRORNOMULTIPLICATION; goto ErrorHandling; }
	actcov++; 
	break;
      }
    }
  }
  if (actcov==0) {
    error=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  }
  
  s->actcov = actcov;

  // transform points using anisotropy and get s->timespacedim
  //if ((param[0][key->totalparam-1]==0.0) && any_time_exception) {
  //  param[0][key->totalparam-1]=1.0;
  //} // see tbm

  s->grid = key->grid && !key->anisotropy;
  GetTrueDim(key->anisotropy, key->timespacedim, s->param[0], 
	     &s->timespacedim, &no_last_comp, start_param, index_dim);
  if ((error=Transform2NoGrid(key, s->param[0], s->timespacedim, 
			      start_param, &(s->x)))!=NOERROR)
      goto ErrorHandling; 
  

  // are methods and parameters fine ?
  for (v=0; v<actcov; v++) {
    cov = &(CovList[covnr[v]]);
    if ((key->Time) && no_last_comp && (cov->isotropic==SPACEISOTROPIC))
      { error= ERRORWITHOUTTIME;goto ErrorHandling;}
    else if ((cov->check!=NULL) &&
	     (error=cov->check(s->param[v], s->timespacedim, AdditiveMpp)))
      goto ErrorHandling;
  }
 
  // determine the minimum area where random field values are to be generated
  if (s->grid) {
    for (d=0; d<key->timespacedim; d++) {
      s->min[d]   = key->x[d][XSTART];
      s->length[d]= key->x[d][XSTEP] * (double) (key->length[d] - 1);
      // x[XSTART] and x[XSTEP] are assumed to be precise, but
      // not x[XEND] 
    }
  } else {
    int ix;
    for (d=0;d<MAXDIM;d++) {s->min[d]=RF_INF; max[d]=RF_NEGINF;}
    for (ix=i=0;i<key->totalpoints;i++,ix+=s->timespacedim) {
      for (d=0; d<s->timespacedim; d++) {
	if (s->x[ix+d] < s->min[d]) s->min[d] = s->x[ix+d];
	if (s->x[ix+d] > max[d]) max[d] = s->x[ix+d];
      }
    }     
    for (d=0;d<MAXDIM;d++) { s->length[d] = max[d] - s->min[d]; }
    for (v=0; v<actcov; v++) {
      if (CovList[covnr[v]].isotropic == FULLISOTROPIC)
	s->param[v][SCALE]=s->param[v][INVSCALE]=1.0;
      else assert(false); // no anisotropic function programmed yet;
      //                     neither considered what the extension should
      //                     look like
    }
  }

  for (i=0; i<actcov; i++) {
    s->addradius[i] = MPP_RADIUS;  // must be set BEFORE the next command!!
    CovList[covnr[i]].add_mpp_scl(s, i); 
    s->MppFct[i] = CovList[covnr[i]].add_mpp_rnd;
  }
  if ((MPP_RADIUS>0.0) && (GENERAL_PRINTLEVEL>=2))
    PRINTF("Note: area has been enlarged by fixed value (%f)", s->addradius[0]);
 
  if (key->anisotropy) return NOERROR_REPEAT;
  else return 0;

 ErrorHandling:
  return error; 
}


void do_addmpp(key_type *key, int m, double *res )
{ 
  int i, d, v;
  double lambda[MAXDIM], min[MAXDIM], max[MAXDIM], 
    factor, average, poisson;
  long segment[MAXDIM+1];
  mpp_storage *s;
  mppmodel model;

  assert(key->active);
  {
    register long i,endfor;
    endfor = key->totalpoints;
    for (i=0;i<endfor;i++) { res[i]=0.0; }
  }

  s =  (mpp_storage*) key->S[m];  
  
  switch (key->distribution) {
  case DISTR_GAUSS : 
    assert(s->actcov==1);
    lambda[0] = ADDMPP_REALISATIONS * s->effectivearea[0];
    factor = sqrt(s->param[0][VARIANCE] / 
		  (ADDMPP_REALISATIONS *s->integralsq[0]));
    average = ADDMPP_REALISATIONS * s->integral[0];
    break;
  case DISTR_POISSON :
    // not done yet -- it is only a fragment
    // on the R level, it is excluded that the program runs 
    // into this part
    //
    // actually, only actcov==1 is possible at the moment!
    //
    // maybe nice to switch MPPFcts back to actcov=1
    assert(s->actcov==1);
    for (v=0; v<s->actcov; v++) {
      lambda[v] = s->param[v][VARIANCE] * s->effectivearea[v]; 
      factor = average = 0.0 /* replace by true value !! */;  
      assert(false);
      if(key->mean!=0) {
	// do not compare mean against VARIANCE since
        // VARIANCEs may differ with v
	//      if (GENERAL_PRINTLEVEL>0)
	ErrorMessage(AdditiveMpp,ERRORVARMEAN); // print it always

	assert(false); // for debugging -- delete this later on

	// finally, the value of MEAN is ignored, and that of VARIANCE
	// is used for simulating the random field
      }
    }
    break;
  default : assert(false);
  }
  
  segment[0] = 1; // only used for s->grid!
  for (d=0; d<key->timespacedim; d++) { 
    assert(d<MAXDIM);
    segment[d+1] = segment[d] * key->length[d];
  }
  for  (i=0; i<key->totalpoints; i++) res[i]=0.0;
  for (v=0; v<s->actcov; v++) {
    for(poisson = -log(1.0 - UNIFORM_RANDOM); 
	poisson < lambda[v];
	poisson -= log(1.0 - UNIFORM_RANDOM)) 
      {   
	(s->MppFct[v])(s, v, min, max, &model );
	if (s->grid) {
	  int start[MAXDIM], end[MAXDIM], resindex, index[MAXDIM],
	    segmentdelta[MAXDIM], dimM1;
	  double coord[MAXDIM], startcoord[MAXDIM];
	  
	  // determine rectangle of indices, where the mpp function
	  // is different from zero
	  for (d=0; d<s->timespacedim; d++) {	 
	    if (min[d]< key->x[d][XSTART]) {start[d]=0;}
	    else start[d] = (int) ( (min[d] - key->x[d][XSTART]) / 
				    key->x[d][XSTEP]);
	    // "end[d] = 1 + *"  since inequalities are "* < end[d]" 
	    // not "* <= end[d]" !
	    end[d] = 1 + (int) ((max[d] - key->x[d][XSTART]) /
				key->x[d][XSTEP]);
	    if (end[d] > key->length[d]) end[d] = key->length[d];	  
	  }	
	  
	  // prepare coord starting value and the segment vectors for res
	  resindex = 0;
	  for (d=0; d<s->timespacedim; d++) {
	    index[d]=start[d];
	    segmentdelta[d] = segment[d] * (end[d]-start[d]);
	    resindex += segment[d] * start[d];
	    coord[d] = startcoord[d] = 
	      key->x[d][XSTART] + (double)start[d] * key->x[d][XSTEP];
	  }
	  
	  // add mpp function to res
	  dimM1 = s->timespacedim -1;
	  d = 0; // only needed for assert(d<dim)
	  while (index[dimM1]<end[dimM1]) {
	    assert(d<s->timespacedim); 
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
	} else { // !s->grid
	  double y[MAXDIM];
	  long j;
	  // the following algorithm can greatly be improved !
	  // but for ease, just the stupid algorithm
	  for (j=i=0; i<key->totalpoints; i++) {
	    for (d=0; d<s->timespacedim; d++,j++) { y[d] = s->x[j]; }
	    res[i] += model(y); 
	  }
	}      
      }
  }
  if (key->distribution==DISTR_GAUSS) {
    for (i=0; i<key->totalpoints;i++) {
      res[i] = factor * (res[i] - average);
    }
  }
}

  










