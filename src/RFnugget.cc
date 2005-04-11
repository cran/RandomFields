
/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 all around the nugget effect -- needs special treatment 

 Copyright (C) 2001 -- 2004 Martin Schlather, 

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
#include <stdio.h>  
#include <stdlib.h>
#include <assert.h>
#include "RFsimu.h"
#include "RFCovFcts.h"


static double NUGGET_TOL_SINGULAR=1e-12;
static double NUGGET_TOL_POINTS=1e-12;

typedef struct nugget_storage {
  double sqrtnugget;
  bool simple;
  int *pos;
} nugget_storage;

void nugget_destruct(void ** S)
{
  if (*S != NULL) {
    nugget_storage *x;
    x = *((nugget_storage**)S); 
    if (x->pos!=NULL) free(x->pos);
    free(*S);
    *S = NULL;
  }
}



/* nugget effect model */
int is_singular(double *C, int dim, bool *singular)
{
  double D[MAXDIM];
  int Xerror, d;
  if ((Xerror=eigenvalues(C, dim, D))!=NOERROR) goto ErrorHandling;
  *singular = false;
  for (d=0; d<dim; d++) 
    if (fabs(D[d])<NUGGET_TOL_SINGULAR) {*singular=true; break;}
  return 0;

 ErrorHandling:
  return Xerror;
}

bool equal(int i, int j, double *ORDERD, int ORDERDIM)
{
  double *x, *y;
  int d;
  x = ORDERD + i * ORDERDIM;
  y = ORDERD + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++) {
      if(fabs(x[d]-y[d]) > fabs(x[d]>y[d] ? x[d] : y[d]) * NUGGET_TOL_POINTS) {
	return false;
      }
  }
  return true;
}

    


// uses global RANDOM !!!
int init_nugget(key_type *key, int m){ 
  int Xerror; 
  nugget_storage *s;
  param_type param;
  double nugget_effect, quotient[MAXCOV];
  unsigned short int actcov;
  int i,covnr[MAXCOV], nonzero_pos;
  int multiply[MAXCOV], TrueDim, v;
  bool singular, no_last_comp;
  double *xx;

  xx=NULL;
  nonzero_pos = -1;
  SET_DESTRUCT(nugget_destruct);
  
  double store_param[TOTAL_PARAM];
  long bytes;
  /* if (key->Time) bytes = key->timespacedim * (key->timespacedim-1); */
  /* else */
  bytes = sizeof(double) * key->timespacedim * key->timespacedim;

  actcov=0;
  for (v=0; v<key->ncov; v++) {
    if ((key->method[v]==Nugget) && (key->left[v])) {
      //Variance==0 is not eliminated anymore !!! -- maybe this could be improved
      //	&& (key->param[v][VARIANCE]>0)) {   
      // do not remove the parenths around METHOD!
      key->left[v] = false;
      assert((key->covnr[v] >= 0) && (key->covnr[v] < currentNrCov));
      assert(key->param[v][VARIANCE] >= 0.0);
      covnr[actcov] = key->covnr[v];
      if (CovList[covnr[actcov]].implemented[Nugget] != IMPLEMENTED) {
      //  formerly:   if (CovList[covnr[actcov]].cov==NULL) {
	Xerror=ERRORNOTDEFINED; goto ErrorHandling;}
      memcpy(param[actcov], key->param[v], sizeof(double) * key->totalparam);
      if (actcov>0) { 
        // actcov>0 not v>0 since next line actcov-1 -- check not for all v
	if ((multiply[actcov-1] = key->op[v-1]) && key->method[v-1] != Nugget){
          if (GENERAL_PRINTLEVEL>0) 
	    PRINTF("severe error - contact author. %d %d %d %d (%s) %d (%s)n",
		   v, key->op[v-1], key->ncov, key->method[v-1],
		   METHODNAMES[key->method[v-1]],
		   Nugget, METHODNAMES[Nugget]);
          Xerror=ERRORMETHODMIX; goto ErrorHandling;
	}
      }

      if (key->anisotropy) {
	bool equal;
	int w, endfor;
	endfor = key->totalparam;
	// if (key->Time && TIMEEXCEPTION) endfor -= key->timespacedim;
	double f_epsilon = 1e-15;
	if (actcov>0) {
	  /* check whether the multiplicative construction are consistent, 
	     i.e. that the SPATIAL part of the anisotropy matrices are 
	     multiplicatives of each other (more precisely, the remaining ones
	     are multiplicatives of the first.) 
	     The time part of the matrix must be of the form (0,..,0,*).
	     A further goal of this part is the collection of additive blocks 
	     that have the same anisotropy matrix structure, in order to minimise
	     the simulation time */
	  assert(nonzero_pos>0);
	  quotient[actcov] = param[actcov][nonzero_pos]/store_param[nonzero_pos];
	  equal = true;
	  for (w=ANISO; w<endfor; w++) {
	    equal &= fabs(store_param[w] * quotient[actcov]
		      - param[actcov][w]) <
              (fabs(param[actcov][w]) + 1.0*(param[actcov][w]==0.0)) *f_epsilon;
	  }
	  if (!equal) { /* even not equal up to a multiplicative constant */
	    if (multiply[actcov-1]) { Xerror=ERRORANISOMIX; goto ErrorHandling; }
	    else { /* ???? nur additive ********* */
	      key->left[v]=true;
	      actcov--;
	    }
	  }
	} else {
	  memcpy(&(store_param[ANISO]), &(param[actcov][ANISO]), bytes);
	  nonzero_pos=ANISO;
	  quotient[0] = 1.0;
	  while ((nonzero_pos<endfor) && (param[actcov][nonzero_pos]==0))
	    nonzero_pos++;
	  if (nonzero_pos>=endfor) { Xerror=ERRORTRIVIAL; goto ErrorHandling; }
	}
      } else {
	assert(fabs(key->param[v][SCALE] * key->param[v][INVSCALE]-1.0)
	       < EPSILON);
      }
      actcov++;
    }
  }
  if (actcov==0) { /* no covariance for the considered method found */
    Xerror=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  }

  if ((key->S[m]=malloc(sizeof(mpp_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (nugget_storage*) key->S[m];
  s->pos = NULL;

  if (key->anisotropy) {
    if ((Xerror = is_singular(&(param[0][ANISO]), key->timespacedim, &singular))
	  !=NOERROR) 
      goto ErrorHandling;
    s->simple = !singular;
    if (singular) {
      int *pos, oldpos, start_param[MAXDIM], index_dim[MAXDIM];
      GetTrueDim(key->anisotropy, key->timespacedim, param[0],
		 &TrueDim, &no_last_comp, start_param, index_dim);
      if ((Xerror=Transform2NoGrid(key, param[0], TrueDim, start_param, &(xx)))
	  !=NOERROR)
	goto ErrorHandling;

      if ((pos = (int*) malloc(sizeof(int) * key->totalpoints))==0){
	Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      ordering(xx, key->totalpoints, TrueDim, pos);
      oldpos = pos[0]; 
      for (i=1 /* ! */; i<key->totalpoints; i++) {
	if (equal(oldpos, pos[i], xx, TrueDim)) pos[i]= -1 - pos[i];
	else oldpos=pos[i];
      }
      s->pos = pos;
    }
  } else {
    s->simple = true;
  }
  for (i=0; i<actcov; i++) assert(CovList[covnr[i]].cov==nugget);
  nugget_effect = CovFct(ZERO,1, covnr, multiply, param, actcov, false);
  s->sqrtnugget = sqrt(nugget_effect);
  if (xx!=NULL) free(NULL);
  
  if (key->anisotropy) return NOERROR_REPEAT;
  else return 0;

 ErrorHandling:
  if (xx!=NULL) free(NULL);
  return Xerror;
}

void do_nugget(key_type *key, bool add, int m, double *res ) {
  double sqrtnugget;
  nugget_storage* s;
  long nx, endfor;

#define RESULT(OP)\
  for (nx=0,endfor=key->totalpoints; nx<endfor; nx++)\
     res[nx] OP GAUSS_RANDOM(sqrtnugget);\

  assert(key->S[m] != NULL);
  s = (nugget_storage*) key->S[m];
  sqrtnugget = s->sqrtnugget;

  if (s->simple) {
    if (add) { RESULT(+=) } else { RESULT(=) } 
  } else {
    int p;
    double dummy;
    assert(s->pos[0]>=0);
    dummy = RF_NAN;
    for (nx=0; nx<key->totalpoints; nx++) {
      if ((p=s->pos[nx])<0) p= -1 - p; // if p<0 then take old variable
      // and -p-1 is the true index 
      else dummy = GAUSS_RANDOM(sqrtnugget);
      if (add) res[p] += dummy;
      else res[p] = dummy;  
    }
  }
}

