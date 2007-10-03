
/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 all around the nugget effect -- needs special treatment 

 Copyright (C) 2001 -- 2006 Martin Schlather, 

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

double NUGGET_TOL=0.0;
void SetParamNugget(int *action, double *nuggettol)
{
  if (*action) {
    NUGGET_TOL = *nuggettol;
    if (NUGGET_TOL < 0) {
      if (GENERAL_PRINTLEVEL>=1) 
	PRINTF("negative tolerance for distance in nugget covariance not allowed; set to zero");
      NUGGET_TOL = 0.0;
    }
  } else {
    *nuggettol = NUGGET_TOL;
  }
}


void nugget_destruct(void ** S)
{
  if (*S != NULL) {
    nugget_storage *x;
    x = *((nugget_storage**)S); 
    if (x->pos!=NULL) free(x->pos);
    if (x->red_field!=NULL) free(x->red_field);
    free(*S);
    *S = NULL;
  }
}

bool equal(int i, int j, double *X, int dim)
{
  param_type p;
  double *x, *y;
  register double dummy, dist;
  int d;
  x = X + i * dim;
  y = X + j * dim;
  for (dist=0.0, d=0; d<dim; d++) {
    dummy = x[d]-y[d];
    dist += dummy * dummy;
  }
  dist = sqrt(dist);
  return nugget(&dist, p)==1.0;
}

// uses global RANDOM !!!
int init_nugget(key_type *key, int m){ 
  methodvalue_type *meth; 
  int Xerror, i, nonzero_pos, v, actcov; 
  nugget_storage *s;
  double nugget_effect;
  covinfo_type *kc, *first=NULL;

  if (key->covFct != CovFct) { Xerror=ERRORNOTPROGRAMMED; goto ErrorHandling;}
  meth = &(key->meth[m]);
  nonzero_pos = -1;
  SET_DESTRUCT(nugget_destruct, m);
  
  actcov=0;
  for (v=0; v<key->ncov; v++) {
    ERRORMODELNUMBER = v;	
    kc = &(key->cov[v]);
    cov_fct *cov;

//    printf("nugget %d %d %d %s \n", v, kc->nr, key->ncov, CovList[kc->nr].name);
    assert((kc->nr >= 0) && (kc->nr < currentNrCov));
    assert(kc->param[VARIANCE] >= 0.0);
    cov = &(CovList[kc->nr]);
    

    if ((kc->method==Nugget) && (kc->left)) {
      meth->covlist[actcov] = v;
      if (cov->type==ISOHYPERMODEL|| cov->type==ANISOHYPERMODEL) {
	  Xerror=ERRORHYPERNOTALLOWED; goto ErrorHandling; }
      if (cov->implemented[Nugget] != IMPLEMENTED) {
	  Xerror=ERRORNOTDEFINED; goto ErrorHandling;} // schliesst hyper mit aus
      assert(cov->cov==nugget);
      if (actcov>0) { 
        // actcov>0 not v>0 since next line actcov-1 -- check not for all v
	if (key->cov[v-1].op) {
	  if (key->cov[v-1].method != Nugget) { // v-1 !!
	    if (GENERAL_PRINTLEVEL>0) 
	       PRINTF("severe error - contact author. %d %d %d %d (%s) %d (%s)n",
		      v, key->cov[v-1].op, key->ncov, key->cov[v-1].method,
		      METHODNAMES[key->cov[v-1].method],
		      Nugget, METHODNAMES[Nugget]);
	    Xerror=ERRORMETHODMIX; goto ErrorHandling;
	  }
	  Xerror=ERRORNOMULTIPLICATION; goto ErrorHandling;
	}
      }
      if (key->anisotropy) {
	bool equal;
	int w, endfor;
	endfor = key->totalparam;
	// if (key->Time && TIMEEXCEPTION) endfor -= key->timespacedim;
	double quotient, f_epsilon = 1e-15;
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
	  quotient = kc->param[nonzero_pos]/first->param[nonzero_pos];
	  equal = quotient != 0.0;
	  for (w=ANISO; w<endfor; w++) {
	    equal &= fabs(first->param[w] * quotient - kc->param[w]) <
              (fabs(kc->param[w]) + 1.0*(kc->param[w]==0.0))
	      * f_epsilon;
	  }
	  if (!equal) { /* even not equal up to a multiplicative constant */
	    if (key->cov[v-1].op) { Xerror=ERRORANISOMIX; goto ErrorHandling; }
	    else { /* ???? nur additive ********* */
		continue;
	    }
	  }
	} else {
	  // needed to determine nonzero_pos above
	  first = kc;
	  nonzero_pos=ANISO;
	  while ((nonzero_pos<endfor) && (first->param[nonzero_pos]==0.0))
	      nonzero_pos++;
	  if (nonzero_pos>=endfor) {
	      Xerror=ERRORLOWRANK; goto ErrorHandling; 
	  }
	}
      }
      actcov++;
      kc->left = false;
    }
  } // v
  ERRORMODELNUMBER = -1;	
  if (actcov==0) { /* no covariance for the considered method found */
    Xerror=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  }
  meth->actcov = actcov;  

  if ((meth->S=malloc(sizeof(mpp_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (nugget_storage*) meth->S;
  s->pos = NULL;
  s->red_field = NULL;
  kc = &(key->cov[meth->covlist[0]]);
  s->simple = key->timespacedim == kc->reduceddim; // fkt, da FULLISOTROPIC
  s->simugrid = kc->simugrid;
  if (key->anisotropy && !s->simple) {
    int *pos, oldpos;
    if (NUGGET_TOL==0.0 && GENERAL_PRINTLEVEL>=1)
      PRINTF("\nThe anisotropy matrix does not have full rank and RFparameters()$nugget.tol equals 0. From a theoretical point of view that's fine, but the simulations will probably be odd. Is this really what you want?\n");
    if ((Xerror=Transform2NoGrid(key, meth->covlist[0])) != NOERROR) 
      goto ErrorHandling;
    if (kc->simugrid) {
      int d, dim=key->timespacedim + 1;
      s->prod_dim[0] = 1; 
      for (d=i=0; d<key->timespacedim; d++, i+=dim) {
        s->reduced_dim[d] =(fabs(kc->param[ANISO + i]) 
			    < NUGGET_TOL)? 1 : key->length[d];
	s->prod_dim[d + 1] = s->prod_dim[d] * s->reduced_dim[d];
      }
      if ((s->red_field=(double *) malloc(sizeof(double) * 
					  s->prod_dim[key->timespacedim])) 
	  == NULL){
	  Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
    } else {
      if ((pos = (int*) malloc(sizeof(int) * key->totalpoints))==0) {
	Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      ordering(kc->x, key->totalpoints, kc->reduceddim, pos);
      oldpos = pos[0];
      for (i=1 /* ! */; i<key->totalpoints; i++) {
	if (equal(oldpos, pos[i], kc->x, kc->reduceddim)) 
	  pos[i]= -1 - pos[i];
	else oldpos=pos[i];
      }
    s->pos = pos;
    }
  }
  nugget_effect = CovFct(ZERO, 1, key->cov, key->meth[m].covlist, actcov, false);
  s->sqrtnugget = sqrt(nugget_effect);
  
  if (key->anisotropy) return NOERROR_REPEAT;
  else return NOERROR;

 ErrorHandling:
  return Xerror;
}

void do_nugget(key_type *key, int m, double *res) {
  methodvalue_type *meth; 
  double sqrtnugget;
  nugget_storage* s;
  long nx, endfor;

  meth = &(key->meth[m]);
  assert(meth->S != NULL);
  s = (nugget_storage*) meth->S;
  sqrtnugget = s->sqrtnugget;

  if (s->simple) {
    for (nx=0, endfor=key->totalpoints; nx<endfor; nx++)
      res[nx] += GAUSS_RANDOM(sqrtnugget);
  } else {
    if (s->simugrid) {
      int d, i, dim, dimM1, index[MAXDIM], *red_dim, *prod_dim;
      long totpnts, idx;
      double *field;
      totpnts = key->totalpoints;
      dim = key->timespacedim;
      dimM1 = dim -1;
      field = s->red_field;
      red_dim = s->reduced_dim;
      prod_dim = s->prod_dim;

      for (i=s->prod_dim[dim] - 1; i>=0; field[i--] = GAUSS_RANDOM(sqrtnugget));
      for (d=0; d<dim; index[d++] = 0);
      for (i=0; i<totpnts; i++) {
//	  printf("i=%d %d %d \n", i, totpnts, key->timespacedim);
	for(idx=d=0; d<dim; d++) idx += (index[d] % red_dim[d]) * prod_dim[d];
	//printf("%d %d %f\n", i, idx, field[idx]);
	res[i] += field[idx];
	d = 0; 
	(index[d])++; 
	while (d < dimM1 && index[d] >= key->length[d]) { 
	  assert(d<dim); // assert(d>=0);
	  index[d] = 0; 
	  d++;   
	  (index[d])++; 
	}
      }
    } else {
      int p;
      double dummy = RF_NAN; // just to avoid warnings from the compiler
      assert(s->pos[0]>=0);
      for (nx=0; nx<key->totalpoints; nx++) {
	if ((p = s->pos[nx]) < 0) p= -1 - p; // if p<0 then take old variable
	// and -p-1 is the true index 
	else dummy = GAUSS_RANDOM(sqrtnugget);
	res[p] += dummy;
      }
    }
  }
}

