
/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 all around the nugget effect -- needs special treatment 

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
	PRINTF("negative tolerance for distance in nugget covariance to allowed;set to zero");
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
  return nugget(&dist, p, 1)==1.0;
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
	  first = kc;
	  nonzero_pos=ANISO;
	  while ((nonzero_pos<endfor) && (first->param[nonzero_pos]==0.0))
	    nonzero_pos++;
	  if (nonzero_pos>=endfor) {
	      Xerror=ERRORTRIVIAL; goto ErrorHandling; 
	  }
	}
      }
      actcov++;
      kc->left = false;
    }
  } // v
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
  kc = &(key->cov[meth->covlist[0]]);
  s->simple = key->timespacedim == kc->truetimespacedim;
  s->simugrid = kc->simugrid;
  if (key->anisotropy && !s->simple) {
    int *pos, oldpos;
    if (NUGGET_TOL==0.0 && GENERAL_PRINTLEVEL>=1)
      PRINTF("\nThe anisotropy matrix does not have full rank and RFparameters()$nugget.tol equals 0. From a theoretical point of view that's fine, but the simulations will probably be odd. Is this really what you want?\n");
    if ((Xerror=Transform2NoGrid(key, meth->covlist[0])) != NOERROR) 
      goto ErrorHandling;
    if (kc->simugrid) {
      int d, dim=key->timespacedim + 1;
      for (d=i=0; i<key->timespacedim; i++, d+=dim)
	s->diag[i] = kc->param[ANISO + d];
    } else {
      if ((pos = (int*) malloc(sizeof(int) * key->totalpoints))==0) {
	Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      ordering(kc->x, key->totalpoints, kc->truetimespacedim, pos);
      oldpos = pos[0];
      for (i=1 /* ! */; i<key->totalpoints; i++) {
	if (equal(oldpos, pos[i], kc->x, kc->truetimespacedim)) 
	  pos[i]= -1 - pos[i];
	else oldpos=pos[i];
      }
    s->pos = pos;
    }
  }
  nugget_effect = CovFct(ZERO, 1, key->cov, key->meth[m].covlist, actcov, false);
  s->sqrtnugget = sqrt(nugget_effect);
  
  if (key->anisotropy) return NOERROR_REPEAT;
  else return 0;

 ErrorHandling:
  return Xerror;
}

void nuggetgrid(double *res, bool *zerodiag, int *length, int k,
		double sqrtnugget, int *p, int *p_ref) {
//  printf("%d %d %d\n", k, *p, *p_ref);
  if (k<0) {
    if (*p_ref < 0) {
      res[(*p)++] = GAUSS_RANDOM(sqrtnugget);
//      printf("res[*p=%d]=%f\n", (*p)-1, res[(*p)-1]);
    } else {
      res[(*p)++] = res[(*p_ref)++];
    }
  } else {
    int i, p_save;
    if (*p_ref<0) {
      p_save = *p;
      nuggetgrid(res, zerodiag, length, k-1, sqrtnugget, p, p_ref);
      for (i=1; i<length[k]; i++) {
	*p_ref = (zerodiag[k]) ? p_save : -1;
//	printf("k=%d ii=%d %d %d\n", k, i, zerodiag[k], *p_ref);
	nuggetgrid(res, zerodiag, length, k-1, sqrtnugget, p, p_ref); 
      }
    } else {
      for (i=0; i<length[k]; i++) {
	nuggetgrid(res, zerodiag, length, k-1, sqrtnugget, p, p_ref);
      }
    }
  }
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
      int p_ref, p ,d;
      bool zerodiag[MAXDIM];
      for (d=0; d<key->timespacedim; d++) 
	zerodiag[d] = (s->diag[d] == 0.0);
      p = 0;
      p_ref = -1;
      nuggetgrid(res, zerodiag, key->length, key->timespacedim-1, sqrtnugget, 
		 &p, &p_ref);
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

