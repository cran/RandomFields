
/* 
 Authors
 Martin Schlather, martin.schlather@cu.lu

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


static Real NUGGET_TOL_SINGULAR=1e-12;
static Real NUGGET_TOL_POINTS=1e-12;

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
Real nugget(Real *x, Real *p, int dim){
  if (*x==0.0) return 1.0;  return 0.0;
}
Real Scalenugget(Real *p, int scaling) { return 1.0; }//or better 0.0 => error?

int is_singular(Real *C, int dim, bool *singular)
{
  double D[MAXDIM];
  int Xerror, d;
  if (Xerror=eigenvalues(C, dim, D)) goto ErrorHandling;
  *singular = false;
  for (d=0; d<dim; d++) 
    if (fabs(D[d])<NUGGET_TOL_SINGULAR) {*singular=true; break;}
  return 0;

 ErrorHandling:
  return Xerror;
}

bool equal(int i, int j, Real *ORDERD, int ORDERDIM)
{
  Real *x, *y;
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
  int v, Xerror; 
  nugget_storage *s;
  param_type param;
  Real nugget_effect, quotient[MAXCOV];
  char actcov;
  int i,covnr[MAXCOV], nonzero_pos;
  int multiply[MAXCOV], TrueDim;
  bool time_exception[MAXCOV], singular, no_last_comp;
  Real *xx;

  xx=NULL;
  SET_DESTRUCT(nugget_destruct);
  FIRSTCHECK_COV_ANISO(Nugget, cov, param, false);

  if ((key->S[m]=malloc(sizeof(mpp_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (nugget_storage*) key->S[m];
  s->pos = NULL;

  if (key->anisotropy) {
    if (Xerror = is_singular(&(param[0][ANISO]), key->timespacedim, &singular)) 
      goto ErrorHandling;
    s->simple = !singular;
    if (singular) {
      int *pos, oldpos, start_param[MAXDIM], index_dim[MAXDIM];
      GetTrueDim(key->anisotropy, key->timespacedim, param[0],
		 &TrueDim, &no_last_comp, start_param, index_dim);
      if (Xerror=Transform2NoGrid(key, param[0], TrueDim, start_param, &(xx)))
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

void do_nugget(key_type *key, bool add, int m, Real *res ) {
  Real sqrtnugget;
  nugget_storage* s;
  long nx, endfor;
#define RESULT(OP)\
  for (nx=0,endfor=key->totalpoints; nx<endfor; nx++)\
     res[nx] OP GAUSS_RANDOM(sqrtnugget);\

  s = (nugget_storage*) key->S[m];
  sqrtnugget = s->sqrtnugget;

  if (s->simple) {
    if (add) { RESULT(+=) } else { RESULT(=) } 
  } else {
    int p;
    Real dummy;
    assert(s->pos[0]>=0);
    for (nx=0; nx<key->totalpoints; nx++) {
      if ((p=s->pos[nx])<0) p= -1 - p; // if p<0 then take old variable
      // and -p-1 is the true index 
      else dummy = GAUSS_RANDOM(sqrtnugget);
      if (add) res[p] += dummy;
      else res[p] = dummy;  
    }
  }
}

void rangenugget(int dim, int *index, Real* range){
  *index = -1;
}
