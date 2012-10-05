
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 all around the nugget effect -- needs special treatment 

 Copyright (C) 2001 -- 2011 Martin Schlather, 

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
 
#include <R_ext/Lapack.h>
#include "RF.h"
#include "primitive.h"

void SetParamNugget(int *action, double *nuggettol, int *nuggetmeth)
{
  nugget_param *gp = &(GLOBAL.nugget);
  if (*action) {
    gp->tol = *nuggettol;
    if (gp->tol < 0) {
      if (PL>=1) 
	PRINTF("negative tolerance for distance in nugget covariance not allowed; set to zero");
      gp->tol = 0.0;
    }
    gp->meth = (bool) *nuggetmeth;
  } else {
    *nuggettol = gp->tol;
    *nuggetmeth = (int) gp->meth;
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
  double *x, *y, v, dummy, dist;
  cov_model cov;
  int d;
  x = X + i * dim;
  y = X + j * dim;
  for (dist=0.0, d=0; d<dim; d++) {
    dummy = x[d]-y[d];
    dist += dummy * dummy;
  }
  dist = sqrt(dist);
  cov.vdim = 1;
  nugget(&dist, &cov, &v);
  return v==1.0;
}

// uses global RANDOM !!!
int init_nugget(method_type *meth){
  cov_model *cov = meth->cov;
  location_type *loc = meth->loc;
  nugget_storage *s;
  globalparam *gp = meth->gp;
  nugget_param* lp = &(gp->nugget);
  int err, i, vdim,
    PL = gp->general.printlevel,
    origdim = loc->timespacedim,
      dim = cov->tsdim,
      dimSq = origdim * origdim;

  SET_DESTRUCT(nugget_destruct);
    if ((meth->S=malloc(sizeof(mpp_storage)))==0){
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (nugget_storage*) meth->S;
  s->pos = NULL;
  s->red_field = NULL;
 
  vdim = cov->calling->vdim;
  if (vdim <= 0) {
    err=ERRORMULTIMISMATCH; goto ErrorHandling;
  }
  cov->vdim = vdim;
  

//  PrintMethodInfo(meth);
//  PrintMethodInfo(meth);

  if ((s->simple = origdim == dim)) {
    double value[MAXNUGGDIM], ivalue[MAXNUGGDIM], dummy[5 * MAXNUGGDIM], *A;
    int j,k, m,
	ndummy = 5 * origdim, err;
    char No = 'N';
//	Upper = 'U';
   
    if (dim > MAXNUGGDIM) {
	ERR("dim larger then MAXNUGGDIM");
    } 
    
    if (meth->cproj != NULL ){
      ERR("projections are not programmed for nugget yet");
    } 
    if (meth->caniso == NULL) {
	s->simple = true;
    } else {    
      A = (double*) malloc(dimSq * sizeof(double));
      // memcpy(A, meth->caniso, sizeof(double) * origdim * origdim);
      // because of numerical errors
      // A =  canio^T * caniso 
      for (k=i=0; i<dimSq; i+=origdim) {
 	for (j=0; j<dimSq; j+=origdim,k++) {
	    A[k] = 0.0;
	    for (m=0; m<origdim; m++) {
		A[k] += meth->caniso[i+m] * meth->caniso[j+m];
	    }
	}
      }
      
      //   printf("%e %e %e %e\n", A[0], A[1], A[2], A[3]);
      //   A[0] = A[1] = A[2] = A[3] = 1.0;
      
      F77_NAME(dgeev)(&No, &No, &origdim, A, &origdim, 
		      value, ivalue, NULL, &origdim, NULL, &origdim,
		      dummy, &ndummy, &err);
      if (err != 0) {
	free(A);
	ERR("dgeev failed in nugget.cc");
      }
      for (i=0; i<origdim; i++) {
//	printf("simple %d %e %e %e  %d\n",
//	       i, value[i], ivalue[i],EIGENVALUE_EPS,
//	       fabs(value[i]) + fabs(ivalue[i]) > EIGENVALUE_EPS);
	if (!(s->simple = fabs(value[i]) + fabs(ivalue[i]) > EIGENVALUE_EPS))
	   break;
      }
      free(A);
    }
  }

  s->simugrid = loc->grid && (meth->type == TypeDiag || meth->type == TypeIso);

//  printf("%d\n", s->simple); assert(false);

  if (!s->simple) {
    int *pos, oldpos;
    if (lp->tol==0.0 && PL>=1)
      PRINTF("\nThe anisotropy matrix does not have full rank and RFparameters()$nugget.tol equals 0. From a theoretical point of view that's fine, but the simulations will probably be odd. Is this really what you want?\n");
    if (s->simugrid) {
      // TypeIso und nicht simple genau dann wenn identisch 0
      int d, dimP1=dim + 1;
      s->prod_dim[0] = 1; 
      if (meth->cproj != NULL ){
	  assert(false); 
	  // only projection possible here, but not programmed yet
      } else {
        for (d=i=0; d<dim; d++, i+=dimP1) {
          s->reduced_dim[d] =
	      fabs(meth->caniso[i]) < lp->tol ? 1 : loc->length[d];
	  s->prod_dim[d + 1] = s->prod_dim[d] * s->reduced_dim[d];
 
//     printf("%d %d %d %d\n", i, s->prod_dim[0], s->prod_dim[1],
//	     s->prod_dim[2]);
	}
      }

      if ((s->red_field=(res_type *) malloc(sizeof(res_type) * cov->vdim *
					  s->prod_dim[dim])) 
	  == NULL){
	  err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
    } else {
      if ((pos = (int*) malloc(sizeof(int) * loc->totalpoints))==0) {
	err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      Transform2NoGridExt(meth, false, true, &(meth->space), &(meth->sptime));
//	assert(false);
      ordering(meth->sptime, loc->totalpoints, dim, pos);
      oldpos = pos[0];
      for (i=1 /* ! */; i<loc->totalpoints; i++) {
	if (equal(oldpos, pos[i], meth->sptime, cov->tsdim)) 
	  pos[i]= -1 - pos[i];
	else oldpos=pos[i];
      }
    s->pos = pos;
    }
  }
  // s->sqrtnugget = sqrt(meth->cvar);
  // printf("%f %f\n", meth->cvar, s->sqrtnugget); assert(false);
  
  return  NOERROR;

 ErrorHandling:
  return err;
}

void do_nugget(method_type *meth, res_type *res) {
  cov_model *cov = meth->cov;
  location_type *loc = meth->loc;
///  double sqrtnugget;
  nugget_storage* s;
  long nx, endfor;
  int v, 
    vdim = cov->vdim;

  assert(meth->S != NULL);
  s = (nugget_storage*) meth->S;
//  sqrtnugget = s->sqrtnugget;

  if (s->simple) {
      for (nx=0, endfor=loc->totalpoints * vdim; nx<endfor; nx++) {
//	  if (nx==0) printf("nug %f ", res[nx]);
	  res[nx] += (res_type) GAUSS_RANDOM(1.0);
//	  if (nx==0)printf("%f \n",  (double) res[nx]);
      }
  } else {
    if (s->simugrid) {
	int d, i, dim, dimM1, index[MAXNUGGDIM], *red_dim, *prod_dim, endfor;
      long totpnts, idx;
      res_type *field;
      totpnts = loc->totalpoints;
      dim = cov->tsdim;
      dimM1 = dim -1;
      field = s->red_field;
      red_dim = s->reduced_dim;
      prod_dim = s->prod_dim;

      //PrintMethodInfo(meth);
      
//      printf("%d %d %d %d\n", dim, s->prod_dim[0], s->prod_dim[1],
//	     s->prod_dim[2]);
      endfor = vdim * s->prod_dim[dim-1];
      for (i=0; i<endfor; i++) {
	  field[i] = (res_type) GAUSS_RANDOM(1.0);
      }
      for (d=0; d<dim; index[d++] = 0);
      for (i=0; i<totpnts; i++) {
	for(idx=d=0; d<dim; d++) idx += (index[d] % red_dim[d]) * prod_dim[d];
	//printf("%d %d %f\n", i, idx, (double) field[idx]);
	for (v=0; v<vdim; v++) {
	  res[i + v] += field[idx + v];
	}
	d = 0; 
	(index[d])++; 
	while (d < dimM1 && index[d] >= loc->length[d]) { 
	  assert(d<dim); // assert(d>=0);
	  index[d] = 0; 
	  d++;   
	  (index[d])++; 
	}
      }
    } else {
      int p;
      double *dummy = (double*) malloc(sizeof(double) * vdim);
      assert(s->pos[0]>=0);
      for (v=0; v<vdim; v++) dummy[v] = RF_NAN; // just to avoid warnings 
      //                                           from the compiler
      for (nx=0; nx<loc->totalpoints; nx++) {
	if ((p = s->pos[nx]) < 0) p = -1 - p; // if p<0 then take old variable
	// and -p-1 is the true index 
	else {
	  for (v=0; v<vdim; v++)
	    dummy[v] = GAUSS_RANDOM(1.0);
	}
	for (v=0; v<vdim; v++)
	  res[p + v] += dummy[v];
      }
      free(dummy);
    }
  }
}

