
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 all around the nugget effect -- needs special treatment 

 Copyright (C) 2001 -- 2014 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
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
#include "randomshape.h"

#define NUGGET_TOL 0
#define NUGGET_VDIM 1
#define NUGGET_PROC_TOL (COMMON_GAUSS + 1)
#define NUGGET_PROC_VDIM (COMMON_GAUSS + 2)

bool equal(cov_model *cov, int i, int j, double *X, int dim)
{
  double *x, *y, v, dummy, dist;
  int d;
  x = X + i * dim;
  y = X + j * dim;
  for (dist=0.0, d=0; d<dim; d++) {
    dummy = x[d]-y[d];
    dist += dummy * dummy;
  }
  dist = sqrt(dist);
  nugget(&dist, cov, &v);
  return v==1.0;
}


/* nugget effect model */
void nugget(double *x, cov_model *cov, double *v) {

  double diag = (*x <= P0(NUGGET_TOL)) ? 1.0 : 0.0;
  int i, endfor,
    vdim   = cov->vdim,
    vdimsq = vdim * vdim;
  
  v[0] = diag;
  for (i = 1; i<vdimsq; v[i++] = diag) {
    endfor = i + vdim;
    for (; i<endfor; v[i++] = 0.0);
  }
}


void covmatrix_nugget(cov_model *cov, double *v) {
  location_type *loc = Loc(cov);
  int i,
    vdim   = cov->vdim,
    n = loc->totalpoints * vdim,
    nP1 = n + 1,
    n2 = n * n;
  for (i=0; i<n2; v[i++]=0.0);
  for (i=0; i<n2; i += nP1) v[i]=1.0;

  //  {  int i,j,k, tot=Loc(cov)->totalpoints; printf("\nnugget %d %d %d\n", n, n2, nP1);
  //   for (k=i=0; i<tot*tot; i+=tot) {
  //     for (j=0; j<tot; j++) printf("%f ", v[k++]);
  //    printf("\n");  }}

  ///APMI(cov);

}

char iscovmatrix_nugget(cov_model VARIABLE_IS_NOT_USED *cov) {  return true; }

void Inversenugget(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) { 
  *v = (*x==1.0) ? 0.0 : RF_INF; //or better 0.0 => error?
}


int check_nugget(cov_model *cov) {
#define nsel 4
  int  err ; // taken[MAX DIM],
  nugget_param *gp  = &(GLOBAL.nugget);

  
  assert(cov->nr == NUGGET);
  ROLE_ASSERT(ROLE_COV || cov->role == ROLE_GAUSS);

  kdefault(cov, NUGGET_TOL, gp->tol);
  if (PisNULL(NUGGET_VDIM)) {

    //printf("nugget vdim %d\n", cov->vdim);

    if (cov->vdim <= 0) cov->vdim = 1;
    kdefault(cov, NUGGET_VDIM, cov->vdim);
  } else {
    cov->vdim = P0INT(NUGGET_VDIM);
  }

  cov->matrix_indep_of_x = true;
  
  if ((err = checkkappas(cov)) != NOERROR) return err;
  
  return NOERROR;
}


void range_nugget(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[NUGGET_TOL] = 0;
  range->max[NUGGET_TOL] = RF_INF;
  range->pmin[NUGGET_TOL] = 0;
  range->pmax[NUGGET_TOL] = 1e-5;
  range->openmin[NUGGET_TOL] = false;
  range->openmax[NUGGET_TOL] = true; 

  range->min[NUGGET_VDIM]  = 1.0;
  range->max[NUGGET_VDIM]  = RF_INF;
  range->pmin[NUGGET_VDIM] = 1.0;
  range->pmax[NUGGET_VDIM] = 10.0;
  range->openmin[NUGGET_VDIM] = false;
  range->openmax[NUGGET_VDIM] = true;
}


int check_nugget_proc(cov_model *cov) {
#define nsel 4
  cov_model *next = cov->sub[0],
    *key = cov->key,
    *sub = (key==NULL) ? next : key,
    *intern;
  int  err, // taken[MAX DIM],
    dim = cov->tsdim;
  //  nugget_param *gp  = &(GLOBAL.nugget);
  //PMI(cov);

  // print("NUGGET PROC !!!! \n");

  ROLE_ASSERT(ROLE_GAUSS);
 
  if (key == NULL) { 
    intern = sub;
    while (intern != NULL && isDollar(intern)) {
      intern = intern->key != NULL ? intern->key : intern->sub[0];
    }
    if (intern->nr != NUGGET) {
      //printf("cov=%s\n", CovList[intern->nr].name);
      //APMI(cov);
      SERR2("'%s' only allows for '%s'", NICK(cov),
	    CovList[NUGGET].nick);
    }
    if (!PisNULL(NUGGET_PROC_TOL)) 
      kdefault(intern, NUGGET_TOL, P0(NUGGET_PROC_TOL));
    if (!PisNULL(NUGGET_PROC_VDIM)) 
      kdefault(intern, NUGGET_VDIM, P0INT(NUGGET_PROC_VDIM));
    if ((err = CHECK(next, dim, dim, PosDefType, KERNEL, SYMMETRIC,
		       SUBMODEL_DEP, ROLE_COV)) != NOERROR) return err;
    if (!PARAMisNULL(intern, NUGGET_TOL))
      kdefault(cov, NUGGET_PROC_TOL, PARAM0(intern, NUGGET_TOL));  
    if (!PARAMisNULL(intern,NUGGET_VDIM))
      kdefault(cov, NUGGET_PROC_VDIM, PARAM0INT(intern, NUGGET_VDIM));
  } else {    // key != NULL  && next != nugget   
    // dann ruft NuggetIntern Nugget auf 
    intern = cov->nr == NUGGET_USER ? sub : cov;
    while (intern != NULL && isAnyDollar(intern)) {
      intern = intern->key != NULL ? intern->key : intern->sub[0];
    }
    if (intern == NULL || intern->nr != NUGGET_INTERN) {
      //PMI(cov);  APMI(intern);
      BUG;
    } else if (intern != cov) {
      //print("****** here\n");
      paramcpy(intern, cov, true, false);
    }

    if (!PisNULL(NUGGET_PROC_TOL)) 
      kdefault(intern, NUGGET_PROC_TOL, P0(NUGGET_PROC_TOL));
    if (!PisNULL(NUGGET_PROC_VDIM)) 
      kdefault(intern, NUGGET_PROC_VDIM, P0INT(NUGGET_PROC_VDIM));
    
    if ((err = CHECK(key, dim, dim, ProcessType, XONLY, CARTESIAN_COORD,
		       SUBMODEL_DEP, ROLE_GAUSS)) != NOERROR) {
      // printf("error nug prox %d\n", err);
      return err;
    }
  } 
   


  cov->vdim = next->vdim;  
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown)
    return ERRORDIM;
  if ((err = check_common_gauss(cov)) != NOERROR) return err;
  cov->vdim = sub->vdim;  
  cov->role = ROLE_GAUSS;
  
  // printf("OK nugget\n");
  EXTRA_STORAGE;

  return NOERROR;
}


void range_nugget_proc(cov_model *cov, range_type *range) {
  range_common_gauss(cov, range);

  range->min[NUGGET_PROC_TOL] = 0;
  range->max[NUGGET_PROC_TOL] = RF_INF;
  range->pmin[NUGGET_PROC_TOL] = 0;
  range->pmax[NUGGET_PROC_TOL] = 1e-5;
  range->openmin[NUGGET_PROC_TOL] = false;
  range->openmax[NUGGET_PROC_TOL] = true; 

  range->min[NUGGET_PROC_VDIM]  = 1.0;
  range->max[NUGGET_PROC_VDIM]  = RF_INF;
  range->pmin[NUGGET_PROC_VDIM] = 1.0;
  range->pmax[NUGGET_PROC_VDIM] = 10.0;
  range->openmin[NUGGET_PROC_VDIM] = false;
  range->openmax[NUGGET_PROC_VDIM] = true;
}

// uses global RANDOM !!!
int init_nugget(cov_model *cov, storage VARIABLE_IS_NOT_USED *S){

  // printf(" CALL OF INIT NUGGET ****************\n");

  location_type 
    *loc=cov->prevloc; 
  if (cov->ownloc!=NULL) {
    LOC_DELETE(&(cov->ownloc));
    //PMI(cov);
    // ERR("unexpected call of nugget");
  }
  cov_model *next = cov->sub[0];
  nugget_storage *s;
  int d, //vdim,
    err = NOERROR,
    origdim = loc->timespacedim,
    dim = cov->tsdim,
    dimSq = origdim * origdim;
  double
    tol = P0(NUGGET_PROC_TOL);
  matrix_type anisotype = TypeMdiag;
  
  ROLE_ASSERT_GAUSS;

  cov->method = Nugget;
  if (cov->Snugget != NULL) NUGGET_DELETE(&(cov->Snugget));
  if ((cov->Snugget = (nugget_storage*) MALLOC(sizeof(nugget_storage)))==NULL){
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (nugget_storage*) cov->Snugget;
  NUGGET_NULL(s);
  s->pos = NULL;
  s->red_field = NULL;

  if (next->nr != NUGGET) {
     err=ERRORCOVFAILED;
     strcpy(ERRORSTRING_OK, CovList[NUGGET].nick);
     strcpy(ERRORSTRING_WRONG, NICK(cov));
     goto ErrorHandling;
  }
  

  if ((s->simple = origdim == dim)) {
    double value[MAXNUGGDIM], ivalue[MAXNUGGDIM], dummy[5 * MAXNUGGDIM], *A;
    int 
      ndummy = 5 * origdim;
    char No = 'N';
//	Upper = 'U';
    //   PMI(cov);
    
    if (loc->caniso == NULL) {
      if (loc->grid) {
	for (d=0; d<dim; d++) 
	  if (fabs(loc->xgr[d][XSTEP]) <=  tol) {
	    s->simple = false;
	    break;
	  }
      }
    } else {    
     if (dim > MAXNUGGDIM) {
	GERR("dim larger then MAXNUGGDIM");
      } 
      anisotype = Type(loc->caniso, loc->cani_nrow, loc->cani_ncol);
      A = (double*) MALLOC(dimSq * sizeof(double));
      AtA(loc->caniso, origdim, origdim, A);
       
      //   print("A %e %e %e %e\n", A[0], A[1], A[2], A[3]);
      //   A[0] = A[1] = A[2] = A[3] = 1.0;
      
      F77_NAME(dgeev)(&No, &No, &origdim, A, &origdim, 
		      value, ivalue, NULL, &origdim, NULL, &origdim,
		      dummy, &ndummy, &err);
      if (err != 0) {
	free(A);
	GERR("dgeev failed in nugget.cc");
      }
      for (d=0; d<origdim; d++) {
	//	print("simple %d %e %e %e  %d\n",
	//	       d, value[d], ivalue[d],EIGENVALUE_EPS,
	//       fabs(value[d]) + fabs(ivalue[d]) > EIGENVALUE_EPS);
	if (!(s->simple = fabs(value[d]) + fabs(ivalue[d]) > EIGENVALUE_EPS))
	   break;
      }
      free(A);
    }
  }

  s->simugrid = loc->grid && isMdiag(anisotype);
  assert(!s->simugrid || loc->caniso == NULL);

//   print("simple %d %d\n", s->simple, s->simugrid); //assert(false);

  if (!s->simple) {
    int *pos, oldpos;
    if (tol==0.0 && PL>=PL_IMPORTANT) {
      PRINTF("\nThe anisotropy matrix does not have full rank and the parameter 'tol' equals 0. From a theoretical point of view that's fine, but the simulations will probably be odd. Is this really what you want?\n");
    }
    if (s->simugrid) {
      // TypeIso und nicht simple genau dann wenn identisch 0
      s->prod_dim[0] = 1; 

      for (d=0; d<dim; d++) {
	s->reduced_dim[d] = fabs(loc->xgr[d][XSTEP]) <= tol ? 1 : loc->length[d];
	s->prod_dim[d + 1] = s->prod_dim[d] * s->reduced_dim[d];
	
      //print("d=%d %d %d %d\n",d,s->prod_dim[0],s->prod_dim[1],s->prod_dim[2]);
      }

      if ((s->red_field=(res_type *) MALLOC(sizeof(res_type) * cov->vdim *
					    s->prod_dim[dim]))
	  == NULL){
	  err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
    } else {
      int i;
      if ((pos = (int*) MALLOC(sizeof(int) * loc->totalpoints)) == NULL) {
	err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      Transform2NoGrid(cov, false, true);
      loc = cov->ownloc;
      ordering(loc->x, loc->totalpoints, dim, pos);

      //      for (i=0; i<loc->totalpoints; i++) {
      //	printf("old %f %f pos[i]=%d new=%f %f\n",
      //	       loc->x[0 + 2 * i], loc->x[1 + 2 * i], pos[i],
      //	       loc->x[0 + 2 * pos[i]], loc->x[1 + 2 * pos[i]]);
      //      }

      oldpos = pos[0];
      for (i=1 /* ! */; i<loc->totalpoints; i++) {
	if (equal(next, oldpos, pos[i], loc->x, cov->tsdim)) {
	  pos[i]= -1 - pos[i];
	  // printf("oldpos %d\n", i);
	} else oldpos=pos[i];
      }
    s->pos = pos;
    }
  }
  
  if ((err = FieldReturn(cov)) != NOERROR) goto ErrorHandling;

 ErrorHandling:
  cov->simu.active = err == NOERROR;

  return err;
}

void do_nugget(cov_model *cov, storage VARIABLE_IS_NOT_USED *S) {
  location_type *loc = Loc(cov);
///  double sqrtnugget;
  nugget_storage* s;
  long nx, endfor;
  int v, 
    vdim = cov->vdim;
  double 
    *res = cov->rf;
  bool loggauss = (bool) P0INT(LOG_GAUSS);

  s = (nugget_storage*) cov->Snugget;
//  sqrtnugget = s->sqrtnugget;

  if (s->simple) {
    for (nx=0, endfor=loc->totalpoints * vdim; nx<endfor; nx++) {
      res[nx] = (res_type) GAUSS_RANDOM(1.0);
    }
  } else {
    if (s->simugrid) {
      int d, i, dim, dimM1, index[MAXNUGGDIM], *red_dim, *prod_dim;
      long totpnts, idx;
      res_type *field;
      totpnts = loc->totalpoints;
      dim = cov->tsdim;
      dimM1 = dim -1;
      field = s->red_field;
      red_dim = s->reduced_dim;
      prod_dim = s->prod_dim;
      endfor = vdim * s->prod_dim[dim]; // not dim-1

      //PrintMethodInfo(meth);
      
      //print("reddim %d %d dim=%d %d %d %d endfor=%d\n", red_dim[0], red_dim[1], dim, s->prod_dim[0], s->prod_dim[1], s->prod_dim[2], endfor);

      for (i=0; i<endfor; i++) {
	field[i] = (res_type) GAUSS_RANDOM(1.0);
      }
      for (d=0; d<dim; index[d++] = 0);
      for (i=0; i<totpnts; i++) {
	for(idx=d=0; d<dim; d++) idx += (index[d] % red_dim[d]) * prod_dim[d];
	//	print("for %d %d %f d=%d %d len=%d %d \n", i, idx, (double) field[idx], index[0], index[1], loc->length[0], loc->length[1]);
	for (v=0; v<vdim; v++) {
	  res[i + v] = field[idx + v];
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
      ALLOC_EXTRA1(dummy, vdim);
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
	  res[p + v] = dummy[v];
      }
    }
  }

  if (loggauss) {
    int i, vdimtot = loc->totalpoints * cov->vdim;
    for (i=0; i<vdimtot; i++) res[i] = exp(res[i]);
  }
}


int struct_nugget(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  //  PMI(cov, "structhyper");
  if (cov->sub[0]->pref[Nugget] == PREF_NONE) {
    return ERRORPREFNONE;
  }
  if (cov->role != ROLE_GAUSS) {
    SERR("type is not Gaussian.");
  }
  return NOERROR;
}
