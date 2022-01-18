
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 all around the nugget effect -- needs special treatment 

 Copyright (C) 2001 -- 2017 Martin Schlather, 

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

#include <stdio.h>  
#include "def.h"
#include <Basic_utils.h>
#include <R_ext/Lapack.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "extern.h"
#include "questions.h"
#include "primitive.h"
#include "Processes.h"
#include "operator.h"


bool equal(model *cov, int i, int j, double *X, int dim) {
  double *x, *y, v, dummy, dist;
  int d;
  x = X + i * dim;
  y = X + j * dim;
  for (dist=0.0, d=0; d<dim; d++) {
    dummy = x[d]-y[d];
    dist += dummy * dummy;
  }
  dist = SQRT(dist);
  nugget(&dist, cov, &v);
  return v==1.0;
}

bool SpatialNugget(model *cov) {
  /* NOTE THAT ONLY THE PRECEDING DOLLAR IS CHECKED, not 
     other ones -- so the latter do not influence the behaviour 
     of the nugget effect 
  */
  
  if (!GLOBAL.internal.allow_duplicated_loc) return true; // all interpreted as
  //                                 spatial nugget, not as measurement error
  assert(cov != NULL);
  // PMI0(cov);
  // printf("isany = %d %d\n", isAnyDollar(cov), PisNULL(DANISO));
  assert(isAnyNugget(cov));
  model *calling = cov->calling;
  if (calling == NULL) return false;
  if (equalsNuggetProc(calling)) calling = calling->calling;
  //printf("calling %.50s %d,%d %d %d\n", NAME(calling), CALLINGNR, MODELNR(calling), GAUSSPROC, CALLINGNR== GAUSSPROC);
 if (calling == NULL) return false;
  if (MODELNR(calling) == GAUSSPROC) calling = calling->calling;
  if (calling == NULL || !isAnyDollar(calling)) return false;

  model *Scale = calling->kappasub[DSCALE],
    *Aniso = calling->kappasub[DANISO],
    *User = calling->kappasub[DAUSER];
  if (!PARAMisNULL(calling, DSCALE) || Scale != NULL)
    ERR0("'scale' does not make sense within a nugget effect. However 'Aniso' does.\nSee the manual.");

  return !(PARAMisNULL(calling, DANISO) && Aniso==NULL &&
	   PARAMisNULL(calling, DAUSER) && User == NULL &&
	   PARAMisNULL(calling, DPROJ));
}

/* nugget effect model */
void nugget(double *x, model *cov, double *v) {
double same = (*x <= P0(NUGGET_TOL)) ? 1.0 : 0.0;
  int i, endfor,
    vdim   = VDIM0,
    vdimsq = vdim * vdim;
  assert(VDIM0 == VDIM1);
  assert(cov->Snugget->spatialnugget);
  //  printf("spatial nugget effect\n");

  v[0] = same;
  for (i = 1; i<vdimsq; v[i++] = same) {
    endfor = i + vdim;
    for (; i<endfor; v[i++] = 0.0);
  }
  //  printf("%10g %10g %d\n", *x, *v, PisNULL(NUGGET_TOL));
}

void nuggetnonstat(double *x, double VARIABLE_IS_NOT_USED *y, model *cov,
		   double *v) {
  // printf("nonstat nugget\n"); 
  int i, endfor,
    dim = ANYDIM, //  OWNXDIM(0)
    vdim   = VDIM0,
    vdimsq = vdim * vdim;
  double same = (*x == 0.0 && y == NULL) || x[dim] == y[dim] ? 1.0 : 0.0;
  assert(VDIM0 == VDIM1);
  assert(!cov->Snugget->spatialnugget);

  v[0] = same;
  for (i = 1; i<vdimsq; v[i++] = same) {
    endfor = i + vdim;
    for (; i<endfor; v[i++] = 0.0);
  }
  //  printf("%10g %10g %d\n", *x, *v, PisNULL(NUGGET_TOL));
}


void covmatrix_nugget(model *cov, double *v) {
  location_type *loc = Loc(cov);
  int 
    vdim   = VDIM0;
  long i,
    n = loc->totalpoints * vdim,
    nP1 = n + 1,
    n2 = n * n;
  
  if (cov->Snugget->spatialnugget) BUG;
  for (i=0; i<n2; v[i++] = 0.0);
  for (i=0; i<n2; i += nP1) v[i]=1.0;
}

char iscovmatrix_nugget(model VARIABLE_IS_NOT_USED *cov) {
  return !cov->Snugget->spatialnugget;
}

void Inversenugget(double VARIABLE_IS_NOT_USED *x, model VARIABLE_IS_NOT_USED *cov, double *v) { 
  *v = 0.0; ///(*x==1.0) ? 0.0 : RF_INF; //or better 0.0 => error?
}


int check_nugget(model *cov) {
#define nsel 4
  int  err ; // taken[MAX DIM],
  nugget_param *gp  = &(GLOBAL.nugget);
 
 
  assert(equalsNugget(COVNR));
  if (!hasAnyEvaluationFrame(cov) && !hasAnyProcessFrame(cov)) ILLEGAL_FRAME;

  kdefault(cov, NUGGET_TOL, gp->tol);
  if (PisNULL(NUGGET_VDIM)) {
    if (VDIM0 <= 0) VDIM0 = VDIM1 = 1;
    kdefault(cov, NUGGET_VDIM, VDIM0);
  } else {
    VDIM0 = VDIM1 = P0INT(NUGGET_VDIM);
  }
  cov->matrix_indep_of_x = true;
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
 
  if (cov->Snugget == NULL) {
    NEW_STORAGE(nugget);
    cov->Snugget->spatialnugget = SpatialNugget(cov);
  }

  if (!GLOBAL.internal.allow_duplicated_loc) {
    for (int i=CircEmbed; i<Nothing; i++)
      cov->pref[i] = (cov->pref[i] > 0) * PREF_BEST;
  } else if (cov->Snugget->spatialnugget) {
    assert(CircEmbed == 0);
    for (int i=CircEmbed; i<Nothing; i++) cov->pref[i] = PREF_NONE;
    cov->pref[Nugget]= cov->pref[Nothing] = PREF_BEST;
  } // else standard, see init.cov.cc

  //  printf("%d\n", cov->Snugget->spatialnugget);  APMI(cov);

  RETURN_NOERROR;
}


void range_nugget(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[NUGGET_TOL] = 0;
  range->max[NUGGET_TOL] = RF_INF;
  range->pmin[NUGGET_TOL] = 0;
  range->pmax[NUGGET_TOL] = 1e-5;
  range->openmin[NUGGET_TOL] = false;
  range->openmax[NUGGET_TOL] = true; 

  range->min[NUGGET_VDIM]  = 1.0;
  range->max[NUGGET_VDIM]  = MAXINT;
  range->pmin[NUGGET_VDIM] = 1.0;
  range->pmax[NUGGET_VDIM] = 10.0;
  range->openmin[NUGGET_VDIM] = false;
  range->openmax[NUGGET_VDIM] = true;
}

Types Typenugget(Types required, model *cov, isotropy_type requ_iso){
  if (cov->Snugget == NULL) {
    NEW_STORAGE(nugget);
    cov->Snugget->spatialnugget = SpatialNugget(cov);
  }
  if (cov->Snugget->spatialnugget || equalsCoordinateSystem(requ_iso) ||
      ( (PisNULL(NUGGET_VDIM) || P0INT(NUGGET_VDIM) == 1) &&
	isSymmetric(requ_iso)))
    return TypeConsistency(required, TcfType); 
    
  return BadType;
}


bool setnugget(model *cov) {
  isotropy_type iso = CONDPREVISO(0);
  if (!isFixed(iso)) return false;

  if (cov->Snugget == NULL) {
    NEW_STORAGE(nugget);
    cov->Snugget->spatialnugget = SpatialNugget(cov);
  }

  if (cov->Snugget->spatialnugget) {
    set_dom(OWN, 0, XONLY);
    set_iso(OWN, 0, IsotropicOf(iso));
  } else {
    set_dom(OWN, 0, KERNEL);
    if (PisNULL(NUGGET_VDIM) || P0INT(NUGGET_VDIM) == 1)
      set_iso(OWN, 0, SymmetricOf(iso));
    else set_iso(OWN, 0, CoordinateSystemOf(iso));
  }
  return true;
}

bool allowedDnugget(model *cov) {
  if (cov->Snugget == NULL) {
    NEW_STORAGE(nugget);
    cov->Snugget->spatialnugget = SpatialNugget(cov);
  }
  bool *D = cov->allowedD;
  for (int i=FIRST_DOMAIN; i<LAST_DOMAINUSER; D[i++]=false);
  D[cov->Snugget->spatialnugget ? XONLY : KERNEL] = true;
  return false;
}
  
bool allowedInugget(model *cov) {
  if (cov->Snugget == NULL) {
    NEW_STORAGE(nugget);
    cov->Snugget->spatialnugget = SpatialNugget(cov);
  }
  bool *I = cov->allowedI;
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  if (cov->Snugget->spatialnugget) {
    I[ISOTROPIC] = I[EARTH_ISOTROPIC] = I[SPHERICAL_ISOTROPIC] = true;
  } else {
     if (PisNULL(NUGGET_VDIM) ||  P0INT(NUGGET_VDIM) == 1)
      I[SYMMETRIC] = I[EARTH_SYMMETRIC] = I[SPHERICAL_SYMMETRIC] = true;
    else
      I[CARTESIAN_COORD] = I[EARTH_COORD] = I[SPHERICAL_COORD] = true;
  }
  return false;
}

int check_nugget_proc(model *cov) {
#define nsel 4
  model *next = cov->sub[0],
    *key = cov->key,
    *sub = (key==NULL) ? next : key,
    *intern;
  int err;

  ASSERT_CARTESIAN;
 
  nugget_storage *s = cov->Snugget;
  if (s == NULL) {
    NEW_STORAGE(nugget);
    s = cov->Snugget;
    s->spatialnugget = SpatialNugget(cov);
  }

  if (key == NULL) { 
    intern = sub;
    while (intern != NULL && isDollar(intern)) {
      intern = intern->key != NULL ? intern->key : intern->sub[0];
    }
    if (!equalsNugget(MODELNR(intern))) {
      SERR2("'%.50s' only allows for '%.50s'", NICK(cov), DefList[NUGGET].nick);
    }
    if (!PisNULL(NUGGET_PROC_TOL)) 
      kdefault(intern, NUGGET_TOL, P0(NUGGET_PROC_TOL));
    if (!PisNULL(NUGGET_PROC_VDIM)) 
      kdefault(intern, NUGGET_VDIM, P0INT(NUGGET_PROC_VDIM));
    assert(equalsCoordinateSystem(OWNISO(0)));
    if ((err = CHECK_THROUGHOUT(next, cov, PosDefType, KERNEL, OWNISO(0),
				SUBMODEL_DEP, EvaluationType))
 	!= NOERROR) RETURN_ERR(err);
   if (!PARAMisNULL(intern, NUGGET_TOL))
      kdefault(cov, NUGGET_PROC_TOL, PARAM0(intern, NUGGET_TOL));  
    if (!PARAMisNULL(intern,NUGGET_VDIM))
      kdefault(cov, NUGGET_PROC_VDIM, PARAM0INT(intern, NUGGET_VDIM));
  } else {    // key != NULL  && next != nugget   
    // dann ruft NuggetIntern Nugget auf 
    intern = COVNR == NUGGET_USER ? sub : cov;
    while (intern != NULL && isAnyDollar(intern)) {
      intern = intern->key != NULL ? intern->key : intern->sub[0];
    }
    if (intern == NULL || MODELNR(intern) != NUGGET_INTERN) {
      //PMI(cov);  APMI(intern);
      BUG;
    } else if (intern != cov) {
      //print("****** here\n");
      paramcpy(intern, cov, true, true, false, false, false);
    }

    if (!PisNULL(NUGGET_PROC_TOL)) 
      kdefault(intern, NUGGET_PROC_TOL, P0(NUGGET_PROC_TOL));
    if (!PisNULL(NUGGET_PROC_VDIM)) 
      kdefault(intern, NUGGET_PROC_VDIM, P0INT(NUGGET_PROC_VDIM));
    
    int dim = ANYDIM;
    if ((err = CHECK(key, dim, dim, ProcessType, XONLY, CARTESIAN_COORD,
		       SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
      // printf("error nug prox %d\n", err);
      RETURN_ERR(err);
    }
  } 

  VDIM0 = next->vdim[0];  
  VDIM1 = next->vdim[1];  
  cov->frame = GaussMethodType;
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}


void range_nugget_proc(model  VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;

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
int init_nugget(model *cov, gen_storage VARIABLE_IS_NOT_USED *S){
#define INTERNAL_ERR -99999
  location_type *loc=PrevLoc(cov);
   int dim = ANYDIM;
 
  if (cov->ownloc!=NULL) {
    LOC_DELETE(&(cov->ownloc));
    //PMI(cov);
    // ERR0("unexpected call of nugget");
  }
  model *next = cov->sub[0];
  nugget_storage *s = cov->Snugget;
  int d, //
    tot = loc->totalpoints,
    *refpos = NULL,
    err = NOERROR,
    //    origdim = loc->timespacedim,
    vdim = VDIM0;
  double
    *A = NULL, 
    *value=NULL, 
    *ivalue=NULL,
    *dummy = NULL,
   tol = P0(NUGGET_PROC_TOL);
  //  matrix_type anisotype = TypeMdiag;
  
  cov->method = Nugget;

  assert(s!=NULL);

 
  if (!equalsNugget(NEXTNR))
    GERR2("'%.50s' was called by '%.50s'", NICK(cov), NICK(next));
  if ((s->datapos = (int*) MALLOC(loc->totalpoints * sizeof(int))) ==NULL) {
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  if (s->spatialnugget) {//spatial nugget; otherwise measurement error
    int *pos;
    //    char No = 'N';
//	Upper = 'U';
    
    if (tol==0.0 && PL>=PL_IMPORTANT) {
      PRINTF("\nThe anisotropy matrix is given and the parameter 'tol' equals 0. From a theoretical point of view that's fine, but the simulations will probably be odd particularly if 'Aniso' does not have full rank. Is this really what you want?\n");
    }
 
    matrix_type anisotype = Type(loc->caniso, loc->cani_nrow, loc->cani_ncol);
    if ( (s->simugrid = loc->grid && isMdiag(anisotype))) { // =  not ==
      // if (loc->caniso != NULL) BUG;
      // TypeIso und nicht simple genau dann wenn identisch 0
      ALLC_NEWINT(Snugget, reduced_dim, dim, reduced_dim);     
      ALLC_NEWINT(Snugget, prod_dim, dim + 1, prod_dim);
      
       prod_dim[0] = 1; 

      for (d=0; d<dim; d++) {
	if (FABS(loc->xgr[d][XSTEP]) > tol) 
	  reduced_dim[d] = loc->xgr[d][XLENGTH];
	else if (FABS(loc->xgr[d][XSTEP]) * loc->xgr[d][XLENGTH] <= tol)
	  reduced_dim[d] = 1;
	else {
	  err = INTERNAL_ERR;
	  warning("'%.50s' is larger than the grid step, but smaller than the whole grid in direction %d. This looks odd.", KNAME(NUGGET_PROC_TOL), d);	
	}
	prod_dim[d + 1] = prod_dim[d] * reduced_dim[d];	
      }
      if (err == NOERROR) {
	if ((s->red_field=(double *) MALLOC(sizeof(double) * vdim *
					    prod_dim[dim])) == NULL ||
	    (s->index = (int*) MALLOC(dim * sizeof(int))) == NULL
	    ){ err=ERRORMEMORYALLOCATION; goto ErrorHandling; }	
      }
    }

    if (!s->simugrid || err != NOERROR) {
      assert(!s->simugrid || err == INTERNAL_ERR);
      int i;
      if ((s->pos = pos = 
	   (int*) MALLOC(sizeof(int) * tot)) == NULL ||
	  (s->index = (int*) MALLOC(sizeof(int) * tot)) == NULL ||
	  (refpos = (int*) MALLOC(sizeof(int) * tot)) == NULL
	  ) {
	err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }

      TransformLoc(cov, false, True, true);

      //      PMI(cov);

      loc = Loc(cov);
      assert(!loc->grid && !loc->Time);
      Ext_ordering(loc->x, loc->totalpoints, dim, pos);
      
      //      for (i=0; i<loc->totalpoints; i++) {
      //	printf("old %10g %10g pos[i]=%d new=%10g %10g\n",
      //	       loc->x[0 + 2 * i], loc->x[1 + 2 * i], pos[i],
      //	       loc->x[0 + 2 * pos[i]], loc->x[1 + 2 * pos[i]]);
      //      }
      
      int oldrefpos,
	last = 0;
      i = 0;
      s->datapos[i] = last; // i in 0...totalpoints-1
      oldrefpos = refpos[last] = i; // last in 0..,[# different locations]
      for (i=1 /* ! */; i<loc->totalpoints; i++) {
	if (equal(next, pos[oldrefpos], pos[i], loc->x, dim)) {
	  s->datapos[i] = s->datapos[oldrefpos];
	} else {
	  bool ok = false;
	  int cur = last; // not last - 1
	  double x1value = loc->x[pos[i] * dim];
	  while(cur >= 0 && 
		(ok = FABS(loc->x[pos[refpos[cur]] * dim] - x1value) < tol) &&
		//                                        "== tol " not possible
		!equal(next, pos[ refpos[cur]], pos[i], loc->x, dim)) cur--;
	  if (cur<0 || !ok) {
	    // first x-coordinate is ordered, so if the distance is too big,
	    // i.e. if the refpos first x-coordinate is too small, there is
	    // no chance to find a close point
	    // so, in any case, a new location has been found
	    last++;
	    s->datapos[i] = last;
	    oldrefpos = refpos[last] = i;
	  } else { // a former reference point found:
	    oldrefpos = refpos[cur];
	    s->datapos[i] = s->datapos[oldrefpos];
	  }
	}
      }
      s->total = last + 1;
      if ((s->red_field=(double *) MALLOC(sizeof(double) * vdim *
					  s->total)) == NULL){
	err=ERRORMEMORYALLOCATION; goto ErrorHandling; 
      }	
    }
  } // ANISO given
  
  if ((err = ReturnOwnField(cov)) != NOERROR) goto ErrorHandling;

 ErrorHandling:
  FREE(A);
  FREE(value);
  FREE(ivalue);
  FREE(dummy);
  FREE(refpos);
  cov->simu.active = err == NOERROR;

  RETURN_ERR(err);
}
 
void do_nugget(model *cov, gen_storage VARIABLE_IS_NOT_USED *S) {
  location_type *loc = Loc(cov);
///  double sqrtnugget;
  nugget_storage* s = (nugget_storage*) cov->Snugget;
  long nx, endfor;
  int 
    vdim = VDIM0;
  double
    *field = s->red_field,
    *res = cov->rf;
  SAVE_GAUSS_TRAFO;
 
 
//  sqrtnugget = s->sqrtnugget;

  //PMI(cov->calling->calling);

  model *intern = cov->key != NULL ? cov->key : cov->sub[0];
  assert(intern != NULL);
  while (intern != NULL && isDollar(intern)) intern = intern->sub[0];
  assert(intern->Snugget != NULL);
  
  if (!intern->Snugget->spatialnugget) {
   for (nx=0, endfor=loc->totalpoints * vdim; nx<endfor; nx++) {
      res[nx] = (double) GAUSS_RANDOM(1.0);
    }
  } else {
    if (s->simugrid) {
      int d, dim, dimM1, *red_dim, *prod_dim;
      long totpnts, idx;
      totpnts = loc->totalpoints;
      dim = OWNTOTALXDIM;
      dimM1 = dim - 1;
      red_dim = s->reduced_dim;
      prod_dim = s->prod_dim;
      endfor = vdim * s->prod_dim[dim]; // not dim-1

      for (int i=0; i<endfor; i++) {
	field[i] = (double) GAUSS_RANDOM(1.0);
      }
      for (d=0; d<dim; s->index[d++] = 0);
      for (int i=0; i<totpnts; i++) {
	for(idx=d=0; d<dim; d++)
	  idx += (s->index[d] % red_dim[d]) * prod_dim[d];
	for (int v=0; v<vdim; v++) {
	  res[i + v] = field[idx + v];
	}
	d = 0; 
	(s->index[d])++; 
	while (d < dimM1 && s->index[d] >= loc->xgr[d][XLENGTH]) { 
	  assert(d<dim); // assert(d>=0);
	  s->index[d] = 0; 
	  d++;   
	  (s->index[d])++; 
	}
      }
    } else {
      int
	pts = Gettotalpoints(cov),
	vdimtotal = s->total * vdim;
      for (int i=0; i<vdimtotal; field[i++] = GAUSS_RANDOM(1.0));
      for (nx=0; nx<pts; nx++) {
	double 
	  *data = field + s->datapos[nx] * vdim,
	  *rp = res + s->pos[nx] * vdim;
	for (int v=0; v<vdim; v++) rp[v] = data[v];
      }
    }
  }
  BOXCOX_INVERSE;
}


int struct_nugget(model *cov, model VARIABLE_IS_NOT_USED **newmodel) {
  //  PMI(cov, "structhyper");
  if (cov->sub[0]->pref[Nugget] == PREF_NONE) {
    RETURN_ERR(ERRORPREFNONE);
  }
  if (!hasGaussMethodFrame(cov)) {
    SERR("type is not Gaussian.");
  }
  RETURN_NOERROR;
}
