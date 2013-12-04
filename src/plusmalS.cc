/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of auxiliary correlation functions 

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is gno error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nonsta     tionary models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2013 Martin Schlather

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
 
#include "RF.h"
#include "Covariance.h"
//#include <R_ext/Lapack.h>
//#include <R_ext/Applic.h>
//#include <R_ext/Utils.h>     
//#include <R_ext/BLAS.h> 





// $
void kappaS(int i, cov_model *cov, int *nr, int *nc){
  switch(i) {
  case DVAR : case DSCALE :
    *nr = *nc = 1;
    break;
  case DANISO :
    *nc = SIZE_NOT_DETERMINED;
    *nr = cov->xdimown;
    break;
  case DALEFT :
    *nr = SIZE_NOT_DETERMINED;
    *nc = cov->xdimown;
    break;
  case DPROJ : 
    *nr = SIZE_NOT_DETERMINED;
    *nc = 1;
    break;
  default : *nr = *nc = -1;
  }
}


#define DOLLAR_SUB 0
// simple transformations as variance, scale, anisotropy matrix, etc.  
void Siso(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[DOLLAR_SUB];
  int i,
    vdimSq = cov->vdim * cov->vdim;
  double y,
    **p = cov->p,
    *aniso=p[DANISO],
    *scale =p[DSCALE],
    var = p[DVAR][0];
  
  //  printf("\nSiso\n");

  y = *x;
  if (aniso != NULL) y = fabs(y * aniso[0]);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] 
      : (y == 0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  // letzteres darf nur passieren wenn dim = 1!!
  COV(&y, next, v);

  //printf("end isoS\n");

  // PMI(cov);
  //  PMI(cov->calling->calling);
  // print("$ %f %f %f %d\n", *x, var, v[0], vdimSq);

  for (i=0; i<vdimSq; i++) v[i] *= var; 

  // print("$ %f %f %f %d\n", *x, var, v[0], vdimSq);
  //  error("von wo wird aufgerufen, dass nan in var?? und loc:x not given?");

}
  

// simple transformations as variance, scale, anisotropy matrix, etc.  
void logSiso(double *x, cov_model *cov, double *v, double *sign){
  cov_model *next = cov->sub[DOLLAR_SUB];
  int i,
    vdimSq = cov->vdim * cov->vdim;
  double y, 
    **p = cov->p,
    *aniso=p[DANISO],
    *scale =p[DSCALE],
    logvar = log(p[DVAR][0]);
  
  y = *x;
  if (aniso != NULL) y = fabs(y * aniso[0]);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] 
      : (y == 0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  LOGCOV(&y, next, v, sign);
  for (i=0; i<vdimSq; i++) v[i] += logvar; 
}
 
void Sstat(double *x, cov_model *cov, double *v){
  logSstat(x, cov, v, NULL);
}

void logSstat(double *x, cov_model *cov, double *v, double *sign){
  cov_model *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO];
  double 
    **p = cov->p,
    var = p[DVAR][0],
    *scale =p[DSCALE], 
    *z = cov->Sdollar->z,
    *aniso=p[DANISO];
  int i,
    nproj = cov->nrow[DPROJ],
    vdimSq = cov->vdim * cov->vdim;
 
  if (nproj > 0) {
    int *proj = (int *)cov->p[DPROJ];
    if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(nproj*sizeof(double));

    if (scale == NULL || scale[0] > 0.0) {
      if (scale == NULL)  for (i=0; i<nproj; i++) z[i] = x[proj[i] - 1];
      else {
	double invscale = 1.0 / scale[0];
	for (i=0; i<nproj; i++) z[i] = invscale * x[proj[i] - 1];
      }
    } else {
      // projection and aniso may not be given at the same time
      for (i=0; i<nproj; i++)
	z[i] = (x[proj[i] - 1] == 0 && scale[0] == 0) ? 0.0 : RF_INF;
    } 
    if (sign==NULL) COV(z, next, v) else LOGCOV(z, next, v, sign);
  } else if (Aniso != NULL) {
    int dim = Aniso->vdim;
    if (z == NULL) 
      z = cov->Sdollar->z =(double*) MALLOC(dim * sizeof(double));
    FCTN(x, Aniso, z);
    if (sign == NULL) COV(z, next, v) else LOGCOV(z, next, v, sign);
  } else if (aniso==NULL && (scale == NULL || scale[0] == 1.0)) {
    if (sign==NULL) COV(x, next, v) else LOGCOV(x, next, v, sign);
  } else {
    int xdimown = cov->xdimown;
    double *xz;
    if (z == NULL) 
      z = cov->Sdollar->z =(double*) MALLOC(xdimown * sizeof(double)); 
    if (aniso!=NULL) {
      int j, k,
	nrow=cov->nrow[DANISO], 
	ncol=cov->ncol[DANISO];
      for (k=i=0; i<ncol; i++) {
	z[i] = 0.0;
	for (j=0; j<nrow; j++) {
	  z[i] += aniso[k++] * x[j];
	}
      }
      xz = z;
    } else xz = x;    
    if (scale != NULL) {
      if (scale[0] > 0.0) {
	double invscale = 1.0 / scale[0];
	for (i=0; i < xdimown; i++) z[i] = invscale * xz[i];
      } else {
	for (i=0; i < xdimown; i++)
	  z[i] = (xz[i] == 0.0 && scale[0] == 0.0) ? 0.0 : RF_INF;
      }
    }
    if (sign==NULL) COV(z, next, v) else LOGCOV(z, next, v, sign);
  }

  if (sign==NULL) {
    for (i=0; i<vdimSq; i++) v[i] *= var; 
  } else {
    double logvar = log(var);
    for (i=0; i<vdimSq; i++) v[i] += logvar; 
  }

}

void Snonstat(double *x, double *y, cov_model *cov, double *v){
  logSnonstat(x, y, cov, v, NULL);
}

void logSnonstat(double *x, double *y, cov_model *cov, double *v, 
		 double *sign){
  cov_model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO];
  double 
    *z1 = cov->Sdollar->z,
    *z2 = cov->Sdollar->z2,
    **p = cov->p, 
    var = p[DVAR][0],
    *scale =p[DSCALE],
    *aniso=p[DANISO];
  int i,
    nproj = cov->nrow[DPROJ],
    vdimSq = cov->vdim * cov->vdim;

  //  PMI(cov, "$$$$");
  
  if (nproj > 0) {
    int *proj = (int *)cov->p[DPROJ];
    if (z1 == NULL) 
      z1 = cov->Sdollar->z = (double*) MALLOC(nproj * sizeof(double));
    if (z2 == NULL) 
      z2 = cov->Sdollar->z2 = (double*) MALLOC(nproj * sizeof(double));
    if (scale==NULL || scale[0] > 0.0) {
      double invscale = scale==NULL ? 1.0 :  1.0 / scale[0];
      for (i=0; i<nproj; i++) {
	z1[i] = invscale * x[proj[i] - 1];
	z2[i] = invscale * y[proj[i] - 1];	
      }
    } else {
      double s = scale[0]; // kann auch negativ sein ...
      for (i=0; i<nproj; i++) {
	z1[i] = (x[proj[i] - 1] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
	z2[i] = (y[proj[i] - 1] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
      }
    }
  } else if (Aniso != NULL) {
    int dim = Aniso->vdim;
    if (z1 == NULL) 
      z1 = cov->Sdollar->z =(double*) MALLOC(dim * sizeof(double));
    if (z2 == NULL) 
      z2 = cov->Sdollar->z2 =(double*) MALLOC(dim * sizeof(double));
    FCTN(x, Aniso, z1);
    FCTN(y, Aniso, z2);
    if (sign == NULL) NONSTATCOV(z1, z2, next, v)
    else LOGNONSTATCOV(z1, z2, next, v, sign);
  } else if (aniso==NULL && (scale==NULL || scale[0] == 1.0)) {
    z1 = x;
    z2 = y;
  } else {
    int xdimown = cov->xdimown;
    double *xz1, *xz2;
    if (z1 == NULL) 
      z1 = cov->Sdollar->z = (double*) MALLOC(xdimown * sizeof(double));
    if (z2 == NULL) 
      z2 = cov->Sdollar->z2 = (double*) MALLOC(xdimown * sizeof(double));
    if (aniso != NULL) {
      int j, k,
	nrow=cov->nrow[DANISO],
	ncol=cov->ncol[DANISO];
      for (k=i=0; i<ncol; i++) {
	z1[i] = z2[i] =0.0;
	for (j=0; j<nrow; j++, k++) {
	  z1[i] += aniso[k] * x[j];
	  z2[i] += aniso[k] * y[j];
	}
      }
      xz1 = z1;
      xz2 = z2;
    } else {
      xz1 = x;
      xz2 = y;
    }
    if (scale != NULL) {
      double s = scale[0];
      if (s > 0.0) {
	double invscale = 1.0 / s;
	for (i=0; i<xdimown; i++) {
	  z1[i] = invscale * xz1[i];
	  z2[i] = invscale * xz2[i];
	}
      } else {
	for (i=0; i<nproj; i++) {
	  z1[i] = (xz1[i] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
	  z2[i] = (xz2[i] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
	}
      }
    }
  }
  
  if (sign == NULL) {
    NONSTATCOV(z1, z2, next, v);
    for (i=0; i<vdimSq; i++) v[i] *= var; 
  } else {
    double logvar = log(var);
    LOGNONSTATCOV(z1, z2, next, v, sign);
    for (i=0; i<vdimSq; i++) v[i] += logvar; 
  }
}


void covmatrixS(cov_model *cov, double *v) {
  location_type *loc = Loc(cov);	      
  cov_model *next = cov->sub[DOLLAR_SUB];
  location_type *locnext = Loc(next);
  assert(locnext != NULL);
  int i, err, tot, totSq,
    dim = loc->timespacedim,
    vdim = cov->vdim;
  double var = cov->p[DVAR][0];
    
  if ((err = alloc_cov(cov, dim, vdim)) != NOERROR) 
    error("memory allocation error in 'covmatrixS'");

  if ((cov->p[DSCALE] != NULL && cov->p[DSCALE][0] != 1.0) || 
       cov->p[DANISO] != NULL || cov->p[DPROJ] != 0) {
     CovarianceMatrix(cov, v); 
    return;
  }

  if (next->xdimprev != next->xdimown) {
    assert(false); // fuehrt zum richtigen Resultat, sollte aber nicht
    // vorkommen!
    CovarianceMatrix(cov, v); 
    return;
  }

  int next_gatter = next->gatternr,
    next_xdim = next->xdimprev;
 
  next->gatternr = cov->gatternr;
  next->xdimprev = cov->xdimprev;
  CovList[next->nr].covmatrix(next, v);//hier wird uU next->totalpoints gesetzt
  next->gatternr = next_gatter;
  next->xdimprev = next_xdim;

  // PMI(cov, "covmatrix S");
  if (locnext==NULL) loc_set(cov, locnext->totalpoints);
  tot = cov->vdim * locnext->totalpoints;
  totSq = tot * tot;
  if (var == 1.0) return;
  for (i=0; i<totSq; v[i++] *= var);
}

char iscovmatrixS(cov_model *cov) {
  cov_model *sub = cov->sub[DOLLAR_SUB];
  return (int) ((cov->p[DSCALE] == NULL || cov->p[DSCALE][0]==1.0) &&
		cov->kappasub[DANISO]==NULL &&
		sub->p[DPROJ] == NULL &&
		sub->p[DANISO] == NULL) * CovList[sub->nr].is_covariance(cov);
}

void DS(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DANISO] == NULL);
  int i,
    vdimSq = cov->vdim * cov->vdim,
    nproj = cov->nrow[DPROJ];
  double y[2], varSc,
    **p = cov->p,
    *scale =p[DSCALE],
    *aniso=p[DANISO],
    spinvscale = 1.0;

  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varSc = p[DVAR][0] * spinvscale;

  if (nproj == 0) {
    y[0] = x[0] * spinvscale; 
    y[1] = (cov->isoown==ISOTROPIC || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    // for (i=0; i<nproj; i++) {
    //   y[i] = spinvscale * x[proj[i] - 1];
  }

  Abl1(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varSc; 
}

void DDS(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DANISO] == NULL);
  int i,
    vdimSq = cov->vdim * cov->vdim,
    nproj = cov->nrow[DPROJ],
    *proj = (int *)cov->p[DPROJ];
  double y[2], varScSq,
    **p = cov->p,
    *scale =p[DSCALE],
    *aniso=p[DANISO],
    spinvscale = 1.0;
  
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScSq = p[DVAR][0] * spinvscale * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (cov->isoown==ISOTROPIC || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i] - 1] * x[proj[i] - 1];
    }
    y[0] = sqrt(y[0]) * spinvscale;
  }
  Abl2(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScSq; 
}


void D3S(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DANISO] == NULL);
  int i,
    vdimSq = cov->vdim * cov->vdim,
    nproj = cov->nrow[DPROJ],
    *proj = (int *)cov->p[DPROJ];
  double y[2], varScS3,
    **p = cov->p,
    *scale =p[DSCALE],
    *aniso=p[DANISO],
    spinvscale = 1.0;
 
  assert(false);
  
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScS3 = p[DVAR][0] * spinvscale * spinvscale * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (cov->isoown==ISOTROPIC || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i] - 1] * x[proj[i] - 1];
    }
    y[0] = sqrt(y[0]) * spinvscale;
  }
  Abl3(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScS3; 
}


void nablahessS(double *x, cov_model *cov, double *v, bool nabla){
  cov_model *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO];
  assert(cov->kappasub[DANISO] == NULL);
  int i, endfor,
    dim = cov->nrow[DANISO],// == ncol == x d i m ?
    xdimown = cov->xdimown,
    nproj = cov->nrow[DPROJ];
  double *xy, *vw,
    *y = cov->Sdollar->y,
    *z = cov->Sdollar->z,
    *w = cov->Sdollar->z2,
    **p = cov->p,
    *scale =p[DSCALE],
    *aniso=p[DANISO],
    var = p[DVAR][0];
  if (nproj != 0) BUG;
  if (Aniso != NULL) BUG;
  if (dim != xdimown) BUG;
	  
  
  if (aniso != NULL) {  
    if (z == NULL) 
      z = cov->Sdollar->z = (double*) MALLOC(sizeof(double) * xdimown);
     if (w == NULL) 
      w = cov->Sdollar->z2 = (double*) MALLOC(sizeof(double) * xdimown);
    xA(x, aniso, xdimown, xdimown, z);
    xy = z;
    vw = w;
  } else {
    xy = x;
    vw = v;
  }

  if (scale != NULL) {
    if (y == NULL) 
      y = cov->Sdollar->y =(double*) MALLOC(sizeof(double) * xdimown);
    assert(scale[0] > 0.0);
    double spinvscale = 1.0 / scale[0];
    var *= spinvscale;
    if (!nabla) var *= spinvscale; // gives spinvscale^2 in case of hess
    for (i=0; i<xdimown; i++) y[i] = xy[i] * spinvscale;
    xy = y;
  }

  endfor = xdimown;
  if (nabla) {
    NABLA(xy, next, vw);
  } else {
    HESSE(xy, next, vw);
    endfor *= xdimown;
  }
     
  if (aniso != NULL) {  
    if (nabla) Ax(aniso, vw, xdimown, xdimown, v);
    else XCXt(aniso, vw, v, xdimown, xdimown);
  }

  for (i=0; i<endfor; i++) v[i] *= var; 
}

void nablaS(double *x, cov_model *cov, double *v){
  nablahessS(x, cov, v, true);
}
void hessS(double *x, cov_model *cov, double *v){
  nablahessS(x, cov, v, false);
}


 
void inverseS(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[DOLLAR_SUB];
  if (cov->kappasub[DANISO] != NULL)
    ERR("inverse can only be calculated if 'Aniso' not an arbitrary function"); 
  int i,
    vdimSq = cov->vdim * cov->vdim,
    nproj = cov->nrow[DPROJ];
  //    *proj = (int *)cov->p[DPROJ];
  double y, 
    s = 1.0,
    **p = cov->p,
    *scale =p[DSCALE],
    *aniso=p[DANISO],
    var = p[DVAR][0];

  if (aniso != NULL) {
    s /= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) s *= scale[0];  
  if (nproj == 0) {
    y= *x / var; // inversion, so variance becomes scale
  } else {
    BUG;  //ERR("nproj is not allowed in invS");
  }
  
  if (CovList[next->nr].inverse == ErrCov) BUG;
  INVERSE(&y, next, v);
 
  for (i=0; i<vdimSq; i++) v[i] *= s; //!

}


//void coinitS(cov_model *cov, localinfotype *li) {
//  cov_model *next = cov->sub[DOLLAR_SUB];
//  if ( CovList[next->nr].coinit == NULL)
//    ERR("# cannot find coinit -- please inform author");
//  CovList[next->nr].coinit(next, li);
//}
//void ieinitS(cov_model *cov, localinfotype *li) {
//  cov_model *next = cov->sub[DOLLAR_SUB];
// 
//  if ( CovList[next->nr].ieinit == NULL)
//    ERR("# cannot find ieinit -- please inform author");
//  CovList[next->nr].ieinit(next, li);
//}

void tbm2S(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[DOLLAR_SUB];
  cov_fct *C = CovList + next->nr; // kein gatternr, da isotrop
  double y[2],  *xy,
    **p=cov->p, 
    *scale =p[DSCALE],
    *aniso = p[DANISO];
  assert(cov->kappasub[DANISO] == NULL);
  
  assert(cov->nrow[DPROJ] == 0);
  if (aniso!=NULL) {
    if (cov->ncol[DANISO]==2) {  // turning layers
      y[0] = x[0] * aniso[0]; // spatial 
      y[1] = x[1] * aniso[3]; // temporal
      assert(aniso[1] == 0.0 && aniso[2] == 0.0);
      if (y[0] == 0.0) C->cov(y, next, v); 
      else C->tbm2(y, next, v);
    } else {
      assert(cov->ncol[DANISO]==1);
      if (cov->nrow[DANISO] == 1) { // turning bands
	y[0] = x[0] * aniso[0]; // purely spatial 
	C->tbm2(y, next, v);
      } else { // turning layers, dimension reduction
	if (p[DANISO][0] == 0.0) {
	  y[0] = x[1] * aniso[1]; // temporal 
	  C->cov(y, next, v); 
	} else {
	  y[0] = x[0] * aniso[0]; // temporal 
	  C->tbm2(y, next, v);
	}
      }
    }
    xy = y;
  } else xy = x;

  if (scale != NULL) {
    double s = scale[0];
    if (s > 0) { 
      double invscale = 1.0 / s;
      if (cov->xdimown == 2){
	y[0] = xy[0] * invscale; // spatial 
	y[1] = xy[1] * invscale; // temporal
	if (y[0] == 0.0) C->cov(y, next, v); 
	else C->tbm2(y, next, v);
      } else {
	y[0] = xy[0] * invscale; // purely spatial 
	C->tbm2(y, next, v);
      }
    } else {
      y[0] = (s < 0.0 || xy[0] != 0.0) ? RF_INF : 0.0;
      if (cov->xdimown == 2)
	y[1] = (s < 0.0 || xy[1] != 0.0) ? RF_INF : 0.0;
      C->tbm2(y, next, v);
    }
  }
  *v *= p[DVAR][0];
}


// TODO : Aniso=Matrix: direkte implementierung in S,
// sodass nicht ueber initS gegangen werden muss, sondern
//  e  < -  e^\top Aniso


int checkS(cov_model *cov) {
  static bool print_warn_Aniso = true;

  // hier kommt unerwartet  ein scale == nan rein ?!!
  cov_model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO],
    *sub = cov->key == NULL ? next : cov->key;

  double  **p = cov->p;
  int i, err,
    xdimown = cov->xdimown,
    xdimNeu = xdimown,
    *proj = ((int *) cov->p[DPROJ]),
    nproj = cov->nrow[DPROJ];
  // bool skipchecks = GLOBAL.general.skipchecks;
  matrix_type type = TypeMany;


  assert(isAnyDollar(cov));
  if (!isDollarProc(cov)) cov->nr = DOLLAR; // wegen nr++ unten !
 
  // cov->q[1] not NULL then A has been given

  // if (cov->method == SpectralTBM && cov->xdimown != next->xdimprev)
  //  return ERRORDIM;

  if (cov->q == NULL && cov->p[DALEFT]!=NULL) {
    if (GLOBAL.warn.Aniso && print_warn_Aniso) {
      print_warn_Aniso = false;
      PRINTF("NOTE! Starting with RandomFields 3.0, the use of 'Aniso' is different from\n the former 'aniso' insofar that 'Aniso' is multiplied from the right by 'x'\n(whereas 'aniso' had been multiplied from the left by 'x').\nSet RFoptions(warn.Aniso=FALSE) to avoid this message.\n");
    }
    // here for the first time
    if (p[DANISO] != NULL) return ERRORANISO_T; 
    int j, k,
      lnrow = cov->nrow[DALEFT],
      lncol = cov->ncol[DALEFT],
      total = lncol * lnrow;
	
    double
      *pA = p[DALEFT]; 
    p[DANISO] = (double*) MALLOC(sizeof(double) * total);
    cov->nrow[DANISO] = lncol;
    cov->ncol[DANISO] = lnrow;
    for (i=k=0; i<lnrow; i++) {
      for (j=i; j<total; j+=lnrow) p[DANISO][k++] = pA[j];
    }
    free(cov->p[DALEFT]);
    cov->p[DALEFT]=NULL;
    cov->ncol[DALEFT] = cov->nrow[DALEFT] = 0;
    cov->q = (double*) CALLOC(1, sizeof(double));
    cov->qlen = 1;
  }
 
  if ((err = checkkappas(cov, false)) != NOERROR) {
    return err;
  }
  kdefault(cov, DVAR, 1.0);

  if (Aniso != NULL) {
    if (p[DANISO] != NULL || p[DPROJ] != NULL || p[DSCALE] != NULL) 
      SERR("if Aniso is an arbitrary function, only 'var' may be given as additional parameter"); 
    if (cov->isoown != SYMMETRIC && cov->isoown != NO_ROTAT_INV) 
      return ERRORANISO;

    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = true;
    if ((err = CHECK(Aniso, cov->tsdim, cov->xdimown, cov->typus, cov->domown,
		       cov->isoown, SUBMODEL_DEP, cov->role)) != NOERROR)
      return err;
    if (cov->key==NULL) {
      if ((err = CHECK(sub, Aniso->vdim, Aniso->vdim, cov->typus,  
			 cov->domown, 
			 cov->isoown, SUBMODEL_DEP, cov->role)) != NOERROR) {
	return err;
      }
    } 
    sub->pref[CircEmbed] = sub->pref[CircEmbedCutoff] = 
      sub->pref[CircEmbedIntrinsic] = sub->pref[TBM] = sub->pref[SpectralTBM] =
      sub->pref[Sequential] = sub->pref[Markov] = sub->pref[Average] =
      sub->pref[Nugget] = sub->pref[RandomCoin] = sub->pref[Hyperplane] 
      = sub->pref[Specific] = PREF_NONE;
 
  } else if (p[DANISO] != NULL) { // aniso given
    int *idx=NULL,
      nrow = cov->nrow[DANISO],
      ncol = cov->ncol[DANISO];
    bool quasidiag;

    idx = (int *) MALLOC((nrow > ncol ? nrow : ncol) * sizeof(int));
    if (nrow==0 || ncol==0) SERR("dimension of the matrix is 0");
    if (p[DPROJ] != NULL) return ERRORANISO_T;
    if (xdimown < nrow) {
      if (PL >= PL_ERRORS) {LPRINT("xdim=%d != nrow=%d\n", xdimown, nrow);}
      SERR("#rows of anisotropy matrix does not match dim. of coordinates");
    }
    if (xdimown != cov->tsdim && nrow != ncol)
      SERR("non-quadratic anisotropy matrices only allowed if dimension of coordinates equals spatio-temporal dimension");

    analyse_matrix(p[DANISO], nrow, ncol, 
		   &(cov->diag),
		   &quasidiag, // &(cov->quasidiag), 
		   idx, // cov->idx
		   &(cov->semiseparatelast),
		   &(cov->separatelast));
    free(idx); idx=NULL;
    type = Type(p[DANISO], nrow, ncol);
    
    cov->full_derivs = cov->rese_derivs 
      = (xdimown == 1 || nrow == ncol) ? 2 : 0;
    cov->loggiven = true;
  
    // printf("hiere\n");
 
    switch (cov->isoown) {
    case ISOTROPIC : 
      if (cov->tsdim != 1) return ERRORANISO;
      break;
    case SPACEISOTROPIC :  
      cov->full_derivs =  cov->rese_derivs = 2;
      if (nrow != 2 || !isMdiag(type))
	SERR("spaceisotropy needs a 2x2 diagonal matrix");
      break;      
    case ZEROSPACEISO : 
      if (!isMdiag(type)) return ERRORANISO;
      break;      
    case VECTORISOTROPIC :
      if (!isMiso(type)) return ERRORANISO; 
      break;
    case SYMMETRIC: 
      break;
    case PREVMODELI : BUG;      
    case NO_ROTAT_INV : 
      if (!isProcess(cov->typus)) return ERRORANISO;
      break;
    default : BUG;
    }

    //  printf(" A hiere\n");
    
    if (!cov->diag) 
       cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = cov->pref[Hyperplane] = PREF_NONE;
    if (cov->isoown != SPACEISOTROPIC && !isMiso(type))
      cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;
   
    if (cov->key==NULL) {
      if ((err = 
	   CHECK(next, ncol, ncol, cov->typus, cov->domown, 
		 ncol==1 && !isProcess(cov->typus) ? ISOTROPIC : cov->isoown, 
		 SUBMODEL_DEP, cov->role))
	  != NOERROR) {
	return err;
      }
      if (next->domown == cov->domown && next->isoown == cov->isoown &&
	  xdimown > 1) next->delflag = DEL_COV - 7;
    } else {
      if ((err =
	   CHECK(cov->key, ncol, ncol, cov->typus, cov->domown, 
		 ncol==1 && !isProcess(cov->typus) ? ISOTROPIC : cov->isoown,
		 SUBMODEL_DEP, cov->role)) != NOERROR) return err;
    }
  
  } else { // p[DANISO] == NULL
    
    int tsdim = cov->tsdim;
    if (nproj > 0) {
      if (cov->ncol[DPROJ] != 1) SERR("proj must be a vector");
      if (cov->xdimprev != cov->xdimown) return ERRORANISO;
      for (i=0; i<nproj; i++) {
	int idx = proj[i] - 1;
	if (idx >= xdimown)
	  SERR2("%d-th value of 'proj' (%d) out of range", i, proj[i]);
      }
      xdimNeu = nproj;
      tsdim = nproj;

      switch (cov->isoown) {
      case ISOTROPIC : return ERRORANISO;
      case SPACEISOTROPIC :  
	if (nproj != 2 || xdimown != 2 || tsdim != 2)
	  SERR("spaceisotropy needs a 2x2 diagonal matrix");
	break;      
      case ZEROSPACEISO : return ERRORANISO; // ginge z.T.; aber kompliziert
	break;      
      case VECTORISOTROPIC : return ERRORANISO; 
      case SYMMETRIC: case NO_ROTAT_INV: break;
      case PREVMODELI : BUG;      
      default : BUG;
      }
 
    }

    if (cov->key==NULL) {
      if ((err = CHECK(next, tsdim, xdimNeu, cov->typus, cov->domown,
		       cov->isoown, SUBMODEL_DEP, cov->role)) != NOERROR) {
	return err;
      }

      if (next->domown == cov->domown &&
	  next->isoown == cov->isoown) // darf nicht nach CHECK als allgemeine Regel ruebergezogen werden, u.a. nicht wegen stat/nicht-stat wechsel !!
	// hier geht es, da domown und isoown nur durchgegeben werden und die Werte      // bereits ein Schritt weiter oben richtig/minimal gesetzt werden.
	next->delflag = DEL_COV - 8;
    } else {
      if ((err = CHECK(cov->key, tsdim, xdimNeu, cov->typus,
			 cov->domown, cov->isoown,
		       SUBMODEL_DEP, cov->role)) != NOERROR) return err;
    }

  } // end no aniso
 
  if (( err = checkkappas(cov, false)) != NOERROR) {
    return err;
  }

   setbackward(cov, sub);
 	
   if (cov->isoown != ISOTROPIC && !isDollarProc(cov)) { // multivariate kann auch xdimNeu == 1 problematisch sein
    cov->nr++;
  }
  
  if (xdimNeu > 1) {
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 0;
  }

  // printf("cov-Nr = %d\n", cov->nr);
  
  // 30.10.11 kommentiert:
  //  cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 
  //      cov->pref[TBM] = cov->pref[SpectralTBM] = 0;
  if ( (p[DANISO]==NULL || isMiso(type)) && p[DPROJ]==NULL) {
    //    print("logspeed %s %f\n", NICK(sub), sub->logspeed);
    cov->logspeed = sub->logspeed * p[DVAR][0];
  }
  //////////////// 

  if (sub->pref[Nugget] > PREF_NONE) cov->pref[Specific] = 100;
  if (nproj == 0) cov->matrix_indep_of_x = sub->matrix_indep_of_x;

  if (p[DPROJ] == NULL && Aniso==NULL && p[DANISO] == NULL) {
    double scale = p[DSCALE] == NULL ? 1.0 : p[DSCALE][0];
    cov->taylorN = sub->taylorN;  
    for (i=0; i<cov->taylorN; i++) {
      cov->taylor[i][TaylorPow] = sub->taylor[i][TaylorPow];
      cov->taylor[i][TaylorConst] = sub->taylor[i][TaylorConst] *
	cov->p[DVAR][0] * pow(scale, -sub->taylor[i][TaylorPow]);   
    }

    cov->tailN = sub->tailN;  
    for (i=0; i<cov->tailN; i++) {
      cov->tail[i][TaylorPow] = sub->tail[i][TaylorPow];
      cov->tail[i][TaylorExpPow] = sub->tail[i][TaylorExpPow];
      cov->tail[i][TaylorConst] = sub->tail[i][TaylorConst] *
	cov->p[DVAR][0] * pow(scale, -sub->tail[i][TaylorPow]);   
      cov->tail[i][TaylorExpConst] = sub->tail[i][TaylorExpConst] *
	pow(scale, -sub->tail[i][TaylorExpPow]);
    }
  } else {
    cov->taylorN = cov->tailN = 0;
  }

  if (cov->Sdollar != NULL && cov->Sdollar->z != NULL)
    DOLLAR_DELETE(&(cov->Sdollar));
  if (cov->Sdollar == NULL) {
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
  } 
  assert(cov->Sdollar->z == NULL);
  
 if (isProcess(cov->typus)) 
    MEMCOPY(cov->pref, PREF_NOTHING, sizeof(pref_shorttype)); 
 
  return NOERROR;
}



bool TypeS(Types required, cov_model *cov) {
  cov_model *sub = cov->key==NULL ? cov->sub[0] : cov->key;

  // printf("\n\n\nxxxxxxxxxxxx %d %s %s\n\n", TypeConsistency(required, sub), TYPENAMES[required], NICK(sub));
  
  // if (required == ProcessType) crash(); //assert(false);

  return (required==TcfType || required==PosDefType || required==NegDefType
	  || required==ShapeType || required==TrendType 
	  || required==ProcessType || required==GaussMethodType)
    && TypeConsistency(required, sub);
}


void spectralS(cov_model *cov, storage *s, double *e) {
  cov_model *next = cov->sub[DOLLAR_SUB];
  double  
    **p = cov->p;
  int d,
    ncol = p[DANISO] == NULL ? cov->tsdim : cov->ncol[DANISO];
  double sube[MAXTBMSPDIM],
    *scale =p[DSCALE],
    invscale = 1.0;
  
  SPECTRAL(next, s, sube); // nicht gatternr

  // Reihenfolge nachfolgend extrem wichtig, da invscale auch bei aniso
  // verwendet wird


  if (scale != NULL) invscale /= scale[0];
  // print("sube %f %f %f %d %d\n", sube[0], sube[1], invscale, cov->tsdim, ncol);
  
  if (p[DANISO] != NULL) {
    int j, k, m,
      nrow = cov->nrow[DANISO],
      total = ncol * nrow;
    double
      *A = p[DANISO]; 
    for (d=0, k=0; d<nrow; d++, k++) {
      e[d] = 0.0;
      for (m=0, j=k; j<total; j+=nrow) {
	e[d] += sube[m++] * A[j] * invscale;
      }
    }
  } else { 
    for (d=0; d<ncol; d++) e[d] = sube[d] * invscale;
  }

}


void rangeS(cov_model *cov, range_type* range){
  int i;

  range->min[DVAR] = 0.0;
  range->max[DVAR] = RF_INF;
  range->pmin[DVAR] = 0.0;
  range->pmax[DVAR] = 100000;
  range->openmin[DVAR] = false;
  range->openmax[DVAR] = true;

  range->min[DSCALE] = 0.0;
  range->max[DSCALE] = RF_INF;
  range->pmin[DSCALE] = 0.0001;
  range->pmax[DSCALE] = 10000;
  range->openmin[DSCALE] = true;
  range->openmax[DSCALE] = true;

  for (i=DANISO; i<= DALEFT; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -10000;
    range->pmax[i] = 10000;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
  
  range->min[DPROJ] = 1;
  range->max[DPROJ] = cov->tsdim;
  range->pmin[DPROJ] = 1;
  range->pmax[DPROJ] =  cov->tsdim;
  range->openmin[DPROJ] = false;
  range->openmax[DPROJ] = false;
}


int structS(cov_model *cov, cov_model **newmodel) {
  //printf(" CALL OF STRUCT $  ****************\n");
  cov_model
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO];
  int err; 
  // cov_model *sub;



  //   printf("%s %d %s\n", ROLENAMES[cov->role], cov->role, TYPENAMES[cov->typus]);
  //  if (cov->role != ROLE_GAUSS) {
    //  APMI(cov->calling);
    //crash();
  //  }
  assert(cov->role == ROLE_GAUSS);

  switch (cov->role) {
  case ROLE_SMITH : 
    
    ASSERT_NEWMODEL_NOT_NULL;
     
    if (cov->p[DANISO] != NULL || cov->p[DALEFT])
       SERR("anisotropy parameter not allowed yet");
    if (cov->p[DPROJ] != NULL) SERR("projections not allowed yet");

    if ((err = STRUCT(next, newmodel)) > NOERROR) return err;

    addModel(newmodel, DOLLAR);
    if (cov->p[DVAR] != NULL) kdefault(*newmodel, DVAR, cov->p[DVAR][0]);
    if (cov->p[DSCALE] != NULL) kdefault(*newmodel, DSCALE, cov->p[DSCALE][0]);
    if (!next->deterministic) SERR("random shapes not programmed yet");
    
    break;
  case ROLE_GAUSS :
    if (isProcess(cov->typus)) {
      cov->nr = DOLLAR_PROC;
      return structSproc(cov, newmodel); // kein STRUCT !!
    } 

    ASSERT_NEWMODEL_NOT_NULL;
   if (newmodel == NULL) 
      SERR1("unexpected call to structure of '%s'", NICK(cov));
    if (cov->key != NULL) COV_DELETE(&(cov->key));

    if (cov->prevloc->distances) 
      SERR("distances do not allow for more sophisticated simulation methods");
    
    if (Aniso!= NULL) {
      SERR("complicated models including arbitrary functions for Aniso cannot be simulated yet");
    } 
    if ((err = STRUCT(next, newmodel)) > NOERROR) return err;

    addModel(newmodel, DOLLAR);
    if (cov->p[DVAR] != NULL) kdefault(*newmodel, DVAR, cov->p[DVAR][0]);
    if (cov->p[DSCALE] != NULL) kdefault(*newmodel, DSCALE, cov->p[DSCALE][0]);
    if (!next->deterministic) SERR("random shapes not programmed yet");
    
    break;
 
  default :
    //  PMI(cov, "structS");
    SERR1("changes in scale/variance not programmed yet for '%s'", 
	  ROLENAMES[cov->role]);      
  }
  
   
  return NOERROR;
}




int initS(cov_model *cov, storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!
  cov_model *next = cov->sub[DOLLAR_SUB];
  //mppinfotype *info = &(s->mppinfo);
  //location_type *loc = cov->prevloc;
  int 
    err = NOERROR;

  if (cov->role == ROLE_BROWNRESNICK || cov->role == ROLE_SMITH ||
      cov->role == ROLE_SCHLATHER || cov->role == ROLE_POISSON
      || cov->role == ROLE_POISSON_GAUSS) {
    double
      var = cov->p[DVAR][0],
      scale = cov->p[DSCALE] == NULL ? 1.0 : cov->p[DSCALE][0],
      // pow_scale = intpow(scale, cov->tsdim),
      *aniso = cov->p[DANISO],
      *anisoleft = cov->p[DALEFT],
      *proj = cov->p[DPROJ];
    cov_model
      *varM = cov->kappasub[DVAR],
      *scaleM = cov->kappasub[DSCALE],
      *anisoM = cov->kappasub[DANISO],
      *anisoleftM = cov->kappasub[DALEFT],
      *projM = cov->kappasub[DPROJ];
    int i,
      nm = cov->mpp.moments,
      dim = cov->tsdim;

    for (i=DVAR; i<=DMAX; i++) 
      if (cov->kappasub[i] != NULL) SERR("not programmed yet");

    if (aniso != NULL || proj != NULL || anisoleft != NULL ||
	anisoM!= NULL || projM!= NULL || anisoleftM!= NULL) 
      SERR("anisotropy and projection not allowed yet in Poisson related models");

    // Achtung INIT_RANDOM ueberschreibt mpp.* !!
    if (varM != NULL) {
      if ((err = INIT_RANDOM(varM, nm == 0 ? 1 : nm, s)) != NOERROR) return err;
      var = varM->mpp.M[1];      
    }
    if (scaleM != NULL) {
      if (dim + nm < 1) SERR("found dimension 0");
      if ((err = INIT_RANDOM(scaleM, dim + nm, s)) != NOERROR)
	return err;
      scale = scaleM->mpp.M[1];      
    }  
    if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) return err;

    if (varM != NULL) {
      for (i=0; i<=nm; i++) {
	cov->mpp.M[i] *= varM->mpp.M[i];
	cov->mpp.Mplus[i] *= varM->mpp.Mplus[i];
      }
    } else {
      double pow_var = 1.0;
      for (i=0; i<=nm; i++, pow_var *= var) {
	cov->mpp.M[i] *= pow_var;
	cov->mpp.Mplus[i] *= pow_var;
      }
    }
    if (scaleM != NULL) {
      for (i=0; i<=nm; i++) {
	cov->mpp.M[i] *= scaleM->mpp.M[i + dim];
	cov->mpp.Mplus[i] *= scaleM->mpp.Mplus[i + dim];
      }
    } else {
      double pow_scale = 1.0;
      for (i=0; i<=nm; i++, pow_scale *= scale) {
	cov->mpp.M[i] *= pow_scale;
	cov->mpp.Mplus[i] *= pow_scale;
      }
    }
 
    cov->mpp.maxheight = next->mpp.maxheight * var;
    //    cov->mpp.refradius *= scale;
    //cov->mpp.refsd *= scale;
    return NOERROR;
  }


  else if (cov->role == ROLE_GAUSS) {
    cov_model 
      *key = cov->key,
      *sub = key == NULL ? next : key;
    assert(sub != NULL);

    assert(key == NULL || ({PMI(cov);false;}));//
      
   
    if ((err=INIT(sub, 0, s)) != NOERROR) return err;

    return NOERROR;
    
  }

  SERR("initiation of scale failed");

}


void doS(cov_model *cov, storage *s){

  if (cov->role == ROLE_BROWNRESNICK || cov->role == ROLE_SMITH ||
      cov->role == ROLE_SCHLATHER || cov->role == ROLE_POISSON
      || cov->role == ROLE_POISSON_GAUSS) {
    cov_model *next = cov->sub[DOLLAR_SUB];
    
    cov_model
      *varM = cov->kappasub[DVAR],
      *scaleM = cov->kappasub[DSCALE];
 
    if (varM != NULL && !varM->deterministic) {
      assert(cov->p[DVAR] != NULL);
      VTLG_R(NULL, varM, cov->p[DVAR]);
    }

    if (scaleM != NULL && !scaleM->deterministic) {
      assert(cov->p[DSCALE] != NULL);
      VTLG_R(NULL, scaleM, cov->p[DSCALE]);
    }
    
    DO(next, s);// nicht gatternr
    cov->mpp.maxheight = next->mpp.maxheight * cov->p[DVAR][0];
    return;
  }
  
  if (cov->role == ROLE_GAUSS) {    
    double 
      *res = cov->rf,
      sd = sqrt(cov->p[DVAR][0]);
    int i,
      totalpoints = Gettotalpoints(cov);
    assert(res != NULL);

    DO(cov->key, s); 

    //    PMI(cov);
    //    printf("totalpoints");

    if (sd != 1.0) for (i=0; i<totalpoints; i++) res[i] *= (res_type) sd;
    return;
  }

  ERR("unknown option in doS "); 
}



int checkplusmal(cov_model *cov) {
  cov_model *sub;
  int i, j, err, dim, xdim, role;
  
  assert(cov->Splus == NULL);  
  dim = cov->tsdim;
  xdim = cov->xdimown;
  role = cov->role;
  


  for (i=0; i<cov->nsub; i++) {

    sub = cov->sub[i];
     
    if (sub == NULL) 
      SERR("+ or *: named submodels are not given in a sequence!");

    Types type = cov->typus;
    domain_type dom = cov->domown;
    isotropy_type iso = cov->isoown;

    //    printf("start\n");
    err = ERRORTYPECONSISTENCY;
    for (j=0; j<=1; j++) { // nur trend als abweichender typus erlaubt
      // 
      //printf("type = %s %s %d\n", TYPENAMES[type], NICK(sub), j);
      if (TypeConsistency(type, sub) &&
	  (err = CHECK(sub, dim, xdim, type, dom, iso, SUBMODEL_DEP, role))
	  == NOERROR) break;
      type = TrendType;
      dom = XONLY;
      iso = NO_ROTAT_INV;
    }
    //    printf("OK\n");

   if (err != NOERROR) {
      //printf("sub %d %s\n", i, NICK(sub));
      //      APMI(cov);
      return err;
    }

 
    if (cov->typus == sub->typus) {
      //printf("cov->iso = %d %d\n", cov->isoown, sub->isoprev);
      setbackward(cov, sub);
    } else {  
      updatepref(cov, sub);
      cov->tbm2num |= sub->tbm2num;
      if (CovList[cov->nr].vdim == SUBMODEL_DEP && 
	  (sub==cov->sub[0] || sub==cov->key)) {
	cov->vdim = sub->vdim;
	cov->vdim2[0] = sub->vdim2[0];
	cov->vdim2[1] = sub->vdim2[1];
      }
      cov->deterministic &= sub->deterministic;
    };

    if (i==0) {
      cov->vdim=sub->vdim; 
      cov->matrix_indep_of_x = sub->matrix_indep_of_x;
    } else {
      cov->matrix_indep_of_x &= sub->matrix_indep_of_x;
      if (cov->vdim != sub->vdim) {
	//      printf("i=%d vdims %d %d\n", i, cov->vdim, sub->vdim);
	SERR("multivariate dimensionality must be equal in the submodels");
      }
    }
    //    printf("i=%d vdims %d %d\n", i, cov->vdim, sub->vdim);
  }

  // !! incorrect  !!
  cov->semiseparatelast = false; 
  cov->separatelast = false; 

  //  PMI(cov, -1);

  return NOERROR;
}





// see private/old.code/operator.cc for some code including variable locations
void select(double *x, cov_model *cov, double *v) {
  int len,
    *element = ((int*) cov->p[SELECT_SUBNR]);
  cov_model *sub = cov->sub[*element];
  if (*element >= cov->nsub) error("select: element out of range");
  COV(x, sub, v);
  if ( (len = cov->nrow[SELECT_SUBNR]) > 1) {
    int i, m,
      vsq = cov->vdim * cov->vdim;
    double *z = cov->Sdollar->z;
    if (z == NULL) z = cov->Sdollar->z =(double*) MALLOC(sizeof(double) * vsq);
    for (i=1; i<len; i++) {
      sub = cov->sub[element[i]];
      COV(x, sub, z);
      for (m=0; m<vsq; m++) v[m] += z[m]; 
    }
  }
}
  

void covmatrix_select(cov_model *cov, double *v) {
  int len = cov->nrow[SELECT_SUBNR];
  
  if (len == 1) {
    int element = ((int*) cov->p[SELECT_SUBNR])[0];
    cov_model *next = cov->sub[element];
    if (element >= cov->nsub) error("select: element out of range");
    CovList[next->nr].covmatrix(next, v);

    // {  int i,j,k, tot=Loc(cov)->totalpoints; printf("\nXcovmat select\n");
    //  for (k=i=0; i<tot*tot; i+=tot) {
    //   for (j=0; j<tot; j++) printf("%f ", v[k++]);
    //  printf("\n");  }}

    // crash();
    //  PMI(next);
      
  }  else StandardCovMatrix(cov, v);
}

char iscovmatrix_select(cov_model VARIABLE_IS_NOT_USED *cov) {  return 2; }

int checkselect(cov_model *cov) {
  int err;

  assert(cov->Splus == NULL);

  kdefault(cov, SELECT_SUBNR, 0);
  if ((err = checkplus(cov)) != NOERROR) return err;

  if ((err = checkkappas(cov)) != NOERROR) return err;

  if (cov->Sdollar != NULL && cov->Sdollar->z != NULL)
    DOLLAR_DELETE(&(cov->Sdollar));
  if (cov->Sdollar == NULL) {
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
  } 
  assert(cov->Sdollar->z == NULL);

  return NOERROR;
}


void rangeselect(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SELECT_SUBNR] = 0;
  range->max[SELECT_SUBNR] = MAXSUB-1;
  range->pmin[SELECT_SUBNR] = 0;
  range->pmax[SELECT_SUBNR] = MAXSUB-1;
  range->openmin[SELECT_SUBNR] = false;
  range->openmax[SELECT_SUBNR] = false;
}




void plusStat(double *x, cov_model *cov, double *v){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  assert(cov->Sdollar != NULL);
  double *z = cov->Sdollar->z;
    if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(sizeof(double) * vsq);
  
  //PMI(cov->calling);
  //print("%d %s %s\n", vsq, NICK(cov->sub[0])),
  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    if (cov->typus == sub->typus) {
      COV(x, sub, z);
      for (m=0; m<vsq; m++) v[m] += z[m]; 
    }

    // if (cov->calling != NULL) printf("stat i=%d %f %f\n", i, v[m], z[m]);
    //  crash(cov);
    //  APMI(cov);
  }
  
  //  printf("plus x=%f %f\n", *x, *v);
}

void plusNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  assert(cov->Sdollar != NULL);
  double *z = cov->Sdollar->z;
  if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(sizeof(double) * vsq);
  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    if (cov->typus == sub->typus) {
      NONSTATCOV(x, y, sub, z);
      for (m=0; m<vsq; m++) v[m] += z[m]; 
    }

    //printf("i=%d %f %f\n", i, v[0], z[0]);
  }
  // printf("plus nonstat x=%f %f\n", *x, *v);
}

void Dplus(double *x, cov_model *cov, double *v){
  cov_model *sub;
  double v1;
  int n = cov->nsub, i;
  *v = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
    if (cov->typus == sub->typus) {
      Abl1(x, sub, &v1);
      (*v) += v1;
    }
  }
}

void DDplus(double *x, cov_model *cov, double *v){
  cov_model *sub;
  double v1;
  int n = cov->nsub, i;
  *v = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
    if (cov->typus == sub->typus) {
      Abl2(x, sub, &v1);
      (*v) += v1;
    }
  }
}


int checkplus(cov_model *cov) {
  int err, i;
  if ((err = checkplusmal(cov)) != NOERROR) {
    return err;
  }
  
  if (cov->domown == STAT_MISMATCH) return ERRORNOVARIOGRAM;
  if (cov->nsub == 0) cov->pref[SpectralTBM] = PREF_NONE;

  if (isPosDef(cov) && cov->domown == XONLY) cov->logspeed = 0.0;
  else if (isNegDef(cov) && cov->domown == XONLY) {
    cov->logspeed = 0.0;
    for (i=0; i<cov->nsub; i++) {
      cov_model *sub = cov->sub[i];
      if (cov->typus == sub->typus) {
	if (ISNAN(sub->logspeed)) {
	  cov->logspeed = RF_NAN;
	  break;
	} else cov->logspeed += sub->logspeed;
      }
    }
  } else cov->logspeed = RF_NAN;

 if (cov->Sdollar != NULL && cov->Sdollar->z != NULL)
    DOLLAR_DELETE(&(cov->Sdollar));
  if (cov->Sdollar == NULL) {
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
  } 
  assert(cov->Sdollar->z == NULL);
  return NOERROR;

  // spectral mit "+" funktioniert, falls alle varianzen gleich sind,
  // d.h. nachfolgend DOLLAR die Varianzen gleich sind oder DOLLAR nicht
  // auftritt; dann zufaellige Auswahl unter "+"
}


bool Typeplus(Types required, cov_model *cov) {
  bool allowed = TypeConsistency(ShapeType, required) || required==TrendType;
    //||required==ProcessType ||required==GaussMethodType; not yet allowed;to do
  if (!allowed) return false;
  int i;
  for (i=0; i<cov->nsub; i++) {
    //printf("type plus %s %s\n",TYPENAMES[required], NICK(cov->sub[i]));
    if (TypeConsistency(required, cov->sub[i])) return true;
  }  
  return false;
}

void spectralplus(cov_model *cov, storage *s, double *e){
  int nr;
  double dummy;
  spec_properties *cs = &(s->spec);
  double *sd_cum = cs->sub_sd_cum;

  nr = cov->nsub - 1;
  dummy = UNIFORM_RANDOM * sd_cum[nr];
  if (ISNA(dummy)) BUG;
  while (nr > 0 && dummy <= sd_cum[nr - 1]) nr--;
  cov_model *sub = cov->sub[nr];
  SPECTRAL(sub, s, e);  // nicht gatternr
}


int structplus(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  int m, err;
  
  
  switch(cov->role) {
  case ROLE_COV :  return NOERROR;
  case ROLE_GAUSS : 
    if (isProcess(cov->typus)) {
      assert(cov->calling != NULL && (isInterface(cov->calling->typus) ||
				      isProcess(cov->calling->typus)));
      assert(cov->nr == PLUS_PROC || cov->nr == PLUS);
      cov->nr = PLUS_PROC;
      return structplusproc(cov, newmodel); // kein STRUCT !!
    } 
    if (cov->Splus != NULL) BUG;
    for (m=0; m<cov->nsub; m++) {
      cov_model *sub = cov->sub[m];
      if ((err = STRUCT(sub, newmodel))  > NOERROR) {
      //	PMI(cov->Splus->keys[m]);
	//	assert(false);
	//printf("end plus\n");
	return err;
      }
    }
    break;
  default :
    SERR2("role '%s' not allowed for '%s'", ROLENAMES[cov->role],
	  NICK(cov));
  }

  return NOERROR;
}

int initplus(cov_model *cov, storage *s){
  int i, err;

  cov->mpp.maxheight = RF_NAN;

  if (cov->role == ROLE_GAUSS) {
    spec_properties *cs = &(s->spec);
    double *sd_cum = cs->sub_sd_cum;
 
    for (i=0; i<cov->nsub; i++) {
      cov_model *sub = cov->Splus == NULL ? cov->sub[i] : cov->Splus->keys[i];

      //print("initplus %d %d initialised=%d\n",i,cov->nsub, sub->initialised);
      //PMI(sub);
      
      if (sub->pref[Nothing] > PREF_NONE) { // to do ??
	// for spectral plus only
	COV(ZERO, sub, sd_cum + i);
	if (i>0) sd_cum[i] += sd_cum[i-1];
      }
      cov->sub[i]->stor = (storage *) MALLOC(sizeof(storage));
      if (!sub->initialised) {
	if ((err = INIT(sub, cov->mpp.moments, s)) != NOERROR) {
	  //  AERR(err);
	  return err;
	}
      }
      sub->simu.active = true;
    }

    // assert(false);
  
    cov->fieldreturn = cov->Splus != NULL;
    cov->origrf = false;
    if (cov->Splus != NULL) cov->rf = cov->Splus->keys[0]->rf;
     
    return NOERROR;
  }

  /*
    pref_type pref;

    for (m=0; m<Forbidden; m++) pref[m] = PREF_BEST;
    if (meth->cov->user[0] == 0 || meth->cov->user[1] == 0) {
 
    // i.e. user defined
    pref_type pref;
    pref[CircEmbed] = pref[TBM] = pref[Direct] = pref[Sequential] 
    = PREF_NONE;
    for (m=0; m<Nothing; m++) { // korrekt auch fuer MaxStable?
    //	  print("%d %d\n", m, pref[m]);
    //	  print("%d %d\n", meth->cov->user[m]);
    if (pref[m] > 0 &&  meth->cov->user[m] > 0) {
    break;
    }
    }
    if (m == Nothing) return ERRORSUBMETHODFAILED;
    }
  */
     
  else if (cov->role == ROLE_COV) {    
    return NOERROR;
  }

  return ERRORFAILED;
}


void doplus(cov_model *cov, storage *s) {
  int i;
  //total = cov->prevloc->totalpoints * cov->vdim;
  //double *res = cov->rf;
  
  if (cov->role == ROLE_GAUSS && cov->method==SpectralTBM) {
    ERR("error in doplus with spectral");
  }
  
  for (i=0; i<cov->nsub; i++) {
    cov_model *sub = cov->Splus==NULL ? cov->sub[i] : cov->Splus->keys[i];
    DO(sub, s);
  }
}




void covmatrix_plus(cov_model *cov, double *v) {
  location_type *loc = Loc(cov);
  //  cov_fct *C = CovList + cov->nr; // nicht gatternr
  int i, 
    totalpoints = loc->totalpoints,
    vdimtot = totalpoints * cov->vdim,
    vdimtotSq = vdimtot * vdimtot,
    nsub = cov->nsub;
  bool is = iscovmatrix_plus(cov) >= 2;
  double *mem = NULL;
  if (is && nsub>1) {
    mem = cov->Sdollar->y;
    if (mem == NULL) 
      mem = cov->Sdollar->y = (double*) MALLOC(sizeof(double) * vdimtotSq);
    is = mem != NULL;
  }
  
  if (is) {
    // APMI(cov);

    //cov_model *sub = cov->sub[0];
    int j;
    if (cov->p[SELECT_SUBNR] == NULL) 
      cov->p[SELECT_SUBNR] = (double*) MALLOC(sizeof(double));
    cov->p[SELECT_SUBNR][0] = 0;
    CovList[SELECT].covmatrix(cov, v);
    for (i=1; i<nsub; i++) {
      if (Loc(cov->sub[i])->totalpoints != totalpoints) BUG;
      cov->p[SELECT_SUBNR][0] = i;
      CovList[SELECT].covmatrix(cov, mem);
      for (j=0; j<vdimtotSq; j++) v[j] += mem[j];
    }
  } else StandardCovMatrix(cov, v);  
}

char iscovmatrix_plus(cov_model *cov) {
  char max=0, is;
  int i,
    nsub = cov->nsub;
  for (i=0; i<nsub; i++) {
    cov_model *sub = cov->sub[i];
    is = CovList[sub->nr].is_covmatrix(sub);
    if (is > max) max = is;
  }
  return max;
}


void malStat(double *x, cov_model *cov, double *v){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  assert(cov->Sdollar != NULL);
  double *z = cov->Sdollar->z;
  if (z == NULL) z = cov->Sdollar->z =(double*) MALLOC(sizeof(double) * vsq);
  
  assert(x[0] >= 0.0 || cov->xdimown > 1);
  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    COV(x, sub, z);
    for (m=0; m<vsq; m++) v[m] *= z[m]; 
  }
}

void logmalStat(double *x, cov_model *cov, double *v, double *sign){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  double *z = cov->Sdollar->z,
    *zsign = cov->Sdollar->z2;
  if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(sizeof(double) * vsq);
  if (zsign == NULL) 
    zsign = cov->Sdollar->z2 = (double*) MALLOC(sizeof(double) * vsq);
  
  assert(x[0] >= 0.0 || cov->xdimown > 1);
  for (m=0; m<vsq; m++) {v[m] = 0.0; sign[m]=1.0;}
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    LOGCOV(x, sub, z, zsign);
    for (m=0; m<vsq; m++) {
      v[m] += z[m]; 
      sign[m] *= zsign[m];
    }
  }
}

void malNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *sub;
  int i, m, nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  assert(cov->Sdollar != NULL);
  double *z = cov->Sdollar->z;
  if (z == NULL) z = cov->Sdollar->z =(double*) MALLOC(sizeof(double) * vsq);
  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    NONSTATCOV(x, y, sub, z);
    for (m=0; m<vsq; m++) v[m] *= z[m]; 
  }
}

void logmalNonStat(double *x, double *y, cov_model *cov, double *v, 
		   double *sign){
  cov_model *sub;
  int i, m, nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  assert(cov->Sdollar != NULL);
  double *z = cov->Sdollar->z,
    *zsign = cov->Sdollar->z2;
  if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(sizeof(double) * vsq);
  if (zsign == NULL) 
    zsign = cov->Sdollar->z2 = (double*) MALLOC(sizeof(double) * vsq);
  for (m=0; m<vsq; m++) {v[m] = 0.0; sign[m]=1.0;}
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    LOGNONSTATCOV(x, y, sub, z, zsign);
    for (m=0; m<vsq; m++) {
      v[m] += z[m]; 
      sign[m] *= zsign[m];
    }
  }
}

void Dmal(double *x, cov_model *cov, double *v){
  cov_model *sub;
  double c[MAXSUB], d[MAXSUB];
  int n = cov->nsub, i;
  for (i=0; i<n; i++) {
    sub = cov->sub[i];
    COV(x, sub, c + i);
    Abl1(x, sub, d + i);
  }
  *v = 0.0;
  for (i=0; i<n; i++) {
    double zw = d[i];
    int j;
    for (j=0; j<n; j++) if (j!=i) zw *= c[j]; 
    (*v) += zw; 
  }
}
int checkmal(cov_model *cov) {
  cov_model *next1 = cov->sub[0];
  cov_model *next2 = cov->sub[1];
  int err;

  if (next2 == NULL) next2 = next1;
  
  if ((err = checkplusmal(cov)) != NOERROR) return err;

  if (cov->domown == STAT_MISMATCH || !isPosDef(cov)) {
    return ERRORNOVARIOGRAM;
  }
  cov->logspeed = cov->domown == XONLY ? 0 : RF_NAN;
  
  if (cov->xdimown >= 2) cov->pref[TBM] = PREF_NONE;
  if (cov->xdimown==2 && cov->nsub == 2 && 
      isAnyDollar(next1) && isAnyDollar(next2)) {
    double *aniso1 = next1->p[DANISO],
      *aniso2 = next2->p[DANISO];
    if (aniso1 != NULL && aniso2 != NULL) {
      if (aniso1[0] == 0.0 && next1->ncol[DANISO] == 1) {
	cov->pref[TBM] = next2->pref[TBM];
      } else if (aniso2[0] == 0.0 && next2->ncol[DANISO] == 1) {
	cov->pref[TBM] = next1->pref[TBM];
      }
    }
  }
 if (cov->Sdollar != NULL && cov->Sdollar->z != NULL)
    DOLLAR_DELETE(&(cov->Sdollar));
  if (cov->Sdollar == NULL) {
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
  } 
  assert(cov->Sdollar->z == NULL);
  return NOERROR;
}

bool Typemal(Types required, cov_model *cov) {
  bool allowed = required==PosDefType || required==ShapeType;
  //||required==NegDefType ||required==TcfType; to do; ist erlaubt unter
  // gewissen Bedingungen
  if (!allowed) return false;
  int i;
  for (i=0; i<cov->nsub; i++) {
    if (!TypeConsistency(required, cov->sub[i])) return false;
  }
  return true;
}


int initmal(cov_model *cov, storage VARIABLE_IS_NOT_USED *s){
//  int err;
//  return err;
  return ERRORFAILED;
  cov->mpp.maxheight = RF_NAN;
}
void domal(cov_model VARIABLE_IS_NOT_USED *cov, storage VARIABLE_IS_NOT_USED *s){
  assert(false);
}


//////////////////////////////////////////////////////////////////////
// mpp-plus


#define PLUS_P 0 // parameter
int  CheckAndSetP(cov_model *cov){   
  int n,
    nsub = cov->nsub;
  double 
    cump = 0.0;
  if (cov->p[PLUS_P] == NULL) {
    assert(cov->nrow[PLUS_P] == 0 && cov->ncol[PLUS_P] == 0);
    cov->p[PLUS_P] = (double*) MALLOC(sizeof(double) * nsub);
    cov->nrow[PLUS_P] = nsub;
    cov->ncol[PLUS_P] = 1;
    for (n=0; n<nsub; n++) cov->p[PLUS_P][n] = 1.0 / (double) nsub;    
  } else {
    cump = 0.0;
    for(n = 0; n<nsub; n++) {
      cump += cov->p[PLUS_P][n];
      //printf("cump %f %f\n",  cump , cov->p[PLUS_P][n]);
      if (cump > 1.0 && n+1<nsub) return ERRORATOMP; 
    }
    if (cump != 1.0) {
      if (nsub == 1) {
	warning("the p-values do not sum up to 1.\nHere only one p-value is given which must be 1.0");
	cov->p[PLUS_P][0] = 1.0;
      } else {
	//printf("%e\n", 1-cov->p[PLUS_P][nsub-1] );
	if (cump < 1.0 && cov->p[PLUS_P][nsub-1]==0)
	  warning("The value of the last component of 'p' is increased."); 
	else SERR("The components of 'p' do not sum up to 1.");
	cov->p[PLUS_P][nsub-1] = 1.0 - (cump - cov->p[PLUS_P][nsub-1]);
      }  
    }
  }
  return NOERROR;
}

void kappamppplus(int i, cov_model *cov, int *nr, int *nc){
  *nr = cov->nsub;
  *nc = i < CovList[cov->nr].kappas ? 1 : -1;
}

void mppplus(double *x, cov_model *cov, double *v) { 
  int i, n,
    vdimSq = cov->vdim * cov->vdim;
  double *z = cov->Sdollar->z;
  if (z == NULL) z = cov->Sdollar->z =(double*) MALLOC(sizeof(double) * vdimSq);
  cov_model *sub;

  if (cov->role == ROLE_COV) {  
    for (i=0; i<vdimSq; i++) v[i] = 0.0;
    for (n=0; n<cov->nsub; n++, sub++) {
      sub = cov->sub[n];
      COV(x, sub, z); // urspruenglich : covlist[sub].cov(x, cov, z); ?!
      for (i=0; i<vdimSq; i++) v[i] += cov->p[PLUS_P][n] * z[i];
    }
  } else {
    assert(cov->role == ROLE_POISSON || cov->role == ROLE_POISSON_GAUSS);
    assert(false);
  }
}

int checkmppplus(cov_model *cov) {
  int err, 
    size = 1;
  cov->maxdim = MAXMPPDIM;
  
  if ((err = checkplusmal(cov)) != NOERROR) {
    // printf("checkmppplus error %d %s\n", err, ERRORSTRING);
    return err;
  }

  if ((err = CheckAndSetP(cov)) != NOERROR) return(err);

  if (cov->q == NULL) {
    if ((cov->q  = (double*) CALLOC(sizeof(double), size)) == NULL)
      return ERRORMEMORYALLOCATION;
    cov->qlen = size;
  }
 if (cov->Sdollar != NULL && cov->Sdollar->z != NULL)
    DOLLAR_DELETE(&(cov->Sdollar));
  if (cov->Sdollar == NULL) {
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
  } 
  assert(cov->Sdollar->z == NULL);
  
  return NOERROR;
}



void rangempplus(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[PLUS_P] = 0.0;
  range->max[PLUS_P] = 1.0;
  range->pmin[PLUS_P] = 0.0;
  range->pmax[PLUS_P] = 1.0;
  range->openmin[PLUS_P] = false;
  range->openmax[PLUS_P] = false;
}

int struct_mppplus(cov_model *cov, cov_model **newmodel){
  int m,
    //nsub = cov->nsub,
    err = NOERROR;

  //if (cov->calling == NULL || cov->calling->ownloc == NULL) 
  // SERR("mppplus does not seem to be used in the right context");
  
  if (cov->role != ROLE_BROWNRESNICK && cov->role != ROLE_SMITH &&
      cov->role != ROLE_SCHLATHER && cov->role != ROLE_POISSON) {
    SERR("method is not based on Poisson point process");
  }

  SERR("'mppplus not programmed yet");
  // Ausnahme: mppplus wird separat behandelt:
  // if (nr == MPPPLUS) return STRUCT(shape, NULL);
 
  if (newmodel == NULL) BUG;  
  if (cov->Splus == NULL) {
    cov->Splus = (plus_storage*) MALLOC(sizeof(plus_storage));
    PLUS_NULL(cov->Splus);
  }

  for (m=0; m<cov->nsub; m++) {
    cov_model *sub = cov->sub[m];
          
    if (cov->Splus->keys[m] != NULL) COV_DELETE(cov->Splus->keys + m);    
    if ((err = covcpy(cov->Splus->keys + m, sub)) != NOERROR) return err;
    if ((err = addShapeFct(cov->Splus->keys + m)) != NOERROR) return err;
    cov->Splus->keys[m]->calling = cov;
  }  
  return NOERROR;
}


int init_mppplus(cov_model *cov, storage *S) {
  cov_model  *sub;
  double M2, M2plus, Eplus, maxheight;
  int n, err;
  ext_bool 
    loggiven = SUBMODEL_DEP,
    fieldreturn = SUBMODEL_DEP;
  pgs_storage *pgs = NULL;
  M2 = M2plus = Eplus = 0.0;
  maxheight = RF_NEGINF;
  
  if (cov->Spgs != NULL) PGS_DELETE(&(cov->Spgs));
  if ((pgs = cov->Spgs = (pgs_storage*) MALLOC(sizeof(pgs_storage))) == NULL) 
    return ERRORMEMORYALLOCATION;
  PGS_NULL(pgs);
  pgs->totalmass = 0.0;
  
  for (n=0; n<cov->nsub; n++) {
    sub = cov->sub[n];
    //if (!sub->mpp.loc_done) 
    //SERR1("submodel %d of '++': the location of the shapes is not defined", 
    // n);
    if ((err = INIT(sub, cov->mpp.moments, S)) != NOERROR) return err;
    //if (!sub->mpp.loc_done) SERR("locations are not initialised");
    if (n==0) {
      loggiven = sub->loggiven;
      fieldreturn = sub->fieldreturn;
    } else {
      if (loggiven != sub->loggiven) loggiven = SUBMODEL_DEP;
      if (fieldreturn != sub->loggiven) fieldreturn = SUBMODEL_DEP;
    }
    pgs->totalmass += sub->Spgs->totalmass * cov->p[PLUS_P][n];
    if (cov->mpp.maxheight > maxheight) maxheight = cov->mpp.maxheight;
    loggiven &= cov->loggiven;

    // Achtung cov->mpp.M2 und Eplus koennten nicht direkt gesetzt
    // werden, da sie vom CovList[sub->nr].Init ueberschrieben werden !!
    
    if (cov->mpp.moments >= 1) {
      Eplus += sub->p[PLUS_P][0] * sub->mpp.Mplus[1]; 
      if (cov->mpp.moments >= 2) {
	M2 += sub->p[PLUS_P][0] * sub->mpp.M[2];
	M2plus += sub->p[PLUS_P][0] * sub->mpp.M[2];
      }
    }
    //assert(sub->mpp.loc_done);
  }
  cov->mpp.maxheight = maxheight;
  //cov->mpp.refsd = RF_NAN;
  //  cov->mpp.refradius = RF_NAN;

  
  if (cov->mpp.moments >= 1) {
    cov->mpp.Mplus[1] = Eplus;
    cov->mpp.M[1] = RF_NAN;
     if (cov->mpp.moments >= 2) {
       cov->mpp.M[2] = M2;
       cov->mpp.Mplus[2] = M2plus;
     }
  }
  
  cov->loggiven = loggiven;
  cov->fieldreturn = fieldreturn;
  //cov->mpp.loc_done = true;       
  cov->origrf = false;
  cov->rf = NULL;

  return NOERROR;
}

void do_mppplus(cov_model *cov, storage *s) {
  cov_model *sub;
  double subselect = UNIFORM_RANDOM;
  int subnr;
  for (subnr=0; (subselect -= cov->sub[subnr]->p[PLUS_P][0]) > 0; subnr++);
  cov->q[0] = (double) subnr;
  sub = cov->sub[subnr];
  
  CovList[sub->nr].Do(sub, s);  // nicht gatternr
  cov->mpp.maxheight = sub->mpp.maxheight;
  cov->fieldreturn = sub->fieldreturn;
  cov->loggiven = sub->loggiven;
  // PrintMPPInfo(cov, s);
}

//////////////////////////////////////////////////////////////////////
// PROCESSES
//////////////////////////////////////////////////////////////////////


int structSproc(cov_model *cov, cov_model **newmodel) {
  //printf(" CALL OF STRUCT $  ****************\n");
  cov_model *key,
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO];
  int dim, err; 
  // cov_model *sub;

  //   printf("%d\n", cov->role); assert(cov->role == ROLE_GAUSS);
  
  assert(isDollarProc(cov));

  switch (cov->role) {
  case ROLE_GAUSS :
    if (newmodel != NULL) 
      SERR1("unexpected call to structure of '%s'", NICK(cov));
    if (cov->key != NULL) COV_DELETE(&(cov->key));

    if (cov->prevloc->distances) 
      SERR("distances do not allow for more sophisticated simulation methods");
    
    if (Aniso!= NULL) {
      SERR("complicated models including arbitrary functions for Aniso cannot be simulated yet");
    } else Transform2NoGrid(cov, true, false); 

    // printf("%f %f %f\n",
	   //	   Loc(cov)->xgr[0][0], Loc(cov)->xgr[0][1], Loc(cov)->xgr[0][2]);
    // assert(false);
 
    if ((err = covcpy(&(cov->key), next)) != NOERROR) return err;

    if (!isGaussProcess(next)) addModel(&(cov->key), GAUSSPROC);
    SetLoc2NewLoc(cov->key, Loc(cov));
    //   APMI(cov);assert(false);
    
    key = cov->key;
    assert(key->calling == cov);    
    
    dim = key->prevloc->timespacedim;
    if ((err = CHECK(key, dim, dim, ProcessType, XONLY, NO_ROTAT_INV,
		     cov->vdim, cov->role)) != NOERROR) return err;
    err = STRUCT(cov->key, NULL);

    //    APMI(cov);

    return err;
  default :
    //  PMI(cov, "structS");
    SERR1("changes in scale/variance not programmed yet for '%s'", 
	  ROLENAMES[cov->role]);      
  }
     
  return NOERROR;
}




int initSproc(cov_model *cov, storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!

  //cov_model *next = cov->sub[DOLLAR_SUB];
  //mppinfotype *info = &(s->mppinfo);
  //  location_type *loc = cov->prevloc;
  int 
    err = NOERROR;
  
  cov_model 
    *key = cov->key;
    //*sub = key == NULL ? next : key;
  location_type *prevloc = cov->prevloc;

  assert(key != NULL);
  
  if ((err = INIT(key, 0, s)) != NOERROR) {
    return err;
  }
  
  key->simu.active = true; 
  assert(s != NULL);

  cov->fieldreturn = true;

  if ((cov->origrf = cov->ownloc != NULL &&
       cov->ownloc->totalpoints != prevloc->totalpoints)) {
    assert(prevloc->grid);
    int dim = prevloc->timespacedim;
    cov->rf = (res_type*) MALLOC(sizeof(res_type) *
				 cov->vdim * 
				 prevloc->totalpoints);
    if (cov->Sdollar != NULL) DOLLAR_DELETE(&(cov->Sdollar));
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
 
    int d,
      *proj = (int *)cov->p[DPROJ],
      bytes = dim * sizeof(int),
      *cumsum = cov->Sdollar->cumsum = (int*) MALLOC(bytes),
      *total = cov->Sdollar->total = (int*) MALLOC(bytes),
      *len = cov->Sdollar->len = (int*) MALLOC(bytes);     
    cov->Sdollar->nx = (int*) MALLOC(bytes); 
    
    for (d=0; d<dim; d++) {
      cumsum[d] = 0;
      len[d] = prevloc->xgr[d][XLENGTH];
    }
    if (proj != NULL) {
      int 
	nproj = cov->nrow[DPROJ];
      cumsum[proj[d] - 1] = 1;
      for (d = 1; d < nproj; d++) {
	cumsum[proj[d] - 1] =
	  cumsum[proj[d - 1] - 1] * prevloc->xgr[d - 1][XLENGTH];
      }
      for (d=0; d<dim;d++) 
	total[d] = cumsum[d] * prevloc->xgr[d][XLENGTH];
    } else {
       int i,
	 iold = 0,
	 nrow = cov->nrow[DANISO],
	 ncol = cov->ncol[DANISO];
      double *A = cov->p[DANISO];
      for (d=0; d<ncol; d++, A += nrow) {
	for (i = 0; i < nrow && A[i] == 0.0; i++);
	if (i == nrow) i = nrow - 1;
	if (d > 0) {
	  cumsum[i] = cumsum[iold] * prevloc->xgr[d - 1][XLENGTH];
	} else { // d ==0
	  cumsum[i] = 1;
	}
	iold = i;
	for (i++; i < nrow; i++) if (A[i] != 0.0) BUG;  // just a check
      }
    }
  } else {
    cov->rf = cov->key->rf;
  }
 
  //  PMI(cov,-1); assert(false);
 
  return NOERROR;
}


void doSproc(cov_model *cov, storage *s){

  if (cov->role == ROLE_BROWNRESNICK || cov->role == ROLE_SMITH ||
      cov->role == ROLE_SCHLATHER || cov->role == ROLE_POISSON
      || cov->role == ROLE_POISSON_GAUSS) {
    cov_model *next = cov->sub[DOLLAR_SUB];
    
    cov_model
      *varM = cov->kappasub[DVAR],
      *scaleM = cov->kappasub[DSCALE];
 
    if (varM != NULL && !varM->deterministic) {
      assert(cov->p[DVAR] != NULL);
      VTLG_R(NULL, varM, cov->p[DVAR]);
    }

    if (scaleM != NULL && !scaleM->deterministic) {
      assert(cov->p[DSCALE] != NULL);
      VTLG_R(NULL, scaleM, cov->p[DSCALE]);
    }
    
    DO(next, s);// nicht gatternr
    cov->mpp.maxheight = next->mpp.maxheight * cov->p[DVAR][0];

  } 

  else if (cov->role == ROLE_GAUSS) {    
    assert(cov->key != NULL);
    double 
      *res = cov->key->rf,
      sd = sqrt(cov->p[DVAR][0]);
    int i,
      totalpoints = Gettotalpoints(cov);

    DO(cov->key, s); 

    if (sd != 1.0) 
      for (i=0; i<totalpoints; i++) {
	res[i] *= (res_type) sd;
    }

  }
  
  else ERR("unknown option in doS "); 

 
  if (cov->origrf) {
    assert(cov->prevloc->grid);
    int i, zaehler, d,
      dim = cov->prevloc->timespacedim,
      *cumsum = cov->Sdollar->cumsum,
      *nx = cov->Sdollar->nx,
      *len = cov->Sdollar->len,
      *total = cov->Sdollar->total;
    assert(cov->key != NULL);
    res_type
      *res = cov->rf,
      *rf = cov->key->rf;
    
    assert(nx != NULL && total != NULL && cumsum != NULL);
    for (d=0; d<dim; d++) {
      nx[d] = 0;
    }
    zaehler = 0;
    i = 0;      
    
    while (true) {
      res[i++] = rf[zaehler];
      // printf("%f %d %d", res[i], zaehler,i);
      d = 0;			
      nx[d]++;			
      zaehler += cumsum[d];	
      while (nx[d] >= len[d]) {	
	nx[d] = 0;		
	zaehler -= total[d];
	if (++d >= dim) break;	
	nx[d]++;			
	zaehler += cumsum[d];					
      }
      if (d >= dim) break;			
    }
  }

}


int checkplusmalproc(cov_model *cov) {
  cov_model *sub;
  int i, err, 
    dim = cov->tsdim, 
    xdim = cov->xdimown,
    role = cov->role;
  Types type = ProcessType; 
  domain_type dom = cov->domown;
  isotropy_type iso = cov->isoown;

  //PMI(cov);
  assert(cov->Splus != NULL);

  for (i=0; i<cov->nsub; i++) {

    sub = cov->Splus->keys[i];
     
    if (sub == NULL) 
      SERR("+ or *: named submodels are not given in a sequence.");

    if (!TypeConsistency(type, sub)) return ERRORTYPECONSISTENCY;
    if ((err= CHECK(sub, dim, xdim, type, dom, iso, SUBMODEL_DEP, role))
	== NOERROR) return err;

    //    printf("OK plusprocess\n");
  }

  return NOERROR;
}



int checkplusproc(cov_model *cov) {
  int err;
  if ((err = checkplusmalproc(cov)) != NOERROR) {
    return err;
  }
  return NOERROR;
}



int structplusproc(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  int m, err;
  assert(cov->nr == PLUS_PROC);

  switch(cov->role) {
  case ROLE_GAUSS : 
    {
      location_type *loc = Loc(cov);
      
      if (cov->Splus == NULL) {
	cov->Splus = (plus_storage*) MALLOC(sizeof(plus_storage));
	PLUS_NULL(cov->Splus);
      }
      
      for (m=0; m<cov->nsub; m++) {

	cov_model *sub = cov->sub[m];

	if (cov->Splus->keys[m] != NULL) COV_DELETE(cov->Splus->keys + m);
	if ((err =  covcpy(cov->Splus->keys + m, sub)) != NOERROR) {
	  return err;
	}
	assert(cov->Splus->keys[m] != NULL);
	assert(cov->Splus->keys[m]->calling == cov);
	
	if (PL >= PL_STRUCTURE) {
	  LPRINT("plus: trying initialisation of submodel #%d (%s).\n", m+1, 
		 NICK(sub));
	}
	
	addModel(cov->Splus->keys + m, GAUSSPROC);
	cov_model *fst = cov; while (fst->calling != NULL) fst = fst->calling; 

	//	assert(false);

	err = CHECK(cov->Splus->keys[m], loc->timespacedim, loc->timespacedim,
		    ProcessType, XONLY, NO_ROTAT_INV, cov->vdim, ROLE_GAUSS);
	if (err != NOERROR) {
	  //
	  return err;
	}
	
	//if (m==1) APMI(cov->Splus->keys[m]);
	if ((cov->Splus->struct_err[m] =
	     err = STRUCT(cov->Splus->keys[m], NULL))  > NOERROR) {
	  //	PMI(cov->Splus->keys[m]);
	  //	assert(false);
	  //printf("end plus\n");
	  //	  AERR(err);
	  return err;
	}
	//AERR(err);
	
	//     printf("structplusmal %d %d\n", m, cov->nsub);
	// PMI(cov->Splus->keys[m]);
      }
  
    
      //assert(false);
      return NOERROR;
    }
    
  default :
    SERR2("role '%s' not allowed for '%s'", ROLENAMES[cov->role],
	  NICK(cov));
  }

  return NOERROR;
}

int initplusproc(cov_model *cov, storage VARIABLE_IS_NOT_USED *s){
  int i, err;

  cov->mpp.maxheight = RF_NAN;
  if (cov->Splus == NULL) BUG;

  if (cov->role == ROLE_GAUSS) {
 
    for (i=0; i<cov->nsub; i++) {
      cov_model *sub = cov->Splus == NULL ? cov->sub[i] : cov->Splus->keys[i];
      assert(cov->sub[i]->stor==NULL);
      cov->sub[i]->stor = (storage *) MALLOC(sizeof(storage));
      if (!sub->initialised) {
	if ((err = INIT(sub, 0, cov->sub[i]->stor)) != NOERROR) {
	  return err;
	}
      }
      sub->simu.active = true;
    }
    cov->simu.active = true;

    // assert(false);
  
    cov->fieldreturn = cov->Splus != NULL;
    cov->origrf = false;
    if (cov->Splus != NULL) cov->rf = cov->Splus->keys[0]->rf;
     
    return NOERROR;
  }
     
  else {
    
  }

  return ERRORFAILED;
}


void doplusproc(cov_model *cov, storage VARIABLE_IS_NOT_USED *s) {
  int m, i,
    total = cov->prevloc->totalpoints * cov->vdim;
  double *res = cov->rf;

  if (cov->role == ROLE_GAUSS && cov->method==SpectralTBM) {
    ERR("error in doplus with spectral");
  }
  assert(cov->Splus != NULL);
  
  for (m=0; m<cov->nsub; m++) {
    cov_model *key = cov->Splus->keys[m],
      *sub = cov->sub[m];
    double *keyrf = key->rf;
    DO(key, sub->stor);
    if (m > 0)
      for(i=0; i<total; i++) res[i] += keyrf[i];
  }
  return;
}


