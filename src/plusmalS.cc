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
 Copyright (C) 2005 -- 2015 Martin Schlather

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
#include <R_ext/Lapack.h>

#include "RF.h"
#include "Operator.h"
#include "variogramAndCo.h"
//#include <R_ext/Applic.h>
//#include <R_ext/Utils.h>     
//#include <R_ext/BLAS.h> 







// $

bool hasVarOnly(cov_model *cov) {
   if (cov == NULL || !isDollar(cov)) BUG;
  if (!PisNULL(DSCALE) && P0(DSCALE) != 1.0) return false;
  if (!PisNULL(DANISO) || !PisNULL(DPROJ)) return false;
  int i,
    kappas = CovList[cov->nr].kappas;
  for (i=0; i<kappas; i++)
    if (cov->kappasub[i] != NULL) return false;
  return true;
}


void kappaS(int i, cov_model *cov, int *nr, int *nc){
  switch(i) {
  case DVAR : case DSCALE :
    *nr = *nc = 1;
    break;
  case DANISO :
    *nr = cov->xdimown;
    *nc = SIZE_NOT_DETERMINED;
    break;
  case DAUSER :
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

// simple transformations as variance, scale, anisotropy matrix, etc.  
void Siso(double *x, cov_model *cov, double *v){

  cov_model *next = cov->sub[DOLLAR_SUB];
 int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim;
  double y,
    *aniso=P(DANISO),
    *scale =P(DSCALE),
    var = P0(DVAR);
  assert(cov->Sdollar->simplevar);
  
  y = *x;
  if (aniso != NULL) y = fabs(y * aniso[0]);

  //  printf("scale = %d %d\n", scale == NULL, aniso == NULL);
  //  printf("scale = %f %f \n", scale[0], *x);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] 
      : (y == 0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  // letzteres darf nur passieren wenn dim = 1!!
  COV(&y, next, v);

  for (i=0; i<vdimSq; i++) v[i] *= var; 
}
  

// simple transformations as variance, scale, anisotropy matrix, etc.  
void logSiso(double *x, cov_model *cov, double *v, double *Sign){
  cov_model *next = cov->sub[DOLLAR_SUB];
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim;
  double y, 
    *aniso=P(DANISO),
    *scale =P(DSCALE),
    logvar = log(P0(DVAR));
   assert(cov->Sdollar->simplevar);

  y = *x;
  if (aniso != NULL) y = fabs(y * aniso[0]);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] 
      : (y == 0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  LOGCOV(&y, next, v, Sign);
  for (i=0; i<vdimSq; i++) v[i] += logvar; 
  
 
}
 
void Sstat(double *x, cov_model *cov, double *v){
  logSstat(x, cov, v, NULL);
}

void logSstat(double *x, cov_model *cov, double *v, double *Sign){
  assert(cov->kappasub[DSCALE] == NULL && 
	 (cov->kappasub[DAUSER]==NULL || 
	  CovList[cov->kappasub[DAUSER]->nr].check==checkAngle));
  cov_model *next = cov->sub[DOLLAR_SUB];
    //  *Aniso = cov->kappasub[DAUSER],
    // *Scale = cov->kappasub[DSCALE];
  double *z1 = NULL,
    *scale =P(DSCALE), 
    *aniso=P(DANISO);
  int i,
    nproj = cov->nrow[DPROJ],
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim;

  if (nproj > 0) {
    int *proj = PINT(DPROJ);
    ALLOC_DOLLAR(z, nproj);
    if (scale == NULL || scale[0] > 0.0) {
      if (scale == NULL)  for (i=0; i<nproj; i++) z[i] = x[proj[i] - 1];
      else {
	double invscale = 1.0 / scale[0];
	for (i=0; i<nproj; i++) {
	  z[i] = invscale * x[proj[i] - 1];
	}
      }
    } else {
      // projection and aniso may not be given at the same time
      for (i=0; i<nproj; i++)
	z[i] = (x[proj[i] - 1] == 0 && scale[0] == 0) ? 0.0 : RF_INF;
    } 
    z1 = z;
    //  } else if (Aniso != NULL) {
    //    int dim = Aniso->vdim[0];
    //    ALLOC_DOLLAR(z, dim);
    //    FCTN(x, Aniso, z);
    //    z1 = z;
    //  } else if (Scale != NULL) {
    //    int dim = Aniso->vdim[0];
    //    ALLOC_DOLLAR(z, dim);
    //    FCTN(x, Aniso, z);
    //    z1 = z;
  } else if (aniso==NULL && (scale == NULL || scale[0] == 1.0)) {
    z1 = x;
  } else {
    int xdimown = cov->xdimown;
    double *xz;
    //    printf("xdimown %d\n", xdimown); BUG;
    ALLOC_DOLLAR(z, xdimown);
    if (aniso!=NULL) {
      xA(x, aniso, cov->nrow[DANISO], cov->ncol[DANISO], z);
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
    z1 = z;
  }  


  double var;
  if (cov->Sdollar->simplevar) {
    var = P0(DVAR);
  } else {
    FCTN(x, cov->kappasub[DVAR], &var);
   }

  if (Sign==NULL) {
    COV(z1, next, v);
    for (i=0; i<vdimSq; i++) v[i] *= var; 
  } else {
    LOGCOV(z1, next, v, Sign);
    var = log(var);
    for (i=0; i<vdimSq; i++) v[i] += var; 
  }

}

void Snonstat(double *x, double *y, cov_model *cov, double *v){
  logSnonstat(x, y, cov, v, NULL);
}

void logSnonstat(double *x, double *y, cov_model *cov, double *v, double *Sign){
  cov_model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DAUSER],
    *Scale = cov->kappasub[DSCALE];
  double *z1, *z2, 
    s1 = RF_NA, s2 = RF_NA, smeanSq=RF_NA,
    *scale =P(DSCALE),
    *aniso=P(DANISO);
  int i,
    nproj = cov->nrow[DPROJ],
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim;
  
  if (nproj > 0) {
    int *proj = PINT(DPROJ);
    ALLOC_DOLLARY(Z1, Z2, nproj);
    z1 = Z1;
    z2 = Z2;
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
    int dim = Aniso->vdim[0];
    ALLOC_DOLLARY(Z1, Z2, dim);
    z1 = Z1;
    z2 = Z2;
    FCTN(x, Aniso, z1);
    FCTN(y, Aniso, z2);
  } else if (Scale != NULL && !isRandom(Scale)) {
    int xdimown = cov->xdimown;
    double s;
    ALLOC_DOLLARY(Z1, Z2, xdimown);
    z1 = Z1;
    z2 = Z2;
    FCTN(x, Scale, &s1);
    FCTN(y, Scale, &s2);
    if (s1 <= 0.0 || s2  <= 0.0)
      ERR1("'%s' must be a positive function", KNAME(DSCALE));
    smeanSq = 0.5 * (s1 * s1 + s2 * s2);
    s = sqrt(smeanSq);
    for (i=0; i<xdimown; i++) {
      z1[i] = x[i] / s;
      z2[i] = y[i] / s;
    }
  } else if (aniso==NULL && (scale==NULL || scale[0] == 1.0)) {
    z1 = x;
    z2 = y;
  } else {
    int xdimown = cov->xdimown;
    double *xz1, *xz2;
    ALLOC_DOLLARY(Z1, Z2, xdimown);
    z1 = Z1;
    z2 = Z2;
    if (aniso != NULL) {
      xA( x, y, aniso,cov->nrow[DANISO], cov->ncol[DANISO], z1, z2);
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
  
  double var;
  if (cov->Sdollar->simplevar) {
    var = P0(DVAR);
    if (Sign != NULL) var = log(var);
  } else {
    double w;
    location_type *loc = Loc(cov);
    int dummy = loc->i_row;
    loc->i_row = loc->i_col;
    FCTN(y, cov->kappasub[DVAR], &w);
    loc->i_row = dummy;
    FCTN(x, cov->kappasub[DVAR], &var);
    var *= w;
    var = Sign == NULL ?  sqrt(var)  : 0.5 * log(var);    
  }

  if (Scale != NULL) {
    if (Sign != NULL) var += 0.5 * log(s1 * s2 / smeanSq);
    else var *= sqrt(s1 * s2 / smeanSq);
  }
 
  if (Sign == NULL) {
    NONSTATCOV(z1, z2, next, v);
    //  printf("S-dim %d %f %f %f %f\n", vdimSq, v[0], v[1], v[2], v[3]);
   for (i=0; i<vdimSq; i++) v[i] *= var; 
   //  printf("--> S-dim %d %f %f %f %f\n", vdimSq, v[0], v[1], v[2], v[3]);
  } else {
    LOGNONSTATCOV(z1, z2, next, v, Sign);
    for (i=0; i<vdimSq; i++) {
      v[i] += var; 
      //printf("%f %f \n", v[i], var);
    }
  }
}


void covmatrixS(cov_model *cov, double *v) {
  location_type *loc = Loc(cov);	      
  cov_model *next = cov->sub[DOLLAR_SUB];
  location_type *locnext = Loc(next);
  assert(locnext != NULL);
  int i, err, tot, totSq,
    dim = loc->timespacedim,
    vdim = cov->vdim[0];
  assert(dim == cov->tsdim);
    

  if (cov->Spgs == NULL && 
      (err = alloc_cov(cov, dim, vdim, vdim)) != NOERROR)
    ERR("memory allocation error in 'covmatrixS'");
  
  if ((!PisNULL(DSCALE) && P0(DSCALE) != 1.0) || 
      !PisNULL(DANISO) || !PisNULL(DPROJ) || 
      cov->kappasub[DSCALE] != NULL ||
      cov->kappasub[DAUSER] != NULL ||
      cov->kappasub[DPROJ] != NULL
      ) {
    CovarianceMatrix(cov, v); 
    return;
  }

  if (next->xdimprev != next->xdimown) {
    BUG; // fuehrt zum richtigen Resultat, sollte aber nicht
    // vorkommen!
    CovarianceMatrix(cov, v); 
    return;
  }

  int next_gatter = next->gatternr,
    next_xdimgatter = next->xdimgatter,
    next_xdim = next->xdimprev;
 
  next->gatternr = cov->gatternr;
  next->xdimprev = cov->xdimprev;
  next->xdimgatter = cov->xdimgatter;
  CovList[next->nr].covmatrix(next, v);//hier wird uU next->totalpoints gesetzt
  next->gatternr = next_gatter;
  next->xdimgatter = next_xdimgatter;
  next->xdimprev = next_xdim;

  tot = cov->vdim[0] * locnext->totalpoints;
  totSq = tot * tot;
  
  if (cov->Sdollar->simplevar) {
    double var = P0(DVAR);
    if (var == 1.0) return;
    for (i=0; i<totSq; v[i++] *= var);
  } else {
    BUG;
   }
}

char iscovmatrixS(cov_model *cov) {
  cov_model *sub = cov->sub[DOLLAR_SUB];
  return (int) ((PisNULL(DSCALE) || P0(DSCALE) ==1.0) &&	
		PARAMisNULL(cov, DAUSER) &&
		PARAMisNULL(cov, DPROJ) &&
		cov->Sdollar->simplevar && // to do
		PARAMisNULL(cov, DANISO)) * CovList[sub->nr].is_covariance(sub);
}

void DS(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[DOLLAR_SUB];
  assert( cov->kappasub[DAUSER] == NULL && cov->kappasub[DSCALE] == NULL);
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim,
    nproj = cov->nrow[DPROJ];
  double y[2], varSc,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;
  assert(cov->Sdollar->simplevar);
  assert(isCartesian(cov->isoown));

  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varSc = P0(DVAR) * spinvscale;

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
  assert(cov->kappasub[DAUSER] == NULL);
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim,
    nproj = cov->nrow[DPROJ],
    *proj = PINT(DPROJ);
  double y[2], varScSq,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;
  
  assert(isCartesian(cov->isoown));
  assert(cov->Sdollar->simplevar);
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScSq = P0(DVAR) * spinvscale * spinvscale;
  
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
  assert(cov->kappasub[DAUSER] == NULL);
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim,
    nproj = cov->nrow[DPROJ],
    *proj = PINT(DPROJ);
  double y[2], varScS3,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;

  assert(isCartesian(cov->isoown));
  assert(cov->Sdollar->simplevar);
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScS3 = P0(DVAR) * spinvscale * spinvscale * spinvscale;
  
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

void D4S(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DAUSER] == NULL);
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim,
    nproj = cov->nrow[DPROJ],
    *proj = PINT(DPROJ);
  double y[2], varScS4,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;

  assert(isCartesian(cov->isoown));
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScS4 = spinvscale * spinvscale;
  varScS4 *= varScS4 * P0(DVAR);
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
  for (i=0; i<vdimSq; i++) v[i] *= varScS4; 
}


void nablahessS(double *x, cov_model *cov, double *v, bool nabla){
  cov_model *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DAUSER];
  if (Aniso != NULL) BUG;
  int i, endfor,
    dim = cov->nrow[DANISO],// == ncol == x d i m ?
    xdimown = cov->xdimown,
    nproj = cov->nrow[DPROJ];
  double *xy, *vw,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    var = P0(DVAR);
  if (nproj != 0) BUG;
  if (dim != xdimown) BUG;
	   
  if (!cov->Sdollar->simplevar) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");

  
  if (aniso != NULL) {  
    ALLOC_DOLLAR(z, xdimown);
    ALLOC_DOLLAR2(w, xdimown);
    xA(x, aniso, xdimown, xdimown, z);
    xy = z;
    vw = w;
  } else {
    xy = x;
    vw = v;
  }

  if (scale != NULL) {
    ALLOC_DOLLAR3(y, xdimown);
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
    else {
      XCXt(aniso, vw, v, xdimown, xdimown); //todo:?reprogramm XCXt with alloc here ?
    }
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
  int i,
    idx[2] = {DAUSER, DPROJ};
  double 
    scale;
  
 if (cov->kappasub[DVAR] != NULL) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");

  for (i=0; i<2; i++) {
    if (cov->kappasub[idx[i]] != NULL) {
      char Msg[100];
      sprintf(Msg, "inverse can only be calculated if '%s' is not an arbitrary function",
  	      KNAME(idx[i])); 
      ERR(Msg);
    }
  }
  if (cov->kappasub[DSCALE] != NULL) {
    double left;
    NONSTATINVERSE(ZERO, cov->kappasub[DSCALE], &left, &scale);
    if (left < 0.0) ERR("scale not defined to be non-negative.");
  } else scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
 
  int
    dim = cov->xdimown,
    nproj = cov->nrow[DPROJ];
  //    *proj = (int *)P(DPROJ];
  double y, 
    s = 1.0,
    *aniso=P(DANISO),
    var = P0(DVAR);

  //PMI(cov); crash();
  if (dim != 1) BUG;

  if (aniso != NULL) {
    if (isMiso(Type(aniso, cov->nrow[DANISO], cov->ncol[DANISO]))) 
      s /= aniso[0];
    else NotProgrammedYet(""); // to do
  }
  s *= scale;  
  if (nproj == 0) {
    y= *x / var; // inversion, so variance becomes scale
  } else {
    BUG;  //ERR("nproj is not allowed in invS");
  }
  
  if (CovList[next->nr].inverse == ErrCov) BUG;
  INVERSE(&y, next, v);
 
  for (i=0; i<dim; i++) v[i] *= s; //!
}


void nonstatinverseS(double *x, cov_model *cov, double *left, double*right,
		     bool log){
  cov_model
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DAUSER],
    *Scale = cov->kappasub[DSCALE];
   int i,
    dim = cov->tsdim,
    nproj = cov->nrow[DPROJ];
  //    *proj = (int *)P(DPROJ];
  double y, 
    s = 1.0,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    var = P0(DVAR);
 
  if (cov->kappasub[DVAR] != NULL) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");


  if (nproj == 0) {
    y= *x / var; // inversion, so variance becomes scale
  } else {
    BUG;  //ERR("nproj is not allowed in invS");
  }
    
  if (CovList[next->nr].nonstat_inverse == ErrInverseNonstat) BUG;
  if (log) {
    NONSTATLOGINVERSE(&y, next, left, right);
  } else {
    NONSTATINVERSE(&y, next, left, right);
  }

 
  if (aniso != NULL) {
    if (isMiso(Type(aniso, cov->nrow[DANISO], cov->ncol[DANISO]))) s/=aniso[0];
    else {
      dollar_storage *S = cov->Sdollar;
      int  
	ncol = cov->ncol[DANISO],
	nrow = cov->nrow[DANISO],
	ncnr = ncol * nrow,
	bytes = ncnr * sizeof(double),
	size = ncol * sizeof(double);
      bool redo;
      if (ncol != nrow) BUG;
      if ((redo = S->save_aniso == NULL)) {
	S->save_aniso = (double *) MALLOC(bytes);
	S->inv_aniso = (double *) MALLOC(bytes);
      }
      ALLOC_DOLLAR4(LR, ncol);
      double *save = S->save_aniso,
	*inv = S->inv_aniso;
      if (!redo) {
	for (i=0; i<ncnr; i++) if ((redo = save[i] != P(DANISO)[i])) break;
      }
      if (redo) {
	MEMCOPY(save, P(DANISO), bytes);
	MEMCOPY(inv, P(DANISO), bytes);
	if (Ext_invertMatrix(inv, nrow) != NOERROR)
	  ERR("inversion of anisotropy matrix failed");
      }
      
      MEMCOPY(LR, right, size);
      xA(LR, inv, nrow, ncol, right);
      
      MEMCOPY(LR, left, size);
      xA(LR, inv, nrow, ncol, left);      
    }    
  }

  if (Aniso != NULL) {
    if (aniso != NULL) BUG;

    if (CovList[Aniso->nr].inverse == ErrInverse) 
      ERR("inverse of anisotropy matrix function unknown");
    int 
      nrow = Aniso->vdim[0],
      ncol = Aniso->vdim[1],
      size = nrow * sizeof(double);
    if (cov->xdimown != ncol || ncol != 1)
      ERR("anisotropy function not of appropriate form");
    ALLOC_DOLLAR4(LR, nrow);
    
    MEMCOPY(LR, right, size);
    INVERSE(LR, Aniso, right);
          
    MEMCOPY(LR, left, size);
    INVERSE(LR, Aniso, left);
  }

  if (Scale != NULL && !isRandom(Scale)) {
    double dummy;
    COV(ZERO, Scale, &dummy);
    s *= dummy;
  } else if (scale != NULL) s *= scale[0];  
  if (s != 1.0) {
    for (i=0; i<dim; i++) {
      left[i] *= s; //!
      right[i] *= s;
    }
  }

}

void nonstatinverseS(double *x, cov_model *cov, double *left, double*right) {
  nonstatinverseS(x, cov, left, right, false);
}

void nonstat_loginverseS(double *x, cov_model *cov, double *left, double*right){
 nonstatinverseS(x, cov, left, right, true);
}

void coinitS(cov_model *cov, localinfotype *li) {
  assert(cov->Sdollar->simplevar);
  cov_model *next = cov->sub[DOLLAR_SUB];
  if ( CovList[next->nr].coinit == NULL)
    ERR("# cannot find coinit -- please inform author");
  CovList[next->nr].coinit(next, li);
}
void ieinitS(cov_model *cov, localinfotype *li) {
   assert(cov->Sdollar->simplevar);
 cov_model *next = cov->sub[DOLLAR_SUB];
  
  if ( CovList[next->nr].ieinit == NULL)
    ERR("# cannot find ieinit -- please inform author");
  CovList[next->nr].ieinit(next, li);
}

void tbm2S(double *x, cov_model *cov, double *v){
   assert(cov->Sdollar->simplevar);
 cov_model *next = cov->sub[DOLLAR_SUB];
  cov_fct *C = CovList + next->nr; // kein gatternr, da isotrop
  double y[2],  *xy,
    *scale =P(DSCALE),
    *aniso = P(DANISO);
  assert(cov->kappasub[DAUSER] == NULL);
  
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
	if (P0(DANISO) == 0.0) {
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
  *v *= P0(DVAR);
}


// TODO : Aniso=Matrix: direkte implementierung in S,
// sodass nicht ueber initS gegangen werden muss, sondern
//  e  < -  e^\top Aniso


int TaylorS(cov_model *cov) {
  cov_model 
    *next = cov->sub[DOLLAR_SUB],
    *sub = cov->key == NULL ? next : cov->key;
  int i;

  if (PisNULL(DPROJ) && PisNULL(DANISO)) {
    double scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
    cov->taylorN = sub->taylorN;  
    for (i=0; i<cov->taylorN; i++) {
      cov->taylor[i][TaylorPow] = sub->taylor[i][TaylorPow];
      cov->taylor[i][TaylorConst] = sub->taylor[i][TaylorConst] *
	P0(DVAR) * pow(scale, -sub->taylor[i][TaylorPow]);   
    }

    cov->tailN = sub->tailN;  
    for (i=0; i<cov->tailN; i++) {
      cov->tail[i][TaylorPow] = sub->tail[i][TaylorPow];
      cov->tail[i][TaylorExpPow] = sub->tail[i][TaylorExpPow];
      cov->tail[i][TaylorConst] = sub->tail[i][TaylorConst] *
	P0(DVAR) * pow(scale, -sub->tail[i][TaylorPow]);   
      cov->tail[i][TaylorExpConst] = sub->tail[i][TaylorExpConst] *
	pow(scale, -sub->tail[i][TaylorExpPow]);
    }
  } else {
    cov->taylorN = cov->tailN = 0;
  }
  return NOERROR;
}

int checkS(cov_model *cov) {

  // hier kommt unerwartet  ein scale == nan rein ?!!
  cov_model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DAUSER],
    *Scale = cov->kappasub[DSCALE],
    *sub = cov->key == NULL ? next : cov->key;
   isotropy_type isoown = cov->isoown;
   int i, err,
    xdimown = cov->xdimown,
     xdimNeu = xdimown,
     // KOMMENTAR NICHT LOESCHEN
     // *proj = PINT(DPROJ), // auf keinen Fall setzen, da Pointer unten neu
   //                        gesetzt wird!!!!
     nproj = cov->nrow[DPROJ];  
  // bool skipchecks = GLOBAL.general.skipchecks;
  matrix_type type = TypeMany;
 
  assert(isAnyDollar(cov));
  if (!isDollarProc(cov)) cov->nr = DOLLAR; // wegen nr++ unten !
  
  // cov->q[1] not NULL then A has been given
  
  // if (cov->method == SpectralTBM && cov->xdimown != next->xdimprev)
  //  return ERRORDIM;
  
  if ((err = checkkappas(cov, false)) != NOERROR) {
    return err;
  }
  kdefault(cov, DVAR, 1.0);

  bool angle = Aniso!=NULL && CovList[Aniso->nr].check == checkAngle;
  if (angle && PisNULL(DANISO) && PisNULL(DAUSER) &&
      cov->xdimown == cov->tsdim) {
    int dim  = cov->tsdim;
    ASSERT_CARTESIAN;
    if (!isCartesian(cov->isoown))  return ERRORANISO;
    if ((err = CHECK(Aniso, dim, dim, ShapeType, XONLY,
		    CARTESIAN_COORD, SUBMODEL_DEP, cov->role)) == NOERROR) {
      if (verysimple(Aniso)) {     
	 PALLOC(DAUSER, dim, dim);
 	 AngleMatrix(Aniso, P(DAUSER));
	 COV_DELETE(cov->kappasub + DAUSER);
	 Aniso = cov->kappasub[DAUSER] = NULL;
	 angle = false;
      }
    }
  }


  if (cov->kappasub[DANISO] != NULL)
    SERR1("'%s' may not be a function.", KNAME(DANISO));
  if (!PisNULL(DAUSER)) {
    if (!isCartesian(cov->isoown))  return ERRORANISO;
  
    //    if (GLOBAL.internal.warn_Aniso) {
    //      PRINTF("NOTE! Starting with RandomFields 3.0, the use of '%s' is different from\nthe former '%s' insofar that '%s' is multiplied from the right by 'x' (i.e. Ax),\nwhereas '%s' had been multiplied from the left by 'x' (i.e. xA).\n", KNAME(DAUSER), KNAME(DANISO), KNAME(DANISO), KNAME(DAUSER));
    //    }
    GLOBAL.internal.warn_Aniso = false;
    // here for the first time
    if (!PisNULL(DANISO)) return ERRORANISO_T; 
    int k,
      lnrow = cov->nrow[DAUSER],
      lncol = cov->ncol[DAUSER];
    long j, 
      total = lncol * lnrow;
	
    double
      *pA = P(DAUSER); 
    PALLOC(DANISO, lncol, lnrow); // !! ACHTUNG col, row gekreuzt
    for (i=k=0; i<lnrow; i++) {
      for (j=i; j<total; j+=lnrow) P(DANISO)[k++] = pA[j];
    }
    PFREE(DAUSER);
  }

  bool simplevar = cov->kappasub[DVAR] == NULL || isRandom(cov->kappasub[DVAR]);
  if (!simplevar) {

    //cov->domown = KERNEL;
    ptwise_type ptt = cov->ptwise_definite;

    isotropy_type isonew = UpgradeToCoordinateSystem(cov->isoown);
    if (isonew == ISO_MISMATCH) SERR("reduced systems not allowed");
    if ((err = CHECK(cov->kappasub[DVAR], cov->tsdim, cov->xdimown, 
		     ShapeType, // only!! -- for pos def use RMprod
		     XONLY, isonew,
		     SCALAR, ROLE_BASE)) != NOERROR) return err;        
    if (cov->kappasub[DVAR]->ptwise_definite != pt_posdef) {
      if (cov->kappasub[DVAR]->ptwise_definite == pt_unknown) {
	if (GLOBAL.internal.warn_negvar && cov->q==NULL) {
	  QALLOC(1);
	  warning("positivity of the variance in '%s' cannot be detected.",
		  NICK(next));
	}
      } else {
	SERR2("positivity of '%s' required. Got '%s'", KNAME(DVAR), 
	      POSITIVITY_NAMES[cov->kappasub[DVAR]->ptwise_definite]);
      }
    }
    cov->ptwise_definite = ptt;
  }

  DOLLAR_STORAGE;
  if (Aniso != NULL) {
    if (!isDollarProc(cov) && !angle && cov->domown != KERNEL) 
      return ERRORFAILED;
    if (!PisNULL(DANISO) || !PisNULL(DPROJ) || !PisNULL(DSCALE) ||
	(Scale!=NULL && !isRandom(Scale)) )
      SERR2("if '%s' is an arbitrary function, only '%s' may be given as additional parameter.", KNAME(DAUSER), KNAME(DVAR));
    if (cov->isoown != SYMMETRIC && !isCoordinateSystem(cov->isoown)) {
      return ERRORANISO;
    }
 
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = true;
    isotropy_type isonew = UpgradeToCoordinateSystem(cov->isoown);
    if (isonew == ISO_MISMATCH) SERR("reduced systems not allowed");
    if ((err = CHECK(Aniso, cov->tsdim, cov->xdimown, ShapeType, XONLY,
		     isonew, SUBMODEL_DEP, cov->role)) != NOERROR) {
      return err;
    }

    if (Aniso->vdim[1] != 1)
      SERR4("'%s' returns a matrix of dimension %d x %d, but dimension %d x 1 is need.",
	    KNAME(DANISO), Aniso->vdim[0], Aniso->vdim[1], cov->xdimown);

    if (cov->key==NULL) {
      if ((err = CHECK(sub, Aniso->vdim[0], Aniso->vdim[0], cov->typus,  
		       cov->domown, 
		       cov->isoown, SUBMODEL_DEP, cov->role)) != NOERROR) {
	return err;
      }
    }
    cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = 
      cov->pref[Hyperplane] = cov->pref[SpectralTBM] = cov->pref[TBM] = 
      PREF_NONE;

    sub->pref[CircEmbed] = sub->pref[CircEmbedCutoff] = 
      sub->pref[CircEmbedIntrinsic] = sub->pref[Sequential] = 
      sub->pref[Specific] = PREF_NONE;
 
  } else if (Scale != NULL && !isRandom(Scale)) {
    if (!isDollarProc(cov) && cov->domown != KERNEL) return ERRORFAILED;
    if (!PisNULL(DANISO) || !PisNULL(DPROJ) || !PisNULL(DSCALE))
      SERR2("if '%s' is an arbitrary function, only '%s' may be given as additional parameter.", KNAME(DSCALE), KNAME(DVAR));
 
   if (cov->isoown != SYMMETRIC && !isCoordinateSystem(cov->isoown)) {
      return ERRORANISO;
    }
  
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = true;
    isotropy_type isonew = UpgradeToCoordinateSystem(cov->isoown);
    if (isonew == ISO_MISMATCH) SERR("reduced systems not allowed");
    if ((err = CHECK(Scale, cov->tsdim, cov->xdimown, ShapeType, XONLY,
		     isonew, SUBMODEL_DEP, cov->role)) != NOERROR) {
      return err;
    }

    if (Scale->vdim[1] != 1 || Scale->vdim[0] != 1)
      SERR3("'%s' must be scalar, not %d x %d",
	    KNAME(DSCALE), Aniso->vdim[0], Aniso->vdim[1]);

    if (cov->key==NULL) {
      if ((err = CHECK(sub, Scale->vdim[0], Scale->vdim[0], cov->typus,  
		       cov->domown, 
		       cov->isoown, SUBMODEL_DEP, cov->role)) != NOERROR) {
	return err;
      }
      if (!isNormalMixture(next)) 
	SERR("scale function only allowed for normal mixtures.");
    }

    cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = 
      cov->pref[Hyperplane] = cov->pref[SpectralTBM] = cov->pref[TBM] = 
      PREF_NONE;

    sub->pref[CircEmbed] = sub->pref[CircEmbedCutoff] = 
      sub->pref[CircEmbedIntrinsic] = sub->pref[Sequential] = 
      sub->pref[Specific] = PREF_NONE;

    } else if (!PisNULL(DANISO)) { // aniso given
    int 
      nrow = cov->nrow[DANISO],
      ncol = cov->ncol[DANISO];

   
    if (nrow==0 || ncol==0) SERR("dimension of the matrix is 0");
    if (!PisNULL(DPROJ)) return ERRORANISO_T;
    if (xdimown < nrow) {
      if (PL >= PL_ERRORS) {LPRINT("xdim=%d != nrow=%d\n", xdimown, nrow);}
      SERR("#rows of anisotropy matrix does not match dim. of coordinates");
    }
    if (xdimown != cov->tsdim && nrow != ncol)
      SERR("non-quadratic anisotropy matrices only allowed if dimension of coordinates equals spatio-temporal dimension");
 
    type = Type(P(DANISO), nrow, ncol);
    
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = true;


    switch (cov->isoown) {
    case ISOTROPIC :   
      if (cov->tsdim != 1) return ERRORANISO;
      cov->full_derivs = cov->rese_derivs = 2;
      break;
    case SPACEISOTROPIC :  
      cov->full_derivs =  cov->rese_derivs = isMdiag(type) ? 2 : 0;
      if (nrow != 2 || !isMdiag(type)) {
	SERR("spaceisotropy needs a 2x2 diagonal matrix");
      }
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
    case CARTESIAN_COORD : case GNOMONIC_PROJ :  case ORTHOGRAPHIC_PROJ:
      if (!isProcess(cov->typus)) return ERRORANISO;
      break;
    default : 
      if (isCartesian(cov->isoown)) { BUG; }
      else return  ERRORANISO;
    }

    if (!isMdiag(type)) 
      cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = cov->pref[Hyperplane] = PREF_NONE;
    if (cov->isoown != SPACEISOTROPIC && !isMiso(type))
      cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;

    if (cov->key==NULL) {
      if ((err = CHECK(next, ncol, ncol, cov->typus, cov->domown, 
		 ncol==1 && !isProcess(cov->typus) ? ISOTROPIC : cov->isoown, 
		 SUBMODEL_DEP, cov->role))
	  != NOERROR) {
 	return err;
      }
      if (next->domown == cov->domown && next->isoown == cov->isoown &&
	  xdimown > 1) next->delflag = DEL_COV - 7;
    } else {
      if ((err = CHECK(cov->key, ncol, ncol, cov->typus, cov->domown, 
		       ncol==1 && 
		       !isProcess(cov->typus) ? ISOTROPIC : cov->isoown,
		       SUBMODEL_DEP, cov->role)) != NOERROR) return err;
    }
  
  } else { // PisNULL(DANISO) 
    
    int tsdim = cov->tsdim;
    if (nproj > 0) {
      if (cov->Sdollar->proj < 0) { // different isoown cause different interpretations; here: restore original value
	PFREE(DPROJ);
	kdefault(cov, DPROJ, cov->Sdollar->proj);
	nproj = 1;
      } 
      if (P0INT(DPROJ) < 0) { // here: interprete "space" and "time"
	//printf("%d\n", nproj);
	if (nproj != 1) {
	  BUG;
	  SERR1("unallowed use of '%s'", KNAME(DPROJ));
	}
	cov->Sdollar->proj = P0INT(DPROJ);
	if (Loc(cov)->Time && tsdim >= 2) {
	  assert(P0INT(DPROJ) == PROJ_TIME || P0INT(DPROJ) == PROJ_SPACE);
	  bool Time = P0INT(DPROJ) == PROJ_TIME;	  
	  PFREE(DPROJ);
	  if (Time) {
	    kdefault(cov, DPROJ, tsdim);
	  } else {
	    nproj = tsdim - 1;
	    PALLOC(DPROJ, nproj, 1);
	    for (i=1; i<=nproj; i++) PINT(DPROJ)[i-1] = i;
	  }
	} else {
	  SERR1("unallowed use of '%s' or model to complicated.", KNAME(DPROJ));
	}
	// PMI(cov, 0);
      }

      if (nproj < tsdim) {	
	for (i = 0; i<Nothing; i++) { if (cov->pref[i] > 0) cov->pref[i] = 1; }
	cov->pref[Specific] = PREF_BEST;
      }


      bool ok = (cov->xdimprev == cov->xdimown) ||
	((cov->calling == NULL || cov->calling->isoown == EARTH_COORDS ||
	  cov->calling->isoown == EARTH_SYMMETRIC) &&
	 isCartesian(cov->isoown) && cov->xdimprev == cov->xdimown - 1);
      if (!ok) {	
	//	printf("!ok %d %d   %d %d\n", cov->xdimprev, cov->xdimown, cov->isoprev == EARTH_COORDS, isCartesian(cov->isoown) );
	//PMI(cov->calling);
	return ERRORANISO;
      }
      for (i=0; i<nproj; i++) {
	int j,
	  idx = PINT(DPROJ)[i];
	if (idx <= 0) 
	  SERR1("only positive values allowed for '%s'", KNAME(DPROJ));
	if (idx > xdimown)
	  SERR3("%d-th value of '%s' (%d) out of range", 
		i, KNAME(DPROJ), PINT(DPROJ)[i]);
	for (j=i+1; j<nproj; j++) {
	  if (PINT(DPROJ)[j] == idx) { BUG;
	    SERR1("values of '%s' must be distinct.", KNAME(DPROJ));
	  }
	}
      }
      xdimNeu = nproj;
      tsdim = nproj;

      switch (cov->isoown) {
      case ISOTROPIC : case EARTH_ISOTROPIC : case SPHERICAL_ISOTROPIC : 
	if (cov->tsdim != 1) return ERRORANISO;
	break;
      case SPACEISOTROPIC : 
	if (nproj > 2) SERR("maximum length of projection vector is 2");
	if (nproj == 2) {
	  if (PINT(DPROJ)[0] >= PINT(DPROJ)[1])
	    SERR1("in case of '%s' projection directions must be ordered",
		  ISONAMES[SPACEISOTROPIC]);
	}
	if (P0INT(DPROJ) == 1) {
	  tsdim = cov->tsdim - 1 ;
	  isoown = ISOTROPIC;
	  xdimNeu = 1;
	} else {
	  assert(P0INT(DPROJ) == 2);
	  tsdim = 1;
	  isoown = ISOTROPIC;
	  xdimNeu = 1;
	}
	break;      
      case ZEROSPACEISO :
	isoown = SYMMETRIC;
	break;      
      case VECTORISOTROPIC : SERR("projection of vectorisotropic fields not programmed yet"); // to do: vdim muss auch reduziert werden ... --- das wird
	// grausam !
	break;
      case SYMMETRIC: case CARTESIAN_COORD:
	break;
      case GNOMONIC_PROJ :  case ORTHOGRAPHIC_PROJ :
	isoown = CARTESIAN_COORD;
	break;
      case PREVMODELI : BUG;      
	break;
      case SPHERICAL_SYMMETRIC : case EARTH_SYMMETRIC :
	if (nproj != 2 || PINT(DPROJ)[0] != 1 || PINT(DPROJ)[1]  != 2)
	  isoown = SYMMETRIC;
	break;
      case SPHERICAL_COORDS : case EARTH_COORDS :
	if (nproj != 2 || PINT(DPROJ)[0] != 1 || PINT(DPROJ)[1]  != 2)
	  isoown = CARTESIAN_COORD;
	break;
      default : 
	if (isCartesian(cov->isoown)) {BUG;}
	else return  ERRORANISO;  // todo
      }
    }
  
    // verhindern, dass die coordinaten-transformation anlaeuft,
    // da aus z.B. Erd-coordinaten durch Projektion (auf die Zeitachse)
    // kartesische werden
     if (cov->key==NULL) {
      if ((err = CHECK_NO_TRAFO(next, tsdim, xdimNeu, cov->typus, cov->domown,
				isoown, 
				cov->vdim[0], // SUBMODEL_DEP; geaendert 20.7.14
				cov->role)) != NOERROR) {
	return err;
      }

      if (next->domown == cov->domown &&
	  next->isoown == cov->isoown) // darf nicht nach C-HECK als allgemeine Regel ruebergezogen werden, u.a. nicht wegen stat/nicht-stat wechsel !!
	// hier geht es, da domown und isoown nur durchgegeben werden und die Werte      // bereits ein Schritt weiter oben richtig/minimal gesetzt werden.
	next->delflag = DEL_COV - 8;
    } else {
      if ((err = CHECK_NO_TRAFO(cov->key, tsdim, xdimNeu, cov->typus,
		       cov->domown, isoown,
		       SUBMODEL_DEP, cov->role)) != NOERROR) return err;
    }

     //     PMI(cov);
  } // end no aniso
 
  if (( err = checkkappas(cov, false)) != NOERROR) {
    return err;
  }
  
  setbackward(cov, sub);

  if ((Aniso != NULL || (Scale != NULL && !isRandom(Scale)) || 
       !PisNULL(DANISO)|| !PisNULL(DPROJ)) && cov->maxdim < cov->xdimown) 
    cov->maxdim = cov->xdimown;
 	
  if (!isAnyIsotropic(cov->isoown) && !isDollarProc(cov)) { // multivariate kann auch xdimNeu == 1 problematisch sein
    cov->nr++;
  }
  
  if (xdimNeu > 1) {
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 0;
  }

  // 30.10.11 kommentiert:
  //  cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 
  //      cov->pref[TBM] = cov->pref[SpectralTBM] = 0;
  if ( (PisNULL(DANISO) || isMiso(type)) && PisNULL(DPROJ)) {
    cov->logspeed = sub->logspeed * P0(DVAR);
  }
  //////////////// 

  if (sub->pref[Nugget] > PREF_NONE) cov->pref[Specific] = 100;
  if (nproj == 0) cov->matrix_indep_of_x = sub->matrix_indep_of_x;

  if ((err = TaylorS(cov)) != NOERROR) return err;  
  cov->Sdollar->simplevar = simplevar;
  cov->Sdollar->isoown = isoown; // for struct Sproc !
  
  if (isProcess(cov->typus)) {
    MEMCOPY(cov->pref, PREF_NOTHING, sizeof(pref_shorttype)); 
  }
 
  if (GLOBAL.coords.coord_system == earth && 
      is_all(isCartesian, CovList + next->nr) &&
      GLOBAL.internal.warn_scale &&
      (PisNULL(DSCALE) || 
       P0(DSCALE) < (strcmp(GLOBAL.coords.newunits[0], "km")== 0 ? 10 : 6.3))) {
    char msg[300];

    sprintf(msg, "value of scale parameter equals '%4.2f',\nwhich is less than 100, although models defined on R^3 are used in the\ncontext of earth coordinates where larger scales are expected.\n(This warning appears only ones per session.)\n", PisNULL(DSCALE) ? 1.0 : P0(DSCALE)); 
    GLOBAL.internal.warn_scale = false;
    warning(msg);
  }
   return NOERROR;
}




void rangeS(cov_model *cov, range_type* range){
  int i;
  bool negdef = isNegDef(cov->typus);
  range->min[DVAR] = negdef ? 0.0 : RF_NEGINF;
  range->max[DVAR] = RF_INF;
  range->pmin[DVAR] = negdef ? 0.0 : -10000;
  range->pmax[DVAR] = 100000;
  range->openmin[DVAR] = !negdef;
  range->openmax[DVAR] = true;

  range->min[DSCALE] = 0.0;
  range->max[DSCALE] = RF_INF;
  range->pmin[DSCALE] = 0.0001;
  range->pmax[DSCALE] = 10000;
  range->openmin[DSCALE] = true;
  range->openmax[DSCALE] = true;

  for (i=DANISO; i<= DAUSER; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -10000;
    range->pmax[i] = 10000;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
  
  range->min[DPROJ] = -2;
  range->max[DPROJ] = cov->tsdim;
  range->pmin[DPROJ] = 1;
  range->pmax[DPROJ] =  cov->tsdim;
  range->openmin[DPROJ] = false;
  range->openmax[DPROJ] = false;
}


bool TypeS(Types required, cov_model *cov, int depth) {
  cov_model *sub = cov->key==NULL ? cov->sub[0] : cov->key;

  // print("  TypeS %s %s :  %d %d\n", TYPENAMES[required], NAME(sub), isShape(required) || isTrend(required) || isProcess(required),  TypeConsistency(required, sub, depth-1));

  return (isShape(required) || isTrend(required) || isProcess(required))
    && TypeConsistency(required, sub, depth-1);
}


void spectralS(cov_model *cov, gen_storage *s, double *e) {
  cov_model *next = cov->sub[DOLLAR_SUB];
  int d,
    ncol = PisNULL(DANISO) ? cov->tsdim : cov->ncol[DANISO];
  double sube[MAXTBMSPDIM],
    *scale =P(DSCALE),
    invscale = 1.0;
  
  SPECTRAL(next, s, sube); // nicht gatternr

  // Reihenfolge nachfolgend extrem wichtig, da invscale auch bei aniso
  // verwendet wird


  if (scale != NULL) invscale /= scale[0];
  
  if (!PisNULL(DANISO)) {
    int 
      nrow = cov->nrow[DANISO];
    long j, k, m,
      total = ncol * nrow;
    double
      *A = P(DANISO); 
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


void ScaleDollarToLoc(cov_model *to, cov_model *from,
		      int VARIABLE_IS_NOT_USED depth) {

  assert(!PARAMisNULL(to, LOC_SCALE));
  assert(isDollar(from));
  assert(!PARAMisNULL(from, DSCALE));
  PARAM(to, LOC_SCALE)[0] = PARAM0(from, DSCALE);
}

bool ScaleOnly(cov_model *cov) {
  return isDollar(cov) && 
    PisNULL(DPROJ) && cov->kappasub[DPROJ] == NULL &&
    PisNULL(DANISO) &&  cov->kappasub[DAUSER] == NULL &&
    (PisNULL(DVAR) || P0(DVAR)==1.0) && cov->kappasub[DVAR] == NULL;
}



int addScales(cov_model **newmodel, double anisoScale, cov_model *scale,
	      double Scale) {
  //  cov_model *calling = (*newmodel)->calling;

  if (anisoScale != 1.0) {
    addModel(newmodel, LOC);
    kdefault(*newmodel, LOC_SCALE, anisoScale);
  }

  if (scale != NULL) {
    if (!isRandom(scale)) SERR("unstationary scale not allowed yet");
    addModel(newmodel, LOC);
    addSetDistr(newmodel, scale->calling, ScaleDollarToLoc, true, MAXINT);
  } else {
    if (Scale != 1.0) {
      addModel(newmodel, LOC);
      kdefault(*newmodel, LOC_SCALE, Scale);
    }
  }  
  return NOERROR;
}

int structS(cov_model *cov, cov_model **newmodel) {
  if (cov->role == ROLE_GAUSS && isProcess(cov->typus)) {
    cov->nr = DOLLAR_PROC;
    return structSproc(cov, newmodel); // kein S-TRUCT(...) !!
  } 

  cov_model *local = NULL,
    *dummy = NULL,
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DAUSER],
    *Scale = cov->kappasub[DSCALE];
  int err = NOERROR; 
  bool generalAniso = (Aniso != NULL && CovList[Aniso->nr].check != checkAngle) 
    || (Scale != NULL && !isRandom(Scale));

   
  ASSERT_NEWMODEL_NOT_NULL;
 
  if (cov->kappasub[DVAR] != NULL) {
    GERR2("Arbitrary functions for '%s' should be replaced by multiplicative models using '%s'", KNAME(DVAR), CovList[PROD].nick);
  } 

  if (generalAniso)
    GERR1("complicated models including arbitrary functions for '%s' cannot be simulated yet", KNAME(DAUSER));

  switch (cov->role) {
  case ROLE_POISSON_GAUSS : 
  case ROLE_SMITH : 
    if (!next->deterministic) GERR("random shapes not programmed yet");
    if (!PisNULL(DANISO)) GERR("anisotropy parameter not allowed yet");
    if (!PisNULL(DPROJ)) GERR("projections not allowed yet");
    if ((err = STRUCT(next, newmodel)) > NOERROR) return err;

    addModel(newmodel, DOLLAR);
    if (!PisNULL(DVAR)) kdefault(*newmodel, DVAR, P0(DVAR));
    if (!PisNULL(DSCALE)) kdefault(*newmodel, DSCALE, P0(DSCALE));
    break;

  case ROLE_MAXSTABLE : // eigentlich nur von RPSmith moeglich !
    if (!next->deterministic) GERR("random shapes not programmed yet");
    if (!PisNULL(DPROJ)) GERR("projections not allowed yet");
    // P(DVAR) hat keine Auswirkungen
    if (!PisNULL(DANISO)) {
      if (!isMonotone(next) || !isIsotropic(next->isoown) ||
	  PisNULL(DANISO) || cov->ncol[DANISO] != cov->nrow[DANISO])
	GERR("anisotropy parameter only allowed for simple models up to now");
    }

    assert(cov->calling != NULL); 
    double anisoScale;
  
    if (!PisNULL(DANISO)) {
      anisoScale = 1.0 / getMinimalAbsEigenValue(P(DANISO), cov->nrow[DANISO]);
      if (!R_FINITE(anisoScale)) 
	GERR("singular anisotropy matrices not allowed");
    } else anisoScale = 1.0;
  
    if (cov->calling->nr == SMITHPROC) {
      if ((err = STRUCT(next, newmodel)) == NOERROR && *newmodel != NULL) {
	assert(	(*newmodel)->calling == cov);
	Types type = 
	  TypeConsistency(PointShapeType, *newmodel, 0) ? PointShapeType :
	  TypeConsistency(RandomType, *newmodel, 0) ? RandomType :
	  TypeConsistency(ShapeType, *newmodel, 0) ? ShapeType :
	  OtherType;
	double scale =  PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
	
	if (type == RandomType) { // random locations given;
	  // so, it must be of pgs type (or standard):
	  
	  if ((err = CHECK_R(*newmodel, cov->tsdim)) != NOERROR) {
	    goto ErrorHandling;
	  }
	  dummy = *newmodel;

	  if ((err=addScales(&dummy, anisoScale, Scale, scale))!=NOERROR){
	    goto ErrorHandling;
	  }	  

	  *newmodel = NULL;     
	  assert(cov->vdim[0] == cov->vdim[1]);
	  if ((err = addPointShape(newmodel, cov, dummy, cov->tsdim, 
				   cov->vdim[0])) != NOERROR) {
	    goto ErrorHandling; 
	  }

	  ASSERT_NEWMODEL_NOT_NULL;
	  (*newmodel)->calling = cov;
	} else { // type not RandomType
	  if (type == PointShapeType && 
	      (err = addScales((*newmodel)->sub + PGS_LOC, anisoScale, Scale, 
			       scale)) != NOERROR) goto ErrorHandling;
	  if ((err = CHECK(*newmodel, cov->tsdim, cov->xdimprev, type, 
			   cov->domprev, cov->isoprev, cov->vdim, 
			   ROLE_MAXSTABLE))
	      != NOERROR) {
	    goto ErrorHandling;
	  }
	  if (type==PointShapeType) {	
	    if ((err = PointShapeLocations(*newmodel, cov)) != NOERROR) 
	      goto ErrorHandling;
	  } else if (type == ShapeType) {
	    dummy = *newmodel;
	    *newmodel = NULL;
	    // suche nach geeigneten locationen
	    if ((err = addScales(&local, anisoScale, Scale, scale))!=NOERROR)
	      goto ErrorHandling;
	    if ((err = addPointShape(newmodel, dummy, NULL, local,
				     cov->tsdim, cov->vdim[0]))
		!= NOERROR) 
	      goto ErrorHandling; 
	  } else BUG;
	} // ! randomtype
      } else { // S TRUCT does not return anything
	int err2;
	if ((err2 = addPointShape(newmodel, cov, NULL, cov->tsdim,
				  cov->vdim[0]))
	    != NOERROR) {
	  if (err == NOERROR) err = err2;
	  goto ErrorHandling; 
	}
	err = NOERROR;
      }

    } else { // not from RPsmith
      BUG;
      // 
     if ((err = STRUCT(next, newmodel)) > NOERROR) return err;
    }
   
    break;
  case ROLE_GAUSS :
    if (cov->key != NULL) COV_DELETE(&(cov->key));

    if (PrevLoc(cov)->distances) 
      GERR("distances do not allow for more sophisticated simulation methods");
    
    if ((err = STRUCT(next, newmodel)) > NOERROR) return err;

    addModel(newmodel, DOLLAR);
    if (!PisNULL(DVAR)) kdefault(*newmodel, DVAR, P0(DVAR));
    if (!PisNULL(DSCALE) ) kdefault(*newmodel, DSCALE, P0(DSCALE));
    if (!next->deterministic) GERR("random shapes not programmed yet");
    
    break;
 
  default :
    GERR2("%s : changes in scale/variance not programmed yet for '%s'", 
	  NICK(cov), ROLENAMES[cov->role]);      
    }
  
  ErrorHandling:
  
    if (dummy != NULL) COV_DELETE(&dummy);
    if (local != NULL) COV_DELETE(&local);

   
  return err;
}




int initS(cov_model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!
  cov_model *next = cov->sub[DOLLAR_SUB],
    *varM = cov->kappasub[DVAR],
    *scaleM = cov->kappasub[DSCALE],
    *anisoM = cov->kappasub[DAUSER],
    *projM = cov->kappasub[DPROJ];
  int 
    vdim = cov->vdim[0],
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    err = NOERROR;
  bool 
    angle = anisoM != NULL && CovList[anisoM->nr].check == checkAngle,
    maxstable = hasExactMaxStableRole(cov);// Realisationsweise 

  if (hasAnyShapeRole(cov)) { // Average  !! ohne maxstable selbst !!
    double
      var[MAXMPPVDIM],
      scale = PisNULL(DSCALE)  ? 1.0 : P0(DSCALE);
    int i,
      dim = cov->tsdim;
    

    if (!PisNULL(DPROJ) || !PisNULL(DANISO) || 
	projM!= NULL || (anisoM != NULL && (!angle || !anisoM->deterministic))){

       SERR("(complicated) anisotropy ond projection not allowed yet in Poisson related models");
     }


    // Achtung I-NIT_RANDOM ueberschreibt mpp.* !!
    if (varM != NULL) {
      int nm_neu = nm == 0 && !maxstable ? 1 : nm;
      if ((err = INIT_RANDOM(varM, nm_neu, s, P(DVAR))) != NOERROR) return err; 
      
      int nmP1 = varM->mpp.moments + 1;
      for (i=0; i<vdim; i++) {
	int idx = i * nmP1;
	var[i] = maxstable ? P0(DVAR) : varM->mpp.mM[idx + 1];      
      }
    } else for (i=0; i<vdim; var[i++] = P0(DVAR));

    if (scaleM != NULL) {
      if (dim + nm < 1) SERR("found dimension <= 0");
      int dim_neu = maxstable ? nm : (dim + nm) < 1 ? 1 : dim + nm; 
    if ((err = INIT_RANDOM(scaleM, dim_neu, s, P(DSCALE)))
	!= NOERROR) return err;
      scale = maxstable ? P0(DSCALE) : scaleM->mpp.mM[1];      
    }
    if ((err = INIT(next, nm, s)) != NOERROR) return err;


    for (i=0; i < nmvdim; i++) {
      cov->mpp.mM[i] = next->mpp.mM[i]; 
      cov->mpp.mMplus[i] = next->mpp.mMplus[i]; 
    }

    if (varM != NULL && !maxstable) {
      for (i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= varM->mpp.mM[i];
	cov->mpp.mMplus[i] *= varM->mpp.mMplus[i];
      }
    } else {
      int j, k;
      double pow_var;
      for (k=j=0; j<vdim; j++) { 
	pow_var = 1.0;
	for (i=0; i <= nm; i++, pow_var *= var[j], k++) {	
	  cov->mpp.mM[k] *= pow_var;
	  cov->mpp.mMplus[k] *= pow_var;
	}
      }
    }
    if (scaleM != NULL && !maxstable) {
      if (scaleM->mpp.moments < dim) SERR("moments can not be calculated.");
      int j, k,
	nmP1 = scaleM->mpp.moments + 1;
      for (k=j=0; j<vdim; j++) { 
	double pow_scale = scaleM->mpp.mM[dim + j * nmP1];
	for (i=0; i <= nm; i++, k++) {
	  cov->mpp.mM[k] *= pow_scale;
	  cov->mpp.mMplus[k] *= pow_scale;
	}
      }
    } else {
      double pow_scale = pow(scale, dim);
      for (i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= pow_scale;
	cov->mpp.mMplus[i] *= pow_scale;
      }
    }

    if (!PisNULL(DANISO)) {      
      if (cov->nrow[DANISO] != cov->ncol[DANISO])
	SERR("only square anisotropic matrices allowed");
      double invdet = fabs(1.0 / getDet(P(DANISO), cov->nrow[DANISO]));
      if (!R_FINITE(invdet))
	SERR("determinant of the anisotropy matrix could not be calculated.");
      
      for (i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= invdet;
	cov->mpp.mMplus[i] *= invdet;
      }

    }

    if (anisoM != NULL) {  
      double invdet;
      if (angle) {
	int 
	  ncol = anisoM->vdim[0],
	  nrow = anisoM->vdim[1];
	if (nrow != ncol) SERR("only square anisotropic matrices allowed");
	double 
	  *diag = PARAM(anisoM, ANGLE_DIAG);
	if (diag != NULL) {
	  invdet = 1.0;
	  for (i=0; i<ncol; i++) invdet /= diag[i];
	} else {
	  invdet = PARAM0(anisoM, ANGLE_RATIO);
	}
      } else {
	SERR("general anisotropy matrices not allowed yet");
      }
      
      invdet = fabs(invdet);      
      if (!R_FINITE(invdet))
	SERR("determinant of the anisotropy matrix could not be calculated.");
      
      for (i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= invdet;
	cov->mpp.mMplus[i] *= invdet;
      }
    }

    
 
    if (R_FINITE(cov->mpp.unnormedmass)) {
      if (vdim > 1) BUG;
      cov->mpp.unnormedmass = next->mpp.unnormedmass * var[0];
    } else {
      for (i=0; i<vdim; i++)
	cov->mpp.maxheights[i] = next->mpp.maxheights[i] * var[i];
    }
 
  }


  else if (cov->role == ROLE_GAUSS) {
    cov_model 
      *key = cov->key,
      *sub = key == NULL ? next : key;
    assert(sub != NULL);

    assert(key == NULL || ({PMI(cov);false;}));//
      
   
    if ((err=INIT(sub, 0, s)) != NOERROR) return err;
   
  }

  else if (cov->role == ROLE_BASE) {
    cov_model 
      *key = cov->key,
      *sub = key == NULL ? next : key;
    assert(sub != NULL);

    assert(key == NULL || ({PMI(cov);false;}));//
      
   
    if ((err=INIT(sub, 0, s)) != NOERROR) return err;
    
  } else { 
    if ((err=INIT(next, 0, s)) != NOERROR) return err; // e.g. from MLE
    // PMI(cov->calling);
    // printf("%s ", NAME(cov->calling));
    // SERR("initiation of scale and/or variance failed");
  }

 if ((err = TaylorS(cov)) != NOERROR) return err;

 return NOERROR;

}


void doS(cov_model *cov, gen_storage *s){
   cov_model
     *varM = cov->kappasub[DVAR],
     *scaleM = cov->kappasub[DSCALE];
   int i,
     vdim = cov->vdim[0];

    if (varM != NULL && !varM->deterministic && !isRandom(varM)) {
      assert(!PisNULL(DVAR));
      DO(varM, s);
    }
    
    if (scaleM != NULL && !scaleM->deterministic && !isRandom(scaleM)) {
      assert(!PisNULL(DSCALE));
      DO(scaleM, s);
    }

  if (hasAnyShapeRole(cov)) {
    cov_model *next = cov->sub[DOLLAR_SUB];
         
    DO(next, s);// nicht gatternr
    for (i=0; i<vdim; i++)
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * P0(DVAR);
    return;
  }  
  
  else if (cov->role == ROLE_GAUSS) {    
    double 
      *res = cov->rf,
      sd = sqrt(P0(DVAR));
    int 
      totalpoints = Gettotalpoints(cov);
    assert(res != NULL);
    if (cov->key == NULL) BUG;

    if (sd != 1.0) for (i=0; i<totalpoints; i++) res[i] *= (double) sd;
    
    return;
  } 

  BUG;
}



int checkplusmal(cov_model *cov) {
  cov_model *sub;
  int i, j, err, 
    vdim[2] = {1, 1},
    dim = cov->tsdim, 
    xdim = cov->xdimown, 
    role = cov->role;
  bool plus = CovList[cov->nr].check == checkplus,
    trend = isTrend(cov->typus);
  //  bool linearmodel = !plus && trend && cov->sub[0]->nr == CONST;
  Types covtype = cov->typus;
  domain_type covdom = trend ? XONLY : cov->domown;
  int trendiso = UpgradeToCoordinateSystem(cov->isoown);
  if (trendiso == ISO_MISMATCH) trendiso = cov->isoown;
  assert(trendiso != ISO_MISMATCH);
  int coviso = trend ? trendiso : cov->isoown;

  assert(cov->Splus == NULL);  

  // printf("\n\n %s %s\n", TYPENAMES[cov->typus], ISONAMES[cov->isoprev]);
  int variants = 1 +
    (int) (!trend && (cov->calling == NULL || !isShape(cov->calling)));


  // PMI(cov->calling);
  cov->matrix_indep_of_x = true;
  for (i=0; i<cov->nsub; i++) { 
    Types type = covtype;
    domain_type dom = covdom;
    int iso = coviso;
    sub = cov->sub[i];

    // PMI(sub); assert(sub->nr != CONST || sub->maxdim > -2);
   //      printf("> %s entering %s\n", NAME(cov), ISONAMES[coviso]);
    //PMI(cov->calling);
    //    assert(covdom == KERNEL);
    //print("stop\n"); if ( iso != EARTH_SYMMETRIC || covdom != KERNEL) return ERRORFAILED;
     
    if (sub == NULL) 
      SERR("+ or *: named submodels are not given in a sequence!");

    if (!plus && type==VariogramType) type=PosDefType;
  
    err = ERRORTYPECONSISTENCY;
    // printf("unten zu 1 !!");
    for (j=0; j<variants; j++) { // nur trend als abweichender typus erlaubt

      // printf(">> %s : %d %d  type =%s   %s \n", NAME(cov), j, variants, TYPENAMES[type], ISONAMES[ iso]);
      // int tt; printf("%d %s\n", tt=TypeConsistency(type, sub, 0), NAME(sub));
 
      if (TypeConsistency(type, sub, 0) &&
	  (err = CHECK(sub, dim, xdim, type, dom, iso, SUBMODEL_DEP, 
		       type == TrendType ? ROLE_BASE : role))
	  == NOERROR) break;
   
      if ((!isNegDef(type) && !isProcess(type)) || isTrend(type)) break;
      //    printf("trying trend %s %d %s %s %s!\n", NAME(sub), j, TYPENAMES[ type], DOMAIN_NAMES[dom], ISONAMES[iso]);     
      type = TrendType;
      dom = XONLY;
      iso = trendiso;
    
      //   
      //    if (j == variants -1) {
      //printf("trying trend %s %d of %d %s!\n", NAME(sub), j, variants, ISONAMES[iso]);
	// PMI(sub);     
	// 
      //      }
    }
    
    if (err != NOERROR) {
      // printf("here  dddd %s %s %d of %d; %d\n", NAME(cov), NAME(sub), j, variants, TypeConsistency(type, sub, 0)); MERR(err);
      return err;
    }

    if (cov->typus == sub->typus) {
      setbackward(cov, sub);
    } else {  
      updatepref(cov, sub);
      cov->tbm2num |= sub->tbm2num;
      cov->vdim[0] = sub->vdim[0];
      cov->vdim[1] = sub->vdim[1];
      cov->deterministic &= sub->deterministic;
    };
    for(j=0; j<2; j++) {
      if (vdim[j] == 1) {
	if (cov->vdim[j] != 1) vdim[j] = cov->vdim[j];
      } else {
	if (cov->vdim[j] != 1 && cov->vdim[j] != vdim[j])
 	SERR4("multivariate dimensionality is different in the submodels (%s is %d-variate; %s is %d-variate)", NICK(cov->sub[0]), cov->vdim[j], NICK(sub), sub->vdim[j]);
      }
    }
    //  if (vdim[1] > vdim[0]) SERR("unclear construction"); // at least currently not allowed 

    cov->matrix_indep_of_x &= sub->matrix_indep_of_x;
  } // i, nsub

  cov->vdim[0] = vdim[0];
  cov->vdim[1] = vdim[1];

  // !! incorrect  !!
  // cov->semiseparatelast = false; 
  //cov->separatelast = false; 

  return NOERROR;
}





// see private/old.code/operator.cc for some code including variable locations
void select(double *x, cov_model *cov, double *v) {
  int len,
    *element = PINT(SELECT_SUBNR);
  cov_model *sub = cov->sub[*element];
  assert(cov->vdim[0] == cov->vdim[1]);
  if (*element >= cov->nsub) ERR("select: element out of range");
  COV(x, sub, v);
  if ( (len = cov->nrow[SELECT_SUBNR]) > 1) {
    int i, m,
      vdim = cov->vdim[0],
      vsq = vdim * vdim;
    ALLOC_EXTRA(z, vsq);
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
    int element = P0INT(SELECT_SUBNR);
    cov_model *next = cov->sub[element];
    if (element >= cov->nsub) ERR("select: element out of range");
    CovList[next->nr].covmatrix(next, v);
      
  }  else StandardCovMatrix(cov, v);
}

char iscovmatrix_select(cov_model VARIABLE_IS_NOT_USED *cov) {  return 2; }

int checkselect(cov_model *cov) {
  int err;

  assert(cov->Splus == NULL);

  if (!isCartesian(cov->isoown)) return ERRORNOTCARTESIAN;
  kdefault(cov, SELECT_SUBNR, 0);
  if ((err = checkplus(cov)) != NOERROR) return err;

  if ((err = checkkappas(cov)) != NOERROR) return err;

  EXTRA_STORAGE;

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
    vdim = cov->vdim[0],
    vsq = vdim * vdim;
  //  assert(cov->vdim[0] == cov->vdim[1]);
  // printf("%d %d\n", cov->vdim[0], cov->vdim[1]);
  ALLOC_EXTRA(z, vsq);

  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];

    //      printf("PLUS %s %s ---> %s %s %d\n", NAME(sub), NICK(sub), TYPENAMES[cov->typus], TYPENAMES[sub->typus], TypeConsistency(cov->typus, sub->typus));

    if (TypeConsistency(cov->typus, sub->typus)) {
      //       pmi(sub, 0);   
      //
        
      FCTN(x, sub, z);

      //printf("%s %d  %s %d\n", NAME(cov), cov->typus, NAME(sub), cov->typus);
      //   assert(cov->typus != 9);

      if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] += z[0]; 
      else for (m=0; m<vsq; m++) v[m] += z[m]; 
     }
    // if (!R_FINITE(z[0])) { PMI(sub); printf("plus i=%d m=%d x=%f %f %f\n", i, m, x[0], v[0], z[0]); }
    // 
  }
  //  APMI(cov);
  //  assert(false);

}

void plusNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim[0] * cov->vdim[1];
  assert(cov->vdim[0] == cov->vdim[1]);
  ALLOC_EXTRA(z, vsq);
  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    if (TypeConsistency(cov->typus, sub->typus)) {
      NONSTATCOV(x, y, sub, z);
      if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] += z[0]; 
      else for (m=0; m<vsq; m++) v[m] += z[m]; 
    }
  }
}

void Dplus(double *x, cov_model *cov, double *v){
  cov_model *sub;
  int n = cov->nsub, i,
    vsq = cov->vdim[0] * cov->vdim[1];
  ALLOC_EXTRA(z, vsq);
  for (int m=0; m<vsq; m++) v[m] = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
    if (TypeConsistency(cov->typus, sub->typus)) {
      Abl1(x, sub, z);
      if (sub->vdim[0] == 1) for (int m=0; m<vsq; m++) v[m] += z[0]; 
      else for (int m=0; m<vsq; m++) v[m] += z[m]; 
    }
  }
}

void DDplus(double *x, cov_model *cov, double *v){
  cov_model *sub;
  int n = cov->nsub, i,
    vsq = cov->vdim[0] * cov->vdim[1];
  ALLOC_EXTRA(z, vsq);
  for (int m=0; m<vsq; m++) v[m] = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
    if (TypeConsistency(cov->typus, sub->typus)) {
      Abl2(x, sub, z);
      if (sub->vdim[0] == 1) for (int m=0; m<vsq; m++) v[m] += z[0]; 
      else for (int m=0; m<vsq; m++) v[m] += z[m]; 
    }
  }
}


int checkplus(cov_model *cov) {
  int err, i;
  if ((err = checkplusmal(cov)) != NOERROR) {
    return err;
  }
  
  if (cov->domown == DOMAIN_MISMATCH) return ERRORNOVARIOGRAM;
  if (cov->nsub == 0) cov->pref[SpectralTBM] = PREF_NONE;

  if (isPosDef(cov) && cov->domown == XONLY) cov->logspeed = 0.0;
  else if (isVariogram(cov) && cov->domown == XONLY) {
    cov->logspeed = 0.0;
    for (i=0; i<cov->nsub; i++) {
      cov_model *sub = cov->sub[i];
      if (TypeConsistency(cov->typus, sub->typus)) {
	if (ISNAN(sub->logspeed)) {
	  cov->logspeed = RF_NA;
	  break;
	} else cov->logspeed += sub->logspeed;
      }
    }
  } else cov->logspeed = RF_NA;

  EXTRA_STORAGE;

  return NOERROR;

  // spectral mit "+" funktioniert, falls alle varianzen gleich sind,
  // d.h. nachfolgend DOLLAR die Varianzen gleich sind oder DOLLAR nicht
  // auftritt; dann zufaellige Auswahl unter "+"
}
 

bool Typeplus(Types required, cov_model *cov, int depth) {
  // assert(false);
  bool allowed = TypeConsistency(ShapeType, required) || isTrend(required);
    //||required==ProcessType ||required==GaussMethodType; not yet allowed;to do

  if (!allowed) return false;
  int i;
  for (i=0; i<cov->nsub; i++) {
    if (TypeConsistency(required, cov->sub[i], depth-1)) return true;
  }  
  return false;
}

void spectralplus(cov_model *cov, gen_storage *s, double *e){
  assert(cov->vdim[0] == 1);
  int nr;
  double dummy;
  spec_properties *cs = &(s->spec);
  double *sd_cum = cs->sub_sd_cum;

  nr = cov->nsub - 1;
  dummy = UNIFORM_RANDOM * sd_cum[nr];
  if (ISNAN(dummy)) BUG;
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
      //assert(cov->calling != NULL && (isInterface(cov->calling->typus) ||
      ///				      isProcess(cov->calling->typus)));
      assert(cov->nr != PLUS_PROC);
      assert(cov->nr == PLUS);
      //cov->nr = PLUS_PROC;
      BUG;
      //return structplusproc(cov, newmodel); // kein S-TRUCT !!
    }
    if (cov->Splus != NULL) BUG;
    for (m=0; m<cov->nsub; m++) {
      cov_model *sub = cov->sub[m];
      if ((err = STRUCT(sub, newmodel))  > NOERROR) {
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


int initplus(cov_model *cov, gen_storage *s){
  int i, err,
    vdim = cov->vdim[0];
  if (cov->vdim[0] != cov->vdim[1]) BUG; // ??

  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_NA;

  if (cov->role == ROLE_GAUSS) {
    spec_properties *cs = &(s->spec);
    double *sd_cum = cs->sub_sd_cum;
 
    if (cov->vdim[0] == 1) {
      for (i=0; i<cov->nsub; i++) {
	cov_model *sub = cov->Splus == NULL ? cov->sub[i] : cov->Splus->keys[i];
      
	if (sub->pref[Nothing] > PREF_NONE) { // to do ??
	  // for spectral plus only
	  COV(ZERO, sub, sd_cum + i);
	  if (i>0) sd_cum[i] += sd_cum[i-1];
	}
	cov->sub[i]->Sgen = (gen_storage *) MALLOC(sizeof(gen_storage));
	if ((err = INIT(sub, cov->mpp.moments, s)) != NOERROR) {
	  //  AERR(err);
	  return err;
	}
	sub->simu.active = true;
      }
    } 
 
    cov->fieldreturn = cov->Splus != NULL;
    cov->origrf = false;
    if (cov->Splus != NULL) cov->rf = cov->Splus->keys[0]->rf;
     
    return NOERROR;
  }

  else if (cov->role == ROLE_COV) {    
    return NOERROR;
  }

  return ERRORFAILED;
}


void doplus(cov_model *cov, gen_storage *s) {
  int i;
  
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
  long i, 
    totalpoints = loc->totalpoints,
    vdimtot = totalpoints * cov->vdim[0],
    vdimtotSq = vdimtot * vdimtot;
  int
    nsub = cov->nsub;
  bool is = iscovmatrix_plus(cov) >= 2;
  double *mem = NULL;
  if (is && nsub > 1) {
    ALLOC_EXTRA2(MEM, vdimtotSq);
    mem = MEM;
    is = mem != NULL;
  }
  
  if (is) {
    int j;
    if (PisNULL(SELECT_SUBNR)) PALLOC(SELECT_SUBNR, 1, 1);
    P(SELECT_SUBNR)[0] = 0;
    CovList[SELECTNR].covmatrix(cov, v);
    for (i=1; i<nsub; i++) {
      if (Loc(cov->sub[i])->totalpoints != totalpoints) BUG;
      P(SELECT_SUBNR)[0] = i;
      CovList[SELECTNR].covmatrix(cov, mem);
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
    vdim = cov->vdim[0],
    vsq = vdim * vdim;
  ALLOC_EXTRA(z, vsq);
  
  //  assert(x[0] >= 0.0 || cov->xdimown > 1);
  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    COV(x, sub, z);
    if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] *= z[0];
    else for (m=0; m<vsq; m++) v[m] *= z[m];  
    //printf("%f %f\n", *x,  *v);
  }
}

void logmalStat(double *x, cov_model *cov, double *v, double *Sign){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vdim = cov->vdim[0],
    vsq = vdim * vdim;
  ALLOC_EXTRA(z, vsq);
  ALLOC_EXTRA2(zSign, vsq);

  assert(cov->vdim[0] == cov->vdim[1]); 
  assert(x[0] >= 0.0 || cov->xdimown > 1);
  for (m=0; m<vsq; m++) {v[m] = 0.0; Sign[m]=1.0;}
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    LOGCOV(x, sub, z, zSign);
    if (sub->vdim[0] == 1) {
      for (m=0; m<vsq; m++) {
	v[m] += z[0]; 
	Sign[m] *= zSign[0];
      }
    } else {
      for (m=0; m<vsq; m++) {
	v[m] += z[m]; 
	Sign[m] *= zSign[m];
      }
    }
  }
}

void malNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *sub;
  int i, m, nsub=cov->nsub,
    vdim = cov->vdim[0],
    vsq = vdim * vdim;
  ALLOC_EXTRA(z, vsq);

  assert(cov->vdim[0] == cov->vdim[1]);
  assert(cov->vdim[0] == cov->vdim[1]);

  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    NONSTATCOV(x, y, sub, z);
    if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] *= z[0];
    else for (m=0; m<vsq; m++) v[m] *= z[m];  
  }
}

void logmalNonStat(double *x, double *y, cov_model *cov, double *v, 
		   double *Sign){
  cov_model *sub;
  int i, m, nsub=cov->nsub,
    vdim = cov->vdim[0],
    vsq = vdim * vdim;
  assert(cov->vdim[0] == cov->vdim[1]);
  ALLOC_EXTRA(z, vsq);
  ALLOC_EXTRA2(zSign, vsq);
  for (m=0; m<vsq; m++) {v[m] = 0.0; Sign[m]=1.0;}
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    LOGNONSTATCOV(x, y, sub, z, zSign);
    if (sub->vdim[0] == 1) {
      for (m=0; m<vsq; m++) {
	v[m] += z[0]; 
	Sign[m] *= zSign[0];
      }
    } else {
      for (m=0; m<vsq; m++) {
	v[m] += z[m]; 
	Sign[m] *= zSign[m];
      }
    }
  }
}

void Dmal(double *x, cov_model *cov, double *v){
  cov_model *sub;
  int n = cov->nsub, i,
    vsq = cov->vdim[0] * cov->vdim[1];
  ALLOC_EXTRA3(c, vsq * MAXSUB);
  ALLOC_EXTRA4(d, vsq * MAXSUB);

  for (i=0; i<n; i++) {
    sub = cov->sub[i];
    COV(x, sub, c + i * vsq);
    Abl1(x, sub, d + i* vsq);
  }
  *v = 0.0;
  for (i=0; i<n; i++) {
    double *zw = d + i *vsq;
    int j;
    for (j=0; j<n; j++) {
      double *cc = c + j * vsq;
      if (j!=i) {
	for (int k=0; k<vsq; k++) zw[j] *= cc[j]; 
      }
    }
    for (int k=0; k<vsq; k++) v[k] += zw[k];
  }
}


int checkmal(cov_model *cov) {
  cov_model *next1 = cov->sub[0];
  cov_model *next2 = cov->sub[1];
  int i, err,
    nsub = cov->nsub;

  if (next2 == NULL) next2 = next1;

  //  printf("cov = %s\n", NAME(next2));

  if ((err = checkplusmal(cov)) != NOERROR) return err;

  //  printf("cov = %s no error \n", NAME(next2));

  bool ok = cov->domown != DOMAIN_MISMATCH &&
    (isTrend(cov->typus) || 
     (isShape(cov->typus) && (!isNegDef(cov->typus) || isPosDef(cov->typus) )));
  if (!ok) return ERRORNOVARIOGRAM;
   
  // to do istcftype und allgemeiner typ zulassen
  if (cov->typus == TrendType) {
    ok = false;
    for (i=0; i<nsub; i++)
      if ((ok = cov->sub[i]->nr == CONST || cov->sub[i]->nr == BIND)) break;
    if (!ok) SERR2("misuse as a trend function. At least one factor must be a constant (including 'NA') or a vector built with '%s(...)' or '%s(...).",
		  CovList[BIND].name,  CovList[BIND].nick);
  }

  cov->logspeed = cov->domown == XONLY ? 0 : RF_NA;
      

  if (cov->xdimown >= 2) cov->pref[TBM] = PREF_NONE;
  if (cov->xdimown==2 && cov->nsub == 2 && 
      isAnyDollar(next1) && isAnyDollar(next2)) {
    double *aniso1 = PARAM(next1, DANISO),
      *aniso2 = PARAM(next2, DANISO);
    if (aniso1 != NULL && aniso2 != NULL) {
      if (aniso1[0] == 0.0 && next1->ncol[DANISO] == 1) {
	cov->pref[TBM] = next2->pref[TBM];
      } else if (aniso2[0] == 0.0 && next2->ncol[DANISO] == 1) {
	cov->pref[TBM] = next1->pref[TBM];
      }
    }
  }
  
  if (cov->ptwise_definite <= last_pt_definite) {
    cov->ptwise_definite = next1->ptwise_definite;
    if (cov->ptwise_definite != pt_zero) {
      for (i=1; i<cov->nsub; i++) {
	cov_model *sub = cov->sub[i];
	if (sub->ptwise_definite == pt_zero) {
	  cov->ptwise_definite = pt_zero;
	  break;
	}
	if (sub->ptwise_definite != pt_posdef) {
	  if (sub->ptwise_definite == pt_negdef) {
	    cov->ptwise_definite = 
	      cov->ptwise_definite == pt_posdef ? pt_negdef
	      : cov->ptwise_definite == pt_negdef ? pt_posdef
	      : pt_indef;
	  } else { // sub = indef
	    cov->ptwise_definite = sub->ptwise_definite;
	    break;
	  }
	}
      }
    }
  }
	  
  EXTRA_STORAGE;

  return NOERROR;
}



bool Typemal(Types required, cov_model *cov, int depth) {
  // printf("%d %d\n", isShape(required), isTrend(required));

  if (!isShape(required) && !isTrend(required)) return false;
  int i;
  for (i=0; i<cov->nsub; i++) {
    //       print("Typemal %s %s %d\n", TYPENAMES[required], NAME(cov->sub[i]), TypeConsistency(required, cov->sub[i], depth-1));
    if (!TypeConsistency(required, cov->sub[i], depth-1)) return false;
  }
  //  print("Typemal OK\n");
  return true;
}



int initmal(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
//  int err;
//  return err;
  return ERRORFAILED;
  int i, 
    vdim = cov->vdim[0];
  if (cov->vdim[0] != cov->vdim[1]) BUG;

  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_NA;
}
void domal(cov_model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  BUG;
}


//////////////////////////////////////////////////////////////////////
// mpp-plus


#define PLUS_P 0 // parameter
int  CheckAndSetP(cov_model *cov){   
  int n,
    nsub = cov->nsub;
  double 
    cump = 0.0;
  if (PisNULL(PLUS_P)) {
    assert(cov->nrow[PLUS_P] == 0 && cov->ncol[PLUS_P] == 0);
    PALLOC(PLUS_P, nsub, 1);
    for (n=0; n<nsub; n++) P(PLUS_P)[n] = 1.0 / (double) nsub;    
  } else {
    cump = 0.0;
    for(n = 0; n<nsub; n++) {
      cump += P(PLUS_P)[n];
      if (cump > 1.0 && n+1<nsub) return ERRORATOMP; 
    }
    if (cump != 1.0) {
      if (nsub == 1) {
	warning("the p-values do not sum up to 1.\nHere only one p-value is given which must be 1.0");
	P(PLUS_P)[0] = 1.0;
      } else {
	if (cump < 1.0 && P(PLUS_P)[nsub-1] == 0) {
	  WARNING1("The value of the last component of '%s' is increased.",
		   KNAME(PLUS_P)); 

	}
	else SERR1("The components of '%s' do not sum up to 1.", KNAME(PLUS_P));
	P(PLUS_P)[nsub-1] = 1.0 - (cump - P(PLUS_P)[nsub-1]);
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
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim;
  cov_model *sub;
  ALLOC_EXTRA(z, vdimSq);

  if (cov->role == ROLE_COV) {  
    for (i=0; i<vdimSq; i++) v[i] = 0.0;
    for (n=0; n<cov->nsub; n++, sub++) {
      sub = cov->sub[n];
      COV(x, sub, z); // urspruenglich : covlist[sub].cov(x, cov, z); ?!
      for (i=0; i<vdimSq; i++) v[i] += P(PLUS_P)[n] * z[i];
    }
  } else {
    assert(hasPoissonRole(cov));
    BUG;
  }
}

int checkmppplus(cov_model *cov) { 
  int err, 
    size = 1;
  cov->maxdim = MAXMPPDIM;
  
  if ((err = checkplusmal(cov)) != NOERROR) {
    return err;
  }

  if ((err = CheckAndSetP(cov)) != NOERROR) return(err);

  if (cov->q == NULL) QALLOC(size);
  EXTRA_STORAGE;
  
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
  
  if (!hasMaxStableRole(cov) && !hasPoissonRole(cov)) {
    SERR("method is not based on Poisson point process");
  }


  return ERRORNOTPROGRAMMEDYET;

  
  // Ausnahme: mppplus wird separat behandelt:
  // if (nr == MPPPLUS) return S TRUCT(shape, NULL);
 
  ASSERT_NEWMODEL_NOT_NULL;
  NEW_STORAGE(plus);
  plus_storage *s = cov->Splus;

  for (m=0; m<cov->nsub; m++) {
    cov_model *sub = cov->sub[m];          
    if (s->keys[m] != NULL) COV_DELETE(s->keys + m);    
    if ((err = covCpy(s->keys + m, sub)) != NOERROR) return err;
    if ((err = addShapeFct(s->keys + m)) != NOERROR) return err;
    s->keys[m]->calling = cov;
  }
  return NOERROR;
}


int init_mppplus(cov_model *cov, gen_storage *S) {
  cov_model  *sub;
  double M2[MAXMPPVDIM], M2plus[MAXMPPVDIM], Eplus[MAXMPPVDIM], 
    maxheight[MAXMPPVDIM];
  int i,n, err,
    vdim = cov->vdim[0];
  if (cov->vdim[0] != cov->vdim[1]) BUG;
  if (vdim > MAXMPPVDIM) BUG;
  ext_bool 
    loggiven = SUBMODEL_DEP,
    fieldreturn = SUBMODEL_DEP;
  pgs_storage *pgs = NULL;
  for (i=0; i<vdim; i++) {
    maxheight[i] = RF_NEGINF;
    M2[i] = M2plus[i] = Eplus[i] = 0.0;
  }
    
  NEW_STORAGE(pgs);
  pgs = cov->Spgs;
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
    pgs->totalmass += sub->Spgs->totalmass * P(PLUS_P)[n];
    for (i=0; i<vdim; i++)   
      if (cov->mpp.maxheights[i] > maxheight[i]) 
	maxheight[i] = cov->mpp.maxheights[i];
    loggiven &= cov->loggiven;

    // Achtung cov->mpp.mM2 und Eplus koennten nicht direkt gesetzt
    // werden, da sie vom CovList[sub->nr].Init ueberschrieben werden !!
    
    if (cov->mpp.moments >= 1) {
      int nmP1 = sub->mpp.moments + 1;
      for (i=0; i<vdim; i++) {
	int idx = i * nmP1;
	Eplus[i] += PARAM0(sub, PLUS_P) * sub->mpp.mMplus[idx + 1]; 
      }
      if (cov->mpp.moments >= 2) {
	for (i=0; i<vdim; i++) {
	  int idx = i * nmP1;
	  M2[i] += PARAM0(sub, PLUS_P)  * sub->mpp.mM[idx + 2];
	  M2plus[i] += PARAM0(sub, PLUS_P) * sub->mpp.mM[idx + 2];
	}
      }
    }
    //assert(sub->mpp.loc_done);
  }
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = maxheight[i];
  //cov->mpp.refsd = RF_NA;
  //  cov->mpp.refradius = RF_NA;

  
  if (cov->mpp.moments >= 1) {
    int nmP1 = cov->mpp.moments + 1;
    for (i=0; i<vdim; i++) {
      int idx = i * nmP1;
      cov->mpp.mMplus[idx + 1] = Eplus[i];
      cov->mpp.mM[idx + 1] = RF_NA;
    }
    if (cov->mpp.moments >= 2) {
      for (i=0; i<vdim; i++) {
	 int idx = i * nmP1;
	 cov->mpp.mM[idx + 2] = M2[i];
	 cov->mpp.mMplus[idx + 2] = M2plus[i];
      }
    }
  }
  
  cov->loggiven = loggiven;
  cov->fieldreturn = fieldreturn;
  //cov->mpp.loc_done = true;       
  cov->origrf = false;
  cov->rf = NULL;

  return NOERROR;
}

void do_mppplus(cov_model *cov, gen_storage *s) {
  cov_model *sub;
  double subselect = UNIFORM_RANDOM;
  int i, subnr,
    vdim = cov->vdim[0];
  assert(cov->vdim[0] == cov->vdim[1]);
  for (subnr=0; (subselect -= PARAM0(cov->sub[subnr], PLUS_P)) > 0; subnr++);
  cov->q[0] = (double) subnr;
  sub = cov->sub[subnr];
  
  CovList[sub->nr].Do(sub, s);  // nicht gatternr
  for (i=0; i<vdim; i++)   
    cov->mpp.maxheights[i] = sub->mpp.maxheights[i];
  cov->fieldreturn = sub->fieldreturn;
  cov->loggiven = sub->loggiven;
}

//////////////////////////////////////////////////////////////////////
// PROCESSES
//////////////////////////////////////////////////////////////////////


int structSproc(cov_model *cov, cov_model **newmodel) {

  //  printf("entering\n");

  cov_model
    *next = cov->sub[DOLLAR_SUB],
    *Scale = cov->kappasub[DSCALE],
    *Aniso = cov->kappasub[DAUSER];
  int dim, err; 
  // cov_model *sub;
   assert(isDollarProc(cov));

   if ((Aniso != NULL && !Aniso->deterministic) ||
       (Scale != NULL && !Scale->deterministic)
       ) {
     SERR1("complicated models including arbitrary functions for '%s' cannot be simulated yet", KNAME(DAUSER));
  }

  switch (cov->role) {
  case ROLE_GAUSS :
    ASSERT_NEWMODEL_NULL;
    if (cov->key != NULL) COV_DELETE(&(cov->key));

    if (PrevLoc(cov)->distances) 
      SERR("distances do not allow for more sophisticated simulation methods");
    
    if (Aniso!= NULL) {
      //A
      // crash();
      Transform2NoGrid(cov, false, true, true); 
       
      location_type *loc = Loc(cov); // do not move before Transform2NoGrid!!
      dim = loc->timespacedim;
      int 
	bytes = dim * sizeof(double);
      long i,
	total = loc->totalpoints;
      
      //      PMI(Aniso);
      //      printf("%d\n", dim);
      
      if (dim != Aniso->vdim[0]) BUG;
      double *v = NULL,
	*x = loc->x;
      assert(x != NULL);
      assert(!loc->grid);
      if ((v = (double*) MALLOC(bytes)) == NULL) return ERRORMEMORYALLOCATION;
      for (i=0; i<total; i++, x+=dim) {
	FCTN(x, Aniso, v);
	MEMCOPY(x, v, bytes);
      }
      FREE(v);
    } else if (Scale != NULL && !isRandom(Scale)) {
      SERR1("Simulation algorithms for arbitrary scale functions do not exist yet -- try using arbitrary function for '%s'", KNAME(DANISO));
    } else {
      //pmi(cov, 0);      printf("next==intern %d\n", next->nr == TBM_PROC_INTERN);

    
      Transform2NoGrid(cov, true, // timesep
		       next->nr == TBM_PROC_INTERN ||
		       next->nr == NUGGET || 
		       next->nr == NUGGET_USER ||
		       next->nr == NUGGET_INTERN 
		       ? false : GRIDEXPAND_AVOID, 
		       true); // involveddollar
     }
    //  APMI(cov);
 
    if ((err = covCpy(&(cov->key), next)) != NOERROR) return err;

  
    if (!isGaussProcess(next)) addModel(&(cov->key), GAUSSPROC);
    SetLoc2NewLoc(cov->key, PLoc(cov));
    
    cov_model *key;
    key = cov->key;
    assert(key->calling == cov);    
    
    dim = PrevLoc(key)->timespacedim;
    if ((err = CHECK_NO_TRAFO(key, dim, dim, ProcessType, XONLY,
			      CoordinateSystemOf(cov->Sdollar->isoown),
			      cov->vdim[0], cov->role)) != NOERROR) {
      return err;
    }
    err = STRUCT(key, NULL);
    //   MERR(err);
    //   APMI(key);
   return err;
  default :
    SERR2("%s: changes in scale/variance not programmed yet for '%s'", 
	  NICK(cov), ROLENAMES[cov->role]);      
  }
     
  return NOERROR;
}




int initSproc(cov_model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!

  //cov_model *next = cov->sub[DOLLAR_SUB];
  //mppinfotype *info = &(s->mppinfo);
  //  location_type *loc = cov->prevloc;
  int 
    err = NOERROR;

  //  assert(false);
  
  cov_model 
    *key = cov->key;
    //*sub = key == NULL ? next : key;
  location_type *prevloc = PrevLoc(cov);

  assert(key != NULL);
  
  if ((err = INIT(key, 0, s)) != NOERROR) {
    return err;
  }
  
  key->simu.active = true; 
  assert(s != NULL);

  cov->fieldreturn = true;

  if ((cov->origrf = cov->ownloc != NULL &&
       Loc(cov)->totalpoints != prevloc->totalpoints)) {
    assert(prevloc->grid);
    int dim = prevloc->timespacedim;
    if (cov->vdim[0] != cov->vdim[1]) BUG;
    cov->rf = (double*) MALLOC(sizeof(double) *
				 cov->vdim[0] * 
				 prevloc->totalpoints);
    DOLLAR_STORAGE;
 
    int d,
      *proj = PINT(DPROJ),
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
      d = 0;
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
      double *A = P(DANISO);
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
  return NOERROR;
}


void doSproc(cov_model *cov, gen_storage *s){
  int i, 
    vdim = cov->vdim[0];

  if (hasMaxStableRole(cov) || hasPoissonRole(cov)) {
    assert(vdim == 1);
    cov_model *next = cov->sub[DOLLAR_SUB];
    
    cov_model
      *varM = cov->kappasub[DVAR],
      *scaleM = cov->kappasub[DSCALE];
 
    assert(cov->vdim[0] == cov->vdim[1]);

    if (varM != NULL && !varM->deterministic) {
      assert(!PisNULL(DVAR) && isRandom(varM));
      VTLG_R(NULL, varM, P(DVAR));
    }

    if (scaleM != NULL && !scaleM->deterministic) { // remote could be deterministic although local ist not
      assert(!PisNULL(DSCALE));
      VTLG_R(NULL, scaleM, P(DSCALE));
    }
    
    DO(next, s);// nicht gatternr
    for (i=0; i<vdim; i++)   
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * P0(DVAR);

  }

  else if (cov->role == ROLE_GAUSS) {    
    assert(cov->key != NULL);
    double 
      *res = cov->key->rf,
      sd = sqrt(P0(DVAR));
    int 
      totalpoints = Gettotalpoints(cov),
      totptsvdim = totalpoints * vdim;
    
    DO(cov->key, s); 
    
    cov_model
      *varM = cov->kappasub[DVAR];
    if (varM != NULL && !isRandom(varM)) {
      ALLOC_NEW(Sdollar, var, totptsvdim, var);
      Fctn(NULL, cov, var);
      for (i=0; i<totptsvdim; i++) {	
	res[i] *= sqrt(var[i]);
      }
    } else if (sd != 1.0) {
      for (i=0; i<totptsvdim; i++) res[i] *= sd;
    }
  }
  
  else BUG;
 
  if (cov->origrf) { // unklar was hier los ist; res bekommt nur einen kleineren  Teil
    if (vdim != 1) BUG;
    assert(PrevLoc(cov)->grid);
    int zaehler, d,
      dim = PrevLoc(cov)->timespacedim,
      *cumsum = cov->Sdollar->cumsum,
      *nx = cov->Sdollar->nx,
      *len = cov->Sdollar->len,
      *total = cov->Sdollar->total;
    assert(cov->key != NULL);
    double
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

  assert(cov->Splus != NULL);
  
  for (i=0; i<cov->nsub; i++) {

    sub = cov->Splus->keys[i];
     
    if (sub == NULL) 
      SERR("named submodels are not given in a sequence.");

    if (!TypeConsistency(type, sub, 0)) {
      return ERRORTYPECONSISTENCY;
    }
    if ((err= CHECK(sub, dim, xdim, type, dom, iso, SUBMODEL_DEP, role))
	!= NOERROR) {
      return err;
    }

   if (i==0) {
      cov->vdim[0]=sub->vdim[0];  // to do: inkonsistent mit vorigen Zeilen !!
      cov->vdim[1]=sub->vdim[1];  // to do: inkonsistent mit vorigen Zeilen !!
    } else {
      if (cov->vdim[0] != sub->vdim[0] || cov->vdim[1] != sub->vdim[1]) {
	SERR("multivariate dimensionality must be equal in the submodels");
      }
    }
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


int structplusmalproc(cov_model *cov, cov_model VARIABLE_IS_NOT_USED**newmodel){
  int m, err,
    dim =  Gettimespacedim(cov);
  //  bool plus = cov->nr == PLUS_PROC ;

 
  switch(cov->role) {
  case ROLE_GAUSS : {
    location_type *loc = Loc(cov);
    NEW_STORAGE(plus);
    plus_storage *s =cov->Splus;
    for (m=0; m<cov->nsub; m++) {
      cov_model *sub = cov->sub[m];
      if (s->keys[m] != NULL) COV_DELETE(s->keys + m);
      if ((err =  covCpy(s->keys + m, sub)) != NOERROR) {
	return err;
      }
      assert(s->keys[m] != NULL);
      assert(s->keys[m]->calling == cov);
      
      if (PL >= PL_STRUCTURE) {
	LPRINT("plus: trying initialisation of submodel #%d (%s).\n", m+1, 
	       NICK(sub));
      }
      
      isotropy_type newiso = UpgradeToCoordinateSystem(cov->isoprev); //!prev
      if (newiso == ISO_MISMATCH) {	  
	SERR2("'%s' for '%s' cannot be upgraded to a coordinate system",
	      NAME(sub), ISONAMES[cov->isoown]);
      }

      addModel(s->keys + m, isTrend(sub->typus) ? TRENDEVAL : GAUSSPROC);
      
      if (isTrend(sub->typus) && sub->Spgs == NULL) {
	// printf("**************** alloc_cov %s\n", NAME(sub));
	if ((err = alloc_cov(sub, dim, sub->vdim[0], sub->vdim[1])) != NOERROR)
	  return err;
      }

      s->keys[m]->calling = cov;      
      if ((err = CHECK(s->keys[m], loc->timespacedim, loc->timespacedim,
		       ProcessType, XONLY,
		       newiso, 
		       cov->vdim,
		       ROLE_GAUSS)) != NOERROR) {
	return err;
      }
     
      if ((s->struct_err[m] =
	   err = STRUCT(s->keys[m], NULL))  > NOERROR) {

	//PMI(s->keys[m]);
	// XERR(err);

	return err;
      }
      
    }
    return NOERROR;
  }
    
  default :
    SERR2("role '%s' not allowed for '%s'", ROLENAMES[cov->role],
	  NICK(cov));
  }

  return NOERROR;
}


int structplusproc(cov_model *cov, cov_model **newmodel){
  assert(cov->nr == PLUS_PROC);
  return structplusmalproc(cov, newmodel);
}


int structmultproc(cov_model *cov, cov_model **newmodel){
  assert(CovList[cov->nr].check == checkmultproc);
  return structplusmalproc(cov, newmodel);
}

int initplusmalproc(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  int i, err,
    vdim = cov->vdim[0];
  assert(cov->vdim[0] == cov->vdim[1]);
  bool plus = cov->nr == PLUS_PROC ;

  for (i=0; i<vdim; i++)   
    cov->mpp.maxheights[i] = RF_NA;
  if (cov->Splus == NULL) BUG;

  if (cov->role == ROLE_GAUSS) {
 
    for (i=0; i<cov->nsub; i++) {
      cov_model *sub = cov->Splus == NULL ? cov->sub[i] : cov->Splus->keys[i];
      if (!plus && 
	  (sub->nr == CONST 
	   //|| CovList[sub[0]->nr].check == checkconstant ||
	   // (isDollar(sub) && CovList[sub->sub[0]->nr].check == checkconstant)
	   ))
	continue;
      assert(cov->sub[i]->Sgen==NULL);
      cov->sub[i]->Sgen = (gen_storage *) MALLOC(sizeof(gen_storage));
      if ((err = INIT(sub, 0, cov->sub[i]->Sgen)) != NOERROR) {
	return err;
      }
      sub->simu.active = true;
    }
    cov->simu.active = true;
    return NOERROR;
  }
    
  else {
    BUG;
  }

    return ERRORFAILED;
}

 
int initplusproc(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  int err;
  if ((err = initplusmalproc(cov, s)) != NOERROR) return err;

  if (cov->role == ROLE_GAUSS) {
    cov->fieldreturn = cov->Splus != NULL;
    cov->origrf = false;
    if (cov->Splus != NULL) cov->rf = cov->Splus->keys[0]->rf;
     
    return NOERROR;
  }
     
  else {
    BUG;
  }

  return ERRORFAILED;
}


void doplusproc(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int m, i,
    total = Loc(cov)->totalpoints * cov->vdim[0];
  double *res = cov->rf;
  assert(cov->vdim[0] == cov->vdim[1]);

  if (cov->role == ROLE_GAUSS && cov->method==SpectralTBM) {
    ERR("error in doplus with spectral");
  }
  assert(cov->Splus != NULL);

  for (m=0; m<cov->nsub; m++) {
    cov_model *key = cov->Splus->keys[m],
      *sub = cov->sub[m];
    double *keyrf = key->rf;
    //    PMIL(key, 2);
    DO(key, sub->Sgen);
    if (m > 0)
      for(i=0; i<total; i++) res[i] += keyrf[i];
  }
  return;
}



#define MULTPROC_COPIES 0
int checkmultproc(cov_model *cov) {
  int err;
  kdefault(cov, MULTPROC_COPIES, GLOBAL.special.multcopies);
  if ((err = checkplusmalproc(cov)) != NOERROR) {
    return err;
  }
  EXTRA_STORAGE;
  return NOERROR;
}


int initmultproc(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  int  err;

  if ((err = initplusmalproc(cov, s)) != NOERROR) {
    BUG;
   return err;
  }

  if (cov->role == ROLE_GAUSS) {
    FieldReturn(cov);
    return NOERROR;
  }
     
  else {
    BUG;
  }

  return ERRORFAILED;
}



void domultproc(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double *res = cov->rf;
  assert(cov->vdim[0] == cov->vdim[1]);
  int  m, i, c, idx,
    vdim = cov->vdim[0],
   //  vdimSq= vdim * vdim,
    total = Loc(cov)->totalpoints,
   totalvdim = total * vdim,
   copies = GLOBAL.special.multcopies,
   factors = 0;

  if (cov->role == ROLE_GAUSS && cov->method==SpectralTBM) {
    ERR("error in do_mult with spectral");
  }

  if (cov->nsub==2 && 
      ((cov->sub[0]->nr==PROD) xor (idx=cov->sub[1]->nr==PROD)) &&
      cov->sub[0]->nr != CONST &&  
      cov->sub[1]->nr != CONST) {
    // koennte noch allgemeiner gemacht werden in verbindung mit CONST; todo
    copies = 1;
    cov->sub[idx]->Sgen->prodproc_random = false;
  }

  ///   APMI(cov);

  assert(cov->Splus != NULL);

  SAVE_GAUSS_TRAFO;
  for (c=0; c<copies; c++) {
    for(i=0; i<totalvdim; res[i++] = 1.0);
    for (m=0; m<cov->nsub; m++) {

      if (PL > PL_REC_DETAILS) {
	PRINTF("\rcopies=%d sub=%d\n", c, m);
      }

      cov_model *key = cov->Splus->keys[m],
	*sub = cov->sub[m];
      double *keyrf = key->rf;
      if (sub->nr == CONST) {
	double 
	  cc = isTrend(sub->typus) ? PARAM0(sub, CONST_C) 
	  : sqrt(PARAM0(sub, CONST_C));
	for(i=0; i<totalvdim; i++) res[i] *= cc;
      }
      /* else {
	bool dollar = isDollar(sub);
	cov_model *Const = dollar ? Const->sub[0] : sub;
	if (CovList[Check->nr].check == checkconstant) {
	  if (dollar) {
	    double var = PARAM0(sub, DVAR);
	    bool random = false;
	    if (sub->kappasub[DVAR] != NULL) {
	      if (random = isRandom(sub->kappasub[DVAR])) {
		Do(sub->kappasub[DVAR], sub->Sgen);
	      } else {
		ALLOC_EXTRA2(VarMem, totalvdim);
		F ct n(NULL, sub->kappasub[DVAR], VarMem);
	      }
	    }
	    if (var != 1.0) {
	      double sd = sqrt(var);
	      for(i=0; i<totalvdim; i++) res[i] *= sd;
	    }
	  } else {
	    ALLOC_EXTRA3(ConstMem, vdimSq);
	  }
	  } */
      else {	  
	factors ++;
	DO(key, sub->Sgen);
	for(i=0; i<totalvdim; i++) {
	  res[i] *= keyrf[i];
	}
      }
    }
    if (factors == 1) return; // no error, just exit
    if (c == 0) {      
      ALLOC_EXTRA(z, totalvdim);
      res = z;
    } else {
      for(i=0; i<totalvdim; i++) cov->rf[i] += res[i];
    }
  }

  double f;
  f = 1 / sqrt((double) copies);
  for(i=0; i<totalvdim; i++) cov->rf[i] *= f;

}



void rangemultproc(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[MULTPROC_COPIES] = 1.0;
  range->max[MULTPROC_COPIES] = RF_INF;
  range->pmin[MULTPROC_COPIES] = 1.0;
  range->pmax[MULTPROC_COPIES] = 1000;
  range->openmin[MULTPROC_COPIES] = false;
  range->openmax[MULTPROC_COPIES] = true;
}


// $power

void PowSstat(double *x, cov_model *cov, double *v){
  logPowSstat(x, cov, v, NULL);
}

void logPowSstat(double *x, cov_model *cov, double *v, double *Sign){
  cov_model *next = cov->sub[POW_SUB];
  double 
    factor,
    var = P0(POWVAR),
    scale =P0(POWSCALE), 
    p = P0(POWPOWER),
    invscale = 1.0 / scale;
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim *vdim,
    xdimown = cov->xdimown;
  assert(cov->vdim[0] == cov->vdim[1]);
  ALLOC_DOLLAR(z, xdimown);

  for (i=0; i < xdimown; i++) z[i] = invscale * x[i];
  if (Sign==NULL) {
    COV(z, next, v);
    factor = var * pow(scale, p);
    for (i=0; i<vdimSq; i++) v[i] *= factor; 
  } else {
    LOGCOV(z, next, v, Sign);
    factor = log(var) + p * log(scale);
    for (i=0; i<vdimSq; i++) v[i] += factor; 
  }
}

void PowSnonstat(double *x, double *y, cov_model *cov, double *v){
  logPowSnonstat(x, y, cov, v, NULL);
}

void logPowSnonstat(double *x, double *y, cov_model *cov, double *v, 
		 double *Sign){
  cov_model *next = cov->sub[POW_SUB];
  double 
    factor,
    var = P0(POWVAR),
    scale =P0(POWSCALE), 
    p = P0(POWPOWER),
     invscale = 1.0 / scale;
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim,
    xdimown = cov->xdimown;
  assert(cov->vdim[0] == cov->vdim[1]);
  ALLOC_DOLLARY(z1, z2, xdimown);

  for (i=0; i<xdimown; i++) {
    z1[i] = invscale * x[i];
    z2[i] = invscale * y[i];
  }

  if (Sign==NULL) {
    NONSTATCOV(z1, z2, next, v);
    factor = var * pow(scale, p);
    for (i=0; i<vdimSq; i++) v[i] *= factor; 
  } else {
    LOGNONSTATCOV(z1, z2, next, v, Sign);
    factor = log(var) + p * log(scale);
    for (i=0; i<vdimSq; i++) v[i] += factor; 
  }
}

 
void inversePowS(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[POW_SUB];
  int i,
    vdim = cov->vdim[0],
    vdimSq = vdim * vdim;
  double y, 
    scale =P0(POWSCALE),
    p = P0(POWPOWER),
    var = P0(POWVAR);
 assert(cov->vdim[0] == cov->vdim[1]);

  y= *x / (var * pow(scale, p)); // inversion, so variance becomes scale
  if (CovList[next->nr].inverse == ErrCov) BUG;
  INVERSE(&y, next, v);
 
  for (i=0; i<vdimSq; i++) v[i] *= scale; 
}


int TaylorPowS(cov_model *cov) {
  if (cov->vdim[0] != 1) SERR("Taylor only known in the unvariate case");
  cov_model 
    *next = cov->sub[POW_SUB];
  int i;
  double scale = PisNULL(POWSCALE) ? 1.0 : P0(POWSCALE);
  cov->taylorN = next->taylorN;  
  for (i=0; i<cov->taylorN; i++) {
    cov->taylor[i][TaylorPow] = next->taylor[i][TaylorPow];
    cov->taylor[i][TaylorConst] = next->taylor[i][TaylorConst] *
      P0(POWVAR) * pow(scale, P0(POWPOWER) - next->taylor[i][TaylorPow]);   
  }
  
  cov->tailN = next->tailN;  
  for (i=0; i<cov->tailN; i++) {
    cov->tail[i][TaylorPow] = next->tail[i][TaylorPow];
    cov->tail[i][TaylorExpPow] = next->tail[i][TaylorExpPow];
    cov->tail[i][TaylorConst] = next->tail[i][TaylorConst] *
      P0(POWVAR) * pow(scale, P0(POWPOWER) - next->tail[i][TaylorPow]);   
    cov->tail[i][TaylorExpConst] = next->tail[i][TaylorExpConst] *
      pow(scale, -next->tail[i][TaylorExpPow]);
  }
  return NOERROR;
}


int checkPowS(cov_model *cov) {
  // hier kommt unerwartet  ein scale == nan rein ?!!
  cov_model 
    *next = cov->sub[POW_SUB];
  int err, 
    tsdim = cov->tsdim,
    xdimown = cov->xdimown,
    xdimNeu = xdimown;

  if (!isCartesian(cov->isoown)) return ERRORNOTCARTESIAN;
    
  kdefault(cov, POWVAR, 1.0);
  kdefault(cov, POWSCALE, 1.0);
  kdefault(cov, POWPOWER, 0.0);
  if ((err = checkkappas(cov)) != NOERROR) {
    return err;
  }
  
  if ((err = CHECK(next, tsdim, xdimNeu, cov->typus, cov->domown,
		   cov->isoown, SUBMODEL_DEP, cov->role)) != NOERROR) {
    return err;
  }

  setbackward(cov, next);
  if ((err = TaylorPowS(cov)) != NOERROR) return err;

  DOLLAR_STORAGE;
  
  return NOERROR;
}



bool TypePowS(Types required, cov_model *cov, int depth) {
  cov_model *next = cov->sub[0];
  return (isShape(required) || isTrend(required) || isProcess(required))
    && TypeConsistency(required, next, depth-1);
}


void rangePowS(cov_model *cov, range_type* range){
  range->min[POWVAR] = 0.0;
  range->max[POWVAR] = RF_INF;
  range->pmin[POWVAR] = 0.0;
  range->pmax[POWVAR] = 100000;
  range->openmin[POWVAR] = false;
  range->openmax[POWVAR] = true;

  range->min[POWSCALE] = 0.0;
  range->max[POWSCALE] = RF_INF;
  range->pmin[POWSCALE] = 0.0001;
  range->pmax[POWSCALE] = 10000;
  range->openmin[POWSCALE] = true;
  range->openmax[POWSCALE] = true;

  range->min[POWPOWER] = RF_NEGINF;
  range->max[POWPOWER] = RF_INF;
  range->pmin[POWPOWER] = -cov->tsdim;
  range->pmax[POWPOWER] = +cov->tsdim;
  range->openmin[POWPOWER] = true;
  range->openmax[POWPOWER] = true;
 }



void PowScaleToLoc(cov_model *to, cov_model *from, int VARIABLE_IS_NOT_USED depth) {
  assert(!PARAMisNULL(to, LOC_SCALE) && !PARAMisNULL(from, POWSCALE));
  PARAM(to, LOC_SCALE)[0] = PARAM0(from, POWSCALE);
}

int structPowS(cov_model *cov, cov_model **newmodel) {
  cov_model
    *next = cov->sub[POW_SUB],
    *scale = cov->kappasub[POWSCALE];
  int err; 

  if (!next->deterministic) SERR("random shapes not programmed yet");

  switch (cov->role) {
  case ROLE_SMITH :  case ROLE_GAUSS :
    ASSERT_NEWMODEL_NOT_NULL;
    
    if ((err = STRUCT(next, newmodel)) > NOERROR) return err;
    
    addModel(newmodel, POWER_DOLLAR);
    if (!PisNULL(POWVAR)) kdefault(*newmodel, POWVAR, P0(POWVAR));
    if (!PisNULL(POWSCALE)) kdefault(*newmodel, POWSCALE, P0(POWSCALE));
    if (!PisNULL(POWPOWER)) kdefault(*newmodel, POWPOWER, P0(POWPOWER));
    
    break;
  case ROLE_MAXSTABLE : {
    ASSERT_NEWMODEL_NOT_NULL;
    
    if ((err = STRUCT(next, newmodel)) > NOERROR) return err;
    
    if (!isRandom(scale)) SERR("unstationary scale not allowed yet");
    addModel(newmodel, LOC);
    addSetDistr(newmodel, scale, PowScaleToLoc, true, MAXINT);
  }
    break;
  default :
    SERR2("'%s': changes in scale/variance not programmed yet for '%s'", 
	  NICK(cov), ROLENAMES[cov->role]);      
  }  
   
  return NOERROR;
}




int initPowS(cov_model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!
  cov_model *next = cov->sub[POW_SUB],
    *varM = cov->kappasub[POWVAR],
    *scaleM = cov->kappasub[POWSCALE];
  int 
    vdim = cov->vdim[0],
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    err = NOERROR;
  bool 
    maxstable = hasExactMaxStableRole(cov);// Realisationsweise 
  assert(cov->vdim[0] == cov->vdim[1]);


  assert(cov->key == NULL || ({PMI(cov);false;}));//   
  
  if (hasAnyShapeRole(cov)) { // !! ohne maxstable selbst !!
    double
      var[MAXMPPVDIM],  
      p = P0(POWPOWER),
      scale = P0(POWSCALE);
    int i,
      intp = (int) p,
      dim = cov->tsdim;

    

    // Achtung I-NIT_RANDOM ueberschreibt mpp.* !!
    if (varM != NULL) {
      int nm_neu = nm == 0 && !maxstable ? 1 : nm;
      if ((err = INIT_RANDOM(varM, nm_neu, s, P(POWVAR))) != NOERROR) 
	return err;
      int nmP1 = varM->mpp.moments + 1;
      for (i=0; i<vdim; i++) {
	int idx = i * nmP1;
	var[i] = maxstable ? P0(DVAR) : varM->mpp.mM[idx + 1];      
      }
    } else for (i=0; i<vdim; var[i++] = P0(POWVAR));

    if (scaleM != NULL) {
      if (p != intp)
	SERR1("random scale can be initialised only for integer values of '%s'",
	     KNAME(POWPOWER));
      if (dim + intp < 0) SERR("negative power cannot be calculated yet");
      int idx_s = maxstable ? nm : dim + nm + intp < 1 ? 1 : dim + nm + intp;
      if ((err = INIT_RANDOM(scaleM, idx_s, s, P(POWSCALE))) != NOERROR) return err;
      assert(scaleM->mpp.moments == 1);
      scale = maxstable ? P0(DSCALE) : scaleM->mpp.mM[1];      
    }
    if ((err = INIT(next, nm, s)) != NOERROR) return err;


    for (i=0; i < nmvdim; i++) {
      cov->mpp.mM[i] = next->mpp.mM[i];
      cov->mpp.mMplus[i] = next->mpp.mMplus[i];
   }


    if (varM != NULL && !maxstable) {
      int j,
	nmP1 = varM->mpp.moments + 1;
      for (j=0; j<vdim; j++) {
	int idx = j * nmP1;
	for (i=0; i <= nm; i++) {
	  cov->mpp.mM[i] *= varM->mpp.mM[idx + i];
	  cov->mpp.mMplus[i] *= varM->mpp.mMplus[idx + i];
	}
      }
    } else {
      int j, k;
      double pow_var;
      for (k=j=0; j<vdim; j++) { 
	pow_var = 1.0;
	for (i=0; i<=nm; i++, pow_var *= var[j], k++) {	
	  cov->mpp.mM[k] *= pow_var;
	  cov->mpp.mMplus[k] *= pow_var;
	}
      }
    }

    if (scaleM != NULL && !maxstable) {
      if (dim + nm * intp < 0 || dim + intp * nm > scaleM->mpp.moments) 
	SERR("moments cannot be calculated");
      assert(scaleM->vdim[0] == 1 && scaleM->vdim[1] == 1 );
      for (i=0; i <= nm; i++) {
	int idx = dim + i * intp;
	cov->mpp.mM[i] *= scaleM->mpp.mM[idx];
	cov->mpp.mMplus[i] *= scaleM->mpp.mMplus[idx];
      }
    } else {
      int j,k;
      double 
	pow_scale,
	pow_s = pow(scale, dim),
	pow_p = pow(scale, p);
      for (k=j=0; j<vdim; j++) { 
	pow_scale = pow_s;
	for (i=0; i <= nm; i++, pow_scale *= pow_p, k++) {
	  cov->mpp.mM[k] *= pow_scale;
	  cov->mpp.mMplus[k] *= pow_scale;
	}
      }
    }
 
    if (R_FINITE(cov->mpp.unnormedmass)) {
      if (vdim > 1) BUG;
      cov->mpp.unnormedmass = next->mpp.unnormedmass * var[0] / pow(scale, p);
    } else {
      double pp = pow(scale, -p);
      for (i=0; i<vdim; i++)   
	cov->mpp.maxheights[i] = next->mpp.maxheights[i] * var[i] * pp;
    }
 
  }


  else if (cov->role == ROLE_GAUSS) {  
    if ((err=INIT(next, 0, s)) != NOERROR) return err;
  }

  else if (cov->role == ROLE_BASE) {
    if ((err=INIT(next, 0, s)) != NOERROR) return err;
    
  }

  else SERR("Initiation of scale and/or variance failed");

 
  if ((err = TaylorPowS(cov)) != NOERROR) return err;

  return NOERROR;
}


void doPowS(cov_model *cov, gen_storage *s){
 
  if (hasAnyShapeRole(cov)) {
    cov_model *next = cov->sub[POW_SUB];
         
    DO(next, s);// nicht gatternr
    double factor = P0(POWVAR) / pow(P0(POWSCALE), P0(POWPOWER));
    int i,
      vdim = cov->vdim[0];
    assert(cov->vdim[0] == cov->vdim[1]);
    for (i=0; i<vdim; i++)   
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * factor;
    return;
  }
  
  BUG;
}
