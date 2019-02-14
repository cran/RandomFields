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
 Copyright (C) 2005 -- 2017 Martin Schlather

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


//#include <Rmath.h>
//#include <stdio.h> 
//#include <R_ext/Lapack.h>
#include "questions.h"
#include "operator.h"
#include "Processes.h"
#include "startGetNset.h"
#include "variogramAndCo.h"
#include "rf_interfaces.h"
#include "shape.h"
#include "primitive.others.h"


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of auxiliary correlation functions 

 Copyright (C) 2017 -- 2017 Martin Schlather

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



// $

bool hasVarOnly(model *cov) {
   if (cov == NULL || !isDollar(cov)) BUG;
  if (!PisNULL(DSCALE) && P0(DSCALE) != 1.0) return false;
  if (!PisNULL(DANISO) || !PisNULL(DPROJ)) return false;
  int i,
    kappas = DefList[COVNR].kappas;
  for (i=0; i<kappas; i++)
    if (cov->kappasub[i] != NULL) return false;
  return true;
}


void kappaS(int i, model *cov, int *nr, int *nc){
  //  PMI(cov);
  switch(i) {
  case DVAR : case DSCALE :
    *nr = *nc = 1;
    break;
  case DANISO :
    *nr = OWNTOTALXDIM;
    *nc = SIZE_NOT_DETERMINED;
    break;
  case DAUSER :
    *nr = SIZE_NOT_DETERMINED;
    *nc = OWNTOTALXDIM;
    break;
  case DPROJ : 
    *nr = SIZE_NOT_DETERMINED;
    *nc = 1;
    break;
  default : *nr = *nc = OUT_OF_RANGE; 
  }
}

// simple transformations as variance, scale, anisotropy matrix, etc.  
void Siso(double *x, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  double y,
    *aniso=P(DANISO),
    *scale =P(DSCALE),
    var = P0(DVAR);
  assert(cov->Sdollar->simplevar);
  
  y = *x;
  if (aniso != NULL) y = FABS(y * aniso[0]);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] : (y==0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  // letzteres darf nur passieren wenn dim = 1!!
  COV(&y, next, v);

  for (i=0; i<vdimSq; i++) v[i] *= var; 
}
  

// simple transformations as variance, scale, anisotropy matrix, etc.  
void logSiso(double *x, model *cov, double *v, double *Sign){
  model *next = cov->sub[DOLLAR_SUB];
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  double y, 
    *aniso=P(DANISO),
    *scale =P(DSCALE),
    logvar = LOG(P0(DVAR));
   assert(cov->Sdollar->simplevar);

  y = *x;
  if (aniso != NULL) y = FABS(y * aniso[0]);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] 
      : (y == 0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  LOGCOV(&y, next, v, Sign);
  for (i=0; i<vdimSq; i++) v[i] += logvar; 
  
 
}
 
void Sstat(double *x, model *cov, double *v){
  logSstat(x, cov, v, NULL);
}




void logSstat(double *x, model *cov, double *v, double *Sign){
  assert(cov->kappasub[DSCALE] == NULL && 
	 (cov->kappasub[DANISO]==NULL || 
	  DefList[MODELNR(cov->kappasub[DANISO])].check==checkAngle)&& 
	 (cov->kappasub[DAUSER]==NULL || 
	  DefList[MODELNR(cov->kappasub[DAUSER])].check==checkAngle));
  model *next = cov->sub[DOLLAR_SUB];
    //  *Aniso = cov->kappasub[DAUSER],
    // *Scale = cov->kappasub[DSCALE];
  double *z1 = NULL,
    *scale =P(DSCALE), 
    *aniso=P(DANISO);
  int i,
    nproj = Nproj,
    vdim = VDIM0,
    vdimSq = vdim * vdim;

  if (nproj > 0) {
    int *proj = PPROJ; 
    TALLOC_X1(z, nproj);
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
    //    A LLOC_DOLLAR(z, dim);
    //    FCTN(x, Aniso, z);
    //    z1 = z;
    //  } else if (Scale != NULL) {
    //    int dim = Aniso->vdim[0];
    //    A LLOC_DOLLAR(z, dim);
    //    FCTN(x, Aniso, z);
    //    z1 = z;
    END_TALLOC_X1;
  } else if (aniso==NULL && (scale == NULL || scale[0] == 1.0)) {
    z1 = x;
  } else {
    int xdimown = OWNTOTALXDIM;
    double *xz;
    TALLOC_X1(z, xdimown);
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
    END_TALLOC_X1
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
    var = LOG(var);
    for (i=0; i<vdimSq; i++) v[i] += var; 
  }

}

void Snonstat(double *x, double *y, model *cov, double *v){
  logSnonstat(x, y, cov, v, NULL);
}

void logSnonstat(double *x, double *y, model *cov, double *v, double *Sign){
  model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE];
  double 
    s1 = RF_NA, s2 = RF_NA, smeanSq=RF_NA,
    *scale =P(DSCALE),
    *aniso=P(DANISO);
  int i,
    //   nr = NA_INTEGER,
    nproj = Nproj,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  bool notdelete = false;
  TALLOC_DOUBLE(z1);
  TALLOC_DOUBLE(z2);

  //  printf("%ld %ld\n", z1, z2);
  
  if (nproj > 0) {
    int *proj = PPROJ;
    TALLOC_GLOBAL_X1(z1, nproj);
    TALLOC_GLOBAL_X2(z2, nproj);
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
    TALLOC_GLOBAL_X1(z1, dim);
    TALLOC_GLOBAL_X2(z2, dim);
    FCTN(x, Aniso, z1);
    FCTN(y, Aniso, z2);
  } else if (Scale != NULL && !isnowRandom(Scale)) {// see also code below variance
    int xdimown = OWNTOTALXDIM;
    double s;    
    //    nr = MODELNR(Scale);

    TALLOC_GLOBAL_X1(z1, xdimown);
    TALLOC_GLOBAL_X2(z2, xdimown);     
    FCTN(x, Scale, &s1);
    FCTN(y, Scale, &s2);
    
    
    if (s1 <= 0.0 || s2  <= 0.0)
      ERR1("'%.50s' must be a positive function", KNAME(DSCALE));
    //  if (x[0]!=1.0) printf("\n%.50s x=%10g,%10g %10g y=%10g,%10g %10g; %d\n", NAME(Scale), x[0], x[1], s1, y[0], y[1], s2, xdimown);
    smeanSq = 0.5 * (s1 * s1 + s2 * s2);
    //  printf("s_x[0]=%10g %10g s=%10g %10g %10g ", x[0], x[1], s1, s2, smeanSq);
    s = SQRT(smeanSq);
    for (i=0; i<xdimown; i++) {
      z1[i] = x[i] / s;
      z2[i] = y[i] / s;
    }
    

      //printf("h/s= %10g\n", SQRT((z1[0] - z2[0]) * (z1[0] - z2[0])  + (z1[1] - z2[1]) * (z1[1] - z2[1])));
     //  if (x[0]!=1.0) printf("s=%10g, %10g %10g; %10g %10g \n", s, z1[0], z1[1], z2[0], z2[1]);
    //    assert(z2[1] < 0.3);
	     
  } else if ((notdelete = aniso==NULL && (scale==NULL || scale[0] == 1.0))) {
    z1 = x;
    z2 = y;
  } else {
    int xdimown = OWNTOTALXDIM;
    double *xz1, *xz2;
    TALLOC_GLOBAL_X1(z1, xdimown);
    TALLOC_GLOBAL_X2(z2, xdimown);
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
    if (Sign != NULL) var = LOG(var);
  } else {
    assert(equalsCoordinateSystem(OWNISO(0)))
    double w;
    FCTN(y, cov->kappasub[DVAR], &w);
    FCTN(x, cov->kappasub[DVAR], &var);
    var *= w;
    var = Sign == NULL ?  SQRT(var)  : 0.5 * LOG(var);    
  }

  if (Scale != NULL) {
     int xdimown = OWNTOTALXDIM;
     double s12 = s1 * s2 /smeanSq,
       dimhalf = 0.5 * xdimown;
     if (Sign == NULL) 
       if (xdimown == 2) var *= s12; // deal most frequent case faster
       else var *= POW(s12, dimhalf);
     else var += dimhalf * LOG(s12); 
  }

  assert(Sign == NULL);
  
  if (Sign == NULL) {
    NONSTATCOV(z1, z2, next, v);
    for (i=0; i<vdimSq; i++) v[i] *= var;
  } else { // Sign != NULL
    LOGNONSTATCOV(z1, z2, next, v, Sign);
    for (i=0; i<vdimSq; i++) v[i] += var;
  }

  if (!notdelete) {
    FREE_TALLOC(z1);
    FREE_TALLOC(z2);
  }
  
}


void covmatrixS(model *cov, double *v) {
  //  location_type *loc = Loc(cov);	      
  model *next = cov->sub[DOLLAR_SUB];
  assert(next != NULL);
  location_type *locnext = Loc(next);
  assert(locnext != NULL);
  int i, tot, totSq,
    dim = GetLoctsdim(cov),
    vdim = VDIM0;
  KEY_type *KT = cov->base;
  assert(dim == PREVTOTALXDIM);

   SPRINTF(KT->error_loc, "'%.50s'", NICK(cov));

   if ((!PisNULL(DSCALE) && P0(DSCALE) != 1.0) || 
       !PisNULL(DANISO) || !PisNULL(DPROJ) || 
       cov->kappasub[DSCALE] != NULL ||
       cov->kappasub[DAUSER] != NULL || cov->kappasub[DANISO] != NULL ||
       cov->kappasub[DPROJ] != NULL ||
       !DefList[NEXTNR].is_covmatrix(next)
       ) {  
    StandardCovMatrix(cov, v); 
    return;
  }
 
  if (cov->Spgs == NULL && alloc_cov(cov, dim, vdim, vdim) != NOERROR)
    XERR(ERRORMEMORYALLOCATION);

  int L = OWNLASTSYSTEM; if (L != LASTSYSTEM(PREVSYSOF(next))) BUG;
  for (int s=0; s<=L; s++) 
    if (XDIM(PREVSYSOF(next), s) != NEXTXDIM(s)) {
      BUG; // fuehrt zum richtigen Resultat, sollte aber nicht
      // vorkommen!
      CovarianceMatrix(cov, v); 
      return;
    }


// save trafos
  Systems_type prev, gatter, own;
  COPYALLSYSTEMS(prev, PREVSYSOF(next), false);
  COPYALLSYSTEMS(gatter, GATTERSYSOF(next), false);
  COPYALLSYSTEMS(own, NEXT, false);

  // next trafos reset by trafos of cov
  COPYALLSYSTEMS(PREVSYSOF(next), PREV, false);
  COPYALLSYSTEMS(GATTERSYSOF(next), GATTER, false);
  COPYALLSYSTEMS(NEXT, OWN, true);
  DefList[NEXTNR].covmatrix(next, v);//hier wird uU next->totalpoints gesetzt

  // restore trafos
  COPYALLSYSTEMS(PREVSYSOF(next), prev, false);
  COPYALLSYSTEMS(GATTERSYSOF(next), gatter, false);
  COPYALLSYSTEMS(NEXT, own, false);

  tot = VDIM0 * locnext->totalpoints;
  totSq = tot * tot;
  
  if (cov->Sdollar->simplevar) {
    double var = P0(DVAR);
    if (var == 1.0) return;
    for (i=0; i<totSq; v[i++] *= var);
  } else {
    BUG;
  }
}

char iscovmatrixS(model *cov) {
  model *sub = cov->sub[DOLLAR_SUB];
  return (int) ((PisNULL(DSCALE) || P0(DSCALE) ==1.0) &&	
		PisNULL( DAUSER) && PisNULL( DANISO) &&
		PisNULL(DPROJ) &&
		cov->Sdollar->simplevar && // to do
		PisNULL(DANISO)) * DefList[SUBNR].is_covmatrix(sub);
}

void DS(double *x, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert( cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL && cov->kappasub[DSCALE] == NULL);
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj;
  double y[2], varSc,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;
  assert(cov->Sdollar->simplevar);
  assert(isCartesian(OWNISO(0)));

  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varSc = P0(DVAR) * spinvscale;

  if (nproj == 0) {
    y[0] = x[0] * spinvscale; 
    y[1] = (equalsIsotropic(OWNISO(0)) || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    assert(equalsSpaceIsotropic(OWNISO(0)));
    //  assert(P0 INT(DPROJ) == PROJ_SPACE); // #define Nproj
    y[0] = x[0] * spinvscale;
    y[1] = RF_NAN;
  }

  Abl1(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varSc; 
}

void DDS(double *x, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL );
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj,
    *proj = PPROJ;
  double y[2], varScSq,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;
  
  assert(isCartesian(OWNISO(0)));
  assert(cov->Sdollar->simplevar);
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScSq = P0(DVAR) * spinvscale * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (equalsIsotropic(OWNISO(0)) || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i] - 1] * x[proj[i] - 1];
    }
    y[0] = SQRT(y[0]) * spinvscale;
  }
  Abl2(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScSq; 
}


void D3S(double *x, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL);
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj,
    *proj = PPROJ;
  double y[2], varScS3,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;

  assert(isCartesian(OWNISO(0)));
  assert(cov->Sdollar->simplevar);
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScS3 = P0(DVAR) * spinvscale * spinvscale * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (OWNISO(0)==ISOTROPIC || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i] - 1] * x[proj[i] - 1];
    }
    y[0] = SQRT(y[0]) * spinvscale;
  }
  Abl3(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScS3; 
}

void D4S(double *x, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL);
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj,
    *proj = PPROJ;
  double y[2], varScS4,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;

  assert(isCartesian(OWNISO(0)));
  if (aniso != NULL) {
    spinvscale *= aniso[0]; 
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScS4 = spinvscale * spinvscale;
  varScS4 *= varScS4 * P0(DVAR);
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (equalsIsotropic(OWNISO(0)) || cov->ncol[DANISO]==1) ? 0.0
     : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i] - 1] * x[proj[i] - 1];
    }
    y[0] = SQRT(y[0]) * spinvscale;
  }
  Abl4(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScS4; 
}


void nablahessS(double *x, model *cov, double *v, bool nabla){
  model *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER] :
    cov->kappasub[DANISO];
  if (Aniso != NULL) BUG;
  int i, endfor,
    dim = cov->nrow[DANISO],// == ncol == x d i m ?
    xdimown = OWNTOTALXDIM,
    nproj = Nproj;
  double *xy, *vw,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    var = P0(DVAR);
  if (nproj != 0) BUG;
  if (dim != xdimown) BUG;
	   
  if (!cov->Sdollar->simplevar) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");


  TALLOC_DOUBLE(z0);
  TALLOC_DOUBLE(w0);
  if (aniso != NULL) {  
    TALLOC_GLOBAL_X1(z0, xdimown);
    TALLOC_GLOBAL_X2(w0, xdimown);
    xA(x, aniso, xdimown, xdimown, z0);
    xy = z0;
    vw = w0;
  } else {
    xy = x;
    vw = v;
  }

  TALLOC_DOUBLE(y0);
  if (scale != NULL) {
    TALLOC_GLOBAL_X3(y0, xdimown);
    assert(scale[0] > 0.0);
    double spinvscale = 1.0 / scale[0];
    var *= spinvscale;
    if (!nabla) var *= spinvscale; // gives spinvscale^2 in case of hess
    for (i=0; i<xdimown; i++) y0[i] = xy[i] * spinvscale;
    xy = y0;
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
    FREE_TALLOC(z0);
    FREE_TALLOC(w0);
  }

  FREE_TALLOC(y0);
  
  for (i=0; i<endfor; i++) v[i] *= var; 
}

void nablaS(double *x, model *cov, double *v){
  nablahessS(x, cov, v, true);
}
void hessS(double *x, model *cov, double *v){
  nablahessS(x, cov, v, false);
}


 

void inverseS(double *x, model *cov, double *v) {
  model *next = cov->sub[DOLLAR_SUB];
  int i,
    idx[3] = {DANISO, DAUSER, DPROJ};
  double 
    scale;
  
  if (cov->kappasub[DVAR] != NULL) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");

  for (i=0; i<3; i++) {
    if (cov->kappasub[idx[i]] != NULL) {
      ERR1( "inverse can only be calculated if '%.20s' is not an arbitrary function",
	    KNAME(idx[i])); 
    }
  }
  if (cov->kappasub[DSCALE] != NULL) {
    double left;
    model *sub = cov->kappasub[DSCALE];
    NONSTATINVERSE(ZERO(sub), sub, &left, &scale);
     if (left < 0.0) ERR("scale not defined to be non-negative.");
  } else scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
 
  int
    dim = PREVTOTALXDIM,
    nproj = Nproj;
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
  
  if (DefList[NEXTNR].inverse == ErrInverse) BUG;
  INVERSE(&y, next, v);
 
  for (i=0; i<dim; i++) v[i] *= s; //!

}


void nonstatinverseS(double *x, model *cov, double *left, double*right,
		     bool log){
  model
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE];
   int i,
    dim = PREVTOTALXDIM,
    nproj = Nproj;
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


  if (DefList[NEXTNR].nonstat_inverse == ErrInverseNonstat) BUG;
  if (log) {
    NONSTATLOGINVERSE(&y, next, left, right);
    // printf("S %f  %f %f %s\n", y, *left, *right, NAME(next));
    
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
	size = ncol * sizeof(double);
      bool redo;
      long bytes = ncnr * sizeof(double);
      if (ncol != nrow) BUG;
      redo = S->save_aniso == NULL;
      ALLC_NEW(Sdollar, save_aniso, ncnr, save_aniso);
      ALLC_NEW(Sdollar, inv_aniso, ncnr, inv_aniso);
      TALLOC_X1(LR, ncol);
      S->done = false;
      int pid;
      if (!redo) {
	for (i=0; i<ncnr;i++) if ((redo = save_aniso[i] != P(DANISO)[i])) break;
	Ext_pid(&pid);
	if (S->pid == 0) {
	  S->pid = pid;	  
	  if (S->pid == pid) { // two processes writing on the same place
	    S->busy = true; // one process doing the job
	  } else if (GLOBAL_UTILS->basic.cores == 1 || !S->busy)
	    ERR("Multiprocessor conflict. Plse inform maintainer & try again");
	} else if (GLOBAL_UTILS->basic.cores == 1 || !S->busy) 
	  ERR("Multiprocessor conflict. Plse inform maintainer & try again");
      }
      if (redo && S->pid == pid) {
	MEMCOPY(save_aniso, P(DANISO), bytes);
	MEMCOPY(inv_aniso, P(DANISO), bytes);
	if (Ext_invertMatrix(inv_aniso, nrow) != NOERROR) XERR(ERRORANISO_INV);
	S->busy = false;
	S->done = true;
      }
      while(S->busy && !S->done) {
	int sleep = (double) nrow * nrow * nrow / 1e8;
	sleep = sleep == 0 ? 1 : sleep;
	if (PL >= PL_BRANCHING) { PRINTF("|");}
	Ext_sleepMicro(&sleep);
      }
      
      MEMCOPY(LR, right, size);
      xA(LR, inv_aniso, nrow, ncol, right);
      
      MEMCOPY(LR, left, size);
      xA(LR, inv_aniso, nrow, ncol, left);      
      END_TALLOC_X1;
    }
  }

  if (Aniso != NULL) {
    if (aniso != NULL) BUG;
    if (DefList[MODELNR(Aniso)].inverse == ErrInverse) XERR(ERRORANISO_INV);
    int 
      nrow = Aniso->vdim[0],
      ncol = Aniso->vdim[1],
      size = nrow * sizeof(double);
    if (dim != ncol || ncol != 1)
      ERR("anisotropy function not of appropriate form");
    TALLOC_X1(LR, nrow);
    
    MEMCOPY(LR, right, size);
    INVERSE(LR, Aniso, right);
          
    MEMCOPY(LR, left, size);
    INVERSE(LR, Aniso, left);
    END_TALLOC_X1;
 }

  if (Scale != NULL && !isnowRandom(Scale)) {
    double dummy;
    COV(ZERO(Scale), Scale, &dummy);
    s *= dummy;
  } else if (scale != NULL) s *= scale[0];  
  if (s != 1.0) {
    for (i=0; i<dim; i++) {
      left[i] *= s; //!
      right[i] *= s;
    }
  }  
}

void nonstatinverseS(double *x, model *cov, double *left, double*right) {
  nonstatinverseS(x, cov, left, right, false);
}

void nonstat_loginverseS(double *x, model *cov, double *left, double*right){
 nonstatinverseS(x, cov, left, right, true);
}

void coinitS(model *cov, localinfotype *li) {
  assert(cov->Sdollar->simplevar);
  model *next = cov->sub[DOLLAR_SUB];
  if ( DefList[NEXTNR].coinit == NULL)
    ERR("# cannot find coinit -- please inform author");
  DefList[NEXTNR].coinit(next, li);
}
void ieinitS(model *cov, localinfotype *li) {
   assert(cov->Sdollar->simplevar);
 model *next = cov->sub[DOLLAR_SUB];
  
  if ( DefList[NEXTNR].ieinit == NULL)
    ERR("# cannot find ieinit -- please inform author");
  DefList[NEXTNR].ieinit(next, li);
}

void tbm2S(double *x, model *cov, double *v){
  assert(cov->Sdollar->simplevar);
 model *next = cov->sub[DOLLAR_SUB];
 // defn *C = DefList + NEXTNR; // kein gatternr, da isotrop
  double y[2],  *xy,
    *scale =P(DSCALE),
    *aniso = P(DANISO);
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL);
 
  assert(Nproj == 0);
  if (aniso!=NULL) {
    if (cov->ncol[DANISO]==2) {  // turning layers
      y[0] = x[0] * aniso[0]; // spatial 
      y[1] = x[1] * aniso[3]; // temporal
      assert(aniso[1] == 0.0 && aniso[2] == 0.0);
      if (y[0] == 0.0) COV(y, next, v) else TBM2CALL(y, next, v);
    } else {
      assert(cov->ncol[DANISO]==1);
      if (cov->nrow[DANISO] == 1) { // turning bands
	y[0] = x[0] * aniso[0]; // purely spatial 
 	TBM2CALL(y, next, v);
      } else { // turning layers, dimension reduction
	if (P0(DANISO) == 0.0) {
	  y[0] = x[1] * aniso[1]; // temporal 
	  COV(y, next, v); 
	} else {
	  y[0] = x[0] * aniso[0]; // temporal 
	  TBM2CALL(y, next, v);
	}
      } 
    }
    xy = y;
  } else xy = x;

   if (scale != NULL) {
    double s = scale[0];
    if (s > 0) { 
      double invscale = 1.0 / s;
      if (OWNTOTALXDIM == 2){
	y[0] = xy[0] * invscale; // spatial 
	y[1] = xy[1] * invscale; // temporal
	if (y[0] == 0.0) COV(y, next, v) else TBM2CALL(y, next, v);
      } else {
	y[0] = xy[0] * invscale; // purely spatial 
	TBM2CALL(y, next, v);
     }
   } else {
      y[0] = (s < 0.0 || xy[0] != 0.0) ? RF_INF : 0.0;
      if (OWNTOTALXDIM == 2)
 	y[1] = (s < 0.0 || xy[1] != 0.0) ? RF_INF : 0.0;
      TBM2CALL(y, next, v);
    }
   } else {
     TBM2CALL(xy, next, v);
   }
  *v *= P0(DVAR);
}


// TODO : Aniso=Matrix: direkte implementierung in S,
// sodass nicht ueber initS gegangen werden muss, sondern
//  e  < -  e^\top Aniso


int TaylorS(model *cov) {
  model 
    *next = cov->sub[DOLLAR_SUB],
    *sub = cov->key == NULL ? next : cov->key;
  int i;

  if (PisNULL(DPROJ) && PisNULL(DANISO)) {
    double scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
    cov->taylorN = sub->taylorN;  
    for (i=0; i<cov->taylorN; i++) {
      cov->taylor[i][TaylorPow] = sub->taylor[i][TaylorPow];
      cov->taylor[i][TaylorConst] = sub->taylor[i][TaylorConst] *
	P0(DVAR) * POW(scale, -sub->taylor[i][TaylorPow]);   
    }

    cov->tailN = sub->tailN;  
    for (i=0; i<cov->tailN; i++) {
      cov->tail[i][TaylorPow] = sub->tail[i][TaylorPow];
      cov->tail[i][TaylorExpPow] = sub->tail[i][TaylorExpPow];
      cov->tail[i][TaylorConst] = sub->tail[i][TaylorConst] *
	P0(DVAR) * POW(scale, -sub->tail[i][TaylorPow]);   
      cov->tail[i][TaylorExpConst] = sub->tail[i][TaylorExpConst] *
	POW(scale, -sub->tail[i][TaylorExpPow]);
    }
  } else {
    cov->taylorN = cov->tailN = 0;
  }
  RETURN_NOERROR;
}

int checkS(model *cov) {
  //  printf("checkS!!\n");

  // hier kommt unerwartet  ein scale == nan rein ?!!
  //  PMI(cov);
  model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE],
    *Var = cov->kappasub[DVAR],
    *sub = cov->key == NULL ? next : cov->key;
   isotropy_type owniso = OWNISO(0);
   int  err,
     //xdimgatter = GATTERTOTALXDIM,
     // KOMMENTAR NICHT LOESCHEN
     // *proj = PPROJ, // auf keinen Fall setzen, da Pointer unten neu
   //                        gesetzt wird!!!!
     last = OWNLASTSYSTEM; 
  // bool skipchecks = GLOBAL.general.skipchecks;
  matrix_type mtype = TypeMany;

  if (Aniso == NULL && Scale == NULL && Var == NULL &&
      PisNULL(DANISO) &&  PisNULL(DAUSER) &&  PisNULL(DSCALE) &&  PisNULL(DVAR))
  //  SERR("no parameter is set.");
 
  assert(isAnyDollar(cov));
  if (!isDollarProc(cov)) set_nr(OWN, DOLLAR); //wegen nr++ unten! NICHT SET_NR!
  
  // cov->q[1] not NULL then A has been given
  
  if ((err = checkkappas(cov, false)) != NOERROR) { // muss weit vorne stehen!
    RETURN_ERR(err);
  }
  kdefault(cov, DVAR, 1.0);
 //  printf("checkS2!!\n");

  if ( (!PisNULL(DSCALE) || Scale != NULL) && OWNLASTSYSTEM != 0)
    SERR1("'%.50s' may be used only if a single coordinate system is envolved", 
	  KNAME(DSCALE));

  bool angle = isAngle(Aniso);
  int dim = OWNLOGDIM(0);
  if (angle && PisNULL(DANISO) && PisNULL(DAUSER) &&
      OWNLASTSYSTEM == 0 && OWNXDIM(0) == dim) {
    ASSERT_CARTESIAN;
    if (!isCartesian(owniso)) RETURN_ERR(ERRORANISO);
    if ((err = CHECK(Aniso, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		     SUBMODEL_DEP, cov->frame)) == NOERROR) {
      if (isverysimple(Aniso)) {     
 	 PALLOC(DAUSER, dim, dim);
 	 AngleMatrix(Aniso, P(DAUSER));
	 COV_DELETE(cov->kappasub + DAUSER);
	 Aniso = cov->kappasub[DAUSER] = NULL;
	 angle = false;
      }
    }
  }

  if (cov->kappasub[DANISO] != NULL) // kann auch weggelassen werden. sollte aufgefangen werden
    SERR1("'%.50s' may not be a function.", KNAME(DANISO));
  if (!PisNULL(DAUSER)) {
    if (!isCartesian(owniso))  RETURN_ERR(ERRORANISO);
    //    if (GLOBAL.internal.warn_Aniso) {
    //      PRINTF("NOTE! Starting with RandomFields 3.0, the use of '%.50s' is different from\nthe former '%.50s' insofar that '%.50s' is multiplied from the right by 'x' (i.e. Ax),\nwhereas '%.50s' had been multiplied from the left by 'x' (i.e. xA).\n", KNAME(DAUSER), KNAME(DANISO), KNAME(DANISO), KNAME(DAUSER));
    //    }
    GLOBAL.internal.warn_Aniso = false;
    // here for the first time
    if (!PisNULL(DANISO)) RETURN_ERR(ERRORANISO_T); 
    int 
      lnrow = cov->nrow[DAUSER],
      lncol = cov->ncol[DAUSER];
    long
      total = lncol * lnrow;
	
    double
      *pA = P(DAUSER); 
    PALLOC(DANISO, lncol, lnrow); // !! ACHTUNG col, row gekreuzt
    for (int i=0, k=0; i<lnrow; i++) {
      for (long j=i; j<total; j+=lnrow) P(DANISO)[k++] = pA[j];
    }
    PFREE(DAUSER);
  }

  bool simplevar = Var == NULL || isnowRandom(Var); // checkkappa s.o.
  if (!simplevar) {
    ptwise_type ptt = cov->ptwise_definite;
    assert(equalsCoordinateSystem(owniso));
    if ((err = CHECK(Var, OWNLOGDIM(0), OWNXDIM(0), 
		     ShapeType, // only!! -- for pos def use RMprod
		     XONLY, owniso,
		     SCALAR, TrendType)) != NOERROR) RETURN_ERR(err);
    if (Var->ptwise_definite != pt_posdef) {
      if (Var->ptwise_definite == pt_unknown) {
	if (GLOBAL.internal.warn_negvar && cov->q==NULL) { // just a flag
	  QALLOC(1);
	  warning("positivity of the variance in '%.50s' cannot be detected.",
		  NICK(next));
	}
      } else {
	SERR2("positivity of '%.50s' required. Got '%.50s'", KNAME(DVAR), 
	      POSITIVITY_NAMES[Var->ptwise_definite]);
      }
    }
    cov->ptwise_definite = ptt;
  }

  ONCE_NEW_STORAGE(dollar);
  EXTRA_STORAGE;
  
  int xdimNeu = OWNXDIM(0);
 //  printf("checkS ff !!\n");

  if (Aniso != NULL) { // CASE 1: Aniso is given by arbitrary fctn
    if (!isDollarProc(cov) && !angle && !isKernel(OWN)) RETURN_ERR(ERRORFAILED);
    if (!PisNULL(DANISO) || !PisNULL(DPROJ) || !PisNULL(DSCALE) ||
	(Scale!=NULL && !isnowRandom(Scale)) )
      SERR2("if '%.50s' is an arbitrary function, only '%.50s' may be given as additional parameter.", KNAME(DAUSER), KNAME(DVAR));

    //assert(equalsCartCoord(owniso) || (angle && equalsSymmetric(owniso)));
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = wahr;
    if ((err = CHECK(Aniso, OWNLOGDIM(0), OWNXDIM(0), ShapeType, XONLY,
		     CARTESIAN_COORD, SUBMODEL_DEP, cov->frame)) != NOERROR) {
       RETURN_ERR(err);
    }

    if (Aniso->vdim[1] != 1)
      SERR3("'%.50s' returns a matrix of dimension %d x %d, but a vector is need.",
	    KNAME(DANISO), Aniso->vdim[0], Aniso->vdim[1]);

    if (cov->key==NULL) {
      ASSERT_QUASIONESYSTEM;
      if ((err = CHECK(sub, Aniso->vdim[0], Aniso->vdim[0], 
		      OWNTYPE(0), OWNDOM(0), owniso, 
		      SUBMODEL_DEP, cov->frame)) != NOERROR) {
	RETURN_ERR(err);
      }
    }
    cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = 
      cov->pref[Hyperplane] = cov->pref[SpectralTBM] = cov->pref[TBM] = 
      PREF_NONE;

    sub->pref[CircEmbed] = sub->pref[CircEmbedCutoff] = 
      sub->pref[CircEmbedIntrinsic] = sub->pref[Sequential] = 
      sub->pref[Specific] = PREF_NONE;
 
  } else if (Scale != NULL && !isnowRandom(Scale)) {
     // CASE 2: Scale is given by arbitrary function
    if (!isDollarProc(cov) && !isKernel(OWN)) RETURN_ERR(ERRORFAILED);
    if (!PisNULL(DANISO) || !PisNULL(DPROJ) || !PisNULL(DSCALE))
       SERR2("if '%.50s' is an arbitrary function, only '%.50s' may be given as additional parameter.", KNAME(DSCALE), KNAME(DVAR));
 
    assert(equalsCartCoord(owniso));
    
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = wahr;
    assert(equalsCartCoord(owniso));

    if ((err = CHECK(Scale, OWNLOGDIM(0), OWNXDIM(0), ShapeType, XONLY,
		     owniso, SCALAR, TrendType)) != NOERROR) {
      RETURN_ERR(err);
    }
    if (Scale->vdim[1] != 1 || Scale->vdim[0] != 1)
      SERR3("'%.50s' must be scalar, not %d x %d",
	    KNAME(DSCALE), Aniso->vdim[0], Aniso->vdim[1]);
    if (cov->key==NULL) {
      ASSERT_QUASIONESYSTEM;
      if ((err = CHECK(sub, OWNLOGDIM(0), OWNXDIM(0),
		       SYSTYPE(OWN, 0), OWNDOM(0), 
		       owniso, SUBMODEL_DEP, cov->frame)) != NOERROR) {
	RETURN_ERR(err);
      }
      if (!isNormalMixture(next)) 
	SERR("scale function only allowed for normal mixtures.");
    }
    
 
    //    APMI(cov);

    cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = 
      cov->pref[Hyperplane] = cov->pref[SpectralTBM] = cov->pref[TBM] = 
      PREF_NONE;

    sub->pref[CircEmbed] = sub->pref[CircEmbedCutoff] = 
      sub->pref[CircEmbedIntrinsic] = sub->pref[Sequential] = 
      sub->pref[Specific] = PREF_NONE;

  } else if (!PisNULL(DANISO)) { // CASE 3, anisotropy matrix is given
    int 
      nrow = cov->nrow[DANISO],
      ncol = cov->ncol[DANISO];

   
    if (nrow==0 || ncol==0) SERR("dimension of the matrix is 0");
    if (!PisNULL(DPROJ)) RETURN_ERR(ERRORANISO_T);

    int xdimown = OWNTOTALXDIM;
    if (xdimown < nrow) {
      if (PL >= PL_ERRORS) {LPRINT("xdim=%d != nrow=%d\n", xdimown, nrow);}
      SERR("#rows of anisotropy matrix does not match dim. of coordinates");
    }
    if (xdimown != OWNLOGDIM(0) && nrow != ncol)
      SERR("non-quadratic anisotropy matrices only allowed if dimension of coordinates equals spatio-temporal dimension");

    // PMI0(cov);
    // printf("nrow=%d %d\n", nrow, ncol);
    
    mtype = Type(P(DANISO), nrow, ncol);
    
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = wahr;

    if (last == 0 || isSpaceIsotropic(OWN)) {
      if (!equalsSpaceIsotropic(OWN)) {
	switch (owniso) {
	case ISOTROPIC :   
	  if (OWNXDIM(0) != 1) RETURN_ERR(ERRORANISO);
	  cov->full_derivs = cov->rese_derivs = 2;
	  break;
	case VECTORISOTROPIC :
	  if (!isMiso(mtype)) RETURN_ERR(ERRORANISO); 
	break;
	case SYMMETRIC: case CARTESIAN_COORD :
	  break;
	case PREVMODEL_I : BUG;      
	case GNOMONIC_PROJ :  case ORTHOGRAPHIC_PROJ:
	  if (!isnowProcess(cov)) RETURN_ERR(ERRORANISO);
	  break;
	default :
	  if (isCartesian(owniso)) { BUG; }
	  else RETURN_ERR(ERRORANISO);
	}
	
	if (!isMiso(mtype)) cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;
	ASSERT_ONESYSTEM;
	err = CHECK(sub, ncol, ncol, SYSTYPE(OWN, 0), OWNDOM(0), 
		    ncol==1 && !isnowProcess(cov)
		    ? IsotropicOf(owniso) : owniso, 
		    SUBMODEL_DEP, cov->frame);
      } else {	
	// spaceisotropic
	/*
    
	  case DOUBLEISOTROPIC :  
	  cov->full_derivs =  cov->rese_derivs = isMdiag(mtype) ? 2 : 0;
	  if (nrow != 2 || !isMdiag(mtype)) {
	  SERR("spaceisotropy needs a 2x2 diagonal matrix");
	  }
	  break;    	  
	*/
	
	cov->full_derivs = cov->rese_derivs = isMdiag(mtype) ? 2 : 0;
	if (nrow != 2 || !isMdiag(mtype)) 
	  SERR("spaceisotropy needs a 2x2 diagonal matrix");
	if (nrow != ncol) BUG;
	COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
	err = CHECK_GEN(sub, SUBMODEL_DEP, SUBMODEL_DEP, cov->frame, false);
      }

      if (err != NOERROR) RETURN_ERR(err);

    } else { // more than 1 coordinate system
      if (!isMdiag(mtype) ) 
	SERR1("If several coordinate systems are envolved, currently '%.50s' can only be a diagonal matrix", KNAME(DANISO));
      
      COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
      if ((err = CHECK_GEN(sub, SUBMODEL_DEP, SUBMODEL_DEP, cov->frame, false)) 
	  != NOERROR) {
	RETURN_ERR(err);
      }
      cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;

  /*

    if (!isMdiag(mtype)) 
      cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = cov->pref[Hyperplane] = PREF_NONE;
    if (owniso != DOUBLEISOTROPIC && !isMiso(mtype))
      cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;

      if ((err = CHECK(suib, ncol, ncol, cov->typus, OWNDOM(0), 
		 ncol==1 && !isProcess(cov->typus) ? ISOTROPIC : owniso, 
		 SUBMODEL_DEP, cov->frame))
	  != NOERROR) {
 	RETURN_ERR(err);
      }
    */
      
    }
 
    if (!isMdiag(mtype)) 
      cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] =
	cov->pref[Hyperplane] = PREF_NONE;
  
  } else if (!PisNULL(DPROJ)) { // geaendert werden: xdim, logdim, iso
    COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
    int xdim = OWNTOTALXDIM,
      p = PINT(DPROJ)[0];  // #define PPROJ

    FREE(PPROJ);
    if (p > 0) {
      Nproj = cov->nrow[DPROJ]; // #define Nproj
      int bytes = sizeof(int) * Nproj; 
      PPROJ = (int*) MALLOC(bytes);
      MEMCOPY(PPROJ, PINT(DPROJ), bytes); // #define PPROJ
    } else {
      int nproj = cov->nrow[DPROJ]; // #define Nproj
      if (nproj != 1)
	SERR1("values 'space' and 'time' in argument '%.50s' do not allow additional values", KNAME(DPROJ));
      if (!(Loc(cov)->Time && xdim >= 2 &&
	  (isCartesian(OWNISO(last))
#if MAXSYSTEMS == 1
	   || (isAnySpherical(PREVISO(0)) && PREVXDIM(0) > 2)
#endif							      
	   )))
	SERR1("unallowed use of '%.50s' or model too complicated.", KNAME(DPROJ));

      bool Space = p == PROJ_SPACE;
      assert(Space || p  == PROJ_TIME);
      if (Space) {
	Nproj = xdim - 1;
	PPROJ = (int*) MALLOC(sizeof(int) * Nproj);
	for (int i=1; i<=Nproj; i++) PPROJ[i-1] = i;
      } else { // Time
	Nproj = 1;
	PPROJ = (int*) MALLOC(sizeof(int) * Nproj);
	PPROJ[0] = xdim;
      }	
    }

    xdimNeu = Nproj;
    int max = xdim;
    if (xdimNeu == 0) SERR("Projection to null space not allowed.");

    if (xdimNeu < xdim) {	
      cov->pref[CircEmbed] = cov->pref[CircEmbedCutoff] =
	cov->pref[CircEmbedIntrinsic] =  cov->pref[Direct] =
	cov->pref[Sequential] = 1;
      cov->pref[TBM] = cov->pref[SpectralTBM] = cov->pref[Average]
	= cov->pref[Nugget]= cov->pref[RandomCoin] = cov->pref[Hyperplane]
	=  PREF_NONE;
      cov->pref[Specific] = PREF_BEST;
    }
     
    /// hier Aenderung ab version 3.2.2 so dass spaceisotropic
    // nicht mehr ueber 2 coordinatensysteme dargestellt werden kann
    int nproj = Nproj;
    for (int i=0; i<nproj; i++) {
      int idx = PPROJ[i];
      if (idx < 1) SERR1("only positive values allowed for '%.50s'", KNAME(DPROJ))
	else if (idx > max) SERR2("value of %.50s[%d] too large", KNAME(DPROJ), i);
	
      for (int j=i+1; j<nproj; j++) {
	int jdx = PPROJ[j];
	if (jdx == idx) SERR1("values of '%.50s' must be distinct.", KNAME(DPROJ));
      }
    }

    set_xdim(PREVSYSOF(sub), 0, xdimNeu);
    set_logdim(PREVSYSOF(sub), 0, xdimNeu);
    if (equalsSpaceIsotropic(OWN)) {
      // printf("%d %d\n", xdimNeu, nproj);
      if (xdimNeu > 2) SERR("maximum length of projection vector is 2");
      if (xdimNeu == 2) {
	if (PPROJ[0] >= PPROJ[1])
	  SERR1("in case of '%.50s' projection directions must be ordered",
		ISO_NAMES[DOUBLEISOTROPIC]);
      } else {
#if MAXSSYSTEMS == 1
	RETURN_ERR(ERRORWRONGISO);
#else	    
	assert(xdimNeu == 1);
	set_iso(PREVSYSOF(sub), 0, ISOTROPIC);
	if (PPROJ[0] == 1) { // space
	  set_logdim(PREVSYSOF(sub), 0, OWNLOGDIM(0) - 1);
	} else { // time
	  // PMI0(cov);
	  // printf("%d %d\n", PPROJ[0], OWNLOGDIM(0));
	  assert(PPROJ[0] == 2);
	  set_logdim(PREVSYSOF(sub), 0, 1);
	}
#endif	    
      }
    } else {// ! spaceisotropic
      switch (OWNISO(0)) {
      case ISOTROPIC : case EARTH_ISOTROPIC : case SPHERICAL_ISOTROPIC : 
	if (OWNXDIM(0) != 1) RETURN_ERR(ERRORANISO);
	break;
      case VECTORISOTROPIC :
	SERR("projection of vectorisotropic fields not programmed yet"); // to do: vdim muss auch reduziert werden ... --- das wird
	// grausam !
	break;
      case SYMMETRIC: case CARTESIAN_COORD:
	break;
      case GNOMONIC_PROJ :  case ORTHOGRAPHIC_PROJ :
	set_iso(PREVSYSOF(sub), 0, CARTESIAN_COORD);
	break;
      case PREVMODEL_I : BUG;      
	break;
      case SPHERICAL_SYMMETRIC : case EARTH_SYMMETRIC :
	if (nproj != 2 || PPROJ[0] != 1 || PPROJ[1]  != 2) {
	  for (int ii=0; ii<nproj; ii++)
	    if (PPROJ[ii] <= 2) APMI0(cov); // RETURN_ERR(ERRORANISO);
	  /// ehemals  owniso = SYMMETRIC; -- ueberall owniso
	  set_iso(PREVSYSOF(sub), 0, SYMMETRIC);
	}
	break;
      case SPHERICAL_COORD : case EARTH_COORD :
	if (nproj != 2 || PPROJ[0] != 1 || PPROJ[1]  != 2) {
	  for (int ii=0; ii<nproj; ii++) // ??
	    if (PPROJ[ii] <= 2) RETURN_ERR(ERRORANISO);
	  set_iso(PREVSYSOF(sub), 0, CARTESIAN_COORD);
	}
	break;
      default : 
	if (isCartesian(OWNISO(0))) {BUG;}
	else return  ERRORANISO;  // todo
      }
    }
    
      
    //    if ((err = CHECK(sub, nproj, nproj, OWNTYPE(0), OWNDOM(0), owniso,
    //		     SUBMODEL_DEP, cov->frame)) != NOERROR)  RETURN_ERR(err); 
    if ((err = CHECK_GEN(sub, 
			 VDIM0, // SUBMODEL_DEP; geaendert 20.7.14
			 VDIM1, cov->frame, true)) != NOERROR) {
      //APMI(cov->calling->calling);
      RETURN_ERR(err);
    }
   
    
  } else { // nproj == 0
 //  printf("checkS  u!!\n");
    COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);

    // verhindern, dass die coordinaten-transformation anlaeuft,
    // da aus z.B. Erd-coordinaten durch Projektion (auf die Zeitachse)
    // kartesische werden

    
    if ((err = CHECK_GEN(sub, 
			 VDIM0, // SUBMODEL_DEP; geaendert 20.7.14
			 VDIM1, cov->frame, true)) != NOERROR) {
      //   printf("\n\n\n\n\nRMS err = %d\n", err);      APMI(sub);
      RETURN_ERR(err);
    }
    
      
    /*   if (cov->key==NULL) {
	 if ((err = CHECK_NO_TRAFO(next, tsdim, xdim Neu, cov->typus, OWNDOM(0),
	 owniso, 
	 VDIM0, // SUBMODEL_DEP; geaendert 20.7.14
	 cov->frame)) != NOERROR) {
	 RETURN_ERR(err);
	 }

	 if (next->domown == OWNDOM(0) &&
	 next->owniso == owniso) // darf nicht nach C-HECK als allgemeine Regel ruebergezogen werden, u.a. nicht wegen stat/nicht-stat wechsel !!
	 // hier geht es, da domown und owniso nur durchgegeben werden und die Werte      // bereits ein Schritt weiter oben richtig/minimal gesetzt werden.
	 } else {
	 if ((err = CHECK_NO_TRAFO(cov->key, tsdim, xdimN eu, cov->typus,
	 OWNDOM(0), owniso,
	 SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
	 }

	 //     PMI(cov);
	 */

  } // end no aniso ((end of CASE 4))
//  printf("checkS 664!!\n");
 
  if (( err = checkkappas(cov, false)) != NOERROR) {
    RETURN_ERR(err);
  }
  
 //  printf("checkS 4!!\n");
  setbackward(cov, sub);
  for (int s=0; s<=last; s++) set_maxdim(OWN, s, OWNLOGDIM(s));

  if ((Aniso != NULL || (Scale != NULL && !isnowRandom(Scale)) || 
       !PisNULL(DANISO) || !PisNULL(DPROJ))) { 
    for (int s=0; s<=last; s++) {
      int max = MAXDIM(SUB, s);
      if (MAXDIM(OWN, s) < max) set_maxdim(OWN, s, max);
    }
  }

  if ((!isAnyIsotropic(owniso) && !isDollarProc(cov)) || last > 0){
    // multivariate kann auch xdimNeu == 1 problematisch sein
    set_nr(OWN, COVNR + 1); // JA NICHT SET_NR!!
  }
  // if (!isAnyIsotropic(owniso) && !isDollarProc(cov)) { // multivariate kann auch xdimNeu == 1 problematisch sein
  //   COVNR++;
  // }
  
  if (xdimNeu > 1) {
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 0;
  }

  // 30.10.11 kommentiert:
  //  cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 
  //      cov->pref[TBM] = cov->pref[SpectralTBM] = 0;
  if ( (PisNULL(DANISO) || isMiso(mtype)) && PisNULL(DPROJ)) {
    cov->logspeed = sub->logspeed * P0(DVAR);
  }
  //////////////// 

  if (sub->pref[Nugget] > PREF_NONE) cov->pref[Specific] = 100;
  if (PisNULL(DPROJ)) cov->matrix_indep_of_x = sub->matrix_indep_of_x;

  if ((err = TaylorS(cov)) != NOERROR) RETURN_ERR(err);  
  cov->Sdollar->simplevar = simplevar;
  ASSERT_ONESYSTEM;
  cov->Sdollar->orig_owniso = owniso; // for struct Sproc !

  if (isnowProcess(cov)) {
    MEMCOPY(cov->pref, PREF_NOTHING, sizeof(pref_shorttype)); 
  }

  if (GLOBAL.coords.coord_system == earth &&
      (PisNULL(DPROJ) ||
       (PPROJ[0] != PROJ_TIME && (Nproj != 1 || PPROJ[0] <= 2))) &&
      isCartesian(DEFSYS(next)) &&//is_all(isCartesian, DefList + NEXTNR) &&
      GLOBAL.internal.warn_scale &&
      (PisNULL(DSCALE) || 
       P0(DSCALE) < (STRCMP(GLOBAL.coords.newunits[0], "km")== 0 ? 10 : 6.3))) {
    GLOBAL.internal.warn_scale = false;
    PRINTF("Value of the scale parameter equals '%4.2f', which is less than 100,\nalthough models defined on R^3 are used in the context of earth coordinates,\nwhere larger scales are expected. (This message appears only ones per session.)\n", PisNULL(DSCALE) ? 1.0 : P0(DSCALE)); 
  }
 
 // printf("checkS! end!\n");

  //PMI0(cov);
  RETURN_NOERROR;
}


bool allowedDS(model *cov) {
  model 
    *Aniso =cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO], // nur cartesian, nur kernel (ausser angle)
    *Daniso = cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE], // && !isRandom(Scale) if (!isDollarProc(cov) && !isKernel(OWN)) RETURN_ERR(ERRORFAILED); NUR cartesian
    *Var = cov->kappasub[DVAR];
  bool angle = isAngle(Aniso) || isAngle(Daniso),
    *D = cov->allowedD;

  if ((Scale != NULL && !isRandom(Scale) && !isDollarProc(cov) ) ||
      (Aniso != NULL && !angle) || // angle erlaubt XONLY!
      (Var != NULL && ! isRandom(Var))) {
    assert(LAST_DOMAINUSER == 1);
    D[XONLY] = false;
    D[KERNEL] = true;
    return false;
  }
  return allowedDstandard(cov);
}


bool allowedIS(model *cov) {
  model 
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],  // nur cartesian, nur kernel (auch angle)
    *Daniso = cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE], // && !isRandom(Scale) if (!isDollarProc(cov) && !isKernel(OWN)) RETURN_ERR(ERRORFAILED); NUR cartesian
    *Var = cov->kappasub[DVAR];
  bool ScaleAniso = (Scale != NULL && !isRandom(Scale)) || Aniso != NULL ||
                  Daniso != NULL,
    //angle = isAngle(Aniso),
    var = Var != NULL && !isRandom(Var);
  bool *I = cov->allowedI;

  bool allowed = allowedIstandard(cov);
  if (allowed) for (int i=FIRST_ISOUSER; i<=LAST_ISOUSER; I[i++] = true);
  //
  //for (int i=0; i<LAST_ISOUSER; i++) printf("%d ", I[i]); BUG;
  
  //  assert(!angle);
  //  printf("%d %d %d\n", Scale != NULL && !isRandom(Scale), Aniso != NULL || Daniso != NULL, !angle);
  //  APMI(cov);

  int nproj = cov->nrow[DPROJ], // #define Nproj
    *proj = PINT(DPROJ);// #define PPROJ
  if (nproj>0) {
    allowed = false;
    assert(SYMMETRIC == 2 + DOUBLEISOTROPIC);
    
    bool hasprev = PREV_INITIALISED;
    isotropy_type previso = hasprev ? PREVISO(0) : ISO_MISMATCH;
    bool anyEarth = hasprev && (isEarth(previso) || isSpherical(previso)),
      time = proj[0] == PROJ_TIME || (nproj == 1 && proj[0] == OWNLOGDIM(0));
    if (!time && anyEarth) {      
      int min = proj[0];  
      for (int i=1; i<nproj; i++) if (proj[i] < min) min = proj[i];
      time = min > 2;
    }
    if (time) {
      cov->IallowedDone = false;
      if (anyEarth && !isAnyIsotropic(previso))
      	set_iso(PREV, 0,
		equalsAnySymmetric(previso) ? SYMMETRIC : CARTESIAN_COORD);
      if (hasprev && (!anyEarth || PREVLOGDIM(0)<=2)) {
	for (int i = 1 + (int) CARTESIAN_COORD; i <= (int) LAST_ISOUSER;
	     I[i++]=false);
      }
    }
  }

  if (ScaleAniso || var) {
     allowed = false;
   int i = (int) FIRST_ISOUSER,
      last = //angle ? SYMMETRIC :
      CARTESIAN_COORD;
    
    while (i <= last && !I[i]) i++;
    I[last] = i <= last;
    for ( ; i < (int) last; I[i++] = false);
    if (ScaleAniso) {
      for (i = 1 + (int)CARTESIAN_COORD; i <= (int) LAST_ISOUSER; I[i++]=false);
      return false; 
    } else {
      i = (int) FIRST_EARTH;
      while (i <= (int) EARTH_COORD && !I[i]) i++;
      I[EARTH_COORD] = i <= EARTH_COORD;
      for ( ; i < (int) EARTH_COORD; I[i++] = false);

      i = (int) FIRST_SPHERICAL;
      while (i <= (int) SPHERICAL_COORD && !I[i]) i++;
      I[SPHERICAL_COORD] = i <= SPHERICAL_COORD;
      for ( ; i < (int) SPHERICAL_COORD; I[i++] = false);
    }
  }

  else if ((!PisNULL(DANISO) /* internal */ && cov->nrow[DANISO] > 1) ||
	   (!PisNULL(DAUSER) && cov->ncol[DAUSER] > 1) || nproj > 0){
    allowed = false;
   int i = (int) FIRST_ISOUSER;
    while (i <= (int) SYMMETRIC && !I[i]) i++;
    I[SYMMETRIC] = i <= SYMMETRIC;
    for ( ; i < (int) SYMMETRIC; I[i++] = false);

    I[EARTH_SYMMETRIC] |= I[EARTH_ISOTROPIC];
    I[EARTH_ISOTROPIC] = false;
    I[SPHERICAL_SYMMETRIC] |= I[SPHERICAL_ISOTROPIC];		       
    I[SPHERICAL_ISOTROPIC] = false;

    if (GetTime(cov)) {      
      I[DOUBLEISOTROPIC] = (!PisNULL(DANISO) && cov->nrow[DANISO] == 2) ||
	(!PisNULL(DAUSER) && cov->ncol[DAUSER] == 2) ||
	nproj > 0; // nicht sehr praezise...
    }   
  }

  return allowed;
}


void rangeS(model *cov, range_type* range){
  int i;
  bool negdef = isnowNegDef(cov);
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
    range->pmin[i] = - 10000;
    range->pmax[i] = 10000;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
  
  range->min[DPROJ] = -2;
  range->max[DPROJ] = GATTERTOTALXDIM;
  range->pmin[DPROJ] = 1;
  range->pmax[DPROJ] =  range->max[DPROJ];
  range->openmin[DPROJ] = false;
  range->openmax[DPROJ] = false;
}


Types TypeS(Types required, model *cov, isotropy_type required_iso){
  bool ok =(COVNR != DOLLAR_PROC &&
	    (isShape(required) || isTrend(required) || equalsRandom(required)))
    || (COVNR == DOLLAR_PROC && isProcess(required));
  //  printf("ok=%d\n", ok);
  if (!ok) return BadType;
 
  model *sub = cov->key==NULL ? cov->sub[0] : cov->key;

  return TypeConsistency(required, sub, required_iso);
}


void spectralS(model *cov, gen_storage *s, double *e) {
  model *next = cov->sub[DOLLAR_SUB];
  int d,
    ncol = PisNULL(DANISO) ? OWNLOGDIM(0) : cov->ncol[DANISO];
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


void ScaleDollarToLoc(model *to, model *from,
		      int VARIABLE_IS_NOT_USED depth) {

  assert(!PARAMisNULL(to, LOC_SCALE));
  assert(isDollar(from));
  assert(!PARAMisNULL(from, DSCALE));
  PARAM(to, LOC_SCALE)[0] = PARAM0(from, DSCALE);
}

bool ScaleOnly(model *cov) {
  return isDollar(cov) && 
    PisNULL(DPROJ) && cov->kappasub[DPROJ] == NULL &&
    PisNULL(DAUSER) &&  cov->kappasub[DAUSER] == NULL &&
    PisNULL(DANISO) &&  cov->kappasub[DANISO] == NULL &&
    (PisNULL(DVAR) || P0(DVAR)==1.0) && cov->kappasub[DVAR] == NULL;
}



int addScales(model **newmodel, model *calling, model *Scale, double scale) {  
  if (scale != 1.0) { 
    addModel(newmodel, LOC, calling);
    kdefault(*newmodel, LOC_SCALE, scale);
  }
  if (Scale != NULL) { // could happen that scale != 1.0 && Scale != NULL
    //                  since scale includes also scale from aniso
    if (!isnowRandom(Scale)) RETURN_ERR_COV(Scale, XERRORNONSTATSCALE);
    addModel(newmodel, LOC, calling);
    addSetDistr(newmodel, Scale->calling, ScaleDollarToLoc, true, MAXINT);
  }
  return NOERROR;
}

int addPGSLocal(model **Key, // struct of shape
		  model *shape,// shape itself
		  model *local_pts,
		  int dim, int vdim, Types frame) {
  // SEE ALSO addPGS in extremes.cc
  /// random coin und smith

  // versucht die automatische Anpassung einer PointShapeType-Funktion;
  // derzeit wird 
  // * PGS und
  // * STANDARD_SHAPE (weiss nicht mehr wieso -> coins?)
  // versucht; inklusive Init
  //
  // ruft t.w.  FillInPts auf
  
#define specific (nPOISSON_SCATTER - 1)
  bool maxstable = hasMaxStableFrame(shape);
  int i,
    err = NOERROR,
    method = GLOBAL.extreme.scatter_method,
    pgs[specific] = {maxstable ? ZHOU : BALLANI, STANDARD_SHAPE}; 
  #define infoN 200
  char info[nPOISSON_SCATTER - 1][LENERRMSG];
  assert(shape != NULL);
  assert(*Key == NULL);

  // to do: strokorbball: raeumlich waere Standard/Uniform besser;
  // praeferenzen programmieren?
  for (i=0; i<specific; i++) {
    if (method != POISSON_SCATTER_ANY && method != i) continue;
    
     if (i > 0) {
      errorMSG(err, info[i-1]);
      //     XERR(err);  // eigentlich muss das hier weg
    }
    //    if (i > 0) XERR(err); assert(i ==0);
    if (*Key != NULL) COV_DELETE(Key);
    assert(shape->calling != NULL);
    addModel(Key, pgs[i], shape->calling);
 
    if ((err = FillInPts(*Key, shape)) != NOERROR) continue;
    if (MODELNR(*Key) != ZHOU) continue;
    model *local, *dummy,
      *cov = *Key;
    if ((err = covcpy(&local, false, local_pts, shape->prevloc, NULL, 
		      true, true, false)) != NOERROR) RETURN_ERR(err);
    SET_CALLING(local, (*Key)->calling);
    dummy = local;
    while (dummy->sub[DOLLAR_SUB] != NULL) dummy = dummy->sub[DOLLAR_SUB];
    if (MODELNR(dummy) != LOC) BUG;
    dummy->sub[DOLLAR_SUB] = *Key;
    SET_CALLING(*Key, dummy);
    
    SET_CALLING(cov, shape->calling);
    SET_CALLING(cov->sub[PGS_FCT], cov);
    SET_CALLING(cov->sub[PGS_LOC], cov);

    assert(cov->sub[PGS_LOC] != NULL && cov->sub[PGS_FCT] != NULL);
    cov->nsub = 2;

     if ((err = CHECK(*Key, dim, dim, PointShapeType, XONLY, 
		     CoordinateSystemOf(ISO(SYSOF(shape), 0)),
		     vdim, frame)) != NOERROR) {
      continue; 
    }
    NEW_COV_STORAGE(cov, gen);

    if ((err = INIT(cov, 1, cov->Sgen)) == NOERROR) break;
  } // for i_pgs

  model *cov = *Key;
  if (err != NOERROR) {
     SERR("error occured when creating the local point-shape fctn");
   // errorstring_type xx = "error occured when creating a local point-shape fctn";
    // errorMSG(err, xx, cov->base, info[i-1]);
    // SERR2("no method found to simulate the given model:\n\tpgs: %.50s\n\tstandard: %.50s)", info[0], info[1]);
  }
  RETURN_NOERROR;
}



int structS(model *cov, model **newmodel) {
  if (isnowProcess(cov)) {
    assert(hasGaussMethodFrame(cov));
    SET_NR(cov, DOLLAR_PROC);
    return structSproc(cov, newmodel); // kein S-TRUCT(...) !!
  } 

  model *local = NULL,
    *dummy = NULL,
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE];
  int err = NOERROR; 
  bool generalAniso = (Aniso != NULL && !isAngle(Aniso)) ||
    (Scale != NULL && !isnowRandom(Scale));

  if (generalAniso)
    GERR1("complicated models including arbitrary functions for '%.50s' cannot be simulated yet", KNAME(DAUSER));
   
  ASSERT_NEWMODEL_NOT_NULL;
 
  if (cov->kappasub[DVAR] != NULL) {
    GERR2("Arbitrary functions for '%.50s' should be replaced by multiplicative models using '%.50s'", KNAME(DVAR), DefList[PROD].nick);
  } 

  switch (cov->frame) {
  case PoissonGaussType :
    BUG;
    break;
  case SmithType : {
    ASSERT_ONESYSTEM;
    ASSERT_NEWMODEL_NOT_NULL; // ?
    if (next->randomkappa) GERR("random shapes not programmed yet");
    if (!PisNULL(DANISO)){err=ERRORANISO; goto ErrorHandling;}
    if (!PisNULL(DANISO)) {
      if (!isMonotone(next) || !isIsotropicXonly(next) ||
	  PisNULL(DANISO) || cov->ncol[DANISO] != cov->nrow[DANISO])
	GERR("anisotropy parameter only allowed for simple models up to now");
    }
    
    if (!PisNULL(DPROJ)) GERR("projections not allowed yet");
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);
   
    assert(*newmodel != NULL);
   
    double anisoScale = 1.0;  
    if (!PisNULL(DANISO)) {
      anisoScale = 1.0 / getMinimalAbsEigenValue(P(DANISO), cov->nrow[DANISO]);
      if (!R_FINITE(anisoScale)) {err=ERRORANISO_INV; goto ErrorHandling;}
      if (anisoScale < 0) {
	err = -anisoScale;
	goto ErrorHandling;
      }
    }

      
    Types type = BadType;
    double scale = PisNULL(DSCALE) ? 1.0 :P0(DSCALE);
    int logdim = OWNLOGDIM(0),
      xdim = OWNXDIM(0);
    // in same cases $ must be set in front of a pointshape function
    // and becomes a pointshape function too, see addpointshapelocal here
     if (isPointShape(*newmodel)) type = PointShapeType;
    else if ((err = CHECK_R(*newmodel, logdim)) == NOERROR) {
      type = RandomType;
      if ((err=addScales(newmodel, *newmodel, Scale, scale * anisoScale))
	  !=NOERROR) goto ErrorHandling;
      SET_CALLING(*newmodel, cov);      
    } 

    if (!equalsRandom(type)) {
      type = type == BadType ? ShapeType : type;
      if ((err = CHECK(*newmodel, logdim, xdim, type, OWNDOM(0), OWNISO(0),
		       cov->vdim, cov->frame)) != NOERROR)
	goto ErrorHandling;

      // crash();
      BUG;
      // ??   addModel(newmodel, DOLLAR); // 2.2.19
      // ??assert( (*newmodel)->calling == cov);

      // ??  if (!PisNULL(DVAR)) kdefault(*newmodel, DVAR, P0(DVAR));

      // ?? double scale = 1.0;
      //    if (!PisNULL(DSCALE)) {
      //      kdefault(*newmodel, DSCALE, P0(DSCALE));
      //      scale = P0(DSCALE);}
 
      if (type == PointShapeType && 
	  (err = addScales((*newmodel)->sub + PGS_LOC, *newmodel,
			   Scale, anisoScale *
			     scale)) != NOERROR) goto ErrorHandling;
      if ((err = CHECK(*newmodel, OWNLOGDIM(0), PREVXDIM(0), type, 
		       PREVDOM(0), PREVISO(0), cov->vdim, 
		       cov->frame))
	  != NOERROR) goto ErrorHandling;
      if (type==PointShapeType) {	
	if ((err = FillInPts(*newmodel, cov)) != NOERROR) goto ErrorHandling;
      } else {
	dummy = *newmodel;
	*newmodel = NULL;
	// suche nach geeigneten locationen
	if ((err = addScales(&local, dummy, Scale, scale * anisoScale))
	    != NOERROR) goto ErrorHandling;
	if ((err = addPGSLocal(newmodel, dummy, local,
			       OWNLOGDIM(0), VDIM0, cov->frame))
	    != NOERROR) goto ErrorHandling; 
      }
    } // ! randomtype
  }
  break;
  
  case SchlatherType: case BrMethodType:
    if (next->randomkappa) GERR("random shapes not programmed yet");
    if (!PisNULL(DPROJ)) GERR("projections not allowed yet");
    // P(DVAR) hat keine Auswirkungen
    if (!PisNULL(DANISO)) {
      if (!isMonotone(next) || !isIsotropicXonly(next) ||
	  PisNULL(DANISO) || cov->ncol[DANISO] != cov->nrow[DANISO])
	GERR("anisotropy parameter only allowed for simple models up to now");
    }
    
    assert(cov->calling != NULL); 
    
    BUG;
    // 
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);   
    break;
    
  case GaussMethodType :
    if (cov->key != NULL) COV_DELETE(&(cov->key));

    if (PrevLoc(cov)->distances) 
      GERR("distances do not allow for more sophisticated simulation methods");
    
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);

    assert(*newmodel != NULL);
    addModel(newmodel, DOLLAR);
    if (!PisNULL(DVAR)) kdefault(*newmodel, DVAR, P0(DVAR));
    if (!PisNULL(DSCALE) ) kdefault(*newmodel, DSCALE, P0(DSCALE));
    if (next->randomkappa) GERR("random shapes not programmed yet");
    
    break;
 
  default :
    BUG;
    GERR2("%.50s : changes in scale/variance not programmed yet for '%.50s'", 
	  NICK(cov), TYPE_NAMES[cov->frame]);      
  }
  
 ErrorHandling:
  
  if (dummy != NULL) COV_DELETE(&dummy);
  if (local != NULL) COV_DELETE(&local);

   
  RETURN_ERR(err);
}




int initS(model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!
  model *next = cov->sub[DOLLAR_SUB],
    *Var = cov->kappasub[DVAR],
    *Scale = cov->kappasub[DSCALE],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *projM = cov->kappasub[DPROJ];
  int 
    vdim = VDIM0,
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    err = NOERROR;
  bool 
    angle = isAngle(Aniso),
    smith = hasSmithFrame(cov);// Realisationsweise 

  if (smith || hasAnyPoissonFrame(cov)) {//Average !!ohne smith selbst!!
    ASSERT_ONESYSTEM;
    double
      var[MAXMPPVDIM],
      scale = PisNULL(DSCALE)  ? 1.0 : P0(DSCALE);
    int dim = OWNLOGDIM(0);
 
    if (!PisNULL(DPROJ) || !PisNULL(DANISO) || 
	projM!= NULL || (Aniso != NULL && (!angle || Aniso->randomkappa))){

      SERR("(complicated) anisotropy ond projection not allowed yet in Poisson related models");
    }
  
     // Achtung I-NIT_RANDOM ueberschreibt mpp.* !!
    if (Var != NULL) {
      int nm_neu = nm == 0 && !smith ? 1 : nm;
      if ((err = INIT_RANDOM(Var, nm_neu, s, P(DVAR))) != NOERROR) RETURN_ERR(err); 
      
      int nmP1 = Var->mpp.moments + 1;
      for (int i=0; i<vdim; i++) {
	int idx = i * nmP1;
	var[i] = smith ? P0(DVAR) : Var->mpp.mM[idx + 1];      
      }
    } else for (int i=0; i<vdim; var[i++] = P0(DVAR));

    if (Scale != NULL) {
      if (dim + nm < 1) SERR("found dimension <= 0");
      int dim_neu = smith ? nm : (dim + nm) < 1 ? 1 : dim + nm; 
      if ((err = INIT_RANDOM(Scale, dim_neu, s, P(DSCALE)))
	  != NOERROR) RETURN_ERR(err);
      scale = smith ? P0(DSCALE) : Scale->mpp.mM[1];      
    }
    if ((err = INIT(next, nm, s)) != NOERROR) RETURN_ERR(err);


    for (int i=0; i < nmvdim; i++) {
      cov->mpp.mM[i] = next->mpp.mM[i]; 
      cov->mpp.mMplus[i] = next->mpp.mMplus[i]; 
    }

    if (Var != NULL && !smith) {
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= Var->mpp.mM[i];
	cov->mpp.mMplus[i] *= Var->mpp.mMplus[i];
      }
    } else {
      int j, k;
      double pow_var;
      for (k=j=0; j<vdim; j++) { 
	pow_var = 1.0;
	for (int i=0; i <= nm; i++, pow_var *= var[j], k++) {	
	  cov->mpp.mM[k] *= pow_var;
	  cov->mpp.mMplus[k] *= pow_var;
	}
      }
    }
    if (Scale != NULL && !smith) {
      if (Scale->mpp.moments < dim) SERR("moments can not be calculated.");
      int j, k,
	nmP1 = Scale->mpp.moments + 1;
      for (k=j=0; j<vdim; j++) { 
	double pow_scale = Scale->mpp.mM[dim + j * nmP1];
	for (int i=0; i <= nm; i++, k++) {
	  cov->mpp.mM[k] *= pow_scale;
	  cov->mpp.mMplus[k] *= pow_scale;
	}
      }
    } else {
      double pow_scale = POW(scale, dim);
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= pow_scale;
	cov->mpp.mMplus[i] *= pow_scale;
      }
    }

    if (!PisNULL(DANISO)) {      
      if (cov->nrow[DANISO] != cov->ncol[DANISO]) RETURN_ERR(ERRORANISO_SQUARE);
      double invdet = FABS(1.0 / getDet(P(DANISO), cov->nrow[DANISO]));
      if (!R_FINITE(invdet)) RETURN_ERR(ERRORANISO_DET);
      
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= invdet;
	cov->mpp.mMplus[i] *= invdet;
      }

    }

    if (Aniso != NULL) {  
      double invdet;
      if (angle) {
	int 
	  ncol = Aniso->vdim[0],
	  nrow = Aniso->vdim[1];
	if (nrow != ncol) RETURN_ERR(ERRORANISO_SQUARE);
	double 
	  *diag = PARAM(Aniso, ANGLE_DIAG);
	if (diag != NULL) {
	  invdet = 1.0;
	  for (int i=0; i<ncol; i++) invdet /= diag[i];
	} else {
	  invdet = PARAM0(Aniso, ANGLE_RATIO);
	}
      } else {
	SERR("only anisotropy matrices basesd on RMangle allowed.");
      }
      
      invdet = FABS(invdet);      
      if (!R_FINITE(invdet)) RETURN_ERR(ERRORANISO_DET);
      
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= invdet;
	cov->mpp.mMplus[i] *= invdet;
      }
    }

    bool unnormed = R_FINITE(next->mpp.unnormedmass);
    if (unnormed) {
      //printf("unnormedmass in $ %.50s\n", NAME(next));
      if (vdim > 1) BUG;
      cov->mpp.unnormedmass = next->mpp.unnormedmass * var[0];
    } else cov->mpp.unnormedmass = RF_NAN;
    
    if (isnowRandom(next)) {
      if (unnormed) for (int i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_INF;
      else RETURN_ERR(ERRORNOTPROGRAMMEDYET);
    } else {
      for (int i=0; i<vdim; i++)
	cov->mpp.maxheights[i] = next->mpp.maxheights[i] * var[i];
    }
  } // hasSmithFrame

  else if (hasGaussMethodFrame(cov)) {
    model 
      *key = cov->key,
      *sub = key == NULL ? next : key;
    assert(sub != NULL);
    assert(key == NULL || ({PMI(cov);false;}));//
   
    if ((err=INIT(sub, 0, s)) != NOERROR) RETURN_ERR(err);
   
  }

  else if (hasAnyEvaluationFrame(cov)) {
    model 
      *key = cov->key,
      *sub = key == NULL ? next : key;
    assert(sub != NULL);
    assert(key == NULL || ({PMI(cov);false;}));//
   
    if ((err=INIT(sub, 0, s)) != NOERROR) RETURN_ERR(err);
    
  } else { 
    if ((err=INIT(next, 0, s)) != NOERROR) RETURN_ERR(err); // e.g. from MLE
  }

  if ((err = TaylorS(cov)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;

}


void doS(model *cov, gen_storage *s){
  model
    *Var = cov->kappasub[DVAR],
    *Scale = cov->kappasub[DSCALE];
  int vdim = VDIM0;
   
  if (Var != NULL) {
    if (isnowRandom(Var)) {
      if (isProcess(Var)) BUG;
      assert(!PisNULL(DVAR));
      DORANDOM(Var, P(DVAR));

      
    } else if (Var->randomkappa) {
      assert(!PisNULL(DVAR));
      DO(Var, s);
    } else BUG;
  }
  
 if (Scale != NULL) {
    if (isnowRandom(Scale)) {
      if (isProcess(Scale)) BUG;
      assert(!PisNULL(DSCALE));
      DORANDOM(Scale, P(DSCALE));
    } else if (Scale->randomkappa) {
      BUG;
      assert(!PisNULL(DSCALE));
      DO(Scale, s);
    } else BUG;
  }

  
  if (hasSmithFrame(cov) || hasAnyPoissonFrame(cov)) {
    model *next = cov->sub[DOLLAR_SUB];
    
    DO(next, s);// nicht gatternr
    for (int i=0; i<vdim; i++)
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * P0(DVAR);
    return;
  }
  
  else if (hasGaussMethodFrame(cov)) {    
    double 
      *res = cov->rf,
      sd = SQRT(P0(DVAR));
    int 
      totalpoints = Gettotalpoints(cov);
    assert(res != NULL);
    if (cov->key == NULL) BUG;

    if (sd != 1.0) for (int i=0; i<totalpoints; i++) res[i] *= (double) sd;
    
    return;
  } 

  BUG;
}

//////////////////////////////////////////////////////////////////////
// PROCESSES
//////////////////////////////////////////////////////////////////////

int structSproc(model *cov, model **newmodel) {
  model
    *next = cov->sub[DOLLAR_SUB],
    *Scale = cov->kappasub[DSCALE],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO];
  int
    dim = GetLoctsdim(cov), 
    newdim = dim,
    err = NOERROR; 
  // model *sub;
  assert(isDollarProc(cov));

  assert(newdim  > 0);
  
  if ((Aniso != NULL && Aniso->randomkappa) ||
      (Scale != NULL && Scale->randomkappa)
      ) {
    SERR1("complicated models including arbitrary functions for '%.50s' cannot be simulated yet", KNAME(DAUSER));
  }
  switch (cov->frame) {
  case GaussMethodType : {
    location_type *loc = Loc(cov); // do not move before TransformLoc!!
    ASSERT_NEWMODEL_NULL;
    if (cov->key != NULL) COV_DELETE(&(cov->key));

    if (PrevLoc(cov)->distances) 
      SERR("distances do not allow for more sophisticated simulation methods");
    
    if (Aniso!= NULL) { // Aniso fctn
      TransformLoc(cov, false, True, true);
      loc = Loc(cov); 
      assert(!loc->grid && !loc->Time);
      newdim = Aniso->vdim[0];
      if (newdim != dim) 
	ERR("change of dimension in struct S not programmed yet");
      
      int bytes = newdim * sizeof(double);
      long total = loc->totalpoints;
      double *v = NULL,
	*x = loc->x,
	*xnew = x; // not correct for newdim != dim;
      //              instead partial_loc_set should be used to reset loc
      assert(x != NULL);
      assert(!loc->grid);
      
      if ((v = (double*) MALLOC(bytes)) == NULL) RETURN_ERR(ERRORMEMORYALLOCATION);
      for (long i=0; i<total; i++, x+=dim, xnew += newdim) {
	FCTN(x, Aniso, v);
	MEMCOPY(xnew, v, bytes);
      }
      FREE(v);
    } else if (Scale != NULL) {
      // check code wether it will still work
      if (isnowRandom(Scale)) {
	SERR("no specific simulation technique available for random scale");
      }
      SERR("no specific simulation technique available for arbitrary scale");
    } else {      
      double scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
      cov->Sdollar->separable = !PisNULL(DPROJ) &&
	(loc->grid || (loc->Time && PINT(DPROJ)[0] < 0)); // #define PPROJ

	
      if (cov->Sdollar->separable) {
	int *proj = PPROJ,
	  nproj = Nproj;
	if (loc->grid || PPROJ[0] == dim) {
	  double *x = (double*) MALLOC(sizeof(double) * dim * 3);
	  for (int i=0; i<nproj; i++) {
	    MEMCOPY(x + i * 3,
		    proj[i] == dim && loc->Time ? loc->T : loc->xgr[proj[i] -1],
		    3 * sizeof(double));
	    x[i * 3 + 1] /= scale;
	  }
	  loc_set(x, NULL, NULL, nproj, nproj, 3, 0, false, true, false, cov);
	  FREE(x);
	} else { // space
	  int i = 0;
	  for ( ; i<nproj; i++) if (PPROJ[i] != i + 1) break;
	  if (i >= nproj) 
	    loc_set(loc->x, NULL, NULL,
		    nproj, nproj,
		    loc->spatialtotalpoints, 0,
		    false, loc->grid, false,
		    cov);
	  else {
	    int pts = Getspatialpoints(cov);
	    long int bytes = sizeof(double) * nproj * pts;
	    double *x = (double*) MALLOC(bytes),
	      *xx = x,
	      *L = loc->x;
	    for (i=0 ; i<pts; i++, L += dim) {	    
	      for (int j=0 ; j<nproj; j++, xx++) *(xx++) = L[proj[j]];
	    }
	    
	    loc_set(x, NULL, NULL,
		    nproj, nproj,
		    loc->spatialtotalpoints, 0,
		    false, loc->grid, false,
		    cov);
	    FREE(x);
	  }
	}
	newdim = GetLoctsdim(cov);
	assert(newdim > 0 );
      } else {
	// falls Aniso-matrix vom projektionstyp und zeitseparabel und
	// zeitkomponente fehlt
	usr_bool gridexpand = NEXTNR == TBM_PROC_INTERN || isAnyNugget(NEXTNR)
	  ? False : GRIDEXPAND_AVOID;
	// LocReduce: ausschiesslich Projektionen werden reduziert
	//	TransformLocReduce(cov, true /* timesep*/, gridexpand, 
	//		   true /* involveddollar */);
	TransformLoc(cov, true /* timesep*/, gridexpand, 
			   true /* involveddollar */);
	newdim = GetLoctsdim(cov);
   	assert(newdim > 0 );
      }
      
    }
 
    if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
    if (!isGaussMethod(cov->key)) addModel(&(cov->key), GAUSSPROC);
    SetLoc2NewLoc(cov->key, PLoc(cov));
    
    model *key;
    key = cov->key;
    assert(key->calling == cov);
    assert(newdim > 0);
    
    ASSERT_ONESYSTEM;
    if ((err = CHECK_NO_TRAFO(key, newdim, newdim, ProcessType, XONLY,
			      CoordinateSystemOf(cov->Sdollar->orig_owniso),
			      VDIM0, GaussMethodType)) != NOERROR) {
      RETURN_ERR(err);
    }
    
    err = STRUCT(key, NULL);

    RETURN_ERR(err);
  }
  default :
    SERR2("%.50s: changes in scale/variance not programmed yet for '%.50s'", 
	  NICK(cov), TYPE_NAMES[cov->frame]);      
  }
     
  RETURN_NOERROR;
}




int initSproc(model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!

  //model *next = cov->sub[DOLLAR_SUB];
  //mppinfotype *info = &(s->mppinfo);
  //  location_type *loc = cov->prevloc;
  model 
    *key = cov->key;
    //*sub = key == NULL ? next : key;
  location_type *prevloc = PrevLoc(cov);

  int 
    prevdim = prevloc->timespacedim,
    //dim = GetLoctsdim(cov),
    err = NOERROR;

  int //zaehler,
    prevtotalpts = prevloc->totalpoints,
    owntotalpts =  Loc(cov)->totalpoints;

      assert(key != NULL);
  
  if ((err = INIT(key, 0, s)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  key->simu.active = true; 
  assert(s != NULL);
  cov->fieldreturn = wahr;
  cov->origrf = prevtotalpts != owntotalpts;
  if (cov->origrf) {
    assert(cov->Sdollar->separable);
    assert(prevtotalpts % owntotalpts == 0);
    assert(!PisNULL(DPROJ) || !PisNULL(DANISO));
    // projection or reducing anisotropy matrix
    if (VDIM0 != VDIM1) BUG;
    cov->rf = (double*) MALLOC(sizeof(double) * VDIM0 * prevloc->totalpoints);
    COND_NEW_STORAGE(dollar, cumsum);
    EXTRA_STORAGE;
 
    int *proj = PPROJ;

    ALLC_NEWINT(Sdollar, cumsum, prevdim + 1, cumsum);
    ALLC_NEWINT(Sdollar, total, prevdim + 1, total);
    ALLC_NEWINT(Sdollar, len, prevdim + 1, len);
    
    if (cov->Sdollar->separable) {
      for (int d=0; d<prevdim; d++) {
	cumsum[d] = 0;
	len[d] = prevloc->xgr[d][XLENGTH];
      }
      if (proj != NULL) {
	int d=0,
	  nproj = Nproj;
	assert(proj[d] > 0);
	cumsum[proj[d] - 1] = 1;
	for (d = 1; d < nproj; d++) {
	  cumsum[proj[d] - 1] = cumsum[proj[d - 1] - 1] * len[d-1];
	}
      } else { // A eine projektionsmatrix
	int i,
	  iold = 0,
	  nrow = cov->nrow[DANISO],
	  ncol = cov->ncol[DANISO];
	double *A = P(DANISO);
	for (int d=0; d<ncol; d++, A += nrow) {
	  for (i = 0; i < nrow && A[i] == 0.0; i++);
	  if (i == nrow) i = nrow - 1;
	  if (d > 0) {
	    cumsum[i] = cumsum[iold] * len[d-1];
	  } else { // d ==0
	    cumsum[i] = 1;
	  }
	  iold = i;
	  for (i++; i < nrow; i++) if (A[i] != 0.0) BUG;  // just a check
	}
      }
      
    } else BUG; /* { // !prevloc->grid spaeter fuer matrizen
		//  die proj-character haben
      if (!prevloc->Time) goto Standard;
      len[0] = prevloc->spatialtotalpoints;
      len[1] = prevloc->T[XLENGTH];
      int nproj = Nproj;
      if (proj[0] != prevdim) { // spatial
	for (d=1; d<nproj; d++) if (proj[d] == prevdim) goto Standard;
	cumsum[0] = 1;
	cumsum[1] = 0;
      } else {
	if (nproj != 1) goto Standard;
	cumsum[0] = 0;
	cumsum[1] = 1;   
      }
      prevdim = 2;
      } */
    
    for (int d=0; d<prevdim; d++) total[d] = cumsum[d] * len[d];
  } else cov->rf = cov->key->rf;

  model *Var = cov->kappasub[DVAR];
  if (Var != NULL) {
    if (isnowRandom(Var) || Var->randomkappa) {
      if (isProcess(Var)) {
	RETURN_ERR(ERRORNOTPROGRAMMEDYET);
      } else {
	if ((err = INIT(Var, 0, s)) != NOERROR)
	  RETURN_ERR(err);	  
      }
    } else {
      assert(cov->Sdollar != NULL);
      int totptsvdim = prevloc->totalpoints * VDIM0;
      cov->Sdollar->sd = (double*) MALLOC(sizeof(double) * totptsvdim);
      double *sd = cov->Sdollar->sd;
      Fctn(NULL, cov, sd);
	for (int i=0; i<totptsvdim; i++) sd[i] = SQRT(sd[i]);
    }
  }
  
  model *Scale = cov->kappasub[DSCALE];
  if (Scale != NULL) {
    RETURN_ERR(ERRORNOTPROGRAMMEDYET);
  }
 
  RETURN_NOERROR;
}


//int zz = 0;
void doSproc(model *cov, gen_storage *s){
  int i, 
    vdim = VDIM0;

  if (hasGaussMethodFrame(cov)) {    
    assert(cov->key != NULL);
    double *res = cov->key->rf;
    int 
      totalpoints = Gettotalpoints(cov),
      totptsvdim = totalpoints * vdim;
    
    DO(cov->key, s); 
    
    model *Var = cov->kappasub[DVAR];

    // printf("her %d\n", zz);
    // if (zz++ == 1) crash();
    
    if (Var != NULL) {
      if (isnowRandom(Var) || Var->randomkappa) {
	if (isProcess(Var)) {
	  XERR(ERRORNOTPROGRAMMEDYET);
	} else {
	  DORANDOM(Var, P(DVAR));
 	}
	double sd = SQRT(P0(DVAR));
	for (i=0; i<totptsvdim; i++) res[i] *= sd;
 	
      } else {
	double *sd = cov->Sdollar->sd;
	assert(sd != NULL);
	for (i=0; i<totptsvdim; i++) res[i] *= sd[i];
      }
    } else {
      assert(!PisNULL(DVAR));
      double sd = SQRT(P0(DVAR)); 
      if (sd != 1.0) {
	for (i=0; i<totptsvdim; i++) res[i] *= sd;
      }
    }
  }
  
  else if (hasMaxStableFrame(cov) || hasAnyPoissonFrame(cov)) {
    BUG; // 2.2.19: darf eigentlich nicht (so) genutzt werden
    assert(vdim == 1);
    model *next = cov->sub[DOLLAR_SUB],
      *Var = cov->kappasub[DVAR],
      *Scale = cov->kappasub[DSCALE];
 
    assert(VDIM0 == VDIM1);

    if (Var != NULL && Var->randomkappa) {
      BUG; // Da muss irgenduas dann auf rf draufmultipliziert werden
      assert(!PisNULL(DVAR) && isnowRandom(Var));
      VTLG_R(NULL, Var, P(DVAR));
    }

    if (Scale != NULL && Scale->randomkappa) {
      // remote could be deterministic although local ist not
      BUG; // otherwise do must depend on the new scale oder man muss
      // interpolieren o.ae.
      assert(!PisNULL(DSCALE));
      VTLG_R(NULL, Scale, P(DSCALE));
    }
    
    DO(next, s);// nicht gatternr

    for (i=0; i<vdim; i++)   
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * P0(DVAR);
  }

  else {
    //PMI(cov);
    BUG;
  }
 
  if (cov->origrf) { // es wird in cov->key->rf nur 
    // der projezierte Teil simuliert. Dieser Teil
    // muss auf das gesamte cov->rf hochgezogen werden
    //if (vdim != 1) BUG;
    assert(cov->Sdollar->separable)
    location_type *prevloc = PrevLoc(cov);
    int zaehler, d,
      // 2 below: one for arbitrary space and one for timesollte 
      prevdim = prevloc->timespacedim,
      dim = prevloc->grid ? prevdim : 2,
      prevtotalpts = prevloc->totalpoints,
      owntotalpts =  Loc(cov)->totalpoints,
      *cumsum = cov->Sdollar->cumsum,
      *len = cov->Sdollar->len,
      *total = cov->Sdollar->total;
    
    if (cov->Sdollar->separable) {
      assert(cov->key != NULL);
      assert(total != NULL && cumsum != NULL);
      TALLOC_L1(nx, prevdim);
      
      for (d=0; d<dim; d++) {
	nx[d] = 0;
      }
      zaehler = 0;
      i = 0;
      
      //   PMI0(cov->calling);   PMI(cov);
      
      for (int v=0; v<vdim; v++) {
	double *res = cov->rf + v * prevtotalpts,
	  *rf = cov->key->rf + v * owntotalpts;
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
      
      END_TALLOC_L1;
    }
  }
}

