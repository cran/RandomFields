/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Definition of auxiliary correlation functions 

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nonsta     tionary models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2011 Martin Schlather

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
#include <assert.h>
#include "RF.h"
#include "Covariance.h"
//#include <R_ext/Lapack.h>
//#include <R_ext/Applic.h>
//#include <R_ext/Utils.h>     
//#include <R_ext/BLAS.h> 


void kdefault(cov_model *cov, int i, double v) {
  cov_fct *C = CovList + cov->nr;
  if (cov->p[i] == NULL) {
    if (C->kappatype[i]==REALSXP) {
      cov->p[i] = (double*) malloc(sizeof(double));
      cov->p[i][0] = v;
    } else if (C->kappatype[i]==INTSXP) {
      cov->p[i] = (double*) malloc(sizeof(int));
      ((int *) (cov->p[i]))[0] = (int) v;
    } else if (C->kappatype[i] == LISTOF + REALSXP) {
      assert(false);
    } else assert(false);
    cov->nrow[i] = cov->ncol[i] = 1;
  } else {
    if (cov->nrow[i] != 1 || cov->ncol[i] != 1) { 
	PRINTF("%d %s %d nrow=%d, ncol=%d\n ", 
	       cov->nr, CovList[cov->nr].name, i, cov->nrow[i], cov->ncol[i]);
	int j; for (j=0; j<4; j++) PRINTF("%f\n", cov->p[i][j]);
      char param_name[100]; 
      strcpy(param_name, CovList[cov->nr].kappanames[i]); 
      PERR("parameter not scalar");
    }
  }
}


void updatepref(cov_model *cov, cov_model *sub) {
  int i;
//  printf("update pref %s\n", CovList[cov->nr].name);
  for (i=0; i<Forbidden; i++) {
//      printf("%d;%d ", cov->pref[i], sub->pref[i]);
      if (sub->pref[i] < cov->pref[i]) {
	  cov->pref[i] = sub->pref[i];
      }
  }
  // printf("\n");
}


void setbackward(cov_model *cov, cov_model *sub) {
  // see also check_co when changing
  assert(cov != NULL);
  assert(sub != NULL);

  if (sub->maxdim < cov->maxdim) cov->maxdim = sub->maxdim;
  cov->normalmix &= sub->normalmix;
  cov->finiterange &= sub->finiterange;
  cov->diag &= sub->diag;
//  cov->quasidiag &= sub->quasidiag;
  cov->separatelast &= sub->separatelast;
  cov->semiseparatelast &= sub->semiseparatelast;  
  if (sub->derivatives < cov->derivatives) 
      cov->derivatives = sub->derivatives;
  
  updatepref(cov, sub);
  cov->tbm2num |= sub->tbm2num;

  cov->hess = (CovList[cov->nr].hess != NULL && sub->hess);
}

int checkkappas(cov_model *cov, bool errornull){
  cov_fct *C = CovList + cov->nr;

  int i,nr, nc,
    kappas= C->kappas,
    *ncol = cov->ncol,
    *nrow = cov->nrow;
  char param_name[PARAMMAXCHAR]; // used in PERR

  for (i=0; i<kappas; i++) {
    strcpy(param_name, C->kappanames[i]);
    if (cov->p[i] == NULL) {
      if (errornull) {
	PERR("parameter unset");
      } else continue;
    }

    C->kappasize(i, cov, &nr, &nc);
    
    if (nc < 0 || nr < 0) PERR("tried to access non-existing parameters ");
    
    if (nc == 1 && ncol[i] != nc) {
      if (nr == 1 && nrow[i] != 1)
	PERR("parameter must be a scalar");
      PERR("parameter must be a vector, not a matrix");
    }
    if (nc > 1 && ncol[i] == 1)
      PERR("parameter must be a matrix");
    
    if ((nc > 0 && ncol[i] != nc) || (nr > 0 && nrow[i] != nr)) {
      // nc==0, nr==0 is coded as "unknown"
      char msg[255];
      sprintf(msg, 
	      "parameter not of the given size: (%d, %d) instead of (%d, %d)", 
	      nrow[i], ncol[i], nr, nc);
      PERR(msg);
    }
    // if nc==0 (nr==0) then this is undetermined.
  }
  return NOERROR;
}

int checkkappas(cov_model *cov){
  return checkkappas(cov, true);
}



// $
void kappaS(int i, cov_model *cov, int *nr, int *nc){
   *nc = (i==DALEFT) ? cov->xdim 
       : (int) (i==DVAR || i==DSCALE); // 0 mean not determined
   *nr = (i == DALEFT) ? 0 : 
       (i == DANISO) ? cov->xdim :
       i <= DMAX ? 1 
       : -1; 
  // nicht cov->tsdim !!
}

  // simple transformations as variance, scale, anisotropy matrix, etc.  
void Siso(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim;
  double 
      var =  cov->p[DVAR][0],
      y = (cov->p[DANISO] == NULL ? 
	   fabs(*x / cov->p[DSCALE][0]) :  fabs(*x * cov->p[DANISO][0]));
  // letzteres darf nur passieren wenn dim = 1!!
  CovList[next->nr].cov(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= var; 

  // printf("Siso %f %f %f\n", x[0], v[0]);

}

void Sstat(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double 
      **p = cov->p,
      var = p[DVAR][0],
      *aniso=p[DANISO];
  int i,
       nproj = cov->nrow[DPROJ],
      vdimSq = cov->vdim * cov->vdim;
 
  if (nproj > 0) {
    int *proj = (int *)cov->p[DPROJ];
    double *z = (double*) malloc(nproj * sizeof(double)),
	invscale = 1.0 / p[DSCALE][0];
    for (i=0; i<nproj; i++) z[i] = invscale * x[proj[i]];
    CovList[next->nr].cov(z, next, v);
    free(z);
  } else { 
    if (aniso==NULL && p[DSCALE][0] == 1.0) {
      CovList[next->nr].cov(x, next, v);
    } else {
      int xdim = cov->xdim;
      double *z = (double*) malloc(xdim * sizeof(double)); 
      if (aniso==NULL) {
        double invscale = 1.0 / p[DSCALE][0];
	for (i=0; i < xdim; i++) z[i] = invscale * x[i];
      } else {
        int j, k,
	    nrow=cov->nrow[DANISO], 
	    ncol=cov->ncol[DANISO];
	for (k=i=0; i<ncol; i++) {
	    z[i] = 0.0;
	    for (j=0; j<nrow; j++) {
		z[i] += aniso[k++] * x[j];
	    }
	}
      }
      CovList[next->nr].cov(z, next, v);    
      free(z);
    }
  }
  for (i=0; i<vdimSq; i++) v[i] *= var; 
}

void Snonstat(double *x, double *y, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double 
      **p = cov->p, 
      var = p[DVAR][0],
      *aniso=p[DANISO];
  int i,
      nproj = cov->nrow[DPROJ],
      vdimSq = cov->vdim * cov->vdim;
  
  if (nproj > 0) {
    int *proj = (int *)cov->p[DPROJ];
    double 
	*z1 = (double*) malloc(nproj * sizeof(double)),
	*z2 = (double*) malloc(nproj * sizeof(double)),
	invscale = 1.0 / p[DSCALE][0];
    for (i=0; i<nproj; i++) {
	z1[i] = invscale * x[proj[i]];
	z2[i] = invscale * y[proj[i]];	
    }
    CovList[next->nr].nonstat_cov(z1, z2, next, v);
    free(z1);
    free(z2);
  } else { 
    if (aniso==NULL && p[DSCALE][0] == 1.0) {
	CovList[next->nr].nonstat_cov(x, y, next, v);
    } else {
      int xdim = cov->xdim;
      double 
	  *z1 = (double*) malloc(xdim * sizeof(double)),
	  *z2 = (double*) malloc(xdim * sizeof(double));
       if (aniso==NULL) {
	  double invscale = 1.0 / p[DSCALE][0];
	  for (i=0; i<xdim; i++) {
	      z1[i] = invscale * x[i];
	      z2[i] = invscale * y[i];
	  }
      } else {
        int j, k,
	    nrow=cov->nrow[DANISO], 
	    ncol=cov->ncol[DANISO];
	for (k=i=0; i<ncol; i++) {
	    z1[i] =  z2[i] =0.0;
	    for (j=0; j<nrow; j++) {
		z1[i] += aniso[k] * x[j];
		z2[i] += aniso[k] * y[j];
	    }
	}
      }
      CovList[next->nr].nonstat_cov(z1, z2, next, v);      
      free(z1);
      free(z2);
    }
  }
  for (i=0; i<vdimSq; i++) v[i] *= var; 
}


void DS(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim,
      nproj = cov->nrow[DPROJ];
   double y[2], 
    **p = cov->p,
    spinvscale = 
       cov->p[DANISO] == NULL ? 1.0 / cov->p[DSCALE][0] : cov->p[DANISO][0], 
    varSc = p[DVAR][0] * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale; 
    y[1] = (cov->isoIn==ISOTROPIC || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * p[DANISO][3]; // temporal; temporal scale
  } else {
      ERR("error in DS; please contact maintainer"); 
      // for (i=0; i<nproj; i++) {
      //   y[i] = spinvscale * x[proj[i]];
  }

  CovList[next->nr].D(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varSc; 
}

void DDS(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim,
    nproj = cov->nrow[DPROJ],
    *proj = (int *)cov->p[DPROJ];
  double y[2], 
    **p = cov->p,
    spinvscale = 
      cov->p[DANISO] == NULL ? 1.0 / cov->p[DSCALE][0] : cov->p[DANISO][0],
    varScSq = p[DVAR][0] * spinvscale * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (cov->isoIn==ISOTROPIC || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * p[DANISO][3]; // temporal; temporal scale
  } else {
    ERR("error in DDS; please contact maintainer"); 
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i]] * x[proj[i]];
    }
    y[0] = sqrt(y[0]) * spinvscale;
  }
  CovList[next->nr].D2(y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScSq; 
}

void nablaS(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim,
      dim = cov->nrow[DANISO],// == ncol == xdim 
      nproj = cov->nrow[DPROJ];
  double *y,  
      **p = cov->p;
  if (nproj != 0)  ERR("error in DS; please contact maintainer");
  y = (double*) malloc(sizeof(double) * dim);
	  
  if (cov->p[DANISO] == NULL) {      
      double spinvscale = 1.0 / cov->p[DSCALE][0],
	  varSc = p[DVAR][0] * spinvscale;
      for (i=0; i<dim; i++) y[i] = x[i] * spinvscale;
      CovList[next->nr].nabla(y, next, v);
      for (i=0; i<vdimSq; i++) v[i] *= varSc; 
  } else {
      double *z,
	  var = p[DVAR][0];
      z = (double*) malloc(sizeof(double) * dim);
      xA(x, p[DANISO], dim, dim, y);
      CovList[next->nr].nabla(y, next, z);
      Ax(p[DANISO], y, dim, dim, v);
      for (i=0; i<dim; i++) v[i] *= var; 
      free(z);
  }
  free(y);
}




void hessS(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim,
      nproj = cov->nrow[DPROJ],
      dim=cov->nrow[DANISO],
      dimsq = dim * dim; // == ncol
  double *y,
      **p = cov->p;
  if (nproj != 0)  ERR("error in DS; please contact maintainer");
  y = (double*) malloc(sizeof(double) * dim);
     
  if (cov->p[DANISO] == NULL) {      
      double spinvscale = 1.0 / cov->p[DSCALE][0],
	  varSc = p[DVAR][0] * spinvscale;
      for (i=0; i<dim; i++) y[i] = x[i] * spinvscale;
      CovList[next->nr].nabla(y, next, v);
      for (i=0; i<vdimSq; i++) v[i] *= varSc; 
  } else {
      double *z,
	  var = p[DVAR][0];
      z = (double*) malloc(sizeof(double) * dimsq);
      xA(x, p[DANISO], dim, dim, y);
      CovList[next->nr].nabla(y, next, z);
      XCXt(p[DANISO], z, v, dim, dim);
      for (i=0; i<dim; i++) v[i] *= var; 
      free(z);
  }
  free(y);
}


void invS(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim,
      nproj = cov->nrow[DPROJ];
//    *proj = (int *)cov->p[DPROJ];
  double y, 
    **p = cov->p,
    invscale = 
      cov->p[DANISO] == NULL ? 1.0 / cov->p[DSCALE][0] : cov->p[DANISO][0], 
    var = p[DVAR][0];
  
  if (nproj == 0) {
    y= *x / var; // inversion, so variance becomes scale
  } else {
    ERR("nproj is not allowed in invS");
  }

  CovList[next->nr].D(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] /= invscale; //!

}


void coinitS(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0]->sub[0];
  if ( CovList[next->nr].coinit == NULL)
    error("# cannot find coinit -- please inform author");
  CovList[next->nr].coinit(next, li);
}
void ieinitS(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0]->sub[0];
 
  if ( CovList[next->nr].ieinit == NULL)
    error("# cannot find ieinit -- please inform author");
  CovList[next->nr].ieinit(next, li);
}
void tbm2S(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  double y[2], 
      **p=cov->p, 
      *aniso = p[DANISO];
  
  assert(cov->nrow[DPROJ] == 0);
  if (aniso==NULL) {
    double invscale = 1.0 / p[DSCALE][0];
    if (cov->xdim == 2){
      y[0] = x[0] * invscale; // spatial 
      y[1] = x[1] * invscale; // temporal
      if (y[0] == 0.0) C->cov(y, next, v); 
      else C->tbm2(y, next, v);
    } else {
      y[0] = x[0] * invscale; // purely spatial 
      C->tbm2(y, next, v);
     }
  } else {
    if (cov->ncol[DANISO]==2) {  // turning layers
      y[0] = x[0] * aniso[0]; // spatial 
      y[1] = x[1] * aniso[3]; // temporal
      assert(aniso[1] == 0.0 && aniso[2] == 0.0);
      if (y[0] == 0.0) C->cov(y, next, v); 
      else C->tbm2(y, next, v);
    } else {
      assert(cov->ncol[DANISO]==1);
      if (cov->nrow[0] == 1) { // turning bands
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
  }
  *v *= p[DVAR][0];
}

int checkS(cov_model *cov) {

// hier kommt unerwartet  ein scale == nan rein ?!!

  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + cov->nr;
  double  **p = cov->p;
  int i, err,
      xdim = cov->xdim,
      nproj = cov->nrow[DPROJ];
  bool skipchecks = GLOBAL.general.skipchecks;

  strcpy(ERROR_LOC, C->name);
  assert(cov->nr >= DOLLAR && cov->nr<=LASTDOLLAR);
  cov->nr = DOLLAR; // wegen nr++ unten !
  cov->manipulating_x = !true; // !! exception as treated directly
  
  // cov->q[1] not NULL then A has been given

  if (cov->q == NULL && cov->p[DALEFT]!=NULL) {
    // here for the first time
    if (p[DANISO] != NULL) ERR("aniso and A may not be given at the same time");
    int j, k,
	nrow = cov->nrow[DALEFT],
        ncol = cov->ncol[DALEFT],
	total = ncol * nrow;
    double
	*pA = p[DALEFT], 
	*pa = p[DANISO];
    p[DANISO] = (double*) malloc(sizeof(double) * total);
    cov->nrow[DANISO] = ncol;
    cov->ncol[DANISO] = nrow;
    for (i=k=0; i<nrow; i++) {
	for (j=i; j<total; j+=nrow) p[DANISO][k++] = pA[j];
    }
    cov->q = (double*) calloc(1, sizeof(double));
    cov->qlen = 1;
  }
 
  kdefault(cov, DVAR, 1.0);
  if (p[DANISO] != NULL) { // aniso given
      int *idx=NULL,
	  nrow = cov->nrow[DANISO],
	  ncol = cov->ncol[DANISO];
    bool quasidiag;
    matrix_type type;

    idx = (int *) malloc((nrow > ncol ? nrow : ncol) * sizeof(int));
    if (nrow==0 || ncol==0) ERR("no dimension of the matrix may be 0");
    if (nproj > 0 || cov->ncol[DPROJ] > 0)
      ERR("the parameters aniso and proj may not be given at the same time");

    analyse_matrix(p[DANISO], nrow, ncol, 
		   &(cov->diag),
		   &quasidiag, // &(cov->quasidiag), 
		   idx, // cov->idx
		   &(cov->semiseparatelast),
		   &(cov->separatelast));
    free(idx); idx=NULL;
    type = Type(p[DANISO], nrow, ncol);
    
    
    cov->derivatives = (cov->xdim == 1 ||  nrow == ncol)  ? 2 : 0;

  
    
    if (xdim < nrow) {
//      || ( xdim > nrow && 
//	  (cov->tsdim != nrow || type!=TypeIso || xdim != 1))) 
      if (PL > 6) PRINTF("xdim=%d != nrow=%d\n", xdim, nrow);
      strcpy(ERRORSTRING,
	      "#rows of anisotropy matrix does not match dim. of coordinates");
      return ERRORM;
    }

    if (cov->isoIn == SPACEISOTROPIC) {
      cov->derivatives = 2;
      if (type != TypeDiag && type != TypeIso)  {
	strcpy(ERRORSTRING, "spaceisotropic needs diagonal matrix");
	return ERRORM + 1;
      }
    }

    if (p[DSCALE] != NULL) {
      strcpy(ERRORSTRING, 
	     "aniso and scale may not be given at the same time");
      return ERRORM + 2;
    }
    if (p[DPROJ] != NULL) {
	strcpy(ERRORSTRING, 
	     "aniso and proj may not be given at the same time");
      return ERRORM + 2;
    }

    if (cov->xdim != cov->tsdim && nrow != ncol) {
      strcpy(ERRORSTRING, "non-quadratic anisotropy matrices only allowed if dimension of coordinates equals spatio-temporal dimension");
      return ERRORM + 3;
    }
    
    if (!cov->diag)  
      cov->pref[RandomCoin] = cov->pref[Hyperplane]=0;//=cov->pref[SpectralTBM]

    if ((err = check2X(next, ncol, ncol, cov->statIn, 
		       ncol==1 ? ISOTROPIC : cov->isoIn, 
		       SUBMODELM))
	!= NOERROR) return err;

    if (next->sub[0]->statIn == cov->statIn &&
	next->sub[0]->isoIn == cov->isoIn &&
	cov->xdim > 1)
	next->tsdim = DEL_COV - 7;

 
    // kdefault(cov, DSCALE, 1.0);
  } else { // aniso == NULL
//      int xdim2 = xdim * xdim,
//      step = xdim + 1;

    //   if (cov->q[1] < 0.0 && (cov->tsdim != xdim || nproj > 0)) {
    //   strcpy(ERRORSTRING, "aniso cannot be forced");
    ///   return ERRORM + 4;
    // }

    int tsdim = cov->tsdim;
    if (nproj > 0) {
	bool *used = (bool*) malloc(xdim * sizeof(bool));
      if (cov->ncol[DPROJ] != 1) ERR("proj must be a vector");
      for (i=0; i<xdim; i++) used[i] = false;
      for (i=0; i<xdim; i++) {
	int idx = ((int *) cov->p[DPROJ])[0];
	if (used[idx] ) ERR("components in proj are given twice");
	used[idx] = true;
      }
      xdim = nproj;
      tsdim = nproj;
      free(used);
    }

    kdefault(cov, DSCALE, 1.0);
     
    //  if (p[DANISO] != NULL) free(p[DANISO]);
    //  p[DANISO] = NULL;
    
    /* 
    p[DANISO] = (double*) calloc(xdim2, sizeof(double));
    ncol = nrow = cov->ncol[DANISO] = cov->nrow[DANISO] = xdim;
    double ani = 1.0 / p[DSCALE][0];
    for (i=0; i<xdim2; i+=step) p[DANISO][i] = ani;
    */

    if ((err = check2X(next, tsdim, xdim, cov->statIn, cov->isoIn,
		       SUBMODELM))
	!= NOERROR) return err;

    if (next->sub[0]->statIn == cov->statIn &&
	next->sub[0]->isoIn == cov->isoIn)
      next->tsdim = DEL_COV - 8;

  }

  // if ((p[DANISO] != NULL && p[DANISO][0] >= 0 && xdim==1) next->tsdim = DEL_COV;
  cov->vdim = next->vdim;
  setbackward(cov, next);
  err = checkkappas(cov, false); // !! unchecked on NULL !!

  if (cov->p[DANISO] != NULL) {
    cov->pref[Nugget] = cov->pref[TBM3] = cov->pref[TBM2] = 0;
  }
  if (err != NOERROR) return err;
	    
  if (cov->xdim > 1) {
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 0;
    cov->nr++;
  }
  
  cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 
      cov->pref[TBM3] = cov->pref[TBM2] = cov->pref[SpectralTBM] = 0;

  return NOERROR;
}

int initspectralS(cov_model *cov) {
  cov_model *next = cov->sub[0];
  //  printf("$ ------> %s %d\n", CovList[next->nr].name, 5);
  if (CovList[next->nr].initspectral == NULL) return ERRORFAILED;
  memcpy(&(next->spec), &(cov->spec), sizeof(spec_covstorage));
  int err=CovList[next->nr].initspectral(next);
  memcpy(&(cov->spec), &(next->spec), sizeof(spec_covstorage));
//  printf("$ ------> %s %d\n", CovList[next->nr].name, err);
  return err;
}


void spectralS(cov_model *cov, spectral_storage *s, double *e) {
    int d;
  cov_model *next = cov->sub[0];
  double *sube;

  if (cov->p[DSCALE] != NULL) {
      int dim = cov->tsdim;
      double 
	  invscale = 1.0 / cov->p[DSCALE][0];
      sube = (double*) malloc(dim * sizeof(double));
      CovList[next->nr].spectral(next, s, sube);
      for (d=0; d<dim; d++)  e[d] = sube[d] * invscale;
  } else {
      int j, k, 
	  nrow = cov->xdim,
	  ncol = next->xdim;
      double 
	  *A = cov->p[DANISO];
      sube = (double*) malloc(ncol * sizeof(double));
      CovList[next->nr].spectral(next, s, sube);
      for (d=0, k=0; d<nrow; d++, k+=nrow) {
	  e[d] = 0.0;
	  for (j=0; j<ncol; j++) {
	      e[d] += sube[j] * A[j + k];
	  }
      }
  }
  free(sube);
}


void rangeS(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;
  int i;

  range->finiterange = true;
  range->maxdim = INFDIM;

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
  range->openmax[DSCALE] = false;

  for (i=DANISO; i<= DALEFT; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -10000;
    range->pmax[i] = 10000;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}


int initS(method_type *meth){
    // am liebsten wuerde ich hier die Koordinaten transformieren;
    // zu grosser Nachteil ist dass GetDiameter nach trafo 
    // grid def nicht mehr ausnutzen kann -- umgehbar?! 
  cov_model *cov=meth->cov;
  globalparam *gp = meth->gp;
  int err, 
      PL = gp->general.printlevel;
      //   ncol = cov->ncol[DANISO],
      //   nrow = cov->nrow[DANISO];
  if (meth->hanging != NULL) {
    if (PL>4) PRINTF("hanging %d\n", meth->hanging->nr);
    if (meth->hanging->nr != ISO2ISO &&	
	meth->hanging->nr != SP2SP &&
	meth->hanging->nr != S2S
      ) {      
      if (PL>5) {
	PRINTF("initS:\n");      
	PrintMethodInfo(meth);
      }
      return ERRORHANGING;
    }
  }

  // setAniso(meth); // schon im S-Teil ist die neue c-aniso gesetzt --
  //                 S-Teil darf somit nicht mit c-aniso aufgerufen werden
  //                 siehe doS, das sofort naechsten Teil aufruft.

  meth->sub[0] = (method_type*) malloc(sizeof(method_type));
  METHOD_NULL(meth->sub[0]);
  meth->nsub = 1;
  method_type *sub = meth->sub[0];
  cpyMethod(meth, sub, false);
  sub->cvar *= cov->p[DVAR][0];
//  int nproj = cov->nrow[DPROJ],
//      mnrow = meth->loc->timespacedim,
//      mncol = nrow;
//  matrix_type type_method;

  cum_dollar(meth, meth->loc->timespacedim, cov, sub);

  sub->cov = cov->sub[0];

//  PrintMethodInfo(sub); assert(false);

  err = initstandard(sub);

  meth->compatible = sub->compatible && (sub->cvar == 1.0);
  meth->domethod = doS;
  meth->S = NULL;
  meth->nr = DOLLAR;
  return err;
}

void doS(method_type *meth, res_type *res){
  int i,
    m=0,
    totalpoints = meth->loc->totalpoints;
  double sd = sqrt(meth->sub[m]->cvar);
  do_meth dometh = meth->sub[m]->domethod;

//   printf("doS %f  %f %s\n", res[0], sd, CovList[meth->cov->nr].name);
  dometh(meth->sub[m], res); 
  // printf("doS  emitte %f\n", res[0]);

  if (sd != 1.0) {
    for (i=0; i<totalpoints; i++) {
	res[i] *= (res_type) sd;
    }
  }
  //  printf("doS  end%f\n", res[0]);
 
}


/* simplifying functions
 turn vector of x into length ||x|| and similar
 reduces C(x,y) to C(x - y)

 needs preceeding analysis of submodels!
 (same for preferred methods; bottom-up-analysis needed)
*/

// #
// keep always the orderung
void iso2iso(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double y=fabs(*x);
  CovList[next->nr].cov(&y, next, v);
}
void spiso2spiso(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  CovList[next->nr].cov(y, next, v);
}
void spacetime2iso(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  CovList[next->nr].cov(&y, next, v);
}

void Stat2iso(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double b = 0.0;
  int i,
    dim=cov->xdim;  

  // printf("dim=%d, %d\n", dim, x[0]);

  for (i=0; i<dim; i++) {
    b += x[i] *x[i];
  }
  b = sqrt(b);
  CovList[next->nr].cov(&b, next, v);
}
void Nonstat2iso(double *x, double *y, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double a, b;
  int dim=cov->xdim, i;  
  for (b=0.0, i=0; i<dim; i++) {
    a = y[i] - x[i];
    b += a * a;
  }
  b = sqrt(b);
  CovList[next->nr].cov(&b, next, v);
}
void Stat2spacetime(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double b, z[2];
  int i,
    dim=cov->xdim - 1;  
  for (b=0.0, i=0; i<dim; i++)  b += x[i] * x[i];
  z[0] = sqrt(b);
  z[1] = fabs(x[dim]);
  CovList[next->nr].cov(z, next, v);
}
void Nonstat2spacetime(double *x, double *y, cov_model *cov, double *v) {
  //assert(false);
  cov_model *next = cov->sub[0];
  double a, b, z[2];
  int dim=cov->xdim - 1, i;  
  for (b=0.0, i=0; i<dim; i++) {
    a = y[i] - x[i];
    b += a * a;
  }
  z[0] = sqrt(b);
  z[1] = fabs(y[dim] - x[dim]);
  CovList[next->nr].cov(z, next, v);
}
void Stat2Stat(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  CovList[next->nr].cov(x, next, v);
}
void Nonstat2Stat(double *x, double *y, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int i,
      dim=cov->xdim;  
  double *z = (double*) malloc(dim * sizeof(double));
  for (i=0; i<dim; i++) z[i] = y[i] - x[i];
  CovList[next->nr].cov(z, next, v);
  free(z);
}
void D_2(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  if (cov->isoIn == ISOTROPIC) {
    double y = fabs(*x);
    PrintModelInfo(cov->calling->calling->calling);  // ok
    PrintModelInfo(cov);  // ok
    
// iso2iso ueberpruefen !!!
// dollar + scale // aniso > 0 und dim = 1
// nach tbm eigentlich auch nicht
//    cov->sub[1000] = cov->sub[1];
    assert(false); // sollte nicht auftreten !!

    CovList[next->nr].D(&y, next, v);
  } else {
    assert(cov->isoIn == SPACEISOTROPIC);
    cov_model *next = cov->sub[0];
    cov_fct *N = CovList + next->nr;
    
    if (N->isotropy==ISOTROPIC) {
      assert(next->xdim == 1);
      double y=sqrt(x[0] * x[0] + x[1] * x[1]);    
      N->D(&y, next, v);
      if (y!=0.0) *v = x[0] / y;
    } else {
      assert(N->isotropy == SPACEISOTROPIC);
      double y[2];
      y[0] = fabs(x[0]);
      y[1] = fabs(x[1]);
      N->D(y, next, v); 
    }
  }
}
void DD_2(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->isoIn == ISOTROPIC) {
    double y = fabs(*x);
    PrintModelInfo(cov->calling->calling->calling);  // ok
    
// iso2iso ueberpruefen !!!
// dollar + scale // aniso > 0 und dim = 1
// nach tbm eigentlich auch nicht
    assert(false); // sollte nicht auftreten !!
    CovList[next->nr].D2(&y, next, v);
  } else {
    assert(cov->isoIn == SPACEISOTROPIC);
    cov_model *next = cov->sub[0];
    cov_fct *N = CovList + next->nr;
    
    if (N->isotropy==ISOTROPIC) {
      assert(next->xdim == 1);
      double w,
	xSq = x[0] * x[0],
	tSq = x[1] * x[1],
	ySq = xSq + tSq,
	y   = sqrt(ySq); 
      
      // (c'(r) * x/r)' = c''(r) * x^2/r^2 + c'(r) [ 1/r - x^2 / r^3]
      N->D2(&y, next, v);
      if (y != 0.0) {
	 N->D(&y, next, &w);
	 w /= y;
	 *v = (*v - w) * xSq / ySq + w;
      }
      else {      
     *v = x[0] / y;
      }
    } else {
      assert(N->isotropy == SPACEISOTROPIC);
      double y[2];
      y[0] = fabs(x[0]);
      y[1] = fabs(x[1]);
      N->D2(y, next, v); 
    }
  }
}


void inv2(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double y=fabs(*x);
  CovList[next->nr].invSq(&y, next, v);
}

// void hyper2


void setdefault(cov_model *cov) {
  // Achilles-Ferse: setdefault wird aufgerufen bevor die Varianten
  // der Modelle und Operatoren feststehen -- die Parameter unten muessen, falls
  // notwendig von Hand in den checks gesetzt werden.

//  cov_model *prev=cov->calling;
  int i;
  cov_fct *C = CovList + cov->nr;
  cov->maxdim = INFDIM;
  cov->normalmix = 
    cov->finiterange = 
    cov->diag = 
//    cov->quasidiag = 
    cov->semiseparatelast = 
    cov->separatelast = true;
//  cov->variogram = false;
//  for(i=0; i<cov->xdim; i++) cov->idx[i] = i;


//  int in=1, op=1;
//  PrintModelList(&in, &op);

  memcpy(cov->pref, C->pref, sizeof(pref_shorttype));  
  for (i=Nothing+1; i<Forbidden; i++) cov->pref[i] = PREF_NONE;
  for (i=0; i<Forbidden; i++) {
      //     printf("%s %s %d\n", C->name, METHODNAMES[i], C->pref[i]);
    if (cov->user[i] < cov->pref[i])
      cov->pref[i] = cov->user[i];
  }


  
//  if (prev != NULL) {
//    cov->tsdim = prev->tsdim;
//    cov->xdim = prev->xdim;
//  }
  cov->tbm2num = C->implemented[TBM2] == NUM_APPROX;

  if (C->isotropy == ISOTROPIC) {
    cov->xdim = 1;
  } else if (C->isotropy == SPACEISOTROPIC) {
    cov->xdim = 2;
  }
}



int check2intern(cov_model *cov){
  bool ANI = false; // GLOBAL.general.aniso;
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + cov->nr,
    *N = CovList + next->nr;
  isotropy_type iso, first_iso, last_iso, iso0;
  stationary_type stat, first_stat, last_stat; 
  int i, err;
  bool skipchecks = GLOBAL.general.skipchecks;

//  printf(">>> %s \n", N->name);

  if (cov->xdim < 1) ERR("dimension less than 1");
  if (PL > 7) PRINTF("checking %s ~ %s\n", C->name, CovList[next->nr].name);

  if (next->nr == GATTER) {
    PRINTF("------------\n");
    PrintModelInfo(cov); // ok
    PRINTF("------------\n");
    assert(false);
  }
  setdefault(cov);
  next->tsdim = cov->tsdim;
  first_iso = last_iso = next->isoIn = N->isotropy;
  first_stat = last_stat = next->statIn = N->stationary;
  
  // erst by check unten
  strcpy(ERROR_LOC, CovList[next->nr].name);
  if (PL > 8) PRINTF("%s\n", ERROR_LOC);
  next->xdim = cov->xdim; // if next is isotropy or spaceisotropic it is 
  //          set to 1 or 2
  setdefault(next);

//  if (CovList[next->nr].primitive) {//29.12.09 geloescht; jetzt in simu.cc l416
//  RUECKGEANDERT 15.1.11 -- darf auf keinen Fall in simu.cc, da sonst der 
// Aufbau gestoert wird!!
  if (CovList[next->nr].primitive &&
      (err = checkkappas(next)) != NOERROR) return err;

  sprintf(ERROR_LOC, "the starting point to `%s'", CovList[next->nr].name);
  if (PL > 8) PRINTF("%s\n", ERROR_LOC);

//  printf("check2intern %s %d %d; %s %d %d\n", 
//	 CovList[cov->nr].name, cov->statIn, cov->isoIn,
//	 N->name, N->stationary, N->isotropy);
	 

  if (first_iso == PREVMODELI) {
    first_iso = ISOTROPIC;
    last_iso = ANISOTROPIC;
  }
  if (last_iso > cov->isoIn) last_iso = cov->isoIn;
  if (last_iso < first_iso) {
    sprintf(ERRORSTRING, 
	    "cannot call non-isotropic model from isotropic one (%s -> %s)", 
	    ISONAMES[(int) cov->isoIn], ISONAMES[(int) next->isoIn]);
    return ERRORM;
  }   
 
  if (first_stat == PREVMODELS) {
    first_stat = STATIONARY;
    last_stat = GENERALISEDCOVARIANCE;
  }
  if (last_stat > cov->statIn) last_stat = cov->statIn;
  if (last_stat < first_stat) {
    // PRINTF("code %d %d\n", prev->stationary, cov->stationary);     
    sprintf(ERRORSTRING, 
	    "cannot call more complex model from less complex one (%s -> %s)",
	    STATNAMES[(int) cov->statIn], STATNAMES[(int) next->statIn]);
    return ERRORM + 1;
 }
  
  err = ERRORNOSTATMATCH;
  next->tsdim = cov->tsdim;
  strcpy(ERROR_LOC, CovList[next->nr].name);
  if (PL > 4)
    PRINTF("%s (stat.start=%d, end=%d)\n", ERROR_LOC, first_stat, last_stat);
    
  for (stat = first_stat; stat <= last_stat; stat++) {
    next->statIn = stat;
    for (iso0=first_iso; iso0 <= last_iso; iso0++) {


//	printf("GEEE %d %d %d iso %d %d %d\n",
//	       stat, first_stat, last_stat, iso0, first_iso, last_iso);

      iso = iso0;
      if (ANI && last_iso >= ANISOTROPIC && first_iso <= SPACEISOTROPIC) {
	// bei mpp darf keine Reduktion stattfinden !
	if (iso == ANISOTROPIC) iso = SPACEISOTROPIC;
	else if (iso == SPACEISOTROPIC) iso = ANISOTROPIC;
      }


      next->derivatives = N->derivs;
 
//      printf("xxx GEEE dsaf %d %d %s %s\n", cov->xdim, N->derivs, N->name,
//	  CovList[next->nr].name);      

      if (cov->xdim == 1) {
	  if ((iso == SPACEISOTROPIC) || (iso == ZEROSPACEISO) || 
	      (iso != ISOTROPIC && (stat == STATIONARY || stat == VARIOGRAM ||
				 stat == IRF) && 
	       cov->vdim<2 && cov->vdim != SUBMODELM)
	      ) continue; 
      }
      
     next->xdim =
	iso == ISOTROPIC ? 1 : iso == SPACEISOTROPIC ? 2 : cov->xdim;
      next->isoIn = iso;
      next->vdim = N->vdim;
//   printf("zz GEEE dsaf %d %d %d %s %s\n",
//	   cov->xdim, N->derivs, next->derivatives, N->name,
//	   CovList[next->nr].name);      
 //      printf("%s GEEE dsaf %d\n", CovList[next->nr].name, err);

      if ((err = CovList[next->nr].check(next)) == NOERROR) break;
//      printf("%s GEEd sd E\n", CovList[next->nr].name);

     }
    if (err == NOERROR) break;
  }
  if (err != NOERROR) return err;
 
  if (next->vdim <= 0) {
//    PrintModelInfo(cov); // Ok
    ERR("something is wrong with the m-dimensionality");
  }
  cov->vdim = next->vdim;

  // printf("bback\n");
  // PrintModelInfo(cov);
    
  sprintf(ERROR_LOC, "Back from `%s':", CovList[next->nr].name);
  if (PL > 8) PRINTF("%s\n", ERROR_LOC);


  if (!skipchecks && (err = check_within_range(next, NAOK_RANGE)) != NOERROR) { 
      return err;
  }

  if (CovList[next->nr].primitive) {
    /* range check of submodel if Primitive !! */

    for (i=0; i<Forbidden; i++) {
      if (next->user[i] < next->pref[i])
	next->pref[i] = next->user[i];
    }
  }
  cov->derivatives = (cov->statIn == STATIONARY || cov->statIn == VARIOGRAM) &&
      ((cov->isoIn == ISOTROPIC) || 
       (cov->isoIn == SPACEISOTROPIC && next->isoIn ==SPACEISOTROPIC)) ? 2 : 0;

   
//  printf("yy GEEE dsaf xdim=%d der=%d n=%d %s; %s %d\n",
//	 cov->xdim, CovList[next->nr].derivs, next->derivatives, N->name,
//	 CovList[cov->nr].name, cov->derivatives);      

  setbackward(cov, next);  
  return NOERROR;
}


int check2(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int  err = ERRORNOSTATMATCH;
  stationary_type stat, first, last;
  
  cov->manipulating_x = !true;  // exception as treated directly
  if (cov->statIn == STATIONARY) {
    first = STATIONARY;
    last = IRF;
/////////////
//    last = STATIONARY;
  } else {
    first = COVARIANCE;
    last = GEN_VARIOGRAM;
  }

  cov->vdim = SUBMODELM;
  cov->isoIn = ANISOTROPIC;
  for (stat=first; stat <= last; stat++) {    
      // printf("%s%s %s [%s %s]\n", CovList[cov->nr].name, 
//	     CovList[cov->sub[0]->nr].name, 
//	     STATNAMES[(int) stat], STATNAMES[(int) first],
//	     STATNAMES[(int) last]);
    cov->statIn = stat;
    if ((err = check2intern(cov)) == NOERROR) break;
//    printf("intern   GEEE dsaf %d \n", err);
    if (PL > 1) {
      PRINTF("error when assuming `%s':\n", STATNAMES[(int) stat]);
      ErrorMessage(Nothing, err);
    }
  }


  
  if (err != NOERROR) return err;

  cov->statIn = next->statIn;
  if (next->statIn == STATIONARY || next->statIn == GENERALISEDCOVARIANCE || 
      next->statIn == VARIOGRAM ||
      next->statIn == IRF) {    
    switch(next->isoIn) {
	case ISOTROPIC :
	  cov->nr = S2ISO; // auch bei xdim == 1, da hier 
	  // nicht S2ISO sondern NonStat2Iso gemeint ist!
	  break;
	case SPACEISOTROPIC : 
	  cov->nr = S2SP;
	  break;
  	case ZEROSPACEISO: case ANISOTROPIC :
	  cov->nr = S2S;
	  break;
	default : assert(false);
    }
  } else {
    if (next->statIn == AUXMATRIX || next->statIn == AUXVECTOR)
      ERR("auxiliary functions are not allowed on top level");
    if (next->statIn != COVARIANCE && next->statIn != GEN_VARIOGRAM)
      ERR("the translation invariance property of the model is unknown");
    cov->tsdim = DEL_COV - 1;
  }
  return NOERROR;
}


int check2X(cov_model *cov, int tsdim, int xdim,
	    stationary_type statprev, isotropy_type isoprev,
	    int vdim
	    ) {
  cov_model *next = cov->sub[0],
    *prev = cov->calling;
  int err;
  
  // long double  sss;

//  PrintModelInfo(cov);

//  printf("%d %d\n", sizeof(long double), sizeof(double)); assert(false);
  
  assert(prev != NULL);
  assert(cov->nr >= GATTER && cov->nr <= LASTGATTER);
  cov->isoIn = isoprev;
  cov->statIn = statprev;
  cov->xdim = xdim;
  cov->tsdim = tsdim;
  cov->vdim = vdim; // e.g. nugget is adapted to required vdim

//  assert(xdim > 1 ||  // deleted 28.12.08 da Aufrufe von ANISO models mit xdim=1fehlerhaft
//	 isoprev == ISOTROPIC || vdim > 1 ||
//	 statprev != STATIONARY && statprev != VARIOGRAM && statprev != IRF);
  sprintf(ERROR_LOC, "#[%s -> %s] (%d %d)",
	  CovList[prev->nr].name,  CovList[next->nr].name, statprev, isoprev);
  if (PL > 4) PRINTF("%s\n", ERROR_LOC);

  if ((err = check2intern(cov)) != NOERROR) {   
      //     printf("2X err=%d\n", err);
      return err;
  }
  
  if (vdim > 0) {
    if (vdim != cov->vdim) ERR("the m-dimensionality does not match");
  } else if (vdim != SUBMODELM)
    ERR("preceding model has unexpected m-dimensionality")
//  printf("Xhere\n");

  if (next->statIn == AUXMATRIX || next->statIn == AUXVECTOR) {
    cov->tsdim = DEL_COV - 2;
    return NOERROR;
  }

  if (isoprev == SPACEISOTROPIC) {
    assert(cov->xdim == 2);
    cov->pref[TBM2] = 0;
  }

  if (statprev < next->statIn) 
    ERR("cannot call more complex models from simpler ones");

  if (statprev == GENERALISEDCOVARIANCE && 
      (next->statIn != GENERALISEDCOVARIANCE && next->statIn != STATIONARY))
    ERR("translation invariance does not match with generalised covariance");

  
  switch(statprev) {
      case STATIONARY : case VARIOGRAM : case IRF : case AUXMATRIX :
      case GENERALISEDCOVARIANCE :	
	switch(isoprev) {
	    case ISOTROPIC :
		cov->nr = ISO2ISO; // iso2iso
	      break;
	    case SPACEISOTROPIC :
	      cov->nr = (next->isoIn == ISOTROPIC) ? SP2ISO : SP2SP;
	      break;	      
	    case ZEROSPACEISO: case ANISOTROPIC :
	      switch (next->isoIn) {
		  case ISOTROPIC :
		    cov->nr = S2ISO;
		    // if (prev->quasidiag) 
		    //    for (i=0; i<xdim; i++) prev->idx[i] = i;
		    
		    break;
		  case SPACEISOTROPIC :
		    cov->nr = S2SP;
		    // if (cov->quasidiag) {
		    // if ((cov->diag = cov->quasidiag = next->idx[1] ==1))
		    //  { for (i=0; i<xdim; i++) prev->idx[i] = i; }
		    break;
		  case ZEROSPACEISO: case ANISOTROPIC :
		    cov->tsdim = DEL_COV - 5;
		    break;
		  default : assert(false);
	      }
	      break;
	    default: 
	      PRINTF("%s %d %d\n", CovList[next->nr].name, statprev, isoprev);
	      assert(false);
	}
	break;
      case AUXVECTOR :
	if (next->statIn != AUXVECTOR) return ERRORAUXVECTOR;
	break;
      case COVARIANCE : case GEN_VARIOGRAM :
	  if (next->statIn == STATIONARY || next->statIn == VARIOGRAM || 
	    next->statIn == IRF || next->statIn==AUXMATRIX){
	    switch(next->isoIn) {
		case ISOTROPIC :
		  cov->nr = S2ISO;
		  // if (prev->quasidiag) 
		  // for (i=0; i<xdim; i++) prev->idx[i] = i;
		  break;
		case SPACEISOTROPIC : 
		  cov->nr = S2SP;
		  // if (prev->quasidiag) {
		  // if ((prev->diag = prev->quasidiag = cov->idx[1] == 1)) {
		  // for (i=0; i<xdim; i++) prev->idx[i] = i;
		  // }
		  //}
		  break;
		case ZEROSPACEISO: case ANISOTROPIC :
		  cov->nr = S2S; 
		  //if (prev->quasidiag) 
		  //  for (i=0; i<xdim; i++) prev->idx[i] = cov->idx[i];
		  
		  break;
		default: assert(false);
	    }
	  } else { // still non-stationary: cannot be simplified
	    cov->tsdim = DEL_COV - 4;
	  }
	  break;

	default : 
	  PrintModelInfo(cov); // ok
	  PRINTF("%d %s\n", prev->statIn, CovList[prev->nr].name);
	  assert(false);
  }
  return NOERROR;
}

void range2(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->finiterange = true;
  range->maxdim = INFDIM;
}
int init2(method_type *meth){
  int err;
  cov_model *newcov=NULL, 
      *cov = meth->cov,
      *next = cov->sub[0];
  meth->hanging = meth->cov;
  
  if (CovList[next->nr].initmethod == NULL)
      return ERRORSUBMETHODFAILED;

  covcpy(&newcov, next, true, true);
  assert(newcov->nr >= GATTER && newcov->nr<=LASTGATTER);
  removeOnly(&newcov);
  setdefault(newcov);
  newcov->statIn = cov->statIn;
  newcov->isoIn = cov->isoIn;
  newcov->xdim = cov->xdim;
//  PrintModelInfo(newcov); assert(false);
  if ((err = CovList[newcov->nr].check(newcov)) != NOERROR) {
      COV_DELETE(&newcov);
      XERR(err);
  }
  DeleteGatter(&newcov);
  // 

// PrintModelInfo(newcov);assert(false);

  meth->sub[0] = (method_type*) malloc(sizeof(method_type));
  METHOD_NULL(meth->sub[0]);
  meth->nsub = 1;
  method_type *sub = meth->sub[0];
  cpyMethod(meth, sub, false);
  sub->cov = newcov;

//  PrintMethodInfo(meth);
  if (meth->gp->general.printlevel >= 5) {
//      PRINTF(" trying initialisation of model `%s' [%d; %ld].\n",
//	     CovList[next->nr].name, next->nr, 
//	     (POINTER) CovList[next->nr].initmethod);
   }

  if ((err = CovList[next->nr].initmethod(sub)) != NOERROR) {
    return err;
  }

  meth->compatible = sub->compatible;
  meth->domethod = do2;
  meth->S = NULL;
  meth->nr = GATTER;

  return NOERROR;
}
void do2(method_type *meth, res_type *res){
    meth->sub[0]->domethod(meth->sub[0], res);   
}
int initspectral2(cov_model *cov) {
  cov_model *next = cov->sub[0];

//  printf("initsp2 %ld %s\n",
//	 CovList[next->nr].initspectral, CovList[next->nr].name);
//  assert(CovList[next->nr].initspectral != NULL);
  if (CovList[next->nr].initspectral == NULL) {
      return ERRORFAILED; 
  }
//  printf("sp err gatter null %s %d\n", CovList[next->nr].name, CovList[next->nr].initspectral(next));
  memcpy(&(next->spec), &(cov->spec), sizeof(spec_covstorage));
  int err=CovList[next->nr].initspectral(next);
  memcpy(&(cov->spec), &(next->spec), sizeof(spec_covstorage));

  return err;
}

void spectral2(cov_model *cov, spectral_storage *s, double *e) {
  cov_model *next = cov->sub[0];
  return CovList[next->nr].spectral(next, s, e);
}

void mppinit2(mpp_storage *s, cov_model *cov) {
    cov_model *next = cov->sub[0];
    CovList[next->nr].mppinit(s, next);
}
void coin2(mpp_storage *s, cov_model *cov){
    cov_model *next = cov->sub[0];
    CovList[next->nr].randomcoin(s, next);  
}
void sd2(mpp_storage *s, cov_model *cov) {
   cov_model *next = cov->sub[0];
   CovList[next->nr].sd(s, next);  
}

res_type mppgetiso2iso(double *x, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  double y=fabs(*x);
  return CovList[next->nr].mppgetstat(&y, next, s);
}
res_type mppgetspiso2spiso(double *x, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  return CovList[next->nr].mppgetstat(y, next, s);
}
res_type mppgetspacetime2iso(double *x, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  return CovList[next->nr].mppgetstat(&y, next, s);
}

res_type mppgetStat2iso(double *x, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  double b = 0.0;
  int i,
    dim=cov->xdim;  
   for (i=0; i<dim; i++) {
    b += x[i] *x[i];
  }
  b = sqrt(b);
  return CovList[next->nr].mppgetstat(&b, next, s);
}
res_type mppgetNonstat2iso(double *x, double *y, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  double a, b;
  int dim=cov->xdim, i;  

//  PrintModelInfo(cov); assert(false);

  for (b=0.0, i=0; i<dim; i++) {
    a = y[i] - x[i];
    b += a * a;
  }
  b = sqrt(b);
  return CovList[next->nr].mppgetstat(&b, next, s);
}
res_type mppgetStat2spacetime(double *x, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  double b, z[2];
  int i,
    dim=cov->xdim - 1;  
  for (b=0.0, i=0; i<dim; i++)  b += x[i] * x[i];
  z[0] = sqrt(b);
  z[1] = fabs(x[dim]);
  return CovList[next->nr].mppgetstat(z, next, s);
}
res_type mppgetNonstat2spacetime(double *x, double *y, cov_model *cov, 
			     mpp_storage *s) {
  //assert(false);
  cov_model *next = cov->sub[0];
  double a, b, z[2];
  int dim=cov->xdim - 1, i;  
  for (b=0.0, i=0; i<dim; i++) {
    a = y[i] - x[i];
    b += a * a;
  }
  z[0] = sqrt(b);
  z[1] = fabs(y[dim] - x[dim]);
  return CovList[next->nr].mppgetstat(z, next, s);
}
res_type mppgetStat2Stat(double *x, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  return CovList[next->nr].mppgetstat(x, next, s);
}
res_type mppgetNonstat2Stat(double *x, double *y, cov_model *cov, mpp_storage *s) {
  cov_model *next = cov->sub[0];
  int i,
      dim=cov->xdim;  
  double *z = (double*) malloc(dim * sizeof(double));
  for (i=0; i<dim; i++) z[i] = y[i] - x[i];
  return CovList[next->nr].mppgetstat(z, next, s);
  free(z);
}





cov_model *covOhneGatter(cov_model *cov) {
  return (cov->nr >= GATTER && cov->nr <= LASTGATTER) ? cov->sub[0] : cov;
}

int checkplusmal(cov_model *cov) {
  cov_model *sub;
  int i, err; // taken[MAX DIM],

  
  cov->manipulating_x = X_PASSED;
  for (i=0; i<cov->nsub; i++) {
    sub=cov->sub[i];
    if (sub == NULL)
      ERR("+ or *: named submodels are not given in a sequence!");

    if ((err = check2X(sub, cov->tsdim, cov->xdim, cov->statIn, cov->isoIn,
		       SUBMODELM)) != NOERROR) {
	// printf("error\n");
	return err;
    }
    if (i==0) cov->vdim=sub->vdim; 
    else if (cov->vdim != sub->vdim) 
      ERR("m-dimensionality must be equal in the submodels ");
    setbackward(cov, sub);//cov->sub[i] may now point to somewhere else!
    if (sub->sub[0]->statIn == cov->statIn &&
	sub->sub[0]->isoIn == cov->isoIn) {
      sub->tsdim = DEL_COV - 6;
    } else {
      if (sub->statIn != cov->statIn || sub->isoIn != cov->isoIn)
	return ERRORFAILED; 
    }
  }

  // !! incorrect  !!
  cov->semiseparatelast = false; // taken[xdim - 1] <= 1;
  cov->separatelast = false; // taken[xdim - 1] <= 1; ?? 
  return NOERROR;
}

void plusStat(double *x, cov_model *cov, double *v){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  double *z = (double*) malloc(sizeof(double) * vsq);
  
  assert(x[0] >= 0.0 || cov->xdim > 1);
  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    CovList[sub->nr].cov(x, sub, z);
//    printf("(%s %f %f) ", CovList[sub->sub[0]->nr].name, *x, *z);
    for (m=0; m<vsq; m++) v[m] += z[m]; 
  }
  free(z);
}

void plusNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  double *z = (double*) malloc(sizeof(double) * vsq);
  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    CovList[sub->nr].nonstat_cov(x, y, sub, z);
    for (m=0; m<vsq; m++) v[m] += z[m]; 
  }
  free(z);
}
void MLEplusStat(double *x, cov_model *cov, double *v){
  cov_model *sub;
  if (ELEMENTNR_PLUS >= 0) {
      sub = cov->sub[ELEMENTNR_PLUS];
      CovList[sub->nr].cov(x, sub, v);
  } else {
      plusStat(x, cov, v); // !!! cov not next 
  }

//  printf("plus %s %d  %f, v=%f\n", CovList[sub->nr].name, ELEMENTNR_PLUS, *x, *v);

}

void MLEplusNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *sub;
  if (ELEMENTNR_PLUS >= 0) {
     sub = cov->sub[ELEMENTNR_PLUS];
     CovList[sub->nr].nonstat_cov(x, y, sub, v);
  } else {
      plusNonStat(x, y, cov, v); // !!! cov not next 
  }
}

void Dplus(double *x, cov_model *cov, double *v){
  cov_model *sub;
  double v1;
  int n = cov->nsub, i;
  *v = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
    CovList[sub->nr].D(x, sub, &v1);
    (*v) += v1;
  }
}

void DDplus(double *x, cov_model *cov, double *v){
  cov_model *sub;
  double v1;
  int n = cov->nsub, i;
  *v = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
    CovList[sub->nr].D2(x, sub, &v1);
    (*v) += v1;
  }
}

int checkplus(cov_model *cov) {
  int err;
  if ((err = checkplusmal(cov)) != NOERROR) {
    return err;
  }
  if (cov->statIn == STAT_MISMATCH) XERR(ERRORNOVARIOGRAM);
  return NOERROR;

  // spectral mit "+" funktioniert, falls alle varianzen gleich sind,
  // d.h. nachfolgend DOLLAR die Varianzen gleich sind oder DOLLAR nicht
  // auftritt; dann zufaellige Auswahl unter "+"
}


int initspectralplus(cov_model *cov) {
// #define spec_ergodic !true
  int i;
  cov_model *sub;
  spec_covstorage *cs = &(cov->spec);
  double *sd_cum = cs->sub_sd_cum;
//  if (spec_ergodic) return ERRORFAILED;
  for (i=0; i<cov->nsub; i++) {
    sub = cov->sub[i];
    CovList[sub->nr].initspectral(sub);
    CovList[sub->nr].cov(ZERO, sub, sd_cum + i);
    if (i>0) sd_cum[i] += sd_cum[i-1];
  }
  return NOERROR;
}
void spectralplus(cov_model *cov, spectral_storage *s, double *e){
  int nr;
  double dummy;
  spec_covstorage *cs = &(cov->spec);
  double *sd_cum = cs->sub_sd_cum;

  nr = cov->nsub - 1;
  dummy = UNIFORM_RANDOM * sd_cum[nr];
  while (nr > 0 && dummy <= sd_cum[nr - 1]) nr--;
  cov_model *next = cov->sub[nr];
  CovList[next->nr].spectral(next, s, e);
}

void rangeplus(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->finiterange = true;
  range->maxdim = INFDIM;
}


void plus_destruct(void ** S) {
  if (*S!=NULL) {
    free(*S);
    *S = NULL;
  }
}

int initplus(method_type *meth){
  cov_model *cov=meth->cov;
  int m, incomp, pos_incomp=-1, err;


  /*
  pref_type pref;

  for (m=0; m<Forbidden; m++) pref[m] = PREF_BEST;
  if (meth->cov->user[0] == 0 || meth->cov->user[1] == 0) {

 
     //PrintModelInfo(meth->cov); assert(false);

      // i.e. user defined
      pref_type pref;
      pref[CircEmbed] = pref[TBM2] = pref[TBM3] = pref[Direct] = pref[Sequential] 
	  = PREF_NONE;
      for (m=0; m<Nothing; m++) { // korrekt auch fuer MaxStable?
//	  printf("%d %d\n", m, pref[m]);
//	  printf("%d %d\n", meth->cov->user[m]);
	  if (pref[m] > 0 &&  meth->cov->user[m] > 0) {
	      break;
	  }
      }
      if (m == Nothing) return ERRORSUBMETHODFAILED;
  }
  */


  SET_DESTRUCT(plus_destruct);
  incomp = 0; 
  meth->nsub = cov->nsub; // could be improved somewhen in future

  for (m=0; m<cov->nsub; m++) {
    // hole m-tes element von "plus"; ersetze "plus" durch "identity" --
    // kann man besser machen
    meth->sub[m] = (method_type*) malloc(sizeof(method_type));
    method_type *sub = meth->sub[m];
    METHOD_NULL(sub);
    cpyMethod(meth, sub, true);
    
    // hier sollte nicht gelinkt, sondern kopiert werden,
    // dann koente obiges pref verwendet werden, um 
    // gewisse Verfahren durch minimumsbildung der prefs
    // auszublenden. Erst dann interessant, wenn User
    // mehrer Moeglichkeiten angeben kann.
    sub->cov = cov->sub[m];
    

    //   printf("initplus:");
    //   PrintModelInfo(sub->cov);
    if (meth->gp->general.printlevel >= 5) {
	PRINTF(" trying initialisation of submodel #%d (initstandard).\n", m);
    }  

    err = initstandard(sub); 

    if (err != NOERROR) return err;
    if (!sub->compatible) {
      pos_incomp = m;
      incomp ++;
    }
  }

  meth->compatible = incomp != 1;
  if (!meth->compatible && pos_incomp > 0) {
    method_type *dummy;
    dummy = meth->sub[pos_incomp];
    meth->sub[pos_incomp] = meth->sub[0];
    meth->sub[0] = dummy;
  } else {
    meth->S =  (incomp > 1) 
      ? (double*) malloc(meth->loc->totalpoints * sizeof(double))
      : NULL;
  }
  char plussign[] = "+";
  meth->nr = getmodelnr(plussign);
  meth->domethod = doplus;
  return NOERROR;
}


void doplus(method_type *meth,  res_type *res){
  int m;
  method_type **sub=meth->sub;
   
  if (meth->S == NULL) {    
      // assert(false);
    for (m=0; m<meth->nsub; m++) {	
	sub[m]->domethod(sub[m], res); 
      }
  } else {
    int i, t=meth->loc->totalpoints;
    m = 0;
    sub[m]->domethod(sub[m], res);   
//    printf("doplus %f\n", (double) res[0]);
    for (m=1; m<meth->nsub; m++) {
//
//	printf("%d doplus comp %d %f\n", m, sub[m]->compatible, (double) res[0]);
      if (sub[m]->compatible) {
        sub[m]->domethod(sub[m], res);   
      } else { 
	for(i=0; i<t; i++)((double *)meth->S)[i] = 0.0;
        sub[m]->domethod(sub[m], (res_type *) meth->S); 
        for(i=0; i<t; i++) res[i] += ((res_type *)meth->S)[i];
//	printf("%d doplus non comp %f\n", m, res[0],  ((double *)meth->S)[0]);
     }
    }
  }
}



void malStat(double *x, cov_model *cov, double *v){
  cov_model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  double *z = (double*) malloc(sizeof(double) * vsq);
  
  assert(x[0] >= 0.0 || cov->xdim > 1);
  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    CovList[sub->nr].cov(x, sub, z);
    for (m=0; m<vsq; m++) v[m] *= z[m]; 
  }
  free(z);
}

void malNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *sub;
  int i, m, nsub=cov->nsub,
    vsq = cov->vdim * cov->vdim;
  double *z = (double*) malloc(sizeof(double) * vsq);
  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    CovList[sub->nr].nonstat_cov(x, y, sub, z);
    for (m=0; m<vsq; m++) v[m] *= z[m]; 
  }
  free(z);
}

int checkmal(cov_model *cov) {
  cov_model *next1 = cov->sub[0];
  cov_model *next2 = cov->sub[1];
  double *aniso1 = next1->p[2], *aniso2 = next2->p[2];
  int err;
  
  if ((err = checkplusmal(cov)) != NOERROR) return err;

  if (cov->statIn == STAT_MISMATCH ||
      (cov->statIn != STATIONARY &&
       cov->statIn != COVARIANCE)
    ) {
    return ERRORNOVARIOGRAM;
  }
  
  cov->pref[TBM2] = 0;
  if (cov->xdim >= 2) cov->pref[TBM3] = 0;
  if (cov->xdim==2 && cov->nsub == 2 && 
      next1->nr >= DOLLAR && next1->nr <= LASTDOLLAR && 
      next2->nr >= DOLLAR && next2->nr <= LASTDOLLAR) {
    if (aniso1[0] == 0.0 && next1->ncol[2] == 1) {
      cov->pref[TBM2] = next2->pref[TBM2];
      cov->pref[TBM3] = next2->pref[TBM3];
    } else if (aniso2[0] == 0.0 && next2->ncol[2] == 1) {
      cov->pref[TBM2] = next1->pref[TBM2];
      cov->pref[TBM3] = next1->pref[TBM3];
    }
  }
  return NOERROR;
}

void Dmal(double *x, cov_model *cov, double *v){
  cov_model *sub;
  cov_fct *C;
  double c[MAXSUB], d[MAXSUB];
  int n = cov->nsub, i;
  for (i=0; i<n; i++) {
    sub = cov->sub[i];
    C = CovList + sub->nr;
    C->cov(x, sub, c + i);
    C->D(x, sub, d + i);
  }
  *v = 0.0;
  for (i=0; i<n; i++) {
    double zw = d[i];
    int j;
    for (j=0; j<n; j++) if (j!=i) zw *= c[j]; 
    (*v) += zw; 
  }
}

void rangemal(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->finiterange = true;
  range->maxdim = INFDIM;
}
int initmal(method_type *meth){
//  int err;
//  return err;
  assert(false);
  meth->domethod = domal;
  char timessign[] = "*";
  meth->nr = getmodelnr(timessign);
  return ERRORFAILED;
}
void domal(method_type *meth, res_type *res){
  assert(false);
}



void location_rules(method_type *meth, pref_type pref) {

  // 1. rules that depend on the the locations and the user's preferences only,
  // but not on the covariance model

  location_type *loc = meth->loc;
  decision_param *dp = &(meth->gp->decision);
//  direct_param *lp = &(meth->gp->direct);
  unsigned int maxmem=500000000;
  int i;

  SimulationType Standard[Nothing] = {
     CircEmbed, CircEmbedIntrinsic, CircEmbedCutoff, SpectralTBM, TBM2, TBM3,
     Direct, Sequential, Markov, Average, Nugget, RandomCoin, Hyperplane
  };
  for (i=0; i<Nothing; i++) {
    pref[Standard[i]] = Nothing - i; 
  }

//  if (loc->totalpoints > max_dirct) pref[Direct] = LOC_PREF_NONE - 0;
  if (dp->stationary_only == DECISION_TRUE) {
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 1;
  }
  if (dp->exactness==DECISION_TRUE) {
    pref[TBM3] = pref[TBM2] = pref[SpectralTBM] = pref[Average] = 
      pref[Markov] = 
      pref[Sequential] = pref[RandomCoin] = pref[Hyperplane] 
      = LOC_PREF_NONE - 2;
  }
  if (loc->grid) {
    if (dp->exactness==DECISION_FALSE && 
	loc->totalpoints * (1 << meth->cov->tsdim) * sizeof(double) >
	maxmem){
      pref[CircEmbed] -= Nothing;
      pref[CircEmbedIntrinsic] -= Nothing;
      pref[CircEmbedCutoff] -= Nothing;
    }
  } else {
    pref[CircEmbed] = pref[CircEmbedIntrinsic] = pref[CircEmbedCutoff] =
      pref[Markov] = LOC_PREF_NONE - 3;
    if (!loc->Time) pref[Sequential] = LOC_PREF_NONE;
  }
}


void mixed_rules(method_type *meth, pref_type locpref, 
		 pref_type pref, int *order) {

  location_type *loc = meth->loc;
  direct_param *lp = &(meth->gp->direct);
  decision_param *dp = &(meth->gp->decision);
//  ce_param *cp = &(meth->gp->ce);
  cov_model *cov = meth->cov;
  int *covpref = cov->pref;
  cov_fct *cf = CovList + cov->nr;
  int i,
    best_dirct=lp->bestvariables,
    max_dirct=lp->maxvariables;

  for (i=0; i<Nothing; i++) {
    pref[i] =  (covpref[i] > PREF_NONE && locpref[i] > LOC_PREF_NONE)
      ? covpref[i] * Nothing + locpref[i]
      : locpref[i] <= LOC_PREF_NONE ?  locpref[i] :  LOC_PREF_NONE - 4;
  }

  if (loc->totalpoints * cov->vdim > max_dirct) pref[Direct] = LOC_PREF_NONE-0;

  if (loc->totalpoints * cov->vdim <= best_dirct && covpref[Direct] == PREF_BEST)
    pref[Direct] = LOC_PREF_BEST;

  if (dp->exactness != DECISION_FALSE && cov->tbm2num)
      pref[TBM2] = LOC_PREF_NONE - 5;

  if (dp->stationary_only==DECISION_CASESPEC && (cf->stationary==STATIONARY))
    pref[CircEmbedIntrinsic] = LOC_PREF_NONE - 6;

  orderingInt(pref, Nothing, 1, order);
}



char PREF_FAILURE[100 * Nothing];
void standard_destruct(void ** S) {
  assert(false);
}
int initstandard(method_type *meth){
  globalparam *gp = meth->gp;
//  decision_param *lp = &(gp->decision);
  pref_type locpref, pref;
  int order[Nothing], i,
      err = ERROROUTOFMETHODLIST,
      errSub = NOERROR;
  SimulationType unimeth;
  int PL = gp->general.printlevel;
  char dummy[100 * Nothing];
  static char FailureMsg[][80] = 
    {"total number of points > maximum value",
     "non-stationary model not allowed",
     "not an exact method as required by user",
     "locations not on a grid",
     "denied by cov. model or by user's pref.",
     "no analytic solutn of Abelian integral",
     "intrinsic fields denied by user"
    };

  assert(meth->destruct==NULL && meth->S==NULL);

  char undefined[] = "undefined";
  meth->nr = getmodelnr(undefined);
  location_rules(meth, locpref);
  mixed_rules(meth, locpref, pref, order);

//  printf("\n **** Entering InitStandard\n");  PrintModelInfo(meth->cov);
//  PrintMethodInfo(meth); assert(false);
 
  if (PL > 4) {
      PRINTF("\n");
    for (i=0; i < Nothing; i++) {
      PRINTF("%-15s:   covprev=%1d   locpref=%4d   totpref=%6d\n", 
	     METHODNAMES[i], meth->cov->pref[i], 
	     (locpref[i] < -9000) ? -999 : locpref[i], pref[i]);
    }
    PRINTF("initstandard %s (%d) [pref value %d]\n", 
	   CovList[meth->cov->nr].name, meth->cov->nr,
	   pref[order[Nothing -1]]);
  }
  
  //  if (pref[order[Nothing -1]] <= 0)
  //   return ERROROUTOFMETHODLIST;
  i = Nothing - 1;
  // printf("hhhhhhhhhhhhhh %d %d %d\n", i, pref[order[i]], PREF_NONE);

  for (i = Nothing - 1; i>=0 && pref[order[i]] > PREF_NONE; i--) {

      //  printf("ggggggggggggggggggggggggg %d\n", i);

    unimeth = (SimulationType) order[i];
    if (PL > 2) PRINTF("%s : pref=%d\n", METHODNAMES[unimeth], pref[unimeth]);
    meth->nr = -unimeth-1;
    meth->compatible = compatible_primitive[unimeth];
    meth->domethod = do_method[unimeth];
    
    //  print("%d\n", unimeth);
    err = init_method[unimeth](meth);  //!!

//    printf("i= %d %d\n", order[i], err);
    //   PrintMethodInfo(meth); assert(false);
 
    if (PL >= 5)
      PRINTF("%d , err =%d [%d, %s]\n", i, err, unimeth, METHODNAMES[unimeth]);
    if (err == NOERROR) {
	if (PL >= 5)  PRINTF("returning to higher level...\n");
      return NOERROR;
    }
    assert(meth->destruct != NULL);
    meth->destruct(&(meth->S));
    meth->destruct = NULL;
  }

  // all failed;
  int nr = meth->cov->nr;
  if (gp->general.printlevel >= 5) {
    PRINTF(" trying direct initialisation of model `%s' [%d; %ld].\n",
   	   CovList[nr].name, nr, CovList[nr].initmethod);
   }  
  
  if (CovList[nr].initmethod == NULL) {
      return ERRORNOTDEFINED;
  }

//  PrintMethodInfo(meth); assert(false);
  if ((errSub = CovList[nr].initmethod(meth)) != NOERROR) {

     if (errSub == ERRORSUBMETHODFAILED || errSub == ERROROUTOFMETHODLIST) {
	
//	 printf("err = %d %d\n", err, errSub);
	 return err;      // not subErr !
    }
    strcpy(PREF_FAILURE, "");
    int p, lp;
#define NMAX 14
    char lpd[255], pd[255], names[NMAX];
    names[NMAX-1] = '\0';
    if (gp->general.printlevel >= 5) {
      for (i=0; i<Nothing; i++) {
        lp = locpref[i] - LOC_PREF_NONE;
	p = pref[i] - LOC_PREF_NONE;
	if (lp>0) {
	    strcpy(lpd, "");
	} else {
	    sprintf(lpd, "%s (locp) ", FailureMsg[-lp]);
	}
	if (p>0 || p == lp) {
	    strcpy(pd, ""); // strcpy(pd, "(method specific failure)");
	} else {
	    sprintf(pd, "%s (pref)", FailureMsg[-p]);
	}
	strncpy(names, METHODNAMES[i], NMAX-1);
	sprintf(dummy, "%s %-13s: %s%s\n", PREF_FAILURE, names, lpd, pd);
	strcpy(PREF_FAILURE, dummy);
      }
    }
    return errSub;
  }
  return NOERROR;
}

void dostandard(method_type *meth, res_type *res){
  assert(false);
  meth->sub[0]->domethod(meth, res); 
}


int init_nothing(method_type *meth) {
  assert(false);
  return ERRORFAILED;
}
void do_nothing(method_type *meth, res_type *res){
  assert(false);
}
