/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Copyright (C) 2006 -- 2011 Martin Schlather

 Definition of correlation functions and derivatives (spectral densities, etc)
 of genuinely anisotropic or non-stationary functions 

 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE may not be used here since these parameters are already 
   considered elsewhere


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
 
#include "RF.h"
#include "Covariance.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

void xA(double *x, double*A, int nrow, int ncol, double *y) {
    int i,j,k;
     for (k=i=0; i<ncol; i++) {
	   y[i] =0.0;
	    for (j=0; j<nrow; j++) {
		y[i] += A[k++] * x[j];
	    }
      }	
}


void Ax(double *A, double*x, int nrow, int ncol, double *y) {
    int i,j,k;
    for (i=0; i<nrow; i++) y[i]=0.0;
    for (k=i=0; i<ncol; i++) { 
	for (j=0; j<nrow; j++) {
	    y[j] += A[k++] * x[i];
	}
    }
}


double XkCXtl(double *X, double *C, int nrow, int dim, int k, int l) {
    double *pX, *pY, scalar, result;
    int i,j,ci,
	size = nrow * dim;
   
  pX = X + k;
  pY = X + l;
  result = 0.0;
  for (ci=0, j=0; j<size; j+=nrow) {
      for (i=0, scalar = 0.0; i<size; i+=nrow) {
	 scalar += pX[i] * C[ci++];
     }
     result += scalar * pY[j];
  }  
  return result;
}


void XCXt(double *X, double *C, double *V, int nrow, int dim) {
    double *pX, *endpX, *dummy, *pdummy, scalar;
    int i, cd, rv, ci, cv,
	size = nrow * dim;
  
  dummy = (double*) malloc(sizeof(double) * nrow * dim);
  if (dummy == NULL) error("memory allocation error in XCXt");
  
  for (pX = X, pdummy=dummy, endpX = pX + nrow;
       pX < endpX; pX++, pdummy++) {
    for (ci=0, cd=0; cd<size; cd+=nrow) {
      for (i=0, scalar = 0.0; i<size; i+=nrow) {
        scalar += pX[i] * C[ci++];
      }
      pdummy[cd] = scalar;
    }
  }

  for (rv = 0; rv<nrow; rv++) {
    for (cv=rv, cd=cv*nrow; cv<size; cd+=nrow) {
      for (scalar=0.0, i=0; i<size; i+=nrow) {
	scalar += dummy[rv + i] * X[cd + i];
      }
      V[rv + cv * nrow] = V[cv + rv * nrow] = scalar;
    }
  } 
  
  free(dummy);
}


double xUy(double *x, double *U, double *y, int dim) {
  double dummy,
    xVy = 0.0;
  int j, d, i,
    dimM1 = dim - 1;
  for (j=d=0; d<dim; d++) {
    j = dim * d;
    for (dummy = 0.0, i=0; i<=d; i++) {
      dummy += x[i] * U[j++];
    }
    for (j += dimM1; i<dim; i++, j+=dim) {
      dummy += x[i] * U[j];
    }
    xVy += dummy * y[d];
  }
  return xVy;
}

double xUxz(double *x, double *U, int dim, double *z) {
  double dummy,
    xVx = 0.0;
  int j, d, i,
    dimM1 = dim - 1;
  for (j=d=0; d<dim; d++) {
    j = dim * d;
    for (dummy = 0.0, i=0; i<=d; i++) {
      dummy += x[i] * U[j++];
    }
    for (j += dimM1; i<dim; i++, j+=dim) { 
      dummy += x[i] * U[j];
    }
    if (z != NULL) z[d] = dummy;
    xVx += dummy * x[d];
  }
  return xVx;
}

double xUx(double *x, double *U, int dim) {
    return xUxz(x, U, dim, NULL);
}

double x_UxPz(double *x, double *U, double *z, int dim) {
// x^t (Ux + z); U dreieckmatrix
  double dummy,
    xVx = 0.0;
  int j, d, i,
    dimM1 = dim - 1;
  for (j=d=0; d<dim; d++) {
    j = dim * d;
    for (dummy = z[d], i=0; i<=d; i++) {
      dummy += x[i] * U[j++];
    }
    for (j += dimM1; i<dim; i++, j+=dim) {
      dummy += x[i] * U[j];
    }
    xVx += dummy * x[d];
  }
  return xVx;
}

double detU(double *C, int dim) {
  /* ACHTUNG!! detU zerstoert !!! */
  int i, info, 
//    job = 10,
    dimP1 = dim + 1,
    dimsq = dim * dim;
  double det = 1.0;

  F77_CALL(dpofa)(C, &dim, &dim, &info); // C i s now cholesky
  if (info != 0) {
    ERR("matrix does not seem to be strictly positive definite");
  }
  for (i=0; i<dimsq; i+=dimP1) det *= C[i];
  return det * det;
}

void det_UpperInv(double *C, double *det, int dim) {
  int i, info, 
    job = 01,
    dimP1 = dim + 1,
    dimsq = dim * dim;
  F77_CALL(dpofa)(C, &dim, &dim, &info); // C i s now cholesky
  if (info != 0) ERR("Inv: dpofa failed -- is matrix positive definite?");

  double Det = 1.0;
  for (i=0; i<dimsq; i+=dimP1) Det *= C[i];
  *det = Det * Det;

  F77_CALL(dpodi)(C, &dim, &dim, det, &job); // C is now Cinv
}



void kappa_ave(int i, int dim, int *nr, int *nc){
  *nr = dim;
  *nc = (i==0) ? dim : (i==1) ? 1 : -1;
}
void kappa_ave1(int i, cov_model *cov, int *nr, int *nc){
  kappa_ave(i, cov->tsdim-1, nr, nc);
}
void kappa_ave2(int i, cov_model *cov, int *nr, int *nc){
  kappa_ave(i, cov->tsdim, nr, nc);
}


void rangeave(int dim, range_arraytype* ra) {
  range_type *range = ra->ranges;
  int i;
  
  for (i=0; i<=1; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -10.0;
    range->pmax[i] = 10.0;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}
void rangeave1(cov_model *cov, range_arraytype* ra){
  rangeave(cov->tsdim - 1, ra);
}
 void rangeave2(cov_model *cov, range_arraytype* ra){
  rangeave(cov->tsdim , ra);
}


int checkave(cov_model *cov, int dim, const char *msg) {
  cov_model *next = cov->sub[0];
  int i, j, err;
//    dimP1 = dim + 1;
  double 
    *A = cov->p[0];
//  char Msg[255];
//  cov_fct *CN = CovList + next->nr;
  assert(cov->xdim == cov->tsdim);

  if (cov->xdim > AveMaxDim) {
      sprintf(ERRORSTRING, 
	      "For technical reasons max. dimension for ave is %d. Got %d.", 
	      AveMaxDim, cov->xdim);
    return ERRORM + 2;
  }

  if (CovList[cov->nr].avef != 0) 
    cov->pref[Average] = cov->user[Average];

  cov->manipulating_x=true;
  if (cov->ncol[0] != dim || cov->nrow[0] != dim) {
    sprintf(ERRORSTRING, "A not %sx%s matrix, but %dx%d (dim=%d)", msg, msg,
	    cov->ncol[0], cov->nrow[0], dim);
    return ERRORM + 2;
  }
  if (cov->ncol[1] != 1 || cov->nrow[1] != dim) {
    sprintf(ERRORSTRING, "z not (%s)-dim vector", msg);
    return ERRORM + 3;
  }

  for (i=0; i<dim; i++)
    for (j=i+1; j<dim; j++)
      if (A[i + j * dim] != A[j + i * dim]) {
	A[j + i * dim] = A[i + j * dim];
	warning("A is not symmetric -- lower part used");
      }
  if (!is_positive_definite(A, dim)) ERR("A is not positive definite");

  if ((err = checkkappas(cov)) != NOERROR) return err;
  if ((err = check2X(next, 1, 1, STATIONARY, ISOTROPIC,
		     UNIVARIATE)) != NOERROR) return err;
  next->tsdim = DEL_COV; // delete GATTER, since non-negativity ensured
  if (!next->normalmix) ERR("sub model is not a normal mixture model");
  if (CovList[next->nr].spectral == NULL)
    ERR("submodel does not have spectral representation");
//  updatepref(cov, next); ## gute idee?
  if (next->pref[SpectralTBM] == PREF_NONE) cov->pref[RandomCoin] = PREF_NONE;
  return NOERROR;
}

int checkave1(cov_model *cov) {
  return checkave(cov, cov->tsdim - 1, "d-1");
}
int checkave2(cov_model *cov) {
  return checkave(cov, cov->tsdim, "d");
}

void ave(double *h, cov_model *cov, double cstart, int dim, double *v) {
  // f = uAu +zu; M = id; 
  cov_model *next = cov->sub[0];
  int i,j,k,d;
    // dimP1 = dim + 1;
  double  detMB, Ah[AveMaxDim],MB[AveMaxDim], Dh[AveMaxDim],
    dummy, c,
    hMh,
    *A = cov->p[0],
    *z = cov->p[1];
 
  c = cstart;
  hMh = 0.0;
  for (k=d=0; d<dim; d++) {
    for (dummy = 0.0, j=0; j<dim; j++) dummy += h[j] * A[k++];
    Ah[d] = dummy;
    c += (Ah[d] + z[d]) * h[d];
    hMh += h[d] * h[d]; 
  }

  for (j=d=0; d<dim; d++) {
    Dh[d] = h[d] /* Mh */ + 2.0 * c * Ah[d]; // 
    for (i=0; i<dim; i++) {
      MB[j++] = 2.0 * Ah[d] * Ah[i];
    }
    MB[j - dim + d] += 1.0;
  }

  det_UpperInv(MB, &detMB, dim);

//  double y = sqrt(c * c - hMh +  0.5 * x_UxPz(Dh, MB, z, dim));
  double y = sqrt(c * c + hMh - 0.5 * xUx(Dh, MB, dim));
  CovList[next->nr].cov(&y, next, v);
  *v /= sqrt(detMB);
}

void ave1(double *x, cov_model *cov, double *v) {
    // example 13
  int dim = cov->tsdim - 1;
  ave(x, cov, x[dim] /* t */, dim, v);
}

void ave2(double *x, cov_model *cov, double *v) {
    // no counterpart exists
  int dim = cov->tsdim;
  ave(x, cov, 0.0 /* t */, dim, v);
}


void sd_standard(mpp_storage *s, cov_model *cov){
  int d;
  double x2;
  for (x2=0.0, d=0; d<s->dim; d++) x2 += s->lensimu[d] * s->lensimu[d];
  x2 = sqrt(x2);
  for (d=0; d<s->dim; d++) s->sdgauss[d] = x2;
}

void sd_ave_stp(mpp_storage *s, cov_model *cov){
  int d;
  double b, alphamin, x2, InvSqrt2a,
    V = s->invscale;
  for (x2=0.0, d=0; d<s->dim; d++) x2 += s->lensimu[d] * s->lensimu[d];
  x2 *= 0.25;
  b = 3 * V * s->c[MPP_EIGEN] * x2 / s->dim;
  alphamin = (4.0 + 4.0 * b - 2.0 * sqrt(4 * b * b + 8.0 * b + 1)) / 3;
  InvSqrt2a = 1.0 / sqrt(2.0 * alphamin * 6 * V);
  for (d=0; d<s->dim; d++) s->sdgauss[d] = InvSqrt2a;
  
//  printf("%f\n", InvSqrt2a);assert(false);

}

cov_model aveGAUSS;
void mppinit_ave(mpp_storage *s, cov_model *cov, int dim) {

  assert(covOhneGatter(cov->sub[0])->tsdim == 1); 

  s->c[DRAWMIX_EXPONENT] = 0.25 * dim;
  COV_NULL(&aveGAUSS);
  aveGAUSS.nr = GAUSS;
  aveGAUSS.tsdim = 1;
  CovList[aveGAUSS.nr].initspectral(&aveGAUSS);
  // s->rsq = -0.5 * s->logapproxzero; //  -0.5 * log(approxzero)
  s->rsq = - s->logapproxzero; //  -0.5 * log(approxzero)
  // rsq: radius^2 ab dem Gaussglocke auf 0 gesetzt wird 
  s->r = s->effectiveRadius = sqrt(s->rsq);

  s->integral   = gaussInt(dim, 1, 1.0, s->effectiveRadius); // 1. Moment
  s->integralsq = gaussInt(dim, 2, 1.0, s->effectiveRadius); // 2. Moment
}
void mppinit_ave1(mpp_storage *s, cov_model *cov) {
  mppinit_ave(s, cov, cov->tsdim-1);
}
void mppinit_ave2(mpp_storage *s, cov_model *cov) {
  mppinit_ave(s, cov, cov->tsdim);
}


void coin_ave(mpp_storage *s, cov_model *cov){
  double spectral[AveMaxDim]; // do never reduce to 1 dimension
  // since CovList.spectral might use all the components 
  static spectral_storage spec = {0.0, 0.0, false, false}; 
  cov_model *next = covOhneGatter(cov->sub[0]);

  s->loginvscale = CovList[next->nr].drawlogmix(next, s);  // V
  s->invscale = exp(s->loginvscale);
  s->c[MPP_EIGEN] = 1.0;

  // printf(" %ld\n", s->location);
  s->location(s, cov); // uniform oder Gauss
  // assert(false);

  
  CovList[GAUSS].spectral(&aveGAUSS, &spec, spectral); // dim == 1 !!!
  s->c[AVE_SPECTRAL] = spectral[0] * sqrt(s->invscale); // dim == 1 !!
  s->c[AVE_VV] = TWOPI * UNIFORM_RANDOM;


}

void ave1_logg(double *u, cov_model *cov, int dim, double *logg, double *sign) {
  // smoothing kernel --- complex call since g could have negative sign
  // on the other hand, for the intensity lambda of the spatial vector,
  // lambda^{-1/2} could be very large causing numerical problems if
  // multivplied with very small g
  int i;
  double r2 = 0.0;

  *logg = *sign = 1.0;
  return;


  // print("%d logg\n", dim);

  for (i=0; i<dim; i++) r2 += u[i] * u[i];
//  print("l %d  logg\n", dim);
  *logg = -0.5 * dim * M_LN_SQRT_PId2 - r2;
//   print("s %d logg\n", dim);
 *sign = 1.0;
// printf("end logg %d %f %f %f\n", dim, *sign, *logg, r2);
}

///double ave1_g(double *u, cov_model *cov, int dim) { // smoothing kernel
//  int i;
//  double r2 = 0.0;
//  for (i=0; i<dim; i++) r2 += u[i] * u[i];
//  return exp(-0.5 * dim * M_LN_SQRT_PId2 - r2);
//}

double ave1_f(double *u, cov_model *cov, int dim) { // f within Y
  int d, j, k;
  double dummy,
    *A = cov->p[0],
    *z = cov->p[1],
    value;
  value = 0.0;
  for (k=d=0; d<dim; d++) {
    for (dummy = z[d], j=0; j<dim; j++) dummy += u[j] * A[k++];
    value += dummy * u[d];
  }
  return value;
}

res_type mppget_ave(double *dux, cov_model *cov, 
		  mpp_storage *s, int dim, double t) {
    // nur stationaer, nicht 
   double f,  logg, sign; // dux[AveMaxDim],  *u = s->u;

  // for (d=0; d < dim; d++) dux[d] = x[d] - u[d];
  ave1_logg(dux, cov, dim, &logg, &sign); 
  f = ave1_f(dux, cov, dim);

  return (res_type) 
      (exp(0.25 * dim * s->loginvscale + s->logInvSqrtDens + logg)
       * M_SQRT2 * cos(s->c[AVE_VV] + s->c[AVE_SPECTRAL] * (f - t))); // Y
}

res_type mppget_ave1(double *x, cov_model *cov, mpp_storage *s) {
  return mppget_ave(x, cov, s, s->dim - 1, x[s->dim-1]);
}

res_type mppget_ave2(double *x, cov_model *cov, mpp_storage *s) {
  return mppget_ave(x, cov, s, s->dim, 0.0);
}


/* coxgauss, cmp with nsst1 !! */
// C = 2 (C + 4 M H M), H = h h^t
// a = t - h M h - zh
// exp(- 0.5 * (h *h + 2 a^2 - mu C mu)) // stimmen die Vorzeichen??
// mu = h - 2 a M h
/* cox, cmp with nsst1 !! */
// coxisham

void GetEu2Dinv(param_type p, double *x, int dim, 
		double *det, double *Eu2Dinv,
		double *newxsq, double *newx, double *z) {
    double t, t2,
	y[AveMaxDim],
	*V = p[0],
	*D= p[1];
    int d,
	dimP1 = dim + 1,
	dimsq = dim * dim;
  t = x[dim];
  t2 = pow(fabs(t), p[2][0]); // standard t^2
  for (d=0; d<dim; d++) {
      y[d] = x[d] - t * V[d];
  }
  
  for (d=0; d<dimsq; d++) {
      Eu2Dinv[d] = t2 * D[d];
  }
  for (d=0; d<dimsq; d+=dimP1)  Eu2Dinv[d] += 1.0; // D + E
  det_UpperInv(Eu2Dinv, det, dim);
//  printf("t=%f, %f %f\n", t, t2, *det);
  *newxsq = xUxz(y, Eu2Dinv, dim, z);
  *newx = sqrt(*newxsq);
}

void cpyUf(double *Eu2Dinv, double factor, int dim, int tsdim, double *v) {
    // Eu2Dinv has dimension dim^2; v dimension tsdim^2
    // Eu2Dinv is copied into the upper left part of v and 
    // multiplied by factor
  int d, i, k, j,
	tsdimsq = tsdim * tsdim;
  
  for (i=0; i<tsdimsq; v[i++] = 0.0);
  for (i=0; i<dim; i++) {
      for (d=i * tsdim, k=i * dim, j=0; j<=i; j++)
	  v[d++] = Eu2Dinv[k++] * factor; 
      for (k += dim - 1; j<dim; j++, k+=dim) { 
	  v[d++] = Eu2Dinv[k] * factor;
      }
  }
}

void addzzT(double *v, double factor, double *z, int dim, int tsdim) {
    // z has dimension dim; v dimension tsdim^2
    // zzT is copied into the upper left part of v after being 
    // multiplied by factor
   
    int i,j,k;
    for (i=0; i<dim; i++) {
	k = i * tsdim;
	for (j=0; j<dim; j++) {
	    v[k++] += z[i] * z[j] * factor;
	}
    }
}


void kappa_cox1(int i, cov_model *cov, int *nr, int *nc){
    switch (i) {
	case 0 :
	    *nc = 1;
	    *nr = cov->tsdim - 1;
	    break;
	case 1 :
	    *nc = *nr = cov->tsdim - 1;
	    break;
	case 2 :
	    *nc = *nr = 1;
	    break;
	default:  *nc = *nr = -1;
    }
}

void cox1(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int dim = cov->tsdim - 1,
      dimsq = dim * dim;
  double det, newx, *Eu2Dinv, newxsq;
 
  Eu2Dinv = (double*) malloc(sizeof(double) * dimsq);
  GetEu2Dinv(cov->p, x, dim, &det, Eu2Dinv, &newxsq, &newx, NULL);
   
  CovList[next->nr].cov(&newx, next, v);
  *v /= sqrt(det);

  free(Eu2Dinv);
}

void coxhess(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int tsdim = cov->tsdim,
      dim = tsdim - 1,
      dimsq = dim * dim;
  double z[AveMaxDim], det, *Eu2Dinv, newx, newxsq, phiD, phiD2;

  Eu2Dinv = (double*) malloc(sizeof(double) * dimsq);
  GetEu2Dinv(cov->p, x, dim, &det, Eu2Dinv, &newxsq, &newx, z);

  CovList[next->nr].D2(&newx, next, &phiD2);  
  if (newxsq == 0.0) {
      cpyUf(Eu2Dinv, phiD2 / sqrt(det), dim, tsdim, v);
      //      printf("%f %f %f\n %f %f %f \n %f %f %f\n Eu2Dinv %f %f %f %f t=%f\n", 
//	 v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], 
//	      Eu2Dinv[0],Eu2Dinv[1],Eu2Dinv[2],Eu2Dinv[3],x[dim]  );
//       assert(Eu2Dinv[0] ==1 );

  } else {
      CovList[next->nr].D(&newx, next, &phiD);
      cpyUf(Eu2Dinv, phiD / (sqrt(det) * newx), dim, tsdim, v);
      addzzT(v, (phiD2 - phiD/newx) / (sqrt(det) * newxsq), z, dim, tsdim);
  }

  free(Eu2Dinv);
}


void coxnabla(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int d,
      tsdim = cov->tsdim,
      dim = tsdim - 1,
      dimsq=dim * dim;
  double z[AveMaxDim], det, newx, newxsq, *Eu2Dinv, phiD, factor;
  
  Eu2Dinv = (double*) malloc(sizeof(double) * dimsq);
  GetEu2Dinv(cov->p, x, dim, &det, Eu2Dinv, &newxsq, &newx, z); 

  if (newxsq == 0.0) {
      for (d=0; d<=dim; d++)  v[d] = 0.0;
  } else {    
    newx = sqrt(newxsq);
    CovList[next->nr].D(&newx, next, &phiD);
    factor = phiD / (det * newx); 
    for (d=0; d<dim; d++) {
	v[d] = factor * z[d];
    }
    for (d=0; d<tsdim; v[d++]=0.0);
  }

  free(Eu2Dinv);
}



int checkcox1(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err, 
       dim = cov->tsdim - 1; 

  cov->manipulating_x=true;
  if (cov->ncol[0] != 1 || cov->nrow[0] != dim) {
//    char Msg[255];
    sprintf(ERRORSTRING, 
	    "V is not given or not a vector of dimension 2 (nrow1=%d ncol2=%d dim=%d)",  
	    cov->nrow[0], cov->ncol[1], dim);
    return ERRORM;
  }

  // is matrix positiv definite?
  if (!is_positive_definite(cov->p[1], dim))
      ERR("D is not (strictly) positive definite.");

  kdefault(cov, 2, 2.0);
  if ((err = checkkappas(cov)) != NOERROR) return err;
  if ((err = check2X(next, dim, 1, STATIONARY, ISOTROPIC,
		     UNIVARIATE)) != NOERROR)
    return err;
  next->tsdim = DEL_COV; // delete GATTER, since non-negativity ensured
  if (!next->normalmix) ERR("submodel is not a normal mixture model");
  if (CovList[next->nr].spectral == NULL)
    ERR("submodel does not have spectral representation");

  updatepref(cov, next);
  if (cov->p[2][0] != 2.0) cov->pref[SpectralTBM] = 0;


  cov->hess = true;

	 
  return NOERROR;
}

void rangecox1(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;

  range->min[0] = RF_NEGINF;
  range->max[0] = RF_INF;
  range->pmin[0] = -100.0;
  range->pmax[0] = +100.0;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = RF_NEGINF;
  range->max[1] = RF_INF;
  range->pmin[1] = -100.0;
  range->pmax[1] = +100.0;
  range->openmin[1] = false;
  range->openmax[1] = false;  

  range->min[2] = 0.0;
  range->max[2] = 2.0;
  range->pmin[2] = 0.1;
  range->pmax[2] = 2.0;
  range->openmin[2] = true;
  range->openmax[2] = false;  
 
  range->maxdim = CoxMaxDim;
}
int initspectralcox1(cov_model *cov) {
  cov_model *next = cov->sub[0];
  if (cov->tsdim != 3) return ERRORFAILED;
  return CovList[next->nr].initspectral(next);
}

void spectralcox1(cov_model *cov, spectral_storage *s, double *e) { 
  cov_model *next = cov->sub[0];
  int d,
    dim = cov->tsdim - 1;
  double t, v[CoxMaxDim],
    *V = cov->p[0],
    rho= cov->p[1][0];
  CovList[next->nr].spectral(next, s, e);
  
  v[0] = rnorm(0.0, INVSQRTTWO);
  v[1] = rho * v[0] + sqrt(1 - rho * rho) * rnorm(0.0, INVSQRTTWO);
 
  for (t = 0.0, d=0; d<dim; d++) {
    t += (v[d] + V[d]) * e[d];
  }
  e[dim] = -t;
}




// GenPS (generalisation of paciore and stein)
// Q = (x-y) Sy (Sx + Sy)^{-1} Sx (x-y) (weicht etwas von Stein ab)
void kappa_stp(int i, cov_model *cov, int *nr, int *nc){
  *nc = (i <= 1) ? cov->tsdim : 1;
  *nr = (i <= 2) ? cov->tsdim : -1;
}
void stp(double *x,  double *y, cov_model *cov, double *v) {
  int d, j, k,
      dim =cov->tsdim,
      dimsq = dim * dim;
  double h[StpMaxDim], 
    Mh[StpMaxDim], hSx[StpMaxDim],
    Syh[StpMaxDim], xi2x, xi2y, 
    detA, hMh, cxy, zh, Q, Amux[StpMaxDim], Amuy[StpMaxDim],
    // Q2, Q3, 
      *Sx, *Sy, *A,
    *Sc = cov->p[0],
    *M = cov->p[1],
    *z = cov->p[2];
  cov_model *next = cov->sub[0],
    *Sf = cov->sub[1],
    *xi2 =cov->sub[2],
    *H =cov->sub[3];

  Sx = (double*) malloc(dimsq * sizeof(double));
  Sy = (double*) malloc(dimsq * sizeof(double));
  A  = (double*) malloc(dimsq * sizeof(double));

  if (Sf != NULL) {
    covfct fct = CovList[Sf->nr].cov;
    fct(x, Sf, Sx); // symmetric, pos definite !!
    fct(y, Sf, Sy);
  } else {
    int bytes = sizeof(double) * dimsq;
    memcpy(Sx, Sc, bytes);
    memcpy(Sy, Sc, bytes);
  }

  if (xi2 != NULL) {
    covfct fct = CovList[xi2->nr].cov;
    fct(x, xi2, &xi2x);
    fct(y, xi2, &xi2y);
  } else {
    xi2x = xi2y = 0.0;
  }

  if (H == NULL) {
    for (k=0, d=0; d<dim; d++) {
      h[d] = x[d] - y[d];
      // printf("%d: %f %f\n", d, x[d], y[d]);
    }
  } else {
    covfct fct = CovList[H->nr].cov;
    double xx[StpMaxDim], yy[StpMaxDim];
    fct(x, H, xx);
    fct(y, H, yy);
    for (k=0, d=0; d<dim; d++) {
      h[d] = xx[d] - yy[d];
      // printf("%d: %f %f\n", d, x[d], y[d]);
    }
  }

  //Q2 = Q3 = 
    zh = hMh = 0.0;
  for (k=0, d=0; d<dim; d++) {
    Mh[d] = hSx[d] = Syh[d] = 0.0;
    for (j=0; j<dim; j++, k++) {
     Mh[d] += h[j] * M[k];
     hSx[d] += h[j] * Sx[k];
     Syh[d] += h[j] * Sy[k]; // uses that S is symmetric
    }
    zh += z[d] * h[d];
    hMh += Mh[d] * h[d];
    //Q2 += Syh[d] * h[d];
    //Q3 += hSx[d] * h[d];
   }
  cxy = xi2x - xi2y - zh;
  //Q2 += (hMh - cxy) * (hMh - cxy); 
  //Q3 += (hMh + cxy) * (hMh + cxy); 

  // printf("%f %f\N", 
  
  for (k=d=0; d<dim; d++) {
    for (j=0; j<dim; j++, k++) {
      A[k] = Sx[k] + Sy[k] + 4.0 * Mh[d] * Mh[j];
    }
    Amux[d] = hSx[d] + 2.0 * (hMh + cxy) * Mh[d]; // uses that M is symmetric
    Amuy[d] = Syh[d] + 2.0 * (hMh - cxy) * Mh[d];
  }

  // printf("Sx=%f %f %f %f\n", Sx[0], Sx[1],  Sx[2],Sx[3]);
  // printf("Sy=%f %f %f %f\n", Sy[0], Sy[1],  Sy[2],Sy[3]);
  // printf("A=%f %f %f %f\n", A[0], A[1],  A[2],A[3]);
//  double AA[StpMaxDim SQ]; memcpy(AA, A, sizeof(double) * StpMaxDim SQ); printf("%f\n", detU(AA, dim));
  det_UpperInv(A, &detA, dim);
  //printf("A=%f %f %f\n", A[0], A[2], A[3]);
//  printf("Ainv=%f %f %f detA=%f %d\n", A[0], A[2], A[3], detA, dim);

  //printf("Q2=%f\n", Q2);
  //Q2 -= xUy(Amuy, A, Amuy, dim);
  // Q3 -= xUy(Amux, A, Amux, dim);
  // printf("Q2=%f %f\n", Q2, xUy(Amuy, A, Amuy, dim));

  Q = cxy * cxy - hMh * hMh + xUy(Amux, A, Amuy, dim);

  //printf("Q=%f Q2=%f Q3=%f (Q2+Q3)/2=%f\n", Q, Q2, Q3, 0.5 *(Q2 + Q3));

  if (Q < 0.0) {
    PRINTF("x=%f,%f y=%f,%f detA=%f\n", 
	   x[0], x[1], y[0], y[1], detA);
    PRINTF("cxy=%4f hMh=%f Amux=%f A[0]=%f Amuy=%f\ndim=%d h=%f,%f hSx=%f,%f, xUy=%f Q=%f\n", 
	  cxy, hMh, Amux[0], A[0], Amuy[0], 
	   dim, h[0], h[1], hSx[0], hSx[1], xUy(Amux, A, Amuy, dim), Q);
    assert(Q >= 0.0);
  }

  Q = sqrt(Q);

  aux_covfct auxcf;
  if ((auxcf = CovList[next->nr].aux_cov) != NULL) 
    auxcf(x, y, Q, next, v);
  else 
    CovList[next->nr].cov(&Q, next, v);

  
  double 
    dx = detU(Sx, dim), 
    dy = detU(Sy, dim);
  
  *v *=  pow(2.0, 0.5 * double(dim)) * pow(dx * dy / (detA * detA), 0.25);
  free(Sx);
  free(Sy);
  free(A);
}



int checkstp(cov_model *cov){
  cov_model 
    *next = cov->sub[0],
    *Sf = cov->sub[1],
    *xi2 =cov->sub[2],
    *H = cov->sub[3];
  int err;
  
  int dim = cov->tsdim;
  
//  printf("%d %d\n", cov->ncol[0], cov->p[0]);
//  assert(false);
  if (cov->xdim > StpMaxDim) {
      sprintf(ERRORSTRING, 
	      "For technical reasons max. dimension for ave is %d. Got %d.", 
	      StpMaxDim, cov->xdim);
    return ERRORMSG;
  }

  cov->manipulating_x=true;
  if (cov->ncol[0] == 0) {  // Sc
    EinheitsMatrix(cov->p, dim);
    cov->ncol[0] = cov->nrow[0] = dim;
  }
  if (cov->ncol[1] == 0) { // M
    EinheitsMatrix(cov->p, dim);
    cov->ncol[1] = cov->nrow[1] = dim;
  }
  if (cov->ncol[2] == 0) { // z
    cov->p[2] = (double*) calloc(dim, sizeof(double));
    cov->ncol[2] = 1;
    cov->nrow[2] = dim;
  }

//  printf("%d\n", dim); assert(false);

  if ((err = checkkappas(cov)) != NOERROR) return err;

  if ((err = check2X(next, 1, 1, COVARIANCE, ANISOTROPIC,
		     UNIVARIATE)) != NOERROR)
    return err;
  if (!next->normalmix) ERR("submodel must be a normal scale mixture model");

  if (Sf != NULL) {
    if ((err = check2X(Sf, dim, dim, AUXMATRIX, ANISOTROPIC,
		       dim)) != NOERROR) 
      return err;
  }
  

  if (xi2 != NULL) {
   if ((err = check2X(xi2, dim, dim, AUXMATRIX, ANISOTROPIC,
		     UNIVARIATE)) != NOERROR)
     return err;
  }

  if (H != NULL) {
   if ((err = check2X(H, dim, dim, AUXVECTOR, ANISOTROPIC, dim)) != NOERROR)
     return err;
  }
  return NOERROR;
}

void rangestp(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;
  int i;

  for (i=0; i<=2; i++) { /* S, M */
    range->min[i] = RF_NEGINF;
    range->max[i]= RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}

cov_model stpGAUSS;
void mppinit_stp(mpp_storage *s, cov_model *cov) {
//  printf("%d %s\n", cov->p[0], CovList[cov->nr].name);
//  assert(false);
  cov_model *Sf = cov->sub[1];

//  warning("dadurch dass U von V abhaengt muss eventuell nachkorrigiert werden");

  s->effectiveRadius = 0.0;
  s->effectivearea = 1.0;

  s->integralpos = RF_NAN;
  s->integral=0.0;  
  s->integralsq=1.0; // falsch!
  s->maxheight = RF_NAN;
  if (Sf != NULL) {
    realfct eigen = CovList[Sf->nr].mineigenvalue;
    assert(eigen != NULL);
    s->c[MPP_EIGEN] = eigen(Sf);
  } else {
    double value[StpMaxDim], ivalue[StpMaxDim], dummy[5 * StpMaxDim], min = RF_INF;
    int i, err,       
      dim = cov->tsdim,
      ndummy = 5 * dim;
    //   char No = 'N';

    //   printf("%d %f %f\n", dim, cov->p[0][0], cov->p[0][1]);

    F77_NAME(dgeev)("No", "No", &dim, cov->p[0], &dim, 
		    value, ivalue, NULL, &dim, NULL, &dim,
		    dummy, &ndummy, &err);
    for (i=0; i<dim; i++) {
      if (min > value[i]) min = value[i];
    }
    s->c[MPP_EIGEN] = min;
  }
  s->c[DRAWMIX_EXPONENT] = 0.25 * cov->tsdim;

  COV_NULL(&stpGAUSS);
  stpGAUSS.nr = GAUSS;
  stpGAUSS.tsdim = 1;
  CovList[stpGAUSS.nr].initspectral(&stpGAUSS);
}

void coin_stp(mpp_storage *s, cov_model *cov){
  cov_model *next = cov->sub[0];
  static spectral_storage spec = {0.0, 0.0, false, false}; 
  double spectral[StpMaxDim];
  
  s->loginvscale = CovList[next->nr].drawlogmix(next, s);  // V
  s->invscale = exp(s->loginvscale);
  s->location(s, cov);


  // PrintModelInfo(&stpGAUSS);
  
  CovList[GAUSS].spectral(&stpGAUSS, &spec, spectral); 
  // return;
  s->c[MPP_SPECTRAL] = spectral[0] * sqrt(s->invscale);  
  s->c[MPP_VV] = TWOPI * UNIFORM_RANDOM;

  if (!(s->c[MPP_SPECTRAL] >= -1e40 && s->c[MPP_SPECTRAL] <1e40)) {
    PRINTF("%f %f %f %f %f\n",
	   spectral,
	   sqrt(s->invscale),
	   s->loginvscale,
	   s->c[DRAWMIX_ALPHA],
	   s->c[MPP_SPECTRAL]);
    assert(false);
  }

  s->effectiveRadius = s->plus + s->relplus / s->invscale;

  // printf("coin: %f %f %f\n", s->u[0], s->u[1], s->logInvSqrtDens);
}

//double *x, double *u, double logInvSqrtDens, double *p, 
//		    double *r, logmixtureweight logmixweight
res_type mppget_stp(double *x, double*u, cov_model *cov, mpp_storage *s){
  // kann um ca. Faktor 2 beschleunigt werden, wenn
  // Sx , logdetU, Hx fuer alle x abgespeichert werden
  // und die Werte dann aus dem Speicher gelesen werden
  // jedoch sehr Speicherintensiv. memcpy braucht man auch nicht

  cov_model *next = cov->sub[0],
    *Sf = cov->sub[1],
    *xi2 =cov->sub[2],
    *H =cov->sub[3];
  int j, k, d, 
      dim= s->dim,
      bytes = sizeof(double) * dim * dim;
  double h[StpMaxDim], hSxh, hSx, xi, Mhd, *Sx,
      //  *u = s->u,
    *Sc = cov->p[0],
    *M = cov->p[1],
    *z = cov->p[2];
  
  Sx= (double*) malloc(bytes);
  if (Sf == NULL) {
    memcpy(Sx, Sc, bytes);
  } else {
    CovList[Sf->nr].cov(x, Sf, Sx); // symmetric, pos definite !!
  }

  if (xi2 == NULL) {
    xi = 0.0;
  } else {
    CovList[xi2->nr].cov(x, xi2, &xi);
  }

  if (H == NULL) {
    for (k=0, d=0; d<dim; d++) {
      h[d] = u[d] - x[d];
      // printf("%d: %f %f\n", d, x[d], y[d]);
    }
  } else {
    covfct fct = CovList[H->nr].cov;
    double xx[StpMaxDim];
    fct(x, H, xx);
    for (k=0, d=0; d<dim; d++) {
      h[d] = u[d] - xx[d];
      // printf("%d: %f %f\n", d, x[d], y[d]);
    }
  }

  hSxh = 0.0;
  for (k=0, d=0; d<dim; d++) {
    Mhd = hSx = 0.0;
    for (j=0; j<dim; j++) {
     Mhd += h[j] * M[k];
     hSx += h[j] * Sx[k++];
     
    }
    xi += Mhd * h[d] + z[d] * h[d];
    hSxh += hSx * h[d];
  }

  double logdetU = log(detU(Sx, dim));
  double exponent =
    0.25 * dim * (2.0 + s->loginvscale - 2.0 * M_LN_SQRT_PI) /* ( 2V/pi)^{d/4} */
    + 0.25 *logdetU                          /* Sx ^1/4 */
    - s->invscale * hSxh             /* exp(-V(U-x) S (U-x)) */
    + CovList[next->nr].logmixweight(x, next, s)/* g*/
    + s->logInvSqrtDens;                /* 1 / sqrt(f) */

  if (!(exponent < 5.0) && PL > 4) {
    if (!(exponent < 6.0)) // could be NA, too
      PRINTF("\n%f logDetU=%f %f lgmx=%f %f alph=%f z=%f",
	     0.25 * dim * (2.0 + s->loginvscale - 2*M_LN_PId2) /* 2 V/pi)^{d/2} */
	     , 0.25 * logdetU                          /* Sx ^1/4 */
	     , -s->invscale * hSxh             /* exp(-V(U-x) S (U-x)) */
	     , CovList[next->nr].logmixweight(x, next, s)/* g*/
	     , s->logInvSqrtDens, 
	     s->c[DRAWMIX_ALPHA]
	     , exponent);
    else PRINTF("!");
  }

  assert(exp(exponent) < 10000000.0);
     
  free(Sx);
  return (res_type)
      (exp(exponent) * cos(s->c[MPP_VV] + s->c[MPP_SPECTRAL] * xi));  /* Y */
}



/* Whittle-Matern or Whittle or Besset */ 
void kappa_rational(int i, cov_model *cov, int *nr, int *nc){
  *nc = (i == 0) ? cov->tsdim : 1;
  *nr = (i==0) ? cov->tsdim : (i==1) ? 2 : -1;
}
double minEigenrational(cov_model *cov) {
  double *a = cov->p[1];
  return (a[0] < a[1]) ? a[0] : a[1];
}
void rational(double *x, cov_model *cov, double *v) {
  int i, k, j, 
    dim = cov->tsdim;
  double nu,
    *A = cov->p[0],
    *a = cov->p[1];
  nu = 0.0;
  for (k=0, i=0; i<dim; i++) {
    double xTC;
    xTC =  0.0;
    for (j=0; j<dim; j++) {
      xTC += x[j] * A[k++];
    }
    nu += xTC * xTC;
  }
  *v = (a[0] + a[1] * nu) / (1 + nu);
}
 
int checkrational(cov_model *cov){
//  pref_type pref = 
//    {5, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 5};
  // CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  if (cov->nrow[1] == 1) {
    double dummy = cov->p[1][0];
    free(cov->p[1]);
    cov->p[1] = (double *) malloc(sizeof(double) * 2);
    cov->p[1][0] = dummy;
    cov->p[1][1] = 0.0;
    cov->nrow[1] = 2;
  }
  return NOERROR;
}

void rangerational(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;

  range->min[0] = RF_NEGINF;
  range->max[0] = RF_INF;
  range->pmin[0] = -1e10;
  range->pmax[0] = 1e10;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = 0.0;
  range->pmax[1] = 10;
  range->openmin[1] = false;
  range->openmax[1] = true;
}


// Sigma(x) = diag>0 + A'xx'A
void kappa_EAxxA(int i, cov_model *cov, int *nr, int *nc){
  *nc = (i == 1) ? cov->tsdim : 1;
  *nr = (i <= 1) ? cov->tsdim : -1;
}
void EAxxA(double *x, cov_model *cov, double *v) {
  int d, k, j, 
    dim = cov->tsdim;
  double xA[EaxxaMaxDim],
    *E = cov->p[0],
    *A = cov->p[1];
  for (k=0, d=0; d<dim; d++) {
    xA[d] =  0.0;
    for (j=0; j<dim; j++) {
      xA[d] += x[j] * A[k++];
    }
  }
  for (k=d=0; d<dim; d++) {
    double xAd = xA[d];
    for (j=0; j<=d; j++) {
      v[k++] = xAd * xA[j];
    }
    v[k-1] += E[d];
    for ( ; j<dim; j++) {
      v[k++] = xAd * xA[j];
    }
  }
}

double minEigenEAxxA(cov_model *cov) {
  double min,
    *E = cov->p[0];
  int i,
    dim = cov->tsdim;
  for (min = RF_INF, i=0; i<dim; i++)
    if (E[i] < min) min = E[i];
  return min;
}
 
int checkEAxxA(cov_model *cov){
//  pref_type pref = 
//    {5, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 5};
  // CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
//  memcpy(cov->pref, pref, sizeof(pref_type));  
    if (cov->xdim > EaxxaMaxDim) {
	sprintf(ERRORSTRING, 
	      "For technical reasons max. dimension for ave is %d. Got %d.", 
	      EaxxaMaxDim, cov->xdim);
    return ERRORMSG;
  }

  cov->vdim = cov->tsdim;
  return NOERROR;
}
 
void rangeEAxxA(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;

  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 0.0001;
  range->pmax[0] = 10;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = RF_NEGINF;
  range->max[1] = RF_INF;
  range->pmin[1] = -1e10;
  range->pmax[1] = 1e10;
  range->openmin[1] = true;
  range->openmax[1] = true;
}




// Sigma(x) = diag>0 + A'xx'A
void kappa_EtAxxA(int i, cov_model *cov, int *nr, int *nc){
  *nc = (i == 1) ? cov->tsdim : 1;
  *nr = (i <= 1) ? cov->tsdim : (i==2) ? 1 : -1;
}
void EtAxxA(double *x, cov_model *cov, double *v) {
  int d, k, j, 
    dim = cov->tsdim,
    time = dim - 1;
  double xDR[EaxxaMaxDim],
    *E = cov->p[0],
    *D = cov->p[1],
    phi = cov->p[2][0],
    c =  cos(phi * x[time]),
    s = sin(phi * x[time]),
    R[9]; assert(dim ==3);
   
//  PrintModelInfo(cov);

 
  R[0] = R[4] = c;
  R[1] = s;
  R[3] = -s;
  R[2] = R[5] = R[6] = R[7] = 0.0;
  R[8] = 1.0;
 
  {
    double xD[EaxxaMaxDim];
    for (k=0, d=0; d<dim; d++) {
      xD[d] =  0.0;
      for (j=0; j<dim; j++) {
	xD[d] += x[j] * D[k++];
      }
    }
    for (k=0, d=0; d<dim; d++) {
      xDR[d] =  0.0;
      for (j=0; j<dim; j++) {
	xDR[d] += xD[j] * R[k++];
      }
    }
    // printf("%f %f\n", c, s);
    // for (d=0; d<dim; d++) assert(xDR[d] == xD[d]);
  }


  for (k=d=0; d<dim; d++) {
    double xDd = xDR[d];
    for (j=0; j<=d; j++) {
      v[k++] = xDd * xDR[j];
    }
    v[k-1] += E[d]; // nur korrekt falls E Vielfaches der EH-Matrix
    for ( ; j<dim; j++) {
      v[k++] = xDd * xDR[j];
    }
  }

//  double w[MAX DIM];
//  EAxxA(x, cov, w);
  // for (d=0; d<dim; d++) {
    // printf("%d %f %f\n", d, w[d], v[d]);
  //  assert(w[d] == v[d]);
  //}

}

double minEigenEtAxxA(cov_model *cov) {
  double min,
    *E = cov->p[0];
  int i,
    dim = cov->tsdim;
  for (min = RF_INF, i=0; i<dim; i++)
    if (E[i] < min) min = E[i];
  return min;
}
 
int checkEtAxxA(cov_model *cov){
//  pref_type pref = 
//    {5, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 5};
  // CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
//  memcpy(cov->pref, pref, sizeof(pref_type));  
  cov->vdim = cov->tsdim;
  return NOERROR;
}
 
void rangeEtAxxA(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;
  int i;

  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 0.0001;
  range->pmax[0] = 10;
  range->openmin[0] = true;
  range->openmax[0] = true;

  for (i=1; i<=2; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}




// Sigma(x) = diag>0 + A'xx'A
void kappa_rotat(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = (i <= 1) ? 1 : -1;
}
void rotat(double *x, cov_model *cov, double *v) {
  int
    dim = cov->tsdim,
    time = dim - 1;
  double
    speed = cov->p[0][0],
    phi = cov->p[0][1],
    absx = sqrt(x[0] * x[0] + x[1] * x[1]);
  *v = (absx == 0.0) ? 0.0
    : speed * (cos(phi * x[time]) * x[0] + sin(phi * x[time]) * x[1]) / absx;
  // printf("%f\n", *v);
}

double minEigenrotat(cov_model *cov) {
  return 0.0;
}
 
int checkrotat(cov_model *cov){
//  if (cov->tsdim != 3) ERR("only 3-d allowed for rotat!");
//  pref_type pref = 
//    {5, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 5};
  // CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
//  memcpy(cov->pref, pref, sizeof(pref_type));  
  cov->vdim = cov->tsdim - 1;
  return NOERROR;
}
 
void rangerotat(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;
  int i;

  for (i=0; i<2; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}



// Sigma(x) = diag>0 + A'xx'A
void kappa_Rotat(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = (i == 0) ?  1 : -1;
}
void Rotat(double *x, cov_model *cov, double *v) {
  int d, k, j, 
    dim = cov->tsdim,
    time = dim - 1;
  double
      phi = cov->p[0][0],
      c =  cos(phi * x[time]),
      s = sin(phi * x[time]),
      R[9]; assert(dim ==3);
   
  R[0] = R[4] = c;
  R[1] = s;
  R[3] = -s;
  R[2] = R[5] = R[6] = R[7] = 0.0;
  R[8] = 1.0;
 
  for (k=0, d=0; d<dim; d++) {
    v[d] =  0.0;
    for (j=0; j<dim; j++) {
      v[d] += x[j] * R[k++];
    }
  }
} 
int checkRotat(cov_model *cov){
//  pref_type pref = 
//    {5, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 5};
  // CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
//  memcpy(cov->pref, pref, sizeof(pref_type));  
  cov->vdim = cov->tsdim;
  return NOERROR;
}

void rangeRotat(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;

  range->min[0] = RF_NEGINF;
  range->max[0] = RF_INF;
  range->pmin[0] = -1e10;
  range->pmax[0] = 1e10;
  range->openmin[0] = true;
  range->openmax[0] = true;
}



/* Whittle-Matern or Whittle or Besset */ 
// statt 0 Parameter: 2 Parameter, M und z fuer xi
void kappaNonStWM(int i, cov_model *cov, int *nr, int *nc){
  *nc =  1;
  *nr = (i==0) ? 1 : -1;
}

void NonStWMQ(double *x, double *y, double sqrtQ, cov_model *cov, double *v){
// check calling functions, like hyperbolic and gneiting if any changings !!
  double loggamma, nuxy, nux, nuy;
  cov_model *nu = cov->sub[0];

  if (nu  == NULL) {
    nuxy = cov->p[0][0];
    loggamma = lgammafn(nuxy);
  } else {
    covfct Cnu = CovList[nu->nr].cov;
    Cnu(x, nu, &nux);
    Cnu(y, nu, &nuy);
    nuxy = 0.5 * (nux + nuy);
    loggamma = 0.5 * (lgammafn(nux) + lgammafn(nuy));
  }
    
  if (sqrtQ == 0.0) {
    *v = 1.0;
    return;
  }
 
  *v = 2.0 * exp(nuxy * log(0.5 * sqrtQ) 
		 - loggamma
		 + log(bessel_k(sqrtQ, nuxy, 2.0)) - sqrtQ);
}

void NonStWM(double *x, double *y, cov_model *cov, double *v){
// check calling functions, like hyperbolic and gneiting if any changings !!
  int d,
    dim = cov->tsdim;
  double Q=0.0;

  for (d=0; d<dim; d++) {
    double delta = x[d] - y[d];
    Q += delta * delta;
  }

  NonStWMQ(x, y, Q, cov, v);
}

int checkNonStWM(cov_model *cov) { 
  cov_model *nu = cov->sub[0];
//    *Q=cov->sub[1];
  int err,
    dim = cov->tsdim;

  cov->manipulating_x=true;
  if (cov->p[0] == NULL) {
    if (nu == NULL)  error("nu is missing");
    cov->p[0] = (double*) malloc(sizeof(double));
    cov->p[0][0] = 1.0; // never used, but avoids errors in checks
    cov->nrow[0] = cov->ncol[0] = 1;
  }

  if ((err = checkkappas(cov)) != NOERROR) return err;
  
  if (nu != NULL) {
    if ((err = check2X(nu, dim, dim, AUXMATRIX, ANISOTROPIC,
		     UNIVARIATE)) != NOERROR) 
      return err;
    if (nu->sub[0]->tsdim != cov->tsdim) ERR("submodel not allowed");
  }
  return NOERROR;
}

void rangeNonStWM(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-2;
  range->pmax[0] = 10.0;
  range->openmin[0] = true;
  range->openmax[0] = true;
}

// using nu^(-1-nu+a)/2 for g and v^-a e^{-1/4v} as density instead of frechet
// the bound 1/3 can be dropped
static double eM025 = exp(-0.25);
double DrawLogMixNonStWM(cov_model *cov, mpp_storage *s) { // inv scale
  // V ~ F in stp
  cov_model *nu = cov->sub[0];
  double
    minnu = (nu == NULL) ? cov->p[0][0] : CovList[nu->nr].mineigenvalue(nu),
    alpha = 1.0 + 0.5 /* 0< . < 1*/ *
    (3.0 * minnu - 2.0 * s->c[DRAWMIX_EXPONENT]);
  if (alpha > 2.0) alpha = 2.0; // original choice
  if (alpha <= 1.0) ERR("minimual nu too low or dimension too high");
  s->c[DRAWMIX_ALPHA] = alpha;
  double beta = s->c[DRAWMIX_BETA],
    p = s->c[DRAWMIX_P],
    logU;
  if (UNIFORM_RANDOM < p){
    s->c[DRAWMIX_ALPHA] = beta;
    logU =  log(UNIFORM_RANDOM * eM025);
    s->c[DRAWMIX_FACTOR] = -0.5 * log(0.25 * p * (beta - 1.0)) + 0.25;
  } else {
    s->c[DRAWMIX_ALPHA] = alpha;
    logU = log(eM025 + UNIFORM_RANDOM * (1.0 - eM025));
    s->c[DRAWMIX_FACTOR] =  -0.5 * log(0.25 * (1.0 - p) * (alpha - 1.0));
  } 
  return  log( -0.25 / logU) / (s->c[DRAWMIX_ALPHA] - 1.0); //=invscale
}

double LogMixWeightNonStWM(double *x, cov_model *cov, mpp_storage *s) {
  // g(v,x) in stp
  double nu,
    alpha = s->c[DRAWMIX_ALPHA],
    loginvscale = s->loginvscale;
  cov_model *Nu = cov->sub[0];
  static double storage = 0.0; 
  
  if (Nu == NULL) 
    nu = cov->p[0][0];
  else 
    CovList[Nu->nr].cov(x, Nu, &nu);

  double z;
  z = - nu  * M_LN2 // in g0  // eine 2 kuerzt sich raus
    + 0.5 * ((1.0 - nu) /*in g0*/ + alpha /*lambda*/ - 2.0 /*fre*/) * loginvscale
    - 0.5 * lgammafn(nu)  // in g0
    + s->c[DRAWMIX_FACTOR] // lambda
    - 0.125 / s->invscale   // g: Frechet
    + 0.125 * pow(s->invscale, 1.0 - alpha); // lambda: frechet

  if (!(z < 7.0)) {
    if (storage != loginvscale) {
      if (PL > 4) 
	PRINTF("alpha=%f, is=%f, cnst=%f logi=%f lgam=%f loga=%f invlogs=%f pow=%f z=%f\n",
	   alpha,s->invscale,
	   (1.0 - nu) * M_LN2 
	  , + ((1.0 - nu) * 0.5 + alpha - 2.0) * loginvscale
	   ,- 0.5 * lgammafn(nu) 
	  , -s->c[DRAWMIX_FACTOR]
	   ,- 0.25 / s->invscale 
	  , + 0.25 * pow(s->invscale, - alpha)
	  , z);
    storage = loginvscale;
    }
    //assert(z < 10.0);
  }

  return z;
				      
}




/* nsst */
/* Tilmann Gneiting's space time models, part I */
void nsst(double *x, cov_model *cov, double *v) {
  cov_model *next1 = cov->sub[0];
  cov_model *next2 = cov->sub[1];
  cov_fct *C1 = CovList + next1->nr,
    *C2 = CovList + next2->nr;
  double v1, v2, psi, phi, y;
  
  C2->cov(ZERO, next2, &v1);
  C2->cov(x + 1, next2, &v2);
  psi = sqrt(1.0 + v1 - v2);  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)

//  printf("%f %f v=%f %f %f p=%f \n", x[0],x[1],  v1, v2, psi, cov->p[0][0]);

  y = x[0] / psi;
  C1->cov(&y, next1, &phi);
  *v = pow(psi, -cov->p[0][0]) * phi;

  //printf("%f  %f;%f %f %f : %f\n", x[0], x[1], phi, psi, cov->p[0][0], *v);
}

void TBM2nsst(double *x, cov_model *cov, double *v) {
  cov_model *next1 = cov->sub[0];
  cov_model *next2 = cov->sub[1];
  cov_fct *C1 = CovList + next1->nr,
    *C2 = CovList + next2->nr;
  double v1, v2, psi, phi, y;

  C2->cov(ZERO, next2, &v1);
  C2->cov(x + 1, next2, &v2);
  psi = 1.0 + v1 - v2;  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  C1->tbm2(&y, next1, &phi);
  *v = pow(psi, -cov->p[0][0]) * phi;
}

void Dnsst(double *x, cov_model *cov, double *v) {
  cov_model *next1 = cov->sub[0];
  cov_model *next2 = cov->sub[1];
  cov_fct *C1 = CovList + next1->nr,
    *C2 = CovList + next2->nr;
  double v1, v2, psi, phi, y;

  C2->cov(ZERO, next2, &v1);
  C2->cov(x + 1, next2, &v2);
  psi = 1.0 + v1 - v2;  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  C1->D(&y, next1, &phi);
  *v = pow(psi, -cov->p[0][0] - 1) * phi;
}

int checknsst(cov_model *cov) {
  cov_model *next1 = cov->sub[0];
  cov_model *next2 = cov->sub[1];
  int err;
      
  cov->manipulating_x=true;
  if (cov->xdim != 2) {
    sprintf(ERRORSTRING, "%s", "reduced dimension must be 2");
    return ERRORM;
  }
  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->normalmix = cov->finiterange = false;
  next1->xdim = 1;
  if ((err = check2X(next1, cov->tsdim, 1, STATIONARY, ISOTROPIC, 
		     UNIVARIATE)) != NOERROR) 
    return err;
  
  if (!next1->normalmix) XERR(ERRORNORMALMIXTURE);
  setbackward(cov, next1);

  cov->normalmix = false;
  next2->tsdim = next2->xdim = 1;
  if ((err = check2X(next2, 1, 1, VARIOGRAM, ISOTROPIC, 
		     UNIVARIATE)) != NOERROR) 
    return err;
  // kein setbackward(cov, next2);
  return NOERROR;
}

void rangensst(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = cov->tsdim - 1;
  range->max[0] = RF_INF;
  range->pmin[0] = cov->tsdim - 1;
  range->pmax[0] = 10.0;
  range->openmin[0] = false;
  range->openmax[0] = true;
}


