/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Kriging methods

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


/*

  PRINTING LEVELS
  ===============
  error messages: 1
  forcation : 2
  minor tracing information : 3--5
  large debugging information: >10

*/

/*
  calculate the transformation of the points only once and store the result in
  a register (indexed by ncov) of key, not of s->. Advantages: points have to
  calculated only once for several tries; if zonal anisotropy then specific 
  methods may be called;  
*/

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
 
#include <string.h>
#include "RF.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R.h>
#include <Rdefines.h>

#define KRIGE_TOLERANCE -1e-10

void CalculateVariogram(double *x, int lx, cov_model *cov, double *vario);
void matmult(double *A, double *B, double * C, int l, int m, int n);
void matmulttransposed(double *A, double *B, double *C, int m, int l, int n);
void InvChol(double *C, int dim); // to replace quote(...solve) in kriging.R

void simpleKriging(double *tgiven, double *x, double *invcov, int *notna,
		   int *Nx, int *Ngiven,
		   int *Dim, int *Rep, double *Mean, double *Krig) {
  // kriging variance is not calculated
    cov_model *cov = STORED_MODEL[MODEL_INTERN];
    int  d, i, j, r, krigi, v,
	err = NOERROR,
	dim = *Dim,
	ngiven = *Ngiven,
	rep = *Rep,
	nx = *Nx,
	vdim = cov->vdim,
	vdimng = vdim * ngiven,
	len_tgiven = dim * ngiven,
	divachtzig = (nx<79) ? 1 : (nx / 79),
	divachtzigM1 = divachtzig - 1;
  
    double xx[MAXSIMUDIM], dummy,
     *M = NULL,
      *pM,
    *dist = NULL;
  // Ngiven : length given = length M vector to be calculated
  // Rep : repetition invcov
  // Nx : number of points in x, length(x) = nx * dim
  if ((M = (double*) malloc(sizeof(double) * vdimng * vdim))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }

  for (krigi=i=0; i<nx; i++) {
    // printf("%d %d %d \n", i, nx, divachtzig);
    if (PL>0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=nx) 
	xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim)
	dist[j] = tgiven[j] - xx[d];
    CovIntern(dist, NULL, true, ngiven, M);

    //   for (r=0; r<ngiven * vdim * vdim; r++) printf("%e ", M[r]); printf("\n");

    krigi = i;
    pM = M;
    for (r=0; r<rep; r++) {
      for (j=v=0; v<vdim; v++, krigi+=nx, pM+=vdimng) {
	dummy = Mean[v];
	for (d=0; d<vdimng; d++) {
	    if (notna[d]) dummy += pM[d] * invcov[d];
//	    printf("%d r=%d v=%d d=%d value=%f %f %f %d \n", 
//		   i, r, v, d, dummy, pM[d], invcov[d], notna[d]);
	}
	Krig[krigi] = dummy;  // krigi=i + v * nx + r * nx * vdim
      }
    }

//    assert(false);
  }
  if (PL > 0) PRINTF("\n");

 ErrorHandling:
  if (M!=NULL) free(M);
  if (dist!=NULL) free(dist);

  if (err != NOERROR) {
      int nKrig;
      nKrig = nx * rep;
      for (i=0; i<nKrig; Krig[i++]=RF_NAN);
  }
}


void simpleKriging2(double *tgiven, double *x, double *data, double *invcov,
		    int *Nx, int *Ngiven, int *Dim, int *Rep, 
		    double *Mean, double *Krig, double *sigma2)
{
  int  d, i, k, j, r, krigi,
    err = NOERROR,
    dim = *Dim,
    ngiven = *Ngiven,
    rep = *Rep,
    nx = *Nx,
    len_tgiven = dim * ngiven,
    divachtzig = (nx<79) ? 1 : (nx / 79),
    divachtzigM1 = divachtzig - 1;
  // Ngiven : length given = length M vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  double xx[MAXSIMUDIM], var, origin[MAXSIMUDIM],
    mean = *Mean,
    *lambda = NULL,
    *M = NULL,
    *dist = NULL;
  if ((M = (double*) malloc(sizeof(double) * ngiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((lambda = (double*) malloc(sizeof(double) * ngiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }  
  for (i=0; i<dim; i++) origin[i]=0.0;
  CovIntern(origin, NULL, true, 1, &var);
  for (krigi=i=0; i<nx; i++) {
    if (PL>0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    CovIntern(dist, NULL, true, ngiven, M); 
    for (d=k=0; k<ngiven; k++) {
      lambda[k] = 0.0;
      for (j=0; j<ngiven; j++) lambda[k] += M[j] * invcov[d++];
    }
    sigma2[i] = var;
    for (j=0; j<ngiven; j++) {
      sigma2[i] -= lambda[j] * M[j];
    }
    if (sigma2[i] < 0) {
	// printf("%e", sigma2[i]);
      assert(sigma2[i] > KRIGE_TOLERANCE);
      sigma2[i] = 0.0;
    }
    for (d=r=0; r<rep; r++, krigi++) {
      Krig[krigi] = mean;
      for (j=0; j<ngiven; j++) Krig[krigi] += lambda[j] * data[d++];
    }
  }
  if (PL > 0) PRINTF("\n");

 ErrorHandling:
  if (M!=NULL) free(M);
  if (lambda!=NULL) free(lambda);
  if (dist!=NULL) free(dist);

  if (err != NOERROR) {
      int nKrig;
      nKrig = nx * rep;
      for (i=0; i<nKrig; Krig[i++]=RF_NAN);
  }
}


void ordinaryKriging(double *tgiven, double *x, double *invcov,
		   int *Nx, int *Ngiven,
		   int *Dim, int *Rep, double *Krig)
{
  // kriging variance is not calculated
  int d, i, j, r,  krigi,
      err = NOERROR,
    dim = *Dim,
    ngiven = *Ngiven,
    ngivenP1 = ngiven + 1,
    rep = *Rep,
    nx = *Nx,
    len_tgiven = dim * ngiven,
    divachtzig = (nx<79) ? 1 : (nx / 79),
    divachtzigM1 = divachtzig - 1;
  // Ngiven : length given = length cov vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  double xx[MAXSIMUDIM],
    *M = NULL,
    *dist =NULL;
  if ((M = (double*) malloc(sizeof(double) * ngivenP1))==NULL) { // + 1
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (krigi=i=0; i<nx; i++) {
//    printf("%d\n", i);
    if (PL>0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
//    printf("Xok\n");
    CovIntern(dist, NULL, true, ngiven, M); 
//    printf("ok1\n");
    M[ngiven] = 1.0; // =
    for (d=r=0; r<rep; r++, krigi++) {
      Krig[krigi] = 0.0;
      for (j=0; j<ngivenP1; j++) Krig[krigi] += M[j] * invcov[d++];
    }
  }
  if (PL > 0) PRINTF("\n");

 ErrorHandling:
  if (M!=NULL) free(M);
  if (dist!=NULL) free(dist);
  
  if (err != NOERROR) {
      int nKrig;
      nKrig = nx * rep;
      for (i=0; i<nKrig; Krig[i++]=RF_NAN);
  }
}


void ordinaryKriging2(double *tgiven, double *x, double *data, double *invcov,
		    int *Nx, int *Ngiven, int *Dim, int *Rep, 
		    double *Krig, double *sigma2)
{
  int d, i, j, k, r, krigi,
      err = NOERROR,
    dim = *Dim,
    ngiven = *Ngiven,
    ngivenP1 = ngiven + 1,
    rep = *Rep,
    nx = *Nx,
    len_tgiven = dim * ngiven,
    divachtzig = (nx<79) ? 1 : (nx / 79),
    divachtzigM1 = divachtzig - 1;
  // Ngiven : length given = length cov vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  double xx[MAXSIMUDIM], origin[MAXSIMUDIM], var,
     *M = NULL,
     *lambda = NULL,
     *dist =NULL;
  if ((M = (double*) malloc(sizeof(double) * ngivenP1))==NULL) { // + 1
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((lambda = (double*) malloc(sizeof(double) * ngivenP1))==NULL) { // + 1
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (i=0; i<dim; i++) origin[i]=0.0;
  CovIntern(origin, NULL, true, 1, &var);
  for (krigi=i=0; i<nx; i++) {
    if (PL>0 && (krigi % divachtzig ==divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    CovIntern(dist, NULL, true, ngiven, M); 
    M[ngiven] = 1.0;
    for (d=k=0; k<ngivenP1; k++) { 
      lambda[k] = 0.0;
      for (j=0; j<ngivenP1; j++) lambda[k] += M[j] * invcov[d++];
    }

    sigma2[i] = var;
    for (j=0; j<ngivenP1; j++) sigma2[i] -= lambda[j] * M[j];
    if (sigma2[i] < 0) {
	// printf("%e", sigma2[i]);
      assert(sigma2[i] > KRIGE_TOLERANCE);
      sigma2[i] = 0.0;
    }
    for (d=r=0; r<rep; r++, krigi++) {
      Krig[krigi] = 0.0;
      for (j=0; j<ngiven; j++) Krig[krigi] += lambda[j] * data[d++]; // < ngiven !!
    }
  }
  if (PL > 0) PRINTF("\n");

 ErrorHandling:
  if (M!=NULL) free(M);
  if (lambda!=NULL) free(lambda);
  if (dist!=NULL) free(dist);

  if (err != NOERROR) {
      int nKrig;
      nKrig = nx * rep;
      for (i=0; i<nKrig; Krig[i++]=RF_NAN);
  }
}


void universalKriging(double *tgiven, double *x, double *invcov, int *Nx,
     int *Ngiven, int *Dim, int *Rep, double *Krig, double *tFmatrix, int *NFCT,
     int *iR, int *jR, SEXP f_expr, SEXP rho)
{
  // kriging variance is not calculated
  int  d, i, j, r,  krigi,
      err = NOERROR,
    nfct = *NFCT,
    ngiven = *Ngiven,
    ngivenPnfct = ngiven + nfct,
    dim = *Dim,
    rep = *Rep,
    nx = *Nx,
    len_tgiven = dim * ngiven,
    divachtzig = (nx<79) ? 1 : (nx / 79),
    divachtzigM1 = divachtzig - 1;
  // Ngiven : length given = length cov vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  double xx[MAXSIMUDIM],
    *covf0 = NULL,
    *dist = NULL;
  if ((covf0 = (double*) malloc(sizeof(double) * ngivenPnfct))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  SEXP f_0;
  PROTECT(f_0 = allocVector(REALSXP,1));
  for (krigi=i=0; i<nx; i++) {
    if (PL>0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    //Initializing covf0
    CovIntern(dist, NULL, true, ngiven, covf0); 
    *iR = i + 1;	//in- and output varingivenle
    for(j=0; j < nfct; j++) { //for each function
	*jR = j + 1; //in- and output varingivenle
	f_0 = eval(f_expr, rho);
	covf0[ngiven+j] = REAL(f_0)[0];
    }
    for (d=r=0; r<rep; r++, krigi++) {
      Krig[krigi] = 0.0;
      for (j=0; j<ngivenPnfct; j++) Krig[krigi] += covf0[j] * invcov[d++];
    }
  }
  UNPROTECT(1);
  if (PL > 0) PRINTF("\n");

 ErrorHandling:
  if (covf0!=NULL) free(covf0);
  if (dist!=NULL) free(dist);
  if (err != NOERROR) {
      int nKrig;
      nKrig = nx * rep;
      for (i=0; i<nKrig; Krig[i++]=RF_NAN);
  }
}


void universalKriging2(double *tgiven, double *x, double *data, double *invcov,
	int *Nx, int *Ngiven, int *Dim, int *Rep, double *Krig, double *sigma2,
	double *tFmatrix, int *NFCT, int *iR, int *jR, SEXP f_expr,
	SEXP rho)
{
 // also calculating kriging variance
 int  d, i, k, j, r, krigi,
     err = NOERROR,
   dim = *Dim,
   ngiven = *Ngiven,
   rep = *Rep,
   nx = *Nx,
   len_tgiven = dim*ngiven,
   nfct = *NFCT,
   divachtzig = (nx < 79) ? 1 : (nx / 79),
   divachtzigM1 = divachtzig - 1;
 // Ngiven: number of points given = length M vector to be calculated
 // Rep: repitition data
 // Nx: number of points in x, length(x) = nx * dim
 // NrowPowers: Number of functions to form the trend
 double  xx[MAXSIMUDIM], var, origin[MAXSIMUDIM],
   *lambda = NULL,
   *M = NULL,
   *dist = NULL,
   *fvector = NULL,
   *X1 = NULL,
   *lambdak = NULL,
   *Qmatrix = NULL,
   *Rvector = NULL,
   *mu = NULL;
 if ((M = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((lambda = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((fvector = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((X1 = (double*) malloc(sizeof(double) * ngiven * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((lambdak = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Qmatrix = (double*) malloc(sizeof(double) * nfct * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Rvector = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((mu = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }

 //Calculating X1 and Q, inverting Q
  matmult(invcov, tFmatrix, X1, ngiven, ngiven, nfct);
  matmulttransposed(tFmatrix, X1, Qmatrix, ngiven, nfct, nfct);
  InvChol(Qmatrix, nfct);

 //Calculating the variance
 for (j=0; j < dim; j++) origin[j]=0.0;
 CovIntern(origin, NULL, true, 1, &var);

 SEXP f_0;
 PROTECT(f_0 = allocVector(REALSXP,1));
 //Kriging for each kriging point
 for (krigi=i=0; i<nx; i++) {
  if (PL > 0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
  //Determining x
  for(d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
  //Calculating cov vector
  for(j=d=0; j < len_tgiven; j++, d= (d+1) % dim) dist[j] = tgiven[j] - xx[d];
  CovIntern(dist, NULL, true, ngiven, M); 
  //Calculating fvector
  *iR = i + 1;   //in- and output varingivenle
  for (j=0; j<nfct; j++)  {
	*jR = j + 1; //in- and output varingivenle
	f_0 = eval(f_expr, rho);
	fvector[j] = REAL(f_0)[0];
  }
  //Solving the kriging equations
  //Calculating lambak
  matmult(invcov, M, lambdak, ngiven, ngiven, 1);
  //Initializing R
  matmulttransposed(tFmatrix, lambdak, Rvector, ngiven, nfct, 1);
  for(k=0; k<nfct; k++) Rvector[k] = Rvector[k] - fvector[k];
  //Calculating mu
  matmult(Qmatrix, Rvector, mu, nfct, nfct, 1);
   //Calculating lambda
  matmult(X1, mu, lambda, ngiven, nfct, 1);
  for(k=0; k<ngiven; k++) lambda[k] = lambdak[k] - lambda[k];
  sigma2[i] = var;
  for (j=0; j<ngiven; j++) sigma2[i] -= lambda[j] * M[j];
  for (j=0; j<nfct; j++) sigma2[i] -= mu[j] * fvector[j];
  if (sigma2[i] < 0) {
      //printf("%e\n", sigma2[i]);
      assert(sigma2[i] > KRIGE_TOLERANCE);
      sigma2[i] = 0.0;
  }
  for (d=r=0; r<rep; r++, krigi++) {
    Krig[krigi] = 0.0;
    for (j=0; j<ngiven; j++) Krig[krigi] += lambda[j] * data[d++];
  }
 }
 UNPROTECT(1);
 if (PL > 0) PRINTF("\n");

 ErrorHandling:
   if (M!=NULL) free(M);
   if (lambda!=NULL) free(lambda);
   if (dist!=NULL) free(dist);
   if (fvector!=NULL) free(fvector);
   if (X1 !=NULL) free(X1);
   if (lambdak !=NULL) free(lambdak);
   if (Qmatrix!=NULL) free(Qmatrix);
   if (Rvector!=NULL) free(Rvector);
   if (mu!=NULL) free(mu);
   if (err != NOERROR) {
       int nKrig;
       nKrig = nx * rep;
       for(i=0; i<nKrig; Krig[i++]=RF_NAN);
   }
}


void universalKriging3(double *tgiven, double *x, double *invcov, int *Nx,
     int *Ngiven, int *Dim, int *Rep, double *Krig, double *tPowers, int *NFCT)
{
  // kriging variance is not calculated
  int d, i, j, r, krigi, 
      err = NOERROR,
    dim = *Dim,
    ngiven = *Ngiven,
    nfct = *NFCT,
    rep = *Rep,
    nx = *Nx,
    len_tgiven = dim * ngiven,
    divachtzig = (nx<79) ? 1 : (nx / 79),
    divachtzigM1 = divachtzig - 1;
  double xx[MAXSIMUDIM],
     *covf0 = NULL,
     *dist = NULL;
  // Ngiven : length given = length cov vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  if ((covf0 = (double*) malloc(sizeof(double) * (ngiven+nfct)))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (krigi=i=0; i<nx; i++) {
    if (PL>0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    //Initializing covf0
    CovIntern(dist, NULL, true, ngiven, covf0); 
    for(j=0; j < nfct; j++) { //for each function
	covf0[ngiven+j] = 1;
	for(d=0; d < dim; d++) //for each component
	  covf0[ngiven+j] *= pow(xx[d], tPowers[j*dim + d]);
    }
    for (d=r=0; r<rep; r++, krigi++) {
      Krig[krigi] = 0.0;
      for (j=0; j<ngiven+nfct; j++) Krig[krigi] += covf0[j] * invcov[d++];
    }
  }
  if (PL > 0) PRINTF("\n");

 ErrorHandling:
  if (covf0!=NULL) free(covf0);
  if (dist!=NULL) free(dist);
  if (err != NOERROR) {
      int nKrig;
      nKrig = nx * rep;
      for (i=0; i<nKrig; Krig[i++]=RF_NAN);
  }
}


void universalKriging4(double *tgiven, double *x, double *data, double *invcov,
	int *Nx, int *Ngiven, int *Dim, int *Rep, double *Krig, double *sigma2,
	double *tPowers, int *NFCT)
{
 // also calculating kriging variance
 int  d, i, k, j, r, krigi,
     err = NOERROR,
   dim = *Dim,
   ngiven = *Ngiven,
   rep = *Rep,
   nx = *Nx,
   len_tgiven = dim*ngiven,
   nfct = *NFCT,
   divachtzig = (nx < 79) ? 1 : (nx / 79),
   divachtzigM1 = divachtzig - 1;
 // Ngiven: number of points given = length cov vector to be calculated
 // Rep: repitition data
 // Nx: number of points in x, length(x) = nx * dim
 // NrowPowers: Number of functions to form the trend
 double xx[MAXSIMUDIM], var, origin[MAXSIMUDIM],
   *lambda = NULL,
   *M = NULL,
   *dist = NULL,
   *Fmatrix = NULL,
   *fvector = NULL,
   *X1 = NULL,
   *lambdak = NULL,
   *Qmatrix = NULL,
   *Rvector = NULL,
   *mu = NULL;
 if ((M = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((lambda = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Fmatrix = (double*) malloc(sizeof(double) * ngiven * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((fvector = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((X1 = (double*) malloc(sizeof(double) * ngiven * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((lambdak = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Qmatrix = (double*) malloc(sizeof(double) * nfct * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Rvector = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((mu = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }

 //Calculate Fmatrix
 double *PFmatrix;
 PFmatrix = Fmatrix;
 for (i=0; i<ngiven; i++) //for each given point
   for(j=0; j < nfct; j++, PFmatrix++) { //for each function
	*PFmatrix = 1;
	for(d=0; d < dim; d++) //for each component
	   *PFmatrix *= pow(tgiven[dim*i + d], tPowers[j*dim + d]);
   }
 //Calculating X1 and Q, inverting Q
 matmult(invcov, Fmatrix, X1, ngiven, ngiven, nfct);
 matmulttransposed(Fmatrix, X1, Qmatrix, ngiven, nfct, nfct);
 //bei variogrammen ist Q i.A. nicht strikt positiv definit
 InvChol(Qmatrix, nfct);

 //Calculating the variance
 for (j=0; j < dim; j++) origin[j]=0.0;
 CovIntern(origin, NULL, true, 1, &var);
 //Kriging for each kriging point
 for (krigi=i=0; i<nx; i++) {
   if (PL > 0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
   //Determining x
   for(d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
   //Calculating cov vector
   for(j=d=0; j < len_tgiven; j++, d= (d+1) % dim) dist[j] = tgiven[j] - xx[d];
   CovIntern(dist, NULL, true, ngiven, M); 
   //Calculating fvector
   for(j=0; j < nfct; j++) { //for each function
	fvector[j] = 1;
	for(d=0; d < dim; d++) //for each component
		fvector[j] *= pow(xx[d], tPowers[j*dim + d]);
   }
  //Solving the kriging equations
  //Calculating lambak
  matmult(invcov, M, lambdak, ngiven, ngiven, 1);
  //Initializing R
  matmulttransposed(Fmatrix, lambdak, Rvector, ngiven, nfct, 1);
  for(k=0; k<nfct; k++) Rvector[k] = Rvector[k] - fvector[k];
  //Calculating mu
  matmult(Qmatrix, Rvector, mu, nfct, nfct, 1);
   //Calculating lambda
  matmult(X1, mu, lambda, ngiven, nfct, 1);
  for(k=0; k<ngiven; k++) lambda[k] = lambdak[k] - lambda[k];
  sigma2[i] = var;
  for (j=0; j<ngiven; j++) sigma2[i] -= lambda[j] * M[j];
  for (j=0; j<nfct; j++) sigma2[i] -= mu[j] * fvector[j];
  if (sigma2[i] < 0) {
      //printf("%e\n", sigma2[i]);
      assert(sigma2[i] > KRIGE_TOLERANCE);
      sigma2[i] = 0.0;
  }
  for (d=r=0; r<rep; r++, krigi++) {
    Krig[krigi] = 0.0;
    for (j=0; j<ngiven; j++) Krig[krigi] += lambda[j] * data[d++];
  }
 }
 if (PL > 0) PRINTF("\n");

 ErrorHandling:
   if (M!=NULL) free(M);
   if (lambda!=NULL) free(lambda);
   if (dist!=NULL) free(dist);
   if (Fmatrix!=NULL) free(Fmatrix);
   if (fvector!=NULL) free(fvector);
   if (X1 !=NULL) free(X1);
   if (lambdak !=NULL) free(lambdak);
   if (Qmatrix!=NULL) free(Qmatrix);
   if (Rvector!=NULL) free(Rvector);
   if (mu!=NULL) free(mu);
   if (err != NOERROR) {
       int nKrig;
       nKrig = nx * rep;
      for(i=0; i<nKrig; Krig[i++]=RF_NAN);
   }
}


void intrinsicKriging(double *tgiven, double *x, double *invcov, int *Nx,
     int *Ngiven, int *Dim, int *Rep, double *Krig, double *tPowers, int *NFCT)
{
  // kriging variance is not calculated
  int  d, i, j, r, krigi, 
      err = NOERROR,
    dim = *Dim,
    ngiven = *Ngiven,
    nfct = *NFCT,
    rep = *Rep,
    nx = *Nx,
    len_tgiven = dim * ngiven,
    divachtzig = (nx<79) ? 1 : (nx / 79),
    divachtzigM1 = divachtzig - 1;
  double xx[MAXSIMUDIM],
     *covf0 = NULL,
     *dist = NULL;
  // Ngiven : length given = length cov vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  if ((covf0 = (double*) malloc(sizeof(double) * (ngiven+nfct)))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (krigi=i=0; i<nx; i++) {
    if (PL>0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    //Initializing covf0
    CalculateVariogram(dist, ngiven, STORED_MODEL[MODEL_INTERN], covf0);
    for(j=0; j<ngiven; j++) covf0[j] = (-1)*covf0[j];
    for(j=0; j < nfct; j++) { //for each function
	covf0[ngiven+j] = 1;
	for(d=0; d < dim; d++) //for each component
	  covf0[ngiven+j] *= pow(xx[d], tPowers[j*dim + d]);
    }
    for (d=r=0; r<rep; r++, krigi++) {
      Krig[krigi] = 0.0;
      for (j=0; j<ngiven+nfct; j++) Krig[krigi] += covf0[j] * invcov[d++];
    }
  }
  if (PL > 0) PRINTF("\n");

 ErrorHandling:
  if (covf0!=NULL) free(covf0);
  if (dist!=NULL) free(dist);
  if (err != NOERROR) {
      int nKrig;
      nKrig = nx * rep;
      for (i=0; i<nKrig; Krig[i++]=RF_NAN);
  }
}

void intrinsicKriging2(double *tgiven, double *x, double *data, double *invcov,
	int *Nx, int *Ngiven, int *Dim, int *Rep, double *Krig, double *sigma2,
	double *tPowers, int *NFCT, SEXP inv_expr, SEXP rho)
{
 // also calculating kriging variance
 int d, i, k, j, r, krigi,
     err = NOERROR,
   one=1,
   dim = *Dim,
   ngiven = *Ngiven,
   rep = *Rep,
   nx = *Nx,
   len_tgiven = dim*ngiven,
   nfct = *NFCT,
   nfct2 = nfct*nfct,
   divachtzig = (nx < 79) ? 1 : (nx / 79),
   divachtzigM1 = divachtzig - 1;
 // Ngiven: number of points given = length cov vector to be calculated
 // Rep: repitition data
 // Nx: number of points in x, length(x) = nx * dim
 // NrowPowers: Number of functions to form the trend
 double xx[MAXSIMUDIM], var, origin[MAXSIMUDIM],
   *lambda = NULL,
   *M = NULL,
   *dist = NULL,
   *Fmatrix = NULL,
   *fvector = NULL,
   *X1 = NULL,
   *lambdak = NULL,
   *Qmatrix = NULL,
   *Rvector = NULL,
   *mu = NULL;
 if ((M = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((lambda = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Fmatrix = (double*) malloc(sizeof(double) * ngiven * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((fvector = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((X1 = (double*) malloc(sizeof(double) * ngiven * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((lambdak = (double*) malloc(sizeof(double) * ngiven))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Qmatrix = (double*) malloc(sizeof(double) * nfct * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((Rvector = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }
 if ((mu = (double*) malloc(sizeof(double) * nfct))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
 }

 //Calculate Fmatrix
 double *PFmatrix;
 PFmatrix = Fmatrix;
 for (i=0; i<ngiven; i++) //for each given point
   for(j=0; j < nfct; j++, PFmatrix++) { //for each function
	*PFmatrix = 1;
	for(d=0; d < dim; d++) //for each component
	   *PFmatrix *= pow(tgiven[dim*i + d], tPowers[j*dim + d]);
   }
 //Calculating X1 and Q, inverting Q
 matmult(invcov, Fmatrix, X1, ngiven, ngiven, nfct);
 matmulttransposed(Fmatrix, X1, Qmatrix, ngiven, nfct, nfct);
 //bei variogrammen ist Q i.A. nicht strikt positiv definit
 SEXP invMatrix;
 PROTECT(invMatrix = allocMatrix(REALSXP, nfct, nfct));
 for(j=0; j < nfct2; j++)  REAL(invMatrix)[j] = Qmatrix[j];
 defineVar(install("Q"), invMatrix, rho);
 invMatrix = eval(inv_expr, rho);
 for(j=0; j < nfct2; j++)  Qmatrix[j] = REAL(invMatrix)[j];
 UNPROTECT(1);

 //Calculating the variance
 for (j=0; j < dim; j++) origin[j]=0.0;
    CalculateVariogram(origin, one, STORED_MODEL[MODEL_INTERN], &var);
    var = (-1)*var;
 //Kriging for each kriging point
 for (krigi=i=0; i<nx; i++) {
   if (PL > 0 && (krigi % divachtzig == divachtzigM1)) PRINTF(".");
   //Determining x
   for(d=0, j=i; d<dim; d++, j+=nx) xx[d]=x[j];
   //Calculating cov vector
   for(j=d=0; j < len_tgiven; j++, d= (d+1) % dim) dist[j] = tgiven[j] - xx[d];
   CalculateVariogram(dist, ngiven, STORED_MODEL[MODEL_INTERN], M);
   for(j=0; j<ngiven; j++) M[j] = (-1)*M[j];
   //Calculating fvector
   for(j=0; j < nfct; j++) { //for each function
	fvector[j] = 1;
	for(d=0; d < dim; d++) //for each component
		fvector[j] *= pow(xx[d], tPowers[j*dim + d]);
   }
  //Solving the kriging equations
  //Calculating lambak
  matmult(invcov, M, lambdak, ngiven, ngiven, 1);
  //Initializing R
  matmulttransposed(Fmatrix, lambdak, Rvector, ngiven, nfct, 1);
  for(k=0; k<nfct; k++) Rvector[k] = Rvector[k] - fvector[k];
  //Calculating mu
  matmult(Qmatrix, Rvector, mu, nfct, nfct, 1);
   //Calculating lambda
  matmult(X1, mu, lambda, ngiven, nfct, 1);
  for(k=0; k<ngiven; k++) lambda[k] = lambdak[k] - lambda[k];
  sigma2[i] = var;
  for (j=0; j<ngiven; j++) sigma2[i] -= lambda[j] * M[j];
  for (j=0; j<nfct; j++) sigma2[i] -= mu[j] * fvector[j];
  if (sigma2[i] < 0) {
      //printf("%e\n", sigma2[i]);
      assert(sigma2[i] > KRIGE_TOLERANCE);
      sigma2[i] = 0.0;
  }
  for (d=r=0; r<rep; r++, krigi++) {
    Krig[krigi] = 0.0;
    for (j=0; j<ngiven; j++) Krig[krigi] += lambda[j] * data[d++];
  }
 }
 if (PL > 0) PRINTF("\n");


 ErrorHandling:
   if (M!=NULL) free(M);
   if (lambda!=NULL) free(lambda);
   if (dist!=NULL) free(dist);
   if (Fmatrix!=NULL) free(Fmatrix);
   if (fvector!=NULL) free(fvector);
   if (X1 !=NULL) free(X1);
   if (lambdak !=NULL) free(lambdak);
   if (Qmatrix!=NULL) free(Qmatrix);
   if (Rvector!=NULL) free(Rvector);
   if (mu!=NULL) free(mu);
   if (err != NOERROR) {
       int nKrig;
       nKrig = nx * rep;
       for(i=0; i<nKrig; Krig[i++]=RF_NAN);
   }
}


void matmult(double *A, double *B, double *C, int l, int m, int n) {
// multiplying an lxm- and an mxn-matrix, saving krigult in C
   int i, j, k;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
	C[i*n + j] = 0;
	for(k=0; k<m; k++)
	   C[i*n +j] = C[i*n+j] + A[i*m+k]*B[k*n+j];
     }
}

void matmulttransposed(double *A, double *B, double *C, int m, int l, int n) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
   int i, j, k;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
	C[i*n + j] = 0;
	for(k=0; k<m; k++)
	   C[i*n +j] += A[k*l+i]*B[k*n+j];
     }
}

void InvChol(double *C, int dim) {
  int i, info, ii, endfor,
    job = 01,
    dimP1 = dim + 1,
    dimsq = dim * dim;
  long ve, ho;
  double Det = 1.0;
  F77_CALL(dpofa)(C, &dim, &dim, &info); // C is now cholesky
  if (info != 0) ERR("Inversion failed, bad functions\n");
  for (i=0; i<dimsq; i+=dimP1) Det *= C[i];
  Det = Det * Det;
  F77_CALL(dpodi)(C, &dim, &dim, &Det, &job); // C is now Cinv
  for (ii=dim, i=0; ii>0; i+=dimP1, ii--) {  // filling lower half
      endfor = i + ii;
      for (ve = i + 1, ho = i + dim; ve < endfor;
	   ve++, ho += dim) {
	C[ve] = C[ho];
      }
  }
}
