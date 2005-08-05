/*
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for unconditional simulation of stationary and isotropic random fields

 Copyright (C) 2001 -- 2005 Martin Schlather, 

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
#include <assert.h>
#include <string.h>
#include "RFsimu.h"

#define KRIGE_TOLERANCE -1e-10

void simpleKriging(double *tgiven, double *x, double *invcov,
		   int *Len_x, int *NN,
		   int *Dim, int *Rep, double *Mean, double *Res)
{
  // kriging variance is not calculated
  int error, len_tgiven, dim, nn, d, i, j, r, rep, resi, len_x, divachtzig, 
    divachtzigM1;
  double *cov, *dist, xx[MAXDIM], mean;
  // NN : length given = length cov vector to be calculated
  // Rep : repetition invcov
  // Len_x : number of points in x, length(x) = len_x * dim
  dim = *Dim;
  nn = *NN;
  mean = *Mean;
  rep = *Rep;
  len_x = *Len_x;
  len_tgiven = dim * nn;
  divachtzigM1 = (divachtzig = (len_x<79) ? 1 : (len_x / 79)) - 1;
  cov = dist =NULL;
  if ((cov = (double*) malloc(sizeof(double) * nn))==NULL) {
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (resi=i=0; i<len_x; i++) {
    // printf("%d %d %d \n", i, len_x, divachtzig);
    if (GENERAL_PRINTLEVEL>0 && (resi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=len_x) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    UncheckedCovFct(dist, NN, cov);
    for (d=r=0; r<rep; r++, resi++) {
      Res[resi] = mean;
      for (j=0; j<nn; j++) Res[resi] += cov[j] * invcov[d++];
    }
  }
  if (GENERAL_PRINTLEVEL > 0) PRINTF("\n");
  free(cov);
  free(dist);
  return;

 ErrorHandling:
  int nRes;
  if (cov!=NULL) free(cov);
  if (dist!=NULL) free(cov);
  nRes = len_x * rep;
  for (i=0; i<nRes; Res[i++]=RF_NAN);
}

void simpleKriging2(double *tgiven, double *x, double *data, double *invcov,
		    int *Len_x, int *NN, int *Dim, int *Rep, 
		    double *Mean, double *Res, double *sigma2)
{
  int error, len_tgiven, dim, nn, d, i, k, j, r, rep, resi, len_x, divachtzig, 
    divachtzigM1, one=1;
  double *cov, *dist, xx[MAXDIM], mean, *lambda, var, origin[MAXDIM];
  // NN : length given = length cov vector to be calculated
  // Rep : repetition data
  // Len_x : number of points in x, length(x) = len_x * dim
  dim = *Dim;
  nn = *NN;
  mean = *Mean;
  rep = *Rep;
  len_x = *Len_x;
  len_tgiven = dim * nn;
  divachtzigM1 = (divachtzig = (len_x<79) ? 1 : (len_x / 79)) - 1;
  lambda = cov = dist = NULL;
  if ((cov = (double*) malloc(sizeof(double) * nn))==NULL) {
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((lambda = (double*) malloc(sizeof(double) * nn))==NULL) {
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }  
  for (i=0; i<dim; i++) origin[i]=0.0;
  UncheckedCovFct(origin, &one, &var);
  for (resi=i=0; i<len_x; i++) {
    if (GENERAL_PRINTLEVEL>0 && (resi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=len_x) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    UncheckedCovFct(dist, NN, cov);
    for (d=k=0; k<nn; k++) {
      lambda[k] = 0.0;
      for (j=0; j<nn; j++) lambda[k] += cov[j] * invcov[d++];
    }
    sigma2[i] = var;
    for (j=0; j<nn; j++) {
      sigma2[i] -= lambda[j] * cov[j];
    }
    if (sigma2[i] < 0) {
	// printf("%e", sigma2[i]);
      assert(sigma2[i] > KRIGE_TOLERANCE);
      sigma2[i] = 0.0;
    }
    for (d=r=0; r<rep; r++, resi++) {
      Res[resi] = mean;
      for (j=0; j<nn; j++) Res[resi] += lambda[j] * data[d++];
    }
  }
  if (GENERAL_PRINTLEVEL > 0) PRINTF("\n");
  free(cov);
  free(lambda);
  free(dist);
  return;

 ErrorHandling:
  int nRes;
  if (cov!=NULL) free(cov);
  if (lambda!=NULL) free(lambda);
  if (dist!=NULL) free(dist);
  nRes = len_x * rep;
  for (i=0; i<nRes; Res[i++]=RF_NAN);
}

void ordinaryKriging(double *tgiven, double *x, double *invcov,
		   int *Len_x, int *NN,
		   int *Dim, int *Rep, double *Res)
{
  // kriging variance is not calculated
  int error, len_tgiven, dim, nn, nnP1, d, i, j, r, rep, resi, len_x, 
    divachtzig, divachtzigM1;
  double *cov, *dist, xx[MAXDIM];
  // NN : length given = length cov vector to be calculated
  // Rep : repetition data
  // Len_x : number of points in x, length(x) = len_x * dim
  dim = *Dim;
  nn = *NN;
  nnP1 = nn + 1;
  rep = *Rep;
  len_x = *Len_x;
  len_tgiven = dim * nn;
  cov = dist =NULL;
  divachtzigM1 = (divachtzig = (len_x<79) ? 1 : (len_x / 79)) - 1;
  if ((cov = (double*) malloc(sizeof(double) * nnP1))==NULL) { // + 1
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (resi=i=0; i<len_x; i++) {
    if (GENERAL_PRINTLEVEL>0 && (resi % divachtzig == divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=len_x) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
   UncheckedCovFct(dist, NN, cov);
   cov[nn] = 1.0; // =
   for (d=r=0; r<rep; r++, resi++) {
      Res[resi] = 0.0;
      for (j=0; j<nnP1; j++) Res[resi] += cov[j] * invcov[d++];
    }
  }
  if (GENERAL_PRINTLEVEL > 0) PRINTF("\n");
  free(cov);
  free(dist);
  return;

 ErrorHandling:
  int nRes;
  if (cov!=NULL) free(cov);
  if (dist!=NULL) free(dist);
  nRes = len_x * rep;
  for (i=0; i<nRes; Res[i++]=RF_NAN);
}

void ordinaryKriging2(double *tgiven, double *x, double *data, double *invcov,
		    int *Len_x, int *NN, int *Dim, int *Rep, 
		    double *Res, double *sigma2)
{
  int error, len_tgiven, dim, nn, nnP1, d, i, j, k, r, rep, resi, len_x,
    divachtzig, divachtzigM1, one=1;
  double *cov, *dist, xx[MAXDIM], *lambda, origin[MAXDIM], var;
  // NN : length given = length cov vector to be calculated
  // Rep : repetition data
  // Len_x : number of points in x, length(x) = len_x * dim
  dim = *Dim;
  nn = *NN;
  nnP1 = nn + 1;
  rep = *Rep;
  len_x = *Len_x;
  len_tgiven = dim * nn;
  lambda = dist =NULL;
  divachtzigM1 = (divachtzig = (len_x<79) ? 1 : (len_x / 79)) - 1;
  if ((cov = (double*) malloc(sizeof(double) * nnP1))==NULL) { // + 1
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((lambda = (double*) malloc(sizeof(double) * nnP1))==NULL) { // + 1
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  if ((dist = (double*) malloc(sizeof(double) * len_tgiven))==NULL) {
    error = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (i=0; i<dim; i++) origin[i]=0.0;
  UncheckedCovFct(origin, &one, &var);
  for (resi=i=0; i<len_x; i++) {
    if (GENERAL_PRINTLEVEL>0 && (resi % divachtzig ==divachtzigM1)) PRINTF(".");
    for (d=0, j=i; d<dim; d++, j+=len_x) xx[d]=x[j];
    for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    UncheckedCovFct(dist, NN, cov);
    cov[nn] = 1.0;
    for (d=k=0; k<nnP1; k++) { 
      lambda[k] = 0.0;
      for (j=0; j<nnP1; j++) lambda[k] += cov[j] * invcov[d++];
    }

    sigma2[i] = var;
    for (j=0; j<nnP1; j++) sigma2[i] -= lambda[j] * cov[j];
    if (sigma2[i] < 0) {
	// printf("%e", sigma2[i]);
      assert(sigma2[i] > KRIGE_TOLERANCE);
      sigma2[i] = 0.0;
    }
    for (d=r=0; r<rep; r++, resi++) {
      Res[resi] = 0.0;
      for (j=0; j<nn; j++) Res[resi] += lambda[j] * data[d++]; // < nn !!
    }
  }
  if (GENERAL_PRINTLEVEL > 0) PRINTF("\n");
  free(cov);
  free(lambda);
  free(dist);
  return;

 ErrorHandling:
  int nRes;
  if (cov!=NULL) free(cov);
  if (lambda!=NULL) free(lambda);
  if (dist!=NULL) free(dist);
  nRes = len_x * rep;
  for (i=0; i<nRes; Res[i++]=RF_NAN);
}

