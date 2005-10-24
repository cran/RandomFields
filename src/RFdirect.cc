/*
 Authors
 Martin Schlather, schlath@hsu-hh.de

 Simulation of a random field by spectral turning bands

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

#include <math.h>
#include <stdio.h>  
#include <stdlib.h>
#include <assert.h>
#include "RFsimu.h"
#include <R_ext/Lapack.h>
//#include <R_ext/Linpack.h>


InversionMethod DIRECTGAUSS_INVERSIONMETHOD=Cholesky;
int DIRECTGAUSS_BESTVARIABLES=800; // see RFsimu.h NOTE: if changed, 
// then also the example in RFparameters must be reconsidered
int DIRECTGAUSS_MAXVARIABLES=4096; // see RFsimu.h


void direct_destruct(void ** S)
{
  if (*S!=NULL) {
    direct_storage *x;
    x = *((direct_storage**)S);    
    if (x->U!=NULL) free(x->U); 
    if (x->G!=NULL) free(x->G);
    free(*S);
    *S = NULL;
  }
}

void SetParamDirectGauss(int *action,int *method, int *bestvariables,
			 int *maxvariables)
{
  if (*action) {
    if ((*method<0) || (*method>=(int) NoFurtherInversionMethod)){
      PRINTF("inversion method out of range; ignored\n");
    }
    else {DIRECTGAUSS_INVERSIONMETHOD= (InversionMethod) *method;}
   DIRECTGAUSS_BESTVARIABLES=*bestvariables;
   DIRECTGAUSS_MAXVARIABLES=*maxvariables;
  } else {
    *method = (int) DIRECTGAUSS_INVERSIONMETHOD;
    *bestvariables = DIRECTGAUSS_BESTVARIABLES;
    *maxvariables = DIRECTGAUSS_MAXVARIABLES;
  }
}  



int init_directGauss(key_type *key, int m) 
{
  long Xerror=NOERROR, totpnts, i;
  double *xx[MAXDIM], *G, *COV, *U, *VT, *work, *D;
  int d, *iwork,
      dim=key->timespacedim; bool freexx=false;
  direct_storage* S;
  methodvalue_type *meth; 

  G=COV=U=VT=work=D=NULL;
  iwork = NULL;
  S=NULL;
  for (d=0; d<dim; d++) xx[d]=NULL;
  meth = &(key->meth[m]);
  SET_DESTRUCT(direct_destruct, m);
  totpnts = key->totalpoints;
  if (totpnts>DIRECTGAUSS_MAXVARIABLES) {
     if (GENERAL_PRINTLEVEL>3)
	 PRINTF("cur. points=%d, max points=%d", 
	    (int) totpnts, (int) DIRECTGAUSS_MAXVARIABLES);
   Xerror=ERRORMETHODNOTALLOWED; goto ErrorHandling;
  }

  Xerror=ERRORMEMORYALLOCATION;  
  if ((COV =(double *) malloc(sizeof(double) * totpnts * totpnts))==NULL)
    goto ErrorHandling;
  if ((U =(double *) malloc(sizeof(double) * totpnts * totpnts))==NULL)
      goto ErrorHandling;
  //for SVD/Chol intermediate results AND  memory space for do_directGauss:
  if ((G = (double *)  malloc(sizeof(double) * (totpnts + 1)))==NULL)
      goto ErrorHandling;
  if ((meth->S = malloc(sizeof(direct_storage))) == NULL)
    goto ErrorHandling;
  S = (direct_storage*) meth->S;
  S->U = S->G = NULL;

  if ((Xerror = FirstCheck_Cov(key, m, true)) !=NOERROR) goto ErrorHandling;

  /* ************************* */
  /* create matrix of explicit */
  /*       x-coordinates       */
  /* ************************* */
  if (key->grid) {
    int index[MAXDIM];
    double step[MAXDIM], p[MAXDIM];
     
    // generate all the grid coordinates exlicitely!
    freexx = true;
    for (d=0; d<dim; d++) xx[d] = NULL;
    for (d=0; d<dim; d++) {
      if ((xx[d]=(double*) malloc(sizeof(double) * totpnts)) == NULL) {
	  Xerror=ERRORMEMORYALLOCATION;  goto ErrorHandling;
      }
      xx[d][0] = p[d] = 0.0;
      index[d]=0;
      step[d]=key->x[d][XSTEP];
    }
    // intialisation of index and p, and defining x[d][i=0]
    // index preciser than the real coordinate values (which are created in p[])
    // index only used to compare with length*
    for (i=1; i<totpnts; i++) {
      d = 0;  //d=dim-1; if x should run the slowest;
      (index[d])++; 
      p[d] += step[d];
      while (index[d] >= key->length[d]) { 
	index[d]=0; 
	p[d]=0.0; 
	d++;   // d--;
	assert(d<dim); // assert(d>=0);
	(index[d])++; 
	p[d] += step[d];
      }
      for (d=0; d<dim; d++) 
	xx[d][i] = p[d];
    }
  } else { // no grid
    if (key->Time) {
      int endfor, t, spatialdim, j;
      double time, step;
      freexx = true;
      spatialdim = key->spatialdim;
      for (d=0; d<dim; d++) { 
	if ((xx[d]=(double*) malloc(sizeof(double)* totpnts))==NULL)
	    goto ErrorHandling;
      }
      endfor = key->length[spatialdim]; 
      step = key->T[XSTEP];
      for (t=0, i=0, time=key->T[XSTART]; t<endfor; t++, time += step) 
	for (j=0; j < key->length[0]; j++, i++) {
	  for (d=0; d<spatialdim; d++) xx[d][i]=key->x[d][j];
	  xx[spatialdim][i] = time;
	}
    } else { 
      freexx = false;
      for (i=0; i<dim; i++) xx[i] = key->x[i];
    }
  }


  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  long j, k, segment;
  int actcov, row, err;
  InversionMethod method;

  actcov = meth->actcov;
  k = 0;
  for (i=0; i<totpnts; i++) {
    double y[MAXDIM];
    k += i;
    for (segment=i* totpnts+i, j=i; j<totpnts; j++) {
      for (d=0; d<dim; d++) 
	y[d] = xx[d][i] - xx[d][j];
      COV[k++] = COV[segment] = 
	  key->covFct(y, dim, key->cov, meth->covlist, actcov, key->anisotropy);
      segment += totpnts;
    }
  }

  /* ********************* */
  /* matrix inversion part */
  /* ********************* */
  method = DIRECTGAUSS_INVERSIONMETHOD;
  switch (method) {
      case Cholesky :

	if (GENERAL_PRINTLEVEL>=3) PRINTF("method=Cholesky\n");
	row=totpnts;
	// dchdc destroys the input matrix; upper half of U contains result!
	memcpy(U, COV, sizeof(double) * row * row);
	// dpotrf	F77_CALL(dpotrf)("Upper", &m, REAL(ans), &m, &i);
	F77_CALL(dpotrf)("Upper", &row, U, &row, &err);
	// F77_NAME(dchdc)(U, &row, &row, G, NULL, &choljob, &err);

	if (err!=NOERROR) {
	  if (GENERAL_PRINTLEVEL>2)
	      PRINTF("Error code F77_CALL(dpotrf) = %d\n", err);
	  Xerror=ERRORDECOMPOSITION;
	} else break;

	// try next method : 
	// most common error: singular matrix 
	
      case SVD :
        int jobint; 
	jobint = 11;
	method = SVD; // necessary if the value of method has been Cholesky.
	//               originally
	if (GENERAL_PRINTLEVEL>=3) PRINTF("method=SVD\n");
	Xerror=ERRORMEMORYALLOCATION;
	if ((VT =(double *) malloc(sizeof(double) * totpnts * totpnts))==NULL)
	    goto ErrorHandling;
	if ((D =(double *) malloc(sizeof(double) * totpnts))==NULL)
	    goto ErrorHandling;
	if ((iwork = (int *) malloc(sizeof(int) * 8 * totpnts))==NULL)
	    goto ErrorHandling;
	row=totpnts;
        // dsvdc destroys the input matrix !!!!!!!!!!!!!!!!!!!!


	// DGESDD (or DGESVD)
      	// dgesdd destroys the input matrix COV;
	// F77_NAME(dsvdc)(COV, &row, &row, &row, D, e, U, &row, V, &row, G,
	//		&jobint, &err);
	double optim_lwork;
	int lwork;
	lwork = -1;
	F77_CALL(dgesdd)("A",  &row,  &row, COV, &row, D, U, &row, VT, &row, 
			 &optim_lwork, &lwork, iwork, &err);
	if ((Xerror=err) != NOERROR) {
	  Xerror=ERRORDECOMPOSITION;
	  goto ErrorHandling;
	}
	lwork = (int) optim_lwork;
	if ((work = (double *) malloc(sizeof(double) * lwork))==NULL)
	    goto ErrorHandling;
	F77_CALL(dgesdd)("A",  &row,  &row, COV, &row, D, U, &row, VT, &row, 
			 work, &lwork, iwork, &err);
	
	if (Xerror==NOERROR && RF_ISNA(D[0]))
	    Xerror=9999;
	if (Xerror!=NOERROR) {
	  if (GENERAL_PRINTLEVEL>2) 
	    PRINTF("Error code F77_CALL(dgesdd) = %d\n", Xerror); 
	  Xerror=ERRORDECOMPOSITION;
	  goto ErrorHandling;
	}
			
	/* calculate SQRT of covariance matrix */
	for (k=0,j=0;j<totpnts;j++) {
	  register double dummy;
	  dummy = sqrt(D[j]);
	  for (i=0;i<totpnts;i++) {
	    U[k++] *= dummy;
	  }
	}
	break;

	default : assert(false);
  } // switch

  if (Xerror!=NOERROR) { goto ErrorHandling; }

  
  S->U=U;
  S->method=method;
  S->G=G;
  free(COV);
  if (D!=NULL) free(D);
  if (work!=NULL) free(work);
  if (iwork!=NULL) free(iwork);
  if (VT!=NULL) free(VT);  
  if (freexx) for (i=0; i<dim; i++) {free(xx[i]);}
  return NOERROR;

 ErrorHandling: 
  if (freexx)
      for (i=0; i<dim; i++) if (xx[i] != NULL) free(xx[i]);
  if (U!=NULL) free(U); 
  if (G!=NULL) free(G); 
  if (COV!=NULL) free(COV);
  if (D!=NULL) free(D);
  if (work!=NULL) free(work);
  if (iwork!=NULL) free(iwork);
  if (VT!=NULL) free(VT);
  return Xerror;
}

void do_directGauss(key_type *key, int m, double *res) 
{  
  direct_storage *S;
  long totpnts, i,j,k;
  double *G,*U;  

  S = (direct_storage*) key->meth[m].S;
  totpnts = key->totalpoints;

  U = S->U;// S^{1/2}
  G = S->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)
  for (i=0; i<totpnts; i++) G[i] = GAUSS_RANDOM(1.0);

  switch (S->method) {
  case Cholesky :
    for (k=0, i=0; i<totpnts; i++, k+=totpnts){
      register double dummy;
      dummy =0.0;
      for (j=0; j<=i; j++){
	dummy += G[j] * U[j+k];
      }
      res[i] += dummy; 
    }
    break;
  case SVD :
    for (i=0; i<totpnts; i++){
      register double dummy;
      dummy = 0.0;
      for (j=0, k=i; j<totpnts; j++, k+=totpnts){
	dummy += U[k] * G[j];
      }
      res[i] += dummy; 
    }
    break;
  default : assert(false);
  }
}

