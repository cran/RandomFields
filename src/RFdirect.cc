/*
 Authors
 Martin Schlather, Martin.Schlather@uni-bayreuth.de

 Simulation of a random field by spectral turning bands

 Copyright (C) 2001 Martin Schlather, 

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

#ifdef MATHLIB_STANDALONE
#include "dsvdc.h"
#include "chol.h"  
#else
#include <R_ext/Applic.h>
#endif


#define CHECK if (0)  

InversionMethod DIRECTGAUSS_INVERSIONMETHOD=Cholesky;
Real DIRECTGAUSS_PRECISION=1E-11;
bool DIRECTGAUSS_CHECKPRECISION=false;
int DIRECTGAUSS_MAXVARIABLES=1000;

typedef struct direct_storage {
  InversionMethod method;
  double *U,*G;
} direct_storage;

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

void SetParamDirectGauss(int *action,int *method,int *checkprecision,
			 Real *requiredprecision, int *maxvariables)
{
  switch (*action) {
  case 0 :
    if ((*method<0) || (*method>=(int) NoFurtherInversionMethod)){
      PRINTF("inversion method out of range; ignored\n");
    }
    else {DIRECTGAUSS_INVERSIONMETHOD= (InversionMethod) *method;}
    if ((*checkprecision<0) || (*checkprecision>1)){
      PRINTF("check precion must be 0 or 1; ignored\n");
    }
    else {
      if ((GENERAL_PRINTLEVEL>=3) && 
	  (DIRECTGAUSS_CHECKPRECISION!=(bool) *checkprecision))
	PRINTF("Note: precision setting works only with Cholesky up to now\n");
      DIRECTGAUSS_CHECKPRECISION=(bool) *checkprecision; 
    }
   DIRECTGAUSS_PRECISION=*requiredprecision;
   if (DIRECTGAUSS_PRECISION<=0) 
     PRINTF("Warning! Non positive value for precision. Algorithm will probably fail.");
   DIRECTGAUSS_MAXVARIABLES=*maxvariables;
   break;
  case 1 :
    *method = (int) DIRECTGAUSS_INVERSIONMETHOD;
    *checkprecision = (int) DIRECTGAUSS_CHECKPRECISION;
    *requiredprecision = DIRECTGAUSS_PRECISION;
    *maxvariables = DIRECTGAUSS_MAXVARIABLES;
    if (GetNotPrint) break;
  case 2 : 
     PRINTF("\nDirect matrix inversion\n=======================\ninversion method=%d\nmaximum number of variables %d\ncheck precision? %d\n",
	    DIRECTGAUSS_INVERSIONMETHOD,DIRECTGAUSS_MAXVARIABLES,
	    DIRECTGAUSS_CHECKPRECISION);
    if (DIRECTGAUSS_CHECKPRECISION) {
      PRINTF("\nprecision=%e(works only with Cholesky up to now)",
	     DIRECTGAUSS_PRECISION);}

   break;
  default : PRINTF(" unknown action\n"); 
  }
}


int CHOLpreciseenough(double *COV, double *U,long nrow) 
{
  double sum;
  int i,j,k;
  long segment,endfor;
  sum = 0;
  endfor = nrow * nrow;
  for (segment=i=0; i<nrow; i++,segment+=nrow) {    
    for (j=segment; j<endfor; j+=nrow) {     
      register Real dummy =0;
      for (k=0; k<=i; k++) {
	dummy += U[segment+k] * U[j+k];
      }
      sum += fabs(COV[i+j]-dummy);
    }
  }
  if (sum>DIRECTGAUSS_PRECISION) {PRINTF("Cholesky, difference=%e\n",sum);} 
  if (GENERAL_PRINTLEVEL>=3) PRINTF("\nChol prec %e\n",sum);
  return sum<DIRECTGAUSS_PRECISION;
}


int init_directGauss(key_type *key) 
{
  double *G,*COV,*U,*V,*e,*D,*xx[MAXDIM];
  long d,i,j,k,halfn, error, job=11,longrow,segment;
  bool first,freexx;
  InversionMethod method;
  covfct cov;



  assert(key->active==false);
  assert((key->covnr>=0) &&(key->covnr<currentNrCov));
  cov = key->cov->cov; /// NOT LOCAL, NO FRACTAL BROWNIAN (YET)
  assert(cov!=NULL); // I do not know an example yet, where cov is unknown, 
  //                    but should be simulated
  assert(key->param[VARIANCE]>=0.0);
  assert(key->param[SILL]>=key->param[VARIANCE]); 
  assert(fabs(key->param[SCALE]*key->param[INVSCALE]-1.0) < EPSILON);
  assert((key->S==NULL) && (key->destruct==NULL));
  G=COV=U=V=e=D=NULL;// necessary!

  freexx = false; // should xx[i] be freed (this is the case if key->grid==TRUE 
  //                 and no eror has occured 

  if ((key->S =malloc(sizeof(direct_storage)))==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  ((direct_storage*)key->S)-> U=NULL;
  ((direct_storage*)key->S)-> G=NULL;
  key->destruct = direct_destruct;
  if(key->cov->check!=NULL) 
    if(error=(key->cov->check(key))) goto ErrorHandling;

  longrow=key->totallength;    
  if (longrow>DIRECTGAUSS_MAXVARIABLES) {
    error=ERRORMETHODNOTALLOWED;goto ErrorHandling;
  }
   
  if ((COV =(double *) malloc(sizeof(double) * longrow * longrow))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  }

  /* calculate covariance matrix */
  if (key->grid) {
    int index[MAXDIM],dimM1=key->dim-1;
    double *yy[MAXDIM],p[MAXDIM];
    for(d=0;d<=dimM1;d++) {
      xx[d]=NULL;
      yy[d]=key->x[d];
   }
    freexx = true;
    
    // generate all the grid coordinates exlicitely!

    for (i=0;i<=dimM1;i++) { 
      if ((xx[i]=(Real*) malloc(sizeof(Real)* longrow))==NULL){
	error=ERRORMEMORYALLOCATION;goto ErrorHandling;
      }
    }
    for (d=0;d<=dimM1;d++) {index[d]=0; xx[d][0]=p[d]=0.0;} 
    // intialisation of index and p, and defining x[d][i=0]
    // index preciser than the real coordinate values (which are created in p[])
    // index only used to compare with length*
    for (i=1;i<key->totallength;i++) {
      d=0;  //d=dimM1; if x should run the slowest;
      index[d]++; p[d]+=yy[d][XSTEP];
      while (index[d]>=key->length[d]) { 
	index[d]=0; p[d]=0.0; 
	d++;   // d--;
	assert(d<=dimM1); // assert(d>=0);
	index[d]++; p[d]+=yy[d][XSTEP];
      }
      for (d=0; d<=dimM1; d++) {xx[d][i]=p[d];}
    }
  } else { 
    for (d=0; d<key->dim; d++) xx[d]=key->x[d];
  }
  k = 0;
  for (i=0;i<longrow;i++) {
    k += i;
    for (segment=i* longrow+i,j=i;j<longrow;j++) {
      register double dx,dy;
      dx =0;
      for (d=0;d<key->dim;d++) {dy=xx[d][i]-xx[d][j]; dx+=dy*dy;}
      COV[k]=COV[segment] = cov(sqrt(dx),key->param);
      k++;
      segment += longrow;
    }
  }
 
  if (key->grid) {for (i=0;i<key->dim;i++) {free(xx[i]);xx[i]=NULL;}}
  if ((U =(double *) malloc(sizeof(double) * longrow * longrow))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  }
  method = DIRECTGAUSS_INVERSIONMETHOD;
  //for SVD intermediate results, and memory space for do_directGauss:
  halfn=(longrow/2)*2;  if (halfn<longrow) { halfn += 2;}
  if ((G = (double *)  malloc(sizeof(double) * halfn))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  } 

  if (GENERAL_PRINTLEVEL>=20) {
    PRINTF("\n COV\n");   
    for (i=0;i<longrow;i++) { for (j=i;j<longrow* longrow;j+=longrow) { 
      PRINTF(" %f ", COV[j]);} PRINTF("\n"); 
    }
  }

  first = true;
  error = 1;
  while (error && (method!=NoFurtherInversionMethod)) {
    if (GENERAL_PRINTLEVEL>=3) PRINTF("method=%d\n",method);
    switch (method) {
    case Cholesky :  
#ifdef MATHLIB_STANDALONE
      chol_(COV,&longrow,&longrow,U,&error);
#else
      { 
	int row,err; row=longrow;
        F77_CALL(chol)(COV,&row,&row,U,&err);
	error = err;
      }
#endif
      if (error==0) {
	if ((DIRECTGAUSS_CHECKPRECISION) && !CHOLpreciseenough(COV,U,longrow)) {
	  error=ERRORPRECISION;
	}
      } else error=ERRORDECOMPOSITION;
      if ((error) && (GENERAL_PRINTLEVEL>=1)) 
	PRINTF("Cholesky decomposition error %d\n",error);
      break;
    case SVD :
	if ((V =(double *) malloc(sizeof(double) * longrow * longrow))==NULL){
	  error=ERRORMEMORYALLOCATION;goto ErrorHandling;
	}
	if ((D =(double *) malloc(sizeof(double) * longrow))==NULL){
	  error=ERRORMEMORYALLOCATION;goto ErrorHandling;
	}
        if ((e = (double *) malloc(sizeof(double) * longrow))==NULL){
	  error=ERRORMEMORYALLOCATION;goto ErrorHandling;
	}
#ifdef MATHLIB_STANDALONE
	dsvdc_(COV,&longrow,&longrow,&longrow,D,e,U,&longrow,V,&longrow,G,
	       &job,&error);
#else
      { 
	int row,err,jobint; row=longrow; jobint=job;
 	F77_CALL(dsvdc)(COV,&row,&row,&row,D,e,U,&row,V,&row,G,&jobint,&err);
	error = err;
      }
#endif
	if (error!=0) {
	  if (GENERAL_PRINTLEVEL>0)PRINTF("ERROR-Code SVD=%d\n",error); 
	  error=ERRORDECOMPOSITION;
	  goto ErrorHandling;
	}
	
	if (GENERAL_PRINTLEVEL>=10) {
	  PRINTF("\n U\n");   
	  for (i=0;i<longrow;i++) { for (j=i;j<longrow* longrow;j+=longrow) {
	    PRINTF(" %f ", U[j]);} PRINTF("\n"); }
	  PRINTF("\n V\n");   
	  for (i=0;i<longrow;i++) { for (j=i;j<longrow* longrow;j+=longrow) { 
	    PRINTF(" %f ", V[j]);} PRINTF("\n"); }
	  PRINTF("\n D\n");   
	  for (i=0;i<longrow;i++) {PRINTF(" %f ", D[i]);}
	}

 	free(e); e=NULL; free(V); V=NULL; // here free(COV) is already possible; 
	//                                   for debugging reasons it is postponed
 
	/* caculate SQRT of covariance matrix */
	for (k=0,j=0;j<longrow;j++) {
	  register double dummy;
	  dummy = sqrt(D[j]);
	  for (i=0;i<longrow;i++) {
	    U[k++] *= dummy;
	  }
	}
	
	if (GENERAL_PRINTLEVEL>=10) { 
	  PRINTF("\nSqrt COV \n");
	  for (i=0;i<longrow;i++) { for (j=i;j<longrow* longrow;j+=longrow) { 
	    PRINTF(" %f ", U[j]);} PRINTF("\n");
	  }
	  PRINTF("\n");
	}
	break;
    default : assert(false);
    }
    if (error) {
      if (first) {method=(InversionMethod) 0;} 
      else {method= (InversionMethod) (1+(int) method);}
      if (method==DIRECTGAUSS_INVERSIONMETHOD) {
	method= (InversionMethod) (1+(int) method);
      }
    }
  }
  if (error) { goto ErrorHandling; }  
  free(COV);
  ((direct_storage*)key->S)->U=U;
  ((direct_storage*)key->S)->method=method;
  ((direct_storage*)key->S)->G=G;
  if (D!=NULL) free(D);
  if (e!=NULL) free(e);
  if (V!=NULL) free(V);
  if (freexx) {for (i=0;i<key->dim;i++) {if (xx[i]!=NULL) free(xx[i]);}}
  return 0;
  
  ErrorHandling:
  //  PRINTF("RFdirect error\n");
  if (COV!=NULL) free(COV);
  if (freexx) {for (i=0;i<key->dim;i++) {if (xx[i]!=NULL) free(xx[i]);}}
  if (U!=NULL) free(U); 
  if (D!=NULL) free(D);
  if (G!=NULL) free(G); 
  if (e!=NULL) free(e);
  if (V!=NULL) free(V);
  assert(key->destruct!=NULL);
  key->destruct(&key->S);
  key->destruct=NULL;
  return error;
}  


void do_directGauss(key_type *key, Real *res END_WITH_GSLRNG) 
{  
  long longrow,halfn,i,j,k;
  double *G,*U;  

#ifdef RF_GSL
  assert(RANDOM!=NULL);
#endif
  assert(key->active);

  longrow=key->totallength;    
  halfn=(longrow/2)*2;  if (halfn<longrow) { halfn += 2;}
  U =  ((direct_storage*)key->S)-> U; // S^{1/2}
  G = ((direct_storage*)key->S)-> G;  // only the memory space is of interest 
  //                                    (stored to avoid allocation errors here)

  for (i=0;i<halfn;) {
    register double UU,VV;
    // PRINTF("i %d\n",i); 
    UU=  SQRTTWO * sqrt(-log(1.0 - UNIFORM_RANDOM));
    VV = TWOPI * UNIFORM_RANDOM; 

    G[i++] = UU * sin(VV);  
    G[i++] = UU * cos(VV); 
  }

 
  switch (((direct_storage*)key->S)->method) {
  case Cholesky :   
    for (k=0,i=0; i<longrow; i++,k+=longrow){
      register double dummy;
      dummy =0.0;
      for (j=0;j<=i;j++){
	dummy += G[j] * U[j+k];	
      }
      res[i] = (Real) dummy;
    } 
    break;
  case SVD :
    for (i=0;i<longrow;i++){
      register double dummy;
      dummy =0;
      for (j=0,k=i;j<longrow;j++,k+=longrow){
	dummy += U[k] * G[j];
      }
      res[i] = (Real) dummy;
    }
    break;
  default : assert(false);
  }

  if (!GENERAL_STORING) { 
    assert(key->destruct!=NULL);
    key->destruct(&key->S);
    key->destruct=NULL;
    key->active=false;}
}









