/*
 Authors
 Martin Schlather, martin.schlather@cu.lu

 Simulation of a random field by spectral turning bands

 Copyright (C) 2001 -- 2004 Martin Schlather, 

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

#ifndef MATHLIB_STANDALONE
#include <R_ext/Applic.h>
#endif


#define CHECK if (0)  

InversionMethod DIRECTGAUSS_INVERSIONMETHOD=Cholesky;
Real DIRECTGAUSS_PRECISION=1E-11;
bool DIRECTGAUSS_CHECKPRECISION=false;
int DIRECTGAUSS_MAXVARIABLES=1800; // see RFsimu.h

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
  if (GENERAL_PRINTLEVEL>=3) PRINTF("\nChol prec %]e\n",sum);
  return sum<DIRECTGAUSS_PRECISION;
}


double *CCOV;
int nn;

void getCov(int *n, double *COV) {
  assert(false);
  int i;
  if (*n==0) { *n=nn; return;}
  for (i=0; i<*n * *n; i++) COV[i]=CCOV[i];
}

int internal_init_directGauss(direct_storage **S, bool grid,
			      int spatialdim, bool Time, double** x, int* length,
			      int totpnts, double *T, CovFctType CovFct,
			      int *covnr, int *op, param_type param,
			      int actcov, bool anisotropy) {
  double *G,*COV,*U,*V,*e,*D;
  Real *xx[MAXDIM];
  long d,i,j,k,halfn, Xerror, job=11,segment;
  bool freexx;
  InversionMethod method;
  int timespacedim;

  timespacedim = spatialdim + (int) Time;
  G=COV=U=V=e=D=NULL;// necessary!
  freexx = false; // should xx[i] be freed (this is the case if key->grid==TRUE 
  //                 and no eror has occured 

  if ((*S = (direct_storage*) malloc(sizeof(direct_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  (*S)->U = (*S)->G = NULL;

  if ((COV =(double *) malloc(sizeof(double) * totpnts * totpnts))==NULL){
    Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;
  }

//////////////////////////
//  if ((CCOV =(double *) malloc(sizeof(double) * totpnts * totpnts))==NULL){
//    Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;
//  }
//  nn = totpnts;
//////////////////////////

  /* calculate covariance matrix */
  if (grid) {
    int index[MAXDIM];
    Real *yy[MAXDIM];
    double p[MAXDIM];
    for(d=0; d<timespacedim; d++) {
      xx[d]=NULL;
      yy[d]=x[d];
    }
    freexx = true;
    
    // generate all the grid coordinates exlicitely!
    for (i=0; i<timespacedim; i++) { 
      if ((xx[i]=(Real*) malloc(sizeof(Real)* totpnts))==NULL){
	Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
    }
    for (d=0; d<timespacedim; d++) {index[d]=0; xx[d][0]=p[d]=0.0;} 
    // intialisation of index and p, and defining x[d][i=0]
    // index preciser than the real coordinate values (which are created in p[])
    // index only used to compare with length*
    for (i=1; i<totpnts; i++) {
      d=0;  //d=timespacedim-1; if x should run the slowest;
      index[d]++; 
      p[d]+=yy[d][XSTEP];
      while (index[d]>=length[d]) { 
	index[d]=0; p[d]=0.0; 
	d++;   // d--;
	assert(d<timespacedim); // assert(d>=0);
	index[d]++; p[d]+=yy[d][XSTEP];
      }
      for (d=0; d<timespacedim; d++) {xx[d][i]=p[d];}
    }
  } else { // no grid
    int t; 
    if (Time) {
      int endfor;
      Real time;
      for(d=0; d<=timespacedim; d++) xx[d]=NULL;
      freexx = true;
      for (i=0; i<timespacedim; i++) { 
	if ((xx[i]=(Real*) malloc(sizeof(Real)* totpnts))==NULL){
	  Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
	}
      }
      endfor = length[spatialdim]; 
      for (t=0, i=0,time=T[XSTART]; t<endfor; t++, time+=T[XSTEP]) 
	for (j=0; j<length[0]; j++,i++) {
	  for (d=0; d<spatialdim; d++) xx[d][i]=x[d][j];
	  xx[spatialdim][i] = time;
	}
    } else { 
      for (d=0; d<spatialdim; d++) xx[d]=x[d];
    }
  }

  k = 0;
  for (i=0;i<totpnts;i++) {
    Real y[MAXDIM];
    k += i;
    for (segment=i* totpnts+i,j=i;j<totpnts;j++) {
      for (d=0; d<timespacedim; d++) {
	y[d]=xx[d][i] - xx[d][j];
      }
//      printf(" %d (%d) ", i, totpnts);

      COV[k]=COV[segment] = CovFct(y, timespacedim, covnr, op, param, 
				   actcov, anisotropy);
      k++;
      segment += totpnts;
    }
  }

//  for (k=0; k<totpnts; k++) {
//    printf(" %d ", k);
//    for (i=0; i<totpnts; i++)
//      if (COV[k + i * totpnts] != COV[i + k * totpnts]) 
//	printf("cov %d %d %f %f\n",
//	       k, i, COV[k + i * totpnts], COV[i + k * totpnts]);
//  }

  if (freexx) {for (i=0; i<timespacedim; i++) {free(xx[i]); xx[i]=NULL;}}
  if ((U =(double *) malloc(sizeof(double) * totpnts * totpnts))==NULL){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  method = DIRECTGAUSS_INVERSIONMETHOD;
  //for SVD intermediate results, and memory space for do_directGauss:
  halfn=(totpnts/2)*2;  if (halfn<totpnts) { halfn += 2;}
  if ((G = (double *)  malloc(sizeof(double) * halfn))==NULL){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 

  Xerror = 1;
  switch (method) {
      case Cholesky :  
	if (GENERAL_PRINTLEVEL>=3) PRINTF("method=Cholesky\n");
	{ 
	  int row,err; row=totpnts;
	  F77_CALL(chol)(COV,&row,&row,U,&err);
	  Xerror = err;
	}
	if (Xerror==NOERROR) {
	  if ((DIRECTGAUSS_CHECKPRECISION) && !CHOLpreciseenough(COV,U,totpnts)) {
	    Xerror=ERRORPRECISION;
	  }
	} else Xerror=ERRORDECOMPOSITION;
	if (Xerror==NOERROR) break;
	else if (GENERAL_PRINTLEVEL>=1)
	  PRINTF("Error code F77_CALL(chol) = %d\n",Xerror);

	// try next method :
	
      case SVD :

///////////////////////////
// {
//   long t2;
//   double dev, max, UV;
//   t2 = totpnts * totpnts;
//   for (i=0; i<t2; i++) CCOV[i]=COV[i];
// }

	method = SVD; // necessary if the value of method has been Cholesky.
	//               originally
	if (GENERAL_PRINTLEVEL>=3) PRINTF("method=SVD\n");
	if ((V =(double *) malloc(sizeof(double) * totpnts * totpnts))==NULL){
	  Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;
	}
	if ((D =(double *) malloc(sizeof(double) * totpnts))==NULL){
	  Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;
	}
	if ((e = (double *) malloc(sizeof(double) * totpnts))==NULL){
	  Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;
	}
	{ 
	  int row,err,jobint; row=totpnts; jobint=job;
	  F77_CALL(dsvdc)(COV,&row,&row,&row,D,e,U,&row,V,&row,G,&jobint,&err);
	  Xerror = err;
	}
	if (Xerror!=NOERROR) {
	  if (GENERAL_PRINTLEVEL>1) 
	    PRINTF("Error code F77_CALL(dsvdc) = %d\n",Xerror); 
	  Xerror=ERRORDECOMPOSITION;
	  goto ErrorHandling;
	}
	
	if (GENERAL_PRINTLEVEL>=2) {
	  long t2;
	  double dev, max, UV;
	  t2 = totpnts * totpnts;
	  dev = 0.0;
	  max = 0;
	  for (i=0; i<t2; i++) {
	    UV = fabs(U[i] - V[i]);
	    if (UV>0.01) printf("%d %f %f %f\n", i, UV, U[i], V[i]);
	    dev += UV;
	    if (max<UV) max=UV;
	  }
	  PRINTF("totpnts %d %d\n", totpnts, t2); 
	  PRINTF("total asymmetry U/V: sum=%f, max=%f\n", dev, max); 
	  if (GENERAL_PRINTLEVEL>=6) { 
	    PRINTF("\n D max=%f, D min=%f \n", D[0], D[totpnts-1]);
	  }
	}
	
	free(e); e=NULL; free(V); V=NULL; // here free(COV) is already possible; 
	//                                  for debugging reasons it is postponed
	
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
  free(COV);
  ((direct_storage*)*S)->U=U;
  ((direct_storage*)*S)->method=method;
  ((direct_storage*)*S)->G=G;
  if (D!=NULL) free(D);
  if (e!=NULL) free(e);
  if (V!=NULL) free(V);
  if (freexx) 
    for (i=0; i<timespacedim; i++) if (xx[i]!=NULL) free(xx[i]);
  return NOERROR;
  
  ErrorHandling:
  if (COV!=NULL) free(COV);
  if (freexx) 
    for (i=0; i<timespacedim; i++) if (xx[i]!=NULL) free(xx[i]);
  if (U!=NULL) free(U); 
  if (D!=NULL) free(D);
  if (G!=NULL) free(G); 
  if (e!=NULL) free(e);
  if (V!=NULL) free(V);
  return Xerror;
}  


int init_directGauss(key_type *key, int m) 
{
  param_type param;
  long Xerror;

  SET_DESTRUCT(direct_destruct);
  FIRSTCHECK_COV(Direct,cov,param); // multiply, covnr, actcov defined
  { 
    int v, timespacedim, start_param[MAXDIM], index_dim[MAXDIM];
    bool no_last_comp;
    cov_fct *cov;
    for (v=0; v<actcov; v++) {
      GetTrueDim(key->anisotropy, key->timespacedim, param[v],
		 &timespacedim, &no_last_comp, start_param, index_dim);
      cov = &(CovList[covnr[v]]);
      if ((key->Time) && no_last_comp && (cov->isotropic==SPACEISOTROPIC))
	{ Xerror= ERRORWITHOUTTIME; goto ErrorHandling;}
      else if ((cov->check!=NULL) &&
	       ((Xerror=cov->check(param[v], timespacedim, Direct)))!=NOERROR)
	goto ErrorHandling;
    }
  }

  if (key->totalpoints>DIRECTGAUSS_MAXVARIABLES) {
    Xerror=ERRORMETHODNOTALLOWED; goto ErrorHandling;
  }
  return(internal_init_directGauss((direct_storage**) (&(key->S[m])),  
				   key->grid, key->spatialdim, key->Time, 
				   key->x, key->length, key->totalpoints,
				   key->T, CovFct,
				   covnr, multiply, param, actcov,
				   key->anisotropy)); 
 ErrorHandling: 
  return Xerror;
}


void internal_do_directGauss(direct_storage *S, bool add, long totpnts, 
			     double *res) {
  long i,j,k;
  double *G,*U;  

  U = S->U;// S^{1/2}
  G = S->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)
  for (i=0; i<totpnts; i++) G[i]=GAUSS_RANDOM(1.0);
#define RESULT(OP)\
  switch (S->method) {\
  case Cholesky :\
    for (k=0,i=0; i<totpnts; i++,k+=totpnts){\
      register double dummy;\
      dummy =0.0;\
      for (j=0;j<=i;j++){\
	dummy += G[j] * U[j+k];\
      }\
      res[i] OP (Real) dummy;\
    }\
    break;\
  case SVD :\
    for (i=0;i<totpnts;i++){\
      register double dummy;\
      dummy =0;\
      for (j=0,k=i;j<totpnts;j++,k+=totpnts){\
	dummy += U[k] * G[j];\
      }\
      res[i] OP (Real) dummy;\
    }\
    break;\
  default : assert(false);\
  }
  
  if (add) {RESULT(+=)} else {RESULT(=)}
}


void do_directGauss(key_type *key, bool add, int m, Real *res) 
{  
  assert(key->active);
  internal_do_directGauss(((direct_storage*) key->S[m]), add, key->totalpoints, 
			  res);
}


void XXinternal_do_directGauss(direct_storage *S, bool add, long totpnts, 
			     double *res) {
  long i,j,k;
  double *G,*U;  


  printf("XX %d \n", totpnts);


  U = S->U;// S^{1/2}
  G = S->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)

  printf("X\n");
  for (i=0; i<totpnts; i++) G[i]=GAUSS_RANDOM(1.0);
  printf("sXX\n");

#define XRESULT(OP)\
  switch (S->method) {\
  case Cholesky :\
  printf("sXX\n");\
    for (k=0,i=0; i<totpnts; i++,k+=totpnts){\
      register double dummy;\
      dummy =0.0;\
      for (j=0;j<=i;j++){\
  	dummy += G[j] * U[j+k];\
      }\
     res[i] OP (Real) dummy;\
    }\
    break;\
  case SVD :\
    for (i=0;i<totpnts;i++){\
      register double dummy;\
      dummy =0;\
      for (j=0,k=i;j<totpnts;j++,k+=totpnts){\
	dummy += U[k] * G[j];\
      }\
      res[i] OP (Real) dummy;\
    }\
    break;\
  default : assert(false);\
  }
  
  if (add) {XRESULT(+=)} else {XRESULT(=)}
}


