/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Simulation of a random field by sequential method

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

#include <math.h>
#include <stdio.h>  
#include <stdlib.h>
#include <assert.h>
#include "RF.h"
#include <R_ext/Lapack.h>
//#include <R_ext/Linpack.h>



bool debugging = true;

void sequential_destruct(void ** S)
{
  if (*S!=NULL) {
    sequential_storage *x;
    x = *((sequential_storage**)S);    
    if (x->U11!=NULL) free(x->U11); 
    if (x->U22!=NULL) free(x->U22); 
    if (x->MuT!=NULL) free(x->MuT); 
    if (x->G!=NULL) free(x->G);
    if (x->Cov21!=NULL) free(x->Cov21);
    if (x->Inv22!=NULL) free(x->Inv22);
    if (x->res0!=NULL) free(x->res0);
    free(*S);
    *S = NULL;
  }
}

void SetParamSequential(int *action, int *maxvariables, int *back, int *initial)
{
  sequ_param  *gp = &(GLOBAL.sequ);
  if (*action) {
    gp->max  = *maxvariables;
    gp->back = *back < 1 ? 1 : *back;
    gp->initial= *initial;
  } else {
    *maxvariables  = gp->max;
    *back  = gp->back;
    *initial = gp->initial; 
   }
}

// (U1, U2) ~ N(0, S) =>
// U1 | U2 ~ N(S_{12} S_{22}^{-1} u2, S_11 - S12 S_22^{-1} S_21)

// start mit S_22; dann nur zeilenweise + Einschwingen

int init_sequential(method_type *meth){
  cov_model *cov = meth->cov;
  cov_fct *C=CovList + cov->nr;
  covfct cf = C->cov;
  location_type *loc = meth->loc;
  globalparam *gp = meth->gp;
  sequ_param* lp = &(gp->sequ);
  long Err=NOERROR, totpnts=0, totpntsSQ, spatialpntsSQ=0, i,
      spatialpntsSQback;
  double 
    *xx = NULL,
    *G = NULL, 
    *U11 = NULL, 
    *COV21 = NULL,
    *MuT = NULL, 
    *U22 = NULL, 
      *Inv22 = NULL;
  res_type
    *res0 = NULL;
  int withoutlast, d, endfor, l, store,
    dim=cov->tsdim,
    spatialdim = dim - 1,
    spatialpnts=0,
    back=0, 
    initial=0,
    PL = gp->general.printlevel;
  sequential_storage* S = NULL;
  bool storing = gp->general.storing;
 
  assert(meth->S == NULL); // ist in simu.cc sichergestellt, schadet aber nicht
  SET_DESTRUCT(sequential_destruct);
  bool Time = loc->Time;

  if (meth->caniso != NULL || meth->cproj!=NULL || meth->cscale != 1.0) {
      Err = ERRORPREVDOLLAR;
      goto ErrorHandling;
  }
   
  if (!loc->grid && !Time) {
    Err = ERRORLASTGRID; 
    goto ErrorHandling;
  }  
  if (C->implemented[Sequential] != IMPLEMENTED) {
    Err=ERRORNOTDEFINED; 
    goto ErrorHandling;
  }
  
  if (cov->vdim > 1) {
      Err=ERRORNOMULTIVARIATE; 
      goto ErrorHandling;   
  }

  spatialpnts = loc->totalpoints / loc->length[spatialdim];

  spatialpntsSQ = spatialpnts * spatialpnts;
  totpnts = (lp->back < 1 ? 1 : lp->back) * spatialpnts;
  if (totpnts > lp->max) {
      sprintf(ERRORSTRING_OK, 
	    "number of points less than RFparameters()$sequential.max (%d)",
	      (int) lp->max);
       sprintf(ERRORSTRING_WRONG,"%d * %d = %ld", 
	       (lp->back < 1 ? 1 : lp->back), 
	       spatialpnts, totpnts);
       Err=ERRORCOVFAILED; goto ErrorHandling;
  }

  Err=ERRORMEMORYALLOCATION;  
  if ((meth->S = malloc(sizeof(sequential_storage))) == NULL)
    goto ErrorHandling;
  S = (sequential_storage*) meth->S;

  back = lp->back;
  if (back < 1) {
    back = lp->max / spatialpnts;
  }
  totpnts = back * spatialpnts;
  totpntsSQ = totpnts * totpnts;
  spatialpntsSQback = totpnts * spatialpnts;

  initial = lp->initial;
  if (initial < 0) initial = back - initial;

  if ((U22 =(double *) malloc(sizeof(double) * totpntsSQ))==NULL)
      goto ErrorHandling;
  if ((Inv22 =(double *) malloc(sizeof(double) * totpntsSQ))==NULL)
      goto ErrorHandling;
  if ((COV21 =(double *) malloc(sizeof(double) * spatialpntsSQback))==NULL)
    goto ErrorHandling;
  if ((U11 =(double *) malloc(sizeof(double) * spatialpntsSQ))==NULL)
    goto ErrorHandling;
  if ((MuT =(double *) malloc(sizeof(double) * spatialpntsSQback))==NULL)
    goto ErrorHandling;
  // memory space for do_sequential:
  if ((G = (double *) malloc(sizeof(double) * totpnts))==NULL)
      goto ErrorHandling;
  if ((res0 = (res_type *) malloc(sizeof(res_type) * 
				(totpnts + spatialpnts * initial))) ==NULL)
      goto ErrorHandling;
  S->U22 = S->U11 = S->MuT = S->G = S->Cov21 = S->Inv22 = NULL;
  S->res0 = NULL;

  if (loc->length[spatialdim] <= back) {
    Err = ERRORTRIVIALTIMEDIRECTION;
    goto ErrorHandling;
  } else Err = NOERROR;

  /* ************************* */
  /* create matrix of explicit */
  /*       x-coordinates       */
  /* ************************* */


  store = loc->length[spatialdim];
  loc->length[spatialdim] = back;
  if (loc->grid) {
     expandgrid(loc->xgr, loc->length, 
		&xx, 
		loc->timespacedim);
     if (!loc->Time) meth->space = meth->sptime;
  } else {
    xtime2x(loc->x, loc->length[0], loc->T, loc->length[spatialdim],
	    &xx, loc->timespacedim);
    assert(loc->Time);
  }
  loc->length[spatialdim] = store;
 
//  PrintModelInfo(cov); assert(false);


  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  long j, k, segment;
  int  row, err, k0, k1, k2;
  double *y;
  y = (double*) malloc(loc->timespacedim * sizeof(double));

  // *** S22

//  printf

  if (PL>4) PRINTF("covariance matrix...\n");
  k = 0;
  for (k0 =i=0; i<totpnts; i++, k0+=dim) {
    k += i;
    for (segment=i* totpnts+i, j=i, k2=k0; j<totpnts; j++) {
      k1 = k0;
      for (d=0; d<dim; d++) {
	y[d] = xx[k1++] - xx[k2++];
//	printf("%f ", y[d]);
      }
      cf(y, cov, U22 + segment);
//      printf("cov=%f\n", U22[segment]);
      U22[k++] = U22[segment];
      segment += totpnts;
    }
  }

  /* Orndung von U22: 
      back ....1
   back
   ...
   1

  Ordnung von COV21
  back+1
  back
  ...
  2
  */

  // *** S11 und S21
  j = k = 0;
  withoutlast = totpnts - spatialpnts;
  for (i=withoutlast * totpnts; i<totpntsSQ; ) {
    j += spatialpnts;
    endfor = j + withoutlast;
    for (; j<endfor; ) {
	COV21[j++] = U22[i++];
    }
    endfor = k + spatialpnts;
    for (; k<endfor; ) U11[k++] = U22[i++];
  }
  

  // *** S21 rest
  k = 0;
  y[spatialdim] = loc->xgr[spatialdim][XSTEP] * back;
  for (k0 =i=0; i<spatialpnts; i++, k0+=dim) { // t_{n+1}
      for (k2=j=0; j<spatialpnts; j++) { // t_1
	k1 = k0;
	for (d=0; d<spatialdim; d++) {
	    y[d] = xx[k1++] - xx[k2++];
	}
	cf(y, cov,  COV21 + k);

	k++; k2++;
    }
    k += withoutlast;
  }

  free(y);
  y = NULL;

//  for (i=0; i<withoutlast * spatialpnts; i++) {
//    printf("%3.4f ", COV21[i]);
//  }
//  assert(false);

/*

source("~/R/RF/RandomFields/tests/source.R")


model <- list(
              "+",
              list("$", aniso=matrix(1:4, ncol=2),
                   list("exp")
                   ),
              list("$", var=0.,
                   list("nugget")
                   )
              )
model <- list("$", aniso=matrix(c(1,0,0,1.1), ncol=2),
                   list("exp")
                   )
#  library(RandomFields, lib"~/R/OLD");model <- list(list(model="exp", aniso=c(1,0,0,1.1), var=1))

RFparameters(Print=5, exactness=FALSE)
x <- c(0,1,0.4) # --> err 66
x <- c(0,1,0.12) # -> U11 neg ## done !!
xx <- seq(x[1],x[2],x[3])
z <- GaussRF(x, x, model=model, gridtri=TRUE, grid=TRUE, Print=6, me="se",
             every=10 ) 
*/

/*
  PRINTF("U11\n");
  for (i=0; i<spatialpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      PRINTF("%3.2f ", U11[i + j * spatialpnts]);
   }
    PRINTF("\n");
  }

  PRINTF("S22 %d %d\n", totpnts,spatialpnts );
  for (i=0; i<totpnts; i++) {
    PRINTF("%2d: ", i);
    for (j=0; j<totpnts; j++) {
      PRINTF("%3.2f ", U22[i + j * totpnts]);
   }
    PRINTF("\n");
  }
*/


  /* ********************* */
  /* matrix compositions  */
  /* ********************* */

  // *** sqrt S 22 
  if (PL>4) PRINTF("Cholesky decomposition of Cov22...\n");
  row=totpnts;
  // dchdc destroys the input matrix; upper half of U22 contains result!
  // dpotrf	F77_CALL(dpotrf)("Upper", &m, REAL(ans), &m, &i);
  //  assert(false);
  F77_CALL(dpotrf)("Upper", &row, U22, &row, &err);
//  assert(false);
  // F77_NAME(dchdc)(U22, &row, &row, G, NULL, &choljob, &err);
  if (err!=NOERROR) {
    if (PL>2)
      PRINTF("Error code F77_CALL(dpotrf) = %d\n", err);
    Err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  } 

  // for (i=0; i<withoutlast * spatialpnts; i++) {
//    PRINTF("%3.4f ", COV21[i]);
//  }
//  assert(false);
  
  // *** inverse of S 22 
  if (PL>4) PRINTF("inverse of Cov22...\n");
  memcpy(Inv22, U22, sizeof(double) * totpntsSQ);
  F77_CALL(dpotri)("Upper", &row, Inv22, &row, &err);
  if (err!=NOERROR) {
    if (PL>2)
      PRINTF("Error code F77_CALL(dpotri) = %d\n", err);
    Err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  }
  for (k=i=0; i<totpnts; i++) {
    k += i;
    for (segment=i* totpnts+i; segment<totpntsSQ; segment += totpnts) {
	Inv22[k++] = Inv22[segment]; 
//	PRINTF("sg %d %f %d \n", k, Inv22[segment] , segment);
    }
  }


/*
  PRINTF("C21\n");
  for (i=0; i<totpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      PRINTF("%3.2f ", COV21[i + j * totpnts]);
   }
    PRINTF("\n");
  }
*/  


  // *** MuT

  if (PL>4) PRINTF("calculating matrix for the mean...\n");
  l = 0;
  for (i=0; i<spatialpntsSQback; i+=totpnts) {
    for (j=0; j<totpntsSQ; j+=totpnts) {
	double dummy;
	dummy = 0.0;
	for (k=0; k<totpnts; k++) dummy += COV21[i + k] * Inv22[j + k];
	MuT[l++] = dummy;
    }
  }

/*
  PRINTF("MuT\n");
  for (i=0; i<totpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      PRINTF("%3.2f ", MuT[i + j * totpnts]);
   }
    PRINTF("\n");
  }
*/
  //  assert(false);
 
  // *** C11
  if (PL>4) PRINTF("calculating cov matrix...\n");
  l = 0;
  for (i=0; i<spatialpntsSQback; i+=totpnts) {
    for (j=0; j<spatialpntsSQback; j+=totpnts) {
	double dummy;
	dummy = 0.0;
	for (k=0; k<totpnts; k++) dummy += MuT[i + k] * COV21[j + k];
	U11[l++] -= dummy;
    }
  }

/*
  PRINTF("U11\n");
  for (i=0; i<spatialpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      PRINTF("%3.2f ", U11[i + j * spatialpnts]);
   }
   PRINTF("\n");
  }
*/


  // *** sqrt C11
  if (PL>4) PRINTF("Cholesky decomposition of Cov11...\n");
  row=spatialpnts;
  // dchdc destroys the input matrix; upper half of U22 contains result!
  // dpotrf	F77_CALL(dpotrf)("Upper", &m, REAL(ans), &m, &i);
  F77_CALL(dpotrf)("Upper", &row, U11, &row, &err);

/*
  for (i=0; i<spatialpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      PRINTF("%3.2f ", U11[i + j * spatialpnts]);
   }
    PRINTF("\n");
  }
*/

  if (err!=NOERROR) {
    if (PL>2)
      PRINTF("U11: Error code F77_CALL(dpotrf) = %d\n", err);
    Err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  }

 ErrorHandling: // and NOERROR...
  if (xx != NULL) free(xx);
  if (S!=NULL) {
      S->totpnts = totpnts;
      S->spatialpnts = spatialpnts;
      S->back = back;
      S->initial = initial;
      S->ntime = loc->length[spatialdim];
  }
  if (!storing && Err!=NOERROR) {
    if (MuT!=NULL) free(MuT);
    if (U11!=NULL) free(U11);
    if (U22!=NULL) free(U22);
    if (G!=NULL) free(G); 
    if (res0!=NULL) free(res0); 
  } else {
    if (S != NULL) {
      S->U22=U22;
      S->U11=U11;
      S->MuT=MuT;
      S->G=G;
      S->res0=res0;
     }
  }
  if (COV21!=NULL) {
      if (S != NULL && debugging)  S->Cov21 = COV21;
      else free(COV21);
  }
  if (Inv22!=NULL) {
      if (S != NULL && debugging)  S->Inv22 = Inv22;
      else free(Inv22);
  }

  return Err;
}


void sequentialpart(res_type *res, long totpnts, int spatialpnts, int ntime,
		    double *U11, double *MuT, double *G) {
  res_type *rp, *oldrp;
  int n, i, k, j, mutj;
  rp = &(res[totpnts]);
  oldrp = &(res[0]);
  for (n=0; n<ntime; n++, rp += spatialpnts, oldrp += spatialpnts) {
    for (i=0; i<spatialpnts; i++) G[i] = GAUSS_RANDOM(1.0);
    for (mutj=k=0, i=0; i<spatialpnts; i++, k+=spatialpnts) {
      res_type dummy;
      double *Uk;
      Uk = &U11[k]; 
      dummy =0.0;
      for (j=0; j<=i; j++) {
	dummy += G[j] * Uk[j];
      }
      for (j=0; j<totpnts; j++) {
	  dummy += MuT[mutj++] * (double) oldrp[j];
      }
      rp[i] = (res_type) dummy; 
    }
  }
}


void do_sequential(method_type *meth, res_type *res) 
{  
  sequential_storage *S;
  long totpnts, i, j, k;
  double *G,*U22, *U11, *MuT;
  res_type *res0; 

  S = (sequential_storage*) meth->S;
  totpnts = S->totpnts;
 
  U22 = S->U22;  // S22^{1/2}
  U11 = S->U11;  // S11^{1/2}
  MuT = S->MuT;
  res0 = S->res0;
  G = S->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)

  // Start
  for (i=0; i<totpnts; i++) G[i] = GAUSS_RANDOM(1.0);
  for (k=0, i=0; i<totpnts; i++, k+=totpnts){
    double dummy, *Uk;
    Uk = &U22[k]; 
    dummy =0.0;
    for (j=0; j<=i; j++){
      dummy += G[j] * Uk[j];
    }
    res0[i] = (res_type) dummy; 
  }
  
  sequentialpart(res0, totpnts, S->spatialpnts, S->initial, U11, MuT, G);
  res0 += S->initial * S->spatialpnts;
  memcpy(res, res0, sizeof(res_type) * totpnts);
  sequentialpart(res, totpnts, S->spatialpnts, S->ntime - S->back, 
		 U11, MuT, G);
}
