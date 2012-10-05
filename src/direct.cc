/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by Cholesky or SVD decomposition

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
 
#include "RF.h"
#include <R_ext/Lapack.h>
//#include <R_ext/Linpack.h>

bool debug=false;

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

void SetParamDirect(int *action,int *method, int *bestvariables,
			 int *maxvariables, double *svdtolerance)
{
  direct_param *gp = &(GLOBAL.direct);
  if (*action) {
    if ((*method<0) || (*method>=(int) NoFurtherInversionMethod)){
      PRINTF("inversion method out of range; ignored\n");
    }
    else {gp->inversionmethod= (InversionMethod) *method;}
    gp->bestvariables = *bestvariables;
    gp->maxvariables  = *maxvariables;
    gp->svdtolerance  = *svdtolerance;
  } else {
    *method = (int) gp->inversionmethod;
    *bestvariables = gp->bestvariables;
    *maxvariables  = gp->maxvariables;
    *svdtolerance  = gp->svdtolerance;
  }
}  


int init_directGauss(method_type *meth) {
  cov_model *cov = meth->cov;
  long err=NOERROR;
  double *xx, 
    *G=NULL, 
    *COV=NULL, 
    *U=NULL, 
    *VT=NULL, 
    *work=NULL,
    *D=NULL, 
    *SICH=NULL;
  int *iwork=NULL, 
    dim=cov->tsdim;
  direct_storage* S=NULL;
  InversionMethod method;
  location_type *loc = meth->loc;
  globalparam *gp = meth->gp;
  direct_param *lp = &(gp->direct);
  bool storing = gp->general.storing;
  nonstat_covfct cf;
  long 
    vdim = cov->vdim,
    locpts = loc->totalpoints,
//    loctot = locpts *dim,
    vdimtot = vdim * locpts,
    vdimSqtot = vdim * vdimtot,
    vdimtotSq = vdimtot * locpts,
    vdimSqtotSq = vdimtot * vdimtot;

  method = lp->inversionmethod;
  assert(meth->S == NULL); // ist in simu.cc sichergestellt, schadet aber nich
  SET_DESTRUCT(direct_destruct);
  if (vdimtot>lp->maxvariables) {
      sprintf(ERRORSTRING_OK, 
	    "number of points less than RFparameters()$direct.maxvariables (%d)",
	     lp->maxvariables);
      sprintf(ERRORSTRING_WRONG,"%ld", vdimtot);
      err=ERRORCOVFAILED; goto ErrorHandling;
  }
  if (meth->caniso != NULL || meth->cproj!=NULL || meth->cscale != 1.0) {
      err = ERRORPREVDOLLAR;
      goto ErrorHandling;
  }
   
  err=ERRORMEMORYALLOCATION;  
  if ((COV =(double *) malloc(sizeof(double) * vdimSqtotSq))==NULL)
    goto ErrorHandling;
  if ((U =(double *) malloc(sizeof(double) * vdimSqtotSq))==NULL)
      goto ErrorHandling;
  //for SVD/Chol intermediate results AND  memory space for do_directGauss:
  if ((G = (double *)  calloc(vdimtot + 1, sizeof(double)))==NULL)
      goto ErrorHandling;
  if ((meth->S = malloc(sizeof(direct_storage))) == NULL)
    goto ErrorHandling;
  S = (direct_storage*) meth->S;
  S->U = S->G = NULL;
  err = NOERROR;

  
  if (cov->statIn != STATIONARY && cov->statIn != COVARIANCE) {
    err = ERRORNOVARIOGRAM; 
    goto ErrorHandling;
  }

//  printf("%d %d %d\n", cov->statIn, STATIONARY, COVARIANCE);
//  PrintModelInfo(cov); assert(false);

  /* ************************* */
  /* create matrix of explicit */
  /*       x-coordinates       */
  /* ************************* */  

  if (loc->grid) {
     expandgrid(loc->xgr, loc->length, 
		&(meth->sptime), 
		loc->timespacedim);
     xx = meth->sptime;
     if (!loc->Time) meth->space = meth->sptime;
  } else {
    if (loc->Time) {
      xtime2x(loc->x, loc->length[0], loc->T, loc->length[dim-1],
	      &(meth->sptime), loc->timespacedim);
      //     printf("direct %d %d %d \n", meth->sptime, meth->space , loc->length[0]); assert(false);
      xx = meth->sptime;

//      for (i=0; i<loc->totalpoints * 3; i++)
//	 printf("%f\n", xx[i]);


    } else xx = loc->x;
  }
     
///  printf("vdimtot=%d %d %d %d %d\n", vdimtot, dim, vdim, locpts, loc->totalpoints);  assert(false);
//  PrintModelInfo(cov);

  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  long k, segment;
  int row, Err;
  double *x, *y;

  k = 0;
  cf = CovList[cov->nr].nonstat_cov;
  x = xx;
  LIST_ELEMENT = 0;
  if (cov->vdim == 1) {
    for (CovMatrixCol=0; CovMatrixCol<locpts; CovMatrixCol++, x+=dim) {
      k += CovMatrixCol;
      y = x;
      segment=CovMatrixCol* vdimtot+CovMatrixCol;
      for (CovMatrixRow=CovMatrixCol;  CovMatrixRow<locpts; 
	   CovMatrixRow++, y+=dim, segment += vdimtot) {
	cf(x, y, cov, COV + segment);
	
//	       PRINTF("%d (%d %d) %f %f %d %f\n", i, k, segment,
//	       x[0], y[0],  dim, COV[segment]);
	
//	     PRINTF("%d (%d %d)   [%f %f ; %f %f]   %d %f\n", i, k, segment,
//	     x[0], y[0], x[1], y[1], dim, COV[segment]);
	
	if (ISNA(COV[segment])) {
	  err=ERRORFAILED; goto ErrorHandling;
	}
	COV[k++] = COV[segment];
      }
    }
  } else {
    int l,n,m;
    double *z, *C;

    z = (double*) malloc(sizeof(double) * vdim * vdim);

    /*
     for (i=0; i<locpts; i++, x+=dim) {
      for (y=x, j=i; j<locpts; j++, y+=dim) {
	cf(x, y, cov, z);
	
	l = k = 0;
	C = COV + vdim * (i + j * vdimtot);
	for (n=0; n<vdim; n++, k += vdimtot - vdim) {
	  for (m=0; m<vdim; m++) {
	    C[k++] = z[l++];
	  }
	}

	if (i != j) {
	  l = k = 0;
	  C = COV + vdim * (j + i * vdimtot);
	  for (n=0; n<vdim; n++) {
	    k = n;
	    for (m=0; m<vdim; m++, k+=vdimtot) {
	      C[k] = z[l++];
	    }
	  }
	}
      }
    }
    */

    //   printf("method %d %d pts=%d\n", method, TCholesky, locpts);
 
    if (method == TCholesky || method == TSVD) {
      for (CovMatrixCol=0; CovMatrixCol<locpts; CovMatrixCol++, x+=dim) {
	for (y=x, CovMatrixRow=CovMatrixCol; CovMatrixRow<locpts; 
	     CovMatrixRow++, y+=dim) {
	  cf(x, y, cov, z);
	  
	  C = COV + vdim * (CovMatrixCol + CovMatrixRow * vdimtot);
	  for (l=n=0; n<vdimSqtot; n+=vdimtot) {
	    int endfor = n + vdim;
	    for (m=n; m<endfor; m++) {
	      C[m] = z[l++];
	    }
	  }
	  
	  if (CovMatrixCol != CovMatrixRow) {
	    C = COV + vdim * (CovMatrixRow + CovMatrixCol * vdimtot);
	    for (l=n=0; n<vdim; n++) {
	      for (m=n; m<vdimSqtot; m+=vdimtot) {
		C[m] = z[l++];
	      }
	    }
	  }
	}
      }
    } else {
      for (CovMatrixCol=0; CovMatrixCol<locpts; CovMatrixCol++, x+=dim) {
	for (y=x, CovMatrixRow=CovMatrixCol; CovMatrixRow<locpts; 
	     CovMatrixRow++, y+=dim) {
	  cf(x, y, cov, z);
	  
	  C = COV + CovMatrixCol + CovMatrixRow * vdimtot;
	  for (l=n=0; n<vdimSqtotSq; n+=vdimtotSq) {
	    int endfor=n + vdimtot;
	    for (m=n; m<endfor; m+=locpts) {
	      C[m] = z[l++];
	    }
	  }
	  
	  if (CovMatrixCol != CovMatrixRow) {
	    C = COV + CovMatrixRow + CovMatrixCol * vdimtot;
	    for (l=m=0; m<vdimtot; m+=locpts) {
	      for (n=m; n<vdimSqtotSq; n+=vdimtotSq) {
		C[n] = z[l++];
	      }
	    }
	  }
	}
      }
    }
    free(z);
  }

//  assert(false);


/*
  for (k=i=0; i<vdimtot; i++){
      for (j=0; j<vdimtot; j++) {
	PRINTF("%3.2f ", COV[k++]);
      }
      PRINTF("\n");
    }

//  assert(false);
*/


  /* ********************** */
  /*  square root of matrix */
  /* ********************** */
  switch (method) {
      case Cholesky : case TCholesky : 
	// only works for strictly positive def. matrices
	if (PL>=3) PRINTF("method to the root=Cholesky\n");
	row=vdimtot;
	// dchdc destroys the input matrix; upper half of U contains result!
	memcpy(U, COV, sizeof(double) * vdimSqtotSq);
	if (debug) {
	    err = ERRORDECOMPOSITION; goto ErrorHandling;
	}   
	F77_CALL(dpotrf)("Upper", &row, U, &row, &Err);


	// F77_NAME(dchdc)(U, &row, &row, G, NULL, &choljob, &err);
	if (Err!=NOERROR) {
	  if (PL>2)
	      PRINTF("Error code F77_CALL(dpotrf) = %d\n", Err);
	  err=ERRORDECOMPOSITION;
	} else break;
	if (lp->svdtolerance <= 0.0) break;
	// try next method : 
	// most common error: singular matrix 
	
      case SVD : case TSVD: // works for any positive semi-definite matrix
	double sum;
	method = SVD; // necessary if the value of method has been Cholesky.
	//               originally
	if (vdimtot>lp->maxvariables * 0.8) {
	  sprintf(ERRORSTRING_OK, 
		  "number of points less than 0.8 * RFparameters()$direct.maxvariables (%d) for SVD",
		  lp->maxvariables);
	  sprintf(ERRORSTRING_WRONG,"%ld", vdimtot);
	  err=ERRORCOVFAILED; goto ErrorHandling;
	}
	if (PL>=3) PRINTF("method to the root=SVD\n");
	err=ERRORMEMORYALLOCATION;
	if ((VT =(double *) malloc(sizeof(double) * vdimSqtotSq))==NULL)
	    goto ErrorHandling;
	if ((D =(double *) malloc(sizeof(double) * vdimtot))==NULL)
	    goto ErrorHandling;
	if ((iwork = (int *) malloc(sizeof(int) * 8 * vdimtot))==NULL)
	    goto ErrorHandling;
	if ((SICH =(double *) malloc(sizeof(double) * vdimSqtotSq))==NULL)
	    goto ErrorHandling;
	memcpy(SICH, COV, sizeof(double) * vdimSqtotSq);
	row=vdimtot;
        // dsvdc destroys the input matrix !!!!!!!!!!!!!!!!!!!!

	// DGESDD (or DGESVD)
      	// dgesdd destroys the input matrix COV;
	// F77_NAME(dsvdc)(COV, &row, &row, &row, D, e, U, &row, V, &row, G,
	//		&jobint /* 11 */ , &err);
	double optim_lwork;
	int lwork;
	lwork = -1;
	F77_CALL(dgesdd)("A", &row, &row, SICH, &row, D, U, &row, VT, &row, 
			 &optim_lwork, &lwork, iwork, &Err);
	if ((err=Err) != NOERROR) {
	  err=ERRORDECOMPOSITION;
	  goto ErrorHandling;
	}
	lwork = (int) optim_lwork;
	if ((work = (double *) malloc(sizeof(double) * lwork))==NULL)
	    goto ErrorHandling;
	F77_CALL(dgesdd)("A",  &row,  &row, SICH, &row, D, U, &row, VT, &row, 
			 work, &lwork, iwork, &Err);
	
	if (err==NOERROR && RF_ISNA(D[0]))
	    err=9999;
	if (err!=NOERROR) {
	  if (PL>2) 
	    PRINTF("Error code F77_CALL(dgesdd) = %d\n", err); 
	  err=ERRORDECOMPOSITION;
	  goto ErrorHandling;
	}
		
	int i,j;
	/* calculate SQRT of covariance matrix */
	for (k=0,j=0;j<vdimtot;j++) {
	  double dummy;
	  dummy = sqrt(D[j]);
	  for (i=0;i<vdimtot;i++) {
	    U[k++] *= dummy;
	  }
	}

	/* check SVD */
	if (lp->svdtolerance >=0) {
	  for (i=0; i<vdimtot; i++) {
	    for (k=i; k<vdimtot; k++) {
	      sum = 0.0;
	      for (j=0; j<vdimSqtotSq; j+=vdimtot) sum += U[i+j] * U[k+j];

	      if (fabs(COV[i * vdimtot + k] - sum) > lp->svdtolerance) {
	        if (PL > 3)
		    PRINTF("difference %e at (%d,%d) between the value (%e) of the covariance matrix and the square of its root (%e).\n", 
			   COV[i * vdimtot +k] - sum, i, k, COV[i* vdimtot +k ], 
			   sum);
		err = ERRORPRECISION;
		goto ErrorHandling;
	      }
	    }
	  }
	}
	break;
      default : assert(false);
  } // switch


 ErrorHandling: // and NOERROR...
  if (S != NULL) S->method = method;
  if (!storing && err!=NOERROR) {
    if (U!=NULL) free(U);
    if (G!=NULL) free(G); 
  } else {
    if (S != NULL) {
      S->U=U;
      S->G=G;
    }
  }
  if (SICH!=NULL) free(SICH);
  if (COV!=NULL) free(COV);
  if (D!=NULL) free(D);
  if (work!=NULL) free(work);
  if (iwork!=NULL) free(iwork);
  if (VT!=NULL) free(VT);
  LIST_ELEMENT = CovMatrixRow = CovMatrixCol = -1;

  return err;
}

void do_directGauss(method_type *meth, res_type *res) 
{  
  location_type *loc = meth->loc;
  direct_storage *S;
  long locpts, vdimtot, i,j,k, vdim;
  double *G,*U, dummy;  
  int m, n;

  S = (direct_storage*) meth->S;
  locpts = loc->totalpoints;
  vdim = meth->cov->vdim;
  vdimtot = locpts * vdim;

//  printf("vdim %d %d\n", meth->cov->vdim, loc->totalpoints);

  U = S->U;// S^{1/2}
  G = S->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)
  for (i=0; i<vdimtot; i++) G[i] = GAUSS_RANDOM(1.0);

  switch (S->method) {
      case TCholesky :
	for (k=0, i=0; i<vdimtot; i++, k+=vdimtot) {
	  double *Uk = U + k; 
	  dummy =0.0;
	  for (j=0; j<=i; j++){
	    dummy += G[j] * Uk[j];
	  }
	  res[i] += (res_type) dummy; 
	}
	break;
      case Cholesky :
	for (k=i=m=0; m<vdim; m++) {
	  for (n=m; n<vdimtot; n+=vdim, i++, k+=vdimtot) {
	    double *Uk = U + k;
	    dummy = 0.0;
	    for (j=0; j<=i; j++){
	      dummy += G[j] * Uk[j];
	    }
	    res[n] += (res_type) dummy; 
	  }
	}
	break;
      case TSVD :
	for (i=0; i<vdimtot; i++){
	  dummy = 0.0;
	  for (j=0, k=i; j<vdimtot; j++, k+=vdimtot){
	    dummy += U[k] * G[j];
	  }
	  res[i] += (res_type) dummy; 
	}
	break;
      case SVD :
	for (i=m=0; m<vdim; m++) {
	  for (n=m; n<vdimtot; n+=vdim, i++) {
	    dummy = 0.0;
	    for (j=0, k=i; j<vdimtot; j++, k+=vdimtot){
	      dummy += U[k] * G[j];
	    }
	    res[n] += (res_type) dummy; 
	  }
	}
	break;
      default : assert(false);
  }
}

