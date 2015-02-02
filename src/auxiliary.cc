//#define DEBUG 1

/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of auxiliary functions

 Copyright (C) 2001 -- 2014 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include <math.h>
#include <unistd.h>
 
#include "auxiliary.h"
//#include <curses.h>
#include "RandomFields.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "RF.h"
#include <Rdefines.h>
//#include <curses.h>


// important check !!
#ifndef SCHLATHERS_MACHINE
#ifdef SHOW_ADDRESSES
SHOW_ADRESSES IS NOT ALLOWED
#endif
#ifdef RANDOMFIELDS_DEBUGGING
RANDOMFIELDS_DEBUGGING IS NOT ALLOWED
#endif
#endif



double intpow(double x, int p) {
  double res = 1.0;
  if (p < 0) {
    p = -p;
    x = 1.0 / x;
  } 
  while (p != 0) {
    if (p % 2 == 1) res *= x;
    x *= x;
    p /= 2;
  }
  return res;
}

void AtA(double *a, int nrow, int ncol, double *A) {
  // A =  a^T %*% a
  int i, k, j, m,
    dimSq = ncol * ncol;
  
  for (k=i=0; i<dimSq; i+=ncol) {
    for (j=0; j<dimSq; j+=ncol, k++) {
      double dummy = 0.0;
      for (m=0; m<nrow; m++) {
	dummy += a[i+m] * a[j+m];
      }
      A[k] = dummy;
    }
  }
}
 

void xA(double *x, double*A, int nrow, int ncol, double *y) {
  int i,j,k;
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(double) * nrow);
  } else {
    for (k=i=0; i<ncol; i++) {
      y[i] = 0.0;
      for (j=0; j<nrow; j++) {
	y[i] += A[k++] * x[j];
      }
    }
  }	
}


void xA(double *x1, double *x2,  double*A, int nrow, int ncol, double *y1,
	double *y2) {
  int i,j,k;
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y1, x1, sizeof(double) * nrow);
    MEMCOPY(y2, x2, sizeof(double) * nrow);
  } else {
    for (k=i=0; i<ncol; i++) {
      y1[i] = y2[i] = 0.0;
      for (j=0; j<nrow; j++) {
	y1[i] += A[k] * x1[j];
	y2[i] += A[k++] * x2[j];
      }
    }
  }	
}


void Ax(double *A, double*x, int nrow, int ncol, double *y) {
  int i,j,k;
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(double) * nrow);
  } else {
    for (i=0; i<nrow; i++) y[i]=0.0;
    for (k=i=0; i<ncol; i++) { 
      for (j=0; j<nrow; j++) {
	y[j] += A[k++] * x[i];
      }
    }
  }
}


void Ax(double *A, double*x1, double*x2, int nrow, int ncol, double *y1,
	double *y2) {
  int i,j,k;
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y1, x1, sizeof(double) * nrow);
    MEMCOPY(y2, x2, sizeof(double) * nrow);
  } else {
    for (i=0; i<nrow; i++) y1[i]=y2[i]=0.0;
    for (k=i=0; i<ncol; i++) { 
      for (j=0; j<nrow; j++) {
	y1[j] += A[k] * x1[i];
	y2[j] += A[k++] * x2[i];
      }
    }
  }
}




void AxResType(double *A, res_type *x, int nrow, int ncol, double *y) {
  int i,j,k;
  for (i=0; i<nrow; i++) y[i]=0.0;
  for (k=i=0; i<ncol; i++) { 
    for (j=0; j<nrow; j++) {
      y[j] += A[k++] * x[i];
    }
  }
}


double XkCXtl(double *X, double *C, int nrow, int dim, int k, int l) {
  // (k-th row of X) * C * (l-th column of X)
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


void XCXt(double *X, double *C, double *V, int nrow, int dim /* dim of C */) {
    double *pX, *endpX, *dummy, *pdummy, scalar;
    int i, cd, rv, ci, cv,
      size = nrow * dim;

  dummy = (double*) MALLOC(sizeof(double) * nrow * dim);
  if (dummy == NULL) ERR("XCXt: memory allocation error in XCXt");
  

  // dummy = XC
  for (pX = X, pdummy=dummy, endpX = pX + nrow;
       pX < endpX; pX++, pdummy++) {
    for (ci=0, cd=0; cd<size; cd+=nrow) {
      for (i=0, scalar = 0.0; i<size; i+=nrow) {
        scalar += pX[i] * C[ci++];
      }
      pdummy[cd] = scalar;
    }
  }

  // V = dummy X^t
  for (rv=0; rv<nrow; rv++) {
    for (cv=rv; cv<nrow; cv++) {
      for (scalar=0.0, i=0; i<size; i+=nrow) {
	scalar += dummy[rv + i] * X[cv + i];
     }
      V[rv + cv * nrow] = V[cv + rv * nrow] = scalar;
    }
  } 

  UNCONDFREE(dummy);
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

double getMinimalAbsEigenValue(double *Aniso, int dim) {
  double dummy, 
    min = RF_INF, 
    *SICH = NULL, 
    *D = NULL, 
    *work = NULL;
  int dd, Err,
    err = NOERROR,
    *iwork = NULL,
     dimSq = dim * dim,
    optim_work = 12 * dim;

  if ((D =(double *) MALLOC(sizeof(double) * dim))==NULL ||
      (work = (double *) MALLOC(sizeof(double) * optim_work))==NULL ||
      (iwork = (int *) MALLOC(sizeof(int) * 8 * dim))==NULL ||
      (SICH =(double *) MALLOC(sizeof(double) * dimSq))==NULL) {
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  MEMCOPY(SICH, Aniso, sizeof(double) * dimSq);
  F77_CALL(dgesdd)("N", &dim, &dim, SICH, &dim, D, NULL, &dim, NULL,
		   &dim, work, &optim_work, iwork, &Err);
  if (Err != 0) GERR("SVD for anisotropy matrix failed.");
  for (dd = 0; dd < dim; dd++) {
    dummy = fabs(D[dd]);
    if (dummy < min) min = dummy;
  } 

 ErrorHandling:
  FREE(D);
  FREE(SICH);
  FREE(work);
  FREE(iwork);
  if (err != NOERROR) XERR(err);

  return min;
}


double getDet(double *Aniso, int dim) {
  double  
    det = 1.0, 
    *SICH = NULL, 
    *D = NULL, 
    *work = NULL;
  int dd, Err,
    err = NOERROR,
    *iwork = NULL,
     dimSq = dim * dim,
    optim_work = 12 * dim;

  if ((D =(double *) MALLOC(sizeof(double) * dim))==NULL ||
      (work = (double *) MALLOC(sizeof(double) * optim_work))==NULL ||
      (iwork = (int *) MALLOC(sizeof(int) * 8 * dim))==NULL ||
      (SICH =(double *) MALLOC(sizeof(double) * dimSq))==NULL) {
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  MEMCOPY(SICH, Aniso, sizeof(double) * dimSq);
  F77_CALL(dgesdd)("N", &dim, &dim, SICH, &dim, D, NULL, &dim, NULL,
		   &dim, work, &optim_work, iwork, &Err);
  if (Err != 0) GERR("SVD for anisotropy matrix failed.");
  for (dd = 0; dd < dim; det *= D[dd++]);

 ErrorHandling:
  FREE(D);
  FREE(SICH);
  FREE(work);
  FREE(iwork);
  if (err != NOERROR) XERR(err);

  return det;
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
    ERR("detU: matrix does not seem to be strictly positive definite");
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
  if (info != 0) {
    ERR("det_UpperInv: dpofa failed -- is matrix positive definite?");
  }

  double Det = 1.0;
  for (i=0; i<dimsq; i+=dimP1) Det *= C[i];
  *det = Det * Det;

  F77_CALL(dpodi)(C, &dim, &dim, det, &job); // C is now Cinv
}



void matmult(double *A, double *B, double *C, int l, int m, int n) {
// multiplying an lxm- and an mxn-matrix, saving krigult in C
   int i, j, k;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
	C[j*l+i] = 0;
	for(k=0; k<m; k++)
	   C[j*l+i] += A[k*l+i]*B[j*m+k];
     }
}

void matmulttransposed(double *A, double *B, double *C, int m, int l, int n) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
   int i, j, k;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
	C[j*l+i] = 0;
	for(k=0; k<m; k++)
	   C[j*l+i] += A[i*m+k]*B[j*m+k];
     }
}



void Xmatmult(double *A, double *B, double *C, int l, int m, int n) {
// multiplying an lxm- and an mxn-matrix, saving krigult in C
  int i, j, k, jl, jm, kl, endfor;
  double dummy;
  for(i=0; i<l; i++) {
    for(jl=i, jm=j=0; j<n; j++, jl+=l, jm+=m) {
       dummy = 0.0;
       endfor = jm + m;
       for(kl=i, k=jm; k<endfor; k++, kl+=l) dummy += A[kl] * B[k]; 
       C[jl] = dummy;
     }
  }
}

void Xmatmulttransposed(double *A, double *B, double *C, int m, int l, int n) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
  int i, j, k, jl, im, jm, endfor, jmk;
  double dummy;
  for(im=i=0; i<l; i++, im+=m) {
    for(jl=i, jm=j=0; j<n; j++, jl+=l, jm+=m) {
      dummy = 0.0;
      endfor = im + m;
      for(jmk=jm, k=im; k<endfor; k++) dummy += A[k] * B[jmk++]; 
      C[jl] = dummy;
    }
  }
}

/*
void InvChol(double *C, int dim) {
  int i, info, ii, endfor,
    job = 01,
    dimP1 = dim + 1,
    dimsq = dim * dim;
  long ve, ho;
  double Det = 1.0;
  F77_CALL(dpofa)(C, &dim, &dim, &info); // C is now cholesky
  if (info != 0) ERR("InvChol: Inversion failed, bad functions\n");
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
*/

int invertMatrix(double *M, int size, int* Methods, int nMeth) {
  // http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F01/f01intro.xml
  int err, m, i,
    method = -1,
    sizeSq = size * size;
  long j;
  double *SICH = NULL;

  if (nMeth > 1) {
    if ((SICH = (double*) MALLOC(sizeof(double) * sizeSq)) == NULL) {
      return ERRORMEMORYALLOCATION;
    }
    memcpy(SICH, M, sizeSq);
  }
  
  for (m=0; m<nMeth; m++) {
    method = Methods[m];
    if (method<0) break;
    if (m>0) {
      memcpy(M, SICH, sizeSq);
    }
    switch(method) {
    case Cholesky : // cholesky
      F77_CALL(dpotrf)("Upper", &size, M, &size, &err);  
       if (err != 0) {
	if (PL >= PL_COV_STRUCTURE) 
	  PRINTF("cholesky decomposition failed; the matrix does not seem to be strictly positive definite\n");
	continue;
      }      
      F77_CALL(dpotri)("U", &size, M, &size, &err);  
      if (err != 0) {
	if (PL >= PL_COV_STRUCTURE) 
	  PRINTF("Cholesky decomposition failed; matrix does not seem to be strictly positive definite\n");
	continue;
      }
      long  i2, i3;
      for (i2=i=0; i<size; i++, i2+=size + 1) {	
	for (i3 = i2 + 1, j = i2 + size; j<sizeSq; j+=size) M[i3++] = M[j];
      }
     break;

    case 1 : {// QR returns transposed of the inverse !!
      int 
	*inc = NULL;
      double *aijmax = NULL,
	*workspaceD = NULL,
	*workspaceU = NULL,
	*workspaceDU = NULL;
      if (//(aijmax = (double*) MALLOC(sizeof(double) * size)) == NULL ||
	  //(inc = (int*) MALLOC(sizeof(int) * size)) == NULL ||
	  (workspaceD = (double*) MALLOC(sizeof(double) * size)) == NULL ||
	  (workspaceU = (double*) MALLOC(sizeof(double) * sizeSq)) == NULL 
	  // || (workspaceDU = (double*) MALLOC(sizeof(double) * size)) == NULL
	  ) {
	err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }


      // PseudoInverse:
      //    F77_CALL(f01blf)(&size, &size, &(GLOBAL.general.matrixtolerance),
       //		       M, &size, aijmax, &irank, inc, workspaceD, 
      //		       workspaceU, &size, workspaceDU, &err);
      F77_CALL(dgeqrf)(&size, &size,
		       M, &size, // aijmax, &irank, inc, workspaceD, 
		       workspaceU, workspaceD,  &size, &err);
      if (err != NOERROR) GERR("QR decomposition failed");
     
     ErrorHandling:
      FREE(aijmax);
      FREE(inc);
      FREE(workspaceD);
      FREE(workspaceU);
      FREE(workspaceDU);
      if (err == NOERROR) break;
      
      break;
    }

    case 2 : {// SVD
      double  optim_lwork, 
	tol = GLOBAL.general.matrixtolerance,
	*VT = NULL,
	*work = NULL,
	*U = NULL,
	*D = NULL;
      int k,  
	lwork = -1,
	*iwork = NULL;
      
      if ((VT =(double *) MALLOC(sizeof(double) * sizeSq))==NULL ||
	  (U =(double *) MALLOC(sizeof(double) * sizeSq))==NULL ||
	  (D =(double *) MALLOC(sizeof(double) * size))==NULL ||
	  (iwork = (int *) MALLOC(sizeof(int) * 8 * size))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling2;
      }

      F77_CALL(dgesdd)("A", &size, &size, M, &size, D, M, &size, VT, &size, 
		       &optim_lwork, &lwork, iwork, &err);
      if (err != NOERROR) {
	err=ERRORDECOMPOSITION;
	goto ErrorHandling2;
      }

      lwork = (int) optim_lwork;
      if ((work = (double *) MALLOC(sizeof(double) * lwork))==NULL) {
	err=ERRORMEMORYALLOCATION;
	goto ErrorHandling2;
      }
      F77_CALL(dgesdd)("A",  &size,  &size, M, &size, D, U, &size, VT, &size, 
		       work, &lwork, iwork, &err);
      
      if (err==NOERROR && ISNAN(D[0])) err=9999;
      if (err!=NOERROR) {
	if (PL>PL_ERRORS) PRINTF("Error code F77_CALL(dgesdd) = %d\n", err); 
	err=ERRORDECOMPOSITION;
	goto ErrorHandling2;
      }
      
      /* calculate SQRT of covariance matrix */
      for (k=0, j=0; j<size; j++) {
	double dummy;
	dummy = fabs(D[j]) < tol ? 0.0 : 1.0 / D[j];
	for (i=0; i<size; i++) U[k++] *= dummy;
      }
      matmult(U, VT, M, size, size, size); // U * V
 
     ErrorHandling2:
      FREE(VT);
      FREE(U);
      FREE(D);
      FREE(iwork);
      FREE(work);
     
      break;
    }
    default : BUG;
      
    } // switch
    if (err==NOERROR) break;
  }
  
  FREE(SICH);
  if (method<0 || m>=nMeth) SERR("matrix inversion failed");
  return NOERROR; //  -method;
}
   
int invertMatrix(double *M, int size) {
  return invertMatrix(M, size, GLOBAL.general.matrix_inversion, MAXINVERSIONS);
}
 
//void solveMatrix(int method, double *M, double *v, double *res, int size) {
//  // alles symmetrische matrizen!!
//  // http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F01/f01intro.xml
//  int i,j, k // err,
//    //sizeSq = size * size
//    ;
//  //  double *U;
// 
//  for (i=k=0; i<size; i++) {
//    double dummy = 0.0;
//    for (j=0; j<size; j++) dummy += M[k++] * v[j];
//    res[i] = dummy;
//  }

  /*
  switch(method) {
  case 0 : // cholesky 
    double dummy;
    for (i=0; i<size; i++) res[i] = 0.0;
    for (i=k=0; i<size; i++) {
      dummy = v[i];
      for (j=0; j<size; j++) res[j] += M[k++] * dummy;
    }     
    break;
  case 1 : // QR
    double dummy;
    for (i=k=0; i++; i<size) {
      dummy = 0.0;
      for (j=0; j++; j<size) dummy += M[k++] * v[j];
      res[i] = dummy;
      }
    break;
  case 2 : // SVD
    double dummy;
    for (i=0; i<size; i++) res[i] = 0.0;
    for (i=k=0; i<size; i++) {
      dummy = v[i];
      for (j=0; j<size; j++) res[j] += M[k++] * dummy;
    } 
    break;
  default : BUG;
  }   
  */
//}


void memory_copy(void *dest, void *src, int bytes) {
  int i, 
    len = bytes / sizeof(int),
    *d = (int*) dest,
    *s = (int *) src;
  if ((len * (int) sizeof(int)) != bytes) {
    ERR("size not a multiple of int");
  }
  for (i=0; i<len; i++) d[i] = s[i];
}


SEXP distInt(SEXP XX, SEXP N, SEXP Genes) {
  int i,j, k, di, diff, *x, *y, ve, ho, endfor,
    *X = INTEGER(XX),
    n = INTEGER(N)[0],
    nP1 = n + 1,
    genes = INTEGER(Genes)[0];
 
  SEXP Dist;
  PROTECT(Dist = allocMatrix(REALSXP, n, n));
  double *dist = REAL(Dist);

  x = y = X;
  for (j=0, i=0;  j<n;  i += nP1, j++, y += genes) {
    dist[i] = 0.0;
    endfor = i + (n - j);
    for (ve = i + 1, ho = i + n, x = y + genes; 
         ve < endfor; 
	 ve++, ho += n) {
      for (di=0.0, k=0; k<genes; k++, x++) {
	diff = *x - y[k];
	di += diff * diff;
      }
      dist[ve] = dist[ho] = sqrt((double) di);
    }
  }
  UNPROTECT(1);
  return Dist;
}

/*

SEXP GetChar(SEXP N, SEXP Choice, SEXP Shorter, SEXP Beep, SEXP Show) {  
    int i, j, 
      start = NA_INTEGER,
      milli = 500,
      len = length(Choice), 
      n = INTEGER(N)[0],
      nline = 100;
  bool shorter = LOGICAL(Shorter)[0],
      piep = LOGICAL(Beep)[0],
      show = LOGICAL(Show)[0];
  char choice[256], *zk,
      backspace = 127,
      eof = 'e'; // crl-D, ^D
  SEXP Zk;
  
  zk = (char*) MALLOC(sizeof(char) * (n+1));
  if (len > 256) len = 256;  
  for (i=0; i<len; i++) {
      choice[i] = CHAR(STRING_ELT(Choice, i))[0];
  }
 
  system("/bin/stty cbreak -echo iuclc"); // or "stty raw" 
  for (j=0; j<n; ) {
   if (j % nline == 0) {
       start = j;
       if (show) {
	   if (j != 0) PRINTF("\n");
	   int m = n-j;
	   if (m > nline) m = nline;
	   for (i=0; i<m; i++) if (i % 20 == 19) PRINTF(":"); else PRINTF(".");
	   for (i=0; i<m; i++) PRINTF("\b");
       }
   }

    zk[j] = getchar();
    if (zk[j] == eof && shorter) break;
    if (zk[j] == backspace && j>start) {
      j--;
      PRINTF("\b \b");
      continue;
    }
    for (i=0; i<len; i++)
      if (zk[j] == choice[i]) break;
    if (i < len) {
	if (show) PRINTF("%c", zk[j], n-1-j);
        j++; 
    } else {
	if (piep) PRINTF("beep does not work.\n"); // beep();
    }
  }
//  beep();
  if (show) {
      PRINTF("\n");
      sleepMilli(&milli);
  }
  system("/bin/stty -cbreak echo -iuclc");
  zk[j] = '\0';
  PROTECT(Zk = allocVector(STRSXP, 1));
  SET_STRING_ELT(Zk, 0, mkChar(zk));
  UNPROTECT(1);      
  return Zk;
}

*/

void strcopyN(char *dest, const char *src, int n) {
  if (n > 1) {
    n--; 
    strncpy(dest, src, n);
  }
  dest[n] = '\0';
}

void I0ML0(double *x, int *n) {
  int i;
  for (i=0; i<*n; i++) x[i] = I0mL0(x[i]);
} 

double I0mL0(double x){
  /* Bessel I_0 - Struve L_0 for non-negative arguments x */
  /* see J. MacLeod, Chebyshev expansions for modified {S}truve and 
     related functions, Mathematics of Computation, 60, 735-747, 1993 */
  static double g2[24] = {
	0.52468736791485599138e0,
	-0.35612460699650586196e0,
	0.20487202864009927687e0,
	-0.10418640520402693629e0,
	0.4634211095548429228e-1,
	-0.1790587192403498630e-1,
	0.597968695481143177e-2,
	-0.171777547693565429e-2,
	0.42204654469171422e-3,
	-0.8796178522094125e-4,
	0.1535434234869223e-4,
	-0.219780769584743e-5,
	0.24820683936666e-6,
	-0.2032706035607e-7,
	0.90984198421e-9,
	0.2561793929e-10,
	-0.710609790e-11,
	0.32716960e-12,
	0.2300215e-13,
	-0.292109e-14,
	-0.3566e-16,
	0.1832e-16,
	-0.10e-18,
	-0.11e-18
  };
  static double g3[24] = {
	2.00326510241160643125e0,
	0.195206851576492081e-2,
	0.38239523569908328e-3,
	0.7534280817054436e-4,
	0.1495957655897078e-4,
	0.299940531210557e-5,
	0.60769604822459e-6,
	0.12399495544506e-6,
	0.2523262552649e-7,
	0.504634857332e-8,
	0.97913236230e-9,
	0.18389115241e-9,
	0.3376309278e-10,
	0.611179703e-11,
	0.108472972e-11,
	0.18861271e-12,
	0.3280345e-13,
	0.565647e-14,
	0.93300e-15,
	0.15881e-15,
	0.2791e-16,
	0.389e-17,
	0.70e-18,
	0.16e-18
  };
  double r, x2, ac;
  int i;
  
  if (x < 0.0) {return RF_NA;}
  if (x < 16.0) {
    r = 0.5 * g2[0];
    ac = acos((6.0 * x - 40.0) / (x + 40.0));
    for (i=1; i<24; i++) {
      r += g2[i] * cos(i * ac);
    }
  } else {
    r = 0.5 * g3[0];
    x2 = x * x;
    ac = acos((800.0 - x2) / (288.0 + x2));
    for (i=1; i<24; i++) {
      r += g3[i] * cos(i * ac);
    }
    r *= T_PI /* 2/pi */ / x;
  }
  return r;
}

void StruveH(double *x, double *nu) {*x=struve(*x, *nu, -1.0, false);}
void StruveL(double *x, double *nu, int * expScaled) 
{
  *x=struve(*x, *nu, 1.0, (bool) *expScaled);
}
double struve(double x, double nu, double factor_Sign, bool expscaled)
{ 
  double Sign, value, epsilon=1e-20;
  double dummy, logx, x1, x2;
  if ((x==0.0) && (nu>-1.0)) return 0.0;
  if (x<=0) return RF_NA; // not programmed yet
  logx = log(0.5 * x);
  x1=1.5;   
  x2=nu+1.5;   
  Sign=1.0;
  if (x2 > 0.0) { 
    dummy = (nu + 1.0) * logx - lgammafn(x1) - lgammafn(x2);
    if (expscaled) dummy -= x;
    value = exp(dummy);
  } else {
    if ( (double) ((int) (x1-0.5)) != x1-0.5 ) return RF_NA;
    value=pow(0.5 * x, nu + 1.0) / (gammafn(x1) * gammafn(x2));
    if (expscaled) value *= exp(-x);
    if ((dummy= value) <0) {
      dummy = -dummy;
      Sign = -1.0;
    }
    dummy = log(dummy);
  }
  logx *= 2.0;
  do {
    if (x2<0) { Sign = -Sign; }    
    dummy += logx - log(x1) - log(fabs(x2));
    value +=  Sign * exp(dummy);
    x1 += 1.0;
    x2 += 1.0;
    Sign = factor_Sign * Sign; 
  } while (exp(dummy) > fabs(value) * epsilon);
  return value;
}

SEXP vectordist(SEXP V, SEXP DIAG){
  bool notdiag = !LOGICAL(DIAG)[0];
  int 
    d, dr, 
    rows = nrows(V), 
    cols = ncols(V),
    rescols = (int) (cols * (cols - 1 + 2 * (int) !notdiag) / 2);
  double *v1, *v2, *end, *Dist,
    *v = REAL(V);
  
  end = v + cols * rows; 
  SEXP DIST;
  
  PROTECT(DIST = allocMatrix(REALSXP, rows, rescols));
  Dist = REAL(DIST);

  for (dr=0, v1=v; v1<end; v1+=rows) { // loop is one to large??
    v2 = v1;
    if (notdiag) {
       v2 += rows;
    }
    for (; v2<end; ) {
      for (d=0; d<rows; v2++) {
	Dist[dr++] = v1[d++] - *v2;
      }
    }
  }
  UNPROTECT(1);
  return DIST;
}

static int ORDERDIM;
static double *ORDERD;
static int *ORDERDINT;

typedef bool (*vergleich)(int, int);
vergleich SMALLER=NULL, GREATER=NULL;

bool smaller(int i, int j)
{
  double *x, *y;
  int d;
  x = ORDERD + i * ORDERDIM;
  y = ORDERD + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
     if (x[d] != y[d]) return x[d] < y[d];
  return false;
}

bool greater(int i, int j)
{
  double *x, *y;
  int d;
  x = ORDERD + i * ORDERDIM;
  y = ORDERD + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
    if (x[d] != y[d]) return x[d] > y[d];
  return false;
}

bool smallerInt(int i, int j)
{
  int *x, *y;
  int d;
  x = ORDERDINT + i * ORDERDIM;
  y = ORDERDINT + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++) {
    if (x[d] != y[d]) {
     return x[d] < y[d];
    }
  }
  return false;
}



bool greaterInt(int i, int j)
{
  int *x, *y;
  int d;
  x = ORDERDINT + i * ORDERDIM;
  y = ORDERDINT + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
    if (x[d] != y[d]) return x[d] > y[d];
  return false;
}

int is_positive_definite(double *C, int dim) {
    int err,
	bytes = sizeof(double) * dim * dim;
  double *test;
  test = (double*) MALLOC(bytes);
  MEMCOPY(test, C, bytes);
  F77_CALL(dpofa)(test, &dim, &dim, &err); 
  UNCONDFREE(test);
  return(err == 0);
}

void order(int *pos, int start, int end) {
  int randpos, pivot, left, right, pivotpos, swap;
  
  if( start < end ) {   
    //Get RNGstate();randpos = start + (int) (UNIFORM_RANDOM * (end-start+1)); PutRNGstate(); // use Get/Put RNGstate with great care !!
    randpos = (int) (0.5 * (start + end));
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    
    pivotpos=start; 
    left = start;
    right=end+1;   
    while (left < right) {      
      //printf("order > %ld start=%d %d left=%d %d %d pivot=%d\n", pos, start, end, left, right, pos[left], pivot);
      while (++left < right && SMALLER(pos[left], pivot)) {
	pivotpos++;
      }
      while (--right > left && GREATER(pos[right], pivot));      
      if (left < right) {
	swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
	pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    order(pos, start, pivotpos-1 );
    order(pos, pivotpos + 1, end );
  }
}

void ordering(double *d, int len, int dim, int *pos) 
{
  /* quicksort algorithm, slightly modified:
     does not sort the data, but d[pos] will be ordered 
     NOTE: pos must have the values 0,1,2,...,start-end !
     (orderdouble is a kind of sorting pos according to
     the variable d)
  */  
  int i;

  for (i=0; i<len;i++) pos[i]=i;
  ORDERD = d;
  ORDERDIM = dim;
  SMALLER = smaller;
  GREATER = greater;
  order(pos, 0, len-1);
}


void Ordering(double *d, int *len, int *dim, int *pos) 
{
  int i;
  for (i=0; i<*len; i++) pos[i]=i;
  ORDERD = d;
  ORDERDIM = *dim;
  SMALLER = smaller;
  GREATER = greater;
  order(pos, 0, *len-1 );
}


void orderingInt(int *d, int len, int dim, int *pos) 
{
  /* quicksort algorithm, slightly modified:
     does not sort the data, but d[pos] will be ordered 
     NOTE: pos must have the values 0,1,2,...,start-end !
     (orderdouble is a kind of sorting pos according to
     the variable d)
  */  
  
 
  int i;
  for (i=0; i<len;i++) pos[i]=i;
  ORDERDINT = d;
  ORDERDIM = dim;
  SMALLER = smallerInt;
  GREATER = greaterInt; 
  order(pos, 0, len-1);
}


int xMatch(char *name, char **list, unsigned int llen)  {
  unsigned int ln, nr;
  // == -1 if no matching function is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  // see also GetModelNr !

  nr=0;
  ln=(unsigned int) strlen(name);
  while ((nr < llen) && strncmp(name, list[nr], ln)) {
    nr++;
  }
  if (nr < llen) { 
    // a matching method is found. Are there other methods that match?
    unsigned int j; 
    if (ln == (unsigned int) strlen(list[nr])) return nr; // exact matching
    j = nr + 1; // if two or more methods have the same name the last one is
    //          taken; stupid here, but nice in GetCovNr
    while (j < llen && strncmp(name, list[j], ln)) j++;
    if (j < llen) {
      if (ln == (unsigned int) strlen(list[j])) return j; //exact matching 
      else return -2; // multiple matching
    }
    return nr;
  } else return -1; // unmatched
}


int CeilIndex(double x, double *cum, int size) {  
   // der kleinste index i so das cum[i] >= x --- sollte das gleiche sein
  // wie searchFirstGreater

  int mitte,
    min = 0,
    max = size - 1;
  while (min < max) {
    mitte = 0.5 * (min + max);
    if (cum[mitte] >= x) max = mitte;
    else min = mitte + 1;
  }
  //  printf("%f < %f <= %f\n", cum[min-1], x, cum[min]);
  assert((min==0) || x > cum[min-1]);
  assert(x <= cum[min] && (min==0 || x > cum[min-1]));


  return min;
} 

/*
int searchFirstGreater(double *v, int len, double z) {
  // assumes the list has non decreasing values
  
  int 
    firsti = 0,
    lasti = len - 1,
    i = len / 2;
  
  if (z > v[lasti]) {
    BUG;
  }
  if (z <= v[firsti]) return firsti;
  
  while (firsti < lasti) {
    if (v[i] < z) {
      firsti = i + 1;
      i = i + (int) ((lasti - i + 1) / 2);
    } else { // v[i] <= z
      lasti = i;
      i = i - (int) ((i - firsti + 1) / 2);
    }
  }
    
  assert((i==0) || z > v[i-1]);
  assert(z <= v[i]);

  //  printf("%f < %f <= %f\n", v[i-1], z, v[i]);//
}
*/


/*
double searchInverse(isofct fct, double start, double *value,
		     double releps) {	       
  while (fct(start, cov) > value) start *= 2.0;    
  while (fct(start, cov) < value) start *= 0.5;
 
  double x = start,
    step = start;
  releps *= start;
  while (step > releps) {
    step *= 0.5;
    if (fct(step, cov) < value) x -= step; else x += step;
  }
}
*/

double searchInverse(covfct fct, cov_model *cov, 
		     double start, double value, double releps) {
  double v;
  fct(&start, cov, &v);
  while (v > value) {start *= 2.0; fct(&start, cov, &v);}
  while (v < value) {start *= 0.5; fct(&start, cov, &v);}
 
  double x = start,
    step = start;
  releps *= step;
  while (step > releps) {
    step *= 0.5;
    fct(&step, cov, &v);
    if (v < value) x -= step; else x += step;
  }
  return x;
}

double searchInverse(covfct fct, cov_model *cov, 
		     double start, double min, double value, double releps) {
  double v;
  assert(start > min);
  fct(&start, cov, &v);
  while (v > value) {start = 2.0 * (start - min) + min; fct(&start, cov, &v);}
  while (v < value) {start = 0.5 * (start - min) + min; fct(&start, cov, &v);}
 
  double x = start,
    step = start - min;
  releps *= step;
  while (step > releps) {
    step *= 0.5;
    fct(&step, cov, &v);
    if (v < value) x -= step; else x += step;
  }
  return x;
}

//double gamma(double x) { BUG; } // use gammafn instead


double incomplete_gamma(double start, double end, double s) {
  // int_start^end t^{s-1} e^{-t} \D t

  // print("incomplete IN s=%f e=%f s=%f\n", start, end, s);

  double
    v = 0.0, 
    w = 0.0;

  if (s <= 1.0) {
    if (start == 0.0) return RF_NA;
  }
  
  double 
    e_start = exp(-start),
    e_end = exp(-end),
    power_start = pow(start, s),      
    power_end = end < RF_INF ? pow(end, s) : 0,
    factor = 1.0; 
  
  
  while (s < 0.0) {
    factor /= s;    
    v +=  factor * (power_end * e_end - power_start * e_start);
    power_start *= start;
    if (end < RF_INF) power_end *= end;
    s += 1.0;
  }
  
  w = pgamma(start, s, 1.0, false, false);  // q, shape, scale, lower, log
  if (R_FINITE(end)) w -= pgamma(end, s, 1.0, false, false);

  //  print("incomplete s=%f e=%f s=%f v=%f g=%f w=%f\n", start, end, s, v, gammafn(s), w);

  return v + gammafn(s) * w * factor;
}


int addressbits(void VARIABLE_IS_NOT_USED *addr) {
#ifndef RANDOMFIELDS_DEBUGGING  
  return 0;
#else
  double x = (intptr_t) addr,
    cut = 1e9;
  x = x - trunc(x / cut) * cut;
  return (int) x;
#endif

}

int Match(char *name, name_type List, int n) {
  // == -1 if no matching name is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  unsigned int ln;
  int Nr;
  Nr=0;
  ln=strlen(name);
  //  print("Match %d %d %s %s %d\n", Nr, n, name, List[Nr], ln);

  while ( Nr < n  && strncmp(name, List[Nr], ln)) {
    Nr++;
  }
  if (Nr < n) { 
    if (ln==strlen(List[Nr])) // exactmatching -- take first -- changed 1/7/07
      return Nr;
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && strncmp(name, List[j], ln)) {j++;}
      if (j<n) {
	if (ln==strlen(List[j])) { // exactmatching -- take first 
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {return MULTIPLEMATCHING;}
  } else return NOMATCHING;
  return Nr;
}

int Match(char *name, const char * List[], int n) {
  // == -1 if no matching name is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  unsigned int ln;
  int Nr;
  Nr=0;
  ln=strlen(name);
  //    print("Matchx %d %d %s %s %d\n", Nr, n, name, List[Nr], ln);

  while ( Nr < n  && strncmp(name, List[Nr], ln)) {
    // print("       %d %d %s %s %d\n", Nr, n, name, List[Nr], ln);
    //   printf("%s\n", List[Nr]);
    Nr++;
  }
  if (Nr < n) { 
    if (ln==strlen(List[Nr])) {// exactmatching -- take first -- changed 1/7/07
      //      print(" found  X    %d %d %s %s %d\n", Nr, n, name, List[Nr], ln);
      return Nr;
    }
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && strncmp(name, List[j], ln)) {j++;}
      if (j<n) {
	if (ln==strlen(List[j])) { // exactmatching -- take first 
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {return MULTIPLEMATCHING;}
  } else return NOMATCHING;

  //  print(" found      %d %d %s %s %d\n", Nr, n, name, List[Nr], ln);
 
  return Nr;
}


//static int ZZ = 0;
double Real(SEXP p, char *name, int idx) {
  char msg[200];
  // if(++ZZ==65724){printf("type=%d %d '%s'\n",ZZ,TYPEOF(p), CHAR(STRING_ELT(p,0)));cov_model *cov;crash(cov);}
  if (p != R_NilValue)
    switch (TYPEOF(p)) {
    case REALSXP :  return REAL(p)[idx];
    case INTSXP : return INTEGER(p)[idx]==NA_INTEGER  
	? RF_NA : (double) INTEGER(p)[idx];
    case LGLSXP : return LOGICAL(p)[idx]==NA_LOGICAL ? RF_NA 
	: (double) LOGICAL(p)[idx];
    default : {}
    }
  // MEMCOPY(msg, p, 300); print("%s\n", msg);
  sprintf(msg, "'%s' cannot be transformed to double! (type=%d)\n",
	  name, TYPEOF(p));  
  //printf("\n>>>> '%s'\n", CHAR(STRING_ELT(p, 0)));
  ERR(msg);
  return RF_NA;  // to avoid warning from compiler
}



void Real(SEXP el,  char *name, double *vec, int maxn) {
  char msg[200];
   int i, j, n;
  if (el == R_NilValue) {
    sprintf(msg,"'%s' cannot be transformed to double.\n", name);
    ERR(msg);
  }
  n = length(el);
  for (j=i=0; i<maxn; i++) {
    vec[i] = Real(el, name, j);
    if (++j >= n) j=0;
  }
  return;
}

int Integer(SEXP p, char *name, int idx, bool nulltoNA) {
  char msg[200];
  if (p != R_NilValue) {
    switch(TYPEOF(p)) {
    case INTSXP : 
      return INTEGER(p)[idx]; 
    case REALSXP : 
      double value;
      value = REAL(p)[idx];      
      if (ISNAN(value)) {
	return NA_INTEGER;
	//sprintf(msg, "%s: NAs not allowed for integer valued parameters", name);
	//	ERR(msg);
      }
      if (value == trunc(value)) return (int) value; 
      else {
	sprintf(msg, "%s: integer value expected", name);
	ERR(msg);
      }
    case LGLSXP :
      return  LOGICAL(p)[idx]==NA_LOGICAL ? NA_INTEGER : (int) LOGICAL(p)[idx];
    default : {}
    }
  } else if (nulltoNA) return NA_INTEGER;
  sprintf(msg, "%s: unmatched type of parameter [type=%d]", name, TYPEOF(p));
  ERR(msg);
  return NA_INTEGER; // compiler warning vermeiden
}

int Integer(SEXP p, char *name, int idx) {
  return Integer(p, name, idx, false);
}


void Integer(SEXP el, char *name, int *vec, int maxn) {
  char msg[200];
  int i, j, n;
  if (el == R_NilValue) {
    sprintf(msg, "'%s' cannot be transformed to integer.\n",name);
    ERR(msg);
  }
  n = length(el);
  for (j=i=0; i<maxn; i++) {
    vec[i] = Integer(el, name, j);
    if (++j >= n) j=0;
  }
}




void Integer2(SEXP el, char *name, int *vec) {
  char msg[200];
  int n;
  if (el == R_NilValue || (n = length(el))==0) {
    sprintf(msg, "'%s' cannot be transformed to integer.\n",name);
    ERR(msg);
  }
 
  vec[0] = Integer(el, name, 0);
  if (n==1) vec[1] = vec[0];
  else {
    vec[1] = Integer(el, name, n-1);
    if (n > 2) {
      int i, 
	v = vec[0] + 1;
      for (i = 1; i<n; i++, v++)
	if (Integer(el, name, i) != v) ERR("not a sequence of numbers"); 
    }
  }
}





bool Logical(SEXP p, char *name, int idx) {
  char msg[200];
  if (p != R_NilValue)
    switch (TYPEOF(p)) {
    case REALSXP: return ISNAN(REAL(p)[idx]) ? NA_LOGICAL : (bool) REAL(p)[idx];
    case INTSXP :
      return INTEGER(p)[idx]==NA_INTEGER ? NA_LOGICAL : (bool) INTEGER(p)[idx];
    case LGLSXP : return LOGICAL(p)[idx];
    default : {}
    }
  sprintf(msg, "'%s' cannot be transformed to logical.\n", name);  
  ERR(msg);
  return NA_LOGICAL;  // to avoid warning from compiler
}


char Char(SEXP el, char *name) {
  char msg[200];
  SEXPTYPE type;
  if (el == R_NilValue) goto ErrorHandling;
  type = TYPEOF(el);
  if (type == CHARSXP) return CHAR(el)[0];
  if (type == STRSXP) {
    if (length(el)==1) {
      if (strlen(CHAR(STRING_ELT(el,0))) == 1)
	return (CHAR(STRING_ELT(el,0)))[0];
      else if (strlen(CHAR(STRING_ELT(el,0))) == 0)
	return '\0';
    }
  }
 
 ErrorHandling:
  sprintf(msg, "'%s' cannot be transformed to character.\n",  name);  
  ERR(msg);
  return 0; // to avoid warning from compiler
}


void String(SEXP el, char *name, char names[MAXUNITS][MAXCHAR]) {
  int i,
    l = length(el);
  char msg[200];
  SEXPTYPE type;  
  if (el == R_NilValue) goto ErrorHandling;
  if (l > MAXUNITS)  {
    ERR1("number of variable names exceeds %d. Take abbreviations?",
	 MAXUNITS);
  }
  type = TYPEOF(el);
  //  printf("type=%d %d %d %d\n", TYPEOF(el), INTSXP, REALSXP, LGLSXP);
  if (type == CHARSXP) {
    for (i=0; i<l; i++) {
      names[i][0] = CHAR(el)[i];
      names[i][1] = '\0';
    }
  } else if (type == STRSXP) {
    for (i=0; i<l; i++) {
      //print("%d %d\n", i, l);
      strcopyN(names[i], CHAR(STRING_ELT(el, i)), MAXCHAR);
    }
  } else goto ErrorHandling;
  return;
 
 ErrorHandling:
  sprintf(msg, "'%s' cannot be transformed to character.\n",  name);  
  ERR(msg);
}



double NonNegInteger(SEXP el, char *name) {
  int num;

  num = INT;
  if (num<0) {
    num=0; 
    char msg[200];
    sprintf(msg,"'%s' which has been negative is set 0.\n",name);
    warning(msg);
  }
  return num;
}

double NonNegReal(SEXP el, char *name) {
  double num;
  num = NUM;
  if (num<0.0) {
    num=0.0; 
    char msg[200];
    sprintf(msg,"%s which has been negative is set 0.\n",name);
    warning(msg);
   }
  return num;
}

double NonPosReal(SEXP el, char *name) {
  double num;
  num = NUM;
  if (num>0.0) {
    num=0.0; 
    char msg[200];
    sprintf(msg,"%s which has been positive is set 0.\n",name);
    warning(msg);
  }
  return num;
}

double PositiveInteger(SEXP el, char *name) {
  int num;
  num = INT;
  if (num<=0) {
    num=0; 
    char msg[200];
    sprintf(msg,"'%s' which has been negative is set 0.\n",name);
    warning(msg);
  }
  return num;
}

double PositiveReal(SEXP el, char *name) {
  double num;
  num = NUM;
  if (num<=0.0) {
    num=0.0; 
    char msg[200];
    sprintf(msg,"%s which has been negative is set 0.\n",name);
    warning(msg);
   }
  return num;
}

int GetName(SEXP el, char *name, const char * List[], int n,
	    int defaultvalue) {
  char msg[1000], dummy[1000];
  int i,
    nM1 = n - 1;

  if (TYPEOF(el) == NILSXP) goto ErrorHandling;
 

  if (TYPEOF(el) == STRSXP) {
    int m = Match((char*) CHAR(STRING_ELT(el, 0)), List, n);

    if (m >= 0) return m; else {
      if (strcmp((char*) CHAR(STRING_ELT(el, 0)), " ") == 0 ||
	  strcmp((char*) CHAR(STRING_ELT(el, 0)), "") == 0) {
	goto ErrorHandling;
      }
    }
  }
  sprintf(dummy, "'%s': unknown value '%s'. Possible values are:", 
	  name, CHAR(STRING_ELT(el, 0)));
  for (i=0; i<nM1; i++) {
    sprintf(msg, "%s '%s',", dummy, List[i]);    
    strcpy(dummy, msg);
  }
  sprintf(msg,"%s '%s'.", dummy, List[i]);  
  ERR(msg);
 
 ErrorHandling:
  if (defaultvalue >= 0) return defaultvalue;
  
  sprintf(msg, "'%s': no value given.", name);
  ERR(msg);

  return 999;// to avoid warning from compiler
}

int GetName(SEXP el, char *name, const char * List[], int n) {
 return GetName(el, name, List, n, -1);
}
