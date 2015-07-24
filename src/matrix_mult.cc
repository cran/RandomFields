/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- Martin Schlather, Reinhard Furrer

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#include <R_ext/Lapack.h>
#include "RF.h"

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



void matmult(double *A, double *B, double *C, int l, int m, int n) {
// multiplying an lxm- and an mxn-matrix, saving lxn result in C
   int i, j, k;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
       double dummy = 0.0;
	for(k=0; k<m; k++) dummy += A[k*l+i]*B[j*m+k];
	C[j*l+i] = dummy;
     }
}

void matmulttransposed(double *A, double *B, double *C, int m, int l, int n) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
   int i, j, k;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
       double dummy = 0.0;
	for(k=0; k<m; k++) dummy += A[i*m+k]*B[j*m+k];
	C[j*l+i] = dummy;
     }
}



void matmult_2ndtransp(double *A, double *B, double *C, int m, int l, int n) {
// multiplying A and t(B) with dim(A)=(m,l) and dim(B)=(n, l),
// saving result in C
  int i, j, k,
    msq = m  * m;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
       double dummy = 0.0;
	for(k=0; k<msq; k+=m) dummy += A[i + k] * B[j + k];
	C[j*l + i] = dummy;
     }
}



void matmult_tt(double *A, double *B, double *C, int m, int l, int n) {
// calculating t(A B) with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
   int i, j, k;
   for(i=0; i<l; i++)
     for(j=0; j<n; j++) {
       double dummy = 0.0;
	for(k=0; k<m; k++) dummy += A[i + k * l] * B[j * m + k];
	C[i * l + j] = dummy;
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

