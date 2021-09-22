/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2001 -- 2017 Martin Schlather, 

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


#include "def.h"
#include "intrinsics.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "xport_import.h"


#if defined SSE2 
int inline DoubleAlign(int X) {
  return (2L + (X * sizeof(double)) / BytesPerBlock) *
    (BytesPerBlock / sizeof(double));
}
#endif


SEXP EstimBRsimuKoeff(SEXP M, SEXP Chol, SEXP N) {
#if defined SSE2
#else
  RFERROR("code only runs under SSE2"); // rest nicht getestet
#endif

  
#define atonce 6
  int
    c = ncols(M),
    c2 = c * c,
     n = INTEGER(N)[0] / atonce;
  if (c != ncols(Chol) || c != nrows(M) || c != nrows(Chol))
    RFERROR("matrices do not match");
  if (c == 0) return R_NilValue;
  if (n * atonce != INTEGER(N)[0])     
    RFERROR1("number of simulations must be divisible by %d", atonce);
  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, c, c));
  double *ans0 = REAL(Ans);
  for (int i=0; i<c2; ans0[i++] = 0.0);

#if defined SSE2 
  int aligned_c = DoubleAlign(c);
  assert( (aligned_c * sizeof(double)) % BytesPerBlock == 0);
#else
  int aligned_c = c;
#endif
  
  double *chol = REAL(Chol),
    *x = (double*) CALLOC(atonce * aligned_c, sizeof(double)),
    *xneu = (double *) CALLOC(atonce * aligned_c, sizeof(double));
  
#if defined SSE2 
    double *X = algn(x),
    *Xneu = algn(xneu);
#else
  double *X = x,
    *Xneu = xneu;
#endif

  if (GLOBAL_UTILS->solve.actual_pivot != PIVOT_NONE &&
     GLOBAL_UTILS->solve.actual_pivot != PIVOT_UNDEFINED)
   RFERROR("pivoting not programmed yet in EstimBRsimuKoeff");

  GetRNGstate();
   
  for (int k=0; k < n; k++) {
    for (int j=0; j<atonce; j++) {
      double *Xneuj = Xneu + j * aligned_c;
      for (int i=0; i<c; i++) Xneuj[i] = GAUSS_RANDOM(1.0);
    }
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) 
#endif
    for (int j=0; j<atonce; j++) {
      double *Xj = X + j * aligned_c,
	*Xneuj = Xneu + j * aligned_c;
      for (int i=0; i<c; i++) // see chol2normNoPivot
	Xj[i] = EXP(Ext_scalarX(Xneuj, chol + c * i, i + 1, SCALAR_AVX));
    }
      
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif
   for (int i=0; i<c; i++) {
     //omp_get_num_procs())
#if defined SSE2 
      double
	*v = X,
	*m = REAL(M) + c * i,
	*end = m + c - doubles;
#define DefMax(I) Double ALIGNED Max##I = ZERODOUBLE;
      DefMax(0); DefMax(1); DefMax(2); DefMax(3); DefMax(4); DefMax(5);
      for ( ; m < end; m += doubles, v += doubles) {
	Double mloc = LOADuDOUBLE(m);
#define MaxX(I)\
	Max##I = MAXDOUBLE(Max##I, MULTDOUBLE(mloc, LOADDOUBLE(v+I*aligned_c)));
	MaxX(0); MaxX(1); MaxX(2); MaxX(3); MaxX(4); MaxX(5);
      }
#define MaxS(I)						\
      double *d##I = (double *) &Max##I,		\
	max##I = d##I[0];				\
      max##I = max##I < d##I[1] ? d##I[1] : max##I;
      MaxS(0); MaxS(1); MaxS(2); MaxS(3); MaxS(4); MaxS(5);
#if defined AVX
#define MaxT(I) max##I = max##I < d##I[2] ? d##I[2] : max##I;	\
      max##I = max##I < d##I[3] ? d##I[3] : max##I;
      MaxT(0); MaxT(1); MaxT(2); MaxT(3); MaxT(4); MaxT(5);
#endif      
      end += doubles;
      for ( ; m<end; m++, v++) {	
	double z;
#define MaxU(I) z = *m * v[I * aligned_c]; max##I = max##I > z ? max##I : z;
	MaxU(0); MaxU(1); MaxU(2); MaxU(3); MaxU(4); MaxU(5);
      }

#define Inv(I) double invmax##I = 1.0 / max##I;
      Inv(0); Inv(1); Inv(2); Inv(3); Inv(4); Inv(5);
      v = X;
      m = REAL(M) + c * i;
      double *ans = ans0 + c * i;
      if (!false) {
	double *im;
#define SetInv(I) Double ALIGNED Invmax##I;			\
	im = (double*) &Invmax##I;				\
	for (int u=0; u<doubles; u++) im[u] = invmax##I;	
	SetInv(0); SetInv(1); SetInv(2); SetInv(3); SetInv(4); SetInv(5);
	end = m + c - doubles;
	for ( ; m < end; m += doubles, v += doubles, ans += doubles) {
	  Double mloc = LOADuDOUBLE(m),
	    a = MULTDOUBLE(Invmax0, LOADDOUBLE(v));
#define Add(I) a = ADDDOUBLE(a, MULTDOUBLE(Invmax##I, LOADDOUBLE(v+I*aligned_c)))
	  Add(1); Add(2); Add(3); Add(4); Add(5);
	  STOREuDOUBLE(ans, ADDDOUBLE(LOADuDOUBLE(ans), MULTDOUBLE(mloc, a) ));
	}
	end += doubles;
#define Mult(I) v[I * aligned_c] * invmax##I
	for ( ; m<end; m++, v++, ans++)
	  *ans += *m * (Mult(0) +Mult(1) +Mult(2) +Mult(3) +Mult(4) +Mult(5));
      } else {     
	for (int j=0; j<c; j++)
	  ans[j] += m[j] * (Mult(0) +Mult(1) +Mult(2) +Mult(3) +Mult(4)+Mult(5));
      }
#else
       double
	max = RF_NEGINF,
	*v = X,
	*m = REAL(M) + c * i;
      for (int j=0; j<c; j++) {
	double dummy = m[j] * v[j];	
	max = dummy > max ? dummy : max;
      }
      double invmax = 1.0 / max;

      double * ans  = ans0 + c * i;
      for (int j=0; j<c; j++) ans[j] += m[j] * v[j] * invmax;
      
#endif    
    }
  }

  double invn = 1.0 / (n * atonce);
  for (int i=0; i<c2; i++) ans0[i] *= invn;
  
  FREE(x); FREE(xneu);
  
  PutRNGstate();// PutRNGstate muss vor UNPROTECT stehen!! 
  UNPROTECT(1);
  return Ans;
}
