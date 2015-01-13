/*
 Authors
 Marco Oesting,
 Martin Schlather, schlather@math.uni-mannheim.de

 Kriging methods

 Copyright (C) 2001 -- 2010 Martin Schlather, 
 Copyright (C) 2011 -- 2014 Marco Oesting & Martin Schlather, 

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
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "RF.h"
#include "variogramAndCo.h"

#define KRIGE_TOLERANCE -1e-10

// Marco: bitte die Definitionen in die RF.h schreiben
void matmult(double *A, double *B, double * C, int l, int m, int n);
void matmulttransposed(double *A, double *B, double *C, int m, int l, int n);
void poly_basis_extern(int *Dim, int *Deg, int *powmatrix);

#define STANDARDKRIGING							\
  double *xx = REAL(X),							\
    *krig = REAL(Krig),							\
    *C = NULL,								\
    *invcov = REAL(Invcov),						\
    *tgiven = REAL(Tgiven);						\
  int *notna = LOGICAL(Notna);						\
  int VARIABLE_IS_NOT_USED d,						\
  i, j, k, r, s, v, krigi,						\
    reg = INTEGER(Reg)[0],						\
    dim = INTEGER(Dim)[0],						\
    ngiven = INTEGER(Ngiven)[0],					\
    rep = INTEGER(Rep)[0],						\
    nx = INTEGER(Nx)[0],						\
    err = NOERROR,							\
    vdim = KEY[reg]->vdim[0],						\
    vdimng = vdim * ngiven,						\
    /*   len_tgiven = dim * ngiven,	*/				\
    divachtzig = (nx<79) ? 1 : (nx / 79),				\
    divachtzigM1 = divachtzig - 1;					\
  bool pr=PL > 0 && GLOBAL.general.pch!='\0' && GLOBAL.general.pch!=' ';\
  if ((C = (double*) MALLOC(sizeof(double) * vdimng * vdim))==NULL	\
      /* || (dist = (double*) MALLOC(sizeof(double) * len_tgiven))==NULL */ \
      ) {								\
    err = ERRORMEMORYALLOCATION;					\
    goto ErrorHandling;							\
  }									\

#define STANDARDKRIGING2						\
  double origin[MAXSIMUDIM],						\
    *var = NULL,							\
    *lambda = NULL,							\
    *data = REAL(Data),							\
    *sigma2 = REAL(Sigma2);						\
  int vari;								\
  STANDARDKRIGING;							\
  if ((lambda = (double*) MALLOC(sizeof(double) * vdimng))==NULL ||	\
      (var = (double*) MALLOC(sizeof(double) * vdim * vdim))==NULL) {	\
    err = ERRORMEMORYALLOCATION;					\
    goto ErrorHandling;							\
  }									\
  for (i=0; i<dim; i++) origin[i] = 0.0;				\

#define STANDARD_END				\
  if (C!=NULL) free(C);				\
  /*  if (dist!=NULL) free(dist);*/		\
  if (err != NOERROR) {				\
    int endforX;				\
    endforX = nx * vdim * rep;			\
    for (i=0; i<endforX; krig[i++]=RF_NA);	\
  }						\
  return NULL;
  
#define STANDARD_END2				\
  if (var!=NULL) free(var);			\
  if (lambda!=NULL) free(lambda);		\
  STANDARD_END					

SEXP simpleKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
		   SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig) {
  // kriging variance is not calculated
  // Ngiven : length given = length C vector to be calculated
  // Rep : repetition invcov
  // Nx : number of points in x, length(x) = nx * dim

  STANDARDKRIGING;
  double dummy;

 
  for (i=0; i<nx; i++, xx+=dim) {
    // print("%d %d %d \n", i, nx, divachtzig);
    if (pr && (i % divachtzig == divachtzigM1)) PRINTF("%c", GLOBAL.general.pch);
    //for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    CovIntern(reg, tgiven, xx, ngiven, 1, C);//! fnktniert nicht mit mixed model! 
    // to do: multivariate: stimmt tgiven - xx???

    //   for (r=0; r<ngiven * vdim * vdim; r++) print("%e ", C[r]); print("\n");
    krigi = i;
    for (v=0; v<vdim; v++, krigi+=nx) {
      for(j=0, r=0; r<rep; r++) {
	dummy = 0.0;
	for (s=v*vdimng, k=0; k<vdimng; k++, s++) {
	  if (notna[k]) dummy += C[s] * invcov[j++]; // memory problem
	  // das kann auch inhaltlich ueberhaupt nicht stimmen. Bitte mit
	  // meinem altem Code abgleichen! DONE!

//	    print("%d r=%d v=%d d=%d value=%f %f %f %d \n", 
//		   i, r, v, d, dummy, pC[d], invcov[d], notna[d]);
	}
	krig[krigi + r*vdim*nx] = dummy;  // krigi=i + v * nx + r * nx * vdim
      }	
    }

//    assert(false);
  }
  if (pr)  PRINTF("\n");

 ErrorHandling:
  STANDARD_END;
}


SEXP simpleKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
		    SEXP Invcov, SEXP Notna,
		   SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
		    SEXP Sigma2) {
 // Ngiven : length given = length M vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim

  STANDARDKRIGING2;

  CovIntern(reg, origin, var);
  for (i=0; i<nx; i++, xx+=dim) {
    //  for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    CovIntern(reg, tgiven, xx, ngiven, 1, C); //! fktniert nicht mit mixed model!
    krigi = i;
    if (pr && (krigi % divachtzig == divachtzigM1))
      PRINTF("%c", GLOBAL.general.pch);    
    for (vari=v=0; v<vdim; v++, krigi += nx, vari += vdim+1) {
      for (k=s=0; k<vdimng; k++) {
        lambda[k] = 0.0;
        if (notna[k])
          for (j=0; j<vdimng; j++) 
	    lambda[k] += C[j+v*vdimng] * invcov[s++];
      }
      sigma2[krigi] = var[vari];
      for (j=0; j<vdimng; j++) {
        sigma2[krigi] -= lambda[j] * C[j+v*vdimng];
      }
      if (sigma2[krigi] < 0) {
	// print("%e", sigma2[krigi]);	
	if (sigma2[krigi] < KRIGE_TOLERANCE) {
	  err =  ERRORKRIGETOL;
	  goto ErrorHandling;
	}
	sigma2[krigi] = 0.0;
      }
      for(d=0, r=0; r<rep; r++) {
        krig[krigi + r*vdim*nx] = 0.0;
        for (k=0; k<vdimng; k++)
          if(notna[k]) krig[krigi + r*vdim*nx] += lambda[k] * data[d++];
      }
    }
  }
  if (pr) PRINTF("\n");

 ErrorHandling:
  STANDARD_END2;
}


SEXP ordinaryKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
		   SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig) {

  // kriging variance is not calculated

  int addv;
  double dummy,
    *F = NULL;

  STANDARDKRIGING;

  if ((F = (double*) MALLOC(sizeof(double) * vdim *vdim))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  
  for (i=0; i<vdim; i++)
    for (j=0; j<vdim; j++)
      F[i+j*vdim] = (i==j) ? 1.0 : 0.0;

  for (i=0; i<nx; i++, xx+=dim) {
    if (pr && (i % divachtzig == divachtzigM1)) PRINTF("%c", GLOBAL.general.pch);
    //for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    CovIntern(reg, tgiven, xx, ngiven, 1, C);//! fnktniert nicht mit mixed model!
    
    //   for (r=0; r<ngiven * vdim * vdim; r++) print("%e ", C[r]); print("\n");
    krigi = i;
    for (v=0; v<vdim; v++, krigi+=nx) {
      for(j=0, r=0; r<rep; r++) {
	dummy = 0.0;
	for (s=v*vdimng, k=0; k<vdimng; k++, s++) {
	  if (notna[k]) dummy += C[s] * invcov[j++];
	}
	for(addv=0; addv<vdim; addv++)
	  dummy += F[v*vdim+addv] * invcov[j++];
	krig[krigi+ r*vdim*nx] = dummy;  // krigi=i + v * nx + r * nx * vdim
      }	
    }

//    assert(false);
  }
  if (pr) PRINTF("\n");

 ErrorHandling:
  if (F!=NULL) free(F);
  STANDARD_END;
}


SEXP ordinaryKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
		    SEXP Invcov, SEXP Notna,
		   SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
		    SEXP Sigma2) {
 // Ngiven : length given = length C vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim

  // ZWINGEND vor STANDARDKRIGING2 !
  int addv;
  double *F = NULL,
     *mu = NULL;

  STANDARDKRIGING2;

  if ((F = (double*) MALLOC(sizeof(double) * vdim * vdim))==NULL || 
      (mu = (double*) MALLOC(sizeof(double) * vdim))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  
  CovIntern(reg, origin, var);
  for(i=0; i<vdim; i++)
    for(j=0; j<vdim; j++)
       F[i+j*vdim] = (i==j) ? 1.0 : 0.0;
 
  for (i=0; i<nx; i++, xx+=dim) {
    krigi = i;    
    if (pr && (krigi % divachtzig == divachtzigM1)) 
      PRINTF("%c", GLOBAL.general.pch);
    //for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    CovIntern(reg, tgiven, xx, ngiven, 1, C);//!!fktniert nicht mit mixed model!
    if (pr && (krigi % divachtzig == divachtzigM1)) 
      PRINTF("%c", GLOBAL.general.pch);    
    for (vari=v=0; v<vdim; v++, krigi += nx, vari += vdim+1) {
      for (k=s=0; k<vdimng; k++) {
        lambda[k] = 0.0;
        if (notna[k])
          for (j=0; j<vdimng; j++) 
	    lambda[k] += C[j+v*vdimng] * invcov[s++];
        for(addv=0; addv<vdim; addv++)
	  lambda[k] += F[addv+v*vdim] * invcov[s++];
      }
      for(k=0; k<vdim; k++) {
        mu[k] = 0.0;
        if (notna[k])
          for (j=0; j<vdimng; j++) 
	    mu[k] += C[j+v*vdimng] * invcov[s++];
        for(addv=0; addv<vdim; addv++)
	  mu[k] += F[addv+v*vdim] * invcov[s++];	
      }
      
      sigma2[krigi] = var[vari];
      for (j=0; j<vdimng; j++) {
        sigma2[krigi] -= lambda[j] * C[j+v*vdimng];
      }
      for(j=0; j<vdim; j++)
	sigma2[krigi] -= mu[j] * F[j+v*vdim];
      if (sigma2[krigi] < 0) {
	// print("%e", sigma2[krigi]);
	if (sigma2[krigi] < KRIGE_TOLERANCE) {
	  err =  ERRORKRIGETOL;
	  goto ErrorHandling;
	}
        sigma2[krigi] = 0.0;
      }
      
      for(d=0, r=0; r<rep; r++) {
        krig[krigi + r*vdim*nx] = 0.0;
        for (k=0; k<vdimng; k++)
          if(notna[k]) krig[krigi + r*vdim*nx] += lambda[k] * data[d++];
      }
    }
  }
  if (pr) PRINTF("\n");

 ErrorHandling:
  if (F!=NULL) free(F);
  if (mu!=NULL) free(mu);
  STANDARD_END2;
}



SEXP universalKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
		      SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
		      SEXP Nfct, SEXP trend_expr, SEXP trend_envir) {
  // kriging variance is not calculated
  // Ngiven : length given = length C vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  int f,
    nfct = INTEGER(Nfct)[0];
  
  double dummy,
     *F = NULL;

  STANDARDKRIGING;
    
  if ((F = (double*) MALLOC(sizeof(double) * nfct * vdim))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }

  SEXP trendargs, trendres;
  PROTECT(trendargs = allocVector(REALSXP,dim));

  for (i=0; i<nx; i++, xx+=dim) {
    krigi = i;
    if (pr && (krigi % divachtzig == divachtzigM1))
      PRINTF("%c", GLOBAL.general.pch);
    //for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    //Initializing C
    CovIntern(reg, tgiven, xx, ngiven, 1, C);
    //install R variables
    for (d=0; d<dim; d++) REAL(trendargs)[d] = xx[d];

    defineVar(install("trendargs"), trendargs, trend_envir);
    PROTECT(trendres = coerceVector(eval(trend_expr, trend_envir),REALSXP));
    for (j=0; j<nfct*vdim; j++) F[j] = REAL(trendres)[j];
    UNPROTECT(1);

    for (v=0; v<vdim; v++, krigi+=nx) {
      for (j=0, r=0; r<rep; r++) {
	dummy = 0.0;
	for (s=v*vdimng, k=0; k<vdimng; k++, s++) {
	  if (notna[k]) dummy += C[s] * invcov[j++];
	}
	for (f=0; f<nfct; f++)
	  dummy += F[v*nfct+f] * invcov[j++];
	krig[krigi + r*vdim*nx] = dummy;  // krigi=i + v * nx + r * nx * vdim
      }
    }
  }
  UNPROTECT(1);
  if (pr) PRINTF("\n");

 ErrorHandling:
  if (F!=NULL) free(F);
  STANDARD_END;
}


SEXP universalKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
		       SEXP Invcov, SEXP Notna,
		       SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
		       SEXP Sigma2, 
		       SEXP Nfct, SEXP trend_expr, SEXP trend_envir) {
 // also calculating kriging variance
  int nfct = INTEGER(Nfct)[0];
 // Ngiven: number of points given = length M vector to be calculated
 // Rep: repitition data
 // Nx: number of points in x, length(x) = nx * dim
  double  
    *Fmatrix = NULL,
    *fvector = NULL,
    *X1 = NULL,
    *lambdak = NULL,
    *Qmatrix = NULL,
    *Rvector = NULL,
    *mu = NULL;
  SEXP trendargs, trendres;  
  
  STANDARDKRIGING2;

  if ((Fmatrix = (double*) MALLOC(sizeof(double) * vdimng * nfct))==NULL||
      (fvector = (double*) MALLOC(sizeof(double) * nfct * vdim))==NULL||
      (X1 = (double*) MALLOC(sizeof(double) * vdimng * nfct))==NULL ||
      (lambdak = (double*) MALLOC(sizeof(double) * vdimng * vdim))==NULL ||
      (Qmatrix = (double*) MALLOC(sizeof(double) * nfct * nfct))==NULL ||
      (Rvector = (double*) MALLOC(sizeof(double) * nfct * vdim))==NULL ||
      (mu = (double*) MALLOC(sizeof(double) * nfct * vdim))==NULL) {
    err=ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  
  PROTECT(trendargs = allocVector(REALSXP,dim));
  //Initialize Fmatrix
  for (s=0, i=0; i<ngiven; i++) {
    for (d=0; d<dim; d++) REAL(trendargs)[d] = tgiven[s++];
    defineVar(install("trendargs"), trendargs, trend_envir);
    PROTECT(trendres = coerceVector(eval(trend_expr, trend_envir),REALSXP));
    for (v=0, k=0; v<vdim; v++)
      for (j=0; j<nfct; j++)
        Fmatrix[j*vdimng + i + v*ngiven] = REAL(trendres)[k++];
    UNPROTECT(1);
  }
  
  //Calculating X1 and Q, inverting Q
  matmult(invcov, Fmatrix, X1, vdimng, vdimng, nfct);
  matmulttransposed(Fmatrix, X1, Qmatrix, vdimng, nfct, nfct);
  if ((err = invertMatrix(Qmatrix, nfct)) > NOERROR) goto ErrorHandling;

  //Calculating the variance
  CovIntern(reg, origin, var);
  
  //Kriging for each kriging point

  for (i=0; i<nx; i++, xx+=dim) {
    krigi = i;
    if (pr && (krigi % divachtzig == divachtzigM1))
      PRINTF("%c", GLOBAL.general.pch);
    //Determining x
    //Calculating C vector
    //for(j=d=0; j<len_tgiven; j++, d= (d+1) % dim) dist[j] = tgiven[j] - xx[d];
    CovIntern(reg, tgiven, xx, ngiven, 1, C); 
    //Calculating fvector
    for (d=0; d<dim; d++) REAL(trendargs)[d] = xx[d];
    defineVar(install("trendargs"), trendargs, trend_envir);
    PROTECT(trendres = coerceVector(eval(trend_expr, trend_envir),REALSXP));   
    for (k=0; k<nfct*vdim; k++) fvector[k] = REAL(trendres)[k];
    UNPROTECT(1);
    //Solving the kriging equations
    //Calculating lambak
    matmult(invcov, C, lambdak, vdimng, vdimng, vdim);
    //Initializing R
    matmulttransposed(Fmatrix, lambdak, Rvector, vdimng, nfct, vdim);
    int endfor = nfct*vdim;
    for(k=0; k<endfor; k++) Rvector[k] -= fvector[k];
    //Calculating mu
    matmult(Qmatrix, Rvector, mu, nfct, nfct, vdim);
    //Calculating lambda
    matmult(X1, mu, lambda, vdimng, nfct, vdim);
    endfor = vdimng*vdim;
    for(k=0; k<endfor; k++) lambda[k] = lambdak[k] - lambda[k];
    for(vari=v=0; v<vdim; v++, vari+=vdim+1, krigi+=nx) {
      sigma2[krigi] = var[vari];
      for (j=0; j<vdimng; j++) {
	sigma2[krigi] -= lambda[j+v*vdimng] * C[j+v*vdimng];
      }
      for (j=0; j<nfct; j++)
	sigma2[krigi] -= mu[j+v*nfct] * fvector[j+v*nfct];
      if (sigma2[krigi] < 0) {
	// print("%e", sigma2[krigi]);
	if (sigma2[krigi] < KRIGE_TOLERANCE) {
	  err=ERRORKRIGETOL;
	  goto ErrorHandling;
	}
	sigma2[krigi] = 0.0;
      }
      for (d=0, r=0; r<rep; r++) {
	krig[krigi + r*vdim*nx] = 0.0;
	for (k=0; k<vdimng; k++) {
	  if(notna[k]) 
	    krig[krigi + r*vdim*nx] += lambda[k+v*vdimng] * data[d++];
	}
      }
    }
  }
  UNPROTECT(1);
  if (pr) PRINTF("\n");
  
 ErrorHandling:
  if (Fmatrix!=NULL) free(Fmatrix);
  if (fvector!=NULL) free(fvector);
  if (X1 !=NULL) free(X1);
  if (lambdak !=NULL) free(lambdak);
  if (Qmatrix!=NULL) free(Qmatrix);
  if (Rvector!=NULL) free(Rvector);
  if (mu!=NULL) free(mu);
  STANDARD_END2;
}



SEXP intrinsicKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
		      SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig,
		      SEXP Polydeg) { 
  // kriging variance is not calculated
 // Ngiven : length given = length C vector to be calculated
  // Rep : repetition data
  // Nx : number of points in x, length(x) = nx * dim
  int  f,
    *polydeg = INTEGER(Polydeg),
    *powmatrix = NULL;

  double dummy,
     *F = NULL;

  STANDARDKRIGING;

  int nfct, vdimnf;
  nfct = binomialcoeff(dim + *polydeg, dim);
  vdimnf = vdim * nfct;

  if ((F = (double*) MALLOC(sizeof(double) * vdimnf * vdim))==NULL ||
      (powmatrix = (int*) MALLOC(sizeof(int) * nfct * dim))==NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }

  poly_basis_extern(&dim, polydeg, powmatrix);
  
  for (i=0; i<nx; i++, xx+=dim) {
    krigi = i;
    if (pr && (krigi % divachtzig == divachtzigM1))
      PRINTF("%c", GLOBAL.general.pch);
    //for (j=d=0; j<len_tgiven; j++, d=(d+1) % dim) dist[j] = tgiven[j] - xx[d];
    //Initializing Cf0
    PseudovariogramIntern(reg, tgiven, xx, ngiven, 1, C);
    for (j=0; j<vdimng*vdim; j++) C[j] = (-1.0)*C[j];
    for (k=0; k<vdimnf*vdim; k++) F[k] = 0.0;
    for (k=0, v=0; v<vdim; v++, k+=vdimnf) {
      for (j=0; j<nfct; j++, k++) {
        F[k] = 1.0;
        for (d=0; d<dim; d++) F[k] *= pow(xx[d],powmatrix[j*dim+d]);
      }
    }
    for (v=0; v<vdim; v++, krigi+=nx) {
      for (j=0, r=0; r<rep; r++) {
	dummy = 0.0;
	for (s=v*vdimng, k=0; k<vdimng; k++, s++) {
	  if (notna[k]) dummy += C[s] * invcov[j++];
	}
	for (f=0; f<vdimnf; f++)
	  dummy += F[v*vdimnf+f] * invcov[j++];
	krig[krigi + r*vdim*nx] = dummy;  // krigi=i + v * nx + r * nx * vdim
      }
    }
  }
  if (pr) PRINTF("\n");

 ErrorHandling:
  if (F!=NULL) free(F);
  if (powmatrix!=NULL) free(powmatrix);
  STANDARD_END;
}

SEXP intrinsicKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
		    SEXP Invcov, SEXP Notna,
		   SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
		    SEXP Sigma2, SEXP Polydeg) {
 // Ngiven: number of points given = length C vector to be calculated
  // Rep: repitition data
  // Nx: number of points in x, length(x) = nx * dim
  // NrowPowers: Number of functions to form the trend

  // Matrix-Inversion .Call("La_dgesv", a, b, tol, PACKAGE = "base")
  // F77_CALL(dgesv)
  // siehe auch static SEXP modLa_dgesv(SEXP A, SEXP Bin, SEXP tolin)
  // oder .Fortran("dqrcf",
  // also calculating kriging variance
  int 
    *powmatrix = NULL,
    *polydeg = INTEGER(Polydeg);

  double 
    *Fmatrix = NULL,
    *fvector = NULL,
    *X1 = NULL,
    *lambdak = NULL,
    *Qmatrix = NULL,
    *Rvector = NULL,
    *mu = NULL;
  
  STANDARDKRIGING2;

  int nfct, vdimnf, endfor;
  nfct = binomialcoeff(dim + *polydeg, dim);
  vdimnf = vdim * nfct;

  if ((Fmatrix = (double*) MALLOC(sizeof(double) * vdimng * vdimnf))==NULL||
      (fvector = (double*) MALLOC(sizeof(double) * vdimnf  * vdim))==NULL||
      (X1 = (double*) MALLOC(sizeof(double) * vdimng * vdimnf))==NULL ||
      (lambdak = (double*) MALLOC(sizeof(double) * vdimng * vdim))==NULL ||
      (Qmatrix = (double*) MALLOC(sizeof(double) * vdimnf * vdimnf))==NULL ||
      (Rvector = (double*) MALLOC(sizeof(double) * vdimnf * vdim))==NULL ||
      (mu = (double*) MALLOC(sizeof(double) * vdimnf * vdim))==NULL ||
      (powmatrix = (int*) MALLOC(sizeof(int) * nfct * dim))==NULL) {
    err=ERRORMEMORYALLOCATION;
   goto ErrorHandling;
  }
  
  poly_basis_extern(&dim, polydeg, powmatrix);

  //Calculate Fmatrix
  endfor = vdimng*vdimnf;
  for(s=0; s<endfor; s++) Fmatrix[s] = 0.0; 
  double *PFmatrix;
  PFmatrix = Fmatrix;
  for (v=0; v<vdim; v++, PFmatrix+=ngiven) {
    // MARTIN: stimmt dies? Unten wird PFmatrix++ ngiven mal aufgerufen
    for (j=0; j<nfct; j++, PFmatrix+=(vdim-1)*ngiven) {//for each function
      for(i=0; i<ngiven; i++, PFmatrix++) { //for each given point
	*PFmatrix = 1.0;
	for(d=0; d<dim; d++) //for each component
	  *PFmatrix *= intpow(tgiven[i*dim + d], powmatrix[j*dim + d]);
      }
    }
  }
  
  
  if (false) {
    printf("in kriging.cc:\n");//
    int fi,fj,fk, size=vdimng;
    for (fk=fi=0; fi<size; fi++) {
      for (fj=0; fj<size; fj++) {
	printf("%3.3f ", invcov[fk++]); // /* check if false */
      }
      printf("\n"); //
    }
    printf("%d %d \n", vdimng, vdimnf);//
 }
  
  
  //Calculating X1 and Q, inverting Q
  matmult(invcov, Fmatrix, X1, vdimng, vdimng, vdimnf);
  matmulttransposed(Fmatrix, X1, Qmatrix, vdimng, vdimnf, vdimnf);
  if ((err = invertMatrix(Qmatrix, vdimnf)) > NOERROR) {     
   //   AERR(err);
    goto ErrorHandling;
  }
  
  /*
    int invertMethod[MAXINVERSION];
    for (i=0; i<MAXINVERSION; i++)
   invertMethod[i] = GLOBAL.general.matrix_inversion[i];
   for (i=0; i<MAXINVERSION; i++) {
   if (invertMethod[i] <= 0) {
   }
   }
  */
  
  
  /*
  //
  //bei variogrammen ist Q i.A. nicht strikt positive definit
  SEXP invMatrix;
  PROTECT(invMatrix = allocMatrix(REALSXP, vdimnf, vdimnf));
  for(j=0; j<vdimnf*vdimnf; j++)  REAL(invMatrix)[j] = Qmatrix[j];
  defineVar(install("Q"), invMatrix, envir);
  invMatrix = eval(inv_expr, envir);
  for(j=0; j<vdimnf*vdimnf; j++)  Qmatrix[j] = REAL(invMatrix)[j];
  UNPROTECT(1);
  */
  

  //Calculating the variance
  PseudovariogramIntern(reg, origin, var);
  endfor = vdim*vdim;
  for(v=0; v<endfor; var[v++] *= -1.0);
  
  //Kriging for each kriging point
  for (i=0; i<nx; i++, xx+=dim) {
    krigi = i;
    if (pr && (krigi % divachtzig == divachtzigM1))
      PRINTF("%c", GLOBAL.general.pch);
    //Determining x
    //Calculating C vector
    //for(j=d=0; j<len_tgiven; j++, d= (d+1) % dim) dist[j] = tgiven[j] - xx[d];
    PseudovariogramIntern(reg, tgiven, xx, ngiven, 1, C);
    endfor = vdimng*vdim;
    for(j=0; j<endfor; j++) C[j] *= -1.0;
    //Calculating fvector
    endfor = vdimnf*vdim;
    for(k=0; k<endfor; k++) fvector[k] = 0.0;
    for(k=0, v=0; v<vdim; v++, k+=vdimnf) {
      for(j=0; j<nfct; j++, k++) { //for each function
	fvector[k] = 1.0;
	for(d=0; d<dim; d++) //for each component
	  fvector[k] *= intpow(xx[d], powmatrix[j*dim + d]);
      }
    }
    
    //Solving the kriging equations
    //Calculating lambak
    matmult(invcov, C, lambdak, vdimng, vdimng, vdim);
    //Initializing R
    matmulttransposed(Fmatrix, lambdak, Rvector, vdimng, vdimnf, vdim);
    endfor = vdimnf*vdim;
    for(k=0; k<endfor; k++) Rvector[k] -= fvector[k];
    //Calculating mu
    matmult(Qmatrix, Rvector, mu, vdimnf, vdimnf, vdim);
    //Calculating lambda
    matmult(X1, mu, lambda, vdimng, vdimnf, vdim);
    endfor = vdimng * vdim;
    for(k=0; k<endfor; k++) lambda[k] = lambdak[k] - lambda[k];
    for(vari=v=0; v<vdim; v++, vari+=vdim+1, krigi+=nx) {
      sigma2[krigi] = var[vari];
      for (k=0; k<vdimng; k++) {
	if(notna[k]) sigma2[krigi] -= lambda[k+v*vdimng] * C[k+v*vdimng];
      }
      for (j=0; j<vdimnf; j++)
	sigma2[krigi] -= mu[j+v*vdimnf] * fvector[j+v*vdimnf];
      if (sigma2[krigi] < 0) {
	//printf("krig %d %d %f\n", krigi, nx, sigma2[krigi]);
	if (sigma2[krigi] < KRIGE_TOLERANCE) {
	  err=ERRORKRIGETOL;
	  goto ErrorHandling;
	} 
	sigma2[krigi] = 0.0;
      }
      for (d=0, r=0; r<rep; r++) {
	krig[krigi + r*vdim*nx] = 0.0;
	for (k=0; k<vdimng; k++)
	  if(notna[k]) krig[krigi+r*vdim*nx] += lambda[k+v*vdimng] * data[d++];
      }
    }
  }
  if (pr) PRINTF("\n");
  
 ErrorHandling:
  if (Fmatrix!=NULL) free(Fmatrix);
  if (fvector!=NULL) free(fvector);
  if (X1 !=NULL) free(X1);
  if (lambdak !=NULL) free(lambdak);
  if (Qmatrix!=NULL) free(Qmatrix);
  if (Rvector!=NULL) free(Rvector);
  if (mu!=NULL) free(mu);
  if (powmatrix!=NULL) free(powmatrix);
  STANDARD_END2;
}



void poly_basis_extern(int *Dim, int *Deg, int *powmatrix) {
  
  int powsum, d, i, j, 
    err = NOERROR,
    *dimi=NULL,
    dim = *Dim, deg = *Deg,
    basislen = binomialcoeff(deg+dim,dim);
     
  dimi = (int *) MALLOC(dim * sizeof(int));
  if (dimi == NULL) {
     err = ERRORMEMORYALLOCATION;
     goto ErrorHandling;
  }
      
  for(d=0; d<dim; dimi[d++] = 0);
  for(j=0; j<basislen; j++) {
    for(d=0; d<dim; d++) powmatrix[j*dim+d] = dimi[d];
    i=0;
    (dimi[i])++;
    powsum = 0;
    for(d=0; d<dim; d++) powsum += dimi[d];
    while(powsum > deg) {
      dimi[i] = 0;
      if (i < dim-1) {
        i++;
	(dimi[i])++;
      }
      powsum = 0;
      for(d=0; d<dim; d++) powsum += dimi[d];
    }
  }
  
  ErrorHandling :
  if (dimi != NULL) free(dimi);   
  dimi = NULL;
  if (err != NOERROR) XERR(err);

  return;
}
