
/*
 Authors
 Martin Schlather, martin.schlather@cu.lu

 Simulation of a random field by turning bands;
 see RFspectral.cc for spectral turning bands

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
//#include <sys/timeb.h>
#include <unistd.h>
#include <assert.h>
#include "RFsimu.h"
#include "RFCovFcts.h"

#define MAXNN 100000000.0 /* max number of points simulated on a line */

typedef struct tbm_lines{
  int lines;           // number of lines simulated
  Real linesimufactor, /*factor by which simulation on the line is finer 
			 than the grid */
    linesimustep;      /* grid lag defined absolutely */
    int every; 
} tbm_lines;

tbm_lines tbm2 =   { 60, 2.0, 0.0, 0};
tbm_lines tbm3d2 = {500, 2.0, 0.0, 0};
tbm_lines tbm3d3 = {500, 2.0, 0.0, 0};
    
#ifdef doubleReal 
ce_param TBMCE = {false, true, TRIVIALSTARTEGY, -1e-7, 1e-3, 3, 0, 0, 0, 0};
#else
ce_param TBMCE = {false, true, TRIVIALSTARTEGY, -1e-5, 1e-3, 3, 0, 0, 0, 0};
#endif

SimulationType TBM_METHOD=CircEmbed;

typedef struct TBM_storage {
  simu1dim linemethod; /* simulation method on the line (up to now only 
		          circulantembedding) */
  SimulationType method;
  int nn[2], m[2], halfm[2], nnhalf, cumm[3]; /*** !!! ***/
  Real linesimuscale, half[MAXDIM], *x;
  int mtot; 
  double *c,*d;
  Real *simuline;
  FFT_storage FFT;
  direct_storage *Direct;
  int simuspatialdim, timespacedim, ce_dim;
  bool grid;
  
} TBM_storage;


double TBM2_INTEGR_PREC=1e-9;
#define NNN 36 /* immmer gerade !! */
Real TBM2integr(int covnr, Real h, Real v, Real *param, int dim, int fct) {
  Real de, ialt, ineu, x, y, z, factor[3] = {0.5, 2.0, 0.5}, 
    xx[3], b[NNN+1], integral, delta, alpha, f, h2, v2, rest; 
  unsigned long int n, k;
  int i, mitte, kk;
    
  h2 = h * h;
  v2 = v * v;
  b[0] = 0.0;
  b[NNN] = 1 - 1e-13;
  mitte = 1;
  b[mitte] = 0.7;
  alpha = 0.19;
  for (i=mitte + 1; i<NNN ; i++) {
    b[i] = exp(0.5 * log(alpha + (1.0 - alpha) * b[i-1]));
  }
  // alpha = 0.7;
  for (i=mitte - 1; i>0 ; i--) {
    b[i] = (double) i / (double) mitte * b[mitte];
  }

  integral = 0.0;
  for (i=0; i<NNN; i++) {
//    printf(" i=%d\n ",i);
    n = 2;
    delta = b[i+1] - b[i];
    de = delta / (double) n;
    ialt = 0.0;
    for (k=0, x=b[i]; k<=2; k++, x+=de) {
      switch(fct) {
	  case 0 : 
	    z = h * x;
	    f = CovList[covnr].cov_tbm3(&z, param, dim);
	    break;
	  case 1 :
	    z = fabs(h * x);
	    f = CovList[covnr].cov(&z, param, dim) +
	      z * CovList[covnr].ableitung(&z, param, dim);
	    break;
	  case 2 :
	    y = x * x * h2;
	    z = sqrt(y + v2); 
	    if (z==0.0) f=1.0; 
	    else f = CovList[covnr].cov(&z, param, dim) +
	      y / z * CovList[covnr].ableitung(&z, param, dim);
	    break;
	  default : break;
	}
      ialt += factor[k]  * x / sqrt(1 - x * x) * f;
//  printf(" %d x=%f %f \n", n, x, ialt * h);
      rest = f * sqrt(1 - b[NNN] * b[NNN]);
    }
    
    while (true) {
      n *= 3;
      de = delta / (double) n;
      ineu = 0.0;
      for (k=1, x=b[i]+de; k<n; k++, x+=de) {
	kk = k % 6;
	if (kk == 0 || kk == 3) continue;
	switch (fct) {
	  case 0 : 
	    z = h * x;
//	    printf(" z=%f\n", z);
	    f = CovList[covnr].cov_tbm3(&z, param, dim);
	    break;
	  case 1 :
	    z = fabs(h * x);
	    f = CovList[covnr].cov(&z, param, dim) +
	      z * CovList[covnr].ableitung(&z, param, dim);
//	    printf("%d %f zz=%f %f %f %f\n", n, x ,z, f, 
//		   CovList[covnr].cov(&z, param, dim),
//		   CovList[covnr].ableitung(&z, param, dim));
	    break;
	  case 2 :
	    y = x * x * h2;
	    z = sqrt(y + v2); 
	    if (z==0.0) f=1.0; 
	    else f = CovList[covnr].cov(&z, param, dim) +
		   y / z * CovList[covnr].ableitung(&z, param, dim);
	    break;
	  default : break;
	}
	if (kk==1 || kk==5) f *= 2.0;
	ineu += x / sqrt(1 - x * x) * f; 
//	printf(" %f %f %f f=%f %f %f %f i=%f\n", 
//	       x, y, z , f, CovList[covnr].cov(&z, param, dim),
//	       y / z,  CovList[covnr].ableitung(&z, param, dim), ineu);
      }
      if ((n>6) && (fabs(ineu - 2 * ialt) * de < TBM2_INTEGR_PREC)) break;
if (n>1000) printf(" %d [%5.3e %5.3e] %e %e %e h=%e v=%e\n",
		   n, 1-b[i], 1-b[i+1], ineu * de, ialt * de,
		  fabs(ineu - 2 * ialt) * de, h, v);
//     assert (n<=10000) ;
      ialt += ineu;
//   assert(false);
    }
//    printf("\n %d [%5.3e %5.3e] %e; y=%f z=%f %e f=%f",
//	   n, 1-b[i], 1-b[i+1], ialt * de / 1.5,
//	   y, z,  fabs(ineu - 2 * ialt) * de, f);
   ialt += ineu;
   integral += ialt * de / 1.5;
  }
//  assert(false);
//printf("%d a=%f h=%f S=%f \n", n, aa, aa * c, integral);
  return integral + rest;
}


double TBM2_INTEGR_TIME_PREC=1e-9;
Real TBM2integr_time(Real aa, Real c, covfct f, Real *param, int dim) {
#define NN 26 /* immmer gerade !! */
  Real h, ialt, ineu, bb=1.0, x, y, z, factor[3]= {0.5, 2.0, 0.5}, 
    xx[3], b[NN+1], integral, delta, alpha; 
  unsigned long int n,k;
  int i, mitte;

  b[NN] = bb;
  b[0] = aa;
  mitte = NN / 2 - 3;
  b[mitte] = 0.5 * aa + 0.5 * bb;
  alpha = 0.65;
  for (i=mitte + 1; i<NN ; i++) {
    b[i] = bb * alpha + (1.0 - alpha) * b[i-1];
  }
  alpha = 0.7;
  for (i=mitte - 1; i>0 ; i--) {
    b[i] = aa * alpha + (1.0 - alpha) * b[i+1];
  }

  integral = 0.0;
  for (i=0; i<NN; i++) {
    n = 2;
    delta = b[i+1] - b[i];
    h = delta / (double) n;
    ialt = 0.0;
    for (k=0, x=b[i]; k<=2; k++, x+=h) {
      y = x * x;
      z = c * sqrt(1.0 - y); 
      ialt += factor[k] * f(&z, param, dim) / y;
//  printf(" %d x=%f %f \n", n, x, ialt * h);
    }
    
    while (true) {
      n *= 3;
      h = delta / (double) n;
      ineu = 0.0;
      for (k=1, x=b[i]+h; k<n; k++, x+=h) {
	switch (k % 6) {
	    case 0: case 3 : break;
	    case 1: case 5 : 
	      y = x * x; 
	      z = c * sqrt(1.0 - y); 
	      ineu += 2.0 * f(&z, param, dim) / y;
	      break;
	    case 2: case 4 : 
	      y = x * x; 
	      z = c * sqrt(1.0 - y); 
	      ineu += f(&z, param, dim) / y;
	      break;
	}
      }
      if ((n>6) && (fabs(ineu - 2 * ialt) * h < TBM2_INTEGR_TIME_PREC)) break;
//    printf(" %d %f \n", n, ineu * h);
      ialt += ineu;
//    assert(false);
    }
//    printf("\n %d [%f %f] %e; %f %f %f", n, b[i], b[i+1], ialt * h / 1.5,
//	  y, z, f(&z, param, dim));
    ialt += ineu;
    integral += ialt * h / 1.5;
  }
//printf("%d a=%f h=%f S=%f \n", n, aa, aa * c, integral);
  return integral;
}


#define ANISOP3 ANISO+3  
Real CovFctTBM2(Real *x, int dim, int *covnr, int *op,
	        param_type param, int ncov, bool anisotropy){
// in dieser und der naechste Funktionen hat $x$ entweder
// 1 Element (Geradenpunkt der turning bands) oder 2 Elemente
// (Punkt auf Flaeche der turning layers)
  int v, vold, w;
  Real zw, result, r;
  result = 0.0;
  v = 0;

  if (dim==2) { // => s->anisotropy! (time-space simulation)
    Real z[2]; //'' !! 2 !
    while (v<ncov) {
      vold = v;
      v++; while ((v<ncov) && op[v-1]) v++; 
      zw = 1.0;
      for (w=vold; w<v; w++) { 
	// there could be multiplication by pure time components (z[0]=0)!
	z[0] = param[w][ANISO] * x[0];
	z[1] = param[w][ANISOP3] * x[1];
	// printf("z=%f %f\n", z[0], z[1]);
	if (CovList[covnr[w]].isotropic==FULLISOTROPIC) {
	  if (z[0]==0.0) {
	    zw *= param[w][VARIANCE] *
	      CovList[covnr[w]].cov(&(z[1]), param[w], dim); 
	  } else { // called only once in the loop
	    z[0] = fabs(z[0]);
	    r = sqrt(z[0] * z[0] + z[1] * z[1]);
	    zw *= param[w][VARIANCE] * z[0] / r *
	      ( CovList[covnr[w]].cov_tbm2(&r, param[w], dim) +
		TBM2integr_time(z[0] / r, r, CovList[covnr[w]].cov,
				param[w], dim) );

	    printf("\nh=%f v=%f, %f %f, c=%f %e %e\n", z[0], z[1], 
	     CovList[covnr[w]].cov_tbm2(&r, param[w], dim)
		   ,
	     TBM2integr_time(z[0] / r, r, CovList[covnr[w]].cov, param[w], dim)
		   ,
	     CovList[covnr[w]].cov_tbm2(&r, param[w], dim) +
	     TBM2integr_time(z[0] / r, r, CovList[covnr[w]].cov, param[w], dim)
		   ,
	     CovList[covnr[w]].cov_tbm2(&r, param[w], dim) +
	     TBM2integr_time(z[0] / r, r, CovList[covnr[w]].cov, param[w], dim)
		   - TBM2integr(covnr[w], z[0], z[1], param[w], dim, 2),
		   TBM2integr(covnr[w], z[0], z[1], param[w], dim, 2)
	);

	  }
	} else { 
	  assert(CovList[covnr[w]].isotropic==SPACEISOTROPIC);
	  if (z[0]==0.0) 
	    zw *=  CovList[covnr[w]].cov(z, param[w], 2);
	  else if (CovList[covnr[w]].cov_tbm2 != NULL) {
	    zw *= CovList[covnr[w]].cov_tbm2(z, param[w],2);
	   } else {
	    assert(CovList[covnr[w]].ableitung != NULL);
	    zw *= TBM2integr(covnr[w], z[0], 0, param[w], 2,
			     CovList[covnr[w]].cov_tbm3 == NULL);
	  }
	  zw *= param[w][VARIANCE];
	}
      }
      result += zw;
    }
  } else {
    assert(dim==1);
    Real z;
    for (v=0; v<ncov; v++){
      // the field must be isotropic, so the x's have been transformed
      // multiplication of covariance functions is not allowed,
      // in contrast to 3dim TBM
    z = x[0] * param[v][ANISO];

    if (z==0.0) zw = 1.0;
    else if (CovList[covnr[v]].cov_tbm2 != NULL) {
      zw = CovList[covnr[v]].cov_tbm2(&z, param[v], dim);
    } else {
      assert(CovList[covnr[v]].ableitung != NULL);
      zw = TBM2integr(covnr[v], z, 0, param[v], 2,
		       CovList[covnr[v]].cov_tbm3 == NULL);
    }
    result += param[v][VARIANCE] * zw;
    }
  }
  return result;
}


Real CovFctTBM3(Real *x, int dim, int *covnr, int *op,
	        param_type param, int ncov, bool anisotropy){
  int v, vold, w, y;
  Real dummy, zw, result;
  Real r, z[2]; //'' !! 2 !
  Real cov[MAXCOV], abl[MAXCOV];
  result = 0.0;
  v = 0;

  if (dim==2) { // => s->anisotropy! 
    while (v<ncov) {
      vold = v;
      v++; while ((v<ncov) && op[v-1]) v++;
      zw = 1.0;
      for (w=vold; w<v; w++) {
	z[0] = param[w][ANISO] * x[0];
	z[1] = param[w][ANISOP3] * x[1];
	if (CovList[covnr[w]].isotropic==FULLISOTROPIC) {
	  r = sqrt(z[0] * z[0] + z[1] * z[1]);
	  if (r==0) abl[w] = 0;
	  else abl[w] = param[w][VARIANCE] * param[w][ANISO] * fabs(z[0]) / r * 
		 CovList[covnr[w]].ableitung(&r, param[w], dim);
	  zw *= (cov[w] =param[w][VARIANCE]*CovList[covnr[w]].cov(&r, param[w], 
								dim));
	} else { 
	  assert(CovList[covnr[w]].isotropic==SPACEISOTROPIC);
	  abl[w] = param[w][VARIANCE] * 
	    param[w][ANISO] * CovList[covnr[w]].ableitung(z, param[w], dim);
	  zw *= (cov[w] =param[w][VARIANCE]*CovList[covnr[w]].cov(z, param[w], 
								dim));
	}
      }
      result += zw;
      if (x[0]!=0.0) {
	// TBM-OP : C + h C' -- C'(0) is often an exception for calculation
	//          but h C' is always zero, so do not consider z[0]=0.0
	dummy = 0.0;
	for (w=vold; w<v; w++) {
	  zw = abl[w];
	  for (y=vold; y<v; y++) if (y!=w) zw *= cov[y];
	  dummy += zw;
	}
	result += fabs(x[0]) * dummy;
      }
    }
  } else {
    assert(dim==1);
    z[1] = 0.0; //in case of (CovList[covnr[v]].isotropic==SPACESOTROPIC) 
     if (ncov==1) {
     z[0] = x[0] * param[0][ANISO];
     return param[v][VARIANCE]*CovList[covnr[0]].cov_tbm3(z, param[0], dim);
    }
    while (v<ncov) {
      vold = v;
      v++; while ((v<ncov) && op[v-1]) v++;
      zw = 1.0;
      for (w=vold; w<v; w++) {
	z[0] = x[0] * param[w][ANISO];
	zw *= (cov[w] = 
	       param[w][VARIANCE] * CovList[covnr[w]].cov(z, param[w], dim)); 
 	abl[w] = param[w][VARIANCE]*CovList[covnr[w]].ableitung(z, param[w],
								dim);
      }
      result += zw;
      if (z[0]!=0.0) {
	dummy = 0.0;
	for (w=vold; w<v; w++) {
	  zw = abl[w];
	  for (y=vold; y<v; y++) if (y!=w) zw *= cov[y];
	  dummy += zw;
	}
	result += fabs(z[0]) * dummy;
      }
    }
  }
  return result;
}

void do_tbmcirculantembedding(key_type *key, int m, Real *res ) {
  TBM_storage *s;
  s = (TBM_storage*) key->S[m];
  internal_do_circ_embed(s->nn, s->m, s->cumm, s->halfm,
			 s->c, s->d, s->mtot, s->ce_dim,
			 &(s->FFT), false, res );
}

void do_tbmdirect(key_type *key, int m, Real *res ) {
  TBM_storage *s;
  s = (TBM_storage*) key->S[m];
  internal_do_directGauss(s->Direct, false, s->mtot, res);
}

void TBM_destruct(void **S) 
{
  if (*S!=NULL) {
   TBM_storage *x;
    x = *((TBM_storage**)S);
    if (x->x!=NULL) free(x->x);
    if (x->c!=NULL) free(x->c);
    if (x->d!=NULL) free(x->d);
    if (x->simuline!=NULL) free(x->simuline);
    if (x->Direct!=NULL) direct_destruct((void**) &(x->Direct));
    FFT_destruct(&(x->FFT));
    free(*S);
    *S = NULL;
  }
}

void SetParamTBMCE( int *action, int *force, Real *tolRe, Real *tolIm,
			int *trials, int *mmin, int *userfft, int *strategy) 
{
  SetParamCE(action, force, tolRe, tolIm, trials, mmin, userfft, strategy,
	     &TBMCE,"TBMCE");
}


void SetParamLines(int *action,int *nLines, Real *linesimufactor, 
		  Real *linesimustep, tbm_lines *tbm, int *every, char *name) {
  switch (*action) {
  case 0 :
    tbm->lines=*nLines;
    if ((*linesimufactor>0.0) ^ (*linesimustep>0.0)) {
      tbm->linesimufactor=*linesimufactor;
      tbm->linesimustep=*linesimustep;    
    } else { 
      if (*linesimufactor!=0.0) PRINTF("\nExactly one of `linesimufactor' and `linesimustep' should be positive.\n"); 
      // else, i.e. both are zero, the old values are kept!
    }
    tbm->every=*every;
    break;
  case 1 :
    *nLines=tbm->lines;
    *linesimufactor=tbm->linesimufactor;
    *linesimustep=tbm->linesimustep;  
    *every=tbm->every;
    if (GetNotPrint)  break;
  case 2:
   PRINTF("\n%s\n====\nLines=%d\nlinesimufactor=%f\nlinesimustep=%f\nevery=%d",
	   name,tbm->lines,tbm->linesimufactor,tbm->linesimustep,tbm->every);
    break;
  default : PRINTF(" unknown action\n"); 
  }
}

void SetParamTBM2(int *action,int *nLines, Real *linesimufactor, 
		  Real *linesimustep, int *every) {
  SetParamLines(action, nLines, linesimufactor, linesimustep, &tbm2, every,
		"TBM2");
}

void SetParamTBM3D2(int *action,int *nLines, Real *linesimufactor, 
		  Real *linesimustep, int *every) {
  SetParamLines(action, nLines, linesimufactor, linesimustep, &tbm3d2, every,
		"TBM3D2");
}

void SetParamTBM3D3(int *action,int *nLines, Real *linesimufactor, 
		  Real *linesimustep, int *every) {
  SetParamLines(action, nLines, linesimufactor, linesimustep, &tbm3d3, every,
		"TBM3D3");
}

void SetParamTBM(int *action, int *tbm_method) {
  switch (*action) {
  case 0 :
    if ((*tbm_method<=(int) Nothing) && (*tbm_method>=0)) {      
      SimulationType m;
      m=(SimulationType) *tbm_method;
      if ((m==CircEmbed) || (m==Direct) || (m==Nothing)){
	TBM_METHOD = m;
	return;
      }
    }
    if (GENERAL_PRINTLEVEL>0) {
      PRINTF("Specified method [%d] is not allowed. Method set to `Nothing'\n",
	     *tbm_method);
    }
    TBM_METHOD = Nothing;
    break;
  case 1 :
    *tbm_method = (int) TBM_METHOD;
    if (GetNotPrint)  break;
  case 2 :
      PRINTF("\nTBM:\n====\nuser defined method=%d\n",TBM_METHOD);
   break;
  default : PRINTF(" unknown action\n"); 
  }
}


// aufpassen auf ( * * 0 )
//               ( * * 0 )
//               ( 0 0 1 )

int init_turningbands(key_type *key, SimulationType method, int m)
{
  param_type param;
  int error, i, d, tbm_dim, v, start_aniso[MAXDIM], index_dim[MAXDIM];
  CovFctType CovFctTBM;
  TBM_storage *s;
  bool no_last_comp;
  Real dist, mindelta, aniso_last_save, 
    quotient[MAXCOV], *directx[2];
  unsigned short int actcov;
  int covnr[MAXCOV], multiply[MAXCOV], nonzero_pos;
  tbm_lines * tbm;
  SimulationType tbm_method;
  

  GENERAL_STORING = true;
  nonzero_pos = -1;
  SET_DESTRUCT(TBM_destruct);

  directx[0] = directx[1] = NULL;
  if ((key->S[m]=malloc(sizeof(TBM_storage)))==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (TBM_storage*)key->S[m];
  s->c = NULL;
  s->d = NULL;
  s->x = NULL;
  s->simuline=NULL;
  s->Direct=NULL;
  FFT_NULL(&(s->FFT));

  // get list of covariance functions that are of interest
  s->method = method;
    
  if (method==TBM2) {  CovFctTBM = CovFctTBM2;  tbm_dim = 2; }
  else {               CovFctTBM = CovFctTBM3;  tbm_dim = 3; }

 
  {    
    /****************************************************************/
    /*            Extraction of matching covariances                */
    /****************************************************************/
    // einfache checks
    // extraktion von covarianzfunktion mit gleichartiger anisotropie
    //    time-komponente wird nicht ge-checked da letzter parameter true
    // bereitstellung von param, quotient, multiply
 
    /* in tbm, time component must be ALWAYS alone, i.e. matrix looks like
     * * 0 
     * * 0
     * * x
   -- this is checked late on
   
   here x may be any constant. In case x is zero 
   (for all matrices (0..actcov-1), the last dimension
   vanishes (performed in the next step)
   
   the asterixes build submatrices for (v=1..actcov-1) that are multiple of the 
   submatrix for v=0, already checked by FIRST_CHECK_
    */

    if (GENERAL_PRINTLEVEL>4) PRINTF("\nextracting matching covariances...");
    Real store_param[TOTAL_PARAM];
    long bytes;
    bytes = sizeof(Real) * key->timespacedim * key->timespacedim;
    actcov=0;
    int v;
    for (v=0; v<key->ncov; v++) {
      if ((key->method[v]==method) && (key->left[v])) {
	//  && (key->param[v][VARIANCE]>0)) { -- 2.1.04
	key->left[v]=false;
	assert((key->covnr[v]>=0) && (key->covnr[v]<currentNrCov));
	assert(key->param[v][VARIANCE]>=0.0);
	covnr[actcov] = key->covnr[v];
	memcpy(param[actcov], key->param[v], sizeof(Real) * key->totalparam);
	if (actcov>0) {
	  if ((multiply[actcov-1] = key->op[v-1]) && key->method[v-1]!=method) {
	    if (GENERAL_PRINTLEVEL>0) 
	      PRINTF("severe error -- contact author, please");
	    error=ERRORMETHODMIX; goto ErrorHandling;
	  }
	}
	if (key->anisotropy) {
	  bool equal;
	  int w, endfor;
	  endfor = key->totalparam;
	  if (key->Time) endfor -= key->timespacedim; // note: time component
	  //   of the anisotropy matrix is not checked, which has the form
	  // (0 0 0 *)^T
	  Real f_epsilon = 1e-15;
	  if (actcov>0) {
	    /* check whether the multiplicative construction are consistent,*/
	    /* i.e. that the SPATIAL part of the anisotropy matrices are    */
	    /* multiplicatives of each other (more precisely, the remaining */
	    /* ones are multiplicatives of the first.)                      */
	    /* The time part of the matrix must be of the form (0,..,0,*).  */
	    /* A further goal of this part is the collection of additive    */
	    /* blocks that have the same anisotropy matrix structure, in    */
	    /* order to minimise the simulation time                        */

	    // note nonzero_pos is set in the else part of "actcov>0", first
	    // nonzero_pos gives the position of the first element in the 
            // first matrix of a multiplicative model that is not zero
	    assert(nonzero_pos > 0);
	    quotient[actcov] = 
	      param[actcov][nonzero_pos] / store_param[nonzero_pos];
	    equal = true;
	    for (w=ANISO; w<endfor; w++) {
	      equal &= fabs(store_param[w] * quotient[actcov] - param[actcov][w])
		< (fabs(param[actcov][w]) + (double)(param[actcov][w]==0.0)) *
		f_epsilon;
	      if (!equal) { /* even not equal up to a multiplicative constant */
		if (multiply[actcov-1]) { 
		  error=ERRORANISOMIX; goto ErrorHandling; 
		} else { /* **** nur additive ********* */
		  key->left[v]=true;
		  actcov--;
		  break;
		}
	      }
	    }
	  } else {
	    memcpy(&(store_param[ANISO]), &(param[actcov][ANISO]), bytes);
	    nonzero_pos=ANISO;
	    quotient[actcov] = 1.0;
	    // endfor: only the spatial components are considered
	    while ((nonzero_pos<endfor) && (param[actcov][nonzero_pos]==0))
	      nonzero_pos++;
	    if (nonzero_pos>=endfor) { error=ERRORTRIVIAL; goto ErrorHandling; }
	  }
	} else {
	  assert(fabs(param[v][SCALE] * param[v][INVSCALE]-1.0) < EPSILON);
	  quotient[actcov] = 1.0;
	}
	actcov++;
      } // v
    }
    if (actcov==0) { /* no covariance for the considered method found */
      if (key->traditional) error=ERRORNOTINITIALIZED;
      else error=NOERROR_ENDOFLIST;
      goto ErrorHandling;
    }
  }


  /****************************************************************/
  /*          investigation of the param structure                */
  /*          and the dimension w.r.t to TBM method               */
  /****************************************************************/
  // determine the reduced dimension of the space
  if (GENERAL_PRINTLEVEL>4) PRINTF("\nchecking parameter structure...");

  // in GetTrueDim only the first anisotropy matrix is investigated
  // so if there is a time component [key->totalparam-1], it 
  // is written temporarily into the first matrix
  aniso_last_save = param[0][key->totalparam-1]; 
  for (v=0; v<actcov; v++) 
    if (param[v][key->totalparam-1]!=0.0) {
      param[0][key->totalparam-1]=1.0;
      break; 
    }
  // no_last_comp: last dim of aniso matrix is identically zero
  GetTrueDim(key->anisotropy, key->timespacedim, 
	     param[0], // note param is modified, having as very last component
	     //           ==1 if any matrix has very last component <>0
	     &s->timespacedim,
	     &no_last_comp, // vanishing for *all* v=0..actcov-1 ?! 
	     start_aniso, index_dim);
  param[0][key->totalparam-1] = aniso_last_save;

  s->simuspatialdim = s->timespacedim;
  if (key->Time && !no_last_comp)  // note: no_last_comp means for all 0..actcov
    s->simuspatialdim--; 
  s->grid = key->grid && !key->anisotropy;
  if ((error=Transform2NoGrid(key, param[0], s->timespacedim,
			      start_aniso, &(s->x))) != NOERROR)
    goto ErrorHandling;



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
 /// ******** den nachfolgenden Teil  mit oben nach aniso_last_save 
/// vereinheitlichen. Notwendigkeit von GetTrueDim??
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

  // determine the *spatial* dimension of the turning band to be simulated
  // check correct dimension
 {
    bool TimeSpaceMix, twodim_ce;
    int w, v;
    
    TimeSpaceMix = false;
    if (twodim_ce = key->Time && !no_last_comp) {
      // note: no_last_comp means for all 0..actcov
      int start, endfor;
      start  = key->totalparam - key->timespacedim;
      endfor = key->totalparam-1;
      for (v=0; v<actcov; v++) {
	/* is it of the form   * * 0
	   .                   * * 0       ?
           .                   * * x 
	*/	
	for (w=start; w<endfor; w++) {
	  if (TimeSpaceMix = param[v][w]!=0.0) {
	    error = ERRORTIMESEPARATE;
	    goto ErrorHandling;
	  }
	}
      }
    }
    s->ce_dim = 1 + (int) (twodim_ce); // the dimension of the turning band
    if (s->timespacedim > twodim_ce + tbm_dim) { 
      error = ERRORWRONGDIM; goto ErrorHandling;
    }
  }


  /****************************************************************/
  /*          investigation of the covariance structure           */
  /****************************************************************/

  if (GENERAL_PRINTLEVEL>4) PRINTF("\nchecking covariance structure...");
  for (v=0; v<actcov; v++) {
    cov_fct *cov;
    cov = &(CovList[covnr[v]]);
    if (((method==TBM3) && (cov->ableitung==NULL)) ||
	((method==TBM2) && 
	 ((cov->cov_tbm2==NULL && cov->ableitung==NULL && !key->anisotropy) ||
	  (key->anisotropy && 
	   (((param[v][nonzero_pos]!=0.0) && cov->cov_tbm2==NULL && 
	     cov->ableitung==NULL) ||
	    ((param[v][nonzero_pos]==0.0) && cov->cov==NULL)
	    )
	    )
	  )
	  )
      )
      { 
	error=ERRORNOTDEFINED; goto ErrorHandling;}
    if ((!key->anisotropy && (CovList[covnr[v]].isotropic!=FULLISOTROPIC)) ||
	(cov->isotropic==ANISOTROPIC) ||
	((cov->isotropic==SPACEISOTROPIC) && (s->ce_dim==1))
      ) { error= ERRORISOTROPICMETHOD; goto ErrorHandling;}
    if ((key->Time) && no_last_comp && (cov->isotropic!=FULLISOTROPIC)) {
      error = ERRORWITHOUTTIME; goto ErrorHandling;} // braucht's net falls
    // GetDim und Trans2nogrid die dimensionen belassen!!
    if ((cov->check!=NULL) && 
	((error=cov->check(param[v], s->timespacedim, method))) != NOERROR)
      goto ErrorHandling;
  }

  if ((s->ce_dim==2) && (method==TBM2)) {
    if (GENERAL_PRINTLEVEL>4) PRINTF("...");
    // TBM2: per multiplication block at most 1 covariance function with
    //       non vanishing spatial components
    // the first matrix must have spatial components;
    // otherwise the the structure is considered as trivial
    int v,actcovM1;    
    // bool spatial_comp;
    actcovM1 = actcov -1; 
    for (v=0; v<actcovM1; ) {
      while ((v<actcovM1) && (!multiply[v])) v++;
      // spatial_comp = param[v][nonzero_pos] != 0.0;
      while ((v<actcovM1) && (multiply[v])) {
	if (param[v+1][nonzero_pos]!=0.0) {
	  error = ERRORTIMESEPARATE; goto ErrorHandling;
	  //if (spatial_comp) { error = ERRORTIMESEPARATE; goto ErrorHandling; }
	  //spatial_comp = true;
	}
	v++;
      }
    } 
  }
 

  /****************************************************************/
  /*          determine length and resolution of the band         */
  /****************************************************************/
 // sort out which finer grid should be used in the spatial domain for
  // low dimensional simulation
  if (GENERAL_PRINTLEVEL>4) PRINTF("\ndetermination of the band length...");
  switch (s->simuspatialdim) {
  case 1 : 
    // for simplicity take for the 1-dim case just the 2-dim parameters
    if (method==TBM2) tbm = &tbm2;
    else { assert(method==TBM3); tbm = &tbm3d2; }
    break;
  case 2 : 
    if (method==TBM2) tbm = &tbm2;
    else { assert(method==TBM3); tbm = &tbm3d2; }
    break;
  case 3 :
    if (method!=TBM3){assert(method==TBM2);error=ERRORFAILED;goto ErrorHandling;}
    tbm = &tbm3d3;
   break;
  default: 
    error=ERRORUNSPECIFIED; goto ErrorHandling;
  }

  // minimum and maximum distance between the points   
  // and set the scale parameters for the band simulation
  if (s->grid) {
    // then, for sure, we have an isotropic definition!
    assert(CovList[covnr[0]].isotropic==FULLISOTROPIC);
    {
      // diameter of the grid
      register Real dummy;
      dist = 0.0; // max distance
      for (d=0; d < s->simuspatialdim; d++) {
	dummy = key->x[d][XSTEP] * (Real) (key->length[d] - 1);
	dist += dummy * dummy;
      }
    }

    // if linesimustep, the fine grid scale is proportional to smallest 
    // distance in the grid. 

    if (tbm->linesimustep>0.0) { s->linesimuscale =  1/tbm->linesimustep; }
    else {
      mindelta = RF_INF;
      for (d=0; d<s->simuspatialdim; d++) {
	if ((key->x[d][XSTEP]<mindelta) && (key->x[d][XSTEP]>0)) 
	  {mindelta=key->x[d][XSTEP];}
      }
      s->linesimuscale = tbm->linesimufactor/mindelta;
    }
    // set the parameters for the line simulation accordingly
    for (v=0; v<actcov; v++) {
      param[v][ANISO] = 1.0 / (s->linesimuscale * param[v][SCALE]);
    }
  } else { // not s->grid
    // s->timespacedim, s->x, s->spatialdim,
    Real min[MAXDIM],max[MAXDIM], dummy; 
    int j, d, ix, jx, endfor;
    for (d=0; d<MAXDIM; d++) {min[d]=RF_INF; max[d]=RF_NEGINF;}
    if (key->grid) { // key->grid     
      Real sxx[ZWEIHOCHMAXDIM * MAXDIM];
      // unsorted, reduced for param[0...0], #=2^Dim, 
      endfor = 1 << key->timespacedim;
      
      // to determine the diameter of the grid determine first 
      // componentwise min and max corner
      GetCornersOfGrid(key, s->timespacedim, start_aniso, param[0], sxx);      
      for (ix=i=0; i<endfor; i++, ix+=s->timespacedim) { 
	for(d=0; d<s->timespacedim; d++) {
	  if (sxx[ix+d]<min[d]) min[d] = sxx[ix+d];
	  if (sxx[ix+d]>max[d]) max[d] = sxx[ix+d];
	}
      }
      
      // determine smallest distance
      GetCornersOfElement(key, s->timespacedim, start_aniso, param[0], sxx);
      GetRangeCornerDistances(key, sxx, s->timespacedim, s->simuspatialdim,
				&mindelta, &dummy);
    } else { // not key->grid
      mindelta = RF_INF; 
      if (key->totalpoints >50000 && tbm->linesimufactor!=0.0) {
	 error=ERRORTOOMANYPOINTS; goto ErrorHandling;
	 /* algorithmus kann verbessert werden, so dass diese Fehlermeldung nicht
            mehr notwendig ist! */ 
      }
      for (ix=i=0; i<key->totalpoints; i++, ix+=s->timespacedim) {
	if ((GENERAL_PRINTLEVEL>4) && ((i % 10000)==0))
	  PRINTF(" %d [%d]\n",i,key->totalpoints);
	// determine componentwise min and max (for the diameter)
	for(d=0; d<s->simuspatialdim; d++){//temporal part need not be considered
	  if (s->x[ix+d]<min[d]) min[d] = s->x[ix+d];
	  if (s->x[ix+d]>max[d]) max[d] = s->x[ix+d];
	}
	if (tbm->linesimustep==0.0) {
	  // determine smallest distance ((may take a lot of time!!))
	  for (jx=j=0; j<i; j++, jx+=s->timespacedim) {
	    register Real diff,dist;
	    for(dist=0.0, d=0; d<s->timespacedim; d++) {
	      diff = s->x[ix+d] - s->x[jx+d]; 
	      dist += diff * diff;
	    }
	    if ((dist>0) && (dist<mindelta)) mindelta=dist; 
	  }
	} // if linesimustep==0.0
      }
      mindelta = sqrt(mindelta);
    }
    // max distance, upperbound
    dist = 0.0;
    for (d=0; d<s->simuspatialdim; d++) {
      dist += (max[d]-min[d])*(max[d]-min[d]);
      s->half[d] = 0.5 * (max[d]+min[d]); // s->half only used if not grid
    }
    for (; d<MAXDIM; d++) s->half[d] = RF_NAN;
    if (tbm->linesimustep==0.0) s->linesimuscale = tbm->linesimufactor/mindelta;
    else s->linesimuscale =  1/tbm->linesimustep;

    // put the SCALE parameters correctly with respect to
    // and the pure time component of the regular case
    for (v=0; v<actcov; v++) param[v][ANISO] = quotient[v] / s->linesimuscale;
    if (s->ce_dim==2) 
      for (v=0; v<actcov; v++) {
	param[v][ANISO+1] = param[v][ANISO+2] = 0;
	param[v][ANISOP3] = param[v][key->totalparam - 1] * key->T[XSTEP];
      }

    if (false) 
    {
      int i,j,k;
      for (k=i=0; i<key->totalpoints; i++) {
	for(j=0; j<s->timespacedim; j++, k++) {
	  printf("%f ", s->x[k]);
	}
	printf("\n");
      }
      for (d=0; d<s->simuspatialdim; d++) 
	printf("min=%f, max=%f\n", min[d], max[d]);
    }

  } // !s->grid

 
  // number of points to be simulated on the "line"
  { 
    int i;
    // nn=number of pixels to be simulated on the lines [==(length of diagonal 
    // through grid that is to be simulated) * (linesimufactor=="precision")]
    // 5 is for safety
    Real dummy;
    dummy = (sqrt(dist) * s->linesimuscale);
    if (dummy > MAXNN) {error=ERRORNN; goto ErrorHandling;}
    s->nn[0]= 5 + (long) dummy; // +5 for safety
    if (s->ce_dim==2) {
      s->nn[1] = key->length[key->timespacedim - 1];
      
      { 
	long t, idx, endtimeindex, spatialpoints;
	double tdouble;
	endtimeindex  =  key->length[key->spatialdim] * s->nn[0];
	spatialpoints = key->totalpoints / key->length[key->spatialdim];
	// s->simuspatialdim: very first time component, +s->timespacedim
	// the following ones
	for (t=0, idx = s->simuspatialdim; t<endtimeindex; t+=s->nn[0]) {
	  tdouble = (double) t;
	  for (i = 0; i < spatialpoints;  i++, idx+= s->timespacedim)  
	    s->x[idx] = tdouble;
	}
      }
    } else s->nn[1]=1; /* safety */
    s->mtot = 1;
    for (i=0; i<s->ce_dim; i++)  s->mtot *= s->nn[i]; 
  }
  s->nnhalf = s->nn[0] / 2;



  tbm_method = key->tbm_method;
  if (s->mtot > DIRECTGAUSS_MAXVARIABLES) {
    tbm_method=CircEmbed;
    if (GENERAL_PRINTLEVEL>1) 
      PRINTF("too many points -- using circulant embedding\n");
  }
  switch (tbm_method) {
      case Direct :
	if (GENERAL_PRINTLEVEL>4)
	  PRINTF("\ninitializing matrix decomposition...");
	s->linemethod = do_tbmdirect;
	for (d=0; d<2; d++) {
	  if ((directx[d]= (Real *) malloc(sizeof(Real) * 3))==0) {
	    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
	  }
	  directx[d][0]=1.0; 
	  directx[d][1]= (double) s->nn[d]; 
	  directx[d][2]=1.0;
	}
	//printkey(key);
	error=internal_init_directGauss(&(s->Direct),
					true, //grid
					s->ce_dim,
					false, // no time component
					directx,
					s->nn,
					s->mtot,
					NULL, // time vector
					CovFctTBM,
					covnr, multiply, param, 
					actcov,
					s->ce_dim>=2 // anisotropy
	  );
	if (error!=NOERROR) goto ErrorHandling;
        // waere auch ein durchschleusen zu circembed moeglich
	// koennte dann aber probleme schaffen, falls weitere Methoden
	// hinzukommen.
	break;
      case CircEmbed : 
	long twoRealmtot;
	if (GENERAL_PRINTLEVEL>4)
	  PRINTF("\ninitializing circulant embedding...");
	s->linemethod = do_tbmcirculantembedding;  
	error=internal_init_circ_embed(UNIT,      // step length ( = 1 )
				       s->ce_dim>=2, // anisotropy
				       covnr, multiply, param, // model
				       s->nn,     // size of the grid
				       s->m, s->cumm, s->halfm,//return variable
				       s->ce_dim, // dimension
				       actcov,    // # of cov.fcts
				       CovFctTBM, // Cov.Fct
				       &TBMCE,    //technical parameters
				       &(s->FFT), // FFT storage
				       &twoRealmtot,//return variable: 
				       //            size of FFT vector
				       &(s->c)); //return variable; FFT transf
	if (error!=NOERROR) goto ErrorHandling; 
	if ((s->d=(double *)malloc(twoRealmtot))==0){
	  error=ERRORMEMORYALLOCATION;goto ErrorHandling;}
	break;
      default : error=ERRORNOTPROGRAMMED; goto ErrorHandling;
  }
  if ((s->simuline=(Real *)malloc(sizeof(Real) * s->mtot))==0){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;}

  if (GENERAL_PRINTLEVEL>4) PRINTF("\ntbm is now initialized.\n");
 
  if (directx[0]!=NULL) free(directx[0]);
  if (directx[1]!=NULL) free(directx[1]);

  // im folgenden anisotropy und nicht traditional, da 
  // oben alle Kovarianzfunktionen zusammengefasst werden mit derselben
  // (An)isotropie-Struktur, hei3t nur bei echter Anisotropie wird TBM2
  // etc. mehrfach aufgerufen!
  if (key->anisotropy) return NOERROR_REPEAT;
  else return NOERROR;
  
 ErrorHandling:
  if (directx[0]!=NULL) free(directx[0]);
  if (directx[1]!=NULL) free(directx[1]);
  return error;
}

int init_turningbands2(key_type *key, int m) {
 return(init_turningbands(key, TBM2, m));
}

int init_turningbands3(key_type *key, int m) {
  return(init_turningbands(key, TBM3, m));
}

void do_turningbands(key_type *key, int m, Real *res )
{ 
  Real phi;
  Real *simuline;
  long n,nn, nnhalf; 
  simu1dim simulate;
  assert(key->active);
  tbm_lines tbm;

  TBM_storage *s;
  s = (TBM_storage*)key->S[m];

  simulate = s->linemethod;
  nn = s->nn[0];
  nnhalf = s->nnhalf; 
  simuline = s->simuline; 
  
  for (n=0; n<key->totalpoints; n++) res[n]=0.0; 
  switch ((int) (s->method==TBM3) + (int) (key->spatialdim>2)) {
      case 0 : tbm = tbm2;
	break;
      case 1 : tbm = tbm3d2;
	break;
      case 2 : tbm = tbm3d3;
	break;
      default: assert(false);
  }
    
  if (GENERAL_PRINTLEVEL>=7) {
    PRINTF("%s : grid=%d % key.spatialdim=%d s.simuspatialdim=%d\n",
	   METHODNAMES[s->method], s->grid, key->spatialdim, s->simuspatialdim);
    if (GENERAL_PRINTLEVEL>=8) printkey(key);
  }
  
  if (s->grid) { // old form, isotropic field
    Real deltax, yoffset,deltay, nxhalf, nyhalf,xoffset,
      linesimuscaleX,linesimuscaleY,linesimuscaleZ;  
    int nx,zaehler,ny;
    nxhalf=((Real)(key->length[0]-1))/2.0; ////
    nyhalf=((Real)(key->length[1]-1))/2.0; ////
    linesimuscaleX = s->linesimuscale * key->x[0][XSTEP]; ////
    if (s->method==TBM2) {// isotropy, TBM2, grid 
      Real deltaphi;
      deltaphi = PI / (Real) tbm.lines;
      phi = deltaphi * UNIFORM_RANDOM; 
      if (key->spatialdim==1) { // dim =1
	for (n=0;n<tbm.lines;n++){
	  if (tbm.every>0 && (n % tbm.every == 0)) PRINTF("%d \n",n);
	  deltax=sin(phi) * linesimuscaleX;
	  simulate(key, m, simuline );
	  xoffset= (double)nnhalf - nxhalf * deltax;
	  zaehler = 0;
	  for(nx=0; nx<key->length[0];  nx++) {
	    assert((xoffset<nn) && (xoffset>=0) );
	    res[zaehler++] += simuline[(long) xoffset];
	    xoffset += deltax;
	  }
	  phi += deltaphi;
	}
      } else { // dim == 2
	linesimuscaleY = s->linesimuscale * key->x[1][XSTEP]; 
	for (n=0;n<tbm.lines;n++){
	  if (tbm.every>0  && (n % tbm.every == 0)) PRINTF("%d \n",n);
	  deltax=sin(phi) * linesimuscaleX;
	  deltay=cos(phi) * linesimuscaleY;
	  simulate(key, m, simuline );
	  /* depending on deltax and deltay,
	     shorter lines might be 
	     simulated! -> improvement in 
	     speed?! */
	  /* then nnhalf must be kept variable, too: */
	  yoffset= (double)nnhalf - nyhalf * deltay - nxhalf * deltax;
	  zaehler = 0;
	  for (ny=0;ny<key->length[1]; ny++) {	  
	    xoffset = yoffset;
	    for(nx=0;nx<key->length[0];  nx++) {
	      assert((xoffset<nn) && (xoffset>=0) );
	      res[zaehler++] += simuline[(long) xoffset];
	      xoffset += deltax;
	    }
	    yoffset += deltay;
	  }
	  phi += deltaphi;
	}
      }
    } else {// isotropy, TBM3, grid
      assert(s->method==TBM3);
      switch (key->spatialdim) {
	  case 1 :
	    for (n=0;n<tbm.lines;n++){
	      if (tbm.every>0  && (n % tbm.every == 0)) PRINTF("%d \n",n);
	      deltax= UNIFORM_RANDOM * linesimuscaleX;
	      simulate(key, m, simuline ); 
	      xoffset=  (double)nnhalf - nxhalf * deltax ;
	      zaehler = 0;
	      for(nx=0;nx<key->length[0];nx++) {
		assert((xoffset<nn) && (xoffset>=0) );
		res[zaehler++] += simuline[(long) xoffset];	  
		if (!(fabs(simuline[(long) xoffset])<10000.0)) 
		  PRINTF("s:%f ",simuline[(long) xoffset]);
		xoffset += deltax;
	      }
	    }
	    break;
	  case 2:
	    linesimuscaleY = s->linesimuscale * key->x[1][XSTEP]; ////
	    for (n=0;n<tbm.lines;n++){
	      if (tbm.every>0  && (n % tbm.every == 0)) 
		PRINTF("%d \n",n);
	      deltax= UNIFORM_RANDOM;// see Martin's tech rep for details
	      deltay= sqrt(1.0-deltax*deltax) * sin(UNIFORM_RANDOM*TWOPI) * 
		linesimuscaleY;
	      deltax *= linesimuscaleX;
	      simulate(key, m, simuline ); 
	      /* depending on deltax and deltay, shorter lines might be 
		 simulated! -> improvement in speed!
	      */
	      yoffset=  (double)nnhalf - nyhalf * deltay - nxhalf * deltax ;
	      zaehler = 0;
	      for (ny=0;ny<key->length[1];ny++) {
		xoffset = yoffset;
		for(nx=0;nx<key->length[0];nx++) {
		  assert((xoffset<nn) && (xoffset>=0) );
		  res[zaehler++] += simuline[(long) xoffset];	  
		  if (!(fabs(simuline[(long) xoffset])<10000.0)) 
		    PRINTF("s:%f ",simuline[(long) xoffset]);
		  xoffset += deltax;
		}
		yoffset += deltay;
	      }
	    }  
	    break;
	  case 3:
	   
	    Real dummy,nzhalf,zoffset,deltaz;
	    int nz;
	    linesimuscaleY = s->linesimuscale * key->x[1][XSTEP]; ////
	    linesimuscaleZ= s->linesimuscale * key->x[2][XSTEP];
	    nzhalf=((Real)(key->length[2]-1))/2.0;
	    for (n=0; n<tbm.lines; n++){
	      if (tbm.every>0  && (n % tbm.every == 0)) PRINTF("%d \n",n);
	      //
	      deltaz = 2.0 * UNIFORM_RANDOM - 1.0;
	      dummy = sqrt(1.0 - deltaz * deltaz);
	      deltaz *= linesimuscaleZ;
	      //
	      deltay= UNIFORM_RANDOM * TWOPI;
	      deltax= cos(deltay) * dummy * linesimuscaleX;
	      deltay= sin(deltay) * dummy * linesimuscaleY;
	      simulate(key, m, simuline );
	      zoffset= (double) nnhalf - nzhalf * deltaz
		- nyhalf * deltay - nxhalf * deltax; 
	      zaehler = 0;    
	      for (nz=0;nz<key->length[2];nz++) {
		yoffset = zoffset;
		for (ny=0;ny<key->length[1];ny++) {
		  xoffset = yoffset;
		  for(nx=0;nx<key->length[0];nx++) {
		    ///
		    if ((xoffset>=nn) || (xoffset<0))  
		      PRINTF(" %d %d %d | %d %d %d \n %f %f %f\n %d %f %f %f\n %d %f\n",
			     nx,ny,nz,key->length[0],key->length[1],
			     key->length[2],
			     deltax,deltay,deltaz,
			     nnhalf,nxhalf,nyhalf,nzhalf,
			     nn,xoffset);
		    ///
		    assert( (xoffset<nn) && (xoffset>=0) );
		    res[zaehler++] += simuline[(long) xoffset];	  
		    if (!(fabs(simuline[(long) xoffset])<10000.0)) 
		      PRINTF("s:%f ",simuline[(long) xoffset]);
		    xoffset += deltax;
		  }
		  yoffset += deltay;
		}
		zoffset += deltaz;
	      }
	    } 
	    break;
	  default : assert(false);
      }
    }
  } else { 
    // not grid, could be time-model!
    // both old and new form included
    //offset=0; for (i=0; i<s->mtot; i++) offset += simuline[i] * 
    //                simuline[i];printf("%f ", offset / (double) s->mtot);

#define TBMLOOP(UNITYVECTOR, OFFSET, INDEX)\
  {\
    for (n=0; n<tbm.lines; n++){\
      if (tbm.every>0  && (n % tbm.every == 0)) PRINTF("%d \n",n); \
      UNITYVECTOR\
      simulate(key, m, simuline); \
      offset=(Real) nnhalf - (OFFSET); \
      for (v = 0, i = 0; i < key->totalpoints; i++) {\
	register long index;\
	index = (long) (offset + INDEX);\
	if ((index>=s->mtot) || (index<0)) {\
	  PRINTF(" index=%d n0=%d n1=%d mtot=%d index=%f; offset=%f x=%f ex=%f t=%f scale=%f half=%f nnhalf=%d v=%d, i=%d\n", \
		 index, s->nn[0], s->nn[1], s->mtot, offset + s->x[v] * ex,\
		 offset, s->x[v], ex, s->x[v+1], linescale, halfx, nnhalf,\
                 v, i);\
	  assert((index<s->mtot) &&  (index>=0));  \
	}\
	res[i] += simuline[index];\
	v += s->timespacedim;\
      }\
    }\
  }

#define TBMST(UNITYVECTOR, OFFSET, INDEX, TIMEINDEX)\
  if (s->ce_dim==1) TBMLOOP(UNITYVECTOR, OFFSET, INDEX)\
  else {assert(s->ce_dim==2); TBMLOOP(UNITYVECTOR, OFFSET, INDEX + TIMEINDEX)};


    Real halfx,halfy,ex,ey,offset,linescale;
    long v;

    halfx = s->half[0];
    halfy = s->half[1];
    linescale = s->linesimuscale;
    int i;
    if (s->method==TBM2) {
      Real deltaphi;
      deltaphi = PI / (Real) tbm.lines;
      phi = deltaphi * UNIFORM_RANDOM; 
      if (s->simuspatialdim==1) { // dim == 1, TBM 2, arbitrary locations
	TBMST({phi += deltaphi;ex=sin(phi) * linescale;}, //UNITYVECTOR
	      halfx * ex,    //OFFSET
	      s->x[v] * ex,  //INDEX
	      s->x[v+1]);    //TIMEINDEX
      } else {  // dim == 2, TBM 2
	TBMST({phi += deltaphi;ex=sin(phi)*linescale; ey=cos(phi)*linescale;},
	      halfx * ex + halfy * ey,
	      s->x[v] * ex + s->x[v+1] * ey,
	      s->x[v+2]);
      }
    } else { // TBM3, not grid, dim=1 or 2
      Real linesimusq, dummy, ez, halfz;
      assert(s->method==TBM3);
      linesimusq = linescale * linescale;
      switch (s->simuspatialdim) {
	  case 1 : //see Martin's techrep f. details
	    TBMST({ex = UNIFORM_RANDOM * linescale;},
		  halfx * ex, 
		  s->x[v] * ex,
		  s->x[v+1]);
	    break;
	  case 2: 
	    TBMST({ez = UNIFORM_RANDOM * linescale;
	           dummy = sqrt(linesimusq - ez * ez);
		   ey= UNIFORM_RANDOM * TWOPI;
		   ex= cos(ey) * dummy;
		   ey= sin(ey) * dummy;},
		  halfx * ex + halfy * ey,
		  s->x[v] * ex + s->x[v+1] * ey,
		  s->x[v+2]);
	    break;
	  case 3:
	    halfz = s->half[2];
	    TBMST({ez = UNIFORM_RANDOM * linescale;
	           dummy = sqrt(linesimusq - ez * ez);
		   ey = UNIFORM_RANDOM * TWOPI;
		   ex = cos(ey) * dummy;
		   ey = sin(ey) * dummy;},
		  halfx * ex + halfy * ey + halfz * ez,
		  s->x[v] * ex + s->x[v+1] * ey + s->x[v+2] * ez,
		  s->x[v+3]);
	    break;
	  default : assert(false);
      }
    }
  }
  {
    register long i;
    register Real InvSqrtNlines;   
    InvSqrtNlines=1.0 / sqrt((Real) tbm.lines);
    for(i=0;i<key->totalpoints;i++) { 
      res[i] *= InvSqrtNlines; 
    }
  }
}


void pokeTBM(int *out, int *in, int *err) {
  int m;
  key_type *keyOut, *keyIn;
  TBM_storage *sOut, *sIn;
  
  *err = 0;
  if ((*out<0) || (*out>=MAXKEYS) || (*in<0) || (*in>=MAXKEYS))
    {*err=ERRORKEYNROUTOFRANGE; goto ErrorHandling;}
  
  keyOut = &(KEY[*out]);
  keyIn =  &(KEY[*in]);
  if (!keyIn->active || !keyOut->active)
    {*err=ERRORNOTINITIALIZED; goto ErrorHandling;}
  if (keyIn->spatialdim != keyOut->spatialdim)
    {*err=ERRORPOKETBMdim; goto ErrorHandling;}
  if (keyIn->n_unimeth != keyOut->n_unimeth)
    {*err=ERRORPOKETBMmeth; goto ErrorHandling;}
  
  for (m=0; m<keyIn->n_unimeth; m++) {
    if (keyIn->unimeth[m] != keyOut->unimeth[m])
      {*err=ERRORPOKETBMmeth; goto ErrorHandling;}
    if ((keyOut->unimeth[m]==TBM2) || (keyOut->unimeth[m]==TBM3)) {
      bool grid;
      sOut = (TBM_storage*) keyOut->S[m];
      sIn = (TBM_storage*) keyIn->S[m];
      if ((sOut->simuspatialdim!=sIn->simuspatialdim) ||
	  (sOut->timespacedim != sIn->timespacedim) ||
	  (sOut->ce_dim != sIn->ce_dim))
	{*err=ERRORPOKETBMdim; goto ErrorHandling;} 
      if (sIn->linesimuscale != sOut->linesimuscale)
	{*err=ERRORPOKETBMparam; goto ErrorHandling;}
      // other parameters should also be tested extensively!!
      grid = sIn->grid;
      TBM_destruct(&(keyIn->S[m]));
      keyIn->S[m] = keyOut->S[m];
      keyOut->S[m]= NULL;
      sIn = (TBM_storage*) keyIn->S[m];
      sIn->grid = grid;
     }
  }
  DeleteKey(out);
  return;
  
 ErrorHandling:
  if (GENERAL_PRINTLEVEL>3) PRINTF("error %d %d..\n", *in, *out);
  DeleteKey(in);
  DeleteKey(out);
  return;
}







