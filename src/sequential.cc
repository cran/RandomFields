/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by sequential method

 Copyright (C) 2001 -- 2014 Martin Schlather, 

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

#include <math.h>
#include <stdio.h>  
#include <stdlib.h>
#include <R_ext/Lapack.h>

#include "RF.h"
#include "shape_processes.h"
#include "Coordinate_systems.h"
//#include <R_ext/Linpack.h>

#define SEQU_MAX (COMMON_GAUSS + 1)
#define SEQU_BACK (COMMON_GAUSS + 2)
#define SEQU_INIT (COMMON_GAUSS + 3)


bool debugging = true;


int check_sequential(cov_model *cov) {
#define nsel 4
  cov_model *next=cov->sub[0];
  int err,
    dim = cov->tsdim; // taken[MAX DIM],
  sequ_param *gp  = &(GLOBAL.sequ);
  location_type *loc = Loc(cov);

  ROLE_ASSERT(ROLE_GAUSS);
  if (!loc->grid && !loc->Time) 
    SERR1("'%s' only possible if at least one direction is a grid", NICK(cov));

  kdefault(cov, SEQU_MAX, gp->max);
  kdefault(cov, SEQU_BACK, gp->back);
  kdefault(cov, SEQU_INIT, gp->initial);
  if ((err = checkkappas(cov)) != NOERROR) return err;
 
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;

  if ((err = CHECK(next, dim, dim, PosDefType, XONLY, SYMMETRIC,
		     SUBMODEL_DEP, ROLE_COV)) != NOERROR) return err;
  if (next->pref[Sequential] == PREF_NONE) return ERRORPREFNONE;
  setbackward(cov, next);

  return NOERROR;
}


void range_sequential(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range) {

  range->min[SEQU_MAX] = 0;
  range->max[SEQU_MAX] = RF_INF;
  range->pmin[SEQU_MAX] = 100;
  range->pmax[SEQU_MAX] = 8000;
  range->openmin[SEQU_MAX] = false;
  range->openmax[SEQU_MAX] = true; 

  range->min[SEQU_BACK] = 0;
  range->max[SEQU_BACK] = RF_INF;
  range->pmin[SEQU_BACK] = 0.1;
  range->pmax[SEQU_BACK] = 10;
  range->openmin[SEQU_BACK] = true;
  range->openmax[SEQU_BACK] = true; 

  range->min[SEQU_INIT] = RF_NEGINF;
  range->max[SEQU_INIT] = RF_INF;
  range->pmin[SEQU_INIT] = -10;
  range->pmax[SEQU_INIT] = 10;
  range->openmin[SEQU_INIT] = false;
  range->openmax[SEQU_INIT] = true; 
}

// (U1, U2) ~ N(0, S) =>
// U1 | U2 ~ N(S_{12} S_{22}^{-1} u2, S_11 - S12 S_22^{-1} S_21)

// start mit S_22; dann nur zeilenweise + Einschwingen

int init_sequential(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  cov_model *next = cov->sub[0];
  // cov_fct *C=CovList + next->gatternr;
  // covfct cf = C->cov;
  location_type *loc = Loc(cov);
  if (loc->distances) return ERRORFAILED;

  int withoutlast, d, endfor, l, 
    err = NOERROR,
    dim = cov->tsdim,
    spatialdim = dim - 1,
    vdim = next->vdim[0],
    max = P0INT(SEQU_MAX),
    back= P0INT(SEQU_BACK), 
    initial= P0INT(SEQU_INIT);
  assert(dim == loc->timespacedim);
  assert(next->vdim[0] == next->vdim[1]);
  if (initial < 0) initial = back - initial;

  double 
    //*caniso = loc->caniso,
    *timecomp = loc->grid ? loc->xgr[spatialdim] : loc->T,
    *xx = NULL,
    *G = NULL, 
    *U11 = NULL, 
    *COV21 = NULL,
    *MuT = NULL, 
    *U22 = NULL, 
    *Inv22 = NULL;
  res_type
    *res0 = NULL;
  sequ_storage* S = NULL;
  long  i, 
    timelength = loc->grid ? loc->xgr[spatialdim][XLENGTH] :loc->T[XLENGTH],
    spatialpnts = loc->totalpoints / timelength,
    totpnts = back * spatialpnts, 
    totpntsSQ =  totpnts * totpnts,
    spatialpntsSQ = spatialpnts * spatialpnts,
    spatialpntsSQback = totpnts * spatialpnts;

  bool 
    storing = GLOBAL.internal.stored_init,
    Time = loc->Time;
 
  if (cov->role == ROLE_COV) {
    return NOERROR;
  }
  
  ROLE_ASSERT_GAUSS;
  cov->method = Sequential;
     
  if (!loc->grid && !Time) 
    GERR("last component must be truely given by a non-trivial grid");

  if (CovList[next->nr].implemented[Sequential] != IMPLEMENTED) {
    err=ERRORNOTDEFINED; 
    goto ErrorHandling;
  }
  
  if (cov->vdim[0] > 1) {
      err=ERRORNOMULTIVARIATE; 
      goto ErrorHandling;   
  }
  
  if (totpnts > max) {
    sprintf(ERRORSTRING_OK, "number of points less than '%s' (=%d)", 
	    KNAME(SEQU_MAX), max);
    sprintf(ERRORSTRING_WRONG,"%d * %ld = %ld", back, spatialpnts, totpnts);
    err=ERRORCOVFAILED; goto ErrorHandling;
  }

 if (timelength <= back) {
    GERR2("the grid in the last direction is too small; use method '%s' instead of '%s'",
	  CovList[DIRECT].nick, CovList[SEQUENTIAL].nick);
  } 
  if (back < 1) back = max / spatialpnts;

  if ((U22 = (double *) MALLOC(sizeof(double) * totpntsSQ))==NULL ||
      (Inv22 = (double *) MALLOC(sizeof(double) * totpntsSQ))==NULL ||
      (COV21 = (double *) MALLOC(sizeof(double) * spatialpntsSQback))==NULL ||
      (U11 = (double *) MALLOC(sizeof(double) * spatialpntsSQ))==NULL ||
      (MuT = (double *) MALLOC(sizeof(double) * spatialpntsSQback))==NULL ||
      (G = (double *) MALLOC(sizeof(double) * totpnts))==NULL ||
      (res0 = (res_type *) MALLOC(sizeof(res_type) * vdim *
				  (totpnts + spatialpnts * initial))) ==NULL) {
    err=ERRORMEMORYALLOCATION;  
    goto ErrorHandling;
  }
  NEW_STORAGE(sequ);
  S = cov->Ssequ;

 
  /* ************************* */
  /* create matrix of explicit */
  /*       x-coordinates       */
  /* ************************* */

 
  if (loc->grid) loc->xgr[spatialdim][XLENGTH] = back; 
  else loc->T[XLENGTH] = back;
  Transform2NoGrid(cov, &xx);
  assert(loc->caniso == NULL);
  if (loc->grid) loc->xgr[spatialdim][XLENGTH]=timelength; 
  else loc->T[XLENGTH] = timelength;
 
  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  long j, k, k0, k1, k2, segment;
  int  row, f77err;
  double *y;
  y = (double*) MALLOC(dim * sizeof(double));

  // *** S22
  if (PL==PL_SUBDETAILS) { LPRINT("covariance matrix...\n"); }
  k = 0;
  for (k0 =i=0; i<totpnts; i++, k0+=dim) {
    k += i;
    for (segment=i* totpnts+i, j=i, k2=k0; j<totpnts; j++) {
      k1 = k0;
      for (d=0; d<dim; d++) {
	y[d] = xx[k1++] - xx[k2++];
      }
      COV(y, next, U22 + segment);
      U22[k++] = U22[segment];
      segment += totpnts;
    }
  }

  /*
  for (k=i=0; i<totpnts; i++) {
    for (d=0; d<totpnts; d++) {
      printf("%f ", U22[k++]); //
    }
    printf("\n"); //
  }
  */
       
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
  y[spatialdim] = timecomp[XSTEP] * back;
  for (k0 =i=0; i<spatialpnts; i++, k0+=dim) { // t_{n+1}
      for (k2=j=0; j<spatialpnts; j++) { // t_1
	k1 = k0;
	for (d=0; d<spatialdim; d++) {
	    y[d] = xx[k1++] - xx[k2++];
	}
	COV(y, next, COV21 + k);	
	k++; k2++;
    }
    k += withoutlast;
  }

  FREE(y);

//  for (i=0; i<withoutlast * spatialpnts; i++) {
//    print("%3.4f ", COV21[i]);
//  }
//  assert(false);

/*
  LPRINT("U11\n");
  for (i=0; i<spatialpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
    }
    LPRINT("\n");
  }

  LPRINT("S22 %d %d\n", totpnts,spatialpnts );
  for (i=0; i<totpnts; i++) {
    LPRINT("%2d: ", i);
    for (j=0; j<totpnts; j++) {
      LPRINT("%3.2f ", U22[i + j * totpnts]);
   }
    LPRINT("\n");
  }
*/


  /* ********************* */
  /* matrix compositions  */
  /* ********************* */

  // *** sqrt S 22 
  if (PL>=PL_STRUCTURE) { LPRINT("Cholesky decomposition of Cov22...\n"); }
  row=totpnts;
  // dchdc destroys the input matrix; upper half of U22 contains result!
  // dpotrf	F77_CALL(dpotrf)("Upper", &m, REAL(ans), &m, &i);
  //  assert(false);
  F77_CALL(dpotrf)("Upper", &row, U22, &row, &f77err);
//  assert(false);
  // F77_NAME(dchdc)(U22, &row, &row, G, NULL, &choljob, &f77err);
  if (f77err!=NOERROR) {
    if (PL>=PL_SUBIMPORTANT) 
      {INDENT; PRINTF("Error code F77_CALL(dpotrf) = %d\n", f77err);}
    err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  } 

  // for (i=0; i<withoutlast * spatialpnts; i++) {
//    LPRINT("%3.4f ", COV21[i]);
//  }
//  assert(false);
  
  // *** inverse of S 22 
  if (PL>=PL_STRUCTURE) { LPRINT("inverse of Cov22...\n"); }
  MEMCOPY(Inv22, U22, sizeof(double) * totpntsSQ);
  F77_CALL(dpotri)("Upper", &row, Inv22, &row, &f77err);
  if (f77err!=NOERROR) {
    if (PL>=PL_SUBIMPORTANT) { INDENT; PRINTF("Error code F77_CALL(dpotri) = %d\n", f77err); }
    err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  }
  for (k=i=0; i<totpnts; i++) {
    k += i;
    for (segment=i* totpnts+i; segment<totpntsSQ; segment += totpnts) {
	Inv22[k++] = Inv22[segment]; 
//	LPRINT("sg %d %f %d \n", k, Inv22[segment] , segment);
    }
  }


/*
  LPRINT("C21\n");
  for (i=0; i<totpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", COV21[i + j * totpnts]);
   }
    LPRINT("\n");
  }
*/  


  // *** MuT

  if (PL>=PL_STRUCTURE) { LPRINT("calculating matrix for the mean...\n"); }
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
  LPRINT("MuT\n");
  for (i=0; i<totpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", MuT[i + j * totpnts]);
   }
    LPRINT("\n");
  }
*/
  //  assert(false);
 
  // *** C11
  if (PL>=PL_STRUCTURE) { LPRINT("calculating cov matrix...\n"); }
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
  LPRINT("U11\n");
  for (i=0; i<spatialpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
   }
   LPRINT("\n");
  }
*/


  // *** sqrt C11
  if (PL>=PL_STRUCTURE) { LPRINT("Cholesky decomposition of Cov11...\n"); }
  row=spatialpnts;
  // dchdc destroys the input matrix; upper half of U22 contains result!
  // dpotrf	F77_CALL(dpotrf)("Upper", &m, REAL(ans), &m, &i);
  F77_CALL(dpotrf)("Upper", &row, U11, &row, &f77err);

/*
  for (i=0; i<spatialpnts; i++) {
    for (j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
   }
    LPRINT("\n");
  }
*/

  if (f77err!=NOERROR) {
    if (PL>=PL_ERRORS) { LPRINT("U11: Error code F77_CALL(dpotrf) = %d\n", f77err); }
    err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  }

  err = FieldReturn(cov);

 ErrorHandling: // and NOERROR...
  FREE(xx);
  if (S!=NULL) {
      S->totpnts = totpnts;
      S->spatialpnts = spatialpnts;
      S->back = back;
      S->initial = initial;
      S->ntime = timecomp[XLENGTH];
  }
  if (!storing && err!=NOERROR) {
    FREE(MuT);
    FREE(U11);
    FREE(U22);
    FREE(G); 
    FREE(res0); 
  } else {
    if (S != NULL) {
      S->U22=U22;
      S->U11=U11;
      S->MuT=MuT;
      S->G=G;
      S->res0=res0;
      //assert(false);
     }
  }
  if (COV21!=NULL) {
      if (S != NULL && debugging)  S->Cov21 = COV21;
      else UNCONDFREE(COV21);
  }
  if (Inv22!=NULL) {
      if (S != NULL && debugging)  S->Inv22 = Inv22;
      else UNCONDFREE(Inv22);
  }

  cov->simu.active = err == NOERROR;

  //printf("active=%d\n", cov->simu.active); assert(false);

  return err;
}


void sequentialpart(res_type *res, long totpnts, int spatialpnts, int ntime,
		    double *U11, double *MuT, double *G) {
  res_type *rp, *oldrp;
  int n, i, k, j, mutj;
  rp = res + totpnts;
  oldrp = res;
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


void do_sequential(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) 
{  
  cov_model *next = cov->sub[0];
  location_type *loc = Loc(cov);
  sequ_storage
    *S = cov->Ssequ;
  assert(S != NULL); 
  
  int vdim = next->vdim[0];
  long  i, j, k,
    totpnts = S->totpnts;
  double *G,*U22, *U11, *MuT;
  res_type *res0,
    *res = cov->rf; 
  bool loggauss = GLOBAL.gauss.loggauss;
 GLOBAL.gauss.loggauss = false;
 
  assert(res != NULL);

  assert(S != NULL); 

 
  //  printf("totpnts %ld %d %d\n", totpnts, S->initial, S->spatialpnts);
  // assert(false);
 
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
  MEMCOPY(res, res0, sizeof(res_type) * totpnts * vdim);
  sequentialpart(res, totpnts, S->spatialpnts, S->ntime - S->back, 
		 U11, MuT, G);

  if (loggauss) {
    long vdimtot = loc->totalpoints * cov->vdim[0];
    for (i=0; i<vdimtot; i++) res[i] = exp(res[i]);
  }

  GLOBAL.gauss.loggauss = loggauss;

}
