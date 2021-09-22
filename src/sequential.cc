/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by sequential method

 Copyright (C) 2001 -- 2017 Martin Schlather, 

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

#include <stdio.h>  
#include "def.h"
#include <Basic_utils.h>
#include <R_ext/Lapack.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "Processes.h"
#include "Coordinate_systems.h"


#define SEQU_BACK (COMMON_GAUSS + 1)
#define SEQU_INIT (COMMON_GAUSS + 2)


bool debugging = true;


int check_sequential(model *cov) {
#define nsel 4
  model *next=cov->sub[0];
  int err,
    dim = ANYDIM; // taken[MAX DIM],
  sequ_param *gp  = &(GLOBAL.sequ);
  location_type *loc = Loc(cov);

  if (!loc->grid && !loc->Time) 
    SERR1("'%.50s' only possible if at least one direction is a grid", NICK(cov));

  kdefault(cov, SEQU_BACK, gp->back);
  kdefault(cov, SEQU_INIT, gp->initial);
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);

  
  if ((err = CHECK(next, dim, dim, PosDefType, XONLY, 
		   SymmetricOf(OWNISO(0)), // todo : eigentlich sollte coordinate of reichen; aber unten alg durchschauen.
		   SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  if (next->pref[Sequential] == PREF_NONE) RETURN_ERR(ERRORPREFNONE);
  setbackward(cov, next);

  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);


  RETURN_NOERROR;
}


void range_sequential(model  VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;

  range->min[SEQU_BACK] = 0;
  range->max[SEQU_BACK] = RF_INF;
  range->pmin[SEQU_BACK] = 0.1;
  range->pmax[SEQU_BACK] = 10;
  range->openmin[SEQU_BACK] = true;
  range->openmax[SEQU_BACK] = true; 

  range->min[SEQU_INIT] = RF_NEGINF;
  range->max[SEQU_INIT] = RF_INF;
  range->pmin[SEQU_INIT] = - 10;
  range->pmax[SEQU_INIT] = 10;
  range->openmin[SEQU_INIT] = false;
  range->openmax[SEQU_INIT] = true; 
}

// (U1, U2) ~ N(0, S) =>
// U1 | U2 ~ N(S_{12} S_{22}^{-1} u2, S_11 - S12 S_22^{-1} S_21)

// start mit S_22; dann nur zeilenweise + Einschwingen

int init_sequential(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  model *next = cov->sub[0];
  location_type *loc = Loc(cov);
  if (loc->distances) RETURN_ERR(ERRORFAILED);

  int withoutlast, d, endfor, l, 
    err = NOERROR,
    dim = ANYDIM,
    spatialdim = dim - 1,
    vdim = next->vdim[0],
    max = GLOBAL.direct.maxvariables,
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
  double
    *res0 = NULL;
  sequ_storage* S = NULL;
  long  
    timelength = loc->grid ? loc->xgr[spatialdim][XLENGTH] :loc->T[XLENGTH],
    spatialpnts = loc->totalpoints / timelength,
    totpnts = back * spatialpnts, 
    totpntsSQ =  totpnts * totpnts,
    spatialpntsSQ = spatialpnts * spatialpnts,
    spatialpntsSQback = totpnts * spatialpnts;

  bool 
    storing = GLOBAL.internal.stored_init,
    Time = loc->Time;
 
  if (hasEvaluationFrame(cov)) {
    RETURN_NOERROR;
  }
  
  cov->method = Sequential;
     
  if (!loc->grid && !Time) 
    GERR("last component must be truely given by a non-trivial grid");

  if (DefList[NEXTNR].implemented[Sequential] != IMPLEMENTED) {
    err=ERRORNOTDEFINED; 
    goto ErrorHandling;
  }
  
  if (VDIM0 > 1) {
      err=ERRORNOMULTIVARIATE; 
      goto ErrorHandling;   
  }
  
  if (totpnts > max) // (totpnts * vdim > max)
    GERR6("'%.50s' valid only if the number of locations is less than '%.50s' (=%d) . Got %d * %ld = %ld.", NICK(cov), direct[DIRECT_MAXVAR_PARAM],
	  max, back, spatialpnts, totpnts);
   
  if (timelength <= back) {
    GERR2("the grid in the last direction is too small; use method '%.50s' instead of '%.50s'",
	  DefList[DIRECT].nick, DefList[SEQUENTIAL].nick);
  } 
  if (back < 1) back = max / spatialpnts;

  if ((U22 = (double *) MALLOC(sizeof(double) * totpntsSQ))==NULL ||
      (Inv22 = (double *) MALLOC(sizeof(double) * totpntsSQ))==NULL ||
      (COV21 = (double *) MALLOC(sizeof(double) * spatialpntsSQback))==NULL ||
      (U11 = (double *) MALLOC(sizeof(double) * spatialpntsSQ))==NULL ||
      (MuT = (double *) MALLOC(sizeof(double) * spatialpntsSQback))==NULL ||
      (G = (double *) MALLOC(sizeof(double) * totpnts))==NULL ||
      (res0 = (double *) MALLOC(sizeof(double) * vdim *
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
  TransformLoc(cov, &xx, false);
  loc = Loc(cov);
  assert(loc->caniso == NULL);
  if (loc->grid) loc->xgr[spatialdim][XLENGTH]=timelength; 
  else loc->T[XLENGTH] = timelength;
 
  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  long k, k0, k1, k2, segment;
  int  row, sub_err,
    icol, irow;
  double *y;
  y = (double*) MALLOC(dim * sizeof(double));

  // *** S22
  // loc->i sollte beim sub eigentlich nie gebraucht werden!!!
  // somit nur zur erinnerung die laufvariablen icol und irow verwendet
  if (PL >= PL_SUBDETAILS) { LPRINT("covariance matrix...\n"); }
  k = 0;
  for (k0 = icol = 0; icol<totpnts; icol++, k0+=dim) {
    k += icol;
    for (segment = icol * (totpnts + 1), irow=icol, k2=k0; 
	 irow < totpnts; irow++) {
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
      printf("%10g ", U22[k++]); //
    }
//    printf("\n");
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
  int jj;
  jj = 0;
  k = 0;
  withoutlast = totpnts - spatialpnts;
  for (int i=withoutlast * totpnts; i<totpntsSQ; ) {
    jj += spatialpnts;
    endfor = jj + withoutlast;
    for (; jj<endfor; ) {
	COV21[jj++] = U22[i++];
    }
    endfor = k + spatialpnts;
    for (; k<endfor; ) U11[k++] = U22[i++];
  }
  

  
  // *** S21 rest
  k = 0;
  y[spatialdim] = timecomp[XSTEP] * back;
  for (k0 = icol = 0; icol<spatialpnts; icol++, k0+=dim) { // t_{n+1}
      for (k2=irow=0; irow<spatialpnts; irow++) { // t_1
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
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
    }
    LPRINT("\n");
  }

  LPRINT("S22 %d %d\n", totpnts,spatialpnts );
  for (i=0; i<totpnts; i++) {
    LPRINT("%2d: ", i);
    for (int j=0; j<totpnts; j++) {
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
  if ((sub_err = Ext_chol(U22, row)) != NOERROR) {
    if (PL>=PL_ERRORS) {
      LPRINT("Error code in cholesky = %d\n", sub_err);
    }
    err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  } 

   
  // *** inverse of S 22 
  if (PL>=PL_STRUCTURE) { LPRINT("inverse of Cov22...\n"); }
  MEMCOPY(Inv22, U22, sizeof(double) * totpntsSQ);
  Ext_chol2inv(Inv22, row);
  k = 0;
  for (int i=0; i<totpnts; i++) {
    k += i;
    for (segment=i* totpnts+i; segment<totpntsSQ; segment += totpnts) {
	Inv22[k++] = Inv22[segment]; 
//	LPRINT("sg %d %10g %d \n", k, Inv22[segment] , segment);
    }
  }


/*
  LPRINT("C21\n");
  for (i=0; i<totpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", COV21[i + j * totpnts]);
   }
    LPRINT("\n");
  }
*/  


  // *** MuT

  if (PL>=PL_STRUCTURE) { LPRINT("calculating matrix for the mean...\n"); }
  l = 0;
  for (int i=0; i<spatialpntsSQback; i+=totpnts) {
    for (int j=0; j<totpntsSQ; j+=totpnts) {
	double dummy;
	dummy = 0.0;
	for (k=0; k<totpnts; k++) dummy += COV21[i + k] * Inv22[j + k];
	MuT[l++] = dummy;
    }
  }

/*
  LPRINT("MuT\n");
  for (i=0; i<totpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", MuT[i + j * totpnts]);
   }
    LPRINT("\n");
  }
*/
  //  assert(false);
 
  // *** C11
  if (PL>=PL_STRUCTURE) { LPRINT("calculating cov matrix...\n"); }
  l = 0;
  for (int i=0; i<spatialpntsSQback; i+=totpnts) {
    for (int j=0; j<spatialpntsSQback; j+=totpnts) {
	double dummy;
	dummy = 0.0;
	for (int kk=0; kk<totpnts; kk++) dummy += MuT[i + kk] * COV21[j + kk];
	U11[l++] -= dummy;
    }
  }

/*
  LPRINT("U11\n");
  for (i=0; i<spatialpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
   }
   LPRINT("\n");
  }
*/


  // *** sqrt C11
  if (PL>=PL_STRUCTURE) { LPRINT("Cholesky decomposition of Cov11...\n"); }
  row=spatialpnts;
  sub_err = Ext_chol(U11, row);
  if (sub_err!=NOERROR) {
    if (PL>=PL_ERRORS) { LPRINT("U11: Error code in chol = %d\n", sub_err); }
    err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  }

  err = ReturnOwnField(cov);

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

  RETURN_ERR(err);
}


void sequentialpart(double *res, long totpnts, int spatialpnts, int ntime,
		    double *U11, double *MuT, double *G) {
  double *rp, *oldrp;
  rp = res + totpnts;
  oldrp = res;
  for (int n=0; n<ntime; n++, rp += spatialpnts, oldrp += spatialpnts) {
    for (int i=0; i<spatialpnts; i++) G[i] = GAUSS_RANDOM(1.0);
    for (int mutj=0, k=0, i=0; i<spatialpnts; i++, k+=spatialpnts) {
      double dummy;
      double *Uk;
      Uk = &U11[k]; 
      dummy =0.0;
      for (int j=0; j<=i; j++) {
	dummy += G[j] * Uk[j];
      }
      for (int j=0; j<totpnts; j++) {
	  dummy += MuT[mutj++] * (double) oldrp[j];
      }
      rp[i] = (double) dummy; 
    }
  }
}


void do_sequential(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) 
{  
  model *next = cov->sub[0];
  sequ_storage
    *S = cov->Ssequ;
  assert(S != NULL); 
  
  int vdim = next->vdim[0];
  long 
    totpnts = S->totpnts;
  double *G,*U22, *U11, *MuT;
  double *res0,
    *res = cov->rf; 
  SAVE_GAUSS_TRAFO;
 
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
  for (int i=0; i<totpnts; i++) G[i] = GAUSS_RANDOM(1.0);
  for (int k=0, i=0; i<totpnts; i++, k+=totpnts){
    double dummy, *Uk;
    Uk = &U22[k]; 
    dummy =0.0;
    for (int j=0; j<=i; j++){
      dummy += G[j] * Uk[j];
    }
    res0[i] = (double) dummy; 
  }
  
  sequentialpart(res0, totpnts, S->spatialpnts, S->initial, U11, MuT, G);
  res0 += S->initial * S->spatialpnts;
  MEMCOPY(res, res0, sizeof(double) * totpnts * vdim);
  sequentialpart(res, totpnts, S->spatialpnts, S->ntime - S->back, 
		 U11, MuT, G);
  BOXCOX_INVERSE;

}
