/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by Cholesky or SVD decomposition

 Copyright (C) 2001 -- 2015 Martin Schlather, 

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
#include "variogramAndCo.h"
//#include <R_ext/Linpack.h>

bool debug=false;


int check_directGauss(cov_model *cov) {
#define nsel 4
  cov_model *next=cov->sub[0];
  location_type *loc = Loc(cov);
  int j, err ; // taken[MAX DIM],
  direct_param *gp  = &(GLOBAL.direct); //
  
  ROLE_ASSERT(ROLE_GAUSS); 

  kdefault(cov, DIRECT_MAXVAR, gp->maxvariables);
  if ((err = checkkappas(cov, false)) != NOERROR) return err;
  if ((cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) &&
      (!loc->distances || cov->xdimprev!=1)) {
    return ERRORDIM;
  } 

  int jj, isotropy[2],  
    endjj = 0;
  if (!isCartesian(cov->isoown)) isotropy[endjj++] = cov->isoown;
  else isotropy[endjj++] = SymmetricOf(cov->isoown);
  Types type[2] = {PosDefType, VariogramType};
  for (jj=0; jj<endjj; jj++) {
     for (j=0; j<=1; j++) {    
       //printf("direct:: %s %s\n", ISONAMES[isotropy[jj]], TYPENAMES[type[j]]);

       //assert(cov->isoown == EARTH_COORD);
     
      if ((err = CHECK(next, cov->tsdim,  cov->xdimprev, 
		       type[j], KERNEL, isotropy[jj],
		       SUBMODEL_DEP, ROLE_COV)) == NOERROR) break;
     }
     if (err == NOERROR) break;
  }
  
  if (err != NOERROR) return err;  
  if (next->pref[Direct] == PREF_NONE) return ERRORPREFNONE;

  setbackward(cov, next);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) return err;
  if ((err = checkkappas(cov)) != NOERROR) return err;

  return NOERROR;
}



void range_direct(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;

  range->min[DIRECT_MAXVAR] = 0;
  range->max[DIRECT_MAXVAR] = MAX_DIRECT_MAXVAR; 
  range->pmin[DIRECT_MAXVAR] = 500;
  range->pmax[DIRECT_MAXVAR] = 5000;
  range->openmin[DIRECT_MAXVAR] = false;
  range->openmax[DIRECT_MAXVAR] = false; 
}


int init_directGauss(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *S) {
  cov_model *next = cov->sub[0];
  double //*xx,
    *Cov=NULL;
  int 
    err = NOERROR,
    maxvariab = P0INT(DIRECT_MAXVAR);
  int dim=cov->tsdim;
  direct_storage* s=NULL;
  location_type *loc = Loc(cov);
  // bool storing = GLOBAL.internal.stored_init; //
  // nonstat_covfct cf;
  long 
    vdim = cov->vdim[0],
    locpts = loc->totalpoints,
//    loctot = locpts *dim,
    vdimtot = vdim * locpts,
    //     vdimSqtot = vdim * vdimtot,
    // vdimtotSq = vdimtot * locpts,
    vdimSqtotSq = vdimtot * vdimtot;


  ROLE_ASSERT_GAUSS;
  EXT_NEW_STORAGE(solve);
  assert(cov->vdim[0] == cov->vdim[1]);

  cov->method = Direct;

  if ((err = alloc_cov(cov, dim, vdim, vdim)) != NOERROR) return err;

  if (vdimtot > maxvariab) {
    GERR4(" '%s' valid only for less than or equal to '%s'=%d data. Got %ld data.",
	  NICK(cov), KNAME(DIRECT_MAXVAR), maxvariab, vdimtot);
  }
  
  //printf("vdim = %d %d %d %d\n", vdim, locpts, vdimtot, vdimSqtotSq); 
  //    PMI(cov);
  
  if ((Cov =(double *) MALLOC(sizeof(double) * vdimSqtotSq))==NULL) {
    err=ERRORMEMORYALLOCATION;  
    goto ErrorHandling;
  }

  NEW_STORAGE(direct);
  s = cov->Sdirect;

  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  
  CovarianceMatrix(next, Cov); 
  assert(R_FINITE(Cov[0]));

  //PMI(cov->calling->calling->calling->calling);   
  if (false) {
    long i,j;
    PRINTF("\n");
    for (i=0; i<locpts * vdim; i++) {
       for (j=0; j<locpts * vdim; j++) {
	 PRINTF("%+2.3e ", Cov[i  + locpts * vdim * j]);
	 assert(R_FINITE( Cov[i  + locpts * vdim * j]));
       }
       PRINTF("\n");
    }
    assert(false); //
  }
  
  if (!isPosDef(next)) {
    if (isVariogram(next)) {
      long i, j, v;
      double min,
	*C = Cov;
      min = RF_INF;
      for (i=0; i< vdimSqtotSq; i++) if (Cov[i] < min) min=Cov[i];
      //       print("Cov %f\n", min);
      // Die Werte der Diagonalbloecke werden erh\"oht:
      for (C=Cov, v=0; v<vdim; v++, C += locpts) { 
	for (i=0; i<locpts; i++, C+=vdimtot) {
	  for (j=0; j<locpts; C[j++] -= min);
	}
      }
    } else {
      err = ERRORNOVARIOGRAM;
      goto ErrorHandling;
    }
  }

  if (false) {
    long i,j,
      endfor = locpts * vdim
      // endfor = 40
      ;
    PRINTF("\n");
    for (i=0; i<endfor; i++) {
       for (j=0; j<endfor; j++) {
	 if (ISNAN(Cov[i  + locpts * vdim * j])) BUG;
	 PRINTF("%+2.2f ", Cov[i  + locpts * vdim * j]);
       }
       PRINTF("\n");
    }
    //   assert(false); 
  }
   
 /* ********************** */
  /*  square root of matrix */
  /* ********************** */
  err = RU_sqrtPosDef(Cov, vdimtot, cov->Ssolve);

  if (err != NOERROR) {
    RU_getErrorString(ERRORSTRING);
    goto ErrorHandling;
  }

  if ((err = FieldReturn(cov)) != NOERROR) goto ErrorHandling; 

  if ( (s->G = (double *) CALLOC(vdimtot + 1, sizeof(double))) == NULL) {
      err=ERRORMEMORYALLOCATION;  
  }

 ErrorHandling: // and NOERROR...
  FREE(Cov);
 
  //printf("init sqrtPosDef emthod %d err=%d\n", cov->Ssolve->method, err);
  return err;
}


void do_directGauss(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *S) {  
  location_type *loc = Loc(cov);
  direct_storage *s = cov->Sdirect;
  long 
    locpts = loc->totalpoints,
    vdim = cov->vdim[0],
    vdimtot = locpts * vdim;
  double 
    *G = NULL,
    //*U = NULL,
    *res = cov->rf;  
  // bool  vdim_close_together = GLOBAL.general.vdim_close_together;

  //printf("do sqrtPosDef emthod %d\n", cov->Ssolve->method);

  SAVE_GAUSS_TRAFO;
  G = s->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)
  for (int i=0; i<vdimtot; i++) G[i] = GAUSS_RANDOM(1.0);  

  //  printf("vdimtot = %d\n", vdimtot);

  RU_sqrtRHS(cov->Ssolve, G, res);

  //res[0]res[1]=res[2]=1;
  //  printf("done vdimtot = %d\n", vdimtot);
  

  //for (int i=0; i<locpts; i++) res[i] = i; print("nonsense");

  BOXCOX_INVERSE;
}

