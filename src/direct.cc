/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by Cholesky or SVD decomposition

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
#include "Coordinate_systems.h"
#include "Processes.h"
#include "variogramAndCo.h"


int check_directGauss(model *cov) {
#define nsel 4
  model *next=cov->sub[0];
  location_type *loc = Loc(cov);
  int j,
    err; // taken[MAX DIM],
  // direct_param *gp  = &(GLOBAL.direct); //

  FRAME_ASSERT_GAUSS_INTERFACE;

     
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  if ((OWNLOGDIM(0) != PREVXDIM(0) || OWNLOGDIM(0) != OWNXDIM(0)) &&
      (!loc->distances || PREVXDIM(0) != 1)) {
    BUG; 
  } 

  err = ERRORFAILED;
  Types type = PosDefType;
  domain_type dom = KERNEL;
  isotropy_type iso = CoordinateSystemOf(OWNISO(0));
  for (j=0; j<=1; j++) {
    if ((err = CHECK(next, OWNLOGDIM(0), OWNXDIM(0), type, dom, iso,
		     SUBMODEL_DEP, GaussMethodType)) == NOERROR) break;
    type = VariogramType;
    iso = SymmetricOf(OWNISO(0));
  }
  if (err != NOERROR) RETURN_ERR(err);
  
  if (next->pref[Direct] == PREF_NONE) RETURN_ERR(ERRORPREFNONE);

  setbackward(cov, next);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}


void range_direct(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;
}

int init_directGauss(model *cov, gen_storage VARIABLE_IS_NOT_USED *S) {
  model *next = cov->sub[0];
  double //*xx,
    *Cov=NULL;
  int
    err = NOERROR,
    maxvariab = GLOBAL.direct.maxvariables;
  direct_storage* s=NULL;
  location_type *loc = Loc(cov);
  long 
    vdim = VDIM0,
    locpts = loc->totalpoints,
//    loctot = locpts *dim,
    vdimtot = vdim * locpts,
    //     vdimSqtot = vdim * vdimtot,
    // vdimtotSq = vdimtot * locpts,
    vdimSqtotSq = vdimtot * vdimtot,
    bytes = sizeof(double) * vdimSqtotSq;
  solve_param sp = GLOBAL_UTILS->solve;


  EXT_NEW_STORAGE(solve);
  assert(VDIM0 == VDIM1);

  cov->method = Direct;

  if ((err = alloc_cov(cov, OWNLOGDIM(0), vdim, vdim)) != NOERROR) RETURN_ERR(err);
 
  //  printf("maxvar %d\n", maxvariab); printf("vdimdot %d\n", vdimtot);
 
  if (vdimtot > maxvariab) {
    GERR5(" '%.50s' valid only for less than or equal to '%.50s'=%d data. Got %ld data. Consider increasing RFoptions(%.50s=)", // %ld ok
	  NICK(cov), direct[DIRECT_MAXVAR_PARAM],maxvariab, vdimtot,
	  direct[DIRECT_MAXVAR_PARAM]);
  }
  
  //  printf("direct vdim = %ld %ld %ld %ld\n", vdim, locpts, vdimtot, vdimSqtotSq); 
  //  PMI(cov->calling);

  if ((Cov =(double *) MALLOC(bytes))==NULL) {
    err=ERRORMEMORYALLOCATION;  
    goto ErrorHandling;
  }
  
  NEW_STORAGE(direct);
  s = cov->Sdirect;

  //  TREE0(cov);  PMI(cov->calling->calling);

  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  CHECKED;
   
  CovarianceMatrix(next, Cov); 
  assert(R_FINITE(Cov[0]));

  // FREE(cov->Spgs->z); printf("zzz\n");
   
  //PMI(cov->calling->calling->calling->calling);   
  if (false) {
    long i,j;
    PRINTF("\n");
    for (i=0; i<locpts * vdim; i++) {
       for (j=0; j<locpts * vdim; j++) {
	 PRINTF("%2.3f ", Cov[i  + locpts * vdim * j]);
	 assert(R_FINITE( Cov[i  + locpts * vdim * j]));
       }
       PRINTF("\n");
    }
    //    assert(false); //
  }


  assert(cov->Spgs->z != NULL);
  
  if (isnowPosDef(next)) {
    if (sp.pivot_check == Nan) sp.pivot_check = False;
    assert(cov->Ssolve != NULL);
    err = Ext_sqrtPosDefFree(Cov, vdimtot, cov->Ssolve, &sp);

    

    if (false) {
      long i,j;
      PRINTF("\n");
      for (i=0; i<locpts * vdim; i++) {
	for (j=0; j<locpts * vdim; j++) {
	  PRINTF("%2.3f ", cov->Ssolve->result[i  + locpts * vdim * j]);
	  assert(R_FINITE( cov->Ssolve->result[i  + locpts * vdim * j]));
	}
	PRINTF("\n");
      }
      //    assert(false); //
    }
   
    

  } else if (isnowVariogram(next)) {
    int r;
    double min = RF_INF,
      *C;
    if (vdim > 1 && GLOBAL.general.vdim_close_together)
      SERR("Simulation of multivariate intrinsic field per RPdirect only possible for 'vdim_close_together=FALSE'");
    for (int i=0; i < vdimSqtotSq; i++) min = Cov[i] < min ? Cov[i] : min;
    if (sp.pivot_check == Nan) sp.pivot_check = False;
    
#define max_rep 5
    for (r = 1; r<=max_rep; r++) {
      // Die Werte der Diagonalbloecke werden erh\"oht:
      double min_2r = (1 << (r-1)) * min;
      //      printf("r=%d %10g %10g bytes=%d %d\n", r, min_2r, min, bytes, bytes / 8);
      // Die Werte nur aif den Diagonalen werden erh\"oht:
      if (r < max_rep) {
	//printf("r < max_re[p vdimlot =%ld bytes=%ld\n", vdimtot, bytes);
	C = (double*) MALLOC(bytes);
	MEMCOPY(C, Cov, bytes);
      } else C = Cov;

     if (false) {
	long i,j,
	  endfor = locpts * vdim
	  // endfor = 40
	  ;
	PRINTF("\n"); 
	for (i=0; i<endfor; i++) {
	  for (j=0; j<endfor; j++) {
	    if (ISNAN(C[i  + locpts * vdim * j])) BUG;
	    PRINTF("%+2.2f ", C[i  + locpts * vdim * j]);
	  }
	  PRINTF("\n");
	}
	//   assert(false); 
      }
    
      
     double *D = C;
      for (int v=0; v<vdim; v++, D += locpts) { 
	for (int i=0; i<locpts; i++, D+=vdimtot) {
	  for (int j=0; j<locpts; D[j++] -= min_2r);
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
	    if (ISNAN(C[i  + locpts * vdim * j])) BUG;
	    PRINTF("%+2.2f ", C[i  + locpts * vdim * j]);
	  }
	  PRINTF("\n");
	}
	//   assert(false); 
      }
      
      if ((err = Ext_sqrtPosDefFree(C, vdimtot, cov->Ssolve, &sp)) == NOERROR)
	break;
    }
    if (r <= max_rep) FREE(Cov);
  } else {
    err = ERRORNOVARIOGRAM;
    goto ErrorHandling;
  }

   
  /* ********************** */
  /*  square root of matrix */
  /* ********************** */

   if (err != NOERROR) {
#ifdef DO_PARALLEL
     GERR1("Calculation of the square root (cholesky decomposition) failed. Rerun with  RFoptions(printlevel=%d)  to see details.", PL_ERRORS-PLoffset);
#else     
     Ext_getErrorString(cov->err_msg);
     goto ErrorHandling;
#endif     
   }

  if ( (s->G = (double *) CALLOC(vdimtot + 1, sizeof(double))) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling; 
  }

  if ((err = ReturnOwnField(cov)) != NOERROR) goto ErrorHandling; 
  
 ErrorHandling: // and NOERROR...
  cov->simu.active = err == NOERROR;
  RETURN_ERR(err);
}


void do_directGauss(model *cov, gen_storage VARIABLE_IS_NOT_USED *S) {  
  location_type *loc = Loc(cov);
  direct_storage *s = cov->Sdirect;
  long 
    locpts = loc->totalpoints,
    vdim = VDIM0,
    vdimtot = locpts * vdim;
  double 
    *G = NULL,
    //*U = NULL,
    *res = cov->rf;  
  // bool  vdim_close_together = GLOBAL.general.vdim_close_together;

  //  PMI(cov);
  
  SAVE_GAUSS_TRAFO;
  G = s->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)

  for (int i=0; i<vdimtot; i++) G[i] = GAUSS_RANDOM(1.0);  

  //  printf("vdimtot = %d\n", vdimtot);

   Ext_sqrtRHS(cov->Ssolve, G, res);

  //res[0]res[1]=res[2]=1;
  //  printf("done vdimtot = %d\n", vdimtot);
  

  //for (int i=0; i<locpts; i++) res[i] = i; print("nonsense");

  BOXCOX_INVERSE;
}

