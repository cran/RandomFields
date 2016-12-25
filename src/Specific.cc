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
//#include <R_ext/Linpack.h>
#include "RF.h"
#include "Operator.h"


int check_specificGauss(cov_model *cov) {
#define nsel 4
  cov_model
    *key = cov->key,
    *next=cov->sub[0];
  int err ; // taken[MAX DIM],
  //  direct_param *gp  = &(GLOBAL.direct); //
  
  ROLE_ASSERT(ROLE_GAUSS);
   
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;
  if (CovList[next->nr].Specific == MISMATCH)
    SERR1("specific method for '%s' not known", NAME(next));

  //    APMI(cov); //assert(false);
  if (key == NULL) {
#define SPEC_TYPES 4
    Types type[SPEC_TYPES] = {PosDefType, PosDefType, VariogramType, TrendType};
    int i,
      max = SymmetricOf(cov->isoown),
      iso[SPEC_TYPES] = {max, max, max, CoordinateSystemOf(cov->isoown)},
      dom[SPEC_TYPES] = {XONLY, KERNEL, XONLY, XONLY};    
    
      // APMI(cov->calling);

    for (i=0; i<SPEC_TYPES; i++) {
      if ((err = CHECK(next, cov->tsdim,  cov->tsdim,
		       type[i], dom[i], iso[i],
		       SUBMODEL_DEP, ROLE_COV)) == NOERROR) break;
    }
    if (err != NOERROR) return err;
    if (next->pref[Specific] == PREF_NONE) return ERRORPREFNONE;
  } else {
    if ((err = CHECK(key, cov->tsdim,  cov->tsdim, ProcessType,
		     XONLY, cov->isoown,
		     SUBMODEL_DEP, ROLE_GAUSS)) != NOERROR) {
      return err;
    }

    //APMI(key);

  }
  cov_model *sub = cov->key == NULL ? next : key;
  setbackward(cov, sub);
  cov->vdim[0] = sub->vdim[0];
  cov->vdim[1] = sub->vdim[1];
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) return err;

  //PMI(cov);

  return NOERROR;
}


int struct_specificGauss(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  cov_model 
    *next = cov->sub[0];
  location_type *loc = PrevLoc(cov);
   int err;
  //   PMI(cov); //assert(false);
  
  if (next->pref[Specific] == PREF_NONE) {
    return ERRORPREFNONE;
  }

  ROLE_ASSERT_GAUSS;
  assert(CovList[next->nr].Specific >= 0);

  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if ((err = covCpy(&(cov->key), next)) != NOERROR) return err;


  if ((err = CHECK(cov->key, next->tsdim, next->xdimprev, next->typus, 
		   next->domprev, next->isoprev, next->vdim,
		   next->role))!= NOERROR) {
    //PMI(cov->key); XERR(err);
    //printf("specific ok\n");
    // crash();
    return err;
  }

  cov->key->nr = CovList[cov->key->nr].Specific ;
  cov->key->role = ROLE_GAUSS;
  cov->key->typus = ProcessType;


  // APMI(cov);

 

  //PMI(cov->key);
  // APMI(cov->key);
  if ((err = STRUCT(cov->key, NULL)) != NOERROR) {
    return err;
  }

  //   APMI(cov->key);

  if ((err = CHECK(cov->key, loc->timespacedim, cov->xdimown, ProcessType,
		   XONLY, CoordinateSystemOf(cov->isoown),
		   cov->vdim, ROLE_GAUSS)) != NOERROR) {
    //PMI(cov->key); XERR(err);
    //printf("specific ok\n");
    // crash();
    return err;
  }

  //APMI(cov->key);
  //  if (next->nr == PLUS) AERR(err);
 
  return err;
}

void range_specificGauss(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){  
  GAUSS_COMMON_RANGE;
}


int init_specificGauss(cov_model *cov, gen_storage *S) {
  cov_model *key = cov->key;
  int err;

  assert(key != NULL);

  if (cov->role == ROLE_COV) {
    return NOERROR;
  }

  ROLE_ASSERT_GAUSS;

  cov->method = Specific;
  if ((err = INIT(key, 0, S)) != NOERROR) return err;

  key->simu.active = true;
  cov->fieldreturn = true;
  cov->origrf = false;
  cov->rf = key->rf;

  return err;
}


void do_specificGauss(cov_model *cov, gen_storage *S) {  
  cov_model *key = cov->key;
  double *res = cov->rf;
  SAVE_GAUSS_TRAFO;
  assert(key != NULL);
  DO(key, S);
  BOXCOX_INVERSE;
}

