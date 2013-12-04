/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by Cholesky or SVD decomposition

 Copyright (C) 2001 -- 2013 Martin Schlather, 

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
 
#include "RF.h"
#include "Covariance.h"
#include <R_ext/Lapack.h>
//#include <R_ext/Linpack.h>


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

  if ((err = check_common_gauss(cov)) != NOERROR) return err;

  //  PMI(cov); //assert(false);
  if (key == NULL) {
    Types type[3] = {PosDefType, NegDefType, TrendType};
    int i,
      iso[3] = {SYMMETRIC, SYMMETRIC, NO_ROTAT_INV};
    
    for (i=0; i<3; i++) {
      if ((err = CHECK(next, cov->tsdim,  cov->tsdim, type[i],
		       cov->domown, iso[i],
		       SUBMODEL_DEP, ROLE_COV)) == NOERROR) break;
    }
    if (err != NOERROR) return err;
    if (next->pref[Specific] == PREF_NONE) return ERRORPREFNONE;
  } else {
    if ((err = CHECK(key, cov->tsdim,  cov->tsdim, ProcessType,
		     XONLY, NO_ROTAT_INV,
		     SUBMODEL_DEP, ROLE_GAUSS)) != NOERROR) {
      return err;
    }
  }
  cov_model *sub = cov->key == NULL ? next : key;
  setbackward(cov, sub);
  cov->vdim = sub->vdim;

  return NOERROR;
}


int struct_specificGauss(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  cov_model 
    *next = cov->sub[0];
  location_type *loc = cov->prevloc;
   int err;
  // PMI(cov); assert(false);
  
  if (next->pref[Specific] == PREF_NONE) {
    return ERRORPREFNONE;
  }

  ROLE_ASSERT_GAUSS;


  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if ((err = covcpy(&(cov->key), next)) != NOERROR) return err;
  cov->key->role = ROLE_GAUSS;
  cov->key->typus = ProcessType;

  //   if (next->nr == PLUS) {PMI(cov); assert(false);}
  //printf("just befor struct\n");

  if ((err = STRUCT(cov->key, NULL) != NOERROR)) return err; ;

  //   APMI(cov->key);

  if ((err = CHECK(cov->key, loc->timespacedim, cov->xdimown, ProcessType,
		   XONLY, NO_ROTAT_INV, cov->vdim, ROLE_GAUSS)) != NOERROR) {
    //PMI(cov->key);
    //printf("specific ok\n");
    // crash();
    return err;
  }

  //APMI(cov->key);
  //  if (next->nr == PLUS) AERR(err);
 
  return err;
}



int init_specificGauss(cov_model *cov, storage *S) {
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


void do_specificGauss(cov_model *cov, storage *S) {  
  cov_model *key = cov->key;
  location_type *loc = Loc(cov);
  bool loggauss = (bool) ((int*) cov->p[LOG_GAUSS])[0];
  double *res = cov->rf;

  assert(key != NULL);
  DO(key, S);
  if (loggauss) {
    int i, vdimtot = loc->totalpoints * cov->vdim;
    for (i=0; i<vdimtot; i++) res[i] = exp(res[i]);
  }
}

