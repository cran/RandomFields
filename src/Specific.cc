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

#include "def.h"
#include <stdio.h>  
#include <Basic_utils.h>
#include <R_ext/Lapack.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "operator.h"
#include "Processes.h"


int check_specificGauss(model *cov) {
#define nsel 4
  model
    *key = cov->key,
    *next= cov->sub[0];
  int err ,
    dim = ANYDIM; // taken[MAX DIM],
  //  direct_param *gp  = &(GLOBAL.direct); //

  FRAME_ASSERT_GAUSS_INTERFACE;
   
  if (DefList[NEXTNR].Specific == MISMATCH || DefList[NEXTNR].Specific ==UNSET)
    SERR1("specific method for '%.50s' not known", NAME(next));

  if (key == NULL) {
#define SPEC_TYPES 4
    Types type[SPEC_TYPES] ={PosDefType, PosDefType, VariogramType, TrendType};

    assert(OWNISO(0) == CoordinateSystemOf(OWNISO(0)));
    isotropy_type
      max = OWNISO(0), // CoordinateSystemOf(OWNISO(0)),
      iso[SPEC_TYPES] = {max, max, SymmetricOf(OWNISO(0)), max};
    domain_type
      dom[SPEC_TYPES] = {XONLY, KERNEL, XONLY, XONLY};
    int n = isAnySpherical(PREVISO(0)) ? 2 : SPEC_TYPES;
      
    ASSERT_ONESYSTEM;

    // printf("\n\n\n\n\n\n\n\n\n\nC HECK %d:::\n", cov->zaehler);
    
    for (int i=0; i<n; i++) {
      if ((err = CHECK(next, dim, dim, type[i], dom[i], iso[i],
		       SUBMODEL_DEP, EvaluationType)) // schwach, damit es
	  // bei plus durchkommt
	  == NOERROR) break;
      // printf("specific i=%d", i);      TREE0(cov);
    }
    if (err != NOERROR) RETURN_ERR(err);
    //PMI0(next);
    if (next->pref[Specific] == PREF_NONE) RETURN_ERR(ERRORPREFNONE);
  } else {
    if ((err = CHECK_PASSTF(key, GaussMethodType, VDIM0, GaussMethodType))
	!= NOERROR) {
      RETURN_ERR(err);
    }
  }
  model *sub = cov->key == NULL ? next : key;
  setbackward(cov, sub);
  VDIM0 = sub->vdim[0];
  VDIM1 = sub->vdim[1];
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  RETURN_NOERROR;
}


int struct_specificGauss(model *cov, model VARIABLE_IS_NOT_USED **newmodel) {
  model 
    *next = cov->sub[0];
  // location_type *loc = PrevLoc(cov);
  int err;
  if (next->pref[Specific] == PREF_NONE) {
    RETURN_ERR(ERRORPREFNONE);
  }

  FRAME_ASSERT_GAUSS_INTERFACE;
  assert(DefList[NEXTNR].Specific >= 0);

  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
  //  printf("copy\n"); PMI(cov);
  if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
  //  printf("copied\n"); PMI(cov);

  COPYALLSYSTEMS(PREVSYSOF(cov->key), PREVSYSOF(next), false);
  cov->key->variant = UNSET; // ???

  // printf("\n\n\n\n\n\n\n\n\nS TRUCT::::\n");
  
  if ((err = CHECK_ONLY(cov->key))!= NOERROR) {
    RETURN_ERR(err);
  }
  //  PMI0(cov->key); 

 
  SET_NR(cov->key, DefList[MODELNR(cov->key)].Specific);
  //  printf("reset modelnr\n");
  // check_passtf hier funktioniert nicht wegen plus, was zuerst
  // struct haben will
  //if ((err = CHECK_PASSTF(cov->key, GaussMethodType, VDIM0, GaussMethodType))
 
  
  cov->key->checked = true;
  
  cov->key->frame = GaussMethodType;
  set_type(PREVSYSOF(cov->key), 0, GaussMethodType);
  set_type(SYSOF(cov->key), 0, GaussMethodType);

  
  if ((err = STRUCT(cov->key, NULL)) != NOERROR) RETURN_ERR(err);

  //  printf("ok\n");

  if ((err = CHECK_PASSTF(cov->key, GaussMethodType, VDIM0, GaussMethodType))
      != NOERROR) {
   //if ((err = CHECK(cov->key, loc->timespacedim, cov->xdimown, ProcessType,
   //		   XONLY, CoordinateSystemOf(OWNISO(0)),
   //		   cov->vdim, GaussMethodType)) != NOERROR) {
    //
    // PMI(cov->key); XERR(err);

    //    printf("err = %d\n", err);
    
    RETURN_ERR(err);
  }

  //  printf("n o err\n");
   RETURN_ERR(NOERROR);
}

void range_specificGauss(model VARIABLE_IS_NOT_USED *cov, range_type *range){  
  GAUSS_COMMON_RANGE;
}


int init_specificGauss(model *cov, gen_storage *S) {
  model *key = cov->key;
  int err;

  assert(key != NULL);

  if (hasEvaluationFrame(cov)) {
    RETURN_NOERROR;
  }

  FRAME_ASSERT_GAUSS_INTERFACE;

  cov->method = Specific;
  //printf("specific : %.50s\n", NAME(key));
  if ((err = INIT(key, 0, S)) != NOERROR) RETURN_ERR(err);

  cov->simu.active = true;
  ReturnOtherField(cov, key);

  RETURN_ERR(err);
}


void do_specificGauss(model *cov, gen_storage *S) {  
  model *key = cov->key;
  double *res = cov->rf;
  SAVE_GAUSS_TRAFO;
  assert(key != NULL);
  DO(key, S);
  BOXCOX_INVERSE;
}

