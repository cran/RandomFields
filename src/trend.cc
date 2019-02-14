
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2011 -- 2015 Marco Oesting & Martin Schlather
               2015 -- 2017 Martin Schlather


 Handling of the different possibilities to pass the trend

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nondomain models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc
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

/*Note:
the parameter 'polycoeff' is hidden in R; it is a vector consisting of the
coefficients of the corresponding trend polynomials:
the first choose(polydeg[1]+d,d) components belong to the first polynomial,
the next choose(polydeg[2]+d,d) to the second one,...
the corresponding monomial functions are of the following order:
1, x, x^2, ..., x^k, y, x y, ... x^(k-1) y, y^2, x y^2..., y^k,
z, x z, x^2 z, ...., x^(k-1) z, y z, x y z, x^2 y z, ..., z^k
*/


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

#include "questions.h"
#include "primitive.others.h"
#include "Processes.h"
#include "shape.h"
#include "Coordinate_systems.h"
#include "rf_interfaces.h"



int binomialcoeff(int n, int k) {
  //programmed as in wikipedia
  int i, res=1;
  
  if((k < 0) || (k > n)) return 0;
  if(k > n-k) k = n-k; //symmetry
  for(i=0; i<k; i++) {
     res *= n-i;
     res /= i+1; //two steps because of integer division
  }
  return res;
}



//////////////////////////////////////////////////////////////////////
//    trend
//////////////////////////////////////////////////////////////////////

void trend(double *x, model *cov, double *v){
  int vdim = VDIM0;

  assert(isnowShape(cov) || isnowTrend(cov));
  if (hasAnyEvaluationFrame(cov)) {
    int vSq = vdim * vdim;
    //    PMI0(cov);
    // BUG;
    for (int i=0; i<vSq; i++) v[i]=0.0;
    return;
  }

  model *musub = cov->kappasub[TREND_MEAN];
  double *mu = P(TREND_MEAN);
  if (musub != NULL) {
    FCTN(x, musub, v);
  } else for (int i=0; i<vdim; i++) v[i] = ISNAN(mu[i]) ? 1.0 : mu[i];
  // 1.0 notwendig fuer likelihood berechnung;
}

void kappatrend(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  *nr = i == TREND_MEAN ? SIZE_NOT_DETERMINED : -1; 
  *nc = 1;
}


bool settrend(model *cov) {
  model *musub = cov->kappasub[TREND_MEAN]; 
  isotropy_type iso = CONDPREVISO(0); 
  if (!isFixed(iso)) return false;
  set_type(OWN, 0, TrendType);
  if (musub == NULL) { // derzeit alles coordinatesystems falls trend,
    // da FCTN dies voraussetzt
    set_iso(OWN, 0, PREVISO(0));
    set_xdim(OWN, 0, PREVXDIM(0));
  } else {
    set_iso(OWN, 0, isCartesian(iso) ? CARTESIAN_COORD
	    : isEarth(iso) ? EARTH_COORD
	    : isSpherical(iso) ? SPHERICAL_COORD
	    : ISO_MISMATCH);
    set_xdim(OWN, 0, PREVXDIM(0));
  }
  return true;
}


bool allowedItrend(model *cov) {
  model *musub = cov->kappasub[TREND_MEAN]; 
 
  bool *I = cov->allowedI;
  if (musub == NULL) return allowedItrue(cov);
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  I[CARTESIAN_COORD] = I[EARTH_COORD] = I[SPHERICAL_COORD] = true;
  return false;
}


int checktrend(model *cov){
  // if (cov->ncol[TREND_LINEAR] > 0 || cov->ncol[TREND_FCT]>0
  //    || cov->ncol[TREND_PARAM_FCT] > 0)
  //  return(ERRORNOTPROGRAMMED);

  model 
    *calling = cov->calling,
    *musub = cov->kappasub[TREND_MEAN];
  int i,  err,
    vdim = 0, 
    logdim = OWNLOGDIM(0);
  //SEXP Rx, fctbody, envir;

  //PMI(cov->calling->calling);
  
  if (musub!=NULL && !equalsCoordinateSystem(PREVISO(0)))
    RETURN_ERR(ERRORWRONGISO);
  assert(calling != NULL);
  if ((musub != NULL) xor PisNULL(TREND_MEAN))
    SERR("no trend argument given");
  bool ok = calling == NULL ||   // passt das so alles??
    ( CALLINGNR == MULT && musub!=NULL && equalsnowTrend(calling) ) ||
    //   (musub == NULL && // warum dies verlangt wird ist nicht klar
    CALLINGNR == TREND_PROC || CALLINGNR == SIMULATE ||
    CALLINGNR == GAUSSPROC || CALLINGNR == SPECIFIC ||
    (CALLINGNR == PLUS &&
     (calling->calling == NULL || isnowProcess(calling->calling) ||
      equalsnowInterface(calling->calling) || equalsnowTrend(calling) ));

  //PMI(cov);
  if (!ok) {
    SERR("trend model may not be used within the given frame.");
  }
  if (!hasTrendFrame(cov) && !hasAnyEvaluationFrame(cov)) {
    //TREE0(cov); PMI(cov);
    ILLEGAL_FRAME;
  }

  if ((cov->matrix_indep_of_x = musub == NULL)) {
    if (calling != NULL && CALLINGNR != PLUS && isnowProcess(calling) && 
	!equalsnowInterface(calling) && CALLINGNR != TREND_PROC) {
      if (OWNLASTSYSTEM > 0) BUG;

      //BUG;
      //set_type(OWN, 0, ShapeType); // 26.9.18
    }
    vdim = cov->nrow[TREND_MEAN];
  } else {
    ASSERT_ONESYSTEM;
    if ((err = CHECK(musub, logdim, OWNXDIM(0), ShapeType, XONLY, OWNISO(0),
		     SUBMODEL_DEP, TrendType)) != NOERROR) RETURN_ERR(err);
    if (isnowRandom(musub)) NotProgrammedYet("mixed effects");
    vdim = musub->vdim[0];
  } 

  if (vdim <= 0) {
    if (calling == NULL || (vdim = calling->vdim[0]) <= 0)
      SERR("multivariate dimension for trend cannot be determined.");
    PALLOC(TREND_MEAN, vdim, 1);
    for(i=0; i<vdim; i++) P(TREND_MEAN)[i] = 0.0;
  }
  VDIM0 = VDIM1 = vdim;
  set_iso(OWN, 0, cov->matrix_indep_of_x ? IsotropicOf(OWNISO(0)) 
	  : CoordinateSystemOf(OWNISO(0)));

  RETURN_NOERROR;

}



void rangetrend(model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  //P(TREND_MEAN]: mu / mean
  
  range->min[TREND_MEAN] = RF_NEGINF;
  range->max[TREND_MEAN] = RF_INF;
  range->pmin[TREND_MEAN] = - 10^10;
  range->pmax[TREND_MEAN] = 10^10;
  range->openmin[TREND_MEAN] = true;
  range->openmax[TREND_MEAN] = true;
}

				
void GetInternalMeanI(model *cov, int vdim, double *mean){
  // assuming that mean[i]=0.0 originally for all i
  int i;
  if (COVNR == TREND) {
    if (cov->ncol[TREND_MEAN]==1) {
      if (cov->nrow[TREND_MEAN] != vdim || cov->kappasub[TREND_MEAN] != NULL) {
	for (i=0; i<vdim; i++) mean[i] = RF_NA;
	return; // only scalar allowed !
      }
      for (i=0; i<vdim; i++) mean[i] += P(TREND_MEAN)[i];
    } 
  } else if (COVNR == CONST && equalsnowTrend(cov)) {
    for (i=0; i<vdim; i++) mean[i] += P(CONST_C)[i];
  } else if (equalsnowTrend(cov)) {
    FCTN(ZERO(cov), cov, mean);//to be improved! to do
  }
  if (COVNR == PLUS || COVNR == TREND) {
    for (i=0; i<cov->nsub; i++) GetInternalMeanI(cov->sub[i], vdim, mean);
  }
}

void GetInternalMean(model *cov, int vdim, double *mean){
  for (int i=0; i<vdim; mean[i++]=0.0);
  GetInternalMeanI(cov, vdim, mean);
}



int checkTrendproc(model *cov) { // auch fuer Trendproc
  model *next = cov->sub[0],
    *musub = cov->kappasub[TREND_MEAN];
  int err;
  ASSERT_ONESYSTEM;
  if ((next != NULL) + (musub != NULL) + !PisNULL(TREND_MEAN) != 1)
    SERR("either 'mu' or a 'sub model' must be given");
  if (musub != NULL) {
    next = cov->sub[0] = musub;
    cov->kappasub[TREND_MEAN] = NULL;
  }
  if (next != NULL) {
     if ((err = CHECK_PASSTF(next, ShapeType, SUBMODEL_DEP, TrendType)) !=
	NOERROR) RETURN_ERR(err);

     //PMI(cov);
     
     assert(isnowShape(next)); // inlcuding trend
    setbackward(cov, next);
    VDIM0 = next->vdim[0]; 
    VDIM1 = next->vdim[1];
  } else {
    VDIM0 = cov->nrow[TREND_MEAN];
    VDIM1 = cov->ncol[TREND_MEAN];
  }

  RETURN_NOERROR;
}
 

int init_Trendproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){// auch fuer Trendproc
  // FRAME_ASSERT_GAUSS;
  if (VDIM0 != 1) NotProgrammedYet("");
  int err;
  if (cov->sub[0] != NULL && (err = check_fctn(cov)) != NOERROR)
    goto ErrorHandling;
  if ((err = ReturnOwnField(cov)) != NOERROR) goto ErrorHandling;
  if (PL>= PL_STRUCTURE) { PRINTF("\n'%.50s' is now initialized.\n", NAME(cov));}

 ErrorHandling:  
  cov->simu.active = err == NOERROR;
  RETURN_ERR(err);
}

void do_Trendproc(model *cov, gen_storage  VARIABLE_IS_NOT_USED *s){
  double
    *res = cov->rf;
  assert(res != NULL);
  errorloc_type errorloc_save;
  KEY_type *KT = cov->base;
  STRCPY(errorloc_save, KT->error_loc);
  SPRINTF(KT->error_loc, "%.50s%.50s", errorloc_save, "add trend model");
  if (cov->sub[0] != NULL) Fctn(NULL, cov, res);
  else {
    location_type *loc = Loc(cov);
    int
      vdim = VDIM0,
      vdimtot = loc->totalpoints * vdim;
    double mu[MAXVDIM];
    assert( cov->nrow[TREND_MEAN] <= MAXVDIM );
    MEMCOPY(mu, P(TREND_MEAN), sizeof(double) * cov->nrow[TREND_MEAN]);
    for (int i=0; i<vdimtot; i++) res[i] = mu[i % vdim];
  }
  STRCPY(KT->error_loc, errorloc_save);
  return; 
}


Types Typetrend(Types required, model *cov, isotropy_type required_iso){
  if (cov->kappasub[TREND_MEAN] == NULL) return required;
  return equalsCoordinateSystem(required_iso) ? required : BadType;
}
