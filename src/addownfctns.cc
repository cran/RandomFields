#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "RF.h"
#include "Covariance.h"
#include <unistd.h>

/*

// see CovFct.cc for examples of possible definitions.

typedef int (*checkfct)(cov_model *cov);
// The function checks whether the intrinsic parameters
// etc, the dimension dim, and the method match or work together.
// if fails with ERR("msg") when inconsistencies detected.

typedef void (*rangefct)(cov_model *cov, range_arraytype* ra);
// RF.h for definition of cov_model and range_arraytype

typedef void (*covfct)(double *x, cov_model *cov, double *v);
// first two parameters are input parameters:
// x : vector of length 1 for ISOTROPIC, of length 2 for SPACEISOTROPIC
//     and of length dim for STATIONARY or NONSTATIONARY
// p : p[KAPPA1], p[KAPPA2], etc 
// IMPORTANT! covfct expect the standard model definition with variance 1 and 
//            scale=1.  (These parameters are set elsewhere.)
// The function returns the function value for a covariance model and
//  - gamma(h) for a variogram model

typedef double (*natscalefct)(cov_model *cov, int scaling);
// all parameters are input parameters:
// scaling : NATSCALE_EXACT, NATSCALE_APPROX, NATSCALE_MLE
// natscalefct returns the scale parameter such that, 
// for x=1, the covariance function value is 0.05. if case
// the scale parameter is not known then the function natscalefct
// should return 0.0 if scaling=NATSCALE_EXACT; 
// if scaling=NATSCALE_APPROX or scaling=NATSCALE_MLE values are return
// as approximation or of interest in MLE of parameters to put
// the parameters p[KAPPA1], etc and p[SCALE] into a somehow orthogonal
// direction.

typedef double (*randommeasure)(cov_model *cov);     
// the function returns a random draw from the spectral measure in the 
// two dimensional spectral turning bands method

nr = IncludePrim(
  char *name,    // name of the model appearing in R
  int kappas,    // number of specific parameters
  checkfct,      // see above
  int isotropic, // values are: ISOTROPIC, SPACEISOTROPIC, STATIONARY,
  //                NONSTATIONARY
  bool variogram,// Is the model a variogramm, e.g. gamma(h)=|h| ?
  //                If so, then the covariance function definition must 
  //                be C(h) = -gamma(h); the derivatives accordingly
  rangefct range // see above
  );

// below, the function parameters can be set also to NULL. Then
// the corresponding methods will not be available
addCov(int nr,                  // the number returned by IncludeModel
       covfct cov,              // see above
       cov_model D,       // the derivative of a ISOTROPIC model
       //                          or the derivative w.r.t. the spatial 
       //                          component in case of a SPACEISOTROPIC model; 
       //                          used in TBM3 method for product models; see 
       //                          also typedef of isofct
       natscalefct naturalscale // see above
  );
addTBM(int nr,                // the number returned by IncludeModel
       cov_model cov_tbm2,       // the solved Abel intregral for TBM2
       cov_model cov_tbm3,       // d(hC(h))/dh -- may become obsolete, 
       //                        since it can be composed from cov and
       //                        derivative. The call of cov_tbm3, however,
       //                        can be much faster
       randommeasure spectral // see above
  );
// other addons exist, but are rarely used 

*/

/* gencauchy */
void gCauchy(double *x, cov_model *cov, double *v){
  double kappa = cov->p[0][0], delta=cov->p[1][0];
  *v = pow(1.0 + pow(*x, kappa), -delta/kappa);
}
double ScalegCauchy(cov_model *cov) {
  double kappa = cov->p[0][0], delta=cov->p[1][0];
  return pow(pow(0.05, -kappa / delta) - 1.0, -1.0 / kappa); 
}
void DgCauchy(double *x, cov_model *cov, double *v){
  double kappa = cov->p[0][0], delta=cov->p[1][0], ha, y=*x;
  if (y ==0.0) {
    *v = ((kappa > 1.0) ? 0.0 : (kappa < 1.0) ? -INFTY : -delta); 
  } else {
    ha=pow(y, kappa - 1.0);
    *v = -delta * ha * pow(1.0 + ha * y, -delta / kappa - 1.0);
  }
}
void DDgCauchy(double *x, cov_model *cov, double *v){
  double kappa = cov->p[0][0], delta=cov->p[1][0], ha, y=*x;
  if (y ==0.0) {
    *v = ((kappa==2.0) ? delta * (delta + 1.0) : INFTY); 
  } else {
    ha=pow(y, kappa);
    *v = delta * ha / (y * y) * (1.0 - kappa + (1.0 + delta) * ha)
      * pow(1.0 + ha, -delta / kappa - 2.0);
  }
}
void checkgCauchy(cov_model *cov){
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 0;
}

void rangegCauchy(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = 0.05;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = 0.05;
  range->pmax[1] = 10.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
}

void addusersfunctions() {
// replace this function by something similar to the code
// found below in the comment
// Only the call of IncludeModel is obligatory.
}

/*
void addusersfunctions() {
  int nr;
  nr=IncludeModel("gencauchy2", 2, checkgCauchy, ISOTROPIC, false, 
                  rangegCauchy);
  addCov(nr,gCauchy, DgCauchy, ScalegCauchy);
  addTBM(nr, NULL, NULL, NULL);
}
*/

