#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <sys/timeb.h>
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include "RFCovFcts.h"
#include "MPPFcts.h"
#include <unistd.h>

/*

// see RFCovFct.cc for examples of possible definitions.

typedef int (*checkfct)(double* p, int timespace dim, 
			SimulationType method);
// The function checks whether the intrinsic parameters p[KAPPAI], p[KAPPAII],
// etc, the dimension dim, and the method match or work together. It
// returns NOERROR or an error code according to error.h

typedef void (*infofct)(double *p, int *maxdim, int *CEbadlybehaved); 
// the function has as input parameter p, the intrinsic parameter
// vector from which p[KAPPAI], p[KAPPAII], etc. is used.
// maxdim is the maximum dimension for which the model is valid; use INFDIM
//        if the model is allowed for all dimensions
// CEbadlybehaved is only 0 (false), 1(true), 2 depending on whether
//        the circulant embedding method frequently fails for the model
// both maxdim and CEbadlybehaved may depend on p[KAPPAI], p[KAPPAII], etc

typedef void (*rangefct)(int dim , int * index, double* range);
// input parameters are
//   dim: the dimension of the random field
//   index: currently called region, starting with 0
// output parameters are
//   index: returns -2 if dimension is not allowed
//                  -1 if end of list has been reached or last region is 
//                     returned
//                  ++(*index) if further regions are available
//   range: a vector of n.kappa blocs of 4 elements, where n.kappa
//          is the number of instrinsic parameters. In each bloc, the first
//          two elements give the -- theoretically funded -- range of the 
//          parameter; the following two elements give the range which is
//          usually not exceeded in practice (worthless or difficult to
//          simulate)

typedef double (*covfct)(double *x, double*p, int dim);
// all parameters are input parameters:
// x : vector of length 1 for FULLISOTROPIC, of length 2 for SPACEISOTROPIC
//     and of length dim for ANISOTROPIC -- currently no ANISOTROPIC model
//     has been programmed yet.
// p : p[KAPPAI], p[KAPPAII], etc 
// dim : currently unused, except for checking
// IMPORTANT! covfct expect the standard model definition with variance 1 and 
//            scale=1. That is, p[VARIANCE], p[SCALE], p[ANISO] may not be 
//            used within covfct. (These parameters are set elsewhere.)
// The function returns the function value for a covariance model and
//  - gamma(h) for a variogram model

typedef double (*isofct)(double*, double*);  
// all parameters are input parameters:
// x : vector of length 1 for FULLISOTROPIC, of length 2 for SPACEISOTROPIC
//     and of length dim for ANISOTROPIC -- currently no ANISOTROPIC model
//     has been programmed yet.
// p : p[KAPPAI], p[KAPPAII], etc 

typedef double (*natscalefct)(double* p, int scaling);
// all parameters are input parameters:
// p : p[KAPPAI], p[KAPPAII], etc 
// scaling : NATSCALE_EXACT, NATSCALE_APPROX, NATSCALE_MLE
// natscalefct returns the scale parameter such that, 
// for x=1, the covariance function value is 0.05. if case
// the scale parameter is not known then the function natscalefct
// should return 0.0 if scaling=NATSCALE_EXACT; 
// if scaling=NATSCALE_APPROX or scaling=NATSCALE_MLE values are return
// as approximation or of interest in MLE of parameters to put
// the parameters p[KAPPAI], etc and p[SCALE] into a somehow orthogonal
// direction.

typedef double (*randommeasure)(double *p);     
// p : p[KAPPAI], p[KAPPAII], etc 
// the function returns a random draw from the spectral measure in the 
// two dimensional spectral turning bands method

nr = IncludeModel(
  char *name,    // name of the model appearing in R
  int kappas,    // number of specific parameters
  checkfct,      // see above
  int isotropic, // values are: FULLISOTROPIC, SPACEISOTROPIC, ANISOTROPIC
  bool variogram,// is the model a variogramm, e.g. gamma(h)=|h|
  //                if so, then the covaiance function definition must 
  //                be C(h) = -\gamma(h); the derivatives accordingly
  infofct info,  // see above
  rangefct range // see above
  );
addCov(int nr,                  // the number returned by IncludeModel
       covfct cov,              // see above
       isofct derivative,       // the derivative of a FULLISOTROPIC model
       //                          or the derivative w.r.t. the spatial 
       //                          component in case of a SPACEISOTROPIC model; 
       //                          used in TBM3 method for product models; see 
       //                          also typedef of isofct
       natscalefct naturalscale // see above
  );
addTBM(int nr,                // the number returned by IncludeModel
       isofct cov_tbm2,       // the solved Abel intregral for TBM2
       isofct cov_tbm3,       // d(hC(h))/dh -- will become obsolete, 
       //                        since it can be composed from cov and
       //                        derivative
       randommeasure spectral // see above
  );
// other addons exist, but are rarely used 

*/


double gCauchy(double *x, double *p, int effectivedim){
  return pow(1.0 + pow(fabs(*x), p[KAPPAI]), -p[KAPPAII]/p[KAPPAI]);
}

double ScalegCauchy(double *p,int scaling) {
  switch(scaling) {
  case NATSCALE_EXACT: case NATSCALE_APPROX:
    return pow(pow(0.05,-p[KAPPAI]/p[KAPPAII])-1.0,-1.0/p[KAPPAI]); 
    break;
  case NATSCALE_MLE: 
    // should be changed! (long tails!)
    return pow(pow(0.05,-p[KAPPAI]/p[KAPPAII])-1.0,-1.0/p[KAPPAI]);
    break;
  default: assert(false);
  }
}

double DgCauchy(double *x, double *p){
  register double ha,y;
  if ((y = fabs(*x))==0.0) 
    return ((p[KAPPAI]>1.0) ? 0.0 : (p[KAPPAI]<1.0) ? -INFTY : -p[KAPPAII]); 
  ha=pow(y, p[KAPPAI] - 1.0);
  return  -  p[KAPPAII] * ha * pow(1.0 + ha * y,-p[KAPPAII] / p[KAPPAI] - 1.0);
}

int checkgCauchy(double *param, int timespacedim, SimulationType method){
  if ((param[KAPPAI]<=0) || (param[KAPPAI]>2.0)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAI]);
    return ERRORCOVFAILED;
    if (param[KAPPAII]>0) return 0;
  }
  if (param[KAPPAII]<=0) {
    strcpy(ERRORSTRING_OK,"0<kappa2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAII]);
    return ERRORCOVFAILED;
    }
  if (method==CircEmbedIntrinsic || method==CircEmbedCutoff)
  {
    if (timespacedim>2) 
    {
      strcpy(ERRORSTRING_OK,"genuine total dim<=2");
      sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
      return ERRORCOVFAILED;
    }
  }
  return 0;
}

static double range_gCauchy[8] = {0, 2, 0.05, 2, 
				  0, RF_INF, 0.05, 10.0};
void rangegCauchy(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_gCauchy, sizeof(double) * 8);
}

void infogCauchy(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}


void addusersfunctions() {
// replace this function by something similar to the code
// found below in the comment
}
/*
void addusersfunctions() {
  int nr;
  nr=IncludeModel("gencauchy2", 2, checkgCauchy, FULLISOTROPIC, false,
		  infogCauchy, rangegCauchy);
  addCov(nr,gCauchy, DgCauchy, ScalegCauchy);
  addTBM(nr, NULL, NULL, NULL);
}
*/

