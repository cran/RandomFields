#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include "RFCovFcts.h"
#include "MPPFcts.h"
#include <unistd.h>

/*

// see RFCovFct.cc for examples of possible definitions.

typedef int (*checkfct)(double* p, int reduceddim, 
			SimulationType method);
// The function checks whether the intrinsic parameters p[KAPPA1], p[KAPPA2],
// etc, the dimension dim, and the method match or work together. It
// returns NOERROR or an error code according to error.h

typedef void (*infofct)(double *p, int *maxdim, int *CEbadlybehaved); 
// the function has as input parameter p, the intrinsic parameter
// vector from which p[KAPPA1], p[KAPPA2], etc. is used.
// maxdim is the maximum dimension for which the model is valid; use INFDIM
//        if the model is allowed for all dimensions
// CEbadlybehaved is only 0 (false), 1(true), 2 depending on whether
//        the circulant embedding method frequently fails for the model
// both maxdim and CEbadlybehaved may depend on p[KAPPA1], p[KAPPA2], etc

typedef void (*rangefct)(int reduceddim , int * index, double* range);
// input parameters are
//   reduceddim: the reduced dimension of the random field, see RFsimu.h
//               for details
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
//      +OPEN in first element if interval left open and left border is integer, 
//      -OPEN in 2nd elemend if interval right open and right border is integer,
//      non-integer border values of an interval are always supposed to be
//      closed.

typedef double (*covfct)(double *x, double *p);
// all parameters are input parameters:
// x : vector of length 1 for ISOTROPIC, of length 2 for SPACEISOTROPIC
//     and of length dim for ANISOTROPIC -- currently no ANISOTROPIC model
//     has been programmed yet.
// p : p[KAPPA1], p[KAPPA2], etc 
// IMPORTANT! covfct expect the standard model definition with variance 1 and 
//            scale=1. That is, p[VARIANCE], p[SCALE], p[ANISO] may not be 
//            used within covfct. (These parameters are set elsewhere.)
// The function returns the function value for a covariance model and
//  - gamma(h) for a variogram model

typedef double (*isofct)(double *x, double *p);  
// all parameters are input parameters:
// x : vector of length 1 for ISOTROPIC, of length 2 for SPACEISOTROPIC
//     and of length dim for ANISOTROPIC -- currently no ANISOTROPIC model
//     has been programmed yet.
// p : p[KAPPA1], p[KAPPA2], etc 
// *isofct is distinguished from *covfct, since it is used for the
// the derivatives of isotropic functions

typedef double (*natscalefct)(double* p, int scaling);
// all parameters are input parameters:
// p : p[KAPPA1], p[KAPPA2], etc 
// scaling : NATSCALE_EXACT, NATSCALE_APPROX, NATSCALE_MLE
// natscalefct returns the scale parameter such that, 
// for x=1, the covariance function value is 0.05. if case
// the scale parameter is not known then the function natscalefct
// should return 0.0 if scaling=NATSCALE_EXACT; 
// if scaling=NATSCALE_APPROX or scaling=NATSCALE_MLE values are return
// as approximation or of interest in MLE of parameters to put
// the parameters p[KAPPA1], etc and p[SCALE] into a somehow orthogonal
// direction.

typedef double (*randommeasure)(double *p);     
// p : p[KAPPA1], p[KAPPA2], etc 
// the function returns a random draw from the spectral measure in the 
// two dimensional spectral turning bands method

nr = IncludeModel(
  char *name,    // name of the model appearing in R
  int kappas,    // number of specific parameters
  checkfct,      // see above
  int isotropic, // values are: ISOTROPIC, SPACEISOTROPIC, ANISOTROPIC
  bool variogram,// Is the model a variogramm, e.g. gamma(h)=|h| ?
  //                If so, then the covariance function definition must 
  //                be C(h) = -gamma(h); the derivatives accordingly
  infofct info,  // see above
  rangefct range // see above
  );

// below, the function parameters can be set also to NULL. Then
// the corresponding methods will not be available
addCov(int nr,                  // the number returned by IncludeModel
       covfct cov,              // see above
       isofct derivative,       // the derivative of a ISOTROPIC model
       //                          or the derivative w.r.t. the spatial 
       //                          component in case of a SPACEISOTROPIC model; 
       //                          used in TBM3 method for product models; see 
       //                          also typedef of isofct
       natscalefct naturalscale // see above
  );
addTBM(int nr,                // the number returned by IncludeModel
       isofct cov_tbm2,       // the solved Abel intregral for TBM2
       isofct cov_tbm3,       // d(hC(h))/dh -- may become obsolete, 
       //                        since it can be composed from cov and
       //                        derivative. The call of cov_tbm3, however,
       //                        can be much faster
       randommeasure spectral // see above
  );
// other addons exist, but are rarely used 

*/


double gCauchy(double *x, double *p){
  return pow(1.0 + pow(fabs(*x), p[KAPPA1]), -p[KAPPA2]/p[KAPPA1]);
}

double ScalegCauchy(double *p,int scaling) {
  switch(scaling) {
  case NATSCALE_EXACT: case NATSCALE_APPROX:
    return pow(pow(0.05,-p[KAPPA1]/p[KAPPA2])-1.0,-1.0/p[KAPPA1]); 
    break;
  case NATSCALE_MLE: 
    // should be changed! (long tails!)
    return pow(pow(0.05,-p[KAPPA1]/p[KAPPA2])-1.0,-1.0/p[KAPPA1]);
    break;
  default: assert(false);
  }
}

double DgCauchy(double *x, double *p){
  double ha,y;
  if ((y = fabs(*x))==0.0) 
    return ((p[KAPPA1]>1.0) ? 0.0 : (p[KAPPA1]<1.0) ? -INFTY : -p[KAPPA2]); 
  ha=pow(y, p[KAPPA1] - 1.0);
  return  -  p[KAPPA2] * ha * pow(1.0 + ha * y,-p[KAPPA2] / p[KAPPA1] - 1.0);
}

int checkgCauchy(double *param, int reduceddim, SimulationType method){
  if (method==CircEmbedIntrinsic || method==CircEmbedCutoff)
  {
    if (reduceddim>2) 
    {
      strcpy(ERRORSTRING_OK,"genuine total dim<=2");
      sprintf(ERRORSTRING_WRONG,"%d",reduceddim);
      return ERRORCOVFAILED;
    }
  }
  return 0;
}

//      +OPEN in first element if interval left open and left border is integer, 
//      -OPEN in second elemend if interval right open and right border is int,
static double range_gCauchy[8] = {OPEN, 2, 0.05, 2, 
				  OPEN, RF_INF, 0.05, 10.0};
void rangegCauchy(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = (reduceddim<=12345) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM; 
  memcpy(range, range_gCauchy, sizeof(double) * 8);
}

void infogCauchy(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
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
		  infogCauchy, rangegCauchy);
  addCov(nr,gCauchy, DgCauchy, ScalegCauchy);
  addTBM(nr, NULL, NULL, NULL);
}
*/

