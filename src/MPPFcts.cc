/*
 Authors 
 Martin Schlather, Martin.Schlather@uni-bayreuth.de 

 *********** this function is still under construction *********

 Copyright (C) 2001 Martin Schlather, 

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

/******************************************************************************
 * 
 *  if deterministic then
 *  functions are always(!) normed to max(f)==1, i.e. for any parameter 
 *  combination.         
 *  ****** or shall they be normed anyhow ??
 * 
 *   the values are finally corrected in do_addmpp by
 *      sqrt(p[VARIANCE] / integral(f^2))  in case of Gaussian random fields
 *                                         (and integral(f) for the mean)     
 *      1.0 / integral(f_+)                  in case of max-stable processes
 *
 **************************************************************************** */

#include <math.h>
//#include "MathlibCC.h"
#include <R_ext/Mathlib.h>
#include <assert.h>
#include "RFsimu.h"
#include "MPPFcts.h"
//#include "RFCovFcts.h"

static int mpp_dim;
static Real mpp_x[MAXDIM];


static Real spherical_Rsq, spherical_adddist;
Real spherical_model(Real *y) 
{
  Real distance;
  int i;
  distance = spherical_adddist;
  for (i=0; i<mpp_dim; i++) {
    register Real dummy;
    dummy = y[i]-mpp_x[i];
    distance += dummy * dummy;
  }
  if (distance>=spherical_Rsq) return 0;
  return 1.0;
}


static Real cone_Rsq, /* square of radius */ 
  cone_rsq, /* square of radius */
  cone_height_socle,
  cone_a, cone_b_socle; /* y = ax + b, for the slope */
Real cone_model(Real *y) 
{
  Real distance;
  int i;
  distance = 0.0;
  for (i=0; i<mpp_dim; i++) {
    register Real dummy;
    dummy = y[i]-mpp_x[i];
    distance += dummy * dummy;
  }
  if (distance>=cone_Rsq) return 0;
  if (distance<=cone_rsq) return cone_height_socle;
  return cone_b_socle + cone_a * sqrt(distance);
}

static Real gauss_invscale, gauss_effectiveradiussq;
Real gauss_model(Real *y)
{
  Real distance;
  int i;
  distance = 0.0;
  for (i=0; i<mpp_dim; i++) {
    register Real dummy;
    dummy = y[i]-mpp_x[i];
    distance += dummy * dummy;
  } 
  if (distance > gauss_effectiveradiussq) return 0;
  return exp(-gauss_invscale * distance);
}


#define SPHERICALRsq 1
#define SPHERICALadddist 2
void generalspherical_init(key_type *key, int balldim)
{
  int d;
  mpp_storage *bs;
  Real gamma;
  bs = (mpp_storage*) key->S;

  bs -> effectiveRadius = 0.5 *  key->param[SCALE];
  if (bs->addradius<=0.0) {bs->addradius=bs -> effectiveRadius;}
  bs -> effectivearea = 1.0;
  for (d=0; d<key->dim; d++) 
    bs->effectivearea *= (bs->length[d] + 2.0 * bs -> effectiveRadius);
  for (d=key->dim; d<balldim; d++) {
    bs -> effectivearea *= (2.0 * bs -> effectiveRadius);
  }   
#ifdef RF_GSL
  gsl_sf_result result; 
  gsl_sf_gamma_impl(1.0+0.5*balldim,&result); // _lngamma_
  gamma = result.val;
#else 
  gamma = gammafn(1.0+0.5*balldim);
#endif
  bs->c[SPHERICALRsq] = bs->effectiveRadius * bs->effectiveRadius;
  bs->integral = bs->integralsq = bs->integralpos =
    pow(SQRTPI * bs->effectiveRadius, balldim) / gamma;
  bs->maxheight = 1.0;
}

void circular_init(key_type *key)
{
  int TRUESPACE=2;
  assert(key!=NULL);
  mpp_storage *bs;
  bs = (mpp_storage*) key->S;
  assert(bs!=NULL);
  generalspherical_init(key,TRUESPACE);
}

void circularMpp(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM)
{ 
  int TRUESPACE = 2;  
  int d;
  assert(bs!=NULL);

  mpp_dim = bs->dim;
  spherical_Rsq = bs->c[SPHERICALRsq];
  *model = spherical_model;

  for (d=0; d<bs->dim; d++) {
    mpp_x[d]= bs->min[d] - bs->addradius + 
      (bs->length[d] + 2.0 * bs->addradius) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - bs->effectiveRadius;
    max[d] = mpp_x[d] + bs->effectiveRadius;
  }

  spherical_adddist = 0.0;
  for (d=bs->dim; d<TRUESPACE; d++) {
    Real dummy;
    dummy = bs->effectiveRadius * UNIFORM_RANDOM;
    spherical_adddist += dummy * dummy;
  }
}

void spherical_init(key_type *key)
{
  int TRUESPACE = 3;  
  mpp_storage *bs;
  bs = (mpp_storage*) key->S;
  generalspherical_init(key,TRUESPACE);
}

void sphericalMpp(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM)
{ 
  int TRUESPACE = 3;  
  int d;
  mpp_dim = bs->dim;
  spherical_Rsq = bs->c[SPHERICALRsq];
 *model = spherical_model;

  for (d=0; d<bs->dim; d++) {
    mpp_x[d]= bs->min[d] - bs->addradius + 
      (bs->length[d] + 2.0 * bs->addradius) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - bs->effectiveRadius;
    max[d] = mpp_x[d] + bs->effectiveRadius;
  }
  spherical_adddist = 0.0;
  for (d=bs->dim; d<TRUESPACE; d++) {
    Real dummy;
    dummy = bs->effectiveRadius * UNIFORM_RANDOM;
    spherical_adddist += dummy * dummy;
  }
}

#define CONErsq 1
#define CONERsq 2
#define CONEHEIGHTSOCLE 3
#define CONEA 4
#define CONEBSOCLE 5
void cone_init(key_type *key)
{
  // ignoring: p[MEAN], p[NUGGET]
  // note : sqrt(p[VARIANCE]) -> height
  //        p[SCALE]          -> R
  // p[KAPPA] = r/R \in [0,1]
  // p[KAPPAX]= socle (height)
  // p[KAPPAXX]=height of cone (without socle)
  // 
  // standard cone has radius 1/2, height of (potential) top: 1
  //           
  int d;
  Real cr,cb,CR,socle,height;
  mpp_storage *bs;
  bs = (mpp_storage*) key->S;
  socle = key->param[KAPPAX];
  height= key->param[KAPPAXX];
  cr= 0.5 * key->param[KAPPA] * key->param[SCALE];
  bs->c[CONErsq] = cr * cr;
  CR  = 0.5 * key->param[SCALE];
  bs->c[CONERsq] = CR * CR;
  bs->c[CONEHEIGHTSOCLE] = socle + height;
  bs->c[CONEA] =  /* a = h / (r - R) */
    height / ((key->param[KAPPA] - 1.0) * 0.5* key->param[SCALE]);
  cb = height / (1.0 - key->param[KAPPA]); /* b = h R / (R - r) */  
  bs->c[CONEBSOCLE] = socle + cb;
  switch (bs->dim) {
  case 2 :
    //
    bs->integral = bs->integralpos = 
      PI * (bs->c[CONERsq] * socle +
	    0.333333333333333333 * height *
	    (bs->c[CONErsq] + bs->c[CONERsq] + cr * CR));
    //
    bs->integralsq 
      =  PI * (bs->c[CONERsq] * socle * socle
	       + bs->c[CONErsq] * height * height 
	       + 0.5 * bs->c[CONEA] * bs->c[CONEA] *
	       (bs->c[CONERsq] *  bs->c[CONERsq] -
		bs->c[CONErsq] * bs->c[CONErsq])
	       + 1.33333333333333333 *  bs->c[CONEA] * cb *
	       (bs->c[CONERsq] * CR - bs->c[CONErsq] * cr)
	       + cb * cb * (bs->c[CONERsq] - bs->c[CONErsq])
	       );
    //
    break;
    // case 3 : bs -> integral = 4 PI / 3 * (R^3 * socle + r^3 * height)
    //                         + 4 PI int_r^R r^2 a x + b) dx
    //      bs -> integralsq = 4 PI / 3 * (R^3 * socle *socle
    //                                    + r^3 * height * height)
    //                     + 4 PI int_r^R r^2 (a x + b)^2 dx
  default : assert(false);
  }

  // here a "mean" effectiveRadius is defined + true radius for simulation;
  // this is important
  bs->effectiveRadius = CR;
  if (bs->addradius<=0.0) {bs->addradius=bs -> effectiveRadius;}

  // the calculation of the effective area can be very complicated
  // if randomised scale is used with upper end point of the
  // distribution ==oo.
  // Then either the simulation is pretty approximate, and probably
  // very bad (if a fixed, small border is used), 
  // or the simulation is very ineffective (if a large fixed border is used)
  // or the formula might be very complicated (if there are adapted borders,
  // and if the distribution density of the scale is weighted by about
  //       (length[1] + 2 scale[1]) *...*(length[dim] + 2 scale[dim])
  //

  bs -> effectivearea = 1.0;
  for (d=0; d<bs->dim; d++) 
     bs->effectivearea *= (bs->length[d] + 2.0 * bs->addradius);

  //
  bs->maxheight = height + socle; //==bs->c[CONEHEIGHTSOCLE]
}



// note. The cone_* values have to be  reset all the time,
// since the function does not know whether there have
// been other simulations of random fields, as cone_init is called
// only at the very fist time.
// Furthermore, the below construction will match with the
// general construction that randomises parameters (e.g. scale mixtures)
void cone(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM)
{ 
  int d;
  mpp_dim = bs->dim;
  cone_Rsq = bs->c[CONERsq];
  cone_rsq = bs->c[CONErsq]; 

  // note that c[CONEHEIGHT]==1 for Gaussian random fields !!
  // However, for max-stable random fields bs->height will have 
  // different values !!
  cone_height_socle = bs->c[CONEHEIGHTSOCLE];
  cone_a            = bs->c[CONEA];
  cone_b_socle      = bs->c[CONEBSOCLE];
  *model = cone_model;
  // logically the current effectiveRadius should be set here
  // no need as here constant scale, so 
  // bs -> effectiveRadius = 0.5 * key->param[SCALE];
  // can be used.  

  // last: random point
  for (d=0; d<mpp_dim; d++) {
    mpp_x[d]= bs->min[d] - bs->addradius + 
      (bs->length[d] + 2.0 * bs->addradius) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - bs->effectiveRadius;
    max[d] = mpp_x[d] + bs->effectiveRadius;
  }
}


int checkcone(key_type *key) {
  switch (key->method) {
  case  AdditiveMpp: 
    if (key->dim!=2) { 
      strcpy(ERRORSTRING_OK,"2 dimensional space");
      sprintf(ERRORSTRING_WRONG,"dim=%d",key->dim);
      return ERRORCOVFAILED;
    }
    break;
  case Nothing :
    break;
  default :   
    strcpy(ERRORSTRING_OK,"method=\"add. MPP (random coins)\"");
    strcpy(ERRORSTRING_WRONG,METHODNAMES[key->method]);
    return ERRORCOVFAILED;
  }
  if ((key->param[KAPPA]  >= 1.0)  ||
      (key->param[KAPPA]  < 0.0)  ||
      (key->param[KAPPAX] < 0.0) ||
      (key->param[KAPPAXX]< 0.0) || 
      ((key->param[KAPPAXX] + key->param[KAPPAX])==0.0)) {
    strcpy(ERRORSTRING_OK,
	   "all kappa_i non-negative, kappa1<1, and kappa2+kappa3>0.0");
    sprintf(ERRORSTRING_WRONG,"c(%f,%f,%f)",key->param[KAPPA],
	    key->param[KAPPAX],key->param[KAPPAXX]);
    return ERRORCOVFAILED;
  }	
  return NOERROR;
}

#define GAUSSINVSCALE 0
#define GAUSSRADIUSSQ 1

Real gaussInt(int d, int xi, Real sigma, Real R) {
  // int_{b(0,R) e^{-a r^2} dr = d b_d int_0^R r^{d-1} e^{-a r^2} dr
  // where a = 2.0 * xi / sigma^2
  // here : 2.0 is a factor so that the mpp function leads to the
  //            gaussian covariance model exp(-x^2)
  //        xi : 1 : integral ()
  //             2 : integral ()^2
  Real a;
  a = 2.0 * xi / (sigma * sigma);
  switch (d) {
  case 1 : 
    Real sqrta;
    sqrta = sqrt(a);
    return SQRTPI / sqrta * (2.0 * pnorm(SQRTTWO * sqrta * R, 0, 1, 1, 0) - 1.0);
  case 2 :
    return PI / a * (1.0 - exp(- a * R * R));
  case 3 :
    Real piDa;
    piDa = PI / a;
    return -2.0 * piDa * R * exp(- a * R * R )
      + piDa * sqrt(piDa) * (2.0 * pnorm(SQRTTWO * sqrt(a) * R, 0, 1,1,0) - 1.0);
  default : assert(false);
  }
}

void gaussmpp_init(key_type *key) 
{
  int d;
  mpp_storage *bs;
  bs = (mpp_storage*) key->S;
  bs -> c[GAUSSINVSCALE] = 2.0 * key->param[INVSCALE] * key->param[INVSCALE]; 
  // e^{-2 x^2} !

  bs -> c[GAUSSRADIUSSQ] =  - 0.5  * log(MPP_APPROXZERO) 
    * key->param[SCALE] * key->param[SCALE];
  bs -> effectiveRadius = sqrt(bs -> c[GAUSSRADIUSSQ]);
  if (bs->addradius<=0.0) {bs->addradius=bs -> effectiveRadius;}

  bs->integral   = gaussInt(key->dim, 1, key->param[SCALE], bs->effectiveRadius);
  bs->integralsq = gaussInt(key->dim, 2, key->param[SCALE], bs->effectiveRadius);

  bs -> integralpos = bs->integral;
  bs -> effectivearea = 1.0;
  for (d=0; d<key->dim; d++) 
     bs -> effectivearea *= (bs->length[d] + 2.0 * bs->addradius);
  bs->maxheight= 1.0;
}

void gaussmpp(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM)
{
  int d;

  gauss_invscale = bs->c[GAUSSINVSCALE];
  gauss_effectiveradiussq =  bs -> c[GAUSSRADIUSSQ];

  *model = gauss_model;  
  mpp_dim = bs->dim;

  for (d=0; d<mpp_dim; d++) {
    mpp_x[d]= bs->min[d] - bs->addradius + 
      (bs->length[d] + 2.0 * bs->addradius) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - bs->effectiveRadius;
    max[d] = mpp_x[d] + bs->effectiveRadius;
  }
}










