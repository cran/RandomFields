/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de 

 *********** this function is still under construction *********

 Copyright (C) 2001 -- 2006 Martin Schlather, 

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

/*
   IF NEW FUNCTION IS DEFINED MAKE SURE THAT NO MORE THAN 5 OR 6
   ADDITIONALLY #defineD CONSTANTS ARE USED -- otherwise increase
   c[MAXCOV][6] in RFsimu.h
 */

// basic shapes are: sphere (2 or 3 dimensions), cone, and Gaussian curve

#include <math.h>
#include <assert.h>
#include "RFsimu.h"
#include "MPPFcts.h"

static int MPP_DIM;
static double MPP_X[MAXDIM];


static double spherical_Rsq, spherical_adddist;
double spherical_model(double *y) 
{
  double distance;
  int i;
  distance = spherical_adddist;
  for (i=0; i<MPP_DIM; i++) {
    register double dummy;
    dummy = y[i]-MPP_X[i];
    distance += dummy * dummy;
  }
  if (distance>=spherical_Rsq) return 0;
  return 1.0;
}


static double cone_Rsq, /* square of radius */ 
  cone_rsq, /* square of radius */
  cone_height_socle,
  cone_a, cone_b_socle; /* y = ax + b, for the slope */
double cone_model(double *y) 
{
  double distance;
  int i;
  distance = 0.0;
  for (i=0; i<MPP_DIM; i++) {
    register double dummy;
    dummy = y[i]-MPP_X[i];
    distance += dummy * dummy;
  }
  if (distance>=cone_Rsq) return 0;
  if (distance<=cone_rsq) return cone_height_socle;
  return cone_b_socle + cone_a * sqrt(distance);
}

static double gauss_effectiveradiussq;
double gauss_model(double *y)
{
  double distance;
  int i;
  distance = 0.0;
  for (i=0; i<MPP_DIM; i++) {
    register double dummy;
    dummy = y[i]-MPP_X[i];
    distance += dummy * dummy;
  } 
  if (distance > gauss_effectiveradiussq) return 0;
  return exp(- 2.0 * distance);
}


#define SPHERICALRsq 1
#define SPHERICALadddist 2
void generalspherical_init(mpp_storage *s, int truedim, param_type param,
			   int balldim) 
{
  int d;
  double gamma;

  assert(s!=NULL);
  s->dim = truedim;
  s->effectiveRadius = 0.5;
  if (s->addradius<=0.0) {s->addradius = s->effectiveRadius;}
  s->effectivearea = 1.0;
  for (d=0; d<truedim; d++) 
    s->effectivearea *= (s->length[d] + 2.0 * s->effectiveRadius);
  for (d=truedim; d<balldim; d++) {
    s->effectivearea *= (2.0 * s->effectiveRadius);
  }   
  gamma = gammafn(1.0 + 0.5 * balldim);
  s->c[SPHERICALRsq] = s->effectiveRadius * s->effectiveRadius;
  s->integral = s->integralsq = s->integralpos =
    pow(SQRTPI * s->effectiveRadius, balldim) / gamma;
  s->maxheight = 1.0;

}

void circular_init(mpp_storage *s, int truedim, param_type param) {
  int TRUESPACE=2;
  generalspherical_init(s, truedim, param, TRUESPACE);
}

void circularMpp(mpp_storage *s, double *min, double *max, mppmodel *model) { 
  int TRUESPACE = 2;  
  int d;
  assert(s!=NULL);

  MPP_DIM = s->dim;
  spherical_Rsq = s->c[SPHERICALRsq];
  *model = spherical_model;

  for (d=0; d<MPP_DIM; d++) {
    MPP_X[d]= s->min[d] - s->addradius + 
      (s->length[d] + 2.0 * s->addradius) * UNIFORM_RANDOM;
    min[d] = MPP_X[d] - s->effectiveRadius;
    max[d] = MPP_X[d] + s->effectiveRadius;
  }

  spherical_adddist = 0.0;
  for (d=MPP_DIM; d<TRUESPACE; d++) {
    double dummy;
    dummy = s->effectiveRadius * UNIFORM_RANDOM;
    spherical_adddist += dummy * dummy;
  }
}

void spherical_init(mpp_storage *s, int truedim, param_type param) {
  int TRUESPACE = 3;  
  generalspherical_init(s, truedim, param, TRUESPACE);
}

void sphericalMpp(mpp_storage *s, double *min, double *max, mppmodel *model) { 
  int TRUESPACE = 3;  
  int d;
  MPP_DIM = s->dim;
spherical_Rsq = s->c[SPHERICALRsq];
 *model = spherical_model;

  for (d=0; d<MPP_DIM; d++) {
    MPP_X[d]= s->min[d] - s->addradius + 
      (s->length[d] + 2.0 * s->addradius) * UNIFORM_RANDOM;
    min[d] = MPP_X[d] - s->effectiveRadius;
    max[d] = MPP_X[d] + s->effectiveRadius;
  }
  spherical_adddist = 0.0;
  for (d=MPP_DIM; d<TRUESPACE; d++) {
    double dummy;
    dummy = s->effectiveRadius * UNIFORM_RANDOM;
    spherical_adddist += dummy * dummy;
  }
}

#define CONErsq 1
#define CONERsq 2
#define CONEHEIGHTSOCLE 3
#define CONEA 4
#define CONEBSOCLE 5
void cone_init(mpp_storage *s, int truedim, param_type param)
{
  // ignoring: p[MEAN], p[NUGGET]
  // note : sqrt(p[VARIANCE])->height
  //        p[scale]         ->R
  // p[KAPPA1] = r/R \in [0,1]
  // p[KAPPA2]= socle (height)
  // p[KAPPA3]=height of cone (without socle)
  // 
  // standard cone has radius 1/2, height of (potential) top: 1
  //           
  int d;
  double cr,cb,CR,socle,height;
  assert(s!=NULL);
  s->dim = truedim;
  socle = param[KAPPA2];
  height= param[KAPPA3];
  cr= 0.5 * param[KAPPA1];
  s->c[CONErsq] = cr * cr;
  CR  = 0.5;
  s->c[CONERsq] = CR * CR;
  s->c[CONEHEIGHTSOCLE] = socle + height;
  s->c[CONEA] =  /* a = h / (r - R) */
    height / ((param[KAPPA1] - 1.0) * 0.5);
  cb = height / (1.0 - param[KAPPA1]); /* b = h R / (R - r) */  
  s->c[CONEBSOCLE] = socle + cb;
  switch (truedim) {
  case 2 :
    //
    s->integral = s->integralpos = 
      PI * (s->c[CONERsq] * socle +
	    0.333333333333333333 * height *
	    (s->c[CONErsq] + s->c[CONERsq] + cr * CR));
    //
    s->integralsq 
      =  PI * (s->c[CONERsq] * socle * socle
	       + s->c[CONErsq] * height * height 
	       + 0.5 * s->c[CONEA] * s->c[CONEA] *
	       (s->c[CONERsq] *  s->c[CONERsq] -
		s->c[CONErsq] * s->c[CONErsq])
	       + 1.33333333333333333 *  s->c[CONEA] * cb *
	       (s->c[CONERsq] * CR - s->c[CONErsq] * cr)
	       + cb * cb * (s->c[CONERsq] - s->c[CONErsq])
	       );
    //
    break;
    // case 3 : s->integral = 4 PI / 3 * (R^3 * socle + r^3 * height)
    //                         + 4 PI int_r^R r^2 a x + b) dx
    //      s->integralsq = 4 PI / 3 * (R^3 * socle *socle
    //                                    + r^3 * height * height)
    //                     + 4 PI int_r^R r^2 (a x + b)^2 dx
  default : assert(false);
  }

  // here a "mean" effectiveRadius is defined + true radius for simulation;
  // this is important
  s->effectiveRadius = CR;
  if (s->addradius<=0.0) {s->addradius = s->effectiveRadius;}

  // the calculation of the effective area can be very complicated
  // if randomised scale is used with upper end point of the
  // distribution ==oo.
  // Then either the simulation is pretty approximate, and probably
  // very bad (if a fixed, small border is used), 
  // or the simulation is very ineffective (if a large fixed border is used)
  // or the formula might be very complicated (if there are adapted borders,
  // and if the distribution density of the scale is weighted by about
  //       (length[1] + 2 scale[1]) *...*(length[truedim] + 2 scale[truedim])
  //

  s->effectivearea = 1.0;
  for (d=0; d<truedim; d++) 
    s->effectivearea*= (s->length[d] + 2.0 * s->addradius);

  //
  s->maxheight = height + socle; //==s->c[CONEHEIGHTSOCLE]
}



// note. The cone_* values have to be  reset all the time,
// since the function does not know whether there have
// been other simulations of random fields, as cone_init is called
// only at the very fist time.
// Furthermore, the below construction will match with the
// general construction that randomises parameters (e.g. scale mixtures)
void cone(mpp_storage *s, double *min, double *max, mppmodel *model)
{ 
  int d;
  MPP_DIM = s->dim;
  cone_Rsq = s->c[CONERsq];
  cone_rsq = s->c[CONErsq]; 

  // note that c[CONEHEIGHT]==1 for Gaussian random fields !!
  // However, for max-stable random fields s->height will have 
  // different values !!
  cone_height_socle = s->c[CONEHEIGHTSOCLE];
  cone_a            = s->c[CONEA];
  cone_b_socle      = s->c[CONEBSOCLE];
  *model = cone_model;
  // logically the current effectiveRadius should be set here
  // no need as here constant scale, so 
  // s->effectiveRadius = 0.5 * s->param[scale];
  // can be used.  

  // last: random point
  for (d=0; d<MPP_DIM; d++) {
    MPP_X[d]= s->min[d] - s->addradius + 
      (s->length[d] + 2.0 * s->addradius) * UNIFORM_RANDOM;
    min[d] = MPP_X[d] - s->effectiveRadius;
    max[d] = MPP_X[d] + s->effectiveRadius;
  }
}


#define GAUSSRADIUSSQ 0
double gaussInt(int d, int xi, double sigma, double R) {
  // int_{b(0,R) e^{-a r^2} dr = d b_d int_0^R r^{d-1} e^{-a r^2} dr
  // where a = 2.0 * xi / sigma^2
  // here : 2.0 is a factor so that the mpp function leads to the
  //            gaussian covariance model exp(-x^2)
  //        xi : 1 : integral ()
  //             2 : integral ()^2
  double a;
  a = 2.0 * xi / (sigma * sigma);
  switch (d) {
  case 1 : 
    double sqrta;
    sqrta = sqrt(a);
    return SQRTPI / sqrta * (2.0 * pnorm(SQRT2 * sqrta * R, 0, 1, 1, 0) - 1.0);
  case 2 :
    return PI / a * (1.0 - exp(- a * R * R));
  case 3 :
    double piDa;
    piDa = PI / a;
    return -2.0 * piDa * R * exp(- a * R * R )
      + piDa * sqrt(piDa) * (2.0 * pnorm(SQRT2 * sqrt(a) * R, 0, 1, 1, 0) - 1.0);
  default : assert(false);
  }
}

void gaussmpp_init(mpp_storage *s, int truedim, param_type param) 
{
  int d;

  s->dim = truedim;
  s->c[GAUSSRADIUSSQ] = -0.5 * log(MPP_APPROXZERO);
  s->effectiveRadius = sqrt(s->c[GAUSSRADIUSSQ]);
  if (s->addradius<=0.0) {s->addradius=s->effectiveRadius;}

  s->integral   = gaussInt(truedim, 1, 1.0, s->effectiveRadius);
  s->integralsq = gaussInt(truedim, 2, 1.0, s->effectiveRadius);

//  printf("%d %f %f %f\n", truedim, s->effectiveRadius, s->integral, s->integralsq);
//  assert(false);

  s->integralpos = s->integral;
  s->effectivearea = 1.0;
  for (d=0; d<truedim; d++) 
    s->effectivearea *= (s->length[d] + 2.0 * s->addradius);
  s->maxheight= 1.0;
}


void gaussmpp(mpp_storage *s, double *min, double *max, mppmodel *model)
{
  int d;

  MPP_DIM = s->dim;
  gauss_effectiveradiussq =  s->c[GAUSSRADIUSSQ];

  *model = gauss_model;  

  for (d=0; d<MPP_DIM; d++) {
    MPP_X[d]= s->min[d] - s->addradius + 
      (s->length[d] + 2.0 * s->addradius) * UNIFORM_RANDOM;
    min[d] = MPP_X[d] - s->effectiveRadius;
    max[d] = MPP_X[d] + s->effectiveRadius;
  }
}











