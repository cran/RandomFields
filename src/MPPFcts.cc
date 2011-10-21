/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de 

 function shapes for the random coin method in MPP.cc
 not included in CovFct.cc since these functions are of 
 particular structure

 Copyright (C) 2001 -- 2011 Martin Schlather, 

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


// basic shapes are: sphere (2 or 3 dimensions), cone, and Gaussian curve

#include <math.h>
 
#include "RF.h"
#include "Covariance.h"


void generalspherical_init(mpp_storage *s, cov_model *cov, int balldim) 
{
  int d, dim= cov->tsdim;
  double gamma,
    radius = 0.5;
  
  assert(s!=NULL);
  s->effectiveRadius = (s->plus < radius  && s->plus>0.0) ? s->plus : radius;
  
  for (d=dim; d<balldim; d++) {
    // uses that s->effectivearea has been set to 1.0 beforehand
    s->effectivearea *= (2.0 * s->effectiveRadius);
  }   
  gamma = gammafn(1.0 + 0.5 * balldim);
  s->r = radius;
  s->rsq = radius * radius;
  s->integral = s->integralsq = s->integralpos =
    pow(SQRTPI * radius, balldim) / gamma;
  s->maxheight = 1.0;
}

void mppinit_circular(mpp_storage *s, cov_model *cov) {
  int TRUESPACE=2;
  generalspherical_init(s, cov, TRUESPACE);
}

void mppinit_spherical(mpp_storage *s, cov_model *cov) {
  int TRUESPACE = 3;  
  generalspherical_init(s, cov, TRUESPACE);
}

void coin_circular(mpp_storage *s, cov_model *cov) { 
  // coin ist deterministisch; nur noch die Position in senkrechter Richtung
  // muss simuliert werden
  int d,
    dim = s->dim,
    TRUESPACE = 2;
  s->spherical_adddist = 0.0;
  s->location(s, cov);
  for (d=dim; d<TRUESPACE; d++) {
    double dummy;
    dummy = s->effectiveRadius * UNIFORM_RANDOM;
    s->spherical_adddist += dummy * dummy;
  }
}

void coin_spherical(mpp_storage *s, cov_model *cov){ 
  // coin is deterministisch; nur noch die Position in senkrechter Richtung
  // muss simuliert werden
  int d, 
    dim = s->dim,
    TRUESPACE = 3;  
  s->spherical_adddist = 0.0;
  s->location(s, cov);
  for (d=dim; d<TRUESPACE; d++) {
    double dummy;
    dummy = s->effectiveRadius * UNIFORM_RANDOM;
    s->spherical_adddist += dummy * dummy;
  }
}

res_type mppget_spherical(double *y, cov_model *cov, mpp_storage *s) {
  double distance = s->spherical_adddist + y[0];
//  int d,
  //   dim = s->dim;
  //for (d=0; d<dim; d++) {
    //   printf("%d %f %f\n", d, y[d], s->u[d]);
  //  double dummy = y[d] - s->u[d];
  // distance += dummy * dummy;
//  }
  return (res_type) ((distance >= s->r) ? 0.0 : 1.0);
}


#define CONER 1
#define CONERsq 2
#define CONEHEIGHTSOCLE 3
#define CONEA 4
#define CONEBSOCLE 5


void mppinit_cone(mpp_storage *s, cov_model *cov)
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
  int dim= cov->tsdim;
  double cr, cb, CR,
    r = cov->p[0][0],
    socle = cov->p[1][0],
    height = cov->p[2][0];
  assert(s!=NULL);
  cr= 0.5 * r;
  s->r = cr;
  s->rsq = cr * cr;
  CR  = 0.5;
  s->c[CONER] = CR;
  s->c[CONERsq] = CR * CR;
  s->c[CONEHEIGHTSOCLE] = socle + height;
  s->c[CONEA] =  /* a = h / (r - R) */
    height / ((r - 1.0) * 0.5);
  cb = height / (1.0 - r); /* b = h R / (R - r) */  
  s->c[CONEBSOCLE] = socle + cb;
  switch (dim) {
  case 2 :
    //
    s->integral = s->integralpos = 
      PI * (s->c[CONERsq] * socle +
	    0.333333333333333333 * height *
	    (s->rsq + s->c[CONERsq] + cr * CR));
    //
    s->integralsq 
      =  PI * (s->c[CONERsq] * socle * socle
	       + s->rsq * height * height 
	       + 0.5 * s->c[CONEA] * s->c[CONEA] *
	       (s->c[CONERsq] *  s->c[CONERsq] -
		s->rsq * s->rsq)
	       + 1.33333333333333333 *  s->c[CONEA] * cb *
	       (s->c[CONERsq] * CR - s->rsq * cr)
	       + cb * cb * (s->c[CONERsq] - s->rsq)
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
  // s->effectiveRadius = CR;

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

  s->effectiveRadius = s->plus <= 0.0 ? CR : s->plus;

  //
  s->maxheight = height + socle; //==s->c[CONEHEIGHTSOCLE]
}



// note. The cone_* values have to be  reset all the time,
// since the function does not know whether there have
// been other simulations of random fields, as cone_init is called
// only at the very fist time.
// Furthermore, the below construction will match with the
// general construction that randomises parameters (e.g. scale mixtures)
void coin_cone(mpp_storage *s, cov_model *cov)
{ 
  // note that c[CONEHEIGHT]==1 for Gaussian random fields !!
  // However, for max-stable random fields s->height will have 
  // different values !!

  // logically the current effectiveRadius should be set here
  // no need as here constant scale, so 
  // s->effectiveRadius = 0.5 * s->param[scale];
  // can be used.  
    s->location(s, cov);

}

/* y = ax + b, for the slope */
res_type mppget_cone(double *y, cov_model *cov, mpp_storage *s) {
  double distance = *y;
//  int d, 
//    dim = s->dim;
//  for (d=0; d<dim; d++) {
//    double dummy = y[d] - s->u[d];
//    distance += dummy * dummy;
//  }
  if (distance >= s->c[CONER]) return 0;
  if (distance <= s->r) return s->c[CONEHEIGHTSOCLE];
  return (res_type) (s->c[CONEBSOCLE] + s->c[CONEA] * distance);
}


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
  default :   
    error("dimension of gauss integral out of range");
    return RF_NAN;

  }
}

void mppinit_Gauss(mpp_storage *s, cov_model *cov){
  int dim= cov->tsdim;
  s->rsq = -0.5 * s->logapproxzero; //  -0.5 * log(approxzero)
  s->r = sqrt(s->rsq);
  if (s->effectiveRadius <= 0.0) s->effectiveRadius = s->r;

  s->integralpos = s->integral = gaussInt(dim, 1, 1.0, s->effectiveRadius);
  s->integralsq = gaussInt(dim, 2, 1.0, s->effectiveRadius);

  s->maxheight= 1.0;

//  printf("%f %f %f %f %f\n", 
//	 s->rsq, s->effectiveRadius, s->integralpos,
//	 s->integralsq);

}


void coin_Gauss(mpp_storage *s, cov_model *cov) { 
  // coin ist deterministisch.
    s->location(s, cov);
}


res_type mppget_Gauss(double *y, cov_model *cov, mpp_storage *s) {
  double distance = *y;
//  int d,
//    dim = s->dim;
//  distance = *y; // 0.0;
  //for (d=0; d<dim; d++) {
  //  double dummy = y[d] - s->u[d];
  //  distance += dummy * dummy;
  //} 
  // printf("dist=%f\n", distance);

 if (distance > s->r) return 0;
 return (res_type) exp(- 2.0 * distance * distance);
}





double gausstestInt(int d, int xi, double sigma, double R) {
  // int_{b(0,R) e^{-a r^2} dr = d b_d int_0^R r^{d-1} e^{-a r^2} dr
  // where a = 2.0 * xi / sigma^2
  // here : 2.0 is a factor so that the mpp function leads to the
  //            gaussian covariance model exp(-x^2)
  //        xi : 1 : integral ()
  //             2 : integral ()^2
  double a;

  // THIS IS NOT CORRECT YET !!!

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
  default :  
    error("dimension of gauss integral out of range");
    return RF_NAN;
  }
}

void mppinit_Gausstest(mpp_storage *s, cov_model *cov){
  int dim= cov->tsdim;
  s->rsq = -0.5 * s->logapproxzero; //  -0.5 * log(approxzero)
  if (s->effectiveRadius <= 0.0) // s->effectiveRadius =  sqrt(s->rsq);
      s->effectiveRadius = 35;

  s->integralpos = s->integral = gausstestInt(dim, 1, 1.0, s->effectiveRadius);
  s->integralsq = gausstestInt(dim, 2, 1.0, s->effectiveRadius);

  s->maxheight= 1.0;
}


void coin_Gausstest(mpp_storage *s, cov_model *cov) { 
  // coin ist deterministisch.
    s->location(s, cov);
}


res_type mppget_Gausstest(double *y, cov_model *cov, mpp_storage *s) {
    // nur stationaer, nicht isotrop
  double distance = 0.0;
  int d,
      dim = s->dim,
      time = dim - 1;
  double deltaT = y[time]; // y[time] - s->u[time];
  for (d=0; d<dim-1; d++) {
      double dummy = y[d] - deltaT; //y[d] - s->u[d] - deltaT;
    distance += dummy * dummy;
  } 
  distance += deltaT * deltaT; 
  if (distance > s->rsq) return 0;
  return (res_type) exp(-distance);
}








double whittleInt(int d, int xi, double sigma, double R) {
  // int_{b(0,R) e^{-a r^2} dr = d b_d int_0^R r^{d-1} e^{-a r^2} dr
  // where a = 2.0 * xi / sigma^2
  // here : 2.0 is a factor so that the mpp function leads to the
  //            whittleian covariance model exp(-x^2)
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
  default : 
    error("dimension of whittle integral out of range");
    return RF_NAN;
  }
}

void mppinit_Whittle(mpp_storage *s, cov_model *cov){
  int dim= cov->tsdim;
  s->rsq = -0.5 * s->logapproxzero; //  -0.5 * log(approxzero)
  s->r = sqrt(s->rsq);
  if (s->effectiveRadius <= 0.0) s->effectiveRadius = s->r;

  s->integralpos = s->integral = whittleInt(dim, 1, 1.0, s->effectiveRadius);
  s->integralsq = whittleInt(dim, 2, 1.0, s->effectiveRadius);

  s->maxheight= 1.0;
}


void coin_Whittle(mpp_storage *s, cov_model *cov) { 
  // coin ist deterministisch.
    s->location(s, cov);
}


res_type mppget_Whittle(double *y, cov_model *cov, mpp_storage *s) {
  double distance;
//  int d,
//    dim = s->dim;
  distance = *y; // 0.0;
//  for (d=0; d<dim; d++) {
//    double dummy = y[d] - s->u[d];
//    distance += dummy * dummy;
//  } 
  if (distance > s->r) return 0;
  return (res_type) exp(- 2.0 * distance * distance);
}








