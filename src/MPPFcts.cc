/*
 Authors 
 Martin Schlather, schlath@hsu-hh.de 

 *********** this function is still under construction *********

 Copyright (C) 2001 -- 2004 Martin Schlather, 

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

static int mpp_dim;
static double mpp_x[MAXDIM];


static double spherical_Rsq, spherical_adddist;
double spherical_model(double *y) 
{
  double distance;
  int i;
  distance = spherical_adddist;
  for (i=0; i<mpp_dim; i++) {
    register double dummy;
    dummy = y[i]-mpp_x[i];
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
  for (i=0; i<mpp_dim; i++) {
    register double dummy;
    dummy = y[i]-mpp_x[i];
    distance += dummy * dummy;
  }
  if (distance>=cone_Rsq) return 0;
  if (distance<=cone_rsq) return cone_height_socle;
  return cone_b_socle + cone_a * sqrt(distance);
}

static double gauss_invscale, gauss_effectiveradiussq;
double gauss_model(double *y)
{
  double distance;
  int i;
  distance = 0.0;
  for (i=0; i<mpp_dim; i++) {
    register double dummy;
    dummy = y[i]-mpp_x[i];
    distance += dummy * dummy;
  } 
  if (distance > gauss_effectiveradiussq) return 0;
  return exp(-gauss_invscale * distance);
}


#define SPHERICALRsq 1
#define SPHERICALadddist 2
void generalspherical_init(mpp_storage* s, int v, int balldim)
{
  int d;
  double gamma;

  assert(s!=NULL);
  s -> effectiveRadius[v] = 0.5 *  s->param[v][SCALE];
  if (s->addradius[v]<=0.0) {s->addradius[v] = s->effectiveRadius[v];}
  s -> effectivearea [v]= 1.0;
  for (d=0; d<s->timespacedim; d++) 
    s->effectivearea[v] *= (s->length[d] + 2.0 * s->effectiveRadius[v]);
  for (d=s->timespacedim; d<balldim; d++) {
    s->effectivearea[v] *= (2.0 * s->effectiveRadius[v]);
  }   
  gamma = gammafn(1.0+0.5*balldim);
  s->c[v][SPHERICALRsq] = s->effectiveRadius[v] * s->effectiveRadius[v];
  s->integral[v] = s->integralsq[v] = s->integralpos[v] =
    pow(SQRTPI * s->effectiveRadius[v], balldim) / gamma;
  s->maxheight[v] = 1.0;
}

void circular_init(mpp_storage *s, int v)
{
  int TRUESPACE=2;
  generalspherical_init(s, v, TRUESPACE);
}

void circularMpp(mpp_storage *s, int v, double *min, double *max,
	  mppmodel *model )
{ 
  int TRUESPACE = 2;  
  int d;
  assert(s!=NULL);

  mpp_dim = s->timespacedim;
  spherical_Rsq = s->c[v][SPHERICALRsq];
  *model = spherical_model;

  for (d=0; d<s->timespacedim; d++) {
    mpp_x[d]= s->min[d] - s->addradius[v] + 
      (s->length[d] + 2.0 * s->addradius[v]) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - s->effectiveRadius[v];
    max[d] = mpp_x[d] + s->effectiveRadius[v];
  }

  spherical_adddist = 0.0;
  for (d=s->timespacedim; d<TRUESPACE; d++) {
    double dummy;
    dummy = s->effectiveRadius[v] * UNIFORM_RANDOM;
    spherical_adddist += dummy * dummy;
  }
}

void spherical_init(mpp_storage *s, int v)
{
  int TRUESPACE = 3;  
  generalspherical_init(s, v, TRUESPACE);
}

void sphericalMpp(mpp_storage *s, int v, double *min, double *max,
	  mppmodel *model )
{ 
  int TRUESPACE = 3;  
  int d;
  mpp_dim = s->timespacedim;
  spherical_Rsq = s->c[v][SPHERICALRsq];
 *model = spherical_model;

  for (d=0; d<s->timespacedim; d++) {
    mpp_x[d]= s->min[d] - s->addradius[v] + 
      (s->length[d] + 2.0 * s->addradius[v]) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - s->effectiveRadius[v];
    max[d] = mpp_x[d] + s->effectiveRadius[v];
  }
  spherical_adddist = 0.0;
  for (d=s->timespacedim; d<TRUESPACE; d++) {
    double dummy;
    dummy = s->effectiveRadius[v] * UNIFORM_RANDOM;
    spherical_adddist += dummy * dummy;
  }
}

#define CONErsq 1
#define CONERsq 2
#define CONEHEIGHTSOCLE 3
#define CONEA 4
#define CONEBSOCLE 5
void cone_init(mpp_storage *s, int v)
{
  // ignoring: p[MEAN], p[NUGGET]
  // note : sqrt(p[VARIANCE]) -> height
  //        p[SCALE]          -> R
  // p[KAPPAI] = r/R \in [0,1]
  // p[KAPPAII]= socle (height)
  // p[KAPPAIII]=height of cone (without socle)
  // 
  // standard cone has radius 1/2, height of (potential) top: 1
  //           
  int d;
  double cr,cb,CR,socle,height;
  assert(s!=NULL);
  socle = s->param[v][KAPPAII];
  height= s->param[v][KAPPAIII];
  cr= 0.5 * s->param[v][KAPPAI] * s->param[v][SCALE];
  s->c[v][CONErsq] = cr * cr;
  CR  = 0.5 * s->param[v][SCALE];
  s->c[v][CONERsq] = CR * CR;
  s->c[v][CONEHEIGHTSOCLE] = socle + height;
  s->c[v][CONEA] =  /* a = h / (r - R) */
    height / ((s->param[v][KAPPAI] - 1.0) * 0.5* s->param[v][SCALE]);
  cb = height / (1.0 - s->param[v][KAPPAI]); /* b = h R / (R - r) */  
  s->c[v][CONEBSOCLE] = socle + cb;
  switch (s->timespacedim) {
  case 2 :
    //
    s->integral[v] = s->integralpos[v] = 
      PI * (s->c[v][CONERsq] * socle +
	    0.333333333333333333 * height *
	    (s->c[v][CONErsq] + s->c[v][CONERsq] + cr * CR));
    //
    s->integralsq[v] 
      =  PI * (s->c[v][CONERsq] * socle * socle
	       + s->c[v][CONErsq] * height * height 
	       + 0.5 * s->c[v][CONEA] * s->c[v][CONEA] *
	       (s->c[v][CONERsq] *  s->c[v][CONERsq] -
		s->c[v][CONErsq] * s->c[v][CONErsq])
	       + 1.33333333333333333 *  s->c[v][CONEA] * cb *
	       (s->c[v][CONERsq] * CR - s->c[v][CONErsq] * cr)
	       + cb * cb * (s->c[v][CONERsq] - s->c[v][CONErsq])
	       );
    //
    break;
    // case 3 : s -> integral = 4 PI / 3 * (R^3 * socle + r^3 * height)
    //                         + 4 PI int_r^R r^2 a x + b) dx
    //      s -> integralsq = 4 PI / 3 * (R^3 * socle *socle
    //                                    + r^3 * height * height)
    //                     + 4 PI int_r^R r^2 (a x + b)^2 dx
  default : assert(false);
  }

  // here a "mean" effectiveRadius is defined + true radius for simulation;
  // this is important
  s->effectiveRadius[v] = CR;
  if (s->addradius[v]<=0.0) {s->addradius[v] = s->effectiveRadius[v];}

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

  s->effectivearea[v] = 1.0;
  for (d=0; d<s->timespacedim; d++) 
     s->effectivearea[v]*= (s->length[d] + 2.0 * s->addradius[v]);

  //
  s->maxheight[v] = height + socle; //==s->c[CONEHEIGHTSOCLE]
}



// note. The cone_* values have to be  reset all the time,
// since the function does not know whether there have
// been other simulations of random fields, as cone_init is called
// only at the very fist time.
// Furthermore, the below construction will match with the
// general construction that randomises parameters (e.g. scale mixtures)
void cone(mpp_storage *s, int v, double *min, double *max,
	  mppmodel *model )
{ 
  int d;
  mpp_dim = s->timespacedim;
  cone_Rsq = s->c[v][CONERsq];
  cone_rsq = s->c[v][CONErsq]; 

  // note that c[CONEHEIGHT]==1 for Gaussian random fields !!
  // However, for max-stable random fields s->height will have 
  // different values !!
  cone_height_socle = s->c[v][CONEHEIGHTSOCLE];
  cone_a            = s->c[v][CONEA];
  cone_b_socle      = s->c[v][CONEBSOCLE];
  *model = cone_model;
  // logically the current effectiveRadius should be set here
  // no need as here constant scale, so 
  // s -> effectiveRadius = 0.5 * s->param[v][SCALE];
  // can be used.  

  // last: random point
  for (d=0; d<mpp_dim; d++) {
    mpp_x[d]= s->min[d] - s->addradius[v] + 
      (s->length[d] + 2.0 * s->addradius[v]) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - s->effectiveRadius[v];
    max[d] = mpp_x[d] + s->effectiveRadius[v];
  }
}
int checkcone(double *param, int timespacedim, SimulationType method) {
  switch (method) {
  case  AdditiveMpp: 
    if (timespacedim!=2) { 
      strcpy(ERRORSTRING_OK,"2 dimensional space");
      sprintf(ERRORSTRING_WRONG,"dim=%d",timespacedim);
      return ERRORCOVFAILED;
    }
    break;
  case Nothing :
    break;
  default :   
    strcpy(ERRORSTRING_OK,"method=\"add. MPP (random coins)\"");
    strcpy(ERRORSTRING_WRONG,METHODNAMES[method]);
    return ERRORCOVFAILED;
  }
  if ((param[KAPPAI]  >= 1.0)  ||
      (param[KAPPAI]  < 0.0)  ||
      (param[KAPPAII] < 0.0) ||
      (param[KAPPAIII]< 0.0) || 
      ((param[KAPPAIII] + param[KAPPAII])==0.0)) {
    strcpy(ERRORSTRING_OK,
	   "all kappa_i non-negative, kappa1<1, and kappa2+kappa3>0.0");
    sprintf(ERRORSTRING_WRONG,"c(%f,%f,%f)",param[KAPPAI],
	    param[KAPPAII],param[KAPPAIII]);
    return ERRORCOVFAILED;
  }	
  return NOERROR;
}
SimulationType methodcone(int spacedim, bool grid){
  return AdditiveMpp;
}
#define RANGE_EPSILON 1E-20
void rangecone(int spatialdim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1;
  double r[12] = {0, 1, 0, 1.0-RANGE_EPSILON,
		0, RF_INF, RANGE_EPSILON, 10.0, 
		0, RF_INF, RANGE_EPSILON, 10.0};
  memcpy(range, r, sizeof(double) * 12);
}
void infocone(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}


#define GAUSSINVSCALE 0
#define GAUSSRADIUSSQ 1

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
      + piDa * sqrt(piDa) * (2.0 * pnorm(SQRT2 * sqrt(a) * R, 0, 1,1,0) - 1.0);
  default : assert(false);
  }
}

void gaussmpp_init(mpp_storage *s, int v) 
{
  int d;
  s -> c[v][GAUSSINVSCALE] = 2.0 * s->param[v][INVSCALE] *
    s->param[v][INVSCALE]; 
  // e^{-2 x^2} !

  s -> c[v][GAUSSRADIUSSQ] =  - 0.5  * log(MPP_APPROXZERO) 
    * s->param[v][SCALE] * s->param[v][SCALE];
  s -> effectiveRadius[v] = sqrt(s -> c[v][GAUSSRADIUSSQ]);
  if (s->addradius[v]<=0.0) {s->addradius[v]=s->effectiveRadius[v];}

  s->integral[v]   = gaussInt(s->timespacedim, 1, s->param[v][SCALE],
			    s->effectiveRadius[v]);
  s->integralsq[v] = gaussInt(s->timespacedim, 2, s->param[v][SCALE], 
			    s->effectiveRadius[v]);

  s->integralpos[v] = s->integral[v];
  s->effectivearea[v] = 1.0;
  for (d=0; d<s->timespacedim; d++) 
     s->effectivearea[v] *= (s->length[d] + 2.0 * s->addradius[v]);
  s->maxheight[v]= 1.0;
}

void gaussmpp(mpp_storage *s, int v,double *min, double *max,
	  mppmodel *model )
{
  int d;

  gauss_invscale = s->c[v][GAUSSINVSCALE];
  gauss_effectiveradiussq =  s -> c[v][GAUSSRADIUSSQ];

  *model = gauss_model;  
  mpp_dim = s->timespacedim;

  for (d=0; d<mpp_dim; d++) {
    mpp_x[d]= s->min[d] - s->addradius[v] + 
      (s->length[d] + 2.0 * s->addradius[v]) * UNIFORM_RANDOM;
    min[d] = mpp_x[d] - s->effectiveRadius[v];
    max[d] = mpp_x[d] + s->effectiveRadius[v];
  }
}










