/*
 Authors 
 Martin Schlather, martin.schlather@cu.lu 

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * when defining your own function, make sure that the covariance function 
   itself allows for an additional nugget effect (spectral measures and tbm 
   operators don't)     
 * VARIANCE, SCALE should not be used

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

/* 
   MAKE SURE IN CHECK OF *SPATIALLY* ISOTROPIC (but not fully 
   isotropic) FUNCTIONS THAT THEY ARE ALWAYS CALLED BY A TIME AND A SPACE 
   COMPONENT, i.e. x is always a vector of two components !!

   dim must always be the true dimension!
*/

#define MINUSINVLOG005 0.3338082006953340674649
#define SQRTINVLOG005 0.5777613700268771079749
#define RANGE_EPSILON 1E-20

#include <math.h>
#include <assert.h>
#include "RFsimu.h"
#include "RFCovFcts.h"

//  {min, max} x {theor, pract}  x length(param) 
static Real range_stable[4] = {0, 2, 0.06, 2};
static Real range_whittle[4]= {0, RF_INF, 1e-2, 10.0};
static Real range_cauchy[4] = {0, RF_INF, 0.09, 10.0};
static Real range_genCauchy[8] = {0, 2, 0.05, 2, 0, RF_INF, 0.05, 10.0};


Real interpolate(Real y, Real *stuetz, int nstuetz, int origin,
		 Real lambda, int delta)
{
  int index,minindex,maxindex,i;
  double weights,sum,diff,a;
  index  = origin + (int) y;
  minindex = index - delta; if (minindex<0) minindex=0;
  maxindex = index + 1 + delta; if (maxindex>nstuetz) maxindex=nstuetz;
  weights = sum = 0.0;
  for (i=minindex;i<maxindex;i++) {    
    diff = y + (Real) (index-i);
    a    = exp( -lambda * diff * diff);
    weights += a;
    sum += a * stuetz[i];  // or  a/stuetz[i]
  }
  return (Real) (weights/sum); // and then   sum/weights       
}

//// NOTE : `*p' may not be changed by any of the functions!

SimulationType methodnugget(int spacedim, bool grid) {return Nugget;}


Real testCov(Real *x,Real *p, int effectivedim){
  Real y;
  y = fabs(*x);  
  if (y==0) return 1.0;
  return (y * M_E < 1) ? (1.0 + 1.0 / log(y)) : 0;
}
SimulationType methodtest(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return SpectralTBM;
    case 2 : return SpectralTBM;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}

/* exponential model */
Real exponential(Real *x,Real *p, int effectivedim){
  return exp(-fabs( *x));
}
Real Scaleexponential(Real *p,int scaling){ return MINUSINVLOG005; } 
Real TBM2exponential(Real *x, Real *p, int effectivedim) 
{
  double y;
  if (*x==0.0) {return 1.0;}
  y = fabs(*x);
  return 1.0 - PIHALF * y * I0mL0(y);
}
Real TBM3exponential(Real *x, Real *p, int effectivedim){
   register Real y;
   y = fabs( *x);
   return (1.0-y)*exp(-y);
}
Real Dexponential(Real *x, Real *p, int effectivedim){
  return - exp(-fabs( *x));
}
Real spectralexponential(Real *p ) { /* see Yaglom ! */
  Real register y;
  y = 1.0 - UNIFORM_RANDOM;
  return sqrt(1.0 / (y * y) - 1.0);
}
// spectral measure, see Lantuejoul, Geostatistical Simulation (2002)
SimulationType methodexponential(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return SpectralTBM;
    case 2 : return SpectralTBM;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangeexponential(int dim, int *index, Real* range){ 
  *index = -1; 
}


/* derivative of exponential: qexponential */
Real qexponential(Real *x,Real *p, int effectivedim){
  register Real y;
  y = exp(-fabs( *x));
  return y * (2.0  - p[KAPPA] * y) / (2.0 - p[KAPPA]);
}
Real Scaleqexponential(Real *p,int scaling){
  return -1.0 / 
    log( (1.0 - sqrt(1.0 - p[KAPPA] * (2.0 - p[KAPPA]) * 0.05)) / p[KAPPA]);
} 
int checkqexponential(Real *param, int timespacedim, SimulationType method){
  if ((param[KAPPA]<0.0) || (param[KAPPA]>1.0)) {
    strcpy(ERRORSTRING_OK,"kappa in [0,1]");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return 0;  
}
Real TBM3Dexponential(Real *x, Real *p, int effectivedim){
   register Real y;
   y = fabs( *x);
   return (y * (2.0  - p[KAPPA] * y) + y * (y * (p[KAPPA] * y - 1.0) * 2.0)) / 
     (2.0 - p[KAPPA]);
}
Real Dqexponential(Real *x,Real *p, int effectivedim){
  register Real y;
  y = exp(-fabs( *x));
  return y * (p[KAPPA] * y - 1.0) * 2.0 / (2.0 - p[KAPPA]);
}
SimulationType methodqexponential(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return Direct;
    case 2 : return Direct;
    default: return Direct;
    }
  }
}
void rangeqexponential(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  Real r[4] = {0, 1, 0, 1};
  memcpy(range, r, sizeof(Real) * 4);
}

/* derivative of exponential: damped cosine */
Real dampedcosine(Real *x, Real*p, int effectivedim){
  register double y;
  y = fabs( *x);
  return exp(-y * p[KAPPA]) * cos(y);
}
Real Scaledampedcosine(Real *p,int scaling){ 
  if (scaling==NATSCALE_EXACT) return 0.0;
  return MINUSINVLOG005; 
} 
Real TBM3dampedcosine(Real *x, Real *p, int effectivedim){
  register double y;
  y = fabs( *x);
  return exp(-p[KAPPA]*y) * ((1.0 - p[KAPPA] * y) * cos(y) - y * sin(y));
}
Real Ddampedcosine(Real *x, Real *p, int effectivedim){
  register double y;
  y = fabs( *x);
  return - exp(-p[KAPPA]*y) * (p[KAPPA] * cos(y) + sin(y));
}
int checkdampedcosine(Real *param, int timespacedim, SimulationType method){
  if ((timespacedim==3) || ((method==TBM3) && (timespacedim<=3))){
    if (param[KAPPA]<1.73205080756889) {
      strcpy(ERRORSTRING_OK,
	     "kappa >= sqrt(3) for 3-dimensional simulations and techniques");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  } else if (timespacedim==2){
    if (param[KAPPA]<1.0) {
      strcpy(ERRORSTRING_OK,"kappa >= 1.0 for 2 dimensions");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  } else if (timespacedim==1){
    if (param[KAPPA]<0.0) {
      strcpy(ERRORSTRING_OK,"kappa >= 0.0 for 1 dimension");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  } else { // timespacedim ==4
    return ERRORNOTPROGRAMMED;
  }
  return 0;  
}
SimulationType methoddampedcosine(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangedampedcosine(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  range[1] = RF_INF;
  range[3] = 10.0;
  switch(dim) {
      case 1 : range[2] = (range[0] = 0.0) + 1e-10; break;
      case 2 : range[2] = (range[0] = 1.0) + 1e-10; break;
      case 3 : range[2] = (range[0] = RF_M_SQRT_3) + 1e-10; break;
      default: *index = -2; range[0] = range[2] = RF_INF;
  }
}


/* Class of "spherical" models  */ 
/* circular model */
Real circular(Real *x, Real *p, int effectivedim)
{
  Real y;
  if ((y=fabs( *x))>1.0) return 0.0; 
  return  1.0 - (2.0 * (y * sqrt(1.0- y * y) + asin(y))) * INVPI;
}
Real Scalecircular(Real *p,int scaling) {return 1.138509531721630274603;}
// spectral measure, see Lantue !! 
SimulationType methodcircular(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return Forbidden;
    }
  } else {
    switch(spacedim) {
    case 1 : return Direct;
    case 2 : return Direct; 
    default: return Forbidden;
    }
  }
}
Real Dcircular(Real *x, Real *p, int dim){
  register Real y;
  if ((y=*x * *x) >= 1.0) {return 0.0;} 
  return -4 * INVPI * sqrt(1 - y);
}
void rangecircular(int dim, int *index, Real* range){
  if (dim<=2) *index=-1; else *index=-2;
}


/* spherical model */ 
Real spherical(Real *x, Real *p, int effectivedim){
  register Real y;
  if ((y= fabs( *x))>1.0) {return 0.0;}
  return (1.0+y*0.5*(y*y-3.0));
}
Real Scalespherical(Real *p,int scaling){ return 1.23243208931941;}
Real TBM2spherical(Real *x, Real *p, int effectivedim){
  register Real y, y2;
  y=fabs( *x); y2=y*y;
  if (effectivedim <= 2) {
    if (y>1.0)
      {return (1.0- 0.75 * y * ((2.0 - y2) * asin(1.0/y) + 
				sqrt(y2 -1.0)));} 
    return (1.0 - 0.375 * PI * y * (2.0 - y2));
  } else {
    assert(false);
  }
}
Real TBM3spherical(Real *x, Real *p, int effectivedim){
  register Real y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return (1.0 + (-3.0 + 2.0 * y * y) * y);
}
Real Dspherical(Real *x, Real *p, int effectivedim){
  register Real y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return  1.5 * (y * y - 1.0);
}
// spectral measure, see Lantue !! 
SimulationType methodspherical(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    case 3 : return CircEmbed;
    default: return Forbidden;
   }
  } else {
    switch(spacedim) {
    case 1 : return TBM2;
    case 2 : return TBM2;
    case 3 : return TBM3;
    default: return Forbidden;
    }
  }
}
void rangespherical(int dim, int *index, Real* range){
  if (dim<=3) *index=-1; else *index=-2;
}


/* power model */ 
Real power(Real *x, Real *p, int effectivedim){
  register Real y;
  if ((y=fabs( *x))>1.0) {return 0.0;}
  return  pow(1.0 - y, p[KAPPA]);
}
Real Scalepower(Real *p,int scaling){ 
  return 1.0 / (1.0 - pow(0.05,1/p[KAPPA]));
}
Real TBM2power(Real *x, Real *p, int effectivedim){
  // only kappa=2 up to now !
  register Real y;
  y=fabs( *x);
  if (y>1.0)
    {return (1.0 - 2.0 * y *(asin(1.0/y) - y + sqrt(y*y-1.0) ));} 
  return (1.0 - y * (PI - 2.0 * y));
}
Real TBM3power(Real *x, Real *p, int effectivedim){
  register Real y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return  (1.0 - y - y * p[KAPPA]) *
    pow(1.0 - y, p[KAPPA]-1.0);
}
Real Dpower(Real *x, Real *p, int effectivedim){
  register Real y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return  -  p[KAPPA] * pow(1.0 - y, p[KAPPA]-1.0);
}
int checkpower(Real *param, int timespacedim, SimulationType method) {
  switch (method) {
  case TBM2 : 
    if (param[KAPPA]!=2) { // not that kappa==2 if fine also for 
      //                              3 dimensions!
      strcpy(ERRORSTRING_OK,"kappa=2 in case of TBM2");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
    break;
  case TBM3 : 
    if (param[KAPPA]<2) { 
      strcpy(ERRORSTRING_OK,"kappa>=2");
      sprintf(ERRORSTRING_WRONG,"method=TBM3 and kappa=%f",
	      param[KAPPA]);
      return ERRORCOVFAILED;
      }
    if (timespacedim<=3) break; //!!
  default: 
    if (param[KAPPA] <  0.5 * (1.0 + timespacedim)) {
      strcpy(ERRORSTRING_OK,"kappa >= (dim+1)/2");
      sprintf(ERRORSTRING_WRONG,"kappa=%f dim=%d",param[KAPPA],
	      timespacedim);
      return ERRORCOVFAILED;
    }
  }
  return 0;
}
SimulationType methodpower(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
// range definition:
// 0: min, theory, 1:max, theory
// 2: min, practically 3:max, practically
void rangepower(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  range[0] = range[2] = 0.5 * ((Real) dim + 1);
  range[1] = RF_INF;
  range[3] = 20.0;
}

/* stable model */
Real stable(Real *x,Real *p, int effectivedim){
  if (*x==0.0) {return 1.0;}
  return  (exp(-pow(fabs( *x),p[KAPPA])));
}
Real Scalestable(Real *p,int scaling){return pow(MINUSINVLOG005,1/p[KAPPA]);}
Real TBM3stable(Real *x, Real *p, int effectivedim){
  register double y;
  y = pow(fabs( *x),p[KAPPA]);
  return exp(-y) * (1 - p[KAPPA] * y);
}
Real Dstable(Real *x, Real *p, int effectivedim){
  register double y,z;
  if ( (z = fabs( *x)) == 0.0) return 0.0; /* WRONG VALUE, but multiplied 
					      with 0 anyway*/
  y = pow(z,p[KAPPA] - 1.0);
  return -  p[KAPPA] * exp(-y * z);
}
int checkstable(Real *param, int timespacedim, SimulationType method) {
  if ((param[KAPPA]<=0) || (param[KAPPA]>2.0)) {
    strcpy(ERRORSTRING_OK,"0<kappa<=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED; 
  }
  return 0;
}
SimulationType methodstable(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangestable(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_stable, sizeof(Real) * 4);
}

/* Whittle-Matern or Whittle or Besset */ 
Real WhittleMatern(Real *x, Real *p, int effectivedim) 
// check calling functions, like hyperbolic and gneiting if any changings !!
{
  static double kappa=RF_INF;
  static double loggamma;
  register double y;
  if ( *x==0.0) {return 1.0;}
  y = fabs( *x); 
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    loggamma = lgammafn(kappa);
  }
  return 2.0 *
    exp(kappa*log(0.5*y)-loggamma + log(bessel_k(y,kappa,2.0))-y);
}
Real ScaleWhittleMatern(Real *p,int scaling){
  // it is working reasonably well if kappa is in [0.001,100]
  // happy to get any better suggestion!!
  static int nstuetz = 19;
  static Real stuetz[]=
  {1.41062516176753e-14, 4.41556861847671e-12, 2.22633601732610e-06, 
   1.58113643548649e-03, 4.22181082102606e-02, 2.25024764696152e-01,
   5.70478218148777e-01, 1.03102016706644e+00, 1.57836638352906e+00,
   2.21866372852304e+00, 2.99573229151620e+00, 3.99852231863082e+00,
   5.36837527567695e+00, 7.30561120838150e+00, 1.00809957038601e+01,
   1.40580075785156e+01, 1.97332533513488e+01, 2.78005149402352e+01,
   3.92400265713477e+01};
  static int stuetzorigin = 11;
  if (scaling==NATSCALE_EXACT) return 0.0;
  return interpolate(log(p[KAPPA]) * INVLOG2, stuetz, nstuetz, stuetzorigin, 
		     1.5, 5);
}

Real TBM2WhittleMatern(Real *x, Real *p, int effectivedim) 
{
  if (p[KAPPA]==0.5) return TBM2exponential(x, p, effectivedim);
  assert(false);
}

Real TBM3WhittleMatern(Real *x, Real *p, int effectivedim)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double kappa=RF_INF;
  static double loggamma;
  register double y,loghalfy;
  if ( *x==0.0) {return 1.0;}
  y = fabs( *x); 
  loghalfy = log(0.5*y);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    loggamma = lgammafn(kappa);
  }
  return  
    (2.0 * exp(kappa * loghalfy -loggamma + log(bessel_k(y, kappa, 2.0)) - y)
     -(4.0 * exp((kappa + 1.0) * loghalfy - loggamma + 
		 log(bessel_k(y, kappa - 1.0, 2.0)) - y)));
}  

#define LOG05 -0.69314718055994528623
Real DWhittleMatern(Real *x, Real *p, int effectivedim)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double kappa=RF_INF;
  static double loggamma;
  register double y;
  if ( *x==0.0) {return 1.0;}
  y = fabs( *x); 
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    loggamma = lgammafn(kappa);
  }
  return -2.0 *exp(kappa * log(0.5*y) - loggamma + 
		   log(bessel_k(y, kappa-1.0, 2.0)) - y);
   
}
Real spectralWhittleMatern(Real *p ) { /* see Yaglom ! */
  return sqrt(pow(1.0 - UNIFORM_RANDOM,-1.0/p[KAPPA]) - 1.0);
}
int checkWhittleMatern(Real *param, int timespacedim, SimulationType method) 
{ 
  static Real spectrallimit=0.17;
  switch(method) {
  case TBM2 : 
    if (param[KAPPA]!=0.5) return ERRORNOTPROGRAMMED;
    break;
  case SpectralTBM : 
    if ((param[KAPPA]>=spectrallimit)) return 0;
    sprintf(ERRORSTRING_OK,
	    "%f<=kappa (the numerical errors are too big for 0<kappa<%f)",
	    spectrallimit,spectrallimit);
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  default :
    if ((param[KAPPA]>0)) return 0;
    strcpy(ERRORSTRING_OK,"0<kappa");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return 0;
}
SimulationType methodWhittleMatern(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return SpectralTBM;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return SpectralTBM;
    case 2 : return SpectralTBM;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangeWhittleMatern(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_whittle, sizeof(Real) * 4);
}


Real hyperbolic(Real *x, Real*p, int effectivedim){ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double kappay;
  Real y;
  if ( *x==0.0) {return 1.0;}
  if (p[KAPPAIII]==0) { // whittle matern
    y = *x * p[KAPPAI];
    return WhittleMatern(&y, p, effectivedim);
  } 
  if (p[KAPPAI]==0) { //cauchy   => KAPPAII < 0 !
    y =  *x/p[KAPPAIII];
    /* note change in sign as KAPPAII<0 */
    return  pow(1+y*y,p[KAPPAII]); 
  }  

  y=fabs( *x); 
  if ((p[KAPPAI]!=kappa) || (p[KAPPAII]!=lambda) || (p[KAPPAIII]!=delta)) {
    kappa=p[KAPPAI]; 
    lambda= p[KAPPAII];
    delta = p[KAPPAIII];
    deltasq = delta * delta;
    kappadelta = kappa * delta;
    logconst = kappadelta - log(bessel_k(kappadelta,lambda,2.0)) 
      - lambda * log(delta);
  }
  y=sqrt(deltasq + y*y);  
  kappay = kappa * y;
  return  
    exp(logconst + lambda * log(y)+log(bessel_k(kappay,lambda,2.0))-kappay);
}
Real TBM3hyperbolic(Real *x, Real*p, int effectivedim)
{ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  Real y;
  double ysq,s,kappas,logs;
  if ( *x==0.0) {return 1.0;}
  if (p[KAPPAIII]==0) { // whittle matern
    y = *x * p[KAPPAI];
    return TBM3WhittleMatern(&y, p, effectivedim);
  } 
  if (p[KAPPAI]==0) { //cauchy
    register Real y,ha;
    y= *x/p[KAPPAIII];
    ha=y*y;
    /* note change in sign as KAPPAII<0 */
    return (1.0+ (1.0+2.0 * p[KAPPAII]) * ha) * 
      pow(1.0+ha,p[KAPPAII]-1.0);
  }
  y=fabs( *x); 
  if ((p[KAPPAI]!=kappa) || (p[KAPPAII]!=lambda) || (p[KAPPAIII]!=delta)) {
    kappa=p[KAPPAI]; 
    lambda= p[KAPPAII];
    delta = p[KAPPAIII];
    deltasq = delta * delta;
    kappadelta = kappa * delta;
    logconst = kappadelta - log(bessel_k(kappadelta,lambda,2.0)) 
      - lambda * log(delta);
  }
  ysq = y * y;
  s=sqrt(deltasq + ysq);
  kappas = kappa * s;
  logs = log(s);  
  return  
    ( exp(logconst + lambda * logs +log(bessel_k(kappas,lambda,2.0))-kappas)
      - ysq*kappa*exp(logconst + (lambda-1.0)*logs 
		      +log(bessel_k(kappas,lambda-1.0,2.0))-kappas)
      );
}
Real Dhyperbolic(Real *x, Real*p, int effectivedim)
{ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  Real y;
  double s,kappas,logs;
  if ( *x==0.0) {return 1.0;}
  if (p[KAPPAIII]==0) { // whittle matern
    y = *x * p[KAPPAI];
    return DWhittleMatern(&y, p, effectivedim);
  } 
  if (p[KAPPAI]==0) { //cauchy
    register Real y,ha;
    y= *x/p[KAPPAIII];
    ha=y*y;
    /* note change in sign as KAPPAII<0 */
    return  2.0 * p[KAPPAII] * fabs(y) * 
      pow(1.0+ha,p[KAPPAII]-1.0);
  }
  y=fabs( *x); 
  if ((p[KAPPAI]!=kappa) || (p[KAPPAII]!=lambda) || (p[KAPPAIII]!=delta)) {
    kappa=p[KAPPAI]; 
    lambda= p[KAPPAII];
    delta = p[KAPPAIII];
    deltasq = delta * delta;
    kappadelta = kappa * delta;
    logconst = kappadelta - log(bessel_k(kappadelta,lambda,2.0)) 
      - lambda * log(delta);
  }
  s=sqrt(deltasq + y * y);
  kappas = kappa * s;
  logs = log(s);  
  return  
    ( 
      - y * kappa*exp(logconst + (lambda-1.0)*logs 
		      +log(bessel_k(kappas,lambda-1.0,2.0))-kappas)
      );
}
int checkhyperbolic(Real *param, int timespacedim, SimulationType method){
  if (param[KAPPAII]>0) {
    if ((param[KAPPAIII]<0) || (param[KAPPAI]<=0)) {
      strcpy(ERRORSTRING_OK,"kappa1>0 and kappa3>=0 if kappa2>0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      param[KAPPAI],param[KAPPAIII],param[KAPPAII]);
      return ERRORCOVFAILED;
    }
  } else if (param[KAPPAII]<0) {
    if ((param[KAPPAIII]<=0) || (param[KAPPAI]<0)) {
      strcpy(ERRORSTRING_OK,"kappa1>=0 and kappa3>0 if kappa2<0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      param[KAPPAI],param[KAPPAIII],param[KAPPAII]);
      return ERRORCOVFAILED;
    }
  } else { // param[KAPPAII]==0.0
    if ((param[KAPPAIII]<=0) || (param[KAPPAI]<=0)) {
      strcpy(ERRORSTRING_OK,"kappa1>0 and kappa3>0 if kappa2=0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      param[KAPPAI],param[KAPPAIII],param[KAPPAII]);
      return ERRORCOVFAILED;
    }
  }
  return 0;
}
SimulationType methodhyperbolic(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangehyperbolic(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  Real r[12] = {0, RF_INF, 0.000001, 10.0,
		RF_NEGINF, RF_INF, -20.0, 20.0,
                0, RF_INF, 0.000001, 10.0
                };
  *index = -1; 
  memcpy(range, r, sizeof(Real) * 12);
}


/* Gneiting's functions */
// alternative to Gaussian */
#define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
#define NumericalScale 0.301187465825
Real Gneiting(Real *x, Real *p, int effectivedim){ 
  register Real y,oneMy8;
  if ((y=fabs(*x * NumericalScale))>1.0) {return(0.0);}
  oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
  return ((1.0+y * ( 8.0 + y * (25.0 + 32.0 *y)))*oneMy8);
}
Real ScaleGneiting(Real *p,int scaling) {return 0.5854160193;}
Real TBM3Gneiting(Real *x, Real *p, int effectivedim){ 
  register Real y,oneMy7;
  if ((y=fabs( *x*Sqrt2TenD47))>1.0) {return 0.0;}  
  oneMy7 = 1.0-y; oneMy7*=oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);
  return 
    (1.0 + y * (7.0  -  y * (5.0 + y * (147.0 + 384.0 * y))))* oneMy7;
}
Real DGneiting(Real *x, Real *p, int effectivedim){ 
  register Real y,oneMy7;
  if ((y=fabs( *x*Sqrt2TenD47))>1.0) {return 0.0;}  
  oneMy7 = 1.0-y; oneMy7*=oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);
  return 
    (-y) * ( 22.0 + y * (154.0 + y * 352.0)) * oneMy7;
}
SimulationType methodGneiting(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangeGneiting(int dim, int *index, Real* range){
  if (dim<=3) *index=-1; else *index=-2;
}


Real genGneiting(Real *x, Real *p, int effectivedim)
{
  register Real y, s;
  if ((y=fabs( *x))>1.0) return 0.0; 
  s = p[KAPPAII] + p[KAPPAI];
  switch ((int) p[KAPPAI]) {
  case 1:
    return  pow ((1.0-y), s) * (1.0 + s*y);   
  case 2:
    return  pow ((1.0-y), s) * 
      (1.0 + y * (s + y * (s*s-1.0)*0.3333333333333333));  
  case 3:
    return   pow ((1.0-y), s) * 
      (1.0 + y * (s 
		  + y * ( 0.2 * (2.0*s*s - 3.0)
			 + y * (s*s-4.0)*s*0.0666666666666666666))) ;
  default : assert(false);  
  }
}
Real TBM3genGneiting(Real *x, Real *p, int effectivedim)
{
  register Real y, s;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  s = p[KAPPAII] + p[KAPPAI];
  switch ((int) p[KAPPAI]) {
  case 1:
    return  pow((1.0 - y), p[KAPPAII]) *
      (1.0 + y * (p[KAPPAII] - y * s * (p[KAPPAII] + 3.0))); 
  case 2:
   return  pow((1.0 - y), (p[KAPPAII]+1.0)) *
      (1.0 + y * ((p[KAPPAII]+1.0)
		- y *(1.0 + 2.0 * s 
		      + y * (s * s - 1.0) * 0.33333333333333 * 
		      (p[KAPPAII] + 5.0))));
  case 3:
    return  pow((1.0 - y), (p[KAPPAII]+2.0)) *
      (1.0 + y * ((p[KAPPAII]+2.0)
		  + y * (0.2 * (-9.0 + s * (-10.0 + s))
			 + y * 0.0666666666666666 * 
			 ((27.0 - s*(7.0 + s*(18.0 + 2.0*s)))
			  - y * (s*s - 4.0) * s * (p[KAPPAII] + 7.0)))));
  default : assert(false);   
  }
}
Real DgenGneiting(Real *x, Real *p, int effectivedim)
{
  register Real y, s;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  s = p[KAPPAII] + p[KAPPA];
  switch ((int) p[KAPPAI]) {
  case 1:
    return  -  pow(1.0 - y, p[KAPPAII]) *
      y *  (p[KAPPAII] + 1.0)  * (p[KAPPAII] + 2.0) ; 
  case 2:
    return -  pow(1.0 - y, s - 1.0) *
      y * (0.333333333333333333 * s * s + 0.66666666666666666667 + s +
	   y * 0.333333333333333333 * (s * s - 1.0) * (p[KAPPAII] + 4.0) );
  case 3:
    return -  pow(1.0 - y, s - 1.0) *
      0.2 * y * (6 + s * (5.0 + s) +
		 y * (-6 + s * (1.0 + s * (4.0 + s)) +
		      - 0.33333333333333333 * y * s * (-12.0 + s * (-4.0 + s * (3.0 + s)))
		      ));
  default : assert(false);   
  }
}
int checkgenGneiting(Real *param, int timespacedim, SimulationType method)
{
  if ((param[KAPPAI] != (double) ((int)param[KAPPAI])) || 
      (param[KAPPAI]<0)) {
      strcpy(ERRORSTRING_OK,"positive integer kappa");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAI]);
      return ERRORCOVFAILED;
  }
  if ((param[KAPPAI] < 1.0) || (param[KAPPAI] > 3.0))
    return ERRORNOTPROGRAMMED;
  if ((method==TBM3) &&
      (param[KAPPAII] < 2.0 + param[KAPPAI])) {
    strcpy(ERRORSTRING_OK,"kappa2 <= 2.0 + kappa1");
    sprintf(ERRORSTRING_WRONG, "method=TBM3, kappa1=%f and kappa2=%f ",
	    param[KAPPAI], param[KAPPAII]);
    return ERRORCOVFAILED;    
  }
  if  (param[KAPPAII] < ( 0.5 * timespacedim + 
				 param[KAPPAI] + 0.5)) {
    strcpy(ERRORSTRING_OK, "kappa2 >= (dim + 2*kappa1 + 1)/2 ");
    sprintf(ERRORSTRING_WRONG, "kappa1=%f, kappa2=%f, and dim=%d ",
	    param[KAPPAI], param[KAPPAII], timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;
}
SimulationType methodgenGneiting(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangegenGneiting(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  range[0] = range[2] = range[1] = range[3] = *index;
  range[4] = range[6] = 0.5 * (Real) (dim + 2 * *index + 1);
  range[5] = RF_INF;
  range[7] = 20.0;
  if ((++(*index)) > 3) *index=-1;
  if (dim>3) *index=-2;
}



/* Gausian model */
Real Gauss(Real *x, Real*p, int effectivedim) {
  return  exp(-  *x * *x);
}
Real ScaleGauss(Real *p,int scaling) {return SQRTINVLOG005;}
Real TBM3Gauss(Real *x, Real*p, int effectivedim) {
  register Real y;
  y =   *x *  *x; 
  return  (1-2.0*y)*exp(-y);
}
Real DGauss(Real *x, Real*p, int effectivedim) {
  register Real y;
  y =   fabs( *x); 
  return  -  2.0 * y * exp(- y * y);
}
Real spectralGauss(Real *p ) {   
  return 2.0  * sqrt(-log(1.0 - UNIFORM_RANDOM));
}
SimulationType methodGauss(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return SpectralTBM;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return SpectralTBM;
    case 2 : return SpectralTBM;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangeGauss(int dim, int *index, Real* range){
  *index=-1;
}


/* Cauchy models */
Real Cauchy(Real *x, Real *p, int effectivedim){
  return  pow(1+ *x * *x,-p[KAPPA]);
}
Real ScaleCauchy(Real *p,int scaling) {
  switch(scaling) {
  case NATSCALE_EXACT: case NATSCALE_APPROX:
    return 1.0/sqrt(pow(0.05,-1/p[KAPPA])-1.0); 
    break;
  case NATSCALE_MLE:
    /// this values should be changed, because of 
    /// possibly long tails!
    return 1.0/sqrt(pow(0.05,-1/p[KAPPA])-1.0); 
    break;
  default: assert(false);
  }
}
Real TBM2Cauchy(Real *x, Real *p, int effectivedim){
  register Real y2, lpy2;
  y2 = *x * *x; 
  lpy2=1.0+y2;
  switch ((int) (p[KAPPA]*2.0 + 0.1)) {
  case 1 : return 1.0 / lpy2;
  case 3 : return  (1.0 - y2)/ (lpy2*lpy2);
  case 5 : return 
	     (1.0-y2*(2.0+0.333333333333333333333*y2))/(lpy2*lpy2*lpy2);
  case 7 : lpy2 *= lpy2;
    return (1.0- y2*(3.0+y2*(1.0+0.2*y2)))/(lpy2 * lpy2);
  default : assert(false);
  }
}
Real TBM3Cauchy(Real *x, Real *p, int effectivedim){
  register Real ha;
  ha= *x *  *x;
  return (1.0+ (1.0-2.0*p[KAPPA])*ha) * pow(1.0+ha,-p[KAPPA]-1.0);
}
Real DCauchy(Real *x, Real *p, int effectivedim){
  register Real y;
  y = fabs( *x);
  return (-2.0 * p[KAPPA] * y) * pow(1.0 + y * y, -p[KAPPA]-1.0);
}
int checkCauchy(Real *param, int timespacedim, SimulationType method){
  switch (method) {
  case TBM2 : 
    if (timespacedim>2) return ERRORCOVFAILED;
    if ((param[KAPPA]!=0.5) && (param[KAPPA]!=1.5) &&
	(param[KAPPA]!=2.5) && (param[KAPPA]!=3.5)) {
      strcpy(ERRORSTRING_OK,"kappa in {0.5, 1.5, 2.5 ,3.5} (TBM2)");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  default : 
    if (param[KAPPA]<=0) {
      strcpy(ERRORSTRING_OK,"0<kappa");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  }
  return 0;
}
SimulationType methodCauchy(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return TBM3;
    default: return TBM3;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangeCauchy(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_cauchy, sizeof(Real) * 4);
}

Real generalisedCauchy(Real *x, Real *p, int effectivedim){
  return pow(1.0 + pow(fabs(*x), p[KAPPAI]), -p[KAPPAII]/p[KAPPAI]);
}
Real ScalegeneralisedCauchy(Real *p,int scaling) {
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
Real TBM3generalisedCauchy(Real *x, Real *p, int effectivedim){
  register Real ha;
  ha=pow(fabs( *x),p[KAPPAI]);
  return 
    (1.0+ (1.0-p[KAPPAII])*ha) * pow(1.0+ha,-p[KAPPAII]/p[KAPPAI]-1.0);
}
Real DgeneralisedCauchy(Real *x, Real *p, int effectivedim){
  register Real ha,y;
  if ((y = fabs(*x))==0.0) return 0.0; // WRONG VALUE, but multiplied with 
                                       // zero anyway 
  ha=pow(y,p[KAPPAI] - 1.0);
  return  -  p[KAPPAII] * ha * pow(1.0 + ha * y,-p[KAPPAII]/p[KAPPAI]-1.0);
}
int checkgeneralisedCauchy(Real *param, int timespacedim, SimulationType method){
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
  return 0;
}
SimulationType methodgeneralisedCauchy(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return TBM3;
    default: return TBM3;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangegeneralisedCauchy(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_genCauchy, sizeof(Real) * 8);
}

Real Cauchytbm(Real *x, Real *p, int effectivedim){
  register Real ha;
  if ( *x==0) {return 1.0;}
  ha=pow(fabs( *x),p[KAPPAI]);
  return 
    (1.0+ (1.0-p[KAPPAII]/p[KAPPAIII])*ha) * pow(1.0+ha,-p[KAPPAII]/p[KAPPAI]-1.0);
}
Real TBM3Cauchytbm(Real *x, Real *p, int effectivedim){
  register Real bg,ha;
  ha=pow(fabs( *x),p[KAPPAI]);
  bg=p[KAPPAII]/p[KAPPAIII];
  return  
    (1 + ha * (1-bg*(1+p[KAPPAI])+(1-p[KAPPAII]) * (1+(1-bg)*ha)))*
    pow(1+ha,-p[KAPPAII]/p[KAPPAI]-2.0);
}
Real DCauchytbm(Real *x, Real *p, int effectivedim){
  register Real y,ha;
  if ((y = fabs(*x)) == 0.0) return 0.0; // WRONG VALUE, but multiplied 
  //                                        zero anyway
  ha=pow(y,p[KAPPAI] - 1.0);
  return  
    p[KAPPAII] *  ha * (-1.0 - p[KAPPAI]/p[KAPPAIII] + 
		       ha * y * (p[KAPPAII]/p[KAPPAIII] - 1.0)) *
    pow(1.0 + ha * y,-p[KAPPAII]/p[KAPPAI]-2.0);
}
int checkCauchytbm(Real *param, int timespacedim, SimulationType method){
  if ((method==TBM3) && (3.0 > param[KAPPAIII])) {
    strcpy(ERRORSTRING_OK,"3 <= kappa3 and method=TBM3");
    sprintf(ERRORSTRING_WRONG,
	    "method=TBM3 and kappa3=%f ",param[KAPPAIII]);
    return ERRORCOVFAILED;    
  }
  if (timespacedim > param[KAPPAIII]) {
    strcpy(ERRORSTRING_OK,"dim <= kappa3");
    sprintf(ERRORSTRING_WRONG,
	    "kappa3=%f < dim=%d ",param[KAPPAIII], timespacedim);
    return ERRORCOVFAILED;    
  }
  // theory: should check whether p[KAPPAIII] is an integer. 
  // practically, it goes through without this restriction !!
  if ((param[KAPPAI]<=0) || (param[KAPPAI]>2.0)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAI]);
    return ERRORCOVFAILED;
  }
  if (param[KAPPAII]<=0) {
    strcpy(ERRORSTRING_OK,"0<kappa2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAII]);
    return ERRORCOVFAILED;
  }
  return 0;
}
SimulationType methodCauchytbm(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return TBM3;
    default: return TBM3;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangeCauchytbm(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_genCauchy, sizeof(Real) * 8);
  range[8] = range[10] = (double) dim;
  range[9] = RF_INF;
  range[11] = (double) dim + 10.0;
}


/* Bessel function */ 
Real Bessel(Real *x,Real *p, int effectivedim){
  static double kappa=RF_INF;
  static double gamma;
  Real y;
  if  ( *x==0.0) {return 1.0;} 
  y = fabs( *x);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    gamma = gammafn(kappa+1.0);
  }
  return  gamma  * pow(0.5 * y,-kappa) * bessel_j(y,kappa);
}
Real spectralBessel(Real *p ) {
  return p[KAPPA]==0 ? 1.0 : 
    (sqrt(1-pow(UNIFORM_RANDOM,1.0/p[KAPPA])));
}
int checkBessel(Real *param,  int timespacedim, SimulationType method){
  // Whenever TBM3Bessel exists, add further check against too small kappa!
  if (param[KAPPA]<0.5*((double) timespacedim) - 1.0) {
    strcpy(ERRORSTRING_OK,"kappa >= dim/2 - 1");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
  if (method==SpectralTBM && param[KAPPA]<0.0) {
    strcpy(ERRORSTRING_OK,"spectral TBM and kappa>=0");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return 0;
}
SimulationType methodBessel(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return SpectralTBM;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return SpectralTBM;
    case 2 : return SpectralTBM;
    case 3 : return Direct;
    default: return Direct;
    }
  }
}
void rangeBessel(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract }
  *index = -1; 
  range[2] = 0.0001 + (range[0] =  0.5 * ((Real) dim - 2.0));
  range[1] = RF_INF;
  range[3] = 10.0;
}


Real wave(Real *x, Real *p, int effectivedim) {
  if ( *x==0.0) { return 1.0;}
  return sin( *x)/ *x;
}
Real Scalewave(Real *p,int scaling) {return 0.302320850755833;}
Real spectralwave(Real *p ) {
  Real x;  x=UNIFORM_RANDOM; return sqrt(1- x * x);
}
SimulationType methodwave(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return SpectralTBM;
    case 2 : return SpectralTBM;
    case 3 : return Direct;
    default: return Direct;
    }
  }
}
void rangewave(int dim, int *index, Real* range){
  if(dim<=3) *index=-1; else *index=-2;
}


Real cubic(Real *x, Real *p, int effectivedim)
{ ///
  Real y, y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return  (1.0 + (((0.75 * y2 - 3.5) * y2 + 8.75) * y - 7) * y2);
}
Real TBM3cubic(Real *x, Real *p, int effectivedim)
{ ///
  Real y, y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return (1.0 + y2 * (-21.0 + y * (35.0 + y2 * (-21.0 + 6.0 * y2))));
}
Real Dcubic(Real *x, Real *p, int effectivedim) 
{ ///
  Real y,y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return  y * (-14.0 + y * (26.25 + y2 * (-17.5 + 5.25 * y2)));
}
Real Scalecubic(Real *p,int scaling) {return 1.44855683156829;}
SimulationType methodcubic(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangecubic(int dim, int *index, Real* range){
  if(dim<=3) *index=-1; else *index=-2;
}


Real penta(Real *x, Real *p, int effectivedim)
{ ///
  Real y,y2;
  if ((y=fabs( *x))>=1.0) { return 0.0; }
  y2=y * y;
  return  
    (1.0 + y2 * (-7.333333333333333 
		 + y2 * (33.0 +
			 y * (-38.5 + 
			      y2 * (16.5 + 
				    y2 * (-5.5 + 
					  y2 * 0.833333333333333))))));
}
Real TBM3penta(Real *x, Real *p, int effectivedim)
{ ///
  Real y,y2;
  if ((y=fabs( *x))>=1.0) { return 0.0; }
  y2 = y * y;
  return 
    (1.0 + y2 * (-22.0 
		 + y2 * (165.0 
			 + y * (-231.0 
				+ y2 * (132.0 
					+ y2 * (-55.0 
						+ 10.0 * y2))))));
}
Real Dpenta(Real *x, Real *p, int effectivedim)
{ ///
  Real y,y2;
  if ((y=fabs( *x))>=1.0) { return 0.0; }
  y2 = y * y;
  return  
    y * (-14.66666666666666667 + 
	 y2 * (132.0 + 
	       y * (-192.5 + 
		    y2 * (115.5 + 
			  y2 * (-49.5 + 
				y2 * 9.16666666666666667)))));
  
 }
Real Scalepenta(Real *p,int scaling) {return 1.6552838957365;}
SimulationType methodpenta(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangepenta(int dim, int *index, Real* range){
 if(dim<=3) *index=-1; else *index=-2;
}



/* Tilmann Gneiting's space time models */

/* following function is unused! -- it has been nsst3*/
Real spacetime3(Real *x,Real *p, int effectivedim){
  // turning bands modified spacetime1 covariance function
  assert(effectivedim==2);
  return spacetime1(x,p,effectivedim) + fabs(x[0]) * Dspacetime1(x,p,effectivedim);
}


Real InvSqrtPsi(Real x, Real a, Real b, int c) {
  Real y;
  y = pow(fabs(x), a);
  switch(c) {
  case 1 : 
    return pow(1.0 + y, - 0.5 * b);
  case 2 : 
    return sqrt( (y + 1.0) / (y/b + 1.0) );
  case 3 :
    return sqrt(-log(b) / log(y + 1.0 / b));
  default: assert(false);
  }
}


/* Tilmann Gneiting's space time models, part I */
Real spacetime1(Real *x, Real *p, int effectivedim){
  Real y, z, invsqrtpsi;
  assert(effectivedim==2);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPAIII], p[KAPPAIV], (int)p[KAPPAV]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPAII]) {
  case 1 : z = stable(&y, p, 1); break;
  case 2 : z = WhittleMatern(&y, p, 1); break;
  case 3 : z = Cauchy(&y, p, 1); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPAVI]) * z;
}
Real TBM2spacetime1(Real *x, Real *p, int effectivedim){
  Real y, z, invsqrtpsi;
  assert(effectivedim==2);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPAIII], p[KAPPAIV], (int)p[KAPPAV]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPAII]) {
  case 1 : assert(false); break;
  case 2 : z = TBM2WhittleMatern(&y, p, 1); break;
  case 3 : z = TBM2Cauchy(&y, p, 1); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPAVI]) * z;
}
Real TBM3spacetime1(Real *x, Real *p, int effectivedim){
  assert(false); // should never happen
  return spacetime1(x, p, effectivedim) +
    Dspacetime1(x, p, effectivedim) * fabs(x[0]);
}
Real Dspacetime1(Real *x, Real *p, int effectivedim){
  Real y, z, invsqrtpsi;
  assert(effectivedim==2);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPAIII], p[KAPPAIV], (int)p[KAPPAV]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPAII]) {
  case 1 : z = Dstable(&y, p, 1); break;
  case 2 : z = DWhittleMatern(&y, p, 1); break;
  case 3 : z = DCauchy(&y, p, 1); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPAVI] + 1.0) * z;
}
int checkspacetime1(Real *param, int timespacedim, SimulationType method) {
  int error;
  // 1 : stable
  // 2 : whittle
  // 3 : cauchy (not generalised)
  // first parameter(s) for phi
  // then choice of phi; then two parameters for psi, then choice of psi
  if (timespacedim<=1) {
    strcpy(ERRORSTRING_OK,"total dim>=2");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  if (param[KAPPAII] != (double)((int) param[KAPPAII])) {
    strcpy(ERRORSTRING_OK,"kappa3 an integer");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPAII]);
    return ERRORCOVFAILED;
  }
  switch((int) param[KAPPAII]) {
  case 1 : 
    if (method==TBM2) {
      strcpy(ERRORSTRING_OK,"TBM2 and kappa2=2,3");
      sprintf(ERRORSTRING_WRONG,"%f", param[KAPPAII]);
      return ERRORCOVFAILED;
    }
    error=checkstable(param, timespacedim, method);
    break;
  case 2 : error=checkWhittleMatern(param, timespacedim, method); break;
  case 3 : error=checkCauchy(param, timespacedim, method); break;
  default : 
    strcpy(ERRORSTRING_OK,"kappa2=1,2,3");
    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPAIII]);
    return ERRORCOVFAILED;
  }
  if (error) return error;
  if ((param[KAPPAIII]<=0) || (param[KAPPAIII]>2)) {
    strcpy(ERRORSTRING_OK,"0<kappa3<=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAIII]);
    return ERRORCOVFAILED;
  }
  if ((param[KAPPAIV]<=0) || (param[KAPPAIV]>1)) {
    strcpy(ERRORSTRING_OK,"0<kappa4<=1");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAIV]);
    return ERRORCOVFAILED;
  } 
  if (param[KAPPAV] != (double)((int) param[KAPPAV])) {
    strcpy(ERRORSTRING_OK,"kappa5 an integer");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPAV]);
    return ERRORCOVFAILED;
  }
  switch((int) param[KAPPAV]) {
  case 1 : case 2: case 3: break;
  default :
    strcpy(ERRORSTRING_OK,"kappa5=1,2,3");
    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPAV]);
    return ERRORCOVFAILED;
  }
  if (timespacedim-1 > param[KAPPAVI]) {
    strcpy(ERRORSTRING_OK,"kappa6>=dim-1");
    sprintf(ERRORSTRING_WRONG,"%f for spatial dim=%d",
	    param[KAPPAVI],timespacedim-1);
    return ERRORCOVFAILED;
  }
  return 0;
}
SimulationType methodspacetime1(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangespacetime1(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract }
  Real *r;
  if ((*index<=0) || (*index>9)) { // see also last line
    int i; for (i=0; i<24; i++) range[i]= RF_NAN; *index=-1; return;
  }
  switch ((*index-1) % 3) {
  case 0 : r=range_stable;  break;
  case 1 : r=range_whittle; break;
  case 2 : r=range_cauchy;  break;
  default: assert(false);
  }
  memcpy(range, r, sizeof(Real) * 4);
    
  range[4] = range[6] =
    range[5] = range[7] = (Real) (1 + (*index - 1) % 3);

  range[8] = 0;
  range[9] = range[11] = 2;
  range[10] = RANGE_EPSILON;

  range[12] = 0;
  range[13] = 1;
  range[14] = RANGE_EPSILON;
  range[15] = 1.0 - RANGE_EPSILON;
  
  range[16] = range[17] = 
    range[18] = range[19] =  (Real) (1 + (*index - 1) / 3);
  
  range[20] = range[22] = (Real) dim;
  range[21] = RF_INF;
  range[23] = (Real) dim + 10.0;

  if ( (++(*index)) > 9) *index=-1;
}


/* Tilmann Gneiting's space time models, part II*/
Real spacetime2(Real *x,Real *p, int effectivedim){
  Real y, z, invsqrtpsi;
  assert(effectivedim==2);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPAIV], p[KAPPAV], (int)p[KAPPAVI]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPAIII]) {
  case 1 : z = generalisedCauchy(&y, p, 1); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPAVII]) * z;
}
Real TBM3spacetime2(Real *x, Real *p, int effectivedim){
  assert(false); // should never happen
  return spacetime2(x, p, effectivedim) +
    Dspacetime2(x, p, effectivedim) * fabs(x[0]);
}
Real Dspacetime2(Real *x, Real *p, int effectivedim){
  Real y, z, invsqrtpsi;
  assert(effectivedim==2);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPAIV], p[KAPPAV], (int)p[KAPPAVI]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPAIII]) {
  case 1 : z = DgeneralisedCauchy(&y, p, 1); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPAV] + 1.0) * z;
}
int checkspacetime2(Real *param, int timespacedim, SimulationType method) {
  int error;
  // 1 : generalisedcauchy 
  // first parameter(s) for phi
  // then choice of phi; then two parameters for psi, then choice of psi
  if (timespacedim<=1) {
    strcpy(ERRORSTRING_OK,"total dim>=2");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  if (param[KAPPAIII] != (double)((int) param[KAPPAIII])) {
    strcpy(ERRORSTRING_OK,"kappa3 an integer");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPAIII]);
    return ERRORCOVFAILED;
  }
  switch((int) param[KAPPAIII]) {
  case 1 : error=checkgeneralisedCauchy(param, timespacedim, method); break;
  default : 
    strcpy(ERRORSTRING_OK,"kappa3=1,2,3");
    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPAIII]);
    return ERRORCOVFAILED;
  }
  if (error) return error;

  if ((param[KAPPAIV]<=0) || (param[KAPPAIV]>2)) {
    strcpy(ERRORSTRING_OK,"0<kappa4<=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAIV]);
    return ERRORCOVFAILED;
  }
  if ((param[KAPPAV]<=0) || (param[KAPPAV]>1)) {
    strcpy(ERRORSTRING_OK,"0<kappa5<=1");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAV]);
    return ERRORCOVFAILED;
  } 
  if (param[KAPPAVI] != (double)((int) param[KAPPAVI])) {
    strcpy(ERRORSTRING_OK,"kappa6 an integer");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPAVI]);
    return ERRORCOVFAILED;
  }
  switch((int) param[KAPPAVI]) {
  case 1 : case 2: case 3: break;
  default :
    strcpy(ERRORSTRING_OK,"kappa6=1,2,3");
    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPAVI]);
    return ERRORCOVFAILED;
  }
  if (timespacedim>param[KAPPAVII]) {
    strcpy(ERRORSTRING_OK,"kappa7>=dim");
    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPAVII]);
  }
  return 0;
}
SimulationType methodspacetime2(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : return TBM3;
    case 2 : return TBM3;
    case 3 : return TBM3;
    default: return Direct;
    }
  }
}
void rangespacetime2(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract }
  Real *r;
  if ((*index<=0) || (*index>3)) { // see also last line
    int i; for (i=0; i<28; i++) range[i]= RF_NAN; *index=-1; return;
  }
  switch ((*index-1) % 1) {
  case 0 : r = range_genCauchy;  break;
  default: assert(false);
  }
  memcpy(range, r, sizeof(Real) * 8);

  range[8] = range[10] =
    range[9] = range[11] = (Real) (1 + (*index - 1) % 1);

  range[12] = 0;
  range[13] = range[15] = 2;
  range[14] = RANGE_EPSILON;

  range[16] = 0;
  range[17] = 1;
  range[18] = RANGE_EPSILON;
  range[19] = 1.0 - RANGE_EPSILON;

  range[20] = range[21] = 
    range[22] = range[23] =  (Real) (1 + (*index - 1) / 1);

  range[24] = range[26] = (Real) dim;
  range[25] = RF_INF;
  range[27] = (Real) dim + 10.0;
  
  if ( (++(*index)) > 3) *index= -1;
}


/* locally defined functions */
// Brownian motion 
Real fractalBrownian(Real *x, Real *p, int effectivdim) {
  return  pow(*x,p[KAPPA]);
}

Real twodimfractalBrownianlocal(Real *x, Real *p, int effectivdim) {
  register Real y;
  if (p[KAPPA]<=1.5) {  
    if ((y=fabs(*x))>1.0) return 0.0;
    return 1.0-pow(y,p[KAPPA])+0.5*p[KAPPA]*(y*y-1.0);
  } else {
    if ((y=fabs(*x))>2.0) return 0.0;
    if (y>1.0) {
      register Real z;
      z = 2.0-y;
      return p[KAPPA] * (2.0-p[KAPPA]) / 18.0 * z * z * z / y;
    } 
    return  
      1.0 - 0.16666666666666666666666666667 * p[KAPPA] * (1+p[KAPPA]) 
       - pow(y,p[KAPPA])  + p[KAPPA] * (5.0 + 2.0 * p[KAPPA])/ 18.0 * y *y;
  }
}

// also corrects the variance and other parameters in param !
Real twodimfractalBrownianS(Real *p) { 
  if (p[KAPPA]<=1.5) return 0.5*p[KAPPA];
  else return p[KAPPA]*(5+2.0*p[KAPPA])/18.0;
}
int check2dfractalBrownian(Real *param, int timespacedim, SimulationType method){
  if ((timespacedim>2) && (method!=Nothing)) {
    strcpy(ERRORSTRING_OK,"total dim<=2");
    sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
    return ERRORCOVFAILED;
  }
 if ((param[KAPPA]<=0.0) || (param[KAPPA]>=2.0)) {
    strcpy(ERRORSTRING_OK,"kappa in (0,2)");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
 if ((method!=Nothing) && (method!=CircEmbedLocal)) {
    strcpy(ERRORSTRING_OK,"method=local CE");
    sprintf(ERRORSTRING_WRONG,"%s",METHODNAMES[method]);
    return ERRORCOVFAILED;
  }
  return 0;  
}
SimulationType method2dfractalBrownian(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbedLocal;
    case 2 : return CircEmbedLocal;
    default: return Forbidden;
    }
  } else {
    switch(spacedim) {
    case 1 : return Direct;
    case 2 : return Direct;
    default: return Forbidden;
    }
  }
}
void range2dfractalBrownian(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  if (dim<=2) *index = -1; else *index=-2;
  Real r[4] = {0, 2, RANGE_EPSILON, 2};
  memcpy(range, r, sizeof(Real) * 4);
}


Real threedimfractalBrownianlocal(Real *x, Real *p, int effectivedim) {
  register Real y;
  if (p[KAPPA] <= 1.0) {
    if ((y=fabs(*x)) > 1.0) return 0.0;
    return 1.0 - pow(y, p[KAPPA]) + 0.5 * p[KAPPA] * (y * y - 1.0);
  } else {
    if ((y=fabs(*x)) > 2.0) return 0.0;
    if (y > 1.0) {
      register Real z;
      z = 2.0 - y;
      return p[KAPPA] * (2.0 - p[KAPPA]) / 18.0 * z * z * z / y;
    } 
    return  
      1.0 - 0.16666666666666666666666666667 * p[KAPPA] * (1+p[KAPPA]) 
       - pow(y,p[KAPPA])  + p[KAPPA] * (5.0 + 2.0 * p[KAPPA])/ 18.0 * y *y;
  }
}
Real threedimfractalBrownianS(Real *p) { 
  if (p[KAPPA] <= 1.0) return 0.5 * p[KAPPA]; 
  else return p[KAPPA] * (5 + 2.0 * p[KAPPA]) / 18.0;
}
int check3dfractalBrownian(Real *param, int timespacedim, SimulationType method){
  if ((timespacedim>3) && (method!=Nothing)) {
    strcpy(ERRORSTRING_OK,"total dim<=3");
    sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
    return ERRORCOVFAILED;
  }
  if ((param[KAPPA]<=0.0) || (param[KAPPA]>=2.0)) {
    strcpy(ERRORSTRING_OK,"kappa in (0,2)");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
  if ((method!=Nothing) && (method!=CircEmbedLocal)) {
    strcpy(ERRORSTRING_OK,"method=local CE");
    sprintf(ERRORSTRING_WRONG,"%s",METHODNAMES[method]);
    return ERRORCOVFAILED;
  }
  return 0;  
}
SimulationType method3dfractalBrownian(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbedLocal;
    case 2 : return CircEmbedLocal;
    case 3 : return CircEmbedLocal;
    default: return Forbidden;
    }
  } else {
    switch(spacedim) {
    case 1 : return Direct;
    case 2 : return Direct;
    case 3 : return Direct;
    default: return Forbidden;
    }
  }
}
void range3dfractalBrownian(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  if (dim<=3) *index = -1; else *index=-2;
  Real r[4] = {0, 2, RANGE_EPSILON, 2};
  memcpy(range, r, sizeof(Real) * 4);
}

Real fractGauss(Real *x, Real *p, int effectivedim){
  register Real y;
  if ((y = fabs(*x))==0.0) return 1.0;
  return 0.5 * (pow(fabs(y + 1.0), p[KAPPA])  
		- 2.0 * pow(y, p[KAPPA]) 
		+ pow(fabs(y - 1.0), p[KAPPA]) 
		);
}
int checkfractGauss(Real *param, int timespacedim, SimulationType method) {
  if ((timespacedim!=1) && (method!=Nothing)) {
    strcpy(ERRORSTRING_OK,"dim=1");
    sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
    return ERRORCOVFAILED;
  }
  if ((param[KAPPA]<=0) || (param[KAPPA]>2.0)) {
    strcpy(ERRORSTRING_OK,"0<kappa<=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return 0;
}
SimulationType methodfractGauss(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    default: return Forbidden;
    }
  } else {
    switch(spacedim) {
    case 1 : return Direct;
    default: return Forbidden;
    }
  }
}
void rangefractGauss(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  if (dim==1) *index = -1; else *index=-2;
  Real r[4] = {0, 2, RANGE_EPSILON, 2};
  memcpy(range, r, sizeof(Real) * 4);
}




/* local-global distinguisher */
Real lgd1(Real *x, Real*p, int effectivedim) {
  Real y;
  if ((y=fabs(*x)) < 1) return 1.0 - p[KAPPAII] / (p[KAPPAI] + p[KAPPAII]) 
	     * exp(p[KAPPAI] * log(y));
  else return p[KAPPAI] / (p[KAPPAI] + p[KAPPAII]) 
	 * exp( -p[KAPPAII] * log(y));
}
Real Scalelgd1(Real *p,int scaling) {
  if (19 * p[KAPPAI] < p[KAPPAII]) 
    return exp( log(0.95 * (p[KAPPAI] + p[KAPPAII]) / p[KAPPAII]) / p[KAPPAI]); 
  else return exp(log(0.05 * (p[KAPPAI] + p[KAPPAII]) / p[KAPPAI])/p[KAPPAII]);
}
Real Dlgd1(Real *x, Real *p, int dim){
  Real y, pp;
  if ( (y=fabs(*x)) == 0) return 0; // falscher Wert, aber sonst gibt NAN-Fehler
  pp = ( (y < 1) ? p[KAPPAI] : -p[KAPPAII] ) - 1.0;
  return - p[KAPPAI] * p[KAPPAII] / (p[KAPPAI] + p[KAPPAII]) * exp(pp * y);
}
/*
Real TBM3(Real *x, Real*p, int effectivedim) {
  register Real y;
}
Real D(Real *x, Real*p, int effectivedim) {
  register Real y;
}
*/
SimulationType methodlgd1(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return CircEmbed;
    case 2 : return CircEmbed;
    default: return CircEmbed;
    }
  } else {
    switch(spacedim) {
    case 1 : 
    case 2 :
    case 3 : 
    default: return Direct;
    }
  }
}
int checklgd1(Real *param, int timespacedim, SimulationType method) {
  if ((timespacedim>2) && (method!=Nothing)) {
    strcpy(ERRORSTRING_OK, "dim<=2");
    sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
    return ERRORCOVFAILED;
  }
  if (((param[KAPPAI]<=0) || (param[KAPPAI]>1.0)) && (timespacedim==1)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=1, dim=1");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAI]);
    return ERRORCOVFAILED;
  }
  if (((param[KAPPAI]<=0) || (param[KAPPAI]>0.5)) && (timespacedim==2)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=1/2, dim=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAI]);
    return ERRORCOVFAILED;
  }
  if (param[KAPPAII]<=0) {
    strcpy(ERRORSTRING_OK,"kappa2>0");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAII]);
    return ERRORCOVFAILED;
  }
  return 0;
}

static Real range_lgd1[8] = {0.0, 1.0, 0.01, 1.0, 
			     0.0, RF_INF, 0.01, 20.0};
void rangelgd1(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  if (dim<=2) *index = -1; else *index=-2;
  memcpy(range, range_lgd1, sizeof(Real) * 8);
  if (dim==2) range[1]=range[3]=0.5;
}


/* FD model */
Real FD(Real *x,Real *p, int effectivedim){
  static Real dold=RF_INF;
  static Real kold, sk;
  Real y, d, k, skP1;
  d = p[KAPPAI];
  y = fabs(*x);
  k = (Real) trunc((double) y);
  if (dold!=d || kold > k) {
    sk = 1;
    kold = 0.0;
  }
  // sign (-1)^k is (kold+d), 16.11.03, checked. 
  for (; kold<k; kold+=1.0) sk =  sk * (kold + d) / (kold + 1.0 - d);
  dold = d;
  kold = k;
  if (k == y) return sk; 
  skP1 = sk * (kold + d) / (kold + 1.0 - d);
  return sk + (y - k) * (skP1 - sk);
}
SimulationType methodFD(int spacedim, bool grid){
 return (grid ? CircEmbed : Direct);
}
int checkFD(Real *param, int timespacedim, SimulationType method) {
  if (param[KAPPAI] < -0.5 || param[KAPPAI] >= 0.5) {
    strcpy(ERRORSTRING_OK,"-0.5 <= kappa < 0.5");
    sprintf(ERRORSTRING_WRONG, "%f",param[KAPPAI]);
    return ERRORCOVFAILED;    
  } 
  if (timespacedim>1) {
    strcpy(ERRORSTRING_OK, "dim=1");
    sprintf(ERRORSTRING_WRONG, "%d", timespacedim);
    return ERRORCOVFAILED;
  } 
  return 0;
}
void rangeFD(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  if (dim==1) *index = -1; else *index=-2;
  range[2] = (range[0] = -0.5) + RANGE_EPSILON;
  range[4] = (range[1] = 0.5) - RANGE_EPSILON;
}


// ---------------------------------------------------------------------



/*
SimulationType method(int spacedim, bool grid){
  if (grid) {
    switch(spacedim) {
    case 1 : return ;
    case 2 : return ;
    default: return ;
    }
  } else {
    switch(spacedim) {
    case 1 : return ;
    case 2 : return ;
    case 3 : return ;
    default: return Direct;
    }
  }
}


void range(int dim, int *index, Real* range){
  //  2 x length(param) x {theor, pract } 
  Real r[] = {};
  memcpy(range, r, sizeof(Real) * );
}


*/
