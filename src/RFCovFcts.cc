/*
 Authors 
 Martin Schlather, Martin.Schlather@uni-bayreuth.de 

 Definition of covariance functions and derivatives (spectral measures, 
 tbm operators)
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * when defining your own function, make sure that the covariance function 
   itself allows for an additional nugget effect (spectral measures and tbm 
   operators don't)     

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



#define MINUSINVLOG005 0.3338082006953340674649
#define SQRTINVLOG005 0.5777613700268771079749
#include <math.h>
#include <Rmath.h>
#include <assert.h>
#include "RFsimu.h"
#include "RFCovFcts.h"


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

//// NOTE `*p' may not be changed by any of the functions!

/* nugget effect model */
Real nugget(Real x,Real *p){
  if (x==0.0) {return p[SILL];} 
  return 0.0;
}
Real Scalenugget(Real *p, int scaling) { return 0.0; } //error!
int checknugget(key_type *key){
  if ((key->param[VARIANCE]!=0.0) &&
      (key->param[NUGGET]!=0.0)) {
    strcpy(ERRORSTRING_OK,
	   "either the variance or the nugget being 0");
    sprintf(ERRORSTRING_WRONG,"variance=%f nugget=%f",
	    key->param[VARIANCE], key->param[NUGGET]);
    return ERRORCOVFAILED;
  }
  key->param[NUGGET] += key->param[VARIANCE];
  key->param[VARIANCE] = 0.0;
  return 0;
}

/* exponential model */
Real exponential(Real x,Real *p){
  if (x==0.0) {return(p[SILL]);}
  return p[VARIANCE]*exp(-fabs(x*p[INVSCALE]));
}
Real Scaleexponential(Real *p,int scaling){ return MINUSINVLOG005; } 
Real TBM3exponential(Real x, Real *p){
   register Real y;
   y = fabs(x*p[INVSCALE]);
   return p[VARIANCE]*(1.0-y)*exp(-y);
}

/* derivative of exponential: qexponential */
Real qexponential(Real x,Real *p){
  register Real y;
  if (x==0.0) {return(p[SILL]);}
  y = exp(-fabs(x*p[INVSCALE]));
  return p[VARIANCE] * ( y * (2.0  - p[KAPPA] * y) / (2.0 - p[KAPPA]) );
}
Real Scaleqexponential(Real *p,int scaling){
  return -1.0 / 
    log( (1.0 - sqrt(1.0 - p[KAPPA] * (2.0 - p[KAPPA]) * 0.05)) / p[KAPPA]);
} 
int checkqexponential(key_type *key){
  if ((key->param[KAPPA]<0.0) || (key->param[KAPPA]>1.0)) {
    strcpy(ERRORSTRING_OK,"kappa in [0,1]");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return 0;  
}

/* derivative of exponential: hole effect */
Real holeeffect(Real x, Real*p){
  register double y;
  if (x==0.0) {return(p[SILL]);}
  y = fabs(x*p[INVSCALE]);
  return p[VARIANCE]*exp(-y * p[KAPPA]) * cos(y);
}
Real Scaleholeeffect(Real *p,int scaling){ 
  if (scaling==NATSCALE_EXACT) return 0.0;
  return MINUSINVLOG005; 
} 
Real TBM3holeeffect(Real x, Real *p){
  register double y;
  y = fabs(x*p[INVSCALE]);
  return p[VARIANCE] * exp(-p[KAPPA]*y) *
    ((1.0 - p[KAPPA] * y) * cos(y) - y * sin(y));
}
int checkholeeffect(key_type *key){
  if ((key->dim==3) || (key->method==TBM3)){
    if (key->param[KAPPA]<1.73205080756889) {
      strcpy(ERRORSTRING_OK,
	     "kappa >= sqrt(3) for 3-dimensional simulations and techniques");
      sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
      return ERRORCOVFAILED;
    }
  } else if (key->dim==2){
    if (key->param[KAPPA]<1.0) {
      strcpy(ERRORSTRING_OK,"kappa >= 1.0 for 2 dimensions");
      sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
      return ERRORCOVFAILED;
    }
  }
  return 0;  
}
 

/* Class of "spherical" models  */ 
/* circular model */
Real circular(Real x, Real *p)
{
  Real y;
  if ((y=fabs(x*p[INVSCALE]))>1.0) return 0.0; 
  if (x==0) {return p[SILL];} 
  return p[VARIANCE]* ( 1.0 - (2.0 * (y * sqrt(1.0- y * y) + asin(y))) / PI);
}
Real Scalecircular(Real *p,int scaling) {return 1.138509531721630274603;}
 
/* spherical model */ 
Real spherical(Real x, Real *p){
  register Real y;
  if (x==0.0) {return p[SILL];} 
  if ((y= fabs(x*p[INVSCALE]))>1.0) {return 0.0;}
  return p[VARIANCE]*(1.0+y*0.5*(y*y-3.0));
}
Real Scalespherical(Real *p,int scaling){ return 1.23243208931941;}
Real TBM2spherical(Real x, Real *p){
  register Real y, y2;
  y=fabs(x*p[INVSCALE]); y2=y*y;
  if (y>1.0)
    {return p[VARIANCE]*(1.0- 0.75 * y * ((2.0 - y2) * asin(1.0/y) + 
					  sqrt(y2 -1.0)));} 
  return p[VARIANCE]*(1.0 - 0.375 * PI * y * (2.0 - y2));
}
Real TBM3spherical(Real x, Real *p){
  register Real y;
  if ((y=fabs(x*p[INVSCALE]))>1.0) {return 0.0;} 
  return p[VARIANCE]*(1.0 + (-3.0 + 2.0 * y * y) * y);
}


/* power model */ 
Real power(Real x, Real *p){
  register Real y;
  if (x==0.0) {return p[SILL];} 
  if ((y=fabs(x*p[INVSCALE]))>1.0) {return 0.0;}
  return p[VARIANCE] * pow(1.0 - y, p[KAPPA]);
}
Real Scalepower(Real *p,int scaling){ 
  return 1.0 / (1.0 - pow(0.05,1/p[KAPPA]));
}
Real TBM2power(Real x, Real *p){
  // only kappa=2 up to now !
  register Real y;
  y=fabs(x*p[INVSCALE]);
  if (y>1.0)
    {return p[VARIANCE]*(1.0 - 2.0 * y *(asin(1.0/y) - y + sqrt(y*y-1.0) ));} 
  return p[VARIANCE]*(1.0 - y * (PI - 2.0 * y));
}
Real TBM3power(Real x, Real *p){
  register Real y;
  if ((y=fabs(x*p[INVSCALE]))>1.0) {return 0.0;} 
  return p[VARIANCE] * (1.0 - y - y * p[KAPPA]) *
    pow(1.0 - y, p[KAPPA]-1.0);
}
int checkpower(key_type *key) {
  switch (key->method) {
  case TBM2 : 
    if (key->param[KAPPA]!=2) { 
      strcpy(ERRORSTRING_OK,"kappa=2 in case of TBM2");
      sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
      return ERRORCOVFAILED;
    }
    break;
  case TBM3 : 
    if (key->param[KAPPA]<2) { 
      strcpy(ERRORSTRING_OK,"kappa>=2");
      sprintf(ERRORSTRING_WRONG,"method=TBM3 and kappa=%f",key->param[KAPPA]);
      return ERRORCOVFAILED;
      }
    break;
  default: 
    if (key->param[KAPPA] <  0.5 * (1.0 + key->dim)) {
      strcpy(ERRORSTRING_OK,"kappa >= (dim+1)/2");
      sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
      return ERRORCOVFAILED;
    }
  }
  return 0;
}


/* stable model */
Real stable(Real x,Real *p){
  if (x==0.0) {return(p[SILL]);}
  return  p[VARIANCE]*(exp(-pow(fabs(x*p[INVSCALE]),p[KAPPA])));
}
Real Scalestable(Real *p,int scaling){return pow(MINUSINVLOG005,1/p[KAPPA]);}
Real TBM3stable(Real x, Real *p){
  register double y;
  y = pow(fabs(x*p[INVSCALE]),p[KAPPA]);
  return p[VARIANCE]*exp(-y) * (1 - p[KAPPA] * y);
}
int checkstable(key_type *key) {
  if ((key->param[KAPPA]>0) && (key->param[KAPPA]<=2.0)) return 0;
  strcpy(ERRORSTRING_OK,"0<kappa<=2");
  sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
  return ERRORCOVFAILED;
}

/* Whittle-Matern or Whittle or Besset */ 
Real WhittleMatern(Real x, Real *p) 
// check calling functions, like hyperbolic and gneiting if any changings !!
{
  static double kappa=RF_INF;
  static double loggamma;
  register double y;
  if (x==0.0) {return p[SILL];}
  y = fabs(x*p[INVSCALE]); 
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
#ifdef RF_GSL
    gsl_sf_result result; 
    gsl_sf_lngamma_impl(kappa,&result);
    loggamma = result.val;
#else
    loggamma = lgammafn(kappa);
#endif
  }
  return p[VARIANCE]*2.0 *
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
Real TBM3WhittleMatern(Real x, Real *p)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double kappa=RF_INF;
  static double loggamma;
  register double y,loghalfy;
  if (x==0.0) {return p[VARIANCE];}
  y = fabs(x*p[INVSCALE]); 
  loghalfy = log(0.5*y);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
#ifdef RF_GSL
    gsl_sf_result result; 
    gsl_sf_lngamma_impl(kappa,&result);
    loggamma = result.val;
#else
    loggamma = lgammafn(kappa);
#endif
  }
  return p[VARIANCE]* 
    (2.0 *exp(kappa*loghalfy -loggamma + log(bessel_k(y,kappa,2.0))-y)
     -(4.0 *exp((kappa+1.0)*loghalfy-loggamma + log(bessel_k(y,kappa-1.0,2.0))-y))
     );
}
Real spectralWhittleMatern(Real *p END_WITH_GSLRNG) { /* see Yaglom ! */
  return  p[INVSCALE]*sqrt(pow(1.0 - UNIFORM_RANDOM,-1.0/p[KAPPA]) -1);
}
int checkWhittleMatern(key_type *key) { 
  static Real spectrallimit=0.17;
  switch(key->method) {
  case SpectralTBM : 
    if ((key->param[KAPPA]>=spectrallimit)) return 0;
    sprintf(ERRORSTRING_OK,
	    "%f<=kappa (the numerical errors are too big for 0<kappa<%f)",
	    spectrallimit,spectrallimit);
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  default :
    if ((key->param[KAPPA]>0)) return 0;
    strcpy(ERRORSTRING_OK,"0<kappa");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  }
}

Real hyperbolic(Real x, Real*p){ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double y,kappay;
  if (x==0.0) {return p[SILL];}
  if (p[KAPPAXX]==0) { // whittle matern
    return WhittleMatern(x * p[KAPPA],p);
  } 
  if (p[KAPPA]==0) { //cauchy   => KAPPAX < 0 !
    y = x*p[INVSCALE]/p[KAPPAXX];
    /* note change in sign as KAPPAX<0 */
    return p[VARIANCE]* pow(1+y*y,p[KAPPAX]); 
  }  
  y=fabs(x*p[INVSCALE]); 
  if ((p[KAPPA]!=kappa) || (p[KAPPAX]!=lambda) || (p[KAPPAXX]!=delta)) {
    kappa=p[KAPPA]; 
    lambda= p[KAPPAX];
    delta = p[KAPPAXX];
    deltasq = delta * delta;
    kappadelta = kappa * delta;
    logconst = kappadelta - log(bessel_k(kappadelta,lambda,2.0)) 
      - lambda * log(delta);
  }
  y=sqrt(deltasq + y*y);  
  kappay = kappa * y;
  return p[VARIANCE]* 
    exp(logconst + lambda * log(y)+log(bessel_k(kappay,lambda,2.0))-kappay);
}
Real TBM3hyperbolic(Real x, Real*p)
{ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double lambdaDkappa;
  static double kappadelta;
  static double logconst;
  double y,ysq,s,kappas,logs;
  if (x==0.0) {return p[VARIANCE];}
  if (p[KAPPAXX]==0) { // whittle matern
    return TBM3WhittleMatern(x * p[KAPPA],p);
  } 
  if (p[KAPPA]==0) { //cauchy
    register Real y,ha;
    y=x*p[INVSCALE]/p[KAPPAXX];
    ha=y*y;
    /* note change in sign as KAPPAX<0 */
    return p[VARIANCE]*(1.0+ (1.0+2.0*p[KAPPAX])*ha) * 
      pow(1.0+ha,p[KAPPAX]-1.0);
  }
  y=fabs(x*p[INVSCALE]); 
  if ((p[KAPPA]!=kappa) || (p[KAPPAX]!=lambda) || (p[KAPPAXX]!=delta)) {
    kappa=p[KAPPA]; 
    lambda= p[KAPPAX];
    delta = p[KAPPAXX];
    deltasq = delta * delta;
    kappadelta = kappa * delta;
    lambdaDkappa= lambda/kappa;
    logconst = kappadelta - log(bessel_k(kappadelta,lambda,2.0)) 
      - lambda * log(delta);
  }
  ysq = y * y;
  s=sqrt(deltasq + ysq);
  kappas = kappa * s;
  logs = log(s);  
  return p[VARIANCE] * 
    ( exp(logconst + lambda * logs +log(bessel_k(kappas,lambda,2.0))-kappas)
      - ysq*kappa*exp(logconst + (lambda-1.0)*logs 
		      +log(bessel_k(kappas,lambda-1.0,2.0))-kappas)
      );
}

int checkhyperbolic(key_type *key){
  if (key->param[KAPPAX]>0) {
    if ((key->param[KAPPAXX]<0) || (key->param[KAPPA]<=0)) {
      strcpy(ERRORSTRING_OK,"kappa1>0 and kappa3>=0 if kappa2>0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      key->param[KAPPA],key->param[KAPPAXX],key->param[KAPPAX]);
      return ERRORCOVFAILED;
    }
  } else if (key->param[KAPPAX]<0) {
    if ((key->param[KAPPAXX]<=0) || (key->param[KAPPA]<0)) {
      strcpy(ERRORSTRING_OK,"kappa1>=0 and kappa3>0 if kappa2<0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      key->param[KAPPA],key->param[KAPPAXX],key->param[KAPPAX]);
      return ERRORCOVFAILED;
    }
  } else { // key->param[KAPPAX]==0.0
    if ((key->param[KAPPAXX]<=0) || (key->param[KAPPA]<=0)) {
      strcpy(ERRORSTRING_OK,"kappa1>0 and kappa3>0 if kappa2=0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      key->param[KAPPA],key->param[KAPPAXX],key->param[KAPPAX]);
      return ERRORCOVFAILED;
    }
  }
  return 0;
}


/* Gneiting's functions */
// alternative to Gaussian */
#define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
Real Gneiting(Real x, Real *p){ 
  register Real y,oneMy8;
  if (x==0.0) {return(p[SILL]);} 
  if ((y=fabs(x*p[INVSCALE]*Sqrt2TenD47))>1.0) {return(0.0);}
  oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
  return p[VARIANCE]*((1.0+y * ( 8.0 + y * (25.0 + 32.0 *y)))*oneMy8);
}
Real ScaleGneiting(Real *p,int scaling) {return 0.5854160193;}
Real TBM3Gneiting(Real x, Real *p){ 
  register Real y,oneMy7;
  if ((y=fabs(x*p[INVSCALE]*Sqrt2TenD47))>1.0) {return 0.0;}  
  oneMy7 = 1.0-y; oneMy7*=oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);
  return p[VARIANCE]*
    (1.0 + y * (7.0  -  y * (5.0 + y * (147.0 + 384.0 * y))))* oneMy7;
}

// differentiable and compact support 
Real Gneitingdiff(Real x, Real *p){ 
  /* three Parameters!  
     p[SCALE]: overall scale 
     p[KAPPA]: whittle-matern shape
     p[KAPPAX]: gneiting scale
   */
  Real y;
  register Real oneMy8;
  if (x==0) {return(p[SILL]);} 
  if ((y=fabs(x*p[INVSCALE]/p[KAPPAX]))>1.0) {return(0.0);} 
  oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
  return (1.0 + y * (8.0 + y * (25.0 + y * 32.0)))*oneMy8*WhittleMatern(x,p); 
}
Real TBM3Gneitingdiff(Real x, Real *p)
{ 
  static double kappa=RF_INF;
  static double loggamma;
  double y,gneiting, dgneiting,loghalfy;
  if (x==0.0) {return p[VARIANCE];} 
  y = fabs(x*p[INVSCALE]);
  { 
    register double oneMz4;
    double z;
    if ((z=y/p[KAPPAX])>=1.0) {return(0.0);} 
    { 
      register double oneMz2;
      oneMz2 = 1.0-z; oneMz2 *= oneMz2;
      oneMz4 = oneMz2 * oneMz2;
      dgneiting = - z/p[KAPPAX] * (22.0 + z * (154.0 + z * 352.0)) * 
	(1.0 -z) * oneMz2 * oneMz4;
    }
    gneiting = (1.0 + z * ( 8.0 + z * (25.0 + 32.0 *z))) * oneMz4 * oneMz4;
  }
  loghalfy = log(0.5*y);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
#ifdef RF_GSL
    gsl_sf_result result; 
    gsl_sf_lngamma_impl(kappa,&result);
    loggamma = result.val;
#else
    loggamma = lgammafn(kappa);
#endif
  }
  return p[VARIANCE]* 
    (2.0 *exp(kappa*loghalfy - loggamma+ log(bessel_k(y,kappa,2.0))-y)
     * (gneiting + y * dgneiting)
     -4.0 *exp((kappa+1.0)*loghalfy -loggamma + log(bessel_k(y,kappa-1.0,2.0))-y)
     * gneiting
     );
}
int checkGneitingdiff(key_type *key) {
  if (key->param[KAPPAX]<=0.0) {
    strcpy(ERRORSTRING_OK,"0<kappa2");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPAX]);
    return ERRORCOVFAILED;
  }
  return checkWhittleMatern(key); 
}

Real genGneiting(Real x, Real *p)
{
  register Real y, s;
 
  if (x==0.0) return p[SILL];
  if ((y=fabs(x*p[INVSCALE]))>1.0) return 0.0; 
  s = p[KAPPAX] + p[KAPPA];
  switch ((int) p[KAPPA]) {
  case 1:
    return p[VARIANCE] * pow ((1.0-y), s) * (1.0 + s*y);   
  case 2:
    return p[VARIANCE] * pow ((1.0-y), s) * 
      (1.0 + y * (s + y * (s*s-1.0)*0.3333333333333333));  
  case 3:
    return p[VARIANCE] *  pow ((1.0-y), s) * 
      (1.0 + y * (s 
		  + y * ( 0.2 * (2.0*s*s - 3.0)
			 + y * (s*s-4.0)*s*0.0666666666666666666))) ;
  default : assert(false);  
  }
}
Real TBM3genGneiting(Real x, Real *p)
{
  register Real y, s;
  if ((y=fabs(x*p[INVSCALE]))>1.0) {return 0.0;} 
  s = p[KAPPAX] + p[KAPPA];
  switch ((int) p[KAPPA]) {
  case 1:
    return p[VARIANCE] * pow((1.0 - y), p[KAPPAX]) *
      (1.0 + y * (p[KAPPAX] - y * s * (p[KAPPAX] + 3.0))); 
  case 2:
   return p[VARIANCE] * pow((1.0 - y), (p[KAPPAX]+1.0)) *
      (1.0 + y * ((p[KAPPAX]+1.0)
		- y *(1.0 + 2.0 * s 
		      + y * (s * s - 1.0) * 0.33333333333333 * 
		      (p[KAPPAX] + 5.0))));
  case 3:
    return p[VARIANCE] * pow((1.0 - y), (p[KAPPAX]+2.0)) *
      (1.0 + y * ((p[KAPPAX]+2.0)
		  + y * (0.2 * (-9.0 + s * (-10.0 + s))
			 + y * 0.0666666666666666 * 
			 ((27.0 - s*(7.0 + s*(18.0 + 2.0*s)))
			  - y * (s*s - 4.0) * s * (p[KAPPAX] + 7.0)))));
  default : assert(false);   
  }
}
int checkgenGneiting(key_type* key)
{
  if ((key->param[KAPPA] != (double) ((int)key->param[KAPPA])) || 
      (key->param[KAPPA]<0)) {
      strcpy(ERRORSTRING_OK,"positive integer kappa");
      sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
      return ERRORCOVFAILED;
  }
  if ((key->param[KAPPA] < 1.0) || (key->param[KAPPA] > 3.0))
    return ERRORNOTPROGRAMMED;
  if ((key->method==TBM3) && (key->param[KAPPAX] < 2.0 + key->param[KAPPA])) {
    strcpy(ERRORSTRING_OK,"kappa2 <= 2.0 + kappa1");
    sprintf(ERRORSTRING_WRONG, "method=TBM3, kappa1=%f and kappa2=%f ",
	    key->param[KAPPA], key->param[KAPPAX]);
    return ERRORCOVFAILED;    
  }
  if  (key->param[KAPPAX] < ( 0.5 * key->dim + key->param[KAPPA] + 0.5)) {
    strcpy(ERRORSTRING_OK, "kappa2 >= (dim + 2*kappa1 + 1)/2 ");
    sprintf(ERRORSTRING_WRONG, "kappa1=%f and kappa2=%f ",
	    key->param[KAPPA], key->param[KAPPAX]);
    return ERRORCOVFAILED;
  }
  return 0;
}



/* Gausian model */
Real Gauss(Real x, Real*p) {
  register Real y;
  if  (x==0.0) {return(p[SILL]);} 
  y = p[INVSCALE] * x;
  return  p[VARIANCE]*exp(- y *y);
}
Real ScaleGauss(Real *p,int scaling) {return SQRTINVLOG005;}
Real TBM3Gauss(Real x, Real*p) {
  register Real y;
  y =  p[INVSCALE] * x; 
  y *= y;
  return  p[VARIANCE]*(1-2.0*y)*exp(-y);
}
Real spectralGauss(Real *p END_WITH_GSLRNG) {   
  return 2.0 * p[INVSCALE] * sqrt(-log(1.0 - UNIFORM_RANDOM));
}


/* Cauchy models */
Real Cauchy(Real x, Real *p){
  register Real y;
  if  (x==0) {return(p[SILL]);} 
  y = x*p[INVSCALE];
  return p[VARIANCE]* pow(1+y*y,-p[KAPPA]);
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
Real TBM2Cauchy(Real x, Real *p){
  register Real y2, lpy2;
  y2 = x*p[INVSCALE]; y2 *= y2; 
  switch ((int) (p[KAPPA]*2.0 + 0.1)) {
  case 1 : return p[VARIANCE] / (1 + y2);
  case 3 : lpy2=1.0+y2; return p[VARIANCE]* (1-y2)/ (lpy2*lpy2);
  case 5 : lpy2=1.0+y2; 
    return p[VARIANCE]*
      (1.0-y2*(2.0+0.333333333333333333333*y2))/(lpy2*lpy2*lpy2);
  case 7 : lpy2=1.0+y2; lpy2*=lpy2;
    return p[VARIANCE]*(1.0- y2*(3.0+y2*(1.0+0.2*y2)))/lpy2;
  default : assert(false);
  }
}
Real TBM3Cauchy(Real x, Real *p){
  register Real y,ha;
  y=x*p[INVSCALE];
  ha=y*y;
  return p[VARIANCE]*(1.0+ (1.0-2.0*p[KAPPA])*ha) * pow(1.0+ha,-p[KAPPA]-1.0);
}
int checkCauchy(key_type *key){
  switch (key->method) {
  case TBM2 : 
    if ((key->param[KAPPA]==0.5) || (key->param[KAPPA]==1.5) ||
	(key->param[KAPPA]==2.5) || (key->param[KAPPA]==3.5)) {return 0;}
    strcpy(ERRORSTRING_OK,"kappa in {0.5, 1.5, 2.5 ,3.5} (TBM2)");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  default : 
    if (key->param[KAPPA]>0) return 0;
    strcpy(ERRORSTRING_OK,"0<kappa");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  }
}

// (1-h^param[2])^(-param[3]/param[2])
Real generalisedCauchy(Real x, Real *p){
  if (x==0) {return(p[SILL]);}
  return p[VARIANCE]*pow(1.0+pow(fabs(x*p[INVSCALE]),p[KAPPA]),
				   -p[KAPPAX]/p[KAPPA]);
}
Real ScalegeneralisedCauchy(Real *p,int scaling) {
  switch(scaling) {
  case NATSCALE_EXACT: case NATSCALE_APPROX:
    return pow(pow(0.05,-p[KAPPA]/p[KAPPAX])-1.0,-1.0/p[KAPPA]); 
    break;
  case NATSCALE_MLE: 
    // should be changed! (long tails!)
    return pow(pow(0.05,-p[KAPPA]/p[KAPPAX])-1.0,-1.0/p[KAPPA]);
    break;
  default: assert(false);
  }
}
Real TBM3generalisedCauchy(Real x, Real *p){
  register Real ha;
  ha=pow(fabs(x*p[INVSCALE]),p[KAPPA]);
  return p[VARIANCE]*
    (1.0+ (1.0-p[KAPPAX])*ha) * pow(1.0+ha,-p[KAPPAX]/p[KAPPA]-1.0);
}
int checkgeneralisedCauchy(key_type *key){
  if ((key->param[KAPPA]<=0) || (key->param[KAPPA]>2.0)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=2");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
    if (key->param[KAPPAX]>0) return 0;
  }
  if (key->param[KAPPAX]<=0) {
    strcpy(ERRORSTRING_OK,"0<kappa2");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPAX]);
    return ERRORCOVFAILED;
    }
  return 0;
}

Real Cauchytbm(Real x, Real *p){
  register Real ha;
  if (x==0) {return(p[SILL]);}
  ha=pow(fabs(x*p[INVSCALE]),p[KAPPA]);
  return p[VARIANCE]*
    (1.0+ (1.0-p[KAPPAX]/p[KAPPAXX])*ha) * pow(1.0+ha,-p[KAPPAX]/p[KAPPA]-1.0);
}
Real TBM3Cauchytbm(Real x, Real *p){
  register Real bg,ha;
  ha=pow(fabs(x*p[INVSCALE]),p[KAPPA]);
  bg=p[KAPPAX]/p[KAPPAXX];
  return p[VARIANCE] * 
    (1 + ha * (1-bg*(1+p[KAPPA])+(1-p[KAPPAX]) * (1+(1-bg)*ha)))*
    pow(1+ha,-p[KAPPAX]/p[KAPPA]-2.0);
}
int checkCauchytbm(key_type *key){
  if ((key->method==TBM3) && (3.0 > key->param[KAPPAXX])) {
    strcpy(ERRORSTRING_OK,"3 <= kappa3 and method=TBM3");
    sprintf(ERRORSTRING_WRONG,
	    "method=TBM3 and kappa3=%f ",key->param[KAPPAXX]);
    return ERRORCOVFAILED;    
  }
  if (key->dim > key->param[KAPPAXX]) {
    strcpy(ERRORSTRING_OK,"dim <= kappa3");
    sprintf(ERRORSTRING_WRONG,
	    "kappa3=%f < dim=%d ",key->param[KAPPAXX],key->dim);
    return ERRORCOVFAILED;    
  }
  // theory: should check whether p[KAPPAXX] is an integer. 
  // But should go through without this restriction
  if ((key->param[KAPPA]<=0) || (key->param[KAPPA]>2.0)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=2");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  }
  if (key->param[KAPPAX]<=0) {
    strcpy(ERRORSTRING_OK,"0<kappa2");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPAX]);
    return ERRORCOVFAILED;
  }
  return 0;
}


/* Bessel function */ 
Real Bessel(Real x,Real *p){
  static double kappa=RF_INF;
  static double gamma;
  Real y;
  if  (x==0.0) {return p[SILL];} 
  y = fabs(x*p[INVSCALE]);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
#ifdef RF_GSL
    gsl_sf_result result; 
    gsl_sf_gamma_impl(p[KAPPA]+1.0,&result); // _lngamma_
    gamma = result.val;
#else 
    gamma = gammafn(kappa+1.0);
#endif
  }
  return p[VARIANCE]* gamma  * pow(0.5 * y,-kappa) * bessel_j(y,kappa);
}
Real spectralBessel(Real *p END_WITH_GSLRNG) {
  return p[KAPPA]==0 ? p[INVSCALE] : 
    (p[INVSCALE] * sqrt(1-pow(UNIFORM_RANDOM,1.0/p[KAPPA])));
}
int checkBessel(key_type *key){
  // Whenever TBM3Bessel exists, add further check against too small kappa!
  if (key->param[KAPPA]<0.5*((double)key->dim)-1.0) {
    strcpy(ERRORSTRING_OK,"kappa >= dim/2 - 1");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return 0;
}


Real expPLUScirc(Real x, Real *p){ 
  Real y,z;
  if (x==0) {return(p[SILL]);}
  z = fabs(x * p[INVSCALE]);
  y = fabs(z/p[KAPPAX]);
  return p[VARIANCE]*
       (p[KAPPA] * exp(-z) + 
	((y>1.0) ? 0.0 : 
	 (1-p[KAPPA])*(1.0-(2.0 *(y * sqrt(1.0 - y * y)+asin(y)))/PI)));
}
Real ScaleexpPLUScirc(Real *p,int scaling) { // pretty heuristic!!!
  Real x,weight;
  if (scaling==NATSCALE_EXACT) return 0.0; //failure!
  x =  -log(0.05/p[KAPPA]);
  weight =  x /  p[KAPPAX]; 
  if (p[KAPPA]>weight) {weight=p[KAPPA];}  else if (weight>1.0) {weight=1.0;}
  return  1.0 /(weight * x +  (1.0-weight) * p[KAPPAX]);
}
int checkexpPLUScirc(key_type *key){
  if ((key->param[KAPPA]<0.0) || (key->param[KAPPA]>1.0)){
    strcpy(ERRORSTRING_OK,"0<=kappa1<=1");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPA]);
    return ERRORCOVFAILED;
  }
  if (key->param[KAPPAX]<=0.0){
    strcpy(ERRORSTRING_OK,"kappa2 > 0");
    sprintf(ERRORSTRING_WRONG,"%f",key->param[KAPPAX]);
    return ERRORCOVFAILED;
  }
  return 0;  
}



Real wave(Real x, Real *p) {
  Real y;
  if (x==0.0) { return p[SILL];}
  y=x*p[INVSCALE];
  return p[VARIANCE]*sin(y)/y;
}
Real Scalewave(Real *p,int scaling) {return 0.302320850755833;}
Real spectralwave(Real *p END_WITH_GSLRNG) {
  Real x;  x=UNIFORM_RANDOM; return (p[INVSCALE] * sqrt(1-x*x));
}


Real cubic(Real x, Real *p) 
{ ///
  Real y,y2;
  if (x==0.0) {return p[SILL];}
  if ((y=fabs(x*p[INVSCALE]))>=1.0) {return 0.0;}
  y2=y * y;
  return p[VARIANCE]* (1.0- (((-0.75 * y2  +3.5) * y2 -8.75) * y + 7)*y2);
}
Real TBM3cubic(Real x, Real *p) 
{ ///
  Real y,y2;
  if ((y=fabs(x*p[INVSCALE]))>=1.0) {return 0.0;}
  y2 = y * y;
  return p[VARIANCE]*(1.0 + y2 * (-21.0 + y * (35.0 + y2 * (-21.0 + 6.0 * y2))));
}
Real Scalecubic(Real *p,int scaling) {return 1.44855683156829;}


Real penta(Real x, Real *p) 
{ ///
  Real y,y2;
  if (x==0.0) { return p[SILL];}
  if ((y=fabs(x*p[INVSCALE]))>=1.0) { return 0.0; }
  y2=y * y;
  return p[VARIANCE]* 
    (1.0 + y2 * (-7.333333333333333 
		 + y2 * (33.0 +
			 y * (-38.5 + 
			      y2 * (16.5 + 
				    y2 * (-5.5 + 
					  y2 * 0.833333333333333))))));
}
Real TBM3penta(Real x, Real *p) 
{ ///
  Real y,y2;
  if ((y=fabs(x*p[INVSCALE]))>=1.0) { return 0.0; }
  y2 = y * y;
  return p[VARIANCE]*
    (1.0 + y2 * (-22.0 
		 + y2 * (165.0 
			 + y * (-231.0 
				+ y2 * (132.0 
					+ y2 * (-55.0 
						+ 10.0 * y2))))));
}
Real Scalepenta(Real *p,int scaling) {return 1.6552838957365;}

// see /home/martin/article/C/RFCovBrownian.cc




