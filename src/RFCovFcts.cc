/*
 Authors 
 Yindeng, Jiang, jiangyindeng@gmail.com
 Martin Schlather, schlath@hsu-hh.de

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * when defining your own function, make sure that the covariance function 
   itself allows for an additional nugget effect (spectral measures and tbm 
   operators don't)     
 * VARIANCE, SCALE may not be used here since these parameters are already 
   considered elsewhere

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2005 Martin Schlather

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

#define RANGE_EPSILON 1E-20
#define UNIT_EPSILON 1E-14

#include <math.h>
#include <assert.h>
#include "RFsimu.h"
#include "RFCovFcts.h"

//  {min, max} x {theor, pract}  x length(param) 
static double range_stable[4] = {0, 2, 0.06, 2};
static double range_whittle[4]= {0, RF_INF, 1e-2, 10.0};
static double range_cauchy[4] = {0, RF_INF, 0.09, 10.0};
static double range_genCauchy[8] = {0, 2, 0.05, 2, 
				    0, RF_INF, 0.05, 10.0};


double interpolate(double y, double *stuetz, int nstuetz, int origin,
		 double lambda, int delta)
{
  int index,minindex,maxindex,i;
  double weights,sum,diff,a;
  index  = origin + (int) y;
  minindex = index - delta; if (minindex<0) minindex=0;
  maxindex = index + 1 + delta; if (maxindex>nstuetz) maxindex=nstuetz;
  weights = sum = 0.0;
  for (i=minindex;i<maxindex;i++) {    
    diff = y + (double) (index-i);
    a    = exp( -lambda * diff * diff);
    weights += a;
    sum += a * stuetz[i];  // or  a/stuetz[i]
  }
  return (double) (weights/sum); // and then   sum/weights       
}

//// NOTE : `*p' may not be changed by any of the functions!
/* constant model */
double constant(double *x, double *p, int effectivedim){
  return 1.0;
}
double TBM2constant(double *x, double *p, int effectivedim) 
{
  return 1.0;
}
double TBM3constant(double *x, double *p, int effectivedim){
   return 1.0;
}
double Dconstant(double *x, double *p, int effectivedim){
  return 0.0;
}
void rangeconstant(int dim, int *index, double* range){ 
  *index = -1; 
}
int checkconstant(double *param, int timespacedim, SimulationType method) {
    return NOERROR;
}
void infoconstant(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM; 
  *CEbadlybehaved=false;
}


/* exponential model */
double exponential(double *x, double *p, int effectivedim){
  return exp(-fabs( *x));
}
double Scaleexponential(double *p, int scaling){ return MINUSINVLOG005; } 
double TBM2exponential(double *x, double *p, int effectivedim) 
{
  double y;
  if (*x==0.0) {return 1.0;}
  y = fabs(*x);
  return 1.0 - PIHALF * y * I0mL0(y);
}
double TBM3exponential(double *x, double *p, int effectivedim){
   register double y;
   y = fabs( *x);
   return (1.0-y)*exp(-y);
}
double Dexponential(double *x, double *p, int effectivedim){
//  printf("\n ******* dexp %f %f\n", *x, - exp(-fabs( *x)));
  return - exp(-fabs( *x));
}
double spectralexponential(double *p ) { /* see Yaglom ! */
  double register y;
  y = 1.0 - UNIFORM_RANDOM;
  return sqrt(1.0 / (y * y) - 1.0);
}
void rangeexponential(int dim, int *index, double* range){ 
  *index = -1; 
}
int checkexponential(double *param, int timespacedim, SimulationType method) {
  if (method==CircEmbedIntrinsic || method==CircEmbedCutoff) {
    if (timespacedim>2) 
    {
      strcpy(ERRORSTRING_OK,"genuine dim<=2");
      sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
      return ERRORCOVFAILED;
    }
  }
  return NOERROR;
}

int hyperexponential(double radius, double *center, double *rx,
		     int dim, bool simulate, 
		     double** Hx, double** Hy, double** Hr)
{
  // lenx : half the length of the rectangle
  // center   : center of the rectangle
  // simulate=false: estimated number of lines returned;
  // simulate=true: number of simulated lines returned;
  // hx, hy : direction of line
  // hr     : distance of the line from the origin
  // rectangular area where center gives the center 
  // 
  // the function expects scale = 1;
  double lambda, phi, lx, ly, *hx, *hy, *hr;
  long i, p, q;
  int k, error;
  
  if (dim==2) {
    // we should be in two dimensions
    // first, we simulate the lines for a rectangle with center (0,0)
    // and half the side length equal to lenx
    lx = rx[0];
    ly = rx[1];
    lambda = TWOPI * radius * 0.5; /* total, integrated, intensity */
    //    0.5 in order to get scale 1
    if (!simulate) return (int) lambda;
    assert(*Hx==NULL);
    assert(*Hy==NULL);
    assert(*Hr==NULL);
    p = (long) rpois(lambda);
    if ((hx=*Hx=(double*) malloc(sizeof(double) * (p + 8 * sizeof(int))))==NULL){
      error=ERRORMEMORYALLOCATION; goto ErrorHandling;
    } 
     // see also bits in RFhyperplan.cc, line 437 about.
    if ((hy=*Hy=(double*) malloc(sizeof(double) * (p + 8 * sizeof(int))))==NULL){
      error=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }  
    if ((hr=*Hr=(double*) malloc(sizeof(double) * (p + 8 * sizeof(int))))==NULL){
      error=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }  

    /* creating the lines; some of the lines are not relevant as they
       do not intersect the rectangle investigated --> k!=4
       (it is checked if all the corners of the rectangle are on one 
       side (?) )
    */
    q=0;
    for(i=0; i<p; i++) {
      phi = UNIFORM_RANDOM * TWOPI;
      hx[q] = cos(phi); 
      hy[q] = sin(phi);    
      hr[q] = UNIFORM_RANDOM * radius;
      k = (hx[q] * (-lx) + hy[q] * (-ly) < hr[q]) +
          (hx[q] * (-lx) + hy[q] * ly < hr[q]) +
          (hx[q] * lx + hy[q] * (-ly) < hr[q]) +
          (hx[q] * lx + hy[q] * ly < hr[q]);
      if (k!=4) { // line inside rectangle, so stored
	// now the simulated line is shifted into the right position 
	hr[q] += center[0] * hx[q] + center[1] * hy[q]; 
	q++; // set pointer for storing to the next element
      }
    }
  } else { 
    // dim = 1  --  not programmed yet
    assert(false);
  }
  return q;

 ErrorHandling:
  PRINTF("error=%d\n", error);
  assert(false);
}
void infoexponential(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved=false;
}


//		      CutOffKnowlege *cutoffknowledge, double *cutoffA,
//		      double *intr_r_theo, double *intr_thresh_theo,
//		      int *intr_gridsize, int *intr_diameter,
//		      double *intr_table) 


/* derivative of exponential: qexponential */
double qexponential(double *x,double *p, int effectivedim){
  register double y;
  y = exp(-fabs( *x));
  return y * (2.0  - p[KAPPA] * y) / (2.0 - p[KAPPA]);
}
double Scaleqexponential(double *p,int scaling){
  return -1.0 / 
    log( (1.0 - sqrt(1.0 - p[KAPPA] * (2.0 - p[KAPPA]) * 0.05)) / p[KAPPA]);
} 
int checkqexponential(double *param, int timespacedim, SimulationType method){
  if ((param[KAPPA]<0.0) || (param[KAPPA]>1.0)) {
    strcpy(ERRORSTRING_OK,"kappa in [0,1]");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return NOERROR;  
}
double TBM3qexponential(double *x, double *p, int effectivedim){
   register double y;
   y = exp(-fabs( *x));
   return (y * (2.0  - p[KAPPA] * y) + 
	   fabs( *x) * y * (y * (p[KAPPA] * y - 1.0) * 2.0)) / 
     (2.0 - p[KAPPA]);
}
double Dqexponential(double *x, double *p, int effectivedim) {
  register double y;
  y = exp(-fabs( *x));
  return y * (p[KAPPA] * y - 1.0) * 2.0 / (2.0 - p[KAPPA]);
}
void rangeqexponential(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  static double r[4] = {0, 1, 0, 1};
  memcpy(range, r, sizeof(double) * 4);
}
void infoqexponential(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = (p[KAPPA]>=0 && p[KAPPA]<=1) ? INFDIM : 0; 
  *CEbadlybehaved=false;
}

/* derivative of exponential: damped cosine */
double dampedcosine(double *x, double*p, int effectivedim){
  register double y;
  y = fabs( *x);
  return exp(-y * p[KAPPA]) * cos(y);
}
double Scaledampedcosine(double *p,int scaling){ 
  if (scaling==NATSCALE_EXACT) return 0.0;
  return MINUSINVLOG005; 
} 
double TBM3dampedcosine(double *x, double *p, int effectivedim){
  register double y;
  y = fabs( *x);
  return exp(-p[KAPPA]*y) * ((1.0 - p[KAPPA] * y) * cos(y) - y * sin(y));
}
double Ddampedcosine(double *x, double *p, int effectivedim){
  register double y;
  y = fabs( *x);
  return - exp(-p[KAPPA]*y) * (p[KAPPA] * cos(y) + sin(y));
}
int checkdampedcosine(double *param, int timespacedim, SimulationType method){
  if (timespacedim==3 || timespacedim<3 && method==TBM3){
    if (param[KAPPA]<1.73205080756889) {
      strcpy(ERRORSTRING_OK,
	     "kappa >= sqrt(3) for 3-dimensional simulations and techniques");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  } else if (timespacedim==2 || timespacedim==1 && method==TBM2){
    if (param[KAPPA]<1.0) {
      strcpy(ERRORSTRING_OK,"kappa >= 1.0 for 2 genuine dimensions");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  } else if (timespacedim==1){
    if (param[KAPPA]<0.0) {
      strcpy(ERRORSTRING_OK,"kappa >= 0.0 for 1 genuine dimension");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
  } else { // timespacedim ==4
    return ERRORNOTPROGRAMMED;
  }
  return 0;  
}
void rangedampedcosine(int dim, int *index, double* range){
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
void infodampedcosine(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim =
    (p[KAPPA] < 0) ? 0 : (p[KAPPA] < 1) ? 1 : (p[KAPPA] < M_SQRT_3) ? 2 : 3;
  *CEbadlybehaved=false;
}


/* Class of "spherical" models  */ 
/* circular model */
double circular(double *x, double *p, int effectivedim)
{
  double y;
  if ((y=fabs( *x))>1.0) return 0.0; 
  return  1.0 - (2.0 * (y * sqrt(1.0- y * y) + asin(y))) * INVPI;
}
double Scalecircular(double *p,int scaling) {return 1.138509531721630274603;}
// spectral measure, see Lantue !! 

double Dcircular(double *x, double *p, int dim){
  register double y;
  if ((y=*x * *x) >= 1.0) {return 0.0;} 
  return -4 * INVPI * sqrt(1 - y);
}
void rangecircular(int dim, int *index, double* range){
  if (dim<=2) *index=-1; else *index=-2;
}
void infocircular(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim= 2;
  *CEbadlybehaved=false;
}
int checkcircular(double *param, int timespacedim, SimulationType method){
  if (timespacedim>2) {
    strcpy(ERRORSTRING_OK,"genuine dimension must be less than or equal to 2");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;  
}



/* spherical model */ 
double spherical(double *x, double *p, int effectivedim){
  register double y;
  if ((y= fabs( *x))>1.0) {return 0.0;}
  return (1.0+y*0.5*(y*y-3.0));
}
double Scalespherical(double *p,int scaling){ return 1.23243208931941;}
double TBM2spherical(double *x, double *p, int effectivedim){
  register double y, y2;
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
double TBM3spherical(double *x, double *p, int effectivedim){
  register double y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return (1.0 + (-3.0 + 2.0 * y * y) * y);
}
double Dspherical(double *x, double *p, int effectivedim){
  register double y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return  1.5 * (y * y - 1.0);
}
void rangespherical(int dim, int *index, double* range){
  if (dim<=3) *index=-1; else *index=-2;
}
void infospherical(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved=false;
}
int checkspherical(double *param, int timespacedim, SimulationType method){
  if (timespacedim > 3) {
    strcpy(ERRORSTRING_OK,"genuine dimension must be less than or equal to 3");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;  
}


/* power model */ 
double power(double *x, double *p, int effectivedim){
  register double y;
  if ((y=fabs( *x))>1.0) {return 0.0;}
  return  pow(1.0 - y, p[KAPPA]);
}
double Scalepower(double *p,int scaling){ 
  return 1.0 / (1.0 - pow(0.05,1/p[KAPPA]));
}
double TBM2power(double *x, double *p, int effectivedim){
  // only kappa=2 up to now !
  register double y;
  y=fabs( *x);
  if (y>1.0)
    {return (1.0 - 2.0 * y *(asin(1.0/y) - y + sqrt(y*y-1.0) ));} 
  return (1.0 - y * (PI - 2.0 * y));
}
double TBM3power(double *x, double *p, int effectivedim){
  register double y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return  (1.0 - y - y * p[KAPPA]) *
    pow(1.0 - y, p[KAPPA]-1.0);
}
double Dpower(double *x, double *p, int effectivedim){
  register double y;
  if ((y=fabs( *x))>1.0) {return 0.0;} 
  return  -  p[KAPPA] * pow(1.0 - y, p[KAPPA]-1.0);
}
int checkpower(double *param, int timespacedim, SimulationType method) {
  int error;
  error = NOERROR;
  switch (method) {
  case TBM2 : 
    if (param[KAPPA]<3/2) { // not that kappa==2 if fine also for 
      //                              3 dimensions!
      strcpy(ERRORSTRING_OK,"kappa>=3/2 in case of TBM2");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVFAILED;
    }
    if (param[KAPPA]!=2) error=ERRORCOVNUMERICAL;
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
      sprintf(ERRORSTRING_WRONG,"kappa=%f genuine dim=%d",param[KAPPA],
	      timespacedim);
      return ERRORCOVFAILED;
    }
  }
  return error;
}
// range definition:
// 0: min, theory, 1:max, theory
// 2: min, practically 3:max, practically
void rangepower(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  range[0] = range[2] = 0.5 * ((double) dim + 1);
  range[1] = RF_INF;
  range[3] = 20.0;
}
void infopower(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = (int) (2.0 * p[KAPPA] - 1.0);
  *CEbadlybehaved=false;
}

/* stable model */
double stable(double *x,double *p, int effectivedim){
  if (*x==0.0) {return 1.0;}
  return  (exp(-pow(fabs( *x), p[KAPPA])));
}
double Scalestable(double *p,int scaling){return pow(MINUSINVLOG005,1/p[KAPPA]);}
double TBM3stable(double *x, double *p, int effectivedim){
  register double y;
  y = pow(fabs( *x),p[KAPPA]);
  return exp(-y) * (1 - p[KAPPA] * y);
}
double Dstable(double *x, double *p, int effectivedim){
  register double y, z;
  if ( (z = fabs( *x)) == 0.0) 
    return ((p[KAPPA] > 1.0) ? 0.0 : (p[KAPPA] < 1.0) ? INFTY : 1.0);
  y = pow(z, p[KAPPA] - 1.0);
  return -p[KAPPA] * y * exp(-y * z);
}
/* stable: second derivative at t=1 */
double DDstable(double *x, double*p, int effectivedim) 
{
  double y, xkappa, z;
  if ( (z = fabs( *x)) == 0.0) return ((p[KAPPA] != 1.0) ? INFTY : 1.0);
  y = pow(z, p[KAPPA] - 2.0);
  xkappa = y *z * z;
  return p[KAPPA] * (1.0 - p[KAPPA] + p[KAPPA] * xkappa) * y * exp(-xkappa);
}

int checkstable(double *param, int timespacedim, SimulationType method) {
  if ((param[KAPPA]<=0) || (param[KAPPA]>2.0)) {
    strcpy(ERRORSTRING_OK,"0<kappa<=2");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED; 
  }
  if (method==CircEmbedIntrinsic || method==CircEmbedCutoff) {
    if (timespacedim>2) 
    {
      strcpy(ERRORSTRING_OK,"genuine total dim<=2");
      sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
      return ERRORCOVFAILED;
    }
//  if (method==CircEmbedCutoff)
//       if (param[KAPPA]>1.0) {
// 	strcpy(ERRORSTRING_OK,"0<kappa<=1");
// 	sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
// 	return ERRORCOVFAILED; 
//       }
  }
  return 0;
}
void rangestable(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_stable, sizeof(double) * 4);
}
void infostable(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = p[KAPPA] == 2.0;
}


/* Whittle-Matern or Whittle or Besset */ 
double WhittleMatern(double *x, double *p, int effectivedim) 
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
    exp(kappa * log(0.5 * y) - loggamma + log(bessel_k(y, kappa, 2.0)) - y);
}
double ScaleWhittleMatern(double *p,int scaling){
  // it is working reasonably well if kappa is in [0.001,100]
  // happy to get any better suggestion!!
  static int nstuetz = 19;
  static double stuetz[]=
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

double TBM2WhittleMatern(double *x, double *p, int effectivedim) 
{
  if (p[KAPPA]==0.5) return TBM2exponential(x, p, effectivedim);
  assert(false);
}

double TBM3WhittleMatern(double *x, double *p, int effectivedim)
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

double DWhittleMatern(double *x, double *p, int effectivedim)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double kappa=RF_INF;
  static double loggamma;
  register double y;
  if ( *x==0.0) 
    return ((p[KAPPA] > 0.5) ? 0.0 : (p[KAPPA] < 0.5) ? INFTY : 1.253314137);
  y = fabs( *x); 
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    loggamma = lgammafn(kappa);
  }
  return -2.0 * exp(kappa * log(0.5 * y) - loggamma + 
		    log(bessel_k(y, kappa - 1.0, 2.0)) - y);
}

double DDWhittleMatern(double *x, double *p, int effectivedim)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double kappa=RF_INF;
  static double gamma;
  register double y;
  if (*x==0.0)  return ((p[KAPPA] > 1.0) ? 0.5 / (p[KAPPA] - 1.0) : INFTY);
  y = fabs( *x); 
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    gamma = gammafn(kappa);
  }
  return pow(0.5 * y , p[KAPPA] - 1.0) / gammafn(p[KAPPA]) *
    (bessel_k(y, p[KAPPA] - 1.0, 1.0) - y * bessel_k(y, p[KAPPA] - 2.0, 1.0)); 
}

double spectralWhittleMatern(double *p ) { /* see Yaglom ! */
  return sqrt(pow(1.0 - UNIFORM_RANDOM, -1.0 / p[KAPPA]) - 1.0);
}

int checkWhittleMatern(double *param, int timespacedim, SimulationType method) { 
  static double spectrallimit=0.17;
  switch(method) {
      case TBM2 : 
	if (param[KAPPA]!=0.5) {
	  strcpy(ERRORSTRING_OK,"1/2");
	  sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
	  return ERRORCOVNUMERICAL;
	}
	break;
      case SpectralTBM : 
	if ((param[KAPPA]>=spectrallimit)) return 0;
	sprintf(ERRORSTRING_OK,
		"%f<=kappa (the numerical errors are too big for 0<kappa<%f)",
		spectrallimit,spectrallimit);
	sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
	return ERRORCOVFAILED;

      case CircEmbedCutoff: case CircEmbedIntrinsic :
	if (timespacedim>2) 
	{
	  strcpy(ERRORSTRING_OK,"genuine total dim<=2");
	  sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
	  return ERRORCOVFAILED;
	}
//     if (method==CircEmbedCutoff)
//       if (2*param[KAPPA]>1.0) {
// 	strcpy(ERRORSTRING_OK,"0<kappa<=0.5");
// 	sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
// 	return ERRORCOVFAILED; 
//       }
	// break;
	// NOTE NO BREAK !
      default :
	  if ((param[KAPPA]<=0)) {
	      strcpy(ERRORSTRING_OK,"0<kappa");
	      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
	      return ERRORCOVFAILED;
	  }
  }
  return NOERROR;
}
void rangeWhittleMatern(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_whittle, sizeof(double) * 4);
}
void infoWhittleMatern(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = p[KAPPA] >= 1.5;
}



double hyperbolic(double *x, double*p, int effectivedim){ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double kappay;
  double y;
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
double TBM3hyperbolic(double *x, double*p, int effectivedim)
{ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double y;
  double ysq,s,kappa_s,logs;
  if ( *x==0.0) {return 1.0;}
  if (p[KAPPAIII]==0) { // whittle matern
    y = *x * p[KAPPAI];
    return TBM3WhittleMatern(&y, p, effectivedim);
  } 
  if (p[KAPPAI]==0) { //cauchy
    register double y,ha;
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
  kappa_s = kappa * s;
  logs = log(s);  
  return  
    ( exp(logconst + lambda * logs +log(bessel_k(kappa_s,lambda,2.0))-kappa_s)
      - ysq*kappa*exp(logconst + (lambda-1.0)*logs 
		      +log(bessel_k(kappa_s,lambda-1.0,2.0))-kappa_s)
      );
}
double Dhyperbolic(double *x, double*p, int effectivedim)
{ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double y;
  double s,kappa_s,logs;
  if ( *x==0.0) {return 1.0;}
  if (p[KAPPAIII]==0) { // whittle matern
    y = *x * p[KAPPAI];
    return DWhittleMatern(&y, p, effectivedim);
  } 
  if (p[KAPPAI]==0) { //cauchy
    register double y,ha;
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
  kappa_s = kappa * s;
  logs = log(s);  
  return  
    ( 
      - y * kappa*exp(logconst + (lambda-1.0)*logs 
		      +log(bessel_k(kappa_s,lambda-1.0,2.0))-kappa_s)
      );
}
int checkhyperbolic(double *param, int timespacedim, SimulationType method){
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
void rangehyperbolic(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  double r[12] = {0, RF_INF, 0.000001, 10.0,
		RF_NEGINF, RF_INF, -20.0, 20.0,
                0, RF_INF, 0.000001, 10.0
                };
  *index = -1; 
  memcpy(range, r, sizeof(double) * 12);
}
void infohyperbolic(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = true;
}


/* Gneiting's functions */
// alternative to Gaussian */
#define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
#define NumericalScale 0.301187465825
double Gneiting(double *x, double *p, int effectivedim){ 
  register double y,oneMy8;
  if ((y=fabs(*x * NumericalScale))>1.0) {return(0.0);}
  oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
  return ((1.0+y * ( 8.0 + y * (25.0 + 32.0 *y)))*oneMy8);
}
double ScaleGneiting(double *p,int scaling) {return 0.5854160193;}
double TBM3Gneiting(double *x, double *p, int effectivedim){ 
  register double y,oneMy7;
  if ((y=fabs( *x * NumericalScale))>1.0) {return 0.0;}  
  oneMy7 = 1.0-y; oneMy7*=oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);
  return 
    (1.0 + y * (7.0  -  y * (5.0 + y * (147.0 + 384.0 * y))))* oneMy7;
}
double DGneiting(double *x, double *p, int effectivedim){ 
  register double y,oneMy7;
  if ((y=fabs( *x * NumericalScale))>1.0) {return 0.0;}  
  oneMy7 = 1.0-y; oneMy7*=oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);
  return 
    (-y) * ( 22.0 + y * (154.0 + y * 352.0)) * oneMy7 * Sqrt2TenD47;
}
void rangeGneiting(int dim, int *index, double* range){
  if (dim<=3) *index=-1; else *index=-2;
}
void infoGneiting(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}
int checkGneiting(double *param, int timespacedim, SimulationType method){
  if (timespacedim > 3) {
    strcpy(ERRORSTRING_OK,"genuine dimension must be less than or equal to 3");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;  
}


double genGneiting(double *x, double *p, int effectivedim)
{
  register double y, s;
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
double TBM3genGneiting(double *x, double *p, int effectivedim)
{
  register double y, s;
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
double DgenGneiting(double *x, double *p, int effectivedim)
{
  register double y, s;
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
int checkgenGneiting(double *param, int timespacedim, SimulationType method)
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
      (param[KAPPAII] < 2.0 + param[KAPPAI])) { // b >= (d + 2a +1)/2 
    strcpy(ERRORSTRING_OK,"kappa2 <= 2.0 + kappa1");
    sprintf(ERRORSTRING_WRONG, "method=TBM3, kappa1=%f and kappa2=%f",
	    param[KAPPAI], param[KAPPAII]);
    return ERRORCOVFAILED;    
  }
  if ((method==TBM2) &&
      (param[KAPPAII] < 1.5 + param[KAPPAI])) { // b >= (d + 2a +1)/2 
    strcpy(ERRORSTRING_OK,"kappa2 <= 2.0 + kappa1");
    sprintf(ERRORSTRING_WRONG, "method=TBM3, kappa1=%f and kappa2=%f",
	    param[KAPPAI], param[KAPPAII]);
    return ERRORCOVFAILED;    
  }
  if  (param[KAPPAII] < ( 0.5 * timespacedim + 
				 param[KAPPAI] + 0.5)) {
    strcpy(ERRORSTRING_OK, "kappa2 >= (dim + 2*kappa1 + 1)/2 ");
    sprintf(ERRORSTRING_WRONG, "kappa1=%f, kappa2=%f, and genuine dim=%d",
	    param[KAPPAI], param[KAPPAII], timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;
}
void rangegenGneiting(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  range[0] = range[2] = range[1] = range[3] = *index;
  range[4] = range[6] = 0.5 * (double) (dim + 2 * *index + 1);
  range[5] = RF_INF;
  range[7] = 20.0;
  if ((++(*index)) > 3) *index=-1;
  if (dim>3) *index=-2;
}
void infogenGneiting(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = (int) (2.0 * p[KAPPAII] - 1.0 - 2.0 * p[KAPPAI]);
  *CEbadlybehaved = false;
}



/* Gausian model */
double Gauss(double *x, double*p, int effectivedim) {
  return  exp(-  *x * *x);
}
double ScaleGauss(double *p,int scaling) {return SQRTINVLOG005;}
double TBM3Gauss(double *x, double*p, int effectivedim) {
  register double y;
  y =   *x *  *x; 
  return  (1-2.0*y)*exp(-y);
}
double DGauss(double *x, double*p, int effectivedim) {
  register double y;
  y =   fabs( *x); 
  return  -  2.0 * y * exp(- y * y);
}
double spectralGauss(double *p ) {   
  return 2.0  * sqrt(-log(1.0 - UNIFORM_RANDOM));
}
void rangeGauss(int dim, int *index, double* range){
  *index=-1;
}
void infoGauss(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}
int checkGauss(double *param, int timespacedim, SimulationType method) {
    return NOERROR;
}


/* Cauchy models */
double Cauchy(double *x, double *p, int effectivedim){
  return  pow(1+ *x * *x,-p[KAPPA]);
}
double ScaleCauchy(double *p,int scaling) {
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
double TBM2Cauchy(double *x, double *p, int effectivedim){
  register double y2, lpy2;
  y2 = *x * *x; 
  lpy2=1.0+y2;
  switch ((int) (p[KAPPA]*2.0 + 0.1)) {
  case 1 : return 1.0 / lpy2;
  case 3 : return (1.0 - y2)/ (lpy2*lpy2);
  case 5 : return (1.0-y2*(2.0+0.333333333333333333333*y2))/(lpy2*lpy2*lpy2);
  case 7 : lpy2 *= lpy2; return (1.0- y2*(3.0+y2*(1.0+0.2*y2)))/(lpy2 * lpy2);
  default : assert(false);
  }
}
double TBM3Cauchy(double *x, double *p, int effectivedim){
  register double ha;
  ha= *x *  *x;
  return (1.0+ (1.0-2.0*p[KAPPA])*ha) * pow(1.0+ha,-p[KAPPA]-1.0);
}
double DCauchy(double *x, double *p, int effectivedim){
  register double y;
  y = fabs( *x);
  return (-2.0 * p[KAPPA] * y) * pow(1.0 + y * y, -p[KAPPA] - 1.0);
}
double DDCauchy(double *x, double *p, int effectivedim){
  register double ha;
  ha = *x * *x;
  return 2.0 * p[KAPPA] * ((2.0 * p[KAPPA] + 1.0) * ha - 1.0) * 
    pow(1.0 + ha, -p[KAPPA] - 2.0);
}
int checkCauchy(double *param, int timespacedim, SimulationType method){
  switch (method) {
  case TBM2 : 
    // not replaced by numerical evaluation due to bad
    // numerical behaviour?!
    if (timespacedim>2) return ERRORCOVFAILED;
    if ((param[KAPPA]!=0.5) && (param[KAPPA]!=1.5) &&
	(param[KAPPA]!=2.5) && (param[KAPPA]!=3.5)) {
      strcpy(ERRORSTRING_OK,"kappa in {0.5, 1.5, 2.5 ,3.5} (TBM2)");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVNUMERICAL;
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
void rangeCauchy(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_cauchy, sizeof(double) * 4);
}
void infoCauchy(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}

double generalisedCauchy(double *x, double *p, int effectivedim){
  return pow(1.0 + pow(fabs(*x), p[KAPPAI]), -p[KAPPAII]/p[KAPPAI]);
}
double ScalegeneralisedCauchy(double *p,int scaling) {
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
double TBM3generalisedCauchy(double *x, double *p, int effectivedim){
  register double ha;
  ha=pow(fabs( *x), p[KAPPAI]);
  return 
    (1.0+ (1.0-p[KAPPAII])*ha) * pow(1.0 + ha, -p[KAPPAII] / p[KAPPAI] - 1.0);
}
double DgeneralisedCauchy(double *x, double *p, int effectivedim){
  register double ha,y;
  if ((y = fabs(*x))==0.0) 
    return ((p[KAPPAI]>1.0) ? 0.0 : (p[KAPPAI]<1.0) ? -INFTY : -p[KAPPAII]); 
  ha=pow(y, p[KAPPAI] - 1.0);
  return  -  p[KAPPAII] * ha * pow(1.0 + ha * y,-p[KAPPAII] / p[KAPPAI] - 1.0);
}
double DDgeneralisedCauchy(double *x, double *p, int effectivedim){
  register double ha,y;
  if ((y = fabs(*x))==0.0) 
    return ((p[KAPPAI]==2.0) ? p[KAPPAII] * (p[KAPPAII] + 1.0) : INFTY); 
  ha=pow(y, p[KAPPAI]);
  return p[KAPPAII] * ha / (y * y) * (1.0 - p[KAPPAI] + (1.0 + p[KAPPAII]) * ha)
	    * pow(1.0 + ha, -p[KAPPAII] / p[KAPPAI] - 2.0);
}
int checkgeneralisedCauchy(double *param, int timespacedim, SimulationType method){
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
//     if (method==CircEmbedCutoff)
//       if (param[KAPPAI]>1.0) {
// 	strcpy(ERRORSTRING_OK,"0<kappa1<=1");
// 	sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAI]);
// 	return ERRORCOVFAILED; 
//       }
  }
  return 0;
}
void rangegeneralisedCauchy(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_genCauchy, sizeof(double) * 8);
}
void infogeneralisedCauchy(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}

double Cauchytbm(double *x, double *p, int effectivedim){
  register double ha;
  if ( *x==0) {return 1.0;}
  ha=pow(fabs( *x),p[KAPPAI]);
  return (1.0 + (1.0 - p[KAPPAII] / p[KAPPAIII]) * ha) * 
      pow(1.0 + ha, -p[KAPPAII] / p[KAPPAI] - 1.0);
}
double TBM3Cauchytbm(double *x, double *p, int effectivedim){
  register double bg,ha;
  ha=pow(fabs( *x),p[KAPPAI]);
  bg=p[KAPPAII]/p[KAPPAIII];
  return  
    (1 + ha * (1-bg*(1+p[KAPPAI])+(1-p[KAPPAII]) * (1+(1-bg)*ha)))*
    pow(1+ha,-p[KAPPAII]/p[KAPPAI]-2.0);
}
double DCauchytbm(double *x, double *p, int effectivedim){
  register double y,ha;
  if ((y = fabs(*x)) == 0.0) return 0.0; // WRONG VALUE, but multiplied 
  //                                        zero anyway
  ha=pow(y,p[KAPPAI] - 1.0);
  return  
    p[KAPPAII] *  ha * (-1.0 - p[KAPPAI]/p[KAPPAIII] + 
		       ha * y * (p[KAPPAII]/p[KAPPAIII] - 1.0)) *
    pow(1.0 + ha * y,-p[KAPPAII]/p[KAPPAI]-2.0);
}
int checkCauchytbm(double *param, int timespacedim, SimulationType method){
  if ((method==TBM3) && (3.0 > param[KAPPAIII])) {
    strcpy(ERRORSTRING_OK,"3 <= kappa3 and method=TBM3");
    sprintf(ERRORSTRING_WRONG,
	    "method=TBM3 and kappa3=%f",param[KAPPAIII]);
    return ERRORCOVFAILED;    
  }
  if ((method==TBM2) && (2.0 > param[KAPPAIII])) {
    strcpy(ERRORSTRING_OK,"2 <= kappa3 and method=TBM2");
    sprintf(ERRORSTRING_WRONG,
	    "method=TBM2 and kappa3=%f", param[KAPPAIII]);
    return ERRORCOVFAILED;    
  }
  if (timespacedim > param[KAPPAIII]) {
    strcpy(ERRORSTRING_OK,"dim <= kappa3");
    sprintf(ERRORSTRING_WRONG,
	    "kappa3=%f and  genuine dim=%d",param[KAPPAIII], timespacedim);
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
void rangeCauchytbm(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = -1; 
  memcpy(range, range_genCauchy, sizeof(double) * 8);
  range[8] = range[10] = (double) dim;
  range[9] = RF_INF;
  range[11] = (double) dim + 10.0;
}
void infoCauchytbm(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}


/* Bessel function */ 
double Bessel(double *x,double *p, int effectivedim){
  static double kappa=RF_INF;
  static double gamma;
  double y;
  if  ( *x==0.0) {return 1.0;} 
  y = fabs( *x);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    gamma = gammafn(kappa+1.0);
  }
  return  gamma  * pow(0.5 * y,-kappa) * bessel_j(y,kappa);
}
double spectralBessel(double *p ) {
  return p[KAPPA]==0 ? 1.0 : 
    (sqrt(1-pow(UNIFORM_RANDOM,1.0/p[KAPPA])));
}
int checkBessel(double *param,  int timespacedim, SimulationType method){
  // Whenever TBM3Bessel exists, add further check against too small kappa!
  if (param[KAPPA]<0.5*((double) timespacedim) - 1.0) {
    strcpy(ERRORSTRING_OK,"kappa >= dim/2 - 1");
    sprintf(ERRORSTRING_WRONG,"kappa=%f and genuine dim=%d",
	    param[KAPPA], timespacedim);
    return ERRORCOVFAILED;
  }
  if (method==SpectralTBM && param[KAPPA]<0.0) {
    strcpy(ERRORSTRING_OK,"spectral TBM and kappa>=0");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return 0;
}
void rangeBessel(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract }
  *index = -1; 
  range[2] = 0.0001 + (range[0] =  0.5 * ((double) dim - 2.0));
  range[1] = RF_INF;
  range[3] = 10.0;
}
void infoBessel(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = (int) (2.0 * p[KAPPA] + 2.0);
  *CEbadlybehaved = 2;
}


double wave(double *x, double *p, int effectivedim) {
  if ( *x==0.0) { return 1.0;}
  return sin( *x)/ *x;
}
double Scalewave(double *p,int scaling) {return 0.302320850755833;}
double spectralwave(double *p ) {
  double x;  x=UNIFORM_RANDOM; return sqrt(1- x * x);
}
void rangewave(int dim, int *index, double* range){
  if(dim<=3) *index=-1; else *index=-2;
}
void infowave(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = 2;
}
int checkwave(double *param, int timespacedim, SimulationType method){
  if (timespacedim > 3) {
    strcpy(ERRORSTRING_OK,"genuine dimension must be less than or equal to 3");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;  
}


double cubic(double *x, double *p, int effectivedim)
{ ///
  double y, y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return  (1.0 + (((0.75 * y2 - 3.5) * y2 + 8.75) * y - 7) * y2);
}
double TBM3cubic(double *x, double *p, int effectivedim)
{ ///
  double y, y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return (1.0 + y2 * (-21.0 + y * (35.0 + y2 * (-21.0 + 6.0 * y2))));
}
double Dcubic(double *x, double *p, int effectivedim) 
{ ///
  double y,y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return  y * (-14.0 + y * (26.25 + y2 * (-17.5 + 5.25 * y2)));
}
double Scalecubic(double *p,int scaling) {return 1.44855683156829;}
void rangecubic(int dim, int *index, double* range){
  if(dim<=3) *index=-1; else *index=-2;
}
void infocubic(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}
int checkcubic(double *param, int timespacedim, SimulationType method){
  if (timespacedim > 3) {
    strcpy(ERRORSTRING_OK,"genuine dimension must be less than or equal to 3");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;  
}

double penta(double *x, double *p, int effectivedim)
{ ///
  double y,y2;
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
double TBM3penta(double *x, double *p, int effectivedim)
{ ///
  double y,y2;
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
double Dpenta(double *x, double *p, int effectivedim)
{ ///
  double y,y2;
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
double Scalepenta(double *p,int scaling) {return 1.6552838957365;}
void rangepenta(int dim, int *index, double* range){
 if(dim<=3) *index=-1; else *index=-2;
}
void infopenta(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}
int checkpenta(double *param, int timespacedim, SimulationType method){
  if (timespacedim > 3) {
    strcpy(ERRORSTRING_OK,"genuine dimension must be less than or equal to 3");
    sprintf(ERRORSTRING_WRONG,"%d", timespacedim);
    return ERRORCOVFAILED;
  }
  return 0;  
}


/* Tilmann Gneiting's space time models */

/* following function is unused! -- it has been nsst3*/
double spacetime3(double *x,double *p, int effectivedim){
  // turning bands modified spacetime1 covariance function
  assert(effectivedim==2);
  return spacetime1(x,p,effectivedim) + fabs(x[0]) * Dspacetime1(x,p,effectivedim);
}


double InvSqrtPsi(double x, double a, double b, int c) {
  double y;
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
double spacetime1(double *x, double *p, int effectivedim){
  double y, z, invsqrtpsi;
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
double TBM2spacetime1(double *x, double *p, int effectivedim){
  double y, z, invsqrtpsi;
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
double TBM3spacetime1(double *x, double *p, int effectivedim){
  assert(false); // should never happen
  return spacetime1(x, p, effectivedim) +
    Dspacetime1(x, p, effectivedim) * fabs(x[0]);
}
double Dspacetime1(double *x, double *p, int effectivedim){
  double y, z, invsqrtpsi;
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
int checkspacetime1(double *param, int timespacedim, SimulationType method) {
  int error;
  // 1 : stable
  // 2 : whittle
  // 3 : cauchy (not generalised)
  // first parameter(s) for phi
  // then choice of phi; then two parameters for psi, then choice of psi
  error = NOERROR;
  if (timespacedim<=1) {
    strcpy(ERRORSTRING_OK,"genuine total dim>=2");
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
    error=checkstable(param, timespacedim, method);
    if (error==0 && method==TBM2)  {
      strcpy(ERRORSTRING_OK,"TBM2 and kappa2=2,3");
      sprintf(ERRORSTRING_WRONG,"%f", param[KAPPAII]);
      error = ERRORCOVNUMERICAL;
    }
    break;
  case 2 : error=checkWhittleMatern(param, timespacedim, method); break;
  case 3 : error=checkCauchy(param, timespacedim, method); break;
  default : 
    strcpy(ERRORSTRING_OK,"kappa2=1,2,3");
    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPAIII]);
    return ERRORCOVFAILED;
  }
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
  if (timespacedim - 1 > param[KAPPAVI]) {
    strcpy(ERRORSTRING_OK,"kappa6>=spatial_dim");
    sprintf(ERRORSTRING_WRONG,"%f for genuine spatial dim=%d",
	    param[KAPPAVI], timespacedim-1);
    return ERRORCOVFAILED;
  }
  return error;
}
void rangespacetime1(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract }
  double *r;
  if ((*index<=0) || (*index>9)) { // see also last line
    int i; for (i=0; i<24; i++) range[i]= RF_NAN; *index=-1; return;
  }
  switch ((*index-1) % 3) {
  case 0 : r=range_stable;  break;
  case 1 : r=range_whittle; break;
  case 2 : r=range_cauchy;  break;
  default: assert(false);
  }
  memcpy(range, r, sizeof(double) * 4);
    
  range[4] = range[6] =
    range[5] = range[7] = (double) (1 + (*index - 1) % 3);

  range[8] = 0;
  range[9] = range[11] = 2;
  range[10] = RANGE_EPSILON;

  range[12] = 0;
  range[13] = 1;
  range[14] = RANGE_EPSILON;
  range[15] = 1.0 - UNIT_EPSILON;
  
  range[16] = range[17] = 
    range[18] = range[19] =  (double) (1 + (*index - 1) / 3);
  
  range[20] = range[22] = (double) dim;
  range[21] = RF_INF;
  range[23] = (double) dim + 10.0;

  if ( (++(*index)) > 9) *index=-1;
}
void infospacetime1(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = p[KAPPAII]==3 || (p[KAPPAII]==2 && p[KAPPAI]>1.5);
}


/* Tilmann Gneiting's space time models, part II*/
double spacetime2(double *x,double *p, int effectivedim){
  double y, z, invsqrtpsi;
  assert(effectivedim==2);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPAIV], p[KAPPAV], (int)p[KAPPAVI]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPAIII]) {
  case 1 : z = generalisedCauchy(&y, p, 1); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPAVII]) * z;
}
double TBM3spacetime2(double *x, double *p, int effectivedim){
  assert(false); // should never happen
  return spacetime2(x, p, effectivedim) +
    Dspacetime2(x, p, effectivedim) * fabs(x[0]);
}
double Dspacetime2(double *x, double *p, int effectivedim){
  double y, z, invsqrtpsi;
  assert(effectivedim==2);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPAIV], p[KAPPAV], (int)p[KAPPAVI]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPAIII]) {
  case 1 : z = DgeneralisedCauchy(&y, p, 1); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPAV] + 1.0) * z;
}
int checkspacetime2(double *param, int timespacedim, SimulationType method) {
  int error;
  // 1 : generalisedcauchy 
  // first parameter(s) for phi
  // then choice of phi; then two parameters for psi, then choice of psi
  error = NOERROR;
  if (timespacedim<=1) {
    strcpy(ERRORSTRING_OK,"genuine total dim>=2");
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
    strcpy(ERRORSTRING_OK,"kappa7 >= dim");
    sprintf(ERRORSTRING_WRONG,"kappa7=%d and genuine dim=%d",
	    (int) param[KAPPAVII], timespacedim);
  }
  return error;
}
void rangespacetime2(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract }
  double *r;
  if ((*index<=0) || (*index>3)) { // see also last line
    int i; for (i=0; i<28; i++) range[i]= RF_NAN; *index=-1; return;
  }
  switch ((*index-1) % 1) {
  case 0 : r = range_genCauchy;  break;
  default: assert(false);
  }
  memcpy(range, r, sizeof(double) * 8);

  range[8] = range[10] =
    range[9] = range[11] = (double) (1 + (*index - 1) % 1);

  range[12] = 0;
  range[13] = range[15] = 2;
  range[14] = RANGE_EPSILON;

  range[16] = 0;
  range[17] = 1;
  range[18] = RANGE_EPSILON;
  range[19] = 1.0 - UNIT_EPSILON;

  range[20] = range[21] = 
    range[22] = range[23] =  (double) (1 + (*index - 1) / 1);

  range[24] = range[26] = (double) dim;
  range[25] = RF_INF;
  range[27] = (double) dim + 10.0;
  
  if ( (++(*index)) > 3) *index= -1;
}
void infospacetime2(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = true;
}

/* locally defined functions */
// Brownian motion 
double fractalBrownian(double*x, double *p, int effectivdim) {
  return - pow(fabs(*x), p[KAPPA]);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
//begin
/* fractalBrownian: first derivative at t=1 */
double DfractalBrownian(double *x, double*p, int effectivedim) 
{
    return -p[KAPPA] * pow(fabs(*x), p[KAPPA] - 1.0);
}
/* fractalBrownian: second derivative at t=1 */
double DDfractalBrownian(double *x, double*p, int effectivedim)  
{
    return -p[KAPPA] * (p[KAPPA]-1.0) * pow(fabs(*x), p[KAPPA] - 2.0);
}

int checkfractalBrownian(double *param, int timespacedim, SimulationType method){
  if ((timespacedim>3) && (method!=Nothing)) {
    strcpy(ERRORSTRING_OK,"genuine total dim<=3");
    sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
    return ERRORCOVFAILED;
  }
  if ((param[KAPPA]<=0.0) || (param[KAPPA]>=2.0)) {
    strcpy(ERRORSTRING_OK,"kappa in (0,2)");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
    return ERRORCOVFAILED;
  }
  if ((method!=Nothing) && (method!=CircEmbedIntrinsic)) {
    strcpy(ERRORSTRING_OK,"method=CircEmbedIntrinsic");
    sprintf(ERRORSTRING_WRONG,"%s",METHODNAMES[method]);
    return ERRORCOVFAILED;
  }
  return NOERROR;  
}
void rangefractalBrownian(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  if (dim<=3) *index = -1; else *index=-2;
  double r[4] = {0.0, 2.0, RANGE_EPSILON, 2.0 - UNIT_EPSILON};
  memcpy(range, r, sizeof(double*) * 4);
}
void infofractalBrownian(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = false;
}


double fractGauss(double *x, double *p, int effectivedim){
  register double y;
  if ((y = fabs(*x))==0.0) return 1.0;
  return 0.5 * (pow(fabs(y + 1.0), p[KAPPA])  
		- 2.0 * pow(y, p[KAPPA]) 
		+ pow(fabs(y - 1.0), p[KAPPA]) 
		);
}
int checkfractGauss(double *param, int timespacedim, SimulationType method) {
  if ((timespacedim!=1)) {
    strcpy(ERRORSTRING_OK,"genuine dim=1");
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
void rangefractGauss(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  if (dim==1) *index = -1; else *index=-2;
  double r[4] = {0, 2, RANGE_EPSILON, 2};
  memcpy(range, r, sizeof(double) * 4);
}
void infofractGauss(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 1;
  *CEbadlybehaved = false;
}




/* local-global distinguisher */
double lgd1(double *x, double*p, int effectivedim) {
  double y;
  if ((y=fabs(*x)) < 1) return 1.0 - p[KAPPAII] / (p[KAPPAI] + p[KAPPAII]) 
	     * exp(p[KAPPAI] * log(y));
  else return p[KAPPAI] / (p[KAPPAI] + p[KAPPAII]) 
	 * exp( -p[KAPPAII] * log(y));
}
double Scalelgd1(double *p,int scaling) {
  if (19 * p[KAPPAI] < p[KAPPAII]) 
    return exp( log(0.95 * (p[KAPPAI] + p[KAPPAII]) / p[KAPPAII]) / p[KAPPAI]); 
  else return exp(log(0.05 * (p[KAPPAI] + p[KAPPAII]) / p[KAPPAI])/p[KAPPAII]);
}
double Dlgd1(double *x, double *p, int dim){
  double y, pp;
  if ( (y=fabs(*x)) == 0) return 0; // falscher Wert, aber sonst gibt NAN-Fehler
  pp = ( (y < 1) ? p[KAPPAI] : -p[KAPPAII] ) - 1.0;
  return - p[KAPPAI] * p[KAPPAII] / (p[KAPPAI] + p[KAPPAII]) * exp(pp * y);
}
int checklgd1(double *param, int timespacedim, SimulationType method) {
  if ((timespacedim>2) && (method!=Nothing)) {
    strcpy(ERRORSTRING_OK, "genuine dim<=2");
    sprintf(ERRORSTRING_WRONG,"%d",timespacedim);
    return ERRORCOVFAILED;
  }
  if (((param[KAPPAI]<=0) || (param[KAPPAI]>1.0)) && (timespacedim==1)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=1, genuine dim=1");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPAI]);
    return ERRORCOVFAILED;
  }
  if (((param[KAPPAI]<=0) || (param[KAPPAI]>0.5)) && 
      (timespacedim==2 || timespacedim==1 && method==TBM2)) {
    strcpy(ERRORSTRING_OK,"0<kappa1<=1/2, genuine dim=2");
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
static double range_lgd1[8] = {0.0, 1.0, 0.01, 1.0, 
			     0.0, RF_INF, 0.01, 20.0};
void rangelgd1(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  if (dim<=2) *index = -1; else *index=-2;
  memcpy(range, range_lgd1, sizeof(double) * 8);
  if (dim==2) range[1]=range[3]=0.5;
}
void infolgd1(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = (p[KAPPA] <= 0) ? 0 : (int) (2.0 * (1.5 - p[KAPPA]));
  *CEbadlybehaved = true;
}


/* FD model */
double FD(double *x,double *p, int effectivedim){
  static double dold=RF_INF;
  static double kold, sk;
  double y, d, k, skP1;
  d = p[KAPPAI];
  y = fabs(*x);
  k = (double) trunc((double) y);
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
int checkFD(double *param, int timespacedim, SimulationType method) {
  if (param[KAPPAI] < -0.5 || param[KAPPAI] >= 0.5) {
    strcpy(ERRORSTRING_OK,"-0.5 <= kappa < 0.5");
    sprintf(ERRORSTRING_WRONG, "%f",param[KAPPAI]);
    return ERRORCOVFAILED;    
  } 
  if (timespacedim>1) {
    strcpy(ERRORSTRING_OK, "genuine dim=1");
    sprintf(ERRORSTRING_WRONG, "%d", timespacedim);
    return ERRORCOVFAILED;
  } 
  return 0;
}
void rangeFD(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  if (dim==1) *index = -1; else *index=-2;
  range[2] = (range[0] = -0.5) + RANGE_EPSILON;
  range[4] = (range[1] = 0.5) - RANGE_EPSILON;
}
void infoFD(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 1;
  *CEbadlybehaved = false;
}



/* iaco cesare model */
double IacoCesare(double *x, double *p, int dim){
    double s;
    int d;
    for (s=0.0, d=dim-2; d>=0; d--) s += x[d] * x[d];
    return pow(1.0 + pow(s, 0.5 * p[KAPPAI]) + 
	       pow(fabs(x[dim-1]),  p[KAPPAII]),
	       - 0.5 * p[KAPPAIII]); 
}
//double ScaleIacoCesare(double *p, int scaling) { return 1.0; } 
void rangeIacoCesare(int dim, int *index, double* range){
    static double range_iacocesare[8]= 
	{1, 2, 1, 2,
	 1, 2, 1, 2};
    memcpy(range, range_iacocesare, sizeof(double) * 8);
    range[8] = range[10] = 0.5 * dim;
    range[9] = RF_INF; range[11] = 10;
    *index = -1;
}
void infoIacoCesare(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = false;
}
int checkIacoCesare(double *param, int timespacedim, SimulationType method) {
  if (method!=Nothing && method!=CircEmbed && method!=Direct)
      return ERRORNOTDEFINED;
  if (timespacedim == 1) return ERRORNOTDEFINED;
  if (param[KAPPAI] < 1 || param[KAPPAI] > 2) {
    strcpy(ERRORSTRING_OK,"1 <= kappa1 <= 2");
    sprintf(ERRORSTRING_WRONG, "%f",param[KAPPAI]);
    return ERRORCOVFAILED;    
  } 
  if (param[KAPPAII] < 1 || param[KAPPAII] > 2) {
    strcpy(ERRORSTRING_OK, "1 <= kappa2 <= 2");
    sprintf(ERRORSTRING_WRONG, "%f", param[KAPPAII]);
    return ERRORCOVFAILED;    
  } 
  if (param[KAPPAIII] < 0.5 * timespacedim) {
    strcpy(ERRORSTRING_OK, "kappa3 >= 0.5 * <genuine time-space dimension>");
    sprintf(ERRORSTRING_WRONG, "%f", param[KAPPAIII]);
    return ERRORCOVFAILED;    
  } 
  return NOERROR;
}



/* iaco cesare model */
double SteinST1(double *x, double *p, int dim){
    // kappa1 : nu
    // kappa2 : tau
    // kappa3-5 : z1-z3
/* 2^(1-nu) / Gamma(nu) [ h^nu K_nu(h) - 2 * tau (x T z) t h^{nu-1} K_{nu-1}(h) /
   (2 nu + d + 1) ]
*/
    register double s;
    double z, logconst;
    int d, time;   
    static double nu=RF_INF;
    static double loggamma;
    time = dim - 1;

    s = x[time] * x[time];
    z = 0.0;
    for (d=0; d<time; d++) {
	s += x[d] * x[d];
	z += x[d] * p[KAPPAIV + d];
    }
    if ( s==0.0 ) {return 1.0;}
    s = sqrt(s);

    if (nu!=p[KAPPAI]) {
	nu=p[KAPPAI];
	loggamma = lgammafn(nu);
    }
    logconst = (nu - 1.0) * log(0.5 * s)  - loggamma;
    return 
	s * exp(logconst + log(bessel_k(s, nu, 2.0)) - s)
	- 2.0 / (2 * nu + p[KAPPAIII]) * p[KAPPAII] * z * x[time] *  
	exp(logconst + log(bessel_k(s, nu - 1.0, 2.0)) - s);
}
// double ScaleSteinST1(double *p, int scaling) { return 1.0; } 
int kappasSteinST1(int dim) {return 2 + dim;}
void rangeSteinST1(int dim, int *index, double* range){
    static double range_steinST1[20]= 
	{1, RF_INF, 1, 10.0,
	 0, 1, 0, 1,
	 -RF_NEGINF, RF_INF, -100000, +100000,
	 -RF_NEGINF, RF_INF, -100000, +100000,
	 -RF_NEGINF, RF_INF, -100000, +100000,
	};
    *index = -1; 
    memcpy(range, range_steinST1, sizeof(double) * 4 * (2 + dim));
}
void infoSteinST1(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = false;
}
int checkSteinST1(double *param, int timespacedim, SimulationType method) {
  if (method!=Nothing && method!=CircEmbed && method!=Direct)
      return ERRORNOTDEFINED;
  if (timespacedim == 1) return ERRORNOTDEFINED;
  if (param[KAPPAI] < 0) {
    strcpy(ERRORSTRING_OK, "kappa1 >= 0");
    sprintf(ERRORSTRING_WRONG, "%f", param[KAPPAII]);
    return ERRORCOVFAILED;    
  } 
  if (param[KAPPAII] < 0 || param[KAPPAII] > 1) {
    strcpy(ERRORSTRING_OK, "0 <= kappa2 <= 1");
    sprintf(ERRORSTRING_WRONG, "%f", param[KAPPAII]);
    return ERRORCOVFAILED;    
  } 
  return checkWhittleMatern(param, timespacedim, method);
}



/* nugget effect model */
double nugget(double *x, double *p, int dim){
  if (*x <= NUGGET_TOL) return 1.0;  return 0.0;
}

double Scalenugget(double *p, int scaling) { return 1.0; }//or better 0.0 => error?
void rangenugget(int dim, int *index, double* range){
  *index = -1;
}
void infonugget(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = false;
}
int checknugget(double *param, int timespacedim, SimulationType method) {
  if (method!=Nothing && method!=CircEmbed && method!=Direct && method!=Nugget)
    return ERRORNOTDEFINED;
  return NOERROR;
}


void infoundefined(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 0;
  *CEbadlybehaved = true;
}
int checkundefined(double *param, int timespacedim, SimulationType method) {
    return ERRORNOTDEFINED;
}
// ---------------------------------------------------------------------



/*

void range(int dim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  double r[] = {};
  memcpy(range, r, sizeof(double) * );
}


*/
