/*
 Authors 
 Yindeng, Jiang, jiangyindeng@gmail.com
 Martin Schlather, martin.schlather@math.uni-goettingen.de

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
 Copyright (C) 2005 -- 2006 Martin Schlather

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

#define UNIT_EPSILON 1E-13

#include <math.h>
#include <assert.h>
#include "RFsimu.h"
#include "RFCovFcts.h"

//  {min, max} x {theor, pract}  x length(param) 
static double range_stable[4] = {OPEN, 2, 0.06, 2};
static double range_whittle[4]= {OPEN, RF_INF, 1e-2, 10.0};
static double range_cauchy[4] = {OPEN, RF_INF, 0.09, 10.0};
static double range_genCauchy[8] = {OPEN, 2, 0.05, 2, 
				    OPEN, RF_INF, 0.05, 10.0};
static double Besselupperbound[Nothing + 1] =
      {80, 80, 80, 80, 80, R_PosInf, 80, NA_REAL, NA_REAL, NA_REAL, NA_REAL, 
       R_PosInf};


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
/* Bessel function */ 
double Bessel(double *x,double *p){
  static double kappa=RF_INF;
  static double gamma;
  double y;
  if  (*x==0.0) {return 1.0;} 
  y = fabs( *x);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    gamma = gammafn(kappa+1.0);
  }
  return  gamma  * pow(0.5 * y,-kappa) * bessel_j(y,kappa);
}
double spectralBessel(double *p ) {
  return p[KAPPA]==0.0 ? 1.0 : 
    (sqrt(1.0 - pow(UNIFORM_RANDOM, 1.0 / p[KAPPA])));
}
int checkBessel(double *param,  int reduceddim, int dim, SimulationType method){
  // Whenever TBM3Bessel exists, add further check against too small kappa!  
  if (param[KAPPA] >= Besselupperbound[method]) {
      sprintf(ERRORSTRING_OK, "%s and kappa < %1.0f", 
	      METHODNAMES[method], Besselupperbound[method]);
      sprintf(ERRORSTRING_WRONG,"%1.3f",param[KAPPA]);
      return ERRORCOVFAILED;
  }

  if (method==SpectralTBM && param[KAPPA]<0.0) {
    strcpy(ERRORSTRING_OK,"spectral TBM and kappa>=0");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPA]);
    return ERRORCOVFAILED;
  }
  return NOERROR;
}
void rangeBessel(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract }
  *index = RANGE_LASTELEMENT; 
  range[2] = 0.0001 + (range[0] =  0.5 * ((double) reduceddim - 2.0));
  range[1] = RF_INF;
  range[3] = range[2] + 10.0;
 }
void infoBessel(double *p, int *maxdim, int *CEbadlybehaved) {
  double dim;
  dim =  (int) (2.0 * p[KAPPA] + 2.0);
  *maxdim = dim > INFDIM ? INFDIM : (int) dim;
  *CEbadlybehaved = 2;
}


/* Cauchy models */
double Cauchy(double *x, double *p){
  return pow(1.0 + *x * *x, -p[KAPPA]);
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
double TBM2Cauchy(double *x, double *p){
  register double y2, lpy2;
  y2 = *x * *x; 
  lpy2 = 1.0 + y2;
  switch ((int) (p[KAPPA] * 2.0 + 0.001)) {// ueber check sichergestellt
  case 1 : return 1.0 / lpy2;
  case 3 : return (1.0 - y2)/ (lpy2 * lpy2);
  case 5 : return (1.0-y2*(2.0+0.333333333333333333333*y2))/(lpy2*lpy2*lpy2);
  case 7 : lpy2 *= lpy2; return (1.0- y2*(3.0+y2*(1.0+0.2*y2)))/(lpy2 * lpy2);
  default : assert(false);
  }
}
double TBM3Cauchy(double *x, double *p){
  register double ha;
  ha= *x * *x;
  return (1.0 + (1.0 - 2.0 * p[KAPPA]) * ha) * pow(1.0 + ha, -p[KAPPA] - 1.0);
}
double DCauchy(double *x, double *p){
  register double y;
  y = fabs( *x);
  return (-2.0 * p[KAPPA] * y) * pow(1.0 + y * y, -p[KAPPA] - 1.0);
}
double DDCauchy(double *x, double *p){
  register double ha;
  ha = *x * *x;
  return 2.0 * p[KAPPA] * ((2.0 * p[KAPPA] + 1.0) * ha - 1.0) * 
    pow(1.0 + ha, -p[KAPPA] - 2.0);
}
int checkCauchy(double *param, int reduceddim, int dim, SimulationType method){
  if (method == TBM2) {
    // not replaced by numerical evaluation due to bad
    // numerical behaviour?!
    if (reduceddim>2) {
      strcpy(ERRORSTRING_OK, "reduced dim <= 2");
      sprintf(ERRORSTRING_WRONG, "%d", reduceddim);  
      return ERRORCOVFAILED;
    }
    if ((param[KAPPA]!=0.5) && (param[KAPPA]!=1.5) &&
	(param[KAPPA]!=2.5) && (param[KAPPA]!=3.5)) {
      strcpy(ERRORSTRING_OK,"kappa in {0.5, 1.5, 2.5 ,3.5}");
      sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
      return ERRORCOVNUMERICAL;
    }
  }
  return NOERROR;
}
void rangeCauchy(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = RANGE_LASTELEMENT; 
  memcpy(range, range_cauchy, sizeof(double) * 4);
}
void infoCauchy(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}

double Cauchytbm(double *x, double *p){
  register double ha;
  if ( *x==0) {return 1.0;}
  ha=pow(fabs( *x),p[KAPPA1]);
  return (1.0 + (1.0 - p[KAPPA2] / p[KAPPA3]) * ha) * 
      pow(1.0 + ha, -p[KAPPA2] / p[KAPPA1] - 1.0);
}
double TBM3Cauchytbm(double *x, double *p){
  register double bg,ha;
  ha=pow(fabs( *x), p[KAPPA1]);
  bg=p[KAPPA2] / p[KAPPA3];
  return  
    (1 + ha * (1-bg*(1+p[KAPPA1])+(1-p[KAPPA2]) * (1+(1-bg)*ha)))*
    pow(1+ha,-p[KAPPA2]/p[KAPPA1]-2.0);
}
double DCauchytbm(double *x, double *p){
  register double y,ha;
  if ((y = fabs(*x)) == 0.0) return 0.0; // WRONG VALUE, but multiplied 
  //                                        by zero anyway
  ha = pow(y, p[KAPPA1] - 1.0);
  return  
    p[KAPPA2] *  ha * (-1.0 - p[KAPPA1]/p[KAPPA3] + 
		       ha * y * (p[KAPPA2]/p[KAPPA3] - 1.0)) *
    pow(1.0 + ha * y,-p[KAPPA2]/p[KAPPA1]-2.0);
}
void rangeCauchytbm(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = RANGE_LASTELEMENT; 
  memcpy(range, range_genCauchy, sizeof(double) * 8);
  range[8] = range[10] = (double) reduceddim;
  range[9] = RF_INF;
  range[11] = (double) reduceddim + 10.0;
}
void infoCauchytbm(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}



/* circular model */
double circular(double *x, double *p)
{
  double y;
  if ((y=fabs(*x)) >= 1.0) return 0.0; 
  return  1.0 - (2.0 * (y * sqrt(1.0 - y * y) + asin(y))) * INVPI;
}
double Scalecircular(double *p,int scaling) {return 1.138509531721630274603;}
// spectral measure, see Lantue !! 

double Dcircular(double *x, double *p){
  register double y;
  if ((y = *x * *x) >= 1.0) {return 0.0;} 
  return -4 * INVPI * sqrt(1.0 - y);
}
void rangecircular(int reduceddim, int *index, double* range){
  *index = (reduceddim<=2) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
}
void infocircular(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim= 2;
  *CEbadlybehaved=false;
}


/* constant model */
double constant(double *x, double *p){
  return 1.0;
}
double TBM2constant(double *x, double *p) 
{
  return 1.0;
}
double TBM3constant(double *x, double *p){
   return 1.0;
}
double Dconstant(double *x, double *p){
  return 0.0;
}
void rangeconstant(int reduceddim, int *index, double* range){ 
  *index = RANGE_LASTELEMENT; 
}
void infoconstant(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM; 
  *CEbadlybehaved=false;
}

/* cone */
// no covariance function defined yet; see MPPFcts.cc for additive method
int checkcone(double *param, int reduceddim, int dim, SimulationType method) {
  switch (method) {
  case  AdditiveMpp: 
    if (reduceddim!=2) { 
      strcpy(ERRORSTRING_OK,"2 dimensional space");
      sprintf(ERRORSTRING_WRONG,"dim=%d",reduceddim);
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
  if (param[KAPPA3] + param[KAPPA2]==0.0) {
    strcpy(ERRORSTRING_OK, "kappa2+kappa3>0.0");
    sprintf(ERRORSTRING_WRONG,"c(%f,%f)", param[KAPPA2],param[KAPPA3]);
    return ERRORCOVFAILED;
  }	
  return NOERROR;
}
void rangecone(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = (reduceddim <= 3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
  double r[12] = {0, 1 - OPEN, 0, 1.0 - UNIT_EPSILON,
		0, RF_INF, UNIT_EPSILON, 10.0, 
		0, RF_INF, UNIT_EPSILON, 10.0};
  memcpy(range, r, sizeof(double) * 12);
}
void infocone(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}



/* cubic */
double cubic(double *x, double *p)
{ ///
  double y, y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return  (1.0 + (((0.75 * y2 - 3.5) * y2 + 8.75) * y - 7) * y2);
}
double TBM3cubic(double *x, double *p)
{ ///
  double y, y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return (1.0 + y2 * (-21.0 + y * (35.0 + y2 * (-21.0 + 6.0 * y2))));
}
double Dcubic(double *x, double *p) 
{ ///
  double y,y2;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;}
  y2 = y * y;
  return  y * (-14.0 + y * (26.25 + y2 * (-17.5 + 5.25 * y2)));
}
double Scalecubic(double *p,int scaling) {return 1.44855683156829;}
void rangecubic(int reduceddim, int *index, double* range){
  *index = (reduceddim<=3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
}
void infocubic(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}


/* cutoff */
// see Hypermodel.cc

/* dagum */
double dagum(double *x, double *p){
    //                     delta
    return 1 - pow((1 + pow(fabs(*x), -p[KAPPA2])), - p[KAPPA1]);
}
double Scaledagum(double *p,int scaling){ 
    return pow(pow(0.95, - 1 / p[KAPPA1]) -1, -1 / p[KAPPA2]);
} 
double Ddagum(double *x, double *p){
    register double y, xd;
    y = fabs(*x);
    xd = pow(y, -p[KAPPA2]);
    return -p[KAPPA1] * p[KAPPA2] * xd / y * pow(1 + xd, -p[KAPPA1] -1);
}
void rangedagum(int reduceddim, int *index, double* range){
    static double range_dagum[8] = {OPEN, 2, 0.01, 2, 
				    OPEN, 1-OPEN, 0.01, 0.99};
    *index = (reduceddim <= 3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
    memcpy(range, range_dagum, sizeof(double) * 8);
}
void infodagum(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved=false;
}


/*  damped cosine -- derivative of exponential:*/
double dampedcosine(double *x, double*p){
  register double y;
  y = fabs(*x);
  return exp(-y * p[KAPPA]) * cos(y);
}
double Scaledampedcosine(double *p,int scaling){ 
  if (scaling==NATSCALE_EXACT) return 0.0;
  return MINUSINVLOG005; 
} 
double TBM3dampedcosine(double *x, double *p){
  register double y;
  y = fabs(*x);
  return exp(-p[KAPPA] * y) * ((1.0 - p[KAPPA] * y) * cos(y) - y * sin(y));
}
double Ddampedcosine(double *x, double *p){
  register double y;
  y = fabs( *x);
  return - exp(-p[KAPPA]*y) * (p[KAPPA] * cos(y) + sin(y));
}
void rangedampedcosine(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = (reduceddim<=3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM; 
  range[1] = RF_INF;
  range[3] = 10.0;
  switch(reduceddim) {
      case 1 : range[2] = (range[0] = 0.0) + 1e-10; break;
      case 2 : range[2] = (range[0] = 1.0) + 1e-10; break;
      case 3 : range[2] = (range[0] = RF_M_SQRT_3) + 1e-10; break;
      default: range[0] = range[2] = RF_INF;
  }
}
void infodampedcosine(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim =
    (p[KAPPA] < 0) ? 0 : (p[KAPPA] < 1) ? 1 : (p[KAPPA] < M_SQRT_3) ? 2 : 3;
  *CEbadlybehaved=false;
}


/* exponential model */
double exponential(double *x, double *p){
  return exp(-fabs( *x));
}
double Scaleexponential(double *p, int scaling){ return MINUSINVLOG005; } 
double TBM2exponential(double *x, double *p) 
{
  double y;
  if (*x==0.0) {return 1.0;}
  y = fabs(*x);
  return 1.0 - PIHALF * y * I0mL0(y);
}
double TBM3exponential(double *x, double *p){
   register double y;
   y = fabs( *x);
   return (1.0-y)*exp(-y);
}
double Dexponential(double *x, double *p){
//  printf("\n ******* dexp %f %f\n", *x, - exp(-fabs( *x)));
  return - exp(-fabs( *x));
}
double DDexponential(double *x, double *p){
  return exp(-fabs( *x));
}
double spectralexponential(double *p) { /* see Yaglom ! */
  double register y;
  y = 1.0 - UNIFORM_RANDOM;
  return sqrt(1.0 / (y * y) - 1.0);
}
void rangeexponential(int reduceddim, int *index, double* range){ 
  *index = RANGE_LASTELEMENT; 
}
int checkexponential(double *param, int reduceddim, int dim, 
		     SimulationType method) {
  if (method==CircEmbedIntrinsic || method==CircEmbedCutoff) {
    if (reduceddim>2) 
    {
      strcpy(ERRORSTRING_OK,"genuine dim<=2");
      sprintf(ERRORSTRING_WRONG,"%d",reduceddim);
      return ERRORCOVFAILED;
    }
  }
  return NOERROR;
}
int hyperexponential(double radius, double *center, double *rx,
		     int reduceddim, bool simulate, 
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
  
  if (reduceddim==2) {
    // we should be in two dimensions
    // first, we simulate the lines for a rectangle with center (0,0)
    // and half the side length equal to lenx
    lx = rx[0];
    ly = rx[1];
    lambda = TWOPI * radius * 0.5; /* total, integrated, intensity */
    //    0.5 in order to get scale 1
    if (!simulate) return lambda < 999999 ? (int) lambda : 999999 ;
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
    // reduceddim = 1  --  not programmed yet
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



// Brownian motion 
double fractalBrownian(double*x, double *p) {
  return - pow(fabs(*x), p[KAPPA]);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
//begin
/* fractalBrownian: first derivative at t=1 */
double DfractalBrownian(double *x, double*p) 
{// FALSE VALUE FOR *x==0 and  p[KAPPA] < 1
    return (*x == 0.0) ? 0.0 : -p[KAPPA] * pow(fabs(*x), p[KAPPA] - 1.0);
}
/* fractalBrownian: second derivative at t=1 */
double DDfractalBrownian(double *x, double*p)  
{// FALSE VALUE FOR *x==0 and  p[KAPPA] < 2
    return (*x == 0.0) ? 0.0 : 
	-p[KAPPA] * (p[KAPPA] - 1.0) * pow(fabs(*x), p[KAPPA] - 2.0);
}
void rangefractalBrownian(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  double r[4] = {OPEN, 2.0, UNIT_EPSILON, 2.0 - UNIT_EPSILON};
  *index = (reduceddim <= 3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
  memcpy(range, r, sizeof(double*) * 4);
}
void infofractalBrownian(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = false;
}



/* FD model */
double FD(double *x,double *p){
  static double dold=RF_INF;
  static double kold, sk;
  double y, d, k, skP1;
  d = p[KAPPA1] * 0.5;
  y = fabs(*x);
  k = (double) trunc((double) y);
  if (dold!=d || kold > k) {
    sk = 1;
    kold = 0.0;
  }
  // sign (-1)^k is (kold+d), 16.11.03, checked. 
  for (; kold<k; kold += 1.0) sk =  sk * (kold + d) / (kold + 1.0 - d);
  dold = d;
  kold = k;
  if (k == y) return sk; 
  skP1 = sk * (kold + d) / (kold + 1.0 - d);
  return sk + (y - k) * (skP1 - sk);
}
void rangeFD(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = (reduceddim==1) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
  range[2] = (range[0] = -1.0) + UNIT_EPSILON;
  range[1] = 1.0 - OPEN;
  range[4] = 1.0 - UNIT_EPSILON;
}
void infoFD(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 1;
  *CEbadlybehaved = false;
}


/* fractgauss */
double fractGauss(double *x, double *p){
  register double y;
  if ((y = fabs(*x)) == 0.0) return 1.0;
  return 0.5 * (pow(fabs(y + 1.0), p[KAPPA])  
		- 2.0 * pow(y, p[KAPPA]) 
		+ pow(fabs(y - 1.0), p[KAPPA]) 
		);
}
void rangefractGauss(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  double r[4] = {OPEN, 2, UNIT_EPSILON, 2};
  *index = (reduceddim==1) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
  memcpy(range, r, sizeof(double) * 4);
}
void infofractGauss(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 1;
  *CEbadlybehaved = false;
}


/* Gausian model */
double Gauss(double *x, double*p) {
  return  exp(- *x * *x);
}
double ScaleGauss(double *p,int scaling) {return SQRTINVLOG005;}
double TBM3Gauss(double *x, double*p) {
  register double y;
  y = *x * *x; 
  return (1 - 2.0 * y) * exp(-y);
}
double DGauss(double *x, double*p) {
  register double y;
  y = fabs(*x); 
  return -2.0 * y * exp(- y * y);
}
double spectralGauss(double *p ) {   
  return 2.0 * sqrt(-log(1.0 - UNIFORM_RANDOM));
}
void rangeGauss(int reduceddim, int *index, double* range){
  *index = RANGE_LASTELEMENT;
}
void infoGauss(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}


/* gencauchy */
double generalisedCauchy(double *x, double *p){
  return pow(1.0 + pow(fabs(*x), p[KAPPA1]), -p[KAPPA2]/p[KAPPA1]);
}
double ScalegeneralisedCauchy(double *p,int scaling) {
  switch(scaling) {
  case NATSCALE_EXACT: case NATSCALE_APPROX:
    return pow(pow(0.05, -p[KAPPA1] / p[KAPPA2]) - 1.0, -1.0 / p[KAPPA1]); 
    break;
  case NATSCALE_MLE: 
    // should be changed! (long tails!)
    return pow(pow(0.05, -p[KAPPA1] / p[KAPPA2]) - 1.0, -1.0 / p[KAPPA1]);
    break;
  default: assert(false);
  }
}
double TBM3generalisedCauchy(double *x, double *p){
  register double ha;
  ha=pow(fabs( *x), p[KAPPA1]);
  return 
    (1.0 + (1.0 - p[KAPPA2]) * ha) * pow(1.0 + ha, -p[KAPPA2] / p[KAPPA1] - 1.0);
}
double DgeneralisedCauchy(double *x, double *p){
  register double ha,y;
  if ((y = fabs(*x))==0.0) 
    return ((p[KAPPA1] > 1.0) ? 0.0 : (p[KAPPA1] < 1.0) ? -INFTY : -p[KAPPA2]); 
  ha=pow(y, p[KAPPA1] - 1.0);
  return  -p[KAPPA2] * ha * pow(1.0 + ha * y, -p[KAPPA2] / p[KAPPA1] - 1.0);
}
double DDgeneralisedCauchy(double *x, double *p){
  register double ha,y;
  if ((y = fabs(*x))==0.0) 
    return ((p[KAPPA1]==2.0) ? p[KAPPA2] * (p[KAPPA2] + 1.0) : INFTY); 
  ha=pow(y, p[KAPPA1]);
  return p[KAPPA2] * ha / (y * y) * (1.0 - p[KAPPA1] + (1.0 + p[KAPPA2]) * ha)
	    * pow(1.0 + ha, -p[KAPPA2] / p[KAPPA1] - 2.0);
}
int checkgeneralisedCauchy(double *param, int reduceddim, int dim, 
			   SimulationType method){
  if (method==CircEmbedIntrinsic || method==CircEmbedCutoff)
  {
    if (reduceddim>2) 
    {
      strcpy(ERRORSTRING_OK,"genuine total dim<=2, currently");
      sprintf(ERRORSTRING_WRONG,"%d",reduceddim);
      return ERRORCOVFAILED;
    }
  }
  return NOERROR;
}
void rangegeneralisedCauchy(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = RANGE_LASTELEMENT; 
  memcpy(range, range_genCauchy, sizeof(double) * 8);
}
void infogeneralisedCauchy(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = 2;
}


/* gengneiting */
double genGneiting(double *x, double *p)
{
  register double y, s;
  if ((y=fabs( *x)) >= 1.0) return 0.0; 
  s = p[KAPPA2] + p[KAPPA1];
  switch ((int) p[KAPPA1]) {
  case 1:
    return  pow ((1.0 - y), s) * (1.0 + s*y);   
  case 2:
    return  pow ((1.0 - y), s) * 
      (1.0 + y * (s + y * (s*s-1.0)*0.3333333333333333));  
  case 3:
    return   pow ((1.0 - y), s) * 
      (1.0 + y * (s 
		  + y * ( 0.2 * (2.0*s*s - 3.0)
			 + y * (s*s-4.0)*s*0.0666666666666666666))) ;
  default : assert(false);  
  }
}
double TBM3genGneiting(double *x, double *p)
{
  register double y, s;
  if ((y=fabs( *x)) >= 1.0) {return 0.0;} 
  s = p[KAPPA2] + p[KAPPA1];
  switch ((int) p[KAPPA1]) {
  case 1:
    return  pow((1.0 - y), p[KAPPA2]) *
      (1.0 + y * (p[KAPPA2] - y * s * (p[KAPPA2] + 3.0))); 
  case 2:
   return  pow((1.0 - y), (p[KAPPA2]+1.0)) *
      (1.0 + y * ((p[KAPPA2]+1.0)
		- y *(1.0 + 2.0 * s 
		      + y * (s * s - 1.0) * 0.33333333333333 * 
		      (p[KAPPA2] + 5.0))));
  case 3:
    return  pow((1.0 - y), (p[KAPPA2]+2.0)) *
      (1.0 + y * (s - 1.0
		  + y * (0.2 * (-9.0 + s * (-10.0 + s))
			 + y * 0.0666666666666666 * 
			 ((27.0 - s*(7.0 + s * (18.0 + 2.0*s)))
			  - y * (s*s - 4.0) * s * (s + 4.0)))));
  default : assert(false);   
  }
}
double DgenGneiting(double *x, double *p)
{
  register double y, s;
  if ((y=fabs(*x)) >= 1.0) {return 0.0;} 
  s = p[KAPPA2] + p[KAPPA1];

  switch ((int) p[KAPPA1]) {
  case 1:
    return  - pow(1.0 - y, s - 1.0) * y * s  * (s + 1.0); 
  case 2:
    return - pow(1.0 - y, s - 1.0) *
      y * (0.333333333333333333 * s * s + 0.66666666666666666667 + s +
	   y * 0.333333333333333333 * (s * s - 1.0) * (p[KAPPA2] + 4.0) );
  case 3:
    return -  pow(1.0 - y, s - 1.0) *
      0.2 * y * (6 + s * (5.0 + s) +
		 y * (-6 + s * (1.0 + s * (4.0 + s)) +
		      0.33333333333333333 * y * s * (-12.0 + s * 
						       (-4.0 + s * (3.0 + s)))
		   ));
  default : assert(false);   
  }
}
void rangegenGneiting(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  range[0] = range[2] = range[1] = range[3] = *index;
  range[4] = range[6] = 0.5 * (double) (reduceddim + 2 * *index + 1);
  range[5] = RF_INF;
  range[7] = 10.0 + range[6];
  if ((++(*index)) > 3) *index = RANGE_LASTELEMENT;
}
void infogenGneiting(double *p, int *maxdim, int *CEbadlybehaved) {
  double dim;
  dim = (2.0 * p[KAPPA2] - 1.0 - 2.0 * p[KAPPA1]);  
  *maxdim = dim > INFDIM ? INFDIM : (int) dim;
  *CEbadlybehaved = false;
}



/* Gneiting's functions -- alternative to Gaussian */
// #define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
#define NumericalScale 0.301187465825
double Gneiting(double *x, double *p){ 
  register double y,oneMy8;
  if ((y=fabs(*x * NumericalScale)) >= 1.0) {return(0.0);}
  oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
  return ((1.0+y * ( 8.0 + y * (25.0 + 32.0 *y)))*oneMy8);
}
double ScaleGneiting(double *p,int scaling) {return 0.5854160193;}
double TBM3Gneiting(double *x, double *p){ 
  register double y,oneMy7;
  if ((y=fabs( *x * NumericalScale)) >= 1.0) {return 0.0;}  
  oneMy7 = 1.0-y; oneMy7*=oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);
  return 
    (1.0 + y * (7.0  -  y * (5.0 + y * (147.0 + 384.0 * y))))* oneMy7;
}
double DGneiting(double *x, double *p){ 
  register double y,oneMy7;
  if ((y=fabs( *x * NumericalScale)) >= 1.0) {return 0.0;}  
  oneMy7 = 1.0-y; oneMy7 *= oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);

/*  double zz,z;
  z= (-y) * ( 22.0 + y * (154.0 + y * 352.0)) * oneMy7;
  p[KAPPA1] = 3; p[KAPPA2] =5.0;
  zz=  DgenGneiting(&y, p);
  printf("%f %f\n",z, zz);
  assert(z==zz);
  assert(false);
*/

  return 
    (-y) * ( 22.0 + y * (154.0 + y * 352.0)) * oneMy7 * NumericalScale;
}
void rangeGneiting(int reduceddim, int *index, double* range){
  *index = (reduceddim<=3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
}
void infoGneiting(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}


/* hyperbolic */
double hyperbolic(double *x, double*p){ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double kappay;
  double y;
  if (*x==0.0) {return 1.0;}
  if (p[KAPPA3]==0) { // whittle matern
    y = *x * p[KAPPA1];
    return WhittleMatern(&y, p);
  } 
  if (p[KAPPA1]==0) { //cauchy   => KAPPA2 < 0 !
    y =  *x / p[KAPPA3];
    /* note change in sign as KAPPA2<0 */
    return pow(1 + y * y, p[KAPPA2]); 
  }  

  y=fabs(*x); 
  if ((p[KAPPA1]!=kappa) || (p[KAPPA2]!=lambda) || (p[KAPPA3]!=delta)) {
    kappa = p[KAPPA1]; 
    lambda= p[KAPPA2];
    delta = p[KAPPA3];
    deltasq = delta * delta;
    kappadelta = kappa * delta;
    logconst = kappadelta - log(bessel_k(kappadelta, lambda, 2.0)) 
      - lambda * log(delta);
  }
  y=sqrt(deltasq + y * y);  
  kappay = kappa * y;
  return exp(logconst + lambda * log(y) + log(bessel_k(kappay, lambda, 2.0)) 
	     - kappay);
}
double TBM3hyperbolic(double *x, double*p)
{ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double y;
  double ysq,s,kappa_s,logs;
  if (*x==0.0) {return 1.0;}
  if (p[KAPPA3]==0) { // whittle matern
    y = *x * p[KAPPA1];
    return TBM3WhittleMatern(&y, p);
  } 
  if (p[KAPPA1]==0) { //cauchy
    register double y,ha;
    y= *x / p[KAPPA3];
    ha=y * y;
    /* note change in sign as KAPPA2<0 */
    return (1.0 + (1.0 + 2.0 * p[KAPPA2]) * ha) * 
      pow(1.0 + ha, p[KAPPA2] - 1.0);
  }
  y=fabs(*x); 
  if ((p[KAPPA1]!=kappa) || (p[KAPPA2]!=lambda) || (p[KAPPA3]!=delta)) {
    kappa = p[KAPPA1]; 
    lambda= p[KAPPA2];
    delta = p[KAPPA3];
    deltasq = delta * delta;
    kappadelta = kappa * delta;
    logconst = kappadelta - log(bessel_k(kappadelta, lambda, 2.0)) 
      - lambda * log(delta);
  }
  ysq = y * y;
  s=sqrt(deltasq + ysq);
  kappa_s = kappa * s;
  logs = log(s);  
  return  
    ( exp(logconst + lambda * logs +log(bessel_k(kappa_s,lambda,2.0))-kappa_s)
      - ysq*kappa*exp(logconst + (lambda-1.0)*logs 
		      +log(bessel_k(kappa_s, lambda - 1.0, 2.0)) - kappa_s)
      );
}
double Dhyperbolic(double *x, double*p)
{ 
  static double kappa = RF_INF;
  static double lambda= RF_INF;
  static double delta = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double y;
  double s,kappa_s,logs;
  if (*x==0.0) {return 1.0;}
  if (p[KAPPA3]==0) { // whittle matern
    y = *x * p[KAPPA1];
    return DWhittleMatern(&y, p);
  } 
  if (p[KAPPA1]==0) { //cauchy
    register double y,ha;
    y= *x / p[KAPPA3];
    ha=y * y;
    /* note change in sign as KAPPA2<0 */
    return  2.0 * p[KAPPA2] * fabs(y) * 
      pow(1.0 + ha, p[KAPPA2] - 1.0);
  }
  y=fabs( *x); 
  if ((p[KAPPA1]!=kappa) || (p[KAPPA2]!=lambda) || (p[KAPPA3]!=delta)) {
    kappa = p[KAPPA1]; 
    lambda= p[KAPPA2];
    delta = p[KAPPA3];
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
      - y * kappa*exp(logconst + (lambda-1.0) * logs 
		      +log(bessel_k(kappa_s,lambda-1.0,2.0))-kappa_s)
      );
}
int checkhyperbolic(double *param, int reduceddim, int dim, 
		    SimulationType method){
  if (param[KAPPA2] >= Besselupperbound[method]) {
      sprintf(ERRORSTRING_OK, "%s and kappa2 < %1.0f", 
	      METHODNAMES[method], Besselupperbound[method]);
      sprintf(ERRORSTRING_WRONG,"%1.3f",param[KAPPA2]);
      return ERRORCOVFAILED;
  }

  if (param[KAPPA2]>0) {
    if ((param[KAPPA3]<0) || (param[KAPPA1]<=0)) {
      strcpy(ERRORSTRING_OK,"kappa1>0 and kappa3>=0 if kappa2>0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      param[KAPPA1],param[KAPPA3],param[KAPPA2]);
      return ERRORCOVFAILED;
    }
  } else if (param[KAPPA2]<0) {
    if ((param[KAPPA3]<=0) || (param[KAPPA1]<0)) {
      strcpy(ERRORSTRING_OK,"kappa1>=0 and kappa3>0 if kappa2<0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      param[KAPPA1],param[KAPPA3],param[KAPPA2]);
      return ERRORCOVFAILED;
    }
  } else { // param[KAPPA2]==0.0
    if ((param[KAPPA3]<=0) || (param[KAPPA1]<=0)) {
      strcpy(ERRORSTRING_OK,"kappa1>0 and kappa3>0 if kappa2=0");
      sprintf(ERRORSTRING_WRONG,"kappa1=%f and kappa3=%f for kappa2=%f",
	      param[KAPPA1],param[KAPPA3],param[KAPPA2]);
      return ERRORCOVFAILED;
    }
  }
  return NOERROR;
}
void rangehyperbolic(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  double r[12] = {0, RF_INF, 0.000001, 10.0,
		  RF_NEGINF, RF_INF, -20.0, 20.0,
		  0, RF_INF, 0.000001, 10.0};
  *index = RANGE_LASTELEMENT; 
  memcpy(range, r, sizeof(double) * 12);
}
void infohyperbolic(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = true;
}



/* iaco cesare model */
double IacoCesare(double *x, double *p){
    return pow(1.0 + pow(x[0], p[KAPPA1]) + pow(fabs(x[1]), p[KAPPA2]),
	       - p[KAPPA3]); 
}
//double ScaleIacoCesare(double *p, int scaling) { return 1.0; } 
void rangeIacoCesare(int reduceddim, int *index, double* range){
    static double range_iacocesare[8]= 
	{1, 2, 1, 2,
	 1, 2, 1, 2};
    memcpy(range, range_iacocesare, sizeof(double) * 8);
    range[8] = range[10] = 0.5 * reduceddim;
    range[9] = RF_INF; range[11] = range[10] + 10.0;
    *index = RANGE_LASTELEMENT;
}
void infoIacoCesare(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = false;
}
int checkIacoCesare(double *param, int reduceddim, int dim, 
		    SimulationType method) {
  if (method!=Nothing && method!=CircEmbed && method!=Direct)
      return ERRORNOTDEFINED;
  if (reduceddim == 1) return ERRORWRONGDIM;
  return NOERROR;
}


/* local-global distinguisher */
double lgd1(double *x, double*p) {
  double y;
  if (*x == 0.0) return 1.0;
  if ((y=fabs(*x)) < 1.0) return 1.0 - p[KAPPA2] / (p[KAPPA1] + p[KAPPA2])
	     * exp(p[KAPPA1] * log(y));
  else return p[KAPPA1] / (p[KAPPA1] + p[KAPPA2]) 
	 * exp( -p[KAPPA2] * log(y));
}
double Scalelgd1(double *p,int scaling) {
  return (19 * p[KAPPA1] < p[KAPPA2])
     ? exp(log(0.95 * (p[KAPPA1] + p[KAPPA2]) / p[KAPPA2]) / p[KAPPA1])
     : exp(log(0.05 * (p[KAPPA1] + p[KAPPA2]) / p[KAPPA1]) / p[KAPPA2]);
}
double Dlgd1(double *x, double *p){
  double y, pp;
  if ( (y=fabs(*x)) == 0.0) return 0.0;// falscher Wert, aber sonst NAN-Fehler
  pp = ( (y < 1.0) ? p[KAPPA1] : -p[KAPPA2] ) - 1.0;
  return - p[KAPPA1] * p[KAPPA2] / (p[KAPPA1] + p[KAPPA2]) * exp(pp * y);
}
void rangelgd1(int reduceddim, int *index, double* range){
  static double range_lgd1[8] =
    {OPEN, 1.0, 0.01, 1.0, 
     OPEN, RF_INF, 0.01, 20.0};
  //  2 x length(param) x {theor, pract } 
  *index = (reduceddim<=2) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
  memcpy(range, range_lgd1, sizeof(double) * 8);
  if (reduceddim==2) range[1] = range[3] = 0.5;
}
void infolgd1(double *p, int *maxdim, int *CEbadlybehaved) {
  double dim;
  dim = 2.0 * (1.5 - p[KAPPA]);
  *maxdim = dim > INFDIM ? INFDIM : (int) dim;
  *CEbadlybehaved = true;
}


/* mastein */
// see Hypermodel.cc


/* nsst */
/* Tilmann Gneiting's space time models, part I */
double InvSqrtPsi(double x, double a, double b, int c) {
  double y;
  if (x == 0.0) return 1.0;
  y = pow(fabs(x), a);
  switch(c) {
  case 1 : 
    return pow(1.0 + y, - 0.5 * b);
  case 2 : 
    return sqrt( (y + 1.0) / (y / b + 1.0) );
  case 3 :
    return sqrt(-log(b) / log(y + 1.0 / b));
  default: assert(false);
  }
}
double spacetime1(double *x, double *p){
  double y, z, invsqrtpsi;
  assert(p[EFFECTIVEDIM] == 2.0);
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPA3], p[KAPPA4], (int) p[KAPPA5]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPA2]) {
   case 1 : z = stable(&y, p); break;  // called with the wrong EFFECTIVEDIM,
       // but not checked by stable, neither by WhittleMatern nor Cauchy
  case 2 : z = WhittleMatern(&y, p); break;
  case 3 : z = Cauchy(&y, p); break;
  default: assert(false);
  }
//  printf("%f %f %f %f %f\n", x[0], x[1], invsqrtpsi, z, pow(invsqrtpsi, p[KAPPA6]) * z);
  return pow(invsqrtpsi, p[KAPPA6]) * z;
}
double TBM2spacetime1(double *x, double *p){
  double y, z, invsqrtpsi;
  invsqrtpsi = InvSqrtPsi(x[1], p[KAPPA3], p[KAPPA4], (int) p[KAPPA5]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPA2]) {
  case 1 : assert(false); break;
  case 2 : z = TBM2WhittleMatern(&y, p); break;
  case 3 : z = TBM2Cauchy(&y, p); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPA6]) * z;
}
double Dspacetime1(double *x, double *p){
  double y, z, invsqrtpsi;
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPA3], p[KAPPA4], (int) p[KAPPA5]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPA2]) {
  case 1 : z = Dstable(&y, p); break;
  case 2 : z = DWhittleMatern(&y, p); break;
  case 3 : z = DCauchy(&y, p); break;
  default: assert(false);
  }
  //printf("%f %f %f %f %f \n", x[0], x[1], invsqrtpsi, z,
//	 pow(invsqrtpsi, p[KAPPA6] + 1.0) * z);
  return pow(invsqrtpsi, p[KAPPA6] + 1.0) * z;
}
int checkspacetime1(double *param, int reduceddim, int dim, 
		    SimulationType method) {
  int error;
  // 1 : stable
  // 2 : whittle
  // 3 : cauchy (not generalised)
  // first parameter(s) for phi
  // then choice of phi; then two parameters for psi, then choice of psi
  error = NOERROR;
  switch((int) param[KAPPA2]) {
  case 1 : 
    error = checkstable(param, reduceddim, dim, method);
    if (error==NOERROR && method==TBM2)  {
      strcpy(ERRORSTRING_OK, "kappa2=2,3");
      sprintf(ERRORSTRING_WRONG,"%f", param[KAPPA2]);
      error = ERRORCOVNUMERICAL;
    }
    break;
  case 2 : 
    error = checkWhittleMatern(param, reduceddim, dim, method);
    break;
  case 3 : error = checkCauchy(param, reduceddim, dim, method); 
    break;
  default : 
      assert(false);
//    strcpy(ERRORSTRING_OK,"kappa2=1,2,3");
//    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPA2]);
//    return ERRORCOVFAILED;
  }
  if (param[KAPPA4]==1 && param[KAPPA5]==3) {
    strcpy(ERRORSTRING_OK,"kappa4<1 if kappa5=3");
    sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA4]);
    return ERRORCOVFAILED;
  } 
  return error;
}
void rangespacetime1(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract }
  double *r;
  if ((*index <= 0) || (*index > 9)) { // see also last line
    int i; 
    for (i=0; i<24; i++) range[i]= RF_NAN; 
    *index = RANGE_LASTELEMENT; 
    return;
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

  range[8] = OPEN;
  range[9] = range[11] = 2;
  range[10] = UNIT_EPSILON;

  range[12] = OPEN;
  range[13] = 1 - OPEN;
  range[14] = UNIT_EPSILON;
  range[15] = 1.0 - UNIT_EPSILON;
  
  range[16] = range[17] = 
    range[18] = range[19] =  (double) (1 + (*index - 1) / 3);
  
  range[20] = range[22] = (double) reduceddim - 1; // spatial reduceddim
  range[21] = RF_INF;
  range[23] = range[22] + 10.0;

  if ( (++(*index)) > 9) *index = RANGE_LASTELEMENT;
}
void infospacetime1(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = p[KAPPA2]==3 || (p[KAPPA2]==2 && p[KAPPA1]>1.5);
}


/* Tilmann Gneiting's space time models, part II*/
double spacetime2(double *x, double *p){
  double y, z, invsqrtpsi;
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPA4], p[KAPPA5], (int)p[KAPPA6]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPA3]) {
  case 1 : z = generalisedCauchy(&y, p); break;
       // called with the wrong EFFECTIVEDIM,
       // but not checked by generalisedCauchy
  default: assert(false);
  }
//  if (fabs(x[0]) < 1e-10) invsqrtpsi=0.;
//  printf("%e %e %e %e %e %d, %d %d %d %d\n", 
//	 x[1], invsqrtpsi, p[KAPPA7], pow(invsqrtpsi, p[KAPPA7]),
//	 R_pow(invsqrtpsi, p[KAPPA7]), 
//	 x[1] == 0.,
//	 invsqrtpsi==1., invsqrtpsi < 1.,invsqrtpsi > 1., p[KAPPA7]>0.);
  return pow(invsqrtpsi, p[KAPPA7]) * z;
}
//double TBM3spacetime2(double *x, double *p){
//  assert(false); // should never happen
//  return spacetime2(x, p, 2) + Dspacetime2(x, p) * fabs(x[0]);
//}
double Dspacetime2(double *x, double *p){
  double y, z, invsqrtpsi;
  invsqrtpsi = InvSqrtPsi(x[1],  p[KAPPA4], p[KAPPA5], (int)p[KAPPA6]);
  y = x[0] * invsqrtpsi; 
  switch((int) p[KAPPA3]) {
  case 1 : z = DgeneralisedCauchy(&y, p); break;
  default: assert(false);
  }
  return pow(invsqrtpsi, p[KAPPA5] + 1.0) * z;
}
int checkspacetime2(double *param, int reduceddim, int dim, 
		    SimulationType method) {
  int error;
  // 1 : generalisedcauchy 
  // first parameter(s) for phi
  // then choice of phi; then two parameters for psi, then choice of psi
  error = NOERROR;
  switch((int) param[KAPPA3]) {
  case 1 : error = checkgeneralisedCauchy(param, reduceddim, dim, method); break;
  default : 
    assert(false);
//    strcpy(ERRORSTRING_OK,"kappa3=1");
//    sprintf(ERRORSTRING_WRONG,"%d",(int) param[KAPPA3]);
//    return ERRORCOVFAILED;
  }
  if (param[KAPPA5]==1 && param[KAPPA6]==3) {
    strcpy(ERRORSTRING_OK, "kappa5<1 if kappa6=3");
    sprintf(ERRORSTRING_WRONG, "%f", param[KAPPA5]);
    return ERRORCOVFAILED;
  }
  return error;
}
void rangespacetime2(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract }
  double *r;
  if ((*index<=0) || (*index>3)) { // see also last line
    int i; 
    for (i=0; i<28; i++) range[i] = RF_NAN; 
    *index = RANGE_LASTELEMENT;
    return;
  }
  switch ((*index-1) % 1) {
  case 0 : r = range_genCauchy;  break;
  default: assert(false);
  }
  memcpy(range, r, sizeof(double) * 8);

  range[8] = range[10] =
    range[9] = range[11] = (double) (1 + (*index - 1) % 1);

  range[12] = OPEN;
  range[13] = range[15] = 2.0;
  range[14] = UNIT_EPSILON;

  range[16] = OPEN;
  range[17] = 1.0;
  range[18] = UNIT_EPSILON;
  range[19] = 1.0 - UNIT_EPSILON;

  range[20] = range[21] = 
    range[22] = range[23] =  (double) (1 + (*index - 1) / 1);

  range[24] = range[26] = (double) (reduceddim - 1); // spatial dim
  range[25] = RF_INF;
  range[27] = range[26] + 10.0;
  
  if ( (++(*index)) > 3) *index = RANGE_LASTELEMENT;
}
void infospacetime2(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = true;
}


/* nugget effect model */
double nugget(double *x, double *p){
  if (*x <= NUGGET_TOL) return 1.0; 
  return 0.0;
}

double Scalenugget(double *p, int scaling) { return 1.0; }//or better 0.0 => error?
void rangenugget(int reduceddim, int *index, double* range){
  *index = RANGE_LASTELEMENT;
}
void infonugget(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = false;
}
int checknugget(double *param, int reduceddim, int dim, SimulationType method) {
  if (method!=Nothing && method!=CircEmbed && method!=Direct && method!=Nugget)
    return ERRORNOTDEFINED;
  return NOERROR;
}



/* penta */
double penta(double *x, double *p)
{ ///
  double y, y2;
  if ((y = fabs(*x)) >= 1.0) return 0.0;
  y2 = y * y;
  return  
    (1.0 + y2 * (-7.333333333333333 
		 + y2 * (33.0 +
			 y * (-38.5 + 
			      y2 * (16.5 + 
				    y2 * (-5.5 + 
					  y2 * 0.833333333333333))))));
}
double TBM3penta(double *x, double *p)
{ ///
  double y, y2;
  if ( (y = fabs(*x)) >= 1.0 ) return 0.0;
  y2 = y * y;
  return 
    (1.0 + y2 * (-22.0 
		 + y2 * (165.0 
			 + y * (-231.0 
				+ y2 * (132.0 
					+ y2 * (-55.0 
						+ 10.0 * y2))))));
}
double Dpenta(double *x, double *p)
{ ///
  double y,y2;
  if ( (y = fabs(*x)) >= 1.0 ) return 0.0;
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
void rangepenta(int reduceddim, int *index, double* range){
  *index = (reduceddim<=3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
}
void infopenta(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = false;
}


/* power model */ 
double power(double *x, double *p){
  register double y;
  if ((y=fabs(*x)) >= 1.0) return 0.0;
  return pow(1.0 - y, p[KAPPA]);
}
double Scalepower(double *p,int scaling){ 
  return 1.0 / (1.0 - pow(0.05, 1.0 / p[KAPPA]));
}
double TBM2power(double *x, double *p){
  // only kappa=2 up to now !
  register double y;
  y = fabs(*x);
  return (y > 1.0) 
      ? (1.0 - 2.0 * y *(asin(1.0 / y) - y + sqrt(y * y - 1.0) ))
      : 1.0 - y * (PI - 2.0 * y);
}
double TBM3power(double *x, double *p){
  register double y;
  if ( (y=fabs(*x)) >= 1.0 ) return 0.0;
  return (1.0 - y - y * p[KAPPA]) * pow(1.0 - y, p[KAPPA] - 1.0);
}
double Dpower(double *x, double *p){
  register double y;
  if ( (y=fabs(*x)) >= 1.0 ) return 0.0;
  return  - p[KAPPA] * pow(1.0 - y, p[KAPPA] - 1.0);
}
int checkpower(double *param, int reduceddim, int dim, SimulationType method) {
  if (method == TBM2 && param[KAPPA]!=2.0) {
    strcpy(ERRORSTRING_OK, "kappa=2");
    sprintf(ERRORSTRING_WRONG,"%f", param[KAPPA]);
    return ERRORCOVNUMERICAL;
  }
  return NOERROR;
}
// range definition:
// 0: min, theory, 1:max, theory
// 2: min, practically 3:max, practically
void rangepower(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = RANGE_LASTELEMENT; 
  range[0] = range[2] = 0.5 * (double) (reduceddim + 1);
  range[1] = RF_INF;
  range[3] = 20.0;
}
void infopower(double *p, int *maxdim, int *CEbadlybehaved) {
  double dim;
  dim = 2.0 * p[KAPPA] - 1.0;
  *maxdim = dim > INFDIM ? INFDIM : (int) dim;
  *CEbadlybehaved=false;
}


/* qexponential -- derivative of exponential */
double qexponential(double *x,double *p){
  register double y;
  y = exp(-fabs(*x));
  return y * (2.0  - p[KAPPA] * y) / (2.0 - p[KAPPA]);
}
double Scaleqexponential(double *p,int scaling){
  return -1.0 / 
    log( (1.0 - sqrt(1.0 - p[KAPPA] * (2.0 - p[KAPPA]) * 0.05)) / p[KAPPA]);
} 
double TBM3qexponential(double *x, double *p) {
   register double y;
   y = exp(-fabs( *x));
   return y * ((2.0  - p[KAPPA] * y) + fabs( *x) * (p[KAPPA] * y - 1.0) * 2.0) /
     (2.0 - p[KAPPA]);
}
double Dqexponential(double *x, double *p) {
  register double y;
  y = exp(-fabs( *x));
  return y * (p[KAPPA] * y - 1.0) * 2.0 / (2.0 - p[KAPPA]);
}
void rangeqexponential(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = RANGE_LASTELEMENT; 
  static double r[4] = {0.0, 1.0, 0.0, 1.0};
  memcpy(range, r, sizeof(double) * 4);
}
void infoqexponential(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = (p[KAPPA]>=0.0 && p[KAPPA]<=1.0) ? INFDIM : 0; 
  *CEbadlybehaved=false;
}



/* spherical model */ 
double spherical(double *x, double *p){
  register double y;
  if ((y=fabs( *x)) >= 1.0) return 0.0;
  return (1.0 + y * 0.5 * (y * y - 3.0));
}
double Scalespherical(double *p,int scaling){ return 1.23243208931941;}
double TBM2spherical(double *x, double *p){
  register double y, y2;
  y  = fabs(*x);
  y2 = y * y;
  return (y>1.0) 
      ? (1.0- 0.75 * y * ((2.0 - y2) * asin(1.0/y) + sqrt(y2 -1.0)))
      : (1.0 - 0.375 * PI * y * (2.0 - y2));
}
double TBM3spherical(double *x, double *p){
  register double y;
  if ((y=fabs(*x)) >= 1.0) return 0.0; 
  return (1.0 + (-3.0 + 2.0 * y * y) * y);
}
double Dspherical(double *x, double *p){
  register double y;
  if ((y=fabs(*x)) >= 1.0) return 0.0;
  return 1.5 * (y * y - 1.0);
}
void rangespherical(int reduceddim, int *index, double* range){
  *index = (reduceddim<=3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
}
void infospherical(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved=false;
}



/* stable model */
double stable(double *x,double *p){
  return (*x==0.0) ? 1.0 : exp(-pow(fabs(*x), p[KAPPA]));
}
double Scalestable(double *p,int scaling){return pow(MINUSINVLOG005,1/p[KAPPA]);}
double TBM3stable(double *x, double *p){
  register double y;
  y = pow(fabs(*x), p[KAPPA]);
  return exp(-y) * (1.0 - p[KAPPA] * y);
}
double Dstable(double *x, double *p){
  register double y, z;
  if ( (z = fabs( *x)) == 0.0) 
    return ((p[KAPPA] > 1.0) ? 0.0 : (p[KAPPA] < 1.0) ? INFTY : 1.0);
  y = pow(z, p[KAPPA] - 1.0);
  return -p[KAPPA] * y * exp(-y * z);
}
/* stable: second derivative at t=1 */
double DDstable(double *x, double*p) 
{
  double y, xkappa, z;
  if ( (z = fabs( *x)) == 0.0) return ((p[KAPPA] != 1.0) ? INFTY : 1.0);
  y = pow(z, p[KAPPA] - 2.0);
  xkappa = y *z * z;
  return p[KAPPA] * (1.0 - p[KAPPA] + p[KAPPA] * xkappa) * y * exp(-xkappa);
}

int checkstable(double *param, int reduceddim, int dim, SimulationType method) {
  if (method==CircEmbedIntrinsic || method==CircEmbedCutoff) {
    if (reduceddim>2) 
    {
      strcpy(ERRORSTRING_OK,"genuine total dim<=2, currently");
      sprintf(ERRORSTRING_WRONG,"%d",reduceddim);
      return ERRORCOVFAILED;
    }
  }
  return NOERROR;
}
void rangestable(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = RANGE_LASTELEMENT; 
  memcpy(range, range_stable, sizeof(double) * 4);
}
void infostable(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = p[KAPPA] == 2.0;
}



/* SPACEISOTROPIC stable model for test purposes only */
double stableX(double *x,double *p){
  double z;
  z = x[0] * x[0] + x[1] * x[1];
  return (z==0.0) ? 1.0 : exp(-pow(z, 0.5 * p[KAPPA]));
}
double DstableX(double *x, double *p){
  register double y, z;
  z = x[0] * x[0] + x[1] * x[1];
  if (z == 0.0) return ((p[KAPPA] > 1.0) ? 0.0 : (p[KAPPA] < 1.0) ? INFTY : 1.0);
  y = pow(z, 0.5 * p[KAPPA] - 1.0);
  return -p[KAPPA] * fabs(x[0]) * y * exp(- y * z);
}
/* END SPACEISOTROPIC stable model for test purposes only */


/* Stein */
// see Hypermodel.cc


/* stein space-time model */
double SteinST1(double *x, double *p){
    // kappa1 : nu
    // kappa2 : tau
    // kappa3-5 : z1-z3
/* 2^(1-nu) / Gamma(nu) [ h^nu K_nu(h) - 2 * tau (x T z) t h^{nu-1} K_{nu-1}(h) /
   (2 nu + d + 1) ]
*/
    register double s;
    double z, logconst;
    int d, time, dim;   
    static double nu=RF_INF;
    static double loggamma;
    dim = (int) p[EFFECTIVEDIM];
    time = dim - 1;

    s = x[time] * x[time];
    z = 0.0;
    for (d=0; d<time; d++) {
	s += x[d] * x[d];
	z += x[d] * p[KAPPA3 + d];
    }
    if ( s==0.0 ) return 1.0;
    s = sqrt(s);

    if (nu != p[KAPPA1]) {
      nu = p[KAPPA1];
      loggamma = lgammafn(nu);
    }
    logconst = (nu - 1.0) * log(0.5 * s)  - loggamma;
    return 
	s * exp(logconst + log(bessel_k(s, nu, 2.0)) - s)
	- 2.0 * z * x[time] * exp(logconst + log(bessel_k(s, nu - 1.0, 2.0)) -s) 
	/ (2.0 * nu + p[KAPPA2])  
	;
}
// double ScaleSteinST1(double *p, int scaling) { return 1.0; } 
int kappasSteinST1(int dim) {return 1 + dim;}
void rangeSteinST1(int reduceddim, int *index, double* range){
    static double range_steinST1[20]= 
	{0, RF_INF, 0, 10.0,
	 RF_NAN, RF_INF, RF_NAN, +100,
	 RF_NEGINF, RF_INF, -100000, +100000,
	 RF_NEGINF, RF_INF, -100000, +100000,
	 RF_NEGINF, RF_INF, -100000, +100000,
	};
    *index = RANGE_LASTELEMENT; 
    memcpy(range, range_steinST1, sizeof(double) * 4 * (2 + reduceddim));
    range[4] = range[6] = reduceddim;
}
void infoSteinST1(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = p[KAPPA] >= 2.5;
;
}
int checkSteinST1(double *param, int reduceddim, int dim, SimulationType method) 
{
  double absz;
  int d;
//  printf("%d %d \n", reduceddim + KAPPA1, KAPPA3); assert(false);
 if (method!=Nothing && method!=CircEmbed && method!=Direct)
      return ERRORNOTDEFINED;
  if (param[KAPPA1] >= Besselupperbound[method]) {
      sprintf(ERRORSTRING_OK, "%s and kappa1 < %1.0f", 
		METHODNAMES[method], Besselupperbound[method]);
      sprintf(ERRORSTRING_WRONG,"%1.3f",param[KAPPA1]);
      return ERRORCOVFAILED;
  }
  if (reduceddim == 1) return ERRORNOTDEFINED;
  for (absz=0.0, d=dim + KAPPA1; d>KAPPA2; d--) {
      absz += param[d] * param[d];
//     printf("%d %f %d %d\n", d, param[d], reduceddim, dim);
  }
//  printf("%f %d\n", absz, GENERAL_SKIPCHECKS);
  if (absz > 1.0 + UNIT_EPSILON && !GENERAL_SKIPCHECKS) {
      strcpy(ERRORSTRING_OK, "||z||^2 <= 1");
      sprintf(ERRORSTRING_WRONG, "%f", absz);
      return ERRORCOVFAILED;    
  }
  return checkWhittleMatern(param, reduceddim, dim, method);
}



/* undefined model -- for technical reasons only */
void infoundefined(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 0;
  *CEbadlybehaved = true;
}
int checkundefined(double *param, int reduceddim, int dim, 
		   SimulationType method) {
    return ERRORNOTDEFINED;
}
int checkOK(double *param, int reduceddim, int dim, SimulationType method){
  return NOERROR;
}


/* wave */
double wave(double *x, double *p) {
  if (*x==0.0) { return 1.0;}
  return sin(*x) / *x;
}
double Scalewave(double *p,int scaling) {return 0.302320850755833;}
double spectralwave(double *p ) {
  double x;  
  x = UNIFORM_RANDOM; 
  return sqrt(1.0 - x * x);
}
void rangewave(int reduceddim, int *index, double* range){
  *index = (reduceddim<=3) ? RANGE_LASTELEMENT : RANGE_INVALIDDIM;
}
void infowave(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = 3;
  *CEbadlybehaved = 2;
}



/* Whittle-Matern or Whittle or Besset */ 
double WhittleMatern(double *x, double *p) 
// check calling functions, like hyperbolic and gneiting if any changings !!
{
  static double kappa=RF_INF;
  static double loggamma;
  register double y;
  if (*x==0.0) return 1.0;
  y = fabs(*x); 
  if (kappa != p[KAPPA]) {
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

double TBM2WhittleMatern(double *x, double *p) 
{
  if (p[KAPPA]==0.5) return TBM2exponential(x, p);
  assert(false);
}

double TBM3WhittleMatern(double *x, double *p)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double kappa=RF_INF;
  static double loggamma;
  register double y,loghalfy;
  if (*x==0.0) return 1.0;
  y = fabs(*x); 
  loghalfy = log(0.5 * y);
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    loggamma = lgammafn(kappa);
  }
  return  
    (2.0 * exp(kappa * loghalfy -loggamma + log(bessel_k(y, kappa, 2.0)) - y)
     - (4.0 * exp((kappa + 1.0) * loghalfy - loggamma + 
		  log(bessel_k(y, kappa - 1.0, 2.0)) - y)));
}  

double DWhittleMatern(double *x, double *p)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double kappa=RF_INF;
  static double loggamma;
  register double y;
  if (*x==0.0) 
    return ((p[KAPPA] > 0.5) ? 0.0 : (p[KAPPA] < 0.5) ? INFTY : 1.253314137);
  y = fabs( *x); 
  if (kappa!=p[KAPPA]) {
    kappa=p[KAPPA];
    loggamma = lgammafn(kappa);
  }
  return -2.0 * exp(kappa * log(0.5 * y) - loggamma + 
		    log(bessel_k(y, kappa - 1.0, 2.0)) - y);
}

double DDWhittleMatern(double *x, double *p)
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

int checkWhittleMatern(double *param, int reduceddim, int dim, 
		       SimulationType method) { 
  static double spectrallimit=0.17;
  if (param[KAPPA] >= Besselupperbound[method]) {
      sprintf(ERRORSTRING_OK, "%s and kappa < %1.0f", 
	      METHODNAMES[method], Besselupperbound[method]);
      sprintf(ERRORSTRING_WRONG,"%1.3f",param[KAPPA]);
      return ERRORCOVFAILED;
  }

  switch(method) {
      case TBM2 : 
	if (param[KAPPA] != 0.5) {
	  strcpy(ERRORSTRING_OK,"1/2");
	  sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
	  return ERRORCOVNUMERICAL;
	}
	break;
      case SpectralTBM : 
        if (param[KAPPA] < spectrallimit) {
	  sprintf(ERRORSTRING_OK,
		  "%f<=kappa (the numerical errors are too big for 0<kappa<%f)",
		  spectrallimit, spectrallimit);
	  sprintf(ERRORSTRING_WRONG,"%f",param[KAPPA]);
	  return ERRORCOVFAILED;
	}
	break;
      case CircEmbedCutoff: case CircEmbedIntrinsic :
	if (reduceddim>2) 
	{
	  strcpy(ERRORSTRING_OK,"genuine total dim<=2, currently");
	  sprintf(ERRORSTRING_WRONG, "%d", reduceddim);
	  return ERRORCOVFAILED;
	}
	break;
      default : {}
  }
  return NOERROR;
}
void rangeWhittleMatern(int reduceddim, int *index, double* range){
  //  2 x length(param) x {theor, pract } 
  *index = RANGE_LASTELEMENT; 
  memcpy(range, range_whittle, sizeof(double) * 4);
}
void infoWhittleMatern(double *p, int *maxdim, int *CEbadlybehaved) {
  *maxdim = INFDIM;
  *CEbadlybehaved = p[KAPPA] >= 2.5;
}


