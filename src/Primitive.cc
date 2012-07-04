/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nonstationary models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc


 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2011 Martin Schlather

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



#include <math.h>
 
#include "RF.h"
#include "primitive.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

//  {min, max} x {theor, pract}  x length(param) 
//static double *range_stable[4] = {0, 2, 0.06, 2};
//double *range_whittle[4]= {0, RF_INF, 1e-2, 10.0};
//double *range_cauchy[4] = {0, RF_INF, 0.09, 10.0};
//static double *range_genCauchy[8] = {0, 2, 0.05, 2, 
//				    0, RF_INF, 0.05, 10.0};



double ScaleOne(cov_model *cov){ return 1.0;} 


double BesselUpperB[Nothing + 1] =
{80, 80, 80, // circulant
 80, 80, R_PosInf, // TBM
 80, 80,    // direct & sequ
 NA_REAL, NA_REAL, NA_REAL,  // GMRF, ave, nugget
 NA_REAL, NA_REAL,   // exotic
 R_PosInf   // Nothing
};

#define LOW_BESSELJ 1e-20
#define LOW_BESSELK 1e-20

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

void Bessel(double *x, cov_model *cov, double *v){
  static double kappaOld=RF_INF;
  static double gamma;
  double kappa = cov->p[0][0], y=*x;

  if  (y <= LOW_BESSELJ) {*v = 1.0; return;} 
  if (kappaOld!=kappa) {
    kappaOld=kappa;
    gamma = gammafn(kappaOld+1.0);
  }
  *v = gamma  * pow(0.5 * y, -kappaOld) * bessel_j(y, kappaOld);
}
int initspectralBessel(cov_model *cov) {
  return (cov->tsdim <= 2) ? NOERROR : ERRORFAILED;
}
void spectralBessel(cov_model *cov, spectral_storage *s, double *e) { 
/* see Yaglom ! */
  // kappa==0.0 ? 1.0 : // not allowed anymore;
	// other wise densityBessel (to define) will not work
   E2(s, sqrt(1.0 - pow(UNIFORM_RANDOM, 1.0 / cov->p[0][0])), e);
}
int checkBessel(cov_model *cov) {
  // Whenever TBM3Bessel exists, add further check against too small kappa! 
  assert(false);
  double kappa = cov->p[0][0];
  int i;
  for (i=0; i<= Nothing; i++) {
    cov->pref[i] *= (ISNA(kappa) || ISNAN(kappa) || kappa < BesselUpperB[i]);
  }
  if (kappa <= 0.0) cov->pref[SpectralTBM] = PREF_NONE;
  return NOERROR;
}
void rangeBessel(cov_model *cov, range_arraytype *ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.5 * ((double) cov->tsdim - 2.0);
  range->max[0] = RF_INF;
  range->pmin[0] = 0.0001 + range->min[0];
  range->pmax[0] = range->pmin[0] + 10.0;
  range->openmin[0] = false;
  range->openmax[0] = true;
  
  int dim = (int) (2.0 * cov->p[0][0] + 2.0);
  range->maxdim = dim >= INFDIM ? INFDIM - 1 : (int) dim;
}


/* Cauchy models */
void Cauchy(double *x, cov_model *cov, double *v){
  double a = cov->p[0][0];
  //printf("%f %f\n", x[0], pow(1.0 + *x * *x, -a));
 *v = pow(1.0 + *x * *x, -a);
}
double ScaleCauchy(cov_model *cov) {
    double a = cov->p[0][0];
    return 1.0 / sqrt(pow(0.05, -1 / a)-1.0); 
}
void TBM2Cauchy(double *x, cov_model *cov, double *v){
  double a = cov->p[0][0], y2, lpy2;
  y2 = *x * *x; 
  lpy2 = 1.0 + y2;
  switch ((int) (a * 2.0 + 0.001)) {// ueber check sichergestellt
      case 1 : *v = 1.0 / lpy2;
      case 3 : *v = (1.0 - y2)/ (lpy2 * lpy2);
      case 5 : *v = (1.0-y2*(2.0+0.333333333333333333333*y2))/(lpy2*lpy2*lpy2);
      case 7 : lpy2 *= lpy2; 
	*v = (1.0- y2*(3.0+y2*(1.0+0.2*y2)))/(lpy2 * lpy2);
      default : assert(false);
  }
}
//void TBM3Cauchy(double *x, cov_model *cov, double *v){
//  double ha = *x * *x, 
//    a = cov->p[0][0];
//  *v = (1.0 + (1.0 - 2.0 * a) * ha) * pow(1.0 + ha, -a - 1.0);
//}
void DCauchy(double *x, cov_model *cov, double *v){
  double y=*x, a = cov->p[0][0];
  *v = (-2.0 * a * y) * pow(1.0 + y * y, -a - 1.0);
}
void DDCauchy(double *x, cov_model *cov, double *v){
  double ha = *x * *x, a = cov->p[0][0];
  *v = 2.0 * a * ((2.0 * a + 1.0) * ha - 1.0) * pow(1.0 + ha, -a - 2.0);
}
void invCauchySq(double*x, cov_model *cov, double *v){
  double
      a = cov->p[0][0];
  if (*x == 0.0) *v = RF_INF;
  else *v = pow(*x, -1.0 / a) - 1.0;
}
int checkCauchy(cov_model *cov){
  double a = cov->p[0][0];
  if (a != 0.5 && a != 1.5 && a != 2.5 && a != 3.5) 
    cov->pref[TBM2] = PREF_NONE;
  return NOERROR;
}
void rangeCauchy(cov_model *cov, range_arraytype *ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 0.09;
  range->pmax[0] = 10.0;
  range->openmin[0] = true;
  range->openmax[0] = true;
}
void coinitCauchy(cov_model *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0; //  q[CUTOFF_A] 
  li->msg[0] = MSGLOCAL_JUSTTRY;
}
double DrawLogMixCauchy(cov_model *cov, mpp_storage *s) { //better GR 3.381.4 ?? !!!!
  // split into parts <1 and >1 ?>
  return log(-log(1.0 -UNIFORM_RANDOM));
}
double LogMixWeightCauchy(double *x, cov_model *cov, mpp_storage *s) {
  double a = cov->p[0][0];
  // da dort 1/w vorgezogen 
  return 0.5 * (a - 1.0) * s->loginvscale - 0.5 * lgammafn(a);
}


/* another Cauchy model */
void Cauchytbm(double *x, cov_model *cov, double *v){
  double ha, a = cov->p[0][0], b = cov->p[1][0], c = cov->p[2][0], y=*x;
  if (y==0) {
    *v = 1.0;
  } else {
    ha=pow(y, a);
    *v = (1.0 + (1.0 - b / c) * ha) * pow(1.0 + ha, -b / a - 1.0);
  }
}
//void TBM3Cauchytbm(double *x, cov_model *cov, double *v){
//  double ha, bg, a = cov->p[0][0], b = cov->p[1][0], c = cov->p[2][0];
//  ha=pow(*x, a);
//  bg=b / c;//  *v = 
//    (1 + ha * (1 - bg * (1 + a) + (1 - b) * (1 + (1 - bg) * ha))) *
//    pow(1 + ha, -b / a - 2.0);
//}
void DCauchytbm(double *x, cov_model *cov, double *v){
  double y= *x, ha, a = cov->p[0][0], b = cov->p[1][0], c = cov->p[2][0];
  if (y == 0.0) *v = 0.0; // WRONG VALUE, but multiplied 
  else {                                  // by zero anyway
    ha = pow(y, a - 1.0);
    *v = b *  ha * (-1.0 - a / c  + ha * y * (b / c - 1.0)) *
      pow(1.0 + ha * y, -b /a - 2.0);
  }
}


void rangeCauchytbm(cov_model *cov, range_arraytype *ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = 0.05;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = 0.05;
  range->pmax[1] = 10.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
 
  range->min[2] = (double) cov->tsdim;
  range->max[2] = RF_INF;
  range->pmin[2] = range->min[2];
  range->pmax[2] = range->pmin[2] + 10.0;
  range->openmin[2] = false;
  range->openmax[2] = true;

  range->maxdim = INFDIM; 
}



/* circular model */
void circular(double *x, cov_model *cov, double *v) {
  double y = *x;
  *v = (y >= 1.0) ? 0.0 
    : 1.0 - (2.0 * (y * sqrt(1.0 - y * y) + asin(y))) * INVPI;
}
//double Scalecircular(cov_model *cov, int scaling) {
//  return 1.138509531721630274603;
//}
// spectral measure, see Lantue !! 
void Dcircular(double *x, cov_model *cov, double *v){
  double y = *x * *x;
  *v = (y >= 1.0) ? 0.0 : -4 * INVPI * sqrt(1.0 - y);
}
void rangecircular(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->maxdim = 2; 
  range->finiterange = true;
}


// constant (matrix); 
void constant(double *x,  cov_model *cov, double *v){
    listoftype *list= (listoftype*) (cov->p[0]);
    int i,j,
	// *ncol = list->ncol,
      *nrow = list->nrow; // anzahl listen elemente
   int vdim = cov->vdim,
	lx = nrow[LIST_ELEMENT] / cov->vdim, //nrow = #pts * vdim
	ly = nrow[LIST_ELEMENT] * lx;
   double **p = list->p,
       *X = (p[LIST_ELEMENT] + 
	     CovMatrixRow + CovMatrixCol * nrow[LIST_ELEMENT]),
	*Y=X,
	*w = v;
//   printf("const %d %d %d listelement=%d %d %d", 
//	  CovMatrixRow, lx, CovMatrixCol, LIST_ELEMENT, nrow[LIST_ELEMENT],
//	  cov->vdim);
    if (CovMatrixRow >= lx || CovMatrixCol >= lx) 
      ERR("constant model: indices out of range");
    if (LIST_ELEMENT >= cov->nrow[0])
      ERR("constant model: list index out of range");
    for (i=0; i<vdim; i++, Y+=ly) {
      for (X=Y, j=0; j<vdim; j++, w++, X+=lx) {
	*w = *X;
      }
    }
}

void constant_nonstat(double *x, double *y, cov_model *cov, double *v){
    constant(x, cov, v);
}

int checkconstant(cov_model *cov) {
  listoftype *list= (listoftype*) (cov->p[0]);
  int info, err, total, i, vdim,
      m = cov->nrow[0],
      *ncol = list->ncol,
      *nrow = list->nrow; // anzahl listen elemente
  double *dummy,
      **p = list->p;

  kdefault(cov, 1, 1); // vdim !
  cov->vdim = vdim = ((int*) cov->p[1])[0];

//  printf("%d \n", cov->p[1][0]);
//  assert(false);

  for (i=0; i<m; i++) {
    if (nrow[i] != ncol[i] || cov->nrow[i] == 0) {
      PRINTF("square matrix in `constant' expected\n");
      return ERRORFAILED;
    }
    if ((int) (nrow[i] / vdim) * vdim != nrow[i]) {
      ERR("size of matrix is not a multiple of the multivariate dimension");
    }

    // check whether positive eigenvalue  
    total = nrow[i] * ncol[i] * sizeof(double); // NICHT weiter vorne
    dummy = (double*) malloc(total);
    memcpy(dummy, p[i], total);

//  printf("vdim=%d %d\n", cov->vdim, total);
//  for (err =0; err < cov->vdim*cov->vdim; err++) printf(" %f ", dummy[err]);
// printf("\n");

    F77_CALL(dpofa)(dummy, &(cov->vdim), &(cov->vdim), &info); // cholesky
    free(dummy);
    if (info != 0) {
	ERR("matrix does not seem to be strictly positive definite");
    }
  }

  err = checkkappas(cov);
  return err;
}
void rangeconstant(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;
  range->min[0] = RF_NEGINF;
  range->max[0] = RF_INF;
  range->pmin[0] = -1e10;
  range->pmax[0] = 1e10;
  range->openmin[0] = true;
  range->openmax[0] = true;
  
  range->min[1] = 1;
  range->max[1] = 9999;
  range->pmin[1] = 1;
  range->pmax[1] = 9999;
  range->openmin[1] = false;
  range->openmax[1] = false;
}



/* cone */
// no covariance function defined yet; see MPPFcts.cc for additive method
int checkcone(cov_model *cov) {
  if (cov->tsdim != 2 || cov->p[2][0] + cov->p[1][0]==0.0) 
    cov->pref[RandomCoin] = PREF_NONE;
  return NOERROR;
}

void rangecone(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 1.0;
  range->pmin[0] = 0.0;
  range->pmax[0] = 1.0 - UNIT_EPSILON;
  range->openmin[0] = false;
  range->openmax[0] = true;

  
  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = UNIT_EPSILON;
  range->pmax[1] = 10.0;
  range->openmin[1] = false;
  range->openmax[1] = true;
 
  range->min[2] = 0.0;
  range->max[2] = RF_INF;
  range->pmin[2] = UNIT_EPSILON;
  range->pmax[2] = 10.0;
  range->openmin[2] = false;
  range->openmax[2] = true;

  range->finiterange = true;
  range->maxdim = 3;
}


/* coxgauss, cmp with nsst1 !! */
// see NonIsoCovFct.cc

/* cubic */
void cubic(double *x, cov_model *cov, double *v) {
  double y=*x, y2=y * y;
  *v = (y >= 1.0) ? 0.0 : (1.0 + (((0.75 * y2 - 3.5) * y2 + 8.75) * y - 7) * y2);
}
void Dcubic(double *x, cov_model *cov, double *v) { 
  double y=*x, y2=y * y;
  *v = (y >= 1.0) ? 0.0 : y * (-14.0 + y * (26.25 + y2 * (-17.5 + 5.25 * y2)));
}
//double Scalecubic(cov_model *cov, int scaling) {return 1.44855683156829;}
void rangecubic(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->maxdim = 3; 
  range->finiterange = true;
}


/* cutoff */
// see Hypermodel.cc

/* dagum */
void dagum(double *x, cov_model *cov, double *v){
  double gamma = cov->p[1][0], 
    beta=cov->p[0][0];
  *v = 1.0 - pow((1 + pow(*x, -beta)), -gamma/beta);
}
double Scaledagum(cov_model *cov){ 
  double gamma = cov->p[1][0], 
    beta=cov->p[0][0];
  return pow(pow(0.95, - beta / gamma ) - 1.0, -1.0 / beta);
} 
void Ddagum(double *x, cov_model *cov, double *v){
  double y=*x, xd, 
    gamma = cov->p[1][0], 
    beta=cov->p[0][0];
  xd = pow(y, -beta);
  *v = -gamma * xd / y * pow(1 + xd, -gamma/ beta -1);
}
void rangedagum(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 1.0;
  range->pmin[0] = 0.01;
  range->pmax[0] = 1.0;
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = 1.0;
  range->pmin[1] = 0.01;
  range->pmax[1] = 1.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
}


/*  damped cosine -- derivative of e xponential:*/
void dampedcosine(double *x, cov_model *cov, double *v){
  double y = *x, lambda = cov->p[0][0];
  *v = exp(-y * lambda) * cos(y);
}
double Scaledampedcosine(cov_model *cov){ 
  return MINUSINVLOG005; 
} 
void Ddampedcosine(double *x, cov_model *cov, double *v){
  double y = *x, lambda = cov->p[0][0];
  *v = - exp(-lambda*y) * (lambda * cos(y) + sin(y));
}
void rangedampedcosine(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  switch(cov->tsdim) {
      case 1 : range->min[0] = 0.0;
      case 2 : range->min[0] = 1.0;
      case 3 : range->min[0] = RF_M_SQRT_3;
      default:  range->min[0] = -RF_INF;
  }
  range->max[0] = RF_INF;
  range->pmin[0] = 0.05;
  range->pmax[0] = range->min[0] + 1e-10;
  range->openmin[0] = false;
  range->openmax[0] = true;

  double lambda = cov->p[0][0];
  range->maxdim =
    (lambda < 0) ? 0 : (lambda < 1) ? 1 : (lambda < M_SQRT_3) ? 2 : 3;
}


/* De Wijsian */
void dewijsian(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0];
  *v = -log(1.0 + pow(*x, alpha));
}
void rangedewijsian(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = UNIT_EPSILON;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false; 
}


/* De Wijsian B */
void DeWijsian(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0],
    range = cov->p[1][0];
  *v = (*x >= range) ? 0.0 : 1.0 -
    log(1.0 + pow(*x, alpha)) / log(1.0 + pow(range, alpha));
}
double ScaleDeWijsian(cov_model *cov){ 
    double range = cov->p[1][0];
    return range; // gibt an, wann die Funktion null ist
} 
void rangeDeWijsian(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 1.0;
  range->pmin[0] = UNIT_EPSILON;
  range->pmax[0] = 1.0;
  range->openmin[0] = true;
  range->openmax[0] = false; 

  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = UNIT_EPSILON;
  range->pmax[1] = 1000;
  range->openmin[1] = true;
  range->openmax[1] = true; 
  range->maxdim = 1;
}


/* exponential model */
void exponential(double *x, cov_model *cov, double *v){
  *v = exp(- *x);
  // printf("exp: %f %f\n", *x, *v);
 }
double Scaleexponential(cov_model *cov){ return MINUSINVLOG005; } 
void TBM2exponential(double *x, cov_model *cov, double *v) 
{
  double y = *x;
  *v = (y==0.0) ?  1.0 : 1.0 - PIHALF * y * I0mL0(y);
}
void Dexponential(double *x, cov_model *cov, double *v){
//  printf("\n ******* dexp %f %f\n", double *- exp(-fabs( *x)));
  *v = - exp(- *x);
}
void DDexponential(double *x, cov_model *cov, double *v){
 *v = exp(-*x);
}
void invexponentialSq(double *x, cov_model *cov, double *v){
 if (*x == 0.0) *v = RF_INF;
 else {
   *v = -log(*x);
   *v *= *v;
 }
}

double densityexponential(double *x, cov_model *cov) {
    double x2, dim12, dim=cov->tsdim;
    int d;
    x2 = x[0] * x[0];
    for (d=1; d<dim; d++) x2 += x[d] * x[d];
    dim12 = 0.5 * ((double) (dim + 1));
    return gammafn(dim12) * pow(M_PI * (1.0 + x2), -dim12);
}

int initspectralexponential(cov_model *cov) {
  spec_covstorage *cs = &(cov->spec);
  if (cov->tsdim <= 2) return NOERROR;
  cs->density = densityexponential;
  return search_metropolis(cov);
}

void spectralexponential(cov_model *cov, spectral_storage *s, double *e) {
  if (cov->tsdim <= 2) {
    double y;
    y = 1.0 - UNIFORM_RANDOM;
    E2(s, sqrt(1.0 / (y * y) - 1.0), e);
  } else {
    metropolis(cov, e);
  }
}

void rangeexponential(cov_model *cov, range_arraytype* ra){ }
int checkexponential(cov_model *cov) {
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  return NOERROR;
}
int hyperexponential(double radius, double *center, double *rx,
		     cov_model *cov, bool simulate, 
		     double ** Hx, double ** Hy, double ** Hr)
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
  long i, p, 
    q = RF_NAN;
  int k, err;
  
  if (cov->tsdim==2) {
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
    if ((hx = *Hx = 
	 (double *) malloc(sizeof(double) * (p + 8 * sizeof(int))))==NULL){
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    } 
     // see also bits in RFhyperplan.cc, line 437 about.
    if ((hy=*Hy=
	 (double *) malloc(sizeof(double) * (p + 8 * sizeof(int))))==NULL){
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }  
    if ((hr=*Hr=
	 (double *) malloc(sizeof(double) * (p + 8 * sizeof(int))))==NULL){
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
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
    // tsdim = 1  --  not programmed yet
    error("impossible dimension in hyperexponential\n");
  }
  return q;

 ErrorHandling:
  PRINTF("error=%d\n", err);
  error("impossible dimension in hyperexponential\n");
  return NA_INTEGER;
}
void coinitExp(cov_model *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0; //  q[CUTOFF_A] 
  li->msg[0] = MSGLOCAL_OK;
}
void ieinitExp(cov_model *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0;
  li->msg[0] = MSGLOCAL_OK;
}
double DrawLogMixExp(cov_model *cov, mpp_storage *s) {
  // GR 3.325: int_-infty^infty exp(x^2/2 - b/x^2) = exp(-\sqrt b)
  double x = GAUSS_RANDOM(1.0);
  return log(1.0 / (x * x));
}
double LogMixWeightExp(double *x, cov_model *cov, mpp_storage *s) {
  return 0.0;
}


// Brownian motion 
void fractalBrownian(double *x, cov_model *cov, double *v) {
  double alpha = cov->p[0][0];
  *v = - pow(*x, alpha);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
//begin
/* fractalBrownian: first derivative at t=1 */
void DfractalBrownian(double *x, cov_model *cov, double *v) 
{// FALSE VALUE FOR *x==0 and  alpha < 1
  double alpha = cov->p[0][0];
  *v = (*x == 0.0) ? 0.0 : -alpha * pow(*x, alpha - 1.0);
}
/* fractalBrownian: second derivative at t=1 */
void DDfractalBrownian(double *x, cov_model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  double alpha = cov->p[0][0];
  *v = (*x == 0.0) ? 0.0 : 
	-alpha * (alpha - 1.0) * pow(*x, alpha - 2.0);
}
void invfractalBrownianSq(double *x, cov_model *cov, double *v) {
  double alpha = cov->p[0][0];
  *v = - pow(*x, 2.0 / alpha);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
void rangefractalBrownian(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = UNIT_EPSILON;
  range->pmax[0] = 2.0 - UNIT_EPSILON;
  range->openmin[0] = true;
  range->openmax[0] = false;
}
void ieinitBrownian(cov_model *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = (cov->tsdim <= 2)
    ? ((cov->p[0][0] <= 1.5) ? 1.0 : 2.0)
    : ((cov->p[0][0] <= 1.0) ? 1.0 : 2.0);
  li->msg[0] = cov->tsdim <= 3 ? MSGLOCAL_OK : MSGLOCAL_JUSTTRY;
}



/* FD model */
void FD(double *x, cov_model *cov, double *v){
  double kappa = cov->p[0][0], y, d, k, skP1;
  static double dold=RF_INF;
  static double kold, sk;
  d = kappa * 0.5;
  y = *x;
  k = trunc(y);
  if (dold!=d || kold > k) {
    sk = 1;
    kold = 0.0;
  }
  // sign (-1)^k is (kold+d), 16.11.03, checked. 
  for (; kold<k; kold += 1.0) sk =  sk * (kold + d) / (kold + 1.0 - d);
  dold = d;
  kold = k;
  if (k == y) {
    *v = sk;
  } else {
    skP1 = sk * (kold + d) / (kold + 1.0 - d);
    *v = sk + (y - k) * (skP1 - sk);
  }
}
void rangeFD(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = -1.0;
  range->max[0] = 1.0;
  range->pmin[0] = range->min[0] + UNIT_EPSILON;
  range->pmax[0] = range->max[0] - UNIT_EPSILON;
  range->openmin[0] = false;
  range->openmax[0] = true;
  range->maxdim = 1;
}



/* fractgauss */
void fractGauss(double *x, cov_model *cov, double *v){
  double y = *x, kappa = cov->p[0][0];
  *v = (y == 0.0) ? 1.0 : 
    0.5 *(pow(y + 1.0, kappa) - 2.0 * pow(y, kappa) + pow(fabs(y - 1.0), kappa));
}
void rangefractGauss(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = UNIT_EPSILON;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;
  range->maxdim = 1;
}


/* Gausian model */
void Gauss(double *x, cov_model *cov, double *v) {
  *v = exp(- *x * *x);
}
double ScaleGauss(cov_model *cov) {return SQRTINVLOG005;}
void DGauss(double *x, cov_model *cov, double *v) {
  double y = *x; 
  *v = -2.0 * y * exp(- y * y);
//  printf("DG %f %f\n", *x, *v);
}
void DDGauss(double *x, cov_model *cov, double *v) {
  double y = *x * *x; 
  *v = (4.0 * y - 2.0)* exp(- y);
//  printf("ZG %f %f\n", *x, *v);
}
void D3Gauss(double *x, cov_model *cov, double *v) {
  double y = *x * *x; 
  *v = *x * (12 - 8 * y) * exp(- y);
//  printf("ZG %f %f\n", *x, *v);
}
void D4Gauss(double *x, cov_model *cov, double *v) {
  double y = *x * *x; 
  *v = ((16 * y - 48) * y + 12) * exp(- y);
//  printf("ZG %f %f\n", *x, *v);
}
void invGaussSq(double *x, cov_model *cov, double *v) {
  if (*x == 0.0) *v = RF_INF;
  else *v = -log(*x);
}
double densityGauss(double *x, cov_model *cov) {  
    int d, dim=cov->tsdim;
    double x2=0.0;
    for (d=0; d<dim; d++) x2 += x[d] * x[d];
    return exp(- 0.25 * x2 - (double) dim * (M_LN2 + M_LN_SQRT_PI));
}
int initspectralGauss(cov_model *cov) {
  spec_covstorage *cs = &(cov->spec);
  if (cov->tsdim <= 2) return NOERROR;
  cs->density = densityGauss;
  return search_metropolis(cov);
}
void spectralGauss(cov_model *cov, spectral_storage *s, double *e) {   
  if (cov->tsdim <= 2) {
    E2(s, 2.0 * sqrt(-log(1.0 - UNIFORM_RANDOM)), e);
  } else {
    metropolis(cov, e);
  }
}
void rangeGauss(cov_model *cov, range_arraytype* ra){ }
double DrawLogMixGauss(cov_model *cov, mpp_storage *s) {
  return 1.0;
}
double LogMixWeightGauss(double *x,cov_model *cov, mpp_storage *s) {
  return 0.0;
}

/* generalised fractal Brownian motion */
void genBrownian(double *x, cov_model *cov, double *v) {
  double alpha = cov->p[0][0], delta =  cov->p[1][0];
  *v = - pow(pow(*x, alpha) + 1.0, delta);
  //this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
void rangegenBrownian(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = UNIT_EPSILON;
  range->pmax[0] = 2.0 - UNIT_EPSILON;
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = 1.0;
  range->pmin[1] = UNIT_EPSILON;
  range->pmax[1] = 1.0 - UNIT_EPSILON;
  range->openmin[1] = true;
  range->openmax[1] = false;
}


/* gencauchy */
void generalisedCauchy(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0];
//  printf("%f\n", alpha);
//  printf("%f\n", beta);
//  printf("%f\n", x[0]);
  *v = pow(1.0 + pow(*x, alpha), -beta/alpha);
  
  // print("gc=%f %f %f %f ", *x, *v, alpha, beta);

  // assert(false);
}
double ScalegeneralisedCauchy(cov_model *cov) {
  double alpha = cov->p[0][0], beta=cov->p[1][0];
  return pow(pow(0.05, -alpha / beta) - 1.0, -1.0 / alpha);
  // MLE works much better with 0.01 then with 0.05
}
void DgeneralisedCauchy(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0], ha, y=*x;
  if (y ==0.0) {
    *v = ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? -INFTY : -beta); 
  } else {
    ha=pow(y, alpha - 1.0);
    *v = -beta * ha * pow(1.0 + ha * y, -beta / alpha - 1.0);
  }
}
void DDgeneralisedCauchy(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0], ha, y=*x;
  if (y ==0.0) {
    *v = ((alpha==2.0) ? beta * (beta + 1.0) : INFTY); 
  } else {
    ha=pow(y, alpha);
    *v = beta * ha / (y * y) * (1.0 - alpha + (1.0 + beta) * ha)
      * pow(1.0 + ha, -beta / alpha - 2.0);
  }
}
void invgeneralisedCauchySq(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0];
//  printf("%f\n", alpha);
//  printf("%f\n", beta);
//  printf("%f\n", x[0]);
  if (*x == 0.0) *v = RF_INF;
  else *v = pow(pow(*x, -alpha/beta) - 1.0, 2.0/alpha);
  // assert(false);
}
int checkgeneralisedCauchy(cov_model *cov){
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  return NOERROR;
}

void rangegeneralisedCauchy(cov_model *cov, range_arraytype* ra) {
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = 0.05;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = 0.05;
  range->pmax[1] = 10.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
}

void coinitgenCauchy(cov_model *cov, localinfotype *li) {
  double thres[2] = {0.5, 1.0}, alpha=cov->p[0][0]; 
  if (alpha <= thres[0]) {
    li->instances = 2;
    li->value[0] = 0.5;
    li->value[1] = 1.0;
    li->msg[0] = li->msg[1] = MSGLOCAL_OK;
  } else {
    li->instances = 1;
    li->value[0] = 1.0; //  q[CUTOFF_A] 
    li->msg[0] = (alpha <= thres[1]) ? MSGLOCAL_OK : MSGLOCAL_JUSTTRY;
  } 
}
void ieinitgenCauchy(cov_model *cov, localinfotype *li) {
  li->instances = 1;
  if (cov->p[0][0] <= 1.0) {
    li->value[0] = 1.0;
    li->msg[0] = MSGLOCAL_OK;
  } else {
    li->value[0] = 1.5;
    li->msg[0] = MSGLOCAL_JUSTTRY;
  }
}



/* gencauchy */
void epsC(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0], eps=cov->p[2][0];
//  printf("%f\n", alpha);
//  printf("%f\n", beta);
//  printf("%f\n", x[0]);
  *v = pow(eps + pow(*x, alpha), -beta/alpha);
  // assert(false);
}
void DepsC(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0], eps=cov->p[2][0], ha, y=*x;
  if (y ==0.0) {
    *v = (eps == 0.0) ? -INFTY : 
      ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? -INFTY : -beta); 
  } else {
    ha=pow(y, alpha - 1.0);
    *v = -beta * ha * pow(eps + ha * y, -beta / alpha - 1.0);
  }
}
void DDepsC(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0], eps=cov->p[2][0], ha, y=*x;
  if (y ==0.0) {
    *v = (eps == 0.0) ? INFTY : ((alpha==2.0) ? beta * (beta + 1.0) : INFTY); 
  } else {
    ha=pow(y, alpha);
    *v = beta * ha / (y * y) * ( (1.0 - alpha) * eps + (1.0 + beta) * ha)
      * pow(eps + ha, -beta / alpha - 2.0);
  }
}
void invepsCSq(double *x, cov_model *cov, double *v){
  double alpha = cov->p[0][0], beta=cov->p[1][0], eps=cov->p[2][0];
  if (*x == 0.0) *v = RF_INF;
  else *v = pow(pow(*x, -alpha/beta) - eps, 2.0/alpha);
}
int checkepsC(cov_model *cov){
  double eps=cov->p[2][0];
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  kdefault(cov, 0, 1.0); 
  kdefault(cov, 1, 1.0); 
  kdefault(cov, 2, 0.0); 
  if (ISNA(eps) || ISNAN(eps) || eps == 0.0) {
    cov->statIn=GENERALISEDCOVARIANCE;
    int i;
    for (i=CircEmbed; i<Nothing; i++) cov->pref[i] = PREF_NONE;
  }
  
  return NOERROR;
}

void rangeepsC(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = 0.05;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = 0.05;
  range->pmax[1] = 10.0;
  range->openmin[1] = true;
  range->openmax[1] = true;

  range->min[2] = 0.0;
  range->max[2] = RF_INF;
  range->pmin[2] = 0.0;
  range->pmax[2] = 10.0;
  range->openmin[2] = false;
  range->openmax[2] = true;
}


/* gengneiting */
void genGneiting(double *x, cov_model *cov, double *v)
{
  int alpha = ((int*)cov->p[0])[0];
  double delta=cov->p[1][0], y=*x, s;
  if (y >= 1.0) {
    *v = 0.0; 
    return;
  }
  s = delta + (double) alpha;
  switch (alpha) {
  case 1:
    *v =  pow ((1.0 - y), s) * (1.0 + s*y);   
    break;
  case 2:
    *v = pow ((1.0 - y), s) * (1.0 + y * (s + y * (s*s-1.0)*0.3333333333333333));
    break;
  case 3:
    *v = pow ((1.0 - y), s) * 
      (1.0 + y * (s + y * ( 0.2 * (2.0*s*s - 3.0)
			    + y * (s*s-4.0)*s*0.0666666666666666666)));
    break;
  default : assert(false);  
  }
}
void DgenGneiting(double *x, cov_model *cov, double *v)
{
  int alpha = ((int*) cov->p[0])[0];
  double delta=cov->p[1][0], y=*x, s;
  if (y >= 1.0) {
    *v = 0.0; 
    return;
  }
  s = delta + (double) alpha;
  switch (alpha) {
      case 1:
	*v =  - pow(1.0 - y, s - 1.0) * y * s  * (s + 1.0); 
	return;
      case 2:
	*v = - pow(1.0 - y, s - 1.0) *
	  y * (0.333333333333333333 * s * s + 0.66666666666666666667 + s +
	       y * 0.333333333333333333 * (s * s - 1.0) * (delta + 4.0) ); 
	return;
      case 3:
	*v = -  pow(1.0 - y, s - 1.0) *
	  0.2 * y * (6 + s * (5.0 + s) +
		     y * (-6 + s * (1.0 + s * (4.0 + s)) +
			  0.33333333333333333 * y * s *(-12.0 + s * 
							(-4.0 + s * (3.0 + s)))
		       ));
	return;
      default : assert(false);   
  }
}
void rangegenGneiting(cov_model *cov, range_arraytype* ra){
  int index;
  range_type *range = ra->ranges + 0;
  ra->n = 3;
  for (index=0; index<ra->n; index++, range++) {
    range->max[0] = range->min[0] =  range->pmin[0] = range->pmax[0] = index;
    range->openmin[0] = false;
    range->openmax[0] = false;
    
    range->min[1] = 0.5 * (double) (cov->tsdim + 2 * index + 1);
    range->max[1] =  RF_INF;
    range->pmin[1] = range->min[1];
    range->pmax[1] = range->pmin[1] + 10.0;
    range->openmin[1] = false;
    range->openmax[1] = true;
    double dim = (2.0 * cov->p[1][0] + 1.0 + 2.0 * ((int *) cov->p[0])[0]);  
    range->maxdim = (dim >= INFDIM) ? INFDIM-1 : (int) dim;
    range->finiterange = true;
  }
}


/* Gneiting's functions -- alternative to Gaussian */
// #define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
#define NumericalScale 0.301187465825
void Gneiting(double *x, cov_model *cov, double *v){ 
  double y=*x * NumericalScale, oneMy8;
  if (y >= 1.0) {
    *v = 0.0;
  } else {
    oneMy8 = 1.0-y; oneMy8*=oneMy8; oneMy8*=oneMy8; oneMy8*=oneMy8;
    *v =((1.0+y * ( 8.0 + y * (25.0 + 32.0 *y)))*oneMy8);
  }
}
//double ScaleGneiting(cov_model *cov, int scaling) {return 0.5854160193;}
void DGneiting(double *x, cov_model *cov, double *v){ 
  double y=*x * NumericalScale, oneMy7;
  if (y >= 1.0) {
    *v = 0.0;
  } else {
    oneMy7 = 1.0-y; oneMy7 *= oneMy7; oneMy7 *= oneMy7 * oneMy7 * (1.0-y);
    *v = (-y) * ( 22.0 + y * (154.0 + y * 352.0)) * oneMy7 * NumericalScale;
  }
}
void rangeGneiting(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->maxdim = 3; 
  range->finiterange = true;
}


/* hyperbolic */
void hyperbolic(double *x, cov_model *cov, double *v){ 
  double kappa = cov->p[0][0], lambda=cov->p[1][0], delta=cov->p[2][0];
  static double kappaOld = RF_INF;
  static double lambdaOld= RF_INF;
  static double deltaOld = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double kappay;
  double y=*x;
  if (y==0.0) { 
    *v = 1.0;
    return;
  }
  if (delta==0) { // whittle matern
      if (kappa > 80) warning("extremely imprecise results due to kappa>80");
    y *= kappa;
    Whittle(&y, cov, v);
  } else if (kappa==0) { //cauchy   => KAPPA2 < 0 !
    y /= delta;
    /* note change in sign as KAPPA2<0 */
    *v = pow(1 + y * y, lambda); 
  } else {
    if ((kappa!=kappaOld) || (lambda!=lambdaOld) || (delta!=deltaOld)) {
      kappaOld = kappa; 
      lambdaOld = lambda;
      deltaOld = delta;
      deltasq = deltaOld * deltaOld;
      kappadelta = kappaOld * deltaOld;
      logconst = kappadelta - log(bessel_k(kappadelta, lambdaOld, 2.0)) 
	- lambda * log(deltaOld);
    }
    y=sqrt(deltasq + y * y);  
    kappay = kappa * y;
    *v = exp(logconst + lambda * log(y) + log(bessel_k(kappay, lambda, 2.0)) 
	     - kappay);
  }
}

void Dhyperbolic(double *x, cov_model *cov, double *v)
{ 
  double kappa = cov->p[0][0], lambda=cov->p[1][0], delta=cov->p[2][0];
  static double kappaOld = RF_INF;
  static double lambdaOld= RF_INF;
  static double deltaOld = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double y = *x;
  double s,kappa_s,logs;
  if (y==0.0) { 
    *v = 1.0;
    return;
  }
  if (delta==0) { // whittle matern
    y *= kappa;
    DWhittle(&y, cov, v);
  } else if (kappa==0) { //cauchy
    double ha;
    y /= delta;
    ha = y * y;
    /* note change in sign as KAPPA2<0 */
    *v = 2.0 * lambda * fabs(y) * pow(1.0 + ha, lambda - 1.0);
  } else {
    if ((kappa!=kappa) || (lambda!=lambda) || (delta!=delta)) {
      kappaOld = kappa; 
      lambdaOld= lambda;
      deltaOld = delta;
      deltasq = deltaOld * deltaOld;
      kappadelta = kappaOld * deltaOld;
      logconst = kappadelta - log(bessel_k(kappadelta,lambdaOld,2.0)) 
	- lambdaOld * log(deltaOld);
    }
    s=sqrt(deltasq + y * y);
    kappa_s = kappa * s;
    logs = log(s);  
    *v = - y * kappa*exp(logconst + (lambda-1.0) * logs 
		      +log(bessel_k(kappa_s,lambda-1.0,2.0))-kappa_s);
  }
}
void invhyperbolicSq(double *x, cov_model *cov, double *v){ 
/*
  double alpha = cov->p[0][0], lambda=cov->p[1][0], delta=cov->p[2][0];
  static double kappaOld = RF_INF;
  static double lambdaOld= RF_INF;
  static double deltaOld = RF_INF;
  static double deltasq;
  static double kappadelta;
  static double logconst;
  double kappay;
  double y=*x;
*/
  ERR("inv hyperbolic not programmed yet");
}
int checkhyperbolic(cov_model *cov){
  double kappa = cov->p[0][0], lambda=cov->p[1][0], delta=cov->p[2][0];
  char msg[255];
  int i;
  for (i=0; i<= Nothing; i++)
    cov->pref[i] *= (ISNA(kappa) || ISNAN(kappa) || kappa < BesselUpperB[i]);
  if (lambda>0) {
    if ((delta<0) || (kappa<=0)) {
      sprintf(msg, "kappa>0 and delta>=0 if lambda>0. Got kappa=%f and delta=%f for lambda=%f", kappa, delta, lambda);
      ERR(msg);
    }
  } else if (lambda<0) {
    if ((delta<=0) || (kappa<0)) {
      sprintf(msg, "kappa>=0 and delta>0 if lambda<0. Got kappa=%f and delta=%f for lambda=%f", kappa, delta, lambda);
      ERR(msg);
    }
  } else { // lambda==0.0
    if ((delta<=0) || (kappa<=0)) {
      sprintf(msg, "kappa>0 and delta>0 if lambda=0. Got kappa=%f and delta=%f for lambda=%f", kappa, delta, lambda);
      ERR(msg);
    }
  }
  return NOERROR;
}
void rangehyperbolic(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 0.000001;
  range->pmax[0] = 10.0;
  range->openmin[0] = false;
  range->openmax[0] = true;

  range->min[1] = RF_NEGINF;
  range->max[1] = RF_INF;
  range->pmin[1] = -20.0;
  range->pmax[1] = 20.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
 
  range->min[2] = 0.0;
  range->max[2] = RF_INF;
  range->pmin[2] = 0.000001;
  range->pmax[2] = 10.0;
  range->openmin[2] = false;
  range->openmax[2] = true;
}



/* iaco cesare model */
void IacoCesare(double *x, cov_model *cov, double *v){
    double kappa = cov->p[0][0], lambda=cov->p[1][0], delta=cov->p[2][0];
  *v = pow(1.0 + pow(x[0], kappa) + pow(x[1], lambda), - delta); 
}
//double ScaleIacoCesare(cov_model *cov, int scaling) { return 1.0; } 
void rangeIacoCesareX(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 1.0;
  range->max[0] = 2.0;
  range->pmin[0] = 1.0;
  range->pmax[0] = 2.0;
  range->openmin[0] = false;
  range->openmax[0] = false;

  range->min[1] = 1.0;
  range->max[1] = 2.0;
  range->pmin[1] = 1.0;
  range->pmax[1] = 2.0;
  range->openmin[1] = false;
  range->openmax[1] = false;
 
  range->min[2] = 0.5 * (double) cov->tsdim;
  range->max[2] = RF_INF;
  range->pmin[2] = range->min[2];
  range->pmax[2] = range->pmin[2] + 10.0;
  range->openmin[2] = false;
  range->openmax[2] = true;
}
void rangeIacoCesare(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = 0.0;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = 2.0;
  range->pmin[1] = 0.0;
  range->pmax[1] = 2.0;
  range->openmin[1] = true;
  range->openmax[1] = false;
 
  range->min[2] = 0.0;
  range->max[2] = RF_INF;
  range->pmin[2] = 0.0;
  range->pmax[2] = 10.0;
  range->openmin[2] = true;
  range->openmax[2] = true;
}
//int checkIacoCesare(cov_model *cov) { }


/* local-global distinguisher */
void lgd1(double *x, cov_model *cov, double *v) {
  double y = *x, kappa = cov->p[0][0], delta=cov->p[1][0];
  *v = (y == 0.0) ? 1.0 : (y < 1.0) 
    ? 1.0 - delta / (kappa + delta) * exp(kappa * log(y))
    : kappa / (kappa + delta) * exp( -delta * log(y));
}
double Scalelgd1(cov_model *cov) {
  double kappa = cov->p[0][0], delta=cov->p[1][0];
  return (19 * kappa < delta)
     ? exp(log(0.95 * (kappa + delta) / delta) / kappa)
     : exp(log(0.05 * (kappa + delta) / kappa) / delta);
}
void Dlgd1(double *x, cov_model *cov, double *v){
  double y=*x, pp, kappa = cov->p[0][0], delta=cov->p[1][0];
  if (y == 0.0) {
    *v = 0.0;// falscher Wert, aber sonst NAN-Fehler#
    return;
  }
  pp = ( (y < 1.0) ? kappa : -delta ) - 1.0;
  *v = - kappa * delta / (kappa + delta) * exp(pp * y);
}
void rangelgd1(cov_model *cov, range_arraytype* ra) {
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = (cov->tsdim==2) ? 0.5 : 1.0;
  range->pmin[0] = 0.01;
  range->pmax[0] = range->max[0];
  range->openmin[0] = true;
  range->openmax[0] = false;

  range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = 0.01;
  range->pmax[1] = 20.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
 
  double dim;
  dim = 2.0 * (1.5 - cov->p[0][0]);
  range->maxdim = dim >= INFDIM ? INFDIM-1 : (int) dim;
}


/* mastein */
// see Hypermodel.cc



/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */ 

/* Whittle-Matern or Whittle or Besset */ 
void WM(double *x, double nu, double *v) 
// check calling functions, like hyperbolic and gneiting if any changings !!
{
  static double nuOld=RF_INF;
  static double loggamma;
  double y = *x;
  
  if (y<=LOW_BESSELK) {
    *v = 1.0;
  } else {
    if (nu != nuOld) {
      nuOld = nu;
      loggamma = lgammafn(nu);
    }
    *v=2.0 * exp(nu * log(0.5 * y) - loggamma + log(bessel_k(y, nu, 2.0)) - y);
//     printf("whittle %e %f %f; %f %f %f; %f \n", y, nu, *v,
//	    nu * log(0.5 * y), - loggamma, log(bessel_k(y, nu, 2.0)),
//	    nu * log(0.5 * y) - loggamma + log(bessel_k(y, nu, 2.0)) - y
//	    );
  }
}


void DWM(double *x, double nu, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double nuOld=RF_INF;
  static double loggamma;
  double y = *x;
  if (y<=LOW_BESSELK) {
    *v = (nu > 0.5) ? 0.0 : (nu < 0.5) ? INFTY : 1.253314137;
  } else {
    if (nu!=nuOld) {
      nuOld = nu;
      loggamma = lgammafn(nu);
    }
    *v = -2.0 * exp(nu * log(0.5 * y) - loggamma + 
		    log(bessel_k(y, nu - 1.0, 2.0)) - y);
  }
//     printf("d1 %f %f %f\n", x[0], *v, INFTY);
}

void DDWM(double *x, double nu, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double nuOld=RF_INF;
  static double gamma;
  double y=*x;

  if (y<=LOW_BESSELK) {
    *v = (nu > 1.0) ? -0.5 / (nu - 1.0) : INFTY;
  } else {
    if (nu!=nuOld) {
      nuOld=nu;
      gamma = gammafn(nu);
    }
    *v = pow(0.5 * y , nu - 1.0) / gamma *
      ( - bessel_k(y, nu - 1.0, 1.0) + y * bessel_k(y, nu - 2.0, 1.0)); 
  }
//  printf("d2 %f %f %f gam=%f y=%f\n", x[0], *v, nu, gamma, y);
}


void D3WM(double *x, double nu, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double nuOld=RF_INF;
  static double gamma;
  double y=*x;

  if (y<=LOW_BESSELK) {
    *v = 0.0;
  } else {
    if (nu!=nuOld) {
      nuOld=nu;
      gamma = gammafn(nu);
    }
    *v = pow(0.5 * y , nu - 1.0) / gamma *
      ( 3.0 * bessel_k(y, nu - 2.0, 1.0) -y * bessel_k(y, nu - 3.0, 1.0)); 
  }
  //printf("d2 %f %f \n", x[0], *v);
}


void D4WM(double *x, double nu, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  static double nuOld=RF_INF;
  static double gamma;
  double y=*x;

  if (y<=LOW_BESSELK) {
      *v = (nu > 2.0) ? 0.75 / ((nu - 1.0) * (nu - 2.0)) : INFTY;
//      printf("LOW! %f %f\n", y, *v);      
  } else {
    if (nu != nuOld) {
      nuOld=nu;
      gamma = gammafn(nu);
    }
    *v = 0.25 * pow(0.5 * y , nu - 3.0) / gamma *
	(+ 6.0 * (nu - 3.0 - y * y) * bessel_k(y, nu - 3.0, 1.0)
	 + y * (3.0  + y * y) * bessel_k(y, nu - 4.0, 1.0)); 
    //  *v = 0.5 * pow(0.5 * y , nu - 2.0) / gamma *
    //   (   3.0 * bessel_k(y, nu - 2.0, 1.0)
//	 - 6.0 * y * bessel_k(y, nu - 3.0, 1.0)
//	 - y * y * bessel_k(y, nu - 4.0, 1.0)); 
  }
  //printf("d2 %f %f \n", x[0], *v);
}

double ScaleWM(cov_model *cov, double nu){
  // it is working reasonably well if nu is in [0.001,100]
  // happy to get any better suggestion!!

//    PrintModelInfo(cov);
    //return 1.0;
//    printf("%ld %d\n", cov, cov->p[0]==NULL); //assert(false);
    

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
  
  return interpolate(log(nu) * INVLOG2, stuetz, nstuetz, stuetzorigin, 
		     1.5, 5);
}

double densityWM(double *x, int dim, double nu) {
  double x2;
  int d;
  x2 = x[0] * x[0];
  for (d=1; d<dim; d++) x2 += x[d] * x[d];
  x2 += 1.0;
  
  return exp(lgammafn(nu + 0.5 * (double) dim)
	     - lgammafn(nu)
	     - (double) dim * M_LN_SQRT_PI
	     - (nu + 0.5 * (double) dim) * log(x2));
}

int checkWM(cov_model *cov, double nu) { 
  static double spectrallimit=0.17;
  int i;
  for (i=0; i<= Nothing; i++)
    cov->pref[i] *= (ISNA(nu) || ISNAN(nu) || nu < BesselUpperB[i]);
  if (nu != 0.5) cov->pref[TBM2] = PREF_NONE;
  if (nu < spectrallimit) cov->pref[SpectralTBM] = PREF_NONE;
  if (cov->tsdim > 2) 
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (nu > 2.5)  cov->pref[CircEmbed] = 2;
  if (nu <= 1.0) cov->derivatives = 0; else if (nu<=2.0) cov->derivatives = 1;
  return NOERROR;
}
void ieinitWM(double nu, localinfotype *li) {
  // intrinsic_rawr
  li->instances = 1;
  if (nu <= 0.5) {
    li->value[0] = 1.0;
    li->msg[0] = MSGLOCAL_OK;
  } else {
    li->value[0] = 1.5;
    li->msg[0] = MSGLOCAL_JUSTTRY;
  }
}

void coinitWM(double nu, localinfotype *li) {
  // cutoff_A
  double thres[2] = {0.25, 0.5}; 
  if (nu <= thres[0]) {
    li->instances = 2;
    li->value[0] = 0.5;
    li->value[1] = 1.0;
    li->msg[0] = li->msg[1] = MSGLOCAL_OK;
  } else {
    li->instances = 1;
    li->value[0] = 1.0; //  q[CUTOFF_A] 
    li->msg[0] = (nu <= thres[1]) ? MSGLOCAL_OK : MSGLOCAL_JUSTTRY;
  } 
}
double DrawLogMixWM(cov_model *cov, mpp_storage *s) { // inv scale
  // V ~ F in PSgen
  return log(-0.25 / log(UNIFORM_RANDOM));
}
double LogMixWeightWM(double nu, mpp_storage *s) {
  // g(v,x) in PSgen
  return (M_LN2  + 0.5 * s->loginvscale) * (1.0 - nu) - 0.5 * lgammafn(nu) ; 
}



/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */ 
#define MATERN_NU_THRES 100
void Matern(double *x, cov_model *cov, double *v) 
// check calling functions, like hyperbolic and gneiting if any changings !!
{
    double y,
	nu = ((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0],
	nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES;

     y = *x * sqrt(2.0 * nuThres);
     WM(&y, nuThres, v);
	
     //   printf("%f %f %f %f\n", nu, nuThres, *x, *v);
    
    if (nu > MATERN_NU_THRES) {
	double w;
	y = *x * INVSQRTTWO;
	Gauss(&y, cov, &w);
	*v = *v * MATERN_NU_THRES / nu + (1 - MATERN_NU_THRES / nu) * w;
    }
//  printf("%f %f %f\n", *x, y, *v);
}

double Scalematern(cov_model *cov) {
    return ScaleWM(cov, ((int*) (cov->p[1]))[0] == 0 
		   ? cov->p[0][0] : 1.0 / cov->p[0][0]);
}

void DMatern(double *x, cov_model *cov, double *v) { 
    double y, scale,
	nu = ((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0],
	nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES;
 
    scale = sqrt(2.0 * nuThres);
    y = *x * scale;
    DWM(&y, nuThres, v);
    *v *= scale;

    if (nu > MATERN_NU_THRES) {
	double w;
	scale = INVSQRTTWO;
	y = *x * scale;
	DGauss(&y, cov, &w);
	w *= scale;
	*v = *v * MATERN_NU_THRES / nu + (1 - MATERN_NU_THRES / nu) * w;
    }
}

void DDMatern(double *x, cov_model *cov, double *v) { 
    double y, scaleSq, scale, 
	nu = ((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0],
	nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES;

    scaleSq = 2.0 * nuThres;
    scale = sqrt(scaleSq);
    y = *x * scale;
    DDWM(&y, nuThres, v);
    *v *= scaleSq;

    if (nu > MATERN_NU_THRES) {
       double w;
       scaleSq = 0.5;
       scale = INVSQRTTWO;
       y = *x * scale;
       DDGauss(&y, cov, &w);
       w *= scaleSq;
       *v = *v * MATERN_NU_THRES / nu + (1 - MATERN_NU_THRES / nu) * w;
  }
}
void D3Matern(double *x, cov_model *cov, double *v) { 
     double y, scaleSq, scale, 
	nu = ((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0],
	nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES;
 
     scaleSq = 2.0 * nuThres;
     scale = sqrt(scaleSq);
     y = *x * scale;
     D3WM(&y, nuThres, v);
     *v *= scaleSq * scale;
  
     if (nu > MATERN_NU_THRES) {
       double w;
       scaleSq = 0.5;
       scale = INVSQRTTWO;
       y = *x * scale;
       D3Gauss(&y, cov, &w);
       w *= scaleSq * scale;
       *v = *v * MATERN_NU_THRES / nu + (1 - MATERN_NU_THRES / nu) * w;
     }
  // printf("ZM %f %f  ", *x, *v);
}
void D4Matern(double *x, cov_model *cov, double *v) { 
    double y, scaleSq, scale, 
	nu = ((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0],
	nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES;
 
    scaleSq = 2.0 * nuThres;
    scale = sqrt(scaleSq);
    y = *x * scale;
    D4WM(&y, nuThres, v);
    *v *= scaleSq * scaleSq;

    if (nu > MATERN_NU_THRES) {
      double w;
      scaleSq = 0.5;
      scale = INVSQRTTWO;
      y = *x * scale;
      D4Gauss(&y, cov, &w);
      w *= scaleSq * scaleSq;
      *v = *v * MATERN_NU_THRES / nu + (1 - MATERN_NU_THRES / nu) * w;
  }
  // printf("ZM %f %f  ", *x, *v);
}

int checkMatern(cov_model *cov) { 
    int i, err;
  kdefault(cov, 1, 0.0);
  if ((err = checkkappas(cov)) != NOERROR) return err;
  for (i=0; i<Forbidden; i++) {
      if (cov->user[i] < cov->pref[i])
	  cov->pref[i] = cov->user[i];
  }
  return checkWM(cov, ((int*) (cov->p[1]))[0] == 0 
		   ? cov->p[0][0] : 1.0 / cov->p[0][0]);
}
void rangeMatern(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-1;
  range->pmax[0] = 20.0;
  range->openmin[0] = false;
  range->openmax[0] = true;

  range->min[1] = 0.0;
  range->max[1] = 1.0;
  range->pmin[1] = 0.0;
  range->pmax[1] = 1.0;
  range->openmin[1] = false;
  range->openmax[1] = false;
}
void ieinitMatern(cov_model *cov, localinfotype *li) {
  ieinitWM(((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0], li);
}
void coinitMatern(cov_model *cov, localinfotype *li) {
  coinitWM(((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0], li);
}
void invMaternSq(double *x, cov_model *cov, double *v) 
// check calling functions, like hyperbolic and gneiting if any changings !!
{
//  double nu = cov->p[0][0];
  
  ERR("inv Matern not programmed yet");
//  printf("%f %f %f\n", *x, y, *v);
}

double densityMatern(double *x, cov_model *cov) {
  double nu = ((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0],
    y = *x / sqrt(2.0 * nu);
  if (nu > 50) error("nu>50 in density of matern class numerically instable");
  return pow(2.0 * nu, 0.5 * cov->tsdim) * densityWM(&y, cov->tsdim, nu);
}
int initspectralMatern(cov_model *cov) {
  spec_covstorage *cs = &(cov->spec);
  if (cov->tsdim <= 2) return NOERROR;
  cs->density = densityMatern;
  return search_metropolis(cov);
}

void spectralMatern(cov_model *cov, spectral_storage *s, double *e) { 
  /* see Yaglom ! */
  if (cov->tsdim <= 2) {
    double nu = ((int*) (cov->p[1]))[0] == 0 ? cov->p[0][0] : 1.0 / cov->p[0][0];
    E2(s, sqrt( 2.0 * nu * (pow(1.0 - UNIFORM_RANDOM, -1.0 / nu) - 1.0) ), e);
  } else {
    metropolis(cov, e);
  }
}



/* nugget effect model */
void nugget(double *x, cov_model *cov, double *v) {
  double diag = (*x <= GLOBAL.nugget.tol) ? 1.0 : 0.0;
  int i, endfor,
    vdim   = cov->vdim,
    vdimsq = vdim * vdim;
  
  v[0] = diag;
  for (i = 1; i<vdimsq; v[i++] = diag) {
    endfor = i + vdim;
    for (; i<endfor; v[i++] = 0.0);
  }
}
double Scalenugget(cov_model *cov) { return 1.0; }//or better 0.0 => error?
void rangenugget(cov_model *cov, range_arraytype* ra){  
  range_type *range = ra->ranges;
  range->min[0]  = 1.0;
  range->max[0]  = RF_INF;
  range->pmin[0] = 1.0;
  range->pmax[0] = 10.0;
  range->openmin[0] = false;
  range->openmax[0] = true;

  range->finiterange = true;
}
int checknugget(cov_model *cov) {
  if (cov->ncol[0]==0 && cov->nrow[0]==0) {
    cov->vdim = cov->calling->vdim;
    if (cov->vdim <= 0) cov->vdim = 1;
    cov->ncol[0] = cov->nrow[0] = 1;
    cov->p[0] = (double*) malloc(sizeof(int));
    ((int*) cov->p[0])[0] = cov->vdim;
  } else {
    if (cov->ncol[0]!=1 || cov->nrow[0]!=1) ERR("parameter must be a scalar");
    cov->vdim =((int*) cov->p[0])[0];
  }
  return NOERROR;
}



// Paciore und Stein 
// see NonIsoCovFct.cc


/* penta */
void penta(double *x, cov_model *cov, double *v)
{ ///
  double y=*x, y2 = y * y;
  *v = (y >= 1.0) ? 0.0 :
    (1.0 + y2 * (-7.333333333333333 
		 + y2 * (33.0 +
			 y * (-38.5 + 
			      y2 * (16.5 + 
				    y2 * (-5.5 + 
					  y2 * 0.833333333333333))))));
}
void Dpenta(double *x, cov_model *cov, double *v)
{ ///
  double y=*x, y2 = y * y;
  *v = (y >= 1.0) ? 0.0 :
    y * (-14.66666666666666667 + 
	 y2 * (132.0 + 
	       y * (-192.5 + 
		    y2 * (115.5 + 
			  y2 * (-49.5 + 
				y2 * 9.16666666666666667)))));
  
}
double Scalepenta(cov_model *cov) {return 1.6552838957365;}
void rangepenta(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->maxdim = 3; 
  range->finiterange = true;
}


/* power model */ 
void power(double *x, cov_model *cov, double *v){
  double kappa = cov->p[0][0], y = *x;
  *v = (y >= 1.0) ? 0.0 : pow(1.0 - y, kappa);
}
//double Scalepower(cov_model *cov){ 
//   double kappa = cov->p[0][0];
//   return 1.0 / (1.0 - pow(0.05, 1.0 / kappa));
//}
void TBM2power(double *x, cov_model *cov, double *v){
  // only kappa=2 up to now !
  double y = *x;
  *v = (y > 1.0)  
    ? (1.0 - 2.0 * y *(asin(1.0 / y) - y + sqrt(y * y - 1.0) ))
    : 1.0 - y * (PI - 2.0 * y);
}
void Dpower(double *x, cov_model *cov, double *v){
  double kappa = cov->p[0][0], y = *x;
  *v = (y >= 1.0) ? 0.0 : - kappa * pow(1.0 - y, kappa - 1.0);
}
int checkpower(cov_model *cov) {
  if (cov->p[0][0] != 2.0) cov->pref[TBM2] = PREF_NONE;
  return NOERROR;
}
// range definition:
// 0: min, theory, 1:max, theory
// 2: min, practically 3:max, practically
void rangepower(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.5 * (double) (cov->tsdim + 1);
  range->max[0] = RF_INF;
  range->pmin[0] = range->min[0];
  range->pmax[0] = 20.0;
  range->openmin[0] = false;
  range->openmax[0] = true;

  double dim = 2.0 * cov->p[0][0] - 1.0;
  range->maxdim = dim >= INFDIM ? INFDIM-1 : (int) dim;
  range->finiterange = true;
}


/* qexponential -- derivative of exponential */
void qexponential(double *x, cov_model *cov, double *v){
  double kappa = cov->p[0][0], y = exp(-*x);
  *v = y * (2.0  - kappa * y) / (2.0 - kappa);
}
double Scaleqexponential(cov_model *cov){
  double kappa = cov->p[0][0];
  return -1.0 / log( (1.0 - sqrt(1.0 - kappa * (2.0 - kappa) * 0.05)) / kappa);
} 
void Dqexponential(double *x, cov_model *cov, double *v) {
  double kappa = cov->p[0][0], y = exp(-*x);
  *v = y * (kappa * y - 1.0) * 2.0 / (2.0 - kappa);
}
void rangeqexponential(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 1.0;
  range->pmin[0] = 0.0;
  range->pmax[0] = 1.0;
  range->openmin[0] = false;
  range->openmax[0] = false;
}
//void infoqexponential(cov_model *cov, info_type* info) {
//  info->maxdim = (kappa>=0.0 && kappa<=1.0) ? INFDIM : 0; 
//}


/* spherical model */ 
void spherical(double *x, cov_model *cov, double *v){
  double y = *x;
  *v = (y >= 1.0) ? 0.0 : 1.0 + y * 0.5 * (y * y - 3.0);
//  printf("%f %f\n", y, *v);
}
// double Scalespherical(cov_model *cov){ return 1.23243208931941;}
void TBM2spherical(double *x, cov_model *cov, double *v){
  double y = *x , y2 = y * y;
  *v = (y>1.0) 
      ? (1.0- 0.75 * y * ((2.0 - y2) * asin(1.0/y) + sqrt(y2 -1.0)))
      : (1.0 - 0.375 * PI * y * (2.0 - y2));
}
void Dspherical(double *x, cov_model *cov, double *v){
  double y = *x;
  *v = (y >= 1.0) ? 0.0 : 1.5 * (y * y - 1.0);
}
void rangespherical(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->maxdim = 3; 
  range->finiterange = true;
}



/* stable model */
void stable(double *x, cov_model *cov, double *v){
  double y = *x, kappa = cov->p[0][0];  
  *v = (y==0.0) ? 1.0 : exp(-pow(y, kappa));
}
double Scalestable(cov_model *cov){
  double kappa = cov->p[0][0];  
  return pow(MINUSINVLOG005, 1.0 / kappa);
}
void Dstable(double *x, cov_model *cov, double *v){
  double z, y = *x, kappa = cov->p[0][0];
  if (y == 0.0) {
    *v = (kappa > 1.0) ? 0.0 : (kappa < 1.0) ? INFTY : 1.0;
  } else {
    z = pow(y, kappa - 1.0);
    *v = -kappa * z * exp(-z * y);
  }
}
/* stable: second derivative at t=1 */
void DDstable(double *x, cov_model *cov, double *v) 
{
  double z, y = *x, kappa = cov->p[0][0], xkappa;
  if (y == 0.0) {
      *v = (kappa == 1.0) ? 1.0 : (kappa == 2.0) ? kappa * (1 - kappa) : INFTY;
  } else {
    z = pow(y, kappa - 2.0);
    xkappa = z * y * y;
    *v = kappa * (1.0 - kappa + kappa * xkappa) * z * exp(-xkappa);
  }
}
void invstableSq(double *x, cov_model *cov, double *v){
  double y = *x, kappa = cov->p[0][0];  
  if (y == 0.0) *v = RF_INF;
  else {
    *v = pow( - log(y), 2.0 / kappa);
  }
}

int checkstable(cov_model *cov) {
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (cov->p[0][0] == 2.0)
    cov->pref[CircEmbed] = 2;
  return NOERROR;
}
  
void rangestable(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = 0.06;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;
}
void coinitstable(cov_model *cov, localinfotype *li) {
  coinitgenCauchy(cov, li);
}
void ieinitstable(cov_model *cov, localinfotype *li) {  
  ieinitgenCauchy(cov, li);
}




/* SPACEISOTROPIC stable model for testing purposes only */
void stableX(double *x, cov_model *cov, double *v){
  double y, kappa = cov->p[0][0];
  y = x[0] * x[0] + x[1] * x[1];
  *v = (y==0.0) ? 1.0 : exp(-pow(y, 0.5 * kappa));
}
void DstableX(double *x, cov_model *cov, double *v){
  double z, y, kappa = cov->p[0][0];
  y = x[0] * x[0] + x[1] * x[1];
  if (y == 0.0) {
    *v = ((kappa > 1.0) ? 0.0 : (kappa < 1.0) ? INFTY : 1.0);
  } else {
    z = pow(y, 0.5 * kappa - 1.0);
    *v = -kappa * x[0] * z * exp(- z * y);
  }
}
/* END SPACEISOTROPIC stable model for testing purposes only */


/* Stein */
// see Hypermodel.cc

/* stein space-time model */
void kappaSteinST1(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = (i == 0) ? 1 : (i==1) ?  cov->tsdim - 1 : -1;
}
void SteinST1(double *x, cov_model *cov, double *v){
/* 2^(1-nu)/Gamma(nu) [ h^nu K_nu(h) - 2 * tau (x T z) t h^{nu-1} K_{nu-1}(h) /
   (2 nu + d + 1) ]

   \|tau z \|<=1 hence \tau z replaced by z !!
*/
  int d,
    dim = cov->tsdim,
    time = dim -1;
  double logconst, y, s,
    nu = cov->p[0][0],
    *z=cov->p[1];
  
  static double nuold=RF_INF;
  static double loggamma;
  static int dimold;

  if (nu != nuold || dimold != dim) {
    nuold = nu;
    dimold = dim;
    loggamma = lgammafn(nu);
  }

  y = 0.0;
  s = x[time] * x[time];
  for (d=0; d<time; d++) {
    s += x[d] * x[d];
    y += x[d] * z[d];
  }

  if ( s==0.0 ) *v = 1.0;
  else {
    s = sqrt(s);
    logconst = (nu - 1.0) * log(0.5 * s)  - loggamma;
    *v =  s * exp(logconst + log(bessel_k(s, nu, 2.0)) - s)
      - 2.0 * y * x[time] * exp(logconst + log(bessel_k(s, nu - 1.0, 2.0)) -s) 
      / (2.0 * nu + time);
  }
}
double densitySteinST1(double *x, cov_model *cov) {
  double x2, wz, dWM, nu = cov->p[0][0], 
    *z=cov->p[1];
  int d,
    dim = cov->tsdim,
    spatialdim = dim - 1;
  static double nuold = RF_INF;
  static int dimold = -1;
  static double constant;
  static double factor;
  if (nu != nuold || dimold != dim) {
    nuold = nu;
    dimold = dim;
    constant = lgammafn(nu) - lgammafn(nu +  0.5 * dim) - dim * M_LN_SQRT_PI;
    factor = nu + dim;
  }

  x2 = x[spatialdim] * x[spatialdim]; // v^2
  wz = 0.0;
  for (d=0; d<spatialdim; d++) {
    x2 += x[d] * x[d];  // w^2 + v^2
    wz += x[d] * z[d];
  }
  dWM = exp(constant - factor * log(x2 + 1.0));
  return (1.0 + 2.0 * wz * x[spatialdim] + x2) * dWM
    / (2.0 * nu + (double) dim + 1.0);
}
int initspectralSteinST1(cov_model *cov) {
  spec_covstorage *cs = &(cov->spec);
  cs->density = densitySteinST1;
  return search_metropolis(cov);
//  return (cov->tsdim == 3) ? NOERROR : ERRORFAILED;
}

void spectralSteinST1(cov_model *cov, spectral_storage *s, double *e){
  metropolis(cov, e);
}

void rangeSteinST1(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;   
  
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-10;
  range->pmax[0] = 10.0;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = RF_NEGINF;
  range->max[1] = RF_INF;
  range->pmin[1] = -10.0;
  range->pmax[1] = 10.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
}
int checkSteinST1(cov_model *cov) {  
  double nu = cov->p[0][0], *z= cov->p[1], absz;
  int d, dim=cov->tsdim-1;

  for (d=0; d<= Nothing; d++) cov->pref[d] *= (nu < BesselUpperB[d]);
  if (nu >= 2.5) cov->pref[CircEmbed] = 2;
  if (dim <= 1) {
    sprintf(ERRORSTRING, "%s", "dim of coordinates must be at least 2");
    return ERRORM;
  }
  	
  for (absz=0.0, d=0;  d<dim; d++)  absz += z[d] * z[d];
  if (ISNA(absz) || ISNAN(absz))
    ERR("currently, components of z cannot be estimated by MLE, so NA's are not allowed");
  if (absz > 1.0 + UNIT_EPSILON && !GLOBAL.general.skipchecks) {
    ERR("||z|| must be less than or equal to 1");
  }
  return NOERROR;
}




/* wave */
void wave(double *x, cov_model *cov, double *v) {
  double y = *x;
  *v = (y==0.0) ? 1.0 : sin(y) / y;
}
double Scalewave(cov_model *cov) {return 0.302320850755833;}
int initspectralwave(cov_model *cov) {
  return (cov->tsdim <= 2) ? NOERROR : ERRORFAILED;
}
void spectralwave(cov_model *cov, spectral_storage *s, double *e) { 
  /* see Yaglom ! */
  double x;  
  x = UNIFORM_RANDOM; 
  E2(s, sqrt(1.0 - x * x), e);
}
void rangewave(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->maxdim = 3;
}





////////////////////////////////////
// Whittlemodel

void Whittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
    WM(x, cov->p[0][0], v);
}

void DWhittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
    DWM(x, cov->p[0][0], v);
}

void DDWhittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
    DDWM(x, cov->p[0][0], v);
}


void D3Whittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
    D3WM(x, cov->p[0][0], v);
}

void D4Whittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
    D4WM(x, cov->p[0][0], v);
}

double ScaleWhittle(cov_model *cov){
    return ScaleWM(cov, cov->p[0][0]);
}

void TBM2Whittle(double *x, cov_model *cov, double *v) 
{
  double nu = cov->p[0][0];
  if (nu == 0.5) TBM2exponential(x, cov, v);
  else assert(false);
}


void invWhittleSq(double *x, cov_model *cov, double *v) 
// check calling functions, like hyperbolic and gneiting if any changings !!
{
//  double nu = cov->p[0][0];
  ERR("inv Matern not programmed yet");
//  printf("%f %f %f\n", *x, y, *v);
}


double densityWhittle(double *x, cov_model *cov) {
   return densityWM(x, cov->tsdim, cov->p[0][0]);
}
int initspectralWhittle(cov_model *cov) {
  spec_covstorage *cs = &(cov->spec);
  if (cov->tsdim <= 2) return NOERROR;
  cs->density = densityWhittle; 
  int err=search_metropolis(cov);
  return err;
}

void spectralWhittle(cov_model *cov, spectral_storage *s, double *e) { 
  /* see Yaglom ! */
  if (cov->tsdim <= 2) {
    double nu = cov->p[0][0];
    E2(s, sqrt(pow(1.0 - UNIFORM_RANDOM, -1.0 / nu) - 1.0), e);
  } else {
    metropolis(cov, e);
  }
}

int checkWhittle(cov_model *cov) { 
  return checkWM(cov, cov->p[0][0]);
}
void rangeWhittle(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-2;
  range->pmax[0] = 10.0;
  range->openmin[0] = true;
  range->openmax[0] = true;
}
void ieinitWhittle(cov_model *cov, localinfotype *li) {
    ieinitWM(cov->p[0][0], li);
}
void coinitWhittle(cov_model *cov, localinfotype *li) {
    coinitWM(cov->p[0][0], li);
}
double LogMixWeightW(double *x, cov_model *cov, mpp_storage *s) {
  // g(v,x) in PSgen
    return LogMixWeightWM(cov->p[0][0], s);  
}








/* Whittle-Matern or Whittle or Besset */ 
/// only lower parts of the matrices, including the diagonals are used when estimating !!

void kappa_biWM(int i, cov_model *cov, int *nr, int *nc){
    *nc = *nr = 1;
    if (i==0 || i==2 || i==4) *nr=2;
}

//
//int zaehler = 0;

/* Whittle-Matern or Whittle or Besset */ 
void biWM2(double *x, cov_model *cov, double *v) {

  double factor, beta, gamma, tsq, t1sq, t2sq, inf, infQ, discr,
      lg[4], z, c[4], nu[4], alpha , a[4], a0sq, a3sq, a1sq,
      dim = (double) cov->tsdim, 
      d2 = dim * 0.5, 
      y = *x,
      *tildenu = cov->p[0],
      tildenu12 = cov->p[1][0],
      *tildea = cov->p[2],
      tildea12 = cov->p[3][0],
      *tildec = cov->p[4],
      rho = cov->p[5][0];
  int i;

  a[0] = 1.0 / tildea[0];
  a[3] = 1.0 / tildea[1];
  a[1] = a[2] = 1.0 / tildea12;
  a0sq = a[0] * a[0];
  a3sq = a[3] * a[3];
  a1sq = a[1] * a[1];

  nu[0] = tildenu[0];
  nu[3] = tildenu[1];
  nu[2] = nu[1] = 0.5 * (nu[0] + nu[3]) / tildenu12;

  lg[0] = lgammafn(nu[0]);
  lg[2] = lg[1] = lgammafn(nu[1]);
  lg[3] = lgammafn(nu[3]);

  alpha = 2 * nu[1] - nu[0] - nu[3];

//  printf("%f %f %f %f\n", tildenu[0],  tildenu[1],  tildenu[2],  tildenu[3]);
  // print("%f\n", alpha); //assert(false);

  // **************** ACHTUNG !!!!! *********
  // nicht gut, sollte in check einmal berechnet werden 
  // dies wiederspricht aber der MLE Maximierung, da dann
  // neu berechnet werden muss; verlg. natsc und cutoff (wo es nicht
  // programmiert ist !!)
  
 
  factor = exp(lgammafn(nu[0] + d2) - lg[0]
	       + lgammafn(nu[3] + d2) - lg[3]
		   + 2.0 * (lg[1]  - lgammafn(nu[1] + d2)
			    + nu[0] * log(a[0]) + nu[3] * log(a[3])
			    - 2.0 * nu[1] * log(a[1]))
	);

//  printf("factor dim %d\n", cov->tsdim);

  // alpha u^2 + beta u + gamma = 0

  gamma = (2.0 * nu[1] + dim) * a0sq * a3sq 
    - (nu[3] + d2) * a0sq * a1sq
    - (nu[0] + d2) * a3sq * a1sq;
      
  beta = (2.0 * nu[1] - nu[3] + d2) * a0sq
      + (2.0 * nu[1] - nu[0] + d2) * a3sq
      - (nu[0] + nu[3] + dim) * a1sq;

//  printf("\nalpha=%f beta=%f gamma=%f\n", alpha, beta, gamma);
     
  if (tildenu12 == 1.0) { // lin. eqn only
    t2sq = (beta == 0.0) ? 0.0 : -gamma / beta;
    if (t2sq < 0.0) t2sq = 0.0;
    t1sq =  t2sq;
  } else { // quadratic eqn.
    discr = beta * beta - 4.0 * alpha * gamma;
    if (discr < 0.0) {
      t1sq = t2sq = 0.0;
    } else {
      discr = sqrt(discr);
      t1sq = (-beta + discr) / (2.0 * alpha);
      if (t1sq < 0.0) t1sq = 0.0;
      t2sq = (-beta - discr) / (2.0 * alpha);
      if (t2sq < 0.0) t2sq = 0.0;	  
    }
  }


  inf = tildenu12 == 1.0 ? 1.0 : RF_INF; // t^2 = infty  ; tildenu[1]>1.0 not
  //                                   allowed by definition
  for (i=0; i<3; i++) {
    tsq = (i == 0) ? 0.0 
      : (i == 1) ? t1sq
      : t2sq;
    
    //      printf("discr=%f a1sq=%f tsq=%f nu1=%f \ndim=%f a0sq=%f nu0=%f a3sq=%f nu3=%f d/2=%f\n",discr,
//	   a1sq, tsq, nu[1], dim, a0sq,  nu[0], a3sq , nu[3], d2);
    
    infQ = pow(a1sq + tsq, 2.0 * nu[1] + dim) /
	  (pow(a0sq + tsq, nu[0] + d2) * pow(a3sq + tsq, nu[3] + d2));

//       printf("i=%d, tsq=%f %f %f %f %e %f\n", 
//	      i, tsq, infQ, discr, beta, alpha, -gamma / beta);
    
 
    if (infQ < inf) inf = infQ;
  }

  c[0] = tildec[0];
  c[3] = tildec[1];
  c[2] = c[1] = rho * sqrt(c[0] *  c[3] * factor * inf);

 
  // print("c "); {int j; for (j=0; j<4; j++) print("%f ", c[j]); printf("\n");}

  //   assert(zaehler++ < 1);

//  printf("%f %f rho=%f factor=%f inf=%f prod=%f\n", c[0], c[3], rho, factor, inf,
//	 rho*sqrt(factor * inf));
    
  for (i=0; i<4; i++) {
    if (i==2) {
      v[2] = v[1];
      continue;
    }
    z = fabs(a[i] * y);
    v[i] = c[i] * ((z <= LOW_BESSELK) 
		   ? 1.0
		   : 2.0 * exp(nu[i] * log(0.5 * z) - lg[i] + 
			       log(bessel_k(z, nu[i], 2.0)) - z)
		   );
    // printf("ayz=%f %f %f\n", a[i], y, z);
  }

  //assert(false);
}



int checkbiWM2(cov_model *cov) { 
  cov->vdim = 2;
  return NOERROR;
}
  

  
void rangebiWM2(cov_model *cov, range_arraytype* ra){
//  kappanames("nu", REALSXP, "a", REALSXP, "c", REALSXP);
  range_type *range = ra->ranges;
//      printf("ij=%d \n", cov->vdim);

  // *nu = cov->p[0], 
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-5;
  range->pmax[0] = 10.0;
  range->openmin[0] = true;
  range->openmax[0] = true;
  
  // *nured12 = cov->p[1], 
  range->min[1] = 0;
  range->max[1] = 1;
  range->pmin[1] = 0.01;
  range->pmax[1] = 1;
  range->openmin[1] = true;
  range->openmax[1] = false;
    
  //   *a = cov->p[2], 
  range->min[2] = 0.0;
  range->max[2] = RF_INF;
  range->pmin[2] = 1e-5;
  range->pmax[2] = 100.0;
  range->openmin[2] = true;
  range->openmax[2] = true;
 
 //   *a12 = cov->p[3], 
  range->min[3] = 0.0;
  range->max[3] = RF_INF;
  range->pmin[3] = 1e-5;
  range->pmax[3] = 100.0;
  range->openmin[3] = true;
  range->openmax[3] = true;
  
  //    *c = cov->p[4]; 
  range->min[4] = 0;
  range->max[4] = RF_INF;
  range->pmin[4] = 0;
  range->pmax[4] = 1000;
  range->openmin[4] = false;
  range->openmax[4] = true;  
  
 //    *rho = cov->p[5]; 
  range->min[5] = -1;
  range->max[5] = 1;
  range->pmin[5] = -1;
  range->pmax[5] = 1;
  range->openmin[5] = false;
  range->openmax[5] = false;    


}






//////////////////////////////////////////////////////////////////////


void kappa_parsbiWM(int i, cov_model *cov, int *nr, int *nc){
    *nc = *nr = 1;
    if (i==0 || i==2) *nr=2;
}

/* Whittle-Matern or Whittle or Besset */ 
void parsbiWM2(double *x, cov_model *cov, double *v) {

  double factor, beta, gamma, tsq, t1sq, t2sq, inf, infQ, discr,
      lg[4], z, c[4], nu[4], alpha , a[4], a0sq, a3sq, a1sq, tildenu12,
      dim = (double) cov->tsdim, 
      d2 = dim * 0.5, 
      y = *x,
      *tildenu = cov->p[0],   
      tildea = cov->p[1][0],
      *tildec = cov->p[2],
      rho = cov->p[3][0];
  int i;

  a[0] = a[3] =  a[1] = a[2] = 1.0 / tildea; 
  a0sq = a[0] * a[0];
  a3sq = a[3] * a[3];
  a1sq = a[1] * a[1];

  tildenu12 =  1.0;
  nu[0] = tildenu[0];
  nu[3] = tildenu[1]; 
  nu[2] = nu[1] = 0.5 * (nu[0] + nu[3]) / tildenu12;

  lg[0] = lgammafn(nu[0]);
  lg[2] = lg[1] = lgammafn(nu[1]);
  lg[3] = lgammafn(nu[3]);

  alpha = 2 * nu[1] - nu[0] - nu[3];

//  printf("%f %f %f %f\n", tildenu[0],  tildenu[1],  tildenu[2],  tildenu[3]);
  // print("%f\n", alpha); //assert(false);

  // **************** ACHTUNG !!!!! *********
  // nicht gut, sollte in check einmal berechnet werden 
  // dies wiederspricht aber der MLE Maximierung, da dann
  // neu berechnet werden muss; verlg. natsc und cutoff (wo es nicht
  // programmiert ist !!)
  
 
  factor = exp(lgammafn(nu[0] + d2) - lg[0]
	       + lgammafn(nu[3] + d2) - lg[3]
		   + 2.0 * (lg[1]  - lgammafn(nu[1] + d2)
			    + nu[0] * log(a[0]) + nu[3] * log(a[3])
			    - 2.0 * nu[1] * log(a[1]))
	);

//  printf("factor dim %d\n", cov->tsdim);

  // alpha u^2 + beta u + gamma = 0

  gamma = (2.0 * nu[1] + dim) * a0sq * a3sq 
    - (nu[3] + d2) * a0sq * a1sq
    - (nu[0] + d2) * a3sq * a1sq;
      
  beta = (2.0 * nu[1] - nu[3] + d2) * a0sq
      + (2.0 * nu[1] - nu[0] + d2) * a3sq
      - (nu[0] + nu[3] + dim) * a1sq;

//  printf("\nalpha=%f beta=%f gamma=%f\n", alpha, beta, gamma);
     
  if (tildenu12 == 1.0) { // lin. eqn only
    t2sq = (beta == 0.0) ? 0.0 : -gamma / beta;
    if (t2sq < 0.0) t2sq = 0.0;
    t1sq =  t2sq;
  } else { // quadratic eqn.
    discr = beta * beta - 4.0 * alpha * gamma;
    if (discr < 0.0) {
      t1sq = t2sq = 0.0;
    } else {
      discr = sqrt(discr);
      t1sq = (-beta + discr) / (2.0 * alpha);
      if (t1sq < 0.0) t1sq = 0.0;
      t2sq = (-beta - discr) / (2.0 * alpha);
      if (t2sq < 0.0) t2sq = 0.0;	  
    }
  }


  inf = tildenu12 == 1.0 ? 1.0 : RF_INF; // t^2 = infty  ; tildenu[1]>1.0 not
  //                                   allowed by definition
  for (i=0; i<3; i++) {
    tsq = (i == 0) ? 0.0 
      : (i == 1) ? t1sq
      : t2sq;
    
    //      printf("discr=%f a1sq=%f tsq=%f nu1=%f \ndim=%f a0sq=%f nu0=%f a3sq=%f nu3=%f d/2=%f\n",discr,
//	   a1sq, tsq, nu[1], dim, a0sq,  nu[0], a3sq , nu[3], d2);
    
    infQ = pow(a1sq + tsq, 2.0 * nu[1] + dim) /
	  (pow(a0sq + tsq, nu[0] + d2) * pow(a3sq + tsq, nu[3] + d2));

//       printf("i=%d, tsq=%f %f %f %f %e %f\n", 
//	      i, tsq, infQ, discr, beta, alpha, -gamma / beta);
    
 
    if (infQ < inf) inf = infQ;
  }

  c[0] = tildec[0];
  c[3] = tildec[1];
  c[2] = c[1] = rho * sqrt(c[0] *  c[3] * factor * inf);

 
  // print("c "); {int j; for (j=0; j<4; j++) print("%f ", c[j]); printf("\n");}

  //   assert(zaehler++ < 1);

//  printf("%f %f rho=%f factor=%f inf=%f prod=%f\n", c[0], c[3], rho, factor, inf,
//	 rho*sqrt(factor * inf));
    
  for (i=0; i<4; i++) {
    if (i==2) {
      v[2] = v[1];
      continue;
    }
    z = fabs(a[i] * y);
    v[i] = c[i] * ((z <= LOW_BESSELK) 
		   ? 1.0
		   : 2.0 * exp(nu[i] * log(0.5 * z) - lg[i] + 
			       log(bessel_k(z, nu[i], 2.0)) - z)
		   );
    // printf("ayz=%f %f %f\n", a[i], y, z);
  }

  //assert(false);
}



int checkparsbiWM2(cov_model *cov) { 
  cov->vdim = 2;
  return NOERROR;
}
  

  
void rangeparsbiWM2(cov_model *cov, range_arraytype* ra){
//  kappanames("nu", REALSXP, "a", REALSXP, "c", REALSXP);
  range_type *range = ra->ranges;
//      printf("ij=%d \n", cov->vdim);

  // *nu = cov->p[0], 
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-5;
  range->pmax[0] = 10.0;
  range->openmin[0] = true;
  range->openmax[0] = true;
      
  //   *a = cov->p[1], 
   range->min[1] = 0.0;
  range->max[1] = RF_INF;
  range->pmin[1] = 1e-5;
  range->pmax[1] = 10.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
  
  //    *c = cov->p[2]; 
  range->min[2] = 0;
  range->max[2] = RF_INF;
  range->pmin[2] = 0;
  range->pmax[2] = 1000;
  range->openmin[2] = false;
  range->openmax[2] = true;  
  
 //    *rho = cov->p[3]; 
  range->min[3] = -1;
  range->max[3] = 1;
  range->pmin[3] = -1;
  range->pmax[3] = 1;
  range->openmin[3] = false;
  range->openmax[3] = false;    


}




void Nonestat(double *x, cov_model *cov, double *v){
    v[0] = RF_INF;
}
void Nonenonstat(double *x, double *y, cov_model *cov, double *v){
    v[0] = RF_INF;
}
void rangeNone(cov_model *cov, range_arraytype* ra){ }
