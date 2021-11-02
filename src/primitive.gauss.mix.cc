/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions that are gaussian scale mixtures

 Copyright (C) 2017 -- 2018 Olga Moreva (bivariate models) & Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/



#include "def.h"
#include <Basic_utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <General_utils.h>

#include <zzz_RandomFieldsUtils.h>
#include "extern.h"
#include "questions.h"
#include "primitive.h"
#include "operator.h"
#include "AutoRandomFields.h"
#include "shape.h"
#include "Processes.h"
//#include "xport_import.h"


//#define LOG2 M_LN2

#define i11 0
#define i21 1
#define i22 2

#define epsilon 1e-15

/* Cauchy models */
#define CAUCHY_GAMMA 0
void Cauchy(double *x, model *cov, double *v){
  double gamma = P0(CAUCHY_GAMMA);
  *v = POW(1.0 + *x * *x, -gamma);
}
void logCauchy(double *x, model *cov, double *v, double *Sign){
  double gamma = P0(CAUCHY_GAMMA);
  *v = LOG(1.0 + *x * *x) * -gamma;
  *Sign = 1.0;
}
void TBM2Cauchy(double *x, model *cov, double *v){
  double gamma = P0(CAUCHY_GAMMA), y2, lpy2;
  y2 = *x * *x; 
  lpy2 = 1.0 + y2;
  switch ((int) (gamma * 2.0 + 0.001)) {// ueber check sichergestellt
  case 1 : *v = 1.0 / lpy2; break;
  case 3 : *v = (1.0 - y2)/ (lpy2 * lpy2); break;
  case 5 : *v = (1.0-y2*(2.0+0.333333333333333333333*y2))/(lpy2*lpy2*lpy2);
    break;
  case 7 : lpy2 *= lpy2; 
    *v = (1.0- y2*(3.0+y2*(1.0+0.2*y2)))/(lpy2 * lpy2);
    break;
  default :
    ERR("TBM2 for cauchy only possible for alpha=0.5 + k; k=0, 1, 2, 3 ");
  }
}
void DCauchy(double *x, model *cov, double *v){
  double y=*x, gamma = P0(CAUCHY_GAMMA);
  *v = (-2.0 * gamma * y) * POW(1.0 + y * y, -gamma - 1.0);
}
void DDCauchy(double *x, model *cov, double *v){
  double ha = *x * *x, gamma = P0(CAUCHY_GAMMA);
  *v = 2.0 * gamma * ((2.0 * gamma + 1.0) * ha - 1.0) * 
    POW(1.0 + ha, -gamma - 2.0);
}
void InverseCauchy(double*x, model *cov, double *v){
  double
    gamma = P0(CAUCHY_GAMMA);
  if (*x == 0.0) *v = RF_INF;
  else *v = SQRT(POW(*x, -1.0 / gamma) - 1.0);
}
int checkCauchy(model VARIABLE_IS_NOT_USED  *cov){
  RETURN_NOERROR;
}
void rangeCauchy(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[CAUCHY_GAMMA] = 0.0;
  range->max[CAUCHY_GAMMA] = RF_INF;
  range->pmin[CAUCHY_GAMMA] = 0.09;
  range->pmax[CAUCHY_GAMMA] = 10.0;
  range->openmin[CAUCHY_GAMMA] = true;
  range->openmax[CAUCHY_GAMMA] = true;
}
void coinitCauchy(model VARIABLE_IS_NOT_USED *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0; //  q[CUTOFF_A] 
  li->msg[0] = MSGLOCAL_JUSTTRY;
}
void DrawMixCauchy(model VARIABLE_IS_NOT_USED *cov, double *random) { //better GR 3.381.4 ?? !!!!
  // split into parts <1 and >1 ?>
  *random = -LOG(1.0 -UNIFORM_RANDOM);
}
double LogMixDensCauchy(double VARIABLE_IS_NOT_USED *x, double logV, model *cov) {
  double gamma = P0(CAUCHY_GAMMA);
  // da dort 1/w vorgezogen 
  return 0.5 * (gamma - 1.0) * logV - 0.5 * lgammafn(gamma);
}


/** another Cauchy model */
#define CTBM_ALPHA 0
#define CTBM_BETA 1
#define CTBM_GAMMA 2
void Cauchytbm(double *x, model *cov, double *v){
  double ha, 
    alpha = P0(CTBM_ALPHA), 
    beta = P0(CTBM_BETA), 
    gamma = P0(CTBM_GAMMA),
    y=*x;
  if (y==0) {
    *v = 1.0;
  } else {
    ha = POW(y, alpha);
    *v = (1.0 + (1.0 - beta / gamma) * ha) * POW(1.0 + ha, -beta / alpha - 1.0);
  }
}

void DCauchytbm(double *x, model *cov, double *v){
  double y= *x, ha, 
    alpha = P0(CTBM_ALPHA), 
    beta = P0(CTBM_BETA),
    gamma = P0(CTBM_GAMMA);
  if (y == 0.0) *v = 0.0; // WRONG VALUE, but multiplied 
  else {                                  // by zero anyway
    ha = POW(y, alpha - 1.0);
    *v = beta *  ha * (-1.0 - alpha / gamma  + ha * y * (beta / gamma - 1.0)) *
      POW(1.0 + ha * y, -beta /alpha - 2.0);
  }
}


void rangeCauchytbm(model *cov, range_type *range){
  range->min[CTBM_ALPHA] = 0.0;
  range->max[CTBM_ALPHA] = 2.0;
  range->pmin[CTBM_ALPHA] = 0.05;
  range->pmax[CTBM_ALPHA] = 2.0;
  range->openmin[CTBM_ALPHA] = true;
  range->openmax[CTBM_ALPHA] = false;

  range->min[CTBM_BETA] = 0.0;
  range->max[CTBM_BETA] = RF_INF;
  range->pmin[CTBM_BETA] = 0.05;
  range->pmax[CTBM_BETA] = 10.0;
  range->openmin[CTBM_BETA] = true;
  range->openmax[CTBM_BETA] = true;
 
  range->min[CTBM_GAMMA] = (double) OWNLOGDIM(0);
  range->max[CTBM_GAMMA] = RF_INF;
  range->pmin[CTBM_GAMMA] = range->min[CTBM_GAMMA];
  range->pmax[CTBM_GAMMA] = range->pmin[CTBM_GAMMA] + 10.0;
  range->openmin[CTBM_GAMMA] = false;
  range->openmax[CTBM_GAMMA] = true;

}



// spatially constant (matrix); 
#define CONSTANT_M 0

void kappaconstant(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc) {
  *nr = *nc = i ==  CONSTANT_M ? 0 : -1;
}

void constant(double VARIABLE_IS_NOT_USED *x, model *cov, double *v){
  int 
    vdimSq = VDIM0 * VDIM1;
  MEMCOPY(v, P(CONSTANT_M), vdimSq * sizeof(double));
}

void nonstatconstant(double VARIABLE_IS_NOT_USED *x, double VARIABLE_IS_NOT_USED *y, model *cov, double *v){
  int 
    vdimSq = VDIM0 * VDIM1;
  MEMCOPY(v, P(CONSTANT_M), vdimSq * sizeof(double));
}

int checkconstant(model *cov) {
  int err; // anzahl listen elemente
 
  if (GLOBAL.internal.warn_constant) {
    GLOBAL.internal.warn_constant = false;
    warning("NOTE that the definition of 'RMconstant' has changed in version 3.0.61. Maybe  'RMfixcov' is the model your are looking for");
  }
 
  VDIM0 = VDIM1 = cov->nrow[CONSTANT_M];
  if (VDIM0 != VDIM1)  RETURN_ERR(ERROR_MATRIX_SQUARE);

  if (equalsVariogram(OWNTYPE(0))) { // ACHTUNG wirklich genau VariogramType
    SERR("strange call");
  }

  if (cov->q != NULL) { 
    return (int) cov->q[0];
  } else QALLOC(1);
  cov->q[0] = NOERROR;
  
  // check whether positive eigenvalue  
  if (!Ext_is_positive_definite(P(CONSTANT_M), VDIM0)) {
    cov->monotone = MONOTONE;
    cov->ptwise_definite = pt_indef;
    if (isnowPosDef(cov)) return cov->q[0] = ERROR_MATRIX_POSDEF;
  } else {
    cov->monotone = COMPLETELY_MON;
    cov->ptwise_definite = pt_posdef;
  }


  cov->matrix_indep_of_x = true;
  int
    vdim = VDIM0,
    vdimP1 = vdim + 1;
  double *p = P(CONSTANT_M);
  for (int i=0; i<vdim; i++) cov->mpp.maxheights[i] = p[i * vdimP1];
  err = checkkappas(cov);

 
  RETURN_ERR(err);
}


void rangeconstant(model VARIABLE_IS_NOT_USED *cov, range_type *range){
 // auch in rangec verwendet!
  int dim = VDIM0;
  if (isnowPosDef(cov)) {
    if (isnowTcf(cov)) {
      range->min[CONSTANT_M] = (int) (dim == 1);
      range->max[CONSTANT_M] = 1;
      range->pmin[CONSTANT_M] = range->min[CONSTANT_M];
      range->pmax[CONSTANT_M] = 1;
      range->openmin[CONSTANT_M] = false;
      range->openmax[CONSTANT_M] = false;    
    } else {
      range->min[CONSTANT_M] = dim == 1 ? 0 : RF_NEGINF;
      range->max[CONSTANT_M] = RF_INF;
      range->pmin[CONSTANT_M] = dim == 1 ? 0 : -1e10;
      range->pmax[CONSTANT_M] = 1e10;
      range->openmin[CONSTANT_M] = dim != 1;
      range->openmax[CONSTANT_M] = true;    
    }
  } else {
    range->min[CONSTANT_M] = RF_NEGINF;
    range->max[CONSTANT_M] = RF_INF;
    range->pmin[CONSTANT_M] = - 1e10;
    range->pmax[CONSTANT_M] = 1e10;
    range->openmin[CONSTANT_M] = true;
    range->openmax[CONSTANT_M] = true;
  }
}




/* exponential model */
void exponential(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
   *v = EXP(- *x);
   //   printf(": %10g %10g\n", *x, *v);
}
void logexponential(double *x, model VARIABLE_IS_NOT_USED *cov, double *v, double *Sign){
  *v = - *x;
  *Sign = 1.0;
 }
void TBM2exponential(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) 
{
  double y = *x;
  *v = 1.0;
  if (y!=0.0) *v = 1.0 - PIHALF * y * Ext_I0mL0(y);
}
void Dexponential(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  *v = - EXP(- *x);
}
void DDexponential(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
 *v = EXP(-*x);
}
void Inverseexponential(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  *v = (*x == 0.0) ? RF_INF : -LOG(*x);
}

void nonstatLogInvExp(double *x, model *cov,
		      double *left, double *right){
  double 
    z = *x <= 0.0 ? - *x : 0.0;  
  int d,
    dim = PREVLOGDIM(0);
  
  for (d=0; d<dim; d++) {
    left[d] = -z;
    right[d] = z;
  }
}


double densityexponential(double *x, model *cov) {
  // to do: ersetzen durch die Familien 

  // spectral density
  int d,
    dim=PREVLOGDIM(0);
  double x2 = 0.0,
    dim12 = 0.5 * ((double) (dim + 1));
  for (d=0; d<dim; d++) x2 += x[d] * x[d];

  //printf("de = %f %d\n",
  
  return gammafn(dim12) * POW(M_PI * (1.0 + x2), -dim12);
}


int initexponential(model *cov, gen_storage *s) {
 int dim = PREVLOGDIM(0);
  double D = (double) dim;

  //  TREE0(cov);
  //  PMI(cov);
  

  if (HAS_SPECTRAL_FRAME(cov)) {
    spec_properties *cs = &(s->spec);
   if (PREVLOGDIM(0) <= 2) RETURN_NOERROR;
    cs->density = densityexponential;
    return search_metropolis(cov, s);
  }

  else if (hasSmithFrame(cov)) {
    //Inverseexponential(&GLOBAL.mpp. about_ zero, NULL, &(cov->mpp.refradius));
    //  R *= GLOBAL.mpp.radius_natscale_factor;
    
    /*
    if (cov->mpp.moments >= 1) {
      int xi, d;
      double i[3], dimfak, Rfactor, sum,
	dimHalf = 0.5 * D;
      dimfak = gammafn(D);
      for (xi=1; xi<=2; xi++) {
	R = xi * cov->mpp.refradius;  //
	// http://de.wikipedia.org/wiki/Kugel
	i[xi] = dimfak * 2 * POW(M_PI / (double) (xi*xi), dimHalf) / 
	  gammafn(dimHalf);
	
	if (R < RF_INF) {
	  // Gradstein 3.351 fuer n=d-1
	  
	  //printf("here\n");
	  for (sum = 1.0, factor=1.0, d=1; d<dim; d++) {
	    factor *= R / D;
	    sum += factor;
	    //	printf("%d %10g %10g %10g\n", d, sum, factor, R);
	}
	  sum *= dimfak * EXP(-R);
	  // printf("%d %10g %10g %10g\n", d, sum, factor, R);
	  i[xi] -= sum;
	}
      }
      cov->mpp.mM[1] = cov->mpp.mMplus[1] = i[1];
      if (cov->mpp.moments >= 2) {
	cov->mpp.mM[2] = cov->mpp.mMplus[2] = i[2];
      }
    }
    */
  
    assert(cov->mpp.maxheights[0] == 1.0);
    if (cov->mpp.moments >= 1) {
       cov->mpp.mM[1]= cov->mpp.mMplus[1] = 
	SurfaceSphere(dim - 1, 1.0) * gammafn(D);
    }
    
  }

  else if (hasRandomFrame(cov)) { RETURN_NOERROR; }
  
  else ILLEGAL_FRAME;

  RETURN_NOERROR;
}
void do_exp(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 
}

void spectralexponential(model *cov, gen_storage *S, double *e) {
  spectral_storage *s = &(S->Sspectral);
  int dim = PREVLOGDIM(0);
  if (dim <= 2) {
    double A = 1.0 - UNIFORM_RANDOM;
    E12(s, dim, SQRT(1.0 / (A * A) - 1.0), e);
    //    printf("dim = %d %f %f %f\n", dim, A, e[0]);
  } else {
    metropolis(cov, S, e);
  }
}

int checkexponential(model *cov) {
  int dim = OWNLOGDIM(0);
  if (dim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (dim != 2) cov->pref[Hyperplane] = PREF_NONE;
  RETURN_NOERROR;
}

int hyperexponential(double radius, double *center, double *rx,
		     model VARIABLE_IS_NOT_USED *cov, bool simulate, 
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
    //
    q=0;
    //   q = RF_NA; assert(false);
  int k,
    Err;
  
 assert(OWNLOGDIM(0)==2);

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
  if ((hx=*Hx=(double *) MALLOC(sizeof(double) * (p + 8 * sizeof(int))))==NULL||
      (hy=*Hy=(double *) MALLOC(sizeof(double) * (p + 8 *sizeof(int))))==NULL||
      (hr=*Hr=(double *) MALLOC(sizeof(double) * (p + 8 * sizeof(int))))==NULL){
    Err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }  
  
  /* creating the lines; some of the lines are not relevant as they
     do not intersect the rectangle investigated --> k!=4
     (it is checked if all the corners of the rectangle are on one 
     side (?) )
  */
   for(i=0; i<p; i++) {
    phi = UNIFORM_RANDOM * TWOPI;
    hx[q] = COS(phi);     hy[q] = SIN(phi);    
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

  return q;

 ErrorHandling:
  return -Err;
}

void coinitExp(model VARIABLE_IS_NOT_USED *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0; //  q[CUTOFF_A] 
  li->msg[0] = MSGLOCAL_OK;
}
void ieinitExp(model VARIABLE_IS_NOT_USED *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0;
  li->msg[0] = MSGLOCAL_OK;
}
void DrawMixExp(model VARIABLE_IS_NOT_USED *cov, double *random) {
  // GR 3.325: int_-infty^infty EXP(x^2/2 - b/x^2) = EXP(-\sqrt b)
  double x = GAUSS_RANDOM(1.0);
  *random = 1.0 / (x * x);
}
double LogMixDensExp(double VARIABLE_IS_NOT_USED *x, double VARIABLE_IS_NOT_USED logV, model VARIABLE_IS_NOT_USED *cov) {
  // todo: check use of mixdens --- likely to programme it now completely differently 
  return 0.0;
}



/* Gausian model */
void Gauss(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = EXP(- *x * *x);
  //printf("x =%10g v=%10g\n", x[0], *v);
}
void logGauss(double *x, model VARIABLE_IS_NOT_USED  *cov, double *v, double *Sign) {
  *v = - *x * *x;
  *Sign = 1.0;
}
void DGauss(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x; 
  *v = -2.0 * y * EXP(- y * y);
}
void DDGauss(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x * *x; 
  *v = (4.0 * y - 2.0)* EXP(- y);
}
void D3Gauss(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x * *x; 
  *v = *x * (12 - 8 * y) * EXP(- y);
}
void D4Gauss(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x * *x; 
  *v = ((16 * y - 48) * y + 12) * EXP(- y);
}
void InverseGauss(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = SQRT(-LOG(*x));
}
void nonstatLogInvGauss(double *x, model VARIABLE_IS_NOT_USED *cov, 
		     double *left, double *right) {
  double 
    z = 0.0;
  if (*x < 0.0) z =SQRT(- *x);  
  int d,
    dim = PREVLOGDIM(0);
  for (d=0; d<dim; d++) {
    left[d] = -z;
    right[d] = z;
  }
}
double densityGauss(double *x, model *cov) {  
    int d, dim=PREVLOGDIM(0);
    double x2=0.0;
    for (d=0; d<dim; d++) x2 += x[d] * x[d];
    return EXP(- 0.25 * x2 - (double) dim * (M_LN2 + M_LN_SQRT_PI));
}
int struct_Gauss(model *cov, model **newmodel) {  
  ASSERT_NEWMODEL_NOT_NULL;

  //printf("dleete next lines\n");  ILLEGAL_FRAME_STRUCT;
  // masse von Rect ist ok (unnormedmass)
  // supportmin / max macht ergebnis auch nicht besser

  switch (cov->frame) {
  case PoissonGaussType :// optimierte density fuer den Gauss-Fall
    double invscale;
    addModel(newmodel, GAUSS, cov);       
    addModel(newmodel, DOLLAR);
    kdefault(*newmodel, DSCALE, INVSQRTTWO);
    addModel(newmodel, TRUNCSUPPORT);
    InverseGauss(&GLOBAL.mpp.about_zero, cov, &invscale);
    kdefault(*newmodel, TRUNC_RADIUS, invscale);
    break;  
  default :
    if (hasSmithFrame(cov)) {  // optimierte density fuer den Gauss-Fall
      // crash();
      addModel(newmodel, GAUSS_DISTR, cov); // to
      kdefault(*newmodel, GAUSS_DISTR_MEAN, 0.0);
      kdefault(*newmodel, GAUSS_DISTR_SD, INVSQRTTWO);
      RETURN_NOERROR;
    } 
    ILLEGAL_FRAME_STRUCT;
  }

  RETURN_NOERROR;
}

double IntUdeU2_intern(int d, double R, double expMR2) {
  // int_0^R u^{d-1} EXP(-u^2) \D u
  if (d==0) return (pnorm(R, 0.0, INVSQRTTWO, true, false) - 0.5)  * SQRTPI;
  else if (d == 1) return 0.5  * (1.0 - expMR2);
  return 0.5 * (expMR2 + (d - 1.0) * IntUdeU2_intern(d - 2, R, expMR2));
}

double IntUdeU2(int d, double R) {
  // int_0^R u^{d-1} EXP(-u^2) \D u
  return IntUdeU2_intern(d, R, EXP(-R*R));
}

int initGauss(model *cov, gen_storage *s) {

  //  if (hasAnyFrame(cov)) RETURN_NOERROR;

  if (HAS_SPECTRAL_FRAME(cov)) {

   spec_properties *cs = &(s->spec);
  
   if (OWNLOGDIM(0) <= 2) RETURN_NOERROR;
    cs->density = densityGauss;
    return search_metropolis(cov, s);
  }

  else if (hasSmithFrame(cov)) {
    // int_{b(0,R) e^{-a r^2} dr = d b_d int_0^R r^{d-1} e^{-a r^2} dr
    // where a = 2.0 * xi / sigma^2
    // here : 2.0 is a factor so that the mpp function leads to the
    //            gaussian covariance model EXP(-x^2)
    //        xi : 1 : integral ()
    //             2 : integral ()^2

   
   double R = RF_INF;
    int 
      dim = OWNLOGDIM(0);
    
    assert(COVNR == GAUSS);
         
  
    assert(cov->mpp.maxheights[0] = 1.0);
    if (cov->mpp.moments >= 1) {
      cov->mpp.mM[1] = cov->mpp.mMplus[1] = 
	SurfaceSphere(dim - 1, 1.0) * IntUdeU2(dim - 1, R);
      int i;
      for (i=2; i<=cov->mpp.moments; i++) {
	cov->mpp.mM[i] = cov->mpp.mM[1] * POW((double) i, -0.5 * dim);
      }
    }    
    //    cov->mpp.maxheights[0] = 1.0;
  }

  else if (hasRandomFrame(cov) || hasAnyPoissonFrame(cov)) { RETURN_NOERROR; }

  else ILLEGAL_FRAME;
 
 
  RETURN_NOERROR;

}

void do_Gauss(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 

}

void spectralGauss(model *cov, gen_storage *S, double *e) {   
  spectral_storage *s = &(S->Sspectral);
  int dim = PREVLOGDIM(0);
  if (dim <= 2) {
    E12(s, dim, 2.0 * SQRT(-LOG(1.0 - UNIFORM_RANDOM)), e);

  } else {
    metropolis(cov, S, e);
  }
}
void DrawMixGauss(model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *random) {
  *random = 1.0;
}
double LogMixDensGauss(double VARIABLE_IS_NOT_USED *x, double VARIABLE_IS_NOT_USED logV, model VARIABLE_IS_NOT_USED *cov) {
  return 0.0;
}

/*
void densGauss(double *x, model *cov, double *v) {
  int 
    factor[MAXMPPDIM+1] = {0, 1 / M_SQRT_PI, INVPI, INVPI / M_SQRT_PI, 
			   INVPI * INVPI},
    dim = cov->tsdim;
    *v = factor[dim] * EXP(- *x * *x);
}
*/

/*
void getMassGauss(double *a, model *cov, double *kappas, double *m) {
  int i, j, k, kold,
    dim = cov->tsdim;
  double val[MAXMPPDIM + 1],
    sqrt2pi = SQRT2 * SQRTPI,
    factor[MAXMPPDIM+1] = {1, 
			   SQRT(2) / M_SQRT_PI, 
			   2 * INVPI,
			   2 * SQRT(2) * INVPI / M_SQRT_PI, 
			   4 * INVPI * INVPI};
  
  val[0] = 1.0;
  for (i=1; i<=dim; i++) 
    val[i] = (2.0 * pnorm(SQRT2 * a[i], 0.0, 1.0, false, false) - 1.0) * M_SQRT_PI;
  for (k=kold=i=0; i<dim; i++) {
    m[k++] = val[i];
    for (j=1; j< kold; j++) m[k++] = val[i] * m[j];
    kold = k;
    pr intf("kold = %d k=%d\n", kold, k);
  }
}
*/

/*
void simuGauss(model *cov, int dim, double *v) {
  int i;
  double 
    dummy;
  *v = 0.0;
  if (dim <= 2) {
    *v = dim == 1 ? FABS(GAUSS_RANDOM(1.0)) : rexp(1.0); 
  } else {
    for (i=0; i<dim; i++) {
      dummy = GAUSS_RANDOM(1.0);
      *v += dummy * dummy;
    }
    *v = SQRT(*v);
  }
}
*/


/* gencauchy */
#define GENC_ALPHA 0
#define GENC_BETA 1
void generalisedCauchy(double *x, model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA);
  *v = POW(1.0 + POW(*x, alpha), -beta/alpha);
}
void DgeneralisedCauchy(double *x, model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA), ha, y=*x;
  if (y ==0.0) {
    *v = ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? RF_NEGINF : -beta); 
  } else {
    ha = POW(y, alpha - 1.0);
    *v = -beta * ha * POW(1.0 + ha * y, -beta / alpha - 1.0);
  }
}
void DDgeneralisedCauchy(double *x, model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA), ha, y=*x;
  if (y ==0.0) {
    *v = ((alpha==1.0) ? beta * (beta + 1.0) 
	  : (alpha==2.0) ?  -beta 
	  : (alpha < 1) ? RF_INF : RF_NEGINF); 
  } else {
    ha = POW(y, alpha);
    *v = beta * ha / (y * y) * (1.0 - alpha + (1.0 + beta) * ha)
      * POW(1.0 + ha, -beta / alpha - 2.0);
  }
}

void D3generalisedCauchy(double *x, model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA), ha, y=*x;
  if (y ==0.0) {
    *v = ((alpha==2.0) ? 0  :  (alpha == 1)? -beta*(beta+1)*(beta+2) :
                                             (alpha < 1)? RF_NEGINF : RF_INF);
  } else {
    ha=POW(y, alpha);
    *v = -beta * ha / (y * y* y) * ( (beta + 1)*(beta + 2)*ha*ha  -
                    (3*beta + alpha + 4)*(alpha - 1)*ha + (alpha-1)*(alpha -2) )
      * POW(1.0 + ha, -beta / alpha - 3.0);
  }
}

void D4generalisedCauchy(double *x, model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA), y=*x;
  if (y ==0.0) {
    *v = ((alpha==2.0) ? 3*beta*(beta + 2)  :  (alpha == 1)?
	  beta*(beta+1)*(beta+2)*(beta + 3) :
	  (alpha < 1)? RF_INF : RF_NEGINF);
  } else {
    int ha=POW(y, alpha),
      haSq = ha * ha;
    *v = beta * ha / (y * y* y *y) *
      ( -(alpha-1) * (alpha-2)*(alpha-3) -
    (alpha - 1) * ( 4*alpha*beta + alpha * (alpha + 7) +
		     6*beta*beta + 22*beta + 18 ) * haSq +
	(alpha - 1) * (alpha * (4 * alpha + 7 * beta + 4) - 11*beta - 18) * ha +
	(beta + 1) * (beta +2) * (beta + 3) * ha * haSq) *
      POW(1.0 + ha, -beta / alpha - 4.0);
  }
}


void loggeneralisedCauchy(double *x, model *cov, double *v, double *Sign){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA);
  *v = LOG(1.0 + POW(*x, alpha)) *  -beta/alpha;
  *Sign = 1.0;
}
void InversegeneralisedCauchy(double *x, model *cov, double *v) {
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA);
  *v = RF_INF;
  if (*x != 0.0) *v = POW(POW(*x, -alpha / beta) - 1.0, 1.0 / alpha);
  // MLE works much better with 0.01 then with 0.05
}
int checkgeneralisedCauchy(model *cov){
  if (OWNLOGDIM(0) > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  cov->monotone = P0(GENC_ALPHA) <= 1.0 ? COMPLETELY_MON : NORMAL_MIXTURE;
  RETURN_NOERROR;
}

void rangegeneralisedCauchy(model *cov, range_type *range) {
  bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));
  
  range->min[GENC_ALPHA] = 0.0;
  range->max[GENC_ALPHA] = tcf ? 1.0 : 2.0;
  range->pmin[GENC_ALPHA] = 0.05;
  range->pmax[GENC_ALPHA] = range->max[GENC_ALPHA];
  range->openmin[GENC_ALPHA] = true;
  range->openmax[GENC_ALPHA] = false;

  range->min[GENC_BETA] = 0.0;
  range->max[GENC_BETA] = RF_INF;
  range->pmin[GENC_BETA] = 0.05;
  range->pmax[GENC_BETA] = 10.0;
  range->openmin[GENC_BETA] = true;
  range->openmax[GENC_BETA] = true;
 
}

void coinitgenCauchy(model *cov, localinfotype *li) {
  double 
    thres[2] = {0.5, 1.0}, 
    alpha=P0(GENC_ALPHA); 
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
void ieinitgenCauchy(model *cov, localinfotype *li) {
  li->instances = 1;
  if (P0(GENC_ALPHA) <= 1.0) {
    li->value[0] = 1.0;
    li->msg[0] = MSGLOCAL_OK;
  } else {
    li->value[0] = 1.5;
    li->msg[0] = MSGLOCAL_JUSTTRY;
  }
}


/*---------------------------------*/


// ************************* bivariate Cauchy


void kappa_biCauchy(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
    *nc = 1;
    *nr = i == BICauchyalpha || i ==  BICauchybeta|| i == BICauchyscale ? 3 :  1;
}

int checkbiCauchy(model *cov) {
  if (OWNLOGDIM(0) > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  RETURN_NOERROR;
}


void rangebiCauchy(model VARIABLE_IS_NOT_USED *cov, range_type *range){

  range->min[BICauchyalpha] = 0.0;
  range->max[BICauchyalpha] = 1;

  range->pmin[BICauchyalpha] = 0.05;
  range->pmax[BICauchyalpha] = 1;

  range->openmin[BICauchyalpha] = true;
  range->openmax[BICauchyalpha] = true;


  range->min[BICauchybeta] = 0.0;
  range->max[BICauchybeta] = RF_INF;
  range->pmin[BICauchybeta] = 0.05;
  range->pmax[BICauchybeta] = 10.0;
  range->openmin[BICauchybeta] = true;
  range->openmax[BICauchybeta] = true;

  range->min[BICauchyscale] = 0.0;
  range->max[BICauchyscale] = RF_INF;
  range->pmin[BICauchyscale] = 1e-2;
  range->pmax[BICauchyscale] = 100.0;
  range->openmin[BICauchyscale] = true;
  range->openmax[BICauchyscale] = true;

  // to do: check rhos

  range->min[BICauchyrho] = - 1;
  range->max[BICauchyrho] = 1;
  range->pmin[BICauchyrho] = - 1;
  range->pmax[BICauchyrho] = 1;
  range->openmin[BICauchyrho] = false;
  range->openmax[BICauchyrho] = false;

}



void coinitbiCauchy(model *cov, localinfotype *li) {
  double
      thres = 1,
      *alpha = P(BICauchyalpha);

  if ( (alpha[0] <= thres) &&  (alpha[1] <= thres) && (alpha[2] <= thres) ) {
      li->instances = 1;
      li->value[0] = CUTOFF_THIRD_CONDITION; //  q[CUTOFF_A]
      li->msg[0] =MSGLOCAL_OK;
    }
  else {
      li->instances = 1;
      li->value[0] = CUTOFF_THIRD_CONDITION ; //  q[CUTOFF_A]
      li->msg[0] = MSGLOCAL_JUSTTRY;
    }
}

void biCauchy (double *x, model *cov, double *v) {
    int i;
    double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA),
            y = *x, z;

    assert(BICauchyalpha == GENC_ALPHA);
    assert(BICauchybeta == GENC_BETA);

     for (i=0; i<3; i++) {
        z = y/P(BICauchyscale)[i];
        P(GENC_ALPHA)[0] = P(BICauchyalpha)[i];
        P(GENC_BETA)[0] = P(BICauchybeta)[i];
        generalisedCauchy(&z, cov, v + i);
    }
     P(BICauchyalpha)[0] = alpha;
     P(BICauchybeta)[0] = beta;
     v[3] = v[2];
     v[1] *= P0(BICauchyrho);
     v[2] = v[1];
}



void DbiCauchy(double *x, model *cov, double *v) {
    int i;
    double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA),
            y = *x, z;

    assert(BICauchyalpha == GENC_ALPHA);
    assert(BICauchybeta == GENC_BETA);

     for (i=0; i<3; i++) {
        z = y/P(BICauchyscale)[i];
        P(GENC_ALPHA)[0] = P(BICauchyalpha)[i];
        P(GENC_BETA)[0] = P(BICauchybeta)[i];
        DgeneralisedCauchy(&z, cov, v + i);
        v[i] /=P(BICauchyscale)[i];
    }
    P(BICauchyalpha)[0] = alpha;
    P(BICauchybeta)[0] = beta;
    v[3] = v[2];
    v[1] *= P0(BICauchyrho);
    v[2] = v[1];
}

void DDbiCauchy(double *x, model *cov, double *v) {
    int i;
    double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA),
            y = *x, z;

    assert(BICauchyalpha == GENC_ALPHA);
    assert(BICauchybeta == GENC_BETA);

     for (i=0; i<3; i++) {
        z = y/P(BICauchyscale)[i];
        P(GENC_ALPHA)[0] = P(BICauchyalpha)[i];
        P(GENC_BETA)[0] = P(BICauchybeta)[i];
        DDgeneralisedCauchy(&z, cov, v + i);
        v[i] /=(P(BICauchyscale)[i]*P(BICauchyscale)[i]);
    }
     P(BICauchyalpha)[0] = alpha;
     P(BICauchybeta)[0] = beta;
     v[3] = v[2];
     v[1] *= P0(BICauchyrho);
     v[2] = v[1];
}

void D3biCauchy(double *x, model *cov, double *v) {
    int i;
    double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA),
            y = *x, z;

    assert(BICauchyalpha == GENC_ALPHA);
    assert(BICauchybeta == GENC_BETA);

     for (i=0; i<3; i++) {
        z = y/P(BICauchyscale)[i];
        P(GENC_ALPHA)[0] = P(BICauchyalpha)[i];
        P(GENC_BETA)[0] = P(BICauchybeta)[i];
        D3generalisedCauchy(&z, cov, v + i);
        v[i] /=(P(BICauchyscale)[i]*P(BICauchyscale)[i]*P(BICauchyscale)[i]);
    }
     P(BICauchyalpha)[0] = alpha;
     P(BICauchybeta)[0] = beta;
     v[3] = v[2];
     v[1] *= P0(BICauchyrho);
     v[2] = v[1];
}

void D4biCauchy(double *x, model *cov, double *v) {
    int i;
    double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA),
            y = *x, z;

    assert(BICauchyalpha == GENC_ALPHA);
    assert(BICauchybeta == GENC_BETA);

     for (i=0; i<3; i++) {
        z = y/P(BICauchyscale)[i];
        P(GENC_ALPHA)[0] = P(BICauchyalpha)[i];
        P(GENC_BETA)[0] = P(BICauchybeta)[i];
        D4generalisedCauchy(&z, cov, v + i);
        double dummy = P(BICauchyscale)[i]*P(BICauchyscale)[i];
        v[i] /=(dummy*dummy);

    }
     P(BICauchyalpha)[0] = alpha;
     P(BICauchybeta)[0] = beta;
     v[3] = v[2];
     v[1] *= P0(BICauchyrho);
     v[2] = v[1];
}


/* epsC -> generalised Cauchy : leading 1 is now an eps */
#define EPS_ALPHA 0
#define EPS_BETA 1
#define EPS_EPS 2
void epsC(double *x, model *cov, double *v){
  double 
    alpha = P0(EPS_ALPHA), 
    beta=P0(EPS_BETA),
    eps=P0(EPS_EPS);
  *v = POW(eps + POW(*x, alpha), -beta/alpha);
 }
void logepsC(double *x, model *cov, double *v, double *Sign){
  double 
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  *v = LOG(eps + POW(*x, alpha)) * -beta/alpha;
  *Sign = 1.0;
 }
void DepsC(double *x, model *cov, double *v){
  double ha, 
    y=*x,
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  if (y ==0.0) {
    *v = (eps == 0.0) ? RF_NEGINF : 
      ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? RF_NEGINF : -beta); 
  } else {
    ha = POW(y, alpha - 1.0);
    *v = -beta * ha * POW(eps + ha * y, -beta / alpha - 1.0);
  }
}
void DDepsC(double *x, model *cov, double *v){
  double ha, 
    y=*x,
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  if (y ==0.0) {
    *v = (eps == 0.0) ? RF_INF : ((alpha==2.0) ? beta * (beta + 1.0) : RF_INF); 
  } else {
    ha = POW(y, alpha);
    *v = beta * ha / (y * y) * ( (1.0 - alpha) * eps + (1.0 + beta) * ha)
      * POW(eps + ha, -beta / alpha - 2.0);
  }
}
void InverseepsC(double *x, model *cov, double *v){
  double 
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  *v = RF_INF;
  if (*x != 0.0) *v = POW(POW(*x, -alpha / beta) - eps, 1.0 / alpha);
}
int checkepsC(model *cov){
  double eps=P0(EPS_ALPHA);
  int i, err;
  if (OWNLOGDIM(0) > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  kdefault(cov, EPS_ALPHA, 1.0); 
  kdefault(cov, EPS_BETA, 1.0); 
  kdefault(cov, EPS_EPS, 0.0); 
  if (ISNAN(eps) || eps == 0.0) {
    //  OWNDOM(0)=GENERALISEDCOVARIANCE; // later
    for (i=CircEmbed; i<Nothing; i++) cov->pref[i] = PREF_NONE;
  }
  
  RETURN_NOERROR;
}

void rangeepsC(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[EPS_ALPHA] = 0.0;
  range->max[EPS_ALPHA] = 2.0;
  range->pmin[EPS_ALPHA] = 0.05;
  range->pmax[EPS_ALPHA] = 2.0;
  range->openmin[EPS_ALPHA] = true;
  range->openmax[EPS_ALPHA] = false;

  range->min[EPS_BETA] = 0.0;
  range->max[EPS_BETA] = RF_INF;
  range->pmin[EPS_BETA] = 0.05;
  range->pmax[EPS_BETA] = 10.0;
  range->openmin[EPS_BETA] = true;
  range->openmax[EPS_BETA] = true;

  range->min[EPS_EPS] = 0.0;
  range->max[EPS_EPS] = RF_INF;
  range->pmin[EPS_EPS] = 0.0;
  range->pmax[EPS_EPS] = 10.0;
  range->openmin[EPS_EPS] = true; // false for generlised covariance
  range->openmax[EPS_EPS] = true;
}



/* hyperbolic */
#define BOLIC_NU 0
#define BOLIC_XI 1
#define BOLIC_DELTA 2
void hyperbolic(double *x, model *cov, double *v){ 
  double Sign;
  loghyperbolic(x, cov, v, &Sign);
  *v = EXP(*v);
}
void loghyperbolic(double *x, model *cov, double *v, double *Sign){ 
  double 
    nu = P0(BOLIC_NU),
    xi=P0(BOLIC_XI), 
    delta=P0(BOLIC_DELTA),
    y=*x;
  *Sign = 1.0;
  if (y==0.0) { 
    *v = 0.0;
    return;
  } else if (y==RF_INF) {
    *v = RF_NEGINF;
    *Sign = 0.0;
    return;
  }
  if (delta==0) { // whittle matern
    if (nu > 80) warning("extremely imprecise results due to nu>80");
    *v = logWM(y * xi, nu, cov->q[WM_LOGGAMMA], 0.0);
  } else if (xi==0) { //cauchy   => NU2 < 0 !
    y /= delta;
    /* note change in Sign as NU2<0 */
    *v = LOG(1.0 + y * y) * 0.5 * nu; 
  } else {
    y=SQRT(delta * delta + y * y);  
    double  bk[MATERN_NU_THRES + 1L],
      xiy = xi * y,
      nuThres = nu <= MATERN_NU_THRES ? nu : MATERN_NU_THRES,
      factor = 1 / SQRT(nuThres);
    if  (xiy <= LOW_MATERN) {*v = 1.0; return;}
    if (xiy == RF_INF)  {*v = 0.0; return;} // also for cosine i.e. nu=-1/2

   *v = QVALUE3 + nu * LOG(y) + LOG(bessel_k_ex(xiy, nu, 2.0, bk)) - xiy;

   if (nu > MATERN_NU_THRES) { // factor!=0.0 && 
     //  printf("UU \n");
     double w,
       g = MATERN_NU_THRES / nu;
     y = 0.5 * xiy * factor;
     Gauss(&y, NULL, &w);
     *v = *v * g + (1.0 - g) * w;
   }
  }
}
void Dhyperbolic(double *x, model *cov, double *v)
{ 
  double 
    nu = P0(BOLIC_NU), 
    xi=P0(BOLIC_XI), 
    delta=P0(BOLIC_DELTA);
   double  
    y = *x;
  if (y==0.0) { 
    *v = 1.0;
    return;
  }
  if (delta==0) { // whittle matern
    *v = xi * Ext_DWM(y * xi, nu, cov->q[WM_LOGGAMMA], 0.0);
    *v *= xi;
  } else if (xi==0) { //cauchy
    double ha;
    y /= delta;
    ha = y * y;
    /* note change in Sign as NU2<0 */
    *v = nu * FABS(y) * POW(1.0 + ha, 0.5 * nu - 1.0) / delta;
  } else {
    double nuThres = nu <= MATERN_NU_THRES ? nu : MATERN_NU_THRES,
      s=SQRT(delta * delta + y * y),
      xi_s = xi * s,
      logs = LOG(s),
      bk[MATERN_NU_THRES + 1L];
    *v = - y * xi * EXP(QVALUE3 + (nuThres-1.0) * logs 
			+ LOG(bessel_k_ex(xi_s, nuThres-1.0, 2.0, bk)) - xi_s);
    if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu,
      factor = 1.0 /  SQRT(nuThres),
      scale = 0.5 * factor;
    y = xi_s * scale;
    DGauss(&y, NULL, &w);
    w *= scale;
    *v = *v * g + (1.0 - g) * w;
    }
  }
}

int inithyperbolic(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double nu = P0(BOLIC_NU),  
    nuThres = nu <= MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    xi=P0(BOLIC_XI), 
    delta=P0(BOLIC_DELTA),
    xd = xi * delta, // xidelta
    bk[MATERN_NU_THRES + 1L];
  
  QVALUE3 = xd - LOG(bessel_k_ex(xd, nuThres, 2.0, bk)) - nuThres * LOG(delta);
  if (nu > MATERN_NU_THRES) { // factor!=0.0 && 
     //  printf("UU \n");
     double w,
       factor = 1 / SQRT(nuThres),
       g = MATERN_NU_THRES / nu,
       y = 0.5 * xd * factor;
     Gauss(&y, NULL, &w);
     QVALUE3 =  QVALUE3 * g + (1.0 - g) * w;
   }

  
  if (!ISNA(delta) && delta == 0.0 && !ISNA(nu)) {
    assert(cov->q + WM_LOGGAMMA == &(QVALUE));
    assert(cov->q + WM_GAMMA == &(QVALUE2));
    cov->q[WM_LOGGAMMA] = lgammafn(nuThres); // QVALUE
    cov->q[WM_GAMMA] = gammafn(nuThres); // QVALUE2
  }
  RETURN_NOERROR;
}

int checkhyperbolic(model *cov){
  double 
    nu = P0(BOLIC_NU),
    xi=P0(BOLIC_XI),
    delta=P0(BOLIC_DELTA);
  int i;
  for (i=0; i<= Nothing; i++)
    cov->pref[i] *= ( ((bool) ISNAN(nu)) || nu < WhittleUpperNu[i]);
  if (nu>0) {
    if ((delta<0) || (xi<=0)) {
      SERR3("xi>0 and delta>=0 if nu>0. Got nu=%10g and delta=%10g for xi=%10g",
	    nu, delta, xi);
    }
  } else if (nu<0) {
    if ((delta<=0) || (xi<0)) {
       SERR3("xi>=0 and delta>0 if nu<0. Got nu=%10g and delta=%10g for xi=%10g",
	    nu, delta, xi);
    }
  } else { // nu==0.0
    if ((delta<=0) || (xi<=0)) {
      SERR3("xi>0 and delta>0 if nu=0. Got nu=%10g and delta=%10g for xi=%10g", 
	    nu, delta, xi);
    }
  }
  if (cov->q == NULL) {
    EXTRA_Q;
    inithyperbolic(cov, NULL);
  }
  RETURN_NOERROR;
}
void rangehyperbolic(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[BOLIC_NU] = RF_NEGINF;
  range->max[BOLIC_NU] = RF_INF;
  range->pmin[BOLIC_NU] = -20.0;
  range->pmax[BOLIC_NU] = 20.0;
  range->openmin[BOLIC_NU] = true;
  range->openmax[BOLIC_NU] = true;

  int i;
  for (i=1; i<=2; i++) { 
    range->min[i] = 0.0;
    range->max[i] = RF_INF;
    range->pmin[i] = 0.000001;
    range->pmax[i] = 10.0;
    range->openmin[i] = false;
    range->openmax[i] = true;
  }
}



/* stable model */
#define STABLE_ALPHA 0
void stable(double *x, model *cov, double *v){
  double y = *x, alpha = P0(STABLE_ALPHA);  
  *v = 1.0;
  if (y!=0.0) *v = EXP(-POW(y, alpha));
}
void logstable(double *x, model *cov, double *v, double *Sign){
  double y = *x, alpha = P0(STABLE_ALPHA);  
  *v = 0.0;
  if (y!=0.0) *v= -POW(y, alpha);
  *Sign = 1.0;
}
void Dstable(double *x, model *cov, double *v){
  double z, y = *x, alpha = P0(STABLE_ALPHA);
  if (y == 0.0) {
    *v = (alpha > 1.0) ? 0.0 : (alpha < 1.0) ? RF_NEGINF : -1.0;
  } else {
    z = POW(y, alpha - 1.0);
    *v = -alpha * z * EXP(-z * y);
  }
}
/* stable: second derivative at t=1 */
void DDstable(double *x, model *cov, double *v) 
{
  double z, y = *x, alpha = P0(STABLE_ALPHA), xalpha;
  if (y == 0.0) {
    // olga please check
    *v = (alpha == 1.0) ? 1.0 : (alpha == 2.0) ? -2.0 : //alpha * (1 - alpha)
	alpha < 1.0 ? RF_INF : RF_NEGINF;
  } else {
    z = POW(y, alpha - 2.0);
    xalpha = z * y * y;
    *v = alpha * (1.0 - alpha + alpha * xalpha) * z * EXP(-xalpha);
  }
}

void D3stable(double *x, model *cov, double *v)
{
  double z, y = *x, alpha = P0(STABLE_ALPHA), xalpha;
  if (y == 0.0) {
    // olga, please check; please note that RF_INF has been incorrect (reserved for integer variables only)
      *v = (alpha == 1.0) ? -1.0 : (alpha == 2.0) ? 0 :
	alpha < 1.0 ? RF_NEGINF : RF_INF;
  } else {
    z = POW(y, alpha - 3.0);
    xalpha = z * y * y * y;
    *v = -alpha * ( 3*alpha*(xalpha -1) + alpha*alpha*
                    (xalpha*xalpha - 3*xalpha + 1) + 2) * z * EXP(-xalpha);
  }
}

void D4stable(double *x, model *cov, double *v)
{
  double z, y = *x, alpha = P0(STABLE_ALPHA), xalpha;
  if (y == 0.0) {
      *v = (alpha == 1.0) ? 1.0 : (alpha == 2.0) ? 0 :
	alpha < 1.0 ? RF_INF : RF_NEGINF;
  } else {
    z = POW(y, alpha - 4.0);
    xalpha = z * y * y * y * y;
    *v = alpha*( alpha*alpha*alpha*(7*xalpha -6*xalpha*xalpha  +
                                    xalpha*xalpha*xalpha - 1 ) +
           +6*alpha*alpha*(-3*xalpha + xalpha*xalpha +1 )+
           11*alpha*(xalpha - 1 ) + 6 ) * z * EXP(-xalpha);

  }
}

void D5stable(double *x, model *cov, double *v)
{
  double z, y = *x, alpha = P0(STABLE_ALPHA), xalpha;
  if (y == 0.0) {
      *v = (alpha == 1.0) ? -1.0 : (alpha == 2.0) ? 0 :
	alpha < 1.0 ? RF_NEGINF : RF_INF;;
  } else {
    z = POW(y, alpha - 5.0);
    xalpha = z * y * y * y * y * y;
    *v = -alpha*( POW(alpha, 4)*(-15*xalpha +25*xalpha*xalpha -
                                 10*POW(xalpha, 3) + POW(xalpha, 4)+ 1 ) +
                  10*alpha*alpha*alpha*(7*xalpha  -6*xalpha*xalpha +
                                        POW(xalpha, 3) - 1 )+
                  35*alpha*alpha*(-3*xalpha + xalpha*xalpha +  1 )+
                  50*alpha*(xalpha -1 ) + 24 )* z * EXP(-xalpha);
  }
}


void Inversestable(double *x, model *cov, double *v){
  double y = *x, alpha = P0(STABLE_ALPHA);  
  *v = y>1.0 ? 0.0 : y == 0.0 ? RF_INF : POW( - LOG(y), 1.0 / alpha);
}
void nonstatLogInversestable(double *x, model *cov,
			     double *left, double *right){
  double 
    alpha = P0(STABLE_ALPHA),
    z = *x <= 0.0 ? POW( - *x, 1.0 / alpha) : 0.0;
  int d,
    dim = OWNLOGDIM(0);
  for (d=0; d<dim; d++) {
    left[d] = -z;
    right[d] = z;
  }
}


int checkstable(model *cov) {
  double alpha = P0(STABLE_ALPHA);
  if (OWNLOGDIM(0) > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (alpha == 2.0)
    cov->pref[CircEmbed] = 2;
  
  cov->monotone = alpha <= 1.0 ? COMPLETELY_MON : NORMAL_MIXTURE;

  RETURN_NOERROR;
}
void rangestable(model *cov, range_type *range){
  bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));
  range->min[STABLE_ALPHA] = 0.0;
  range->max[STABLE_ALPHA] = tcf ? 1.0 : 2.0;
  range->pmin[STABLE_ALPHA] = 0.06;
  range->pmax[STABLE_ALPHA] =  range->max[STABLE_ALPHA];
  range->openmin[STABLE_ALPHA] = true;
  range->openmax[STABLE_ALPHA] = false;
}
void coinitstable(model *cov, localinfotype *li) {
  coinitgenCauchy(cov, li);
}
void ieinitstable(model *cov, localinfotype *li) {  
  ieinitgenCauchy(cov, li);
}


/* DOUBLEISOTROPIC stable model for testing purposes only */
void stableX(double *x, model *cov, double *v){
  double y, alpha = P0(STABLE_ALPHA);
  y = x[0] * x[0] + x[1] * x[1];
  *v =  1.0;
  if (y!=0.0) *v = EXP(-POW(y, 0.5 * alpha));
}
void DstableX(double *x, model *cov, double *v){
  double z, y, alpha = P0(STABLE_ALPHA);
  y = x[0] * x[0] + x[1] * x[1];
  if (y == 0.0) {
    *v = ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? RF_INF : 1.0);
  } else {
    z = POW(y, 0.5 * alpha - 1.0);
    *v = -alpha * x[0] * z * EXP(- z * y);
  }
}
/* END DOUBLEISOTROPIC stable model for testing purposes only */



// ************************* bivariate powered exponential or bivariate stable
void kappa_biStable(int i, model VARIABLE_IS_NOT_USED *cov, int *nr,
                    int *nc){
  *nc = 1;
  *nr  = ( i == BIStablealpha || i == BIStablescale ) ? 3 :
           (i == BIStablecdiag || i == BIStablealphadiag ) ? 2 :
           (i == BIStablerho || i == BIStablebetared || i == BIStablerhored ) ?
                                                               1 : -1;
}


int checkbiStable(model *cov) {
  int err;
  gen_storage s;
  gen_NULL(&s);
  s.check = true;

  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  //   bistable_storage *S = cov->Sbistable;

  if (cov->Sbistable == NULL) {
    ONCE_NEW_STORAGE(bistable);
    bistable_storage *S = cov->Sbistable;
    S->alphadiag_given = !PisNULL(BIStablealphadiag);
    S->rhored_given = !PisNULL(BIStablerhored);
    // printf("check: alphadiag_given = %d\n", S->alphadiag_given);
  }

  if ((err=initbiStable(cov, &s)) != NOERROR) RETURN_ERR(err);
  VDIM0 = VDIM1 = 2;

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  //PMI(cov);

  RETURN_NOERROR;
}


// alphadiag always out
sortsofparam sortof_bistable(model *cov, int k,
			     int  VARIABLE_IS_NOT_USED *row,
			     int  VARIABLE_IS_NOT_USED *col,
			     sort_origin origin) {
  bistable_storage *S = cov->Sbistable;
  if (S == NULL) return UNKNOWNPARAM;
  switch(k) {
  case BIStablealphadiag : case BIStablebetared :
    return S->alphadiag_given || origin != original ? ANYPARAM : IGNOREPARAM;
  case BIStablescale: return SCALEPARAM;
  case BIStablecdiag : return VARPARAM;
  case BIStablerho :
    return S->rhored_given || origin == mle_conform ? IGNOREPARAM : ONLYRETURN;
  case BIStablerhored :
    return S->rhored_given || origin != original ? ANYPARAM : IGNOREPARAM;
  case BIStablealpha : return S->alphadiag_given || origin == mle_conform
      ? IGNOREPARAM : ONLYRETURN;
  default : BUG;
  }
}



sortsofparam sortof_bistable_INisOUT(model *cov, int k,
			     int  VARIABLE_IS_NOT_USED *row,
			     int  VARIABLE_IS_NOT_USED *col) {
  bistable_storage *S = cov->Sbistable;
  if (S == NULL) return UNKNOWNPARAM;
  switch(k) {
  case BIStablealphadiag : return !S->alphadiag_given ? ONLYMLE : ANYPARAM;
  case BIStablescale: return SCALEPARAM;
  case BIStablecdiag : return VARPARAM;
  case BIStablerho : return S->rhored_given ? ONLYRETURN : IGNOREPARAM;
  case BIStablerhored : return S->rhored_given ? ONLYMLE : ANYPARAM;
  case BIStablebetared : return !S->alphadiag_given ? ONLYMLE : ANYPARAM;
  case BIStablealpha : return !S->alphadiag_given ? ONLYRETURN : IGNOREPARAM;
  default : BUG;
  }
}


void rangebiStable(model VARIABLE_IS_NOT_USED *cov, range_type *range){


  range->min[BIStablealpha] = 0.0;
  range->max[BIStablealpha] = 1;
  range->pmin[BIStablealpha] = 0.06;
  range->pmax[BIStablealpha] = 1;
  range->openmin[BIStablealpha] = true;
  range->openmax[BIStablealpha] = false;


  range->min[BIStablealphadiag] = 0.0;
  range->max[BIStablealphadiag] = 1;
  range->pmin[BIStablealphadiag] = 0.06;
  range->pmax[BIStablealphadiag] = 1;
  range->openmin[BIStablealphadiag] = true;
  range->openmax[BIStablealphadiag] = false;

  range->min[BIStablebetared] = 0.0;
  range->max[BIStablebetared] = 1.0;
  range->pmin[BIStablebetared] = 0.0;
  range->pmax[BIStablebetared] =1.0;
  range->openmin[BIStablebetared] = false;
  range->openmax[BIStablebetared] = false;

  range->min[BIStablescale] = 0.0;
  range->max[BIStablescale] = RF_INF;
  range->pmin[BIStablescale] = 1e-2;
  range->pmax[BIStablescale] = 100.0;
  range->openmin[BIStablescale] = true;
  range->openmax[BIStablescale] = true;

  //    *c = P0(BIc];
  // to do: check rhos
  
  //  int dim = cov->tsdim;
  //  double *scale = P(BIStablescale);

  range->min[BIStablerho] = - 1.0;
  range->max[BIStablerho] = 1.0;
  range->pmin[BIStablerho] = - 1.0;
  range->pmax[BIStablerho] = 1.0;
  range->openmin[BIStablerho] = false;
  range->openmax[BIStablerho] = false;

  range->min[BIStablerhored] = - 1.0;
  range->max[BIStablerhored] = 1.0;
  range->pmin[BIStablerhored] = - 1.0;
  range->pmax[BIStablerhored] = 1.0;
  range->openmin[BIStablerhored] = false;
  range->openmax[BIStablerhored] = false;

  range->min[BIStablecdiag] = 0;
  range->max[BIStablecdiag] = RF_INF;
  range->pmin[BIStablecdiag] = 0.001;
  range->pmax[BIStablecdiag] = 1000;
  range->openmin[BIStablecdiag] = false;
  range->openmax[BIStablecdiag] = true;

}

void biStablePolynome(double r, double alpha, double a, int dim, double *v ) {
  double x = POW(a*r, alpha);
  
  if (dim == 1) {
    *v = alpha*x - alpha + 1;
  }
  if ( dim == 2 || dim == 3 ) {
    *v = alpha*alpha*x*x - 3*alpha*alpha*x + 4*alpha*x + alpha*alpha -
        4*alpha + 3;
  }
  
}

void biStableUnderInfLog(double r, double *alpha, double *a,
                         int dim, double *res ) {
  
  double x = POW(a[i11]*r, alpha[i11]),
    y = POW(a[i21]*r, alpha[i21]),
    z = POW(a[i22]*r, alpha[i22]);
  double p1, p2, p3;
  biStablePolynome(r, alpha[i11], a[i11], dim, &p1 );
  biStablePolynome(r, alpha[i21], a[i21], dim, &p2 );
  biStablePolynome(r, alpha[i22], a[i22], dim, &p3 );
  
  if (r == 0) {
    *res = 0;
  }
  else {
          *res =  (alpha[i11] + alpha[i22] - 2*alpha[i21])*LOG(r) +
        (2*y - x - z) + LOG(p1*p3/(p2*p2));
  }
  
}

void biStableInterval(double *alpha, double *a, int dim,
                      double *left, double *right ) {
  double middle = 1,
      fmiddle, fright, fleft,
      scalesratio1, scalesratio2, firstmin = 1.0/epsilon;
  *left = *right = middle;
  
  // check first if POW(a[0]/a[1], alpha[0] > 11 or
  // check first if POW(a[2]/a[1], alpha[2] > 11
  
  scalesratio1 = POW(a[0]/a[1], alpha[0]);
  scalesratio2 = POW(a[2]/a[1], alpha[2]);
   
  if ( (scalesratio1 >= 11) || (scalesratio2 >= 11) ) {
    biStableUnderInfLog( 1.0/(POW(2, 1/alpha[1])*a[1] ), alpha, a, dim, &firstmin);
    if (EXP(firstmin) < epsilon ) {
      *left = 0;
      *right = 0;
      return;
    }
  }
  
  biStableUnderInfLog(middle, alpha, a, dim, &fmiddle);
  
  
  
  // if the value at 1/(POW(2, 1/alpha[1])*a[1] ) is smaller that the value at 1
  // start searching around 1/(POW(2, 1/alpha[1])*a[1] )
  if ( fmiddle > firstmin ) {
    middle = 1.0/(POW(2, 1/alpha[1])*a[1] );
    *left = *right = middle;
    fmiddle = firstmin;
    
  }
  fright = fleft = fmiddle;
  
  while (( fmiddle >= MIN(fleft, fright))  &&
         (EXP(MIN(fleft, MIN(fright, fmiddle))) > epsilon ) )  {
    
    if ( fleft <= fmiddle )  {
      middle = *left;
      fmiddle = fleft;
      *left = *left/2;
    }
    
    if (fright <= fmiddle) {
      middle = *right;
      fmiddle = fright;
      *right = *right*2;
    }
    
    biStableUnderInfLog(*right, alpha, a, dim, &fright );
    biStableUnderInfLog(*left, alpha, a, dim,  &fleft );
  }
  
  
  if (EXP(MIN(fleft, MIN(fright, fmiddle))) <= epsilon ) {
    *left = 0;
    *right = 0;
  }
}

/*
 * Golden search algorithm from Numerical Recipes in C,
 * book by William H. Press
 *
 * (c, c)  - interval for searching
 *  *alpha, *a - array of parameters of biStable model
 *  *rhomax - maximum allowable rho for *alpha, *a, the result of
 * the searching
 * dim - dimension of the model
 * */


#define GOLDENR 0.61803399
#define GOLDENC (1.0 -GOLDENR)
#define GOLDENTOL 1e-9

#define GOLDENSHIFT2(a, b, c) (a) = (b); (b) = (c);
#define GOLDENSHIFT3(a, b, c, d) (a) = (b); (b) = (c); (c) = (d);
void biStableMinRho(model *cov, double *alpha, double *a, double ax, double cx, double *rhomax) {
  
  double bx = ax + (cx - ax)*GOLDENC,
  f1, f2, x0, x1, x2, x3, dummy,
    logconst = LOG(alpha[i11]*alpha[i22]/(alpha[i21]*alpha[i21])*
                             POW(a[i11], alpha[i11])*POW(a[i22], alpha[i22])/
                            POW(a[i21], 2*alpha[i21]));
  x0 = ax;
  x3 = cx;
  if ( FABS(cx - bx) > FABS(bx - ax) ) {
    x1 = bx;
    x2 = bx + GOLDENC*(cx - bx);
  } else {
    x2 = bx;
    x1 = bx - GOLDENC*(bx - ax);
  }
  int dim = OWNLOGDIM(0);
  biStableUnderInfLog(x1, alpha, a, dim, &f1);
  biStableUnderInfLog(x2, alpha, a, dim, &f2);
  
  while ( FABS(x3 - x0) > GOLDENTOL*(FABS(x1) + FABS(x2) ) ) {
    if (f2 < f1) {
      GOLDENSHIFT3(x0, x1, x2, GOLDENR*x1 + GOLDENC*x3 )
          biStableUnderInfLog(x2, alpha, a, dim, &dummy);
      GOLDENSHIFT2(f1, f2, dummy )
    } else {
      GOLDENSHIFT3(x3, x2, x1, GOLDENR*x2 + GOLDENC*x0 )
          biStableUnderInfLog(x1, alpha, a, dim, &dummy);
      GOLDENSHIFT2(f2, f1, dummy )
    }
  }
  

  
  if (f1 < f2) {
    *rhomax  = SQRT(EXP(f1+logconst));
  } else {
    *rhomax  = SQRT(EXP(f2+logconst));
  }
  *rhomax = MIN(*rhomax, 1);
}


int initbiStable(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double a[3],
    *cdiag = P(BIStablecdiag),
    *rho = P(BIStablerho),
    *rhored = P(BIStablerhored),
    // *cdiag = P(BIStablecdiag),
    *alpha = P(BIStablealpha),
    betared = RF_NAN,
    *alphadiag,
    *scale = P(BIStablescale),
    betaa, betac, //beta_red in [0, 1]; beta_red*betaa + betac in [betac, betaa + betac]
    rhomax = -2,
    left = 0,
    right = 0;
  int dim = OWNLOGDIM(0);
  bool notallowed = false;
  bistable_storage *S = cov->Sbistable;
  assert(S != NULL);

  if (PisNULL(BIStablescale)) {
    PALLOC(BIStablescale, 3, 1);
    scale = P(BIStablescale);
    for (int i = 0; i < 3; scale[i++] = 1.0);
   }
 // set a = 1/s
  a[i11] = 1.0/scale[i11];
  a[i21] = 1.0/scale[i21];
  a[i22] = 1.0/scale[i22];


 // default variance is 1
  if (PisNULL(BIStablecdiag)) {
    PALLOC(BIStablecdiag, 2, 1);
    cdiag = P(BIStablecdiag);
    cdiag[0] = 1;
    cdiag[1] = 1; 
   }

 // alpha diagonal are set => we are in likelohood estimation
 // we map alphadiag and beta red to alpha[i11], alpha[i21], alpha[i22]
 // alpha[i21] = betared*betaa + betac

 // beta = betaa + (2 - betaa)*beta_red
  if (S->alphadiag_given) {
    alphadiag = P(BIStablealphadiag);
    betaa = MAX(alphadiag[0], alphadiag[1]);
    betac = 2 - betaa;
    if (!PisNULL(BIStablebetared)) {
      betared = P0(BIStablebetared);
    }

  if (s->check && !PisNULL(BIStablealpha) ) {
    // if both alphas and alphadiag are given and check is true,
    // then check if they alphas and alphadiag are equal (or close to equal)
    alpha = P(BIStablealpha);

    if (cov->nrow[BIStablealpha] != 3 || cov->ncol[BIStablealpha] != 1)
       QERRC(BIStablealpha, "must be a 3 x 1 vector");
    if (FABS(alpha[i11] - alphadiag[0]) > alpha[i11] * epsilon ||
        FABS(alpha[i22] - alphadiag[1]) > alpha[i22] * epsilon ||
        FABS(alpha[i21] - (betaa + betac*betared) )
        > alpha[i21] * epsilon) {
       //       PMI(cov);
      //printf("\n\n\n --- FABS(alpha[i21] - (betared*betaa + betac) )  = %10g \n",             FABS(alpha[i21] - (betaa + betared*betac) ) );
      //     printf("\n\n\n --- alpha[i21] * epsilon  = %10g \n",             alpha[i21] * epsilon );


      QERRC2(BIStablealpha, "does not match '%.50s' and '%.50s'.",
             BIStablealphadiag, BIStablebetared);
      }
   }

   // allocate memory for the alphas
  if (PisNULL(BIStablealpha)) {
     PALLOC(BIStablealpha, 3, 1);
    }
    alpha = P(BIStablealpha);

   // fill the alphas
   alpha[i11] = alphadiag[0];
   alpha[i22] = alphadiag[1];
   alpha[i21] = betaa + betared*betac;

   // These combinations of smoothness parameters are not allowed
   // if user mode, return an Error and stop
   // if likelihood mode, set correlation to 0

   //notallowed = ( alpha[i21] < MAX(alpha[i11], alpha[i22]) ) ||

  } else  {
   //alphas are given

   if ( !PisNULL(BIStablealpha) ) {
     // all alphas are given, s->check => parameters are set by user,
     // not likelihood,  or they are NA
     // we map alphas to alphadiag and betared

     if (PisNULL(BIStablealphadiag)) PALLOC(BIStablealphadiag, 2, 1);
     alphadiag = P(BIStablealphadiag);

     if (PisNULL(BIStablebetared)) PALLOC(BIStablebetared, 1, 1);

     alphadiag[0] = alpha[i11];
     //      alpha[i21];
     alphadiag[1] = alpha[i22];

     betaa = MAX(alphadiag[0], alphadiag[1]);
     betac = 2 - MAX(alphadiag[0], alphadiag[1]);

     P(BIStablebetared)[0] = (alpha[i21] - betaa)/betac ;
     betared  = P(BIStablebetared)[0];

    }

   if (PisNULL(BIStablealpha)) {
     // PMI(cov);
     QERRC2(BIStablealpha, "either '%.50s' or '%.50s' must be set", BIStablealpha,
           BIStablealphadiag);
    }
 }
   notallowed =
     ( ( alpha[i11] == alpha[i21] ) &&  ( alpha[i22] == alpha[i21] ) &&
       ( POW(a[i21], alpha[i11]) < 0.5*POW(a[i11], alpha[i11])+0.5*
         POW(a[i22], alpha[i11]) ) ) ||
     ( ( alpha[i11] == alpha[i21] ) &&  ( alpha[i11] > alpha[i22] ) &&
       ( a[i21] <= POW(0.5, 1/alpha[i11])*a[i11] ) ) ||
     ( ( alpha[i22] == alpha[i21] ) &&  ( alpha[i22] > alpha[i11] ) &&
       ( a[i21] <= POW(0.5, 1/alpha[i22])*a[i22] ) )||
      ( alpha[i11] > alpha[i21] ) ||  ( alpha[i22] > alpha[i21] ) ;

   if(  notallowed && !s->check ) {
     {rhomax = 0;}
    }

   if(  notallowed && s->check ) {
     QERRC(BIStablealpha, "This combination of smoothness parameters is not allowed.");
    }

  // alphas are set

 // find maximum allowable rho
  
   if (rhomax != 0) {
   biStableInterval(alpha, a, dim, &left, &right );
   if ( (right == 0) && (left == 0)   ) {
     rhomax = 0;
    } else {
     //     PMI(cov);
     //     printf("\nDear Olga, please check the PMI output; if the output looks reasonable you should check your code in the next line.\n");
     biStableMinRho(cov, alpha, a, left, right, &rhomax);
     //     printf("done\n");
    }
  }

 // if rho is given, then it is set by user. Check if it is correct or not
  // ERR("Olga, there must be a call of type \"if (rhored_given)...\"");
  // if (!PisNULL(BIStablerhored) ){
    if ( S->rhored_given ){
         
   //if rhored is given, we are in the likelihood estimation
      if (PisNULL(BIStablerho)) PALLOC(BIStablerho, 1, 1);
      rho = P(BIStablerho);
      *rho = P(BIStablerho)[0] = (*rhored)*rhomax;
      
    } else if (!PisNULL(BIStablerho) && s->check  ) {
      
      if ( PisNULL(BIStablerhored) ) PALLOC(BIStablerhored, 1, 1);
      rhored = P(BIStablerhored);
      *rhored = (*rho)/rhomax;
      
    } else if (!PisNULL(BIStablerho) && s->check  ) {
      if (FABS(*rho) > rhomax ) {
        QERRC(BIStablealpha,
              "The value of cross-correlation parameter rho is too high.");
      }
      
      if ( PisNULL(BIStablerhored) ) PALLOC(BIStablerhored, 1, 1);
      rhored = P(BIStablerhored);
      *rhored = P(BIStablerhored)[0]= (*rho)/rhomax;
      
  }
    
    // printf("end init: alphadiag_given = %d\n", S->alphadiag_given);
  cov->initialised = true;
  RETURN_NOERROR;
}

void coinitbiStable(model *cov, localinfotype *li) {
  double
      thres = 1,
      *alpha = P(BIStablealpha);

  if ( ( alpha[0] <= thres ) &&  ( alpha[1] <= thres ) &&
       ( alpha[2] <= thres ) ) {
      li->instances = 1;
      li->value[0] = 1; //  q[CUTOFF_A]
      li->msg[0] =MSGLOCAL_OK;
    }
  else {
      li->instances = 1;
      li->value[0] = CUTOFF_THIRD_CONDITION ; //  q[CUTOFF_A]
      li->msg[0] = MSGLOCAL_JUSTTRY;
    }
}

void biStable (double *x, model *cov, double *v) {
  int i;
  double
    alpha = P0(STABLE_ALPHA),
    *cdiag = P(BIStablecdiag),
    y = *x, z;
  
  assert(BIStablealpha == STABLE_ALPHA);
  
  /*   biStable_storage *S = cov->SbiStable;
    assert(S != NULL);
    */
  
  
  for (i=0; i<3; i++) {
    z = y/P(BIStablescale)[i];
    P(STABLE_ALPHA)[0] = P(BIStablealpha)[i];
    stable(&z, cov, v + i);
  }
  P(BIStablealpha)[0] = alpha;
  
  v[0] = v[0]*cdiag[0];
  v[3] = v[2]*cdiag[1];
  v[1] *= P0(BIStablerho)*SQRT(cdiag[0]*cdiag[1]);
  v[2] = v[1];
  
  // printf("(%4.3f, %4.3f; %4.3e %4.3e %4.3e %4.3e)\t", x[0], x[1], v[0], v[1], v[2], v[3]);
  //  PMI(cov->calling->calling);  
  // if (!R_FINITE(v[0]) || !R_FINITE(v[1]) || !R_FINITE(v[2]) || !R_FINITE(v[3])) { PMI(cov); printf("(%4.3f, %4.3f; %4.3e %4.3e %4.3e %4.3e)\t", x[0], x[1], v[0], v[1], v[2], v[3]); BUG; }
  
}



void DbiStable(double *x, model *cov, double *v) {
  int i;
  double
    alpha = P0(STABLE_ALPHA),
    *cdiag = P(BIStablecdiag),
    y = *x, z;
  assert(BIStablealpha == STABLE_ALPHA);
  
  /*  biStable_storage *S = cov->SbiStable;
    assert(S != NULL);
   */
  
  for (i=0; i<3; i++) {
    z = y/P(BIStablescale)[i];
    P(STABLE_ALPHA)[0] = P(BIStablealpha)[i];
    Dstable(&z, cov,  v + i);
    v[i] /= P(BIStablescale)[i];
  }
  P(BIStablealpha)[0] = alpha;
  
  v[0] = v[0]*cdiag[0];
  v[3] = v[2]*cdiag[1];
  v[1] *= P0(BIStablerho)*SQRT(cdiag[0]*cdiag[1]);
  v[2] = v[1];
}


void DDbiStable(double *x, model *cov, double *v) {
  int i;
  double
    alpha = P0(STABLE_ALPHA),
    *cdiag = P(BIStablecdiag),
    y = *x, z;
  assert(BIStablealpha == STABLE_ALPHA);

  /*  biStable_storage *S = cov->SbiStable;
    assert(S != NULL);
    */

  for (i=0; i<3; i++) {
      z = y/P(BIStablescale)[i];
      P(STABLE_ALPHA)[0] = P(BIStablealpha)[i];
      DDstable(&z, cov,  v + i);
      v[i] /= P(BIStablescale)[i]*P(BIStablescale)[i];
    }
  P(BIStablealpha)[0] = alpha;

  v[0] = v[0]*cdiag[0];
  v[3] = v[2]*cdiag[1];
  v[1] *= P0(BIStablerho)*SQRT(cdiag[0]*cdiag[1]);
  v[2] = v[1];
}
void D3biStable(double *x, model *cov, double *v) {
  int i;
  double
    alpha = P0(STABLE_ALPHA),
    *cdiag = P(BIStablecdiag),
    y = *x, z;
  assert(BIStablealpha == STABLE_ALPHA);

  /*   biStable_storage *S = cov->SbiStable;
    assert(S != NULL);
   */

  for (i=0; i<3; i++) {
      z = y/P(BIStablescale)[i];
      P(STABLE_ALPHA)[0] = P(BIStablealpha)[i];
      D3stable(&z, cov, v + i);
      v[i] /= P(BIStablescale)[i]*P(BIStablescale)[i]*P(BIStablescale)[i];
    }
  P(BIStablealpha)[0] = alpha;

  v[0] = v[0]*cdiag[0];
  v[3] = v[2]*cdiag[1];
  v[1] *= P0(BIStablerho)*SQRT(cdiag[0]*cdiag[1]);
  v[2] = v[1];
}
void D4biStable(double *x, model *cov, double *v) {
  int i;
  double
    alpha = P0(STABLE_ALPHA),
    *cdiag = P(BIStablecdiag),
    y = *x, z;
  //   assert(BIStablealpha == STABLE_ALPHA);
  
  /*    biStable_storage *S = cov->SbiStable;
    assert(S != NULL);
  */
  //    assert(cov->initialised);
  
  
  for (i=0; i<3; i++) {
    z = y/P(BIStablescale)[i];
    P(STABLE_ALPHA)[0] = P(BIStablealpha)[i];
    D4stable(&z, cov, v + i);
    double dummy =P(BIStablescale)[i]*P(BIStablescale)[i];
    v[i] /= dummy*dummy;
  }
  P(BIStablealpha)[0] = alpha;
  
  v[0] = v[0]*cdiag[0];
  v[3] = v[2]*cdiag[1];
  v[1] *= P0(BIStablerho)*SQRT(cdiag[0]*cdiag[1]);
  v[2] = v[1];
}


