/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nondomain models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc


 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2014 Martin Schlather

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



#include <math.h>
 
#include "RF.h"
#include "primitive.h"
#include "Operator.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "Coordinate_systems.h"


#define i11 0
#define i21 1
#define i22 2

#define epsilon 1e-15


double BesselUpperB[Nothing + 1] =
{80, 80, 80, // circulant
 80, RF_INF, // TBM
 80, 80,    // direct & sequ
 RF_NA, RF_NA, RF_NA,  // GMRF, ave, nugget
 RF_NA, RF_NA, RF_NA,  // exotic
 RF_INF   // Nothing
};


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




/* BCW */
#define BCW_ALPHA 0
#define BCW_BETA 1
#define BCW_EPS 1e-7
#define BCW_TAYLOR_ZETA \
  (- LOG2 * (1.0 + 0.5 * zetalog2 * (1.0 + ONETHIRD * zetalog2)))
#define BCW_CAUCHY pow(1.0 + pow(*x, alpha), zeta)
void bcw(double *x, cov_model *cov, double *v){
  double alpha = P0(BCW_ALPHA), beta=P0(BCW_BETA),
    zeta = beta / alpha,
    absZeta = fabs(zeta);
  if (absZeta > BCW_EPS) 
    *v = BCW_CAUCHY / (1.0 - pow(2.0, zeta));
  else {
    double dewijs = log(1.0 + pow(*x, alpha)),
      zetadewijs = zeta * dewijs,
      zetalog2 = zeta * LOG2
      ;
    if (fabs(zetadewijs) > BCW_EPS)
      *v = BCW_CAUCHY / (zeta * BCW_TAYLOR_ZETA);
    else {
      *v = dewijs * (1.0 + 0.5 * zetadewijs * (1.0 + ONETHIRD *zetadewijs))
	/ BCW_TAYLOR_ZETA;
    }
  } 
}


void Dbcw(double *x, cov_model *cov, double *v){
  double ha,
    alpha = P0(BCW_ALPHA), 
    beta=P0(BCW_BETA),
    zeta=beta / alpha,
    absZeta = fabs(zeta),
    y=*x;
  
  if (y ==0.0) {
    *v = ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? -INFTY : alpha); 
  } else {
    ha=pow(y, alpha - 1.0);
    *v = alpha * ha * pow(1.0 + ha * y, zeta - 1.0);
  }
  
  if (absZeta > BCW_EPS) *v *= zeta / (1.0 - pow(2.0, zeta));
  else {
    double zetalog2 = zeta * LOG2;
    *v /= BCW_TAYLOR_ZETA;
  }
}


void DDbcw(double *x, cov_model *cov, double *v){
  double ha,
    alpha = P0(BCW_ALPHA), 
    beta=P0(BCW_BETA), 
    zeta = beta / alpha,
    absZeta = fabs(zeta),
     y=*x;

  if (y == 0.0) {
    *v = ((alpha==2.0) ? -beta * (1.0 - beta) : INFTY); 
  } else {
    ha = pow(y, alpha);
    *v = -alpha * ha / (y * y) * (1.0 - alpha + (1.0 - beta) * ha)
      * pow(1.0 + ha, beta / alpha - 2.0);
  }
  if (absZeta > BCW_EPS) *v *= zeta / (1.0 - pow(2.0, zeta));
  else {
    double zetalog2 = zeta * LOG2;
    *v /= BCW_TAYLOR_ZETA;
  }
}


void Inversebcw(double *x, cov_model *cov, double *v) {  
  double 
    alpha = P0(BCW_ALPHA), 
    beta=P0(BCW_BETA),
    zeta = beta / alpha;
  if (*x == 0.0) {
    *v = beta < 0.0 ? RF_INF : 0.0; 
    return;
  }
  if (zeta != 0.0)
    *v = pow(pow(*x * fabs(pow(2.0, zeta) - 1.0) + 1.0, 1.0/zeta) - 1.0,
	       1.0/alpha); 
  else 
    *v =  pow(exp(*x * LOG2) - 1.0, 1.0 / alpha);   
}
int checkbcw(cov_model *cov) {
  double
    alpha = P0(BCW_ALPHA), 
    beta=P0(BCW_BETA);
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  cov->logspeed = beta > 0 ? RF_INF : beta < 0.0 ? 0 : alpha * INVLOG2;
  return NOERROR;
}

void rangebcw(cov_model *cov, range_type *range) {
  bool tcf = isTcf(cov->typus) || cov->isoown == SPHERICAL_ISOTROPIC;
  bool posdef = isPosDef(cov->typus);
  range->min[BCW_ALPHA] = 0.0;
  range->max[BCW_ALPHA] = tcf ? 1.0 : 2.0;
  range->pmin[BCW_ALPHA] = 0.05;
  range->pmax[BCW_ALPHA] = range->max[BCW_ALPHA];
  range->openmin[BCW_ALPHA] = true;
  range->openmax[BCW_ALPHA] = false;

  range->min[BCW_BETA] = RF_NEGINF;
  range->max[BCW_BETA] = posdef ? 0.0 : 2.0;
  range->pmin[BCW_BETA] = -10.0;
  range->pmax[BCW_BETA] = 2.0;
  range->openmin[BCW_BETA] = true;
  range->openmax[BCW_BETA] = posdef;
}


void coinitbcw(cov_model *cov, localinfotype *li) {
  double
    beta=P0(BCW_BETA);
  if (beta < 0) coinitgenCauchy(cov, li);
  else {
    li->instances = 0;
  }
}
void ieinitbcw(cov_model *cov, localinfotype *li) {
  if (P0(BCW_BETA) < 0) ieinitgenCauchy(cov, li);
  else {
    ieinitBrownian(cov, li); // to do: can nicht passen!!
    // todo: nachrechnen, dass OK!
  }
}



#define LOW_BESSELJ 1e-20
#define LOW_BESSELK 1e-20
#define BESSEL_NU 0
void Bessel(double *x, cov_model *cov, double *v){

  static double nuOld=RF_INF;
  static double gamma;
  double nu = P0(BESSEL_NU), 
    y=*x;
  if  (y <= LOW_BESSELJ) {*v = 1.0; return;}
  if (y == RF_INF)  {*v = 0.0; return;} // also for cosine i.e. nu=-1/2
  if (nuOld!=nu) {
    nuOld=nu;
    gamma = gammafn(nuOld+1.0);
  }
  *v = gamma  * pow(2.0 / y, nuOld) * bessel_j(y, nuOld);
}
int initBessel(cov_model *cov, gen_storage 
	       VARIABLE_IS_NOT_USED *s) {
  ASSERT_GAUSS_METHOD(SpectralTBM);
  return NOERROR;
}
int checkBessel(cov_model *cov) {
  // Whenever TBM3Bessel exists, add further check against too small nu! 
  double nu = P0(BESSEL_NU);
  int i;
  double dim = (2.0 * P0(BESSEL_NU) + 2.0);

  for (i=0; i<= Nothing; i++)
    cov->pref[i] *= (ISNAN(nu) || nu < BesselUpperB[i]);
  if (cov->tsdim>2) cov->pref[SpectralTBM] = PREF_NONE; // do to d > 2 !
  cov->maxdim = (ISNAN(dim) || dim >= INFDIM) ? INFDIM : (int) dim;

  return NOERROR;
}
void spectralBessel(cov_model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  double 
    nu =  P0(BESSEL_NU);
/* see Yaglom ! */
  // nu==0.0 ? 1.0 : // not allowed anymore;
	// other wise densityBessel (to define) will not work
  if (nu >= 0.0) {
    E12(s, cov->tsdim, nu > 0 ? sqrt(1.0 - pow(UNIFORM_RANDOM, 1 / nu)) : 1, e);
  } else {
    double A;
    assert(cov->tsdim == 1);
    if (nu == -0.5) A = 1.0;
    else { // siehe private/bessel.pdf, p. 6, Remarks
      // http://homepage.tudelft.nl/11r49/documents/wi4006/bessel.pdf
      while (true) {
	A = 1.0 - pow(UNIFORM_RANDOM, 1.0 / ( P0(BESSEL_NU) + 0.5));
	if (UNIFORM_RANDOM <= pow(1 + A, nu - 0.5)) break;
      }
    }
    E1(s, A, e);
  }
}
void rangeBessel(cov_model *cov, range_type *range){
  range->min[BESSEL_NU] = 0.5 * ((double) cov->tsdim - 2.0);
  range->max[BESSEL_NU] = RF_INF;
  range->pmin[BESSEL_NU] = 0.0001 + range->min[BESSEL_NU];
  range->pmax[BESSEL_NU] = range->pmin[BESSEL_NU] + 10.0;
  range->openmin[BESSEL_NU] = false;
  range->openmax[BESSEL_NU] = true;
}


/* Cauchy models */
#define CAUCHY_GAMMA 0
void Cauchy(double *x, cov_model *cov, double *v){
  double gamma = P0(CAUCHY_GAMMA);
  *v = pow(1.0 + *x * *x, -gamma);
}
void logCauchy(double *x, cov_model *cov, double *v, double *Sign){
  double gamma = P0(CAUCHY_GAMMA);
  *v = log(1.0 + *x * *x) * -gamma;
  *Sign = 1.0;
}
void TBM2Cauchy(double *x, cov_model *cov, double *v){
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
void DCauchy(double *x, cov_model *cov, double *v){
  double y=*x, gamma = P0(CAUCHY_GAMMA);
  *v = (-2.0 * gamma * y) * pow(1.0 + y * y, -gamma - 1.0);
}
void DDCauchy(double *x, cov_model *cov, double *v){
  double ha = *x * *x, gamma = P0(CAUCHY_GAMMA);
  *v = 2.0 * gamma * ((2.0 * gamma + 1.0) * ha - 1.0) * 
    pow(1.0 + ha, -gamma - 2.0);
}
void InverseCauchy(double*x, cov_model *cov, double *v){
  double
    gamma = P0(CAUCHY_GAMMA);
  if (*x == 0.0) *v = RF_INF;
  else *v = sqrt(pow(*x, -1.0 / gamma) - 1.0);
}
int checkCauchy(cov_model VARIABLE_IS_NOT_USED  *cov){
  return NOERROR;
}
void rangeCauchy(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[CAUCHY_GAMMA] = 0.0;
  range->max[CAUCHY_GAMMA] = RF_INF;
  range->pmin[CAUCHY_GAMMA] = 0.09;
  range->pmax[CAUCHY_GAMMA] = 10.0;
  range->openmin[CAUCHY_GAMMA] = true;
  range->openmax[CAUCHY_GAMMA] = true;
}
void coinitCauchy(cov_model VARIABLE_IS_NOT_USED *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0; //  q[CUTOFF_A] 
  li->msg[0] = MSGLOCAL_JUSTTRY;
}
void DrawMixCauchy(cov_model VARIABLE_IS_NOT_USED *cov, double *random) { //better GR 3.381.4 ?? !!!!
  // split into parts <1 and >1 ?>
  *random = -log(1.0 -UNIFORM_RANDOM);
}
double LogMixDensCauchy(double VARIABLE_IS_NOT_USED *x, double logV, cov_model *cov) {
  double gamma = P0(CAUCHY_GAMMA);
  // da dort 1/w vorgezogen 
  return 0.5 * (gamma - 1.0) * logV - 0.5 * lgammafn(gamma);
}


/* another Cauchy model */
#define CTBM_ALPHA 0
#define CTBM_BETA 1
#define CTBM_GAMMA 2
void Cauchytbm(double *x, cov_model *cov, double *v){
  double ha, 
    alpha = P0(CTBM_ALPHA), 
    beta = P0(CTBM_BETA), 
    gamma = P0(CTBM_GAMMA),
    y=*x;
  if (y==0) {
    *v = 1.0;
  } else {
    ha=pow(y, alpha);
    *v = (1.0 + (1.0 - beta / gamma) * ha) * pow(1.0 + ha, -beta / alpha - 1.0);
  }
}

void DCauchytbm(double *x, cov_model *cov, double *v){
  double y= *x, ha, 
    alpha = P0(CTBM_ALPHA), 
    beta = P0(CTBM_BETA),
    gamma = P0(CTBM_GAMMA);
  if (y == 0.0) *v = 0.0; // WRONG VALUE, but multiplied 
  else {                                  // by zero anyway
    ha = pow(y, alpha - 1.0);
    *v = beta *  ha * (-1.0 - alpha / gamma  + ha * y * (beta / gamma - 1.0)) *
      pow(1.0 + ha * y, -beta /alpha - 2.0);
  }
}


void rangeCauchytbm(cov_model *cov, range_type *range){
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
 
  range->min[CTBM_GAMMA] = (double) cov->tsdim;
  range->max[CTBM_GAMMA] = RF_INF;
  range->pmin[CTBM_GAMMA] = range->min[CTBM_GAMMA];
  range->pmax[CTBM_GAMMA] = range->pmin[CTBM_GAMMA] + 10.0;
  range->openmin[CTBM_GAMMA] = false;
  range->openmax[CTBM_GAMMA] = true;

}



/* circular model */
void circular(double *x, cov_model VARIABLE_IS_NOT_USED  *cov, double *v) {
  double y = *x;
  *v = (y >= 1.0) ? 0.0 
    : 1.0 - (2.0 * (y * sqrt(1.0 - y * y) + asin(y))) * INVPI;
}
void Dcircular(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x * *x;
  *v = (y >= 1.0) ? 0.0 : -4 * INVPI * sqrt(1.0 - y);
}
int structCircSph(cov_model *cov, cov_model **newmodel, int dim) { 
  ASSERT_NEWMODEL_NOT_NULL;

  switch (cov->role) {
  case ROLE_POISSON_GAUSS :
    {
      addModel(newmodel, BALL, cov);      
      addModel(newmodel, DOLLAR);
      addModelKappa(*newmodel, DSCALE, SCALESPHERICAL);
      kdefault((*newmodel)->kappasub[DSCALE], SPHERIC_SPACEDIM,
	       (double) cov->tsdim);
      kdefault((*newmodel)->kappasub[DSCALE], SPHERIC_BALLDIM, (double) dim);
     } 
    break;
  case ROLE_POISSON : 
    return addUnifModel(cov, 1.0, newmodel); // to do: felix
  case ROLE_MAXSTABLE : 
    return addUnifModel(cov, 1.0, newmodel);
  default:
    BUG;
  }
  return NOERROR;
}

int structcircular(cov_model *cov, cov_model **newmodel) {
  return structCircSph(cov, newmodel, 2);
}


// spatially constant (matrix); 
#define CONSTANT_M 0

void kappaconstant(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc) {
  *nr = *nc = i ==  CONSTANT_M ? 0 : -1;
}

void constant(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v){
  int 
    vdimSq = cov->vdim[0] * cov->vdim[1];
  MEMCOPY(v, P(CONSTANT_M), vdimSq * sizeof(double));
}

void nonstatconstant(double VARIABLE_IS_NOT_USED *x, double VARIABLE_IS_NOT_USED *y, cov_model *cov, double *v){
  int 
    vdimSq = cov->vdim[0] * cov->vdim[1];
  MEMCOPY(v, P(CONSTANT_M), vdimSq * sizeof(double));
}

int checkconstant(cov_model *cov) {
  int info, err; // anzahl listen elemente

 
  if (GLOBAL.internal.warn_constant) {
    GLOBAL.internal.warn_constant = false;
    warning("NOTE THAT THE DEFINITION OF 'RMconstant' HAS CHANGED. MAYBE, 'RMfixcov' IS THE CORRECT MODEL");
  }
 
  cov->vdim[0] = cov->vdim[1] = cov->nrow[CONSTANT_M];
  if (cov->vdim[0] != cov->vdim[1])  return ERROR_MATRIX_SQUARE;

  if (cov->typus == VariogramType) { // ACHTUNG wirklich genau VariogramType
    SERR("strange call");
  }

  
  if (cov->q != NULL) {
    return (int) cov->q[0];
  } else QALLOC(1);

  cov->q[0] = NOERROR;

  cov->monotone = COMPLETELY_MON;
  cov->ptwise_definite = pt_posdef;
   
  // check whether positive eigenvalue  
  long total = cov->nrow[CONSTANT_M] * cov->ncol[CONSTANT_M] * sizeof(double);
  double *dummy = (double*) MALLOC(total);
  MEMCOPY(dummy, P(CONSTANT_M), total);
      
  F77_CALL(dpofa)(dummy, cov->nrow + CONSTANT_M, cov->ncol + CONSTANT_M, 
		  &info); // cholesky
  
  FREE(dummy);
  if (info != 0) {
    if (isPosDef(cov)) return cov->q[0] = ERROR_MATRIX_POSDEF;
    cov->monotone = MONOTONE;
    cov->ptwise_definite = pt_indef;
  }
    

  cov->matrix_indep_of_x = true;
  cov->mpp.maxheights[0] = RF_NA;
  err = checkkappas(cov);

 
  return err;
}

void rangeconstant(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){

  range->min[CONSTANT_M] = RF_NEGINF;
  range->max[CONSTANT_M] = RF_INF;
  range->pmin[CONSTANT_M] = -1e10;
  range->pmax[CONSTANT_M] = 1e10;
  range->openmin[CONSTANT_M] = true;
  range->openmax[CONSTANT_M] = true;

 }






/* Covariate models */
#define COVARIATE_C 0
#define COVARIATE_X 1
#define COVARIATE_T 2
#define COVARIATE_GRID 3
#define COVARIATE_FACTOR 4
void kappa_covariate(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
    *nc = *nr = i <= COVARIATE_T ? 0 
      : i <= COVARIATE_FACTOR ? 1 
      : -1;
}

int cutidx(double Idx, int len) {
  int idx = (int) round(Idx);
  if (idx < 0) idx = 0;
  if (idx >= len) idx = len - 1;
  return idx;
}
void covariate(double *x, cov_model *cov, double *v){
  location_type *loc = Loc(cov);
  int nr, i, d, idx,
    vdim = cov->vdim[0],    
    ntot = loc->totalpoints;
  double X[2],
    *p = P(COVARIATE_C),
    f = P0(COVARIATE_FACTOR);

  assert(cov->vdim[1] == 1);
  if (loc->i_col == I_COL_NA) { 
    // interpolate: here nearest neighbour/voronoi
   int 
     cummul = 1.0,
     tsdim = cov->xdimprev;
   if (loc->grid) {
     nr = 0;
     
     for (d = 0; d<tsdim; d++) {
       int len = loc->xgr[d][XLENGTH];
       double
	 step = loc->xgr[d][XSTEP];      
          
       if (d > 1 || !isAnySpherical(cov->isoown)) {
	 idx = cutidx((x[d] - loc->xgr[d][XSTART]) / step, len);
       } else {

	 if (d == 0) { // to do: report technique?
	   double full, half, y[2];
	   int idx2,
	     dim = 2;
	   for (i=0; i<dim; i++) y[i] = loc->xgr[i][XSTART];
	   if (isSpherical(cov->isoown)) {
	     full = M_2_PI;
	     half = M_PI;
	     if (GLOBAL.coords.polar_coord) NotProgrammedYet("");
	   } else if (isEarth(cov->isoown)) {
	     full = 360.0;
	     half = 180.0;
	   } else BUG;
	   statmod2(y, full, half, X);

	   idx = cutidx((x[0] - X[0]) / step, len);
	   double X0 = X[0] + full * (2.0 * (double) (x[0] > X[0]) - 1);
	   idx2 = cutidx((x[0] - X0) / step, len);
	   if (fabs(x[0] - (X0 + idx2 * step)) <
	       fabs(x[0] - (X[0] + idx * step))) idx = idx2;
	 } else { 
	   assert(d==1);
	   idx = cutidx((x[d] - X[1]) / step, len);
	 }
       }
      
       nr += cummul * idx;
       cummul *= len;  
     } 
          
   } else { 
      cov_model *next = cov->sub[0];
      double distance,
	mindist = RF_INF,
	*given = loc->x;
      nr = 0;
      for (i=0; i<ntot; i++, given+=tsdim) {
	NONSTATCOV(x, given, next, &distance);
	if (distance < mindist) {
	  mindist = distance;
	  nr = i;
	}	
      }
    }
  } else { // should happen only in a matrix context!
    nr = loc->i_row;
  }
  
  if (GLOBAL.general.vdim_close_together) {
    p += nr * vdim;
    for (i=0; i<vdim; i++, p++) v[i] = *p * f;
  } else {
    p += nr;
    for (i=0; i<vdim; i++, p+=ntot) v[i] = *p * f;
  }
}

int checkcovariate(cov_model *cov){
  cov_model *next = cov->sub[0];
  int err, 
    data = cov->nrow[COVARIATE_C] * cov->ncol[COVARIATE_C];
  location_type *loc = Loc(cov);
  coord_sys_enum ccs = GLOBAL.coords.coord_system;

  //printf("%s\n",  ISONAMES[cov->isoown]);
  if ((ccs == cartesian && cov->isoown != CARTESIAN_COORD) ||
      (ccs == earth && cov->isoown != EARTH_COORD) ||
      (ccs == sphere && cov->isoown != SPHERICAL_COORD))
    SERR2("'%s' not the global coordinate system ('%s')",
	  ISONAMES[cov->isoown], COORD_SYS_NAMES[GLOBAL.coords.coord_system]);

  kdefault(cov, COVARIATE_FACTOR, 1.0);
  if ((!PisNULL(COVARIATE_X) || !PisNULL(COVARIATE_T)) && cov->ownloc == NULL) {
    if (cov->ncol[COVARIATE_T] != 3 && cov->ncol[COVARIATE_T] != 0)
      SERR1("'%s' should have three components", KNAME(COVARIATE_T));
    if (cov->ncol[COVARIATE_X] != loc->xdimOZ) 
      SERR1("number of columns of '%s' do not equal the dimension",
	    KNAME(COVARIATE_X));
    if ((err = loc_set(P(COVARIATE_X), NULL, P(COVARIATE_T), 
		       loc->spatialdim, loc->xdimOZ, 
		       cov->nrow[COVARIATE_X], 0, cov->ncol[COVARIATE_T]==3, 
		       PINT(COVARIATE_GRID), 
		       false, &(cov->ownloc))) != NOERROR) {
      if (cov->ownloc != NULL) BUG;
      return err;
    }
    loc = Loc(cov); // zwingend nochmal!
  } else {
    //   kdefault(cov, COVARIATE_X, 0);
    // kdefault(cov, COVARIATE_T, 0);
    // kdefault(cov, COVARIATE_GRID, 0);
  }

  int ntot = loc->totalpoints,
    vdim = data / ntot;
  if (vdim * ntot != data) 
    SERR("number of data not a multiple of the number of locations");
  cov->vdim[0] = vdim;
  cov->vdim[1] = 1;

  if (next == NULL) {
    QALLOC(1);
    addModel(cov, 0, TRAFO);
    next = cov->sub[0];
    kdefault(next, TRAFO_ISO, ISOTROPIC);
  }
 
  assert(cov->q != NULL);
  PARAMINT(next, TRAFO_ISO)[0] =
    isCartesian(cov->isoprev) ? ISOTROPIC
    : isSpherical(cov->isoprev) ? SPHERICAL_ISOTROPIC
    : isEarth(cov->isoprev) ? EARTH_ISOTROPIC 
    : ISO_MISMATCH;

   if ((err = CHECK(next, cov->tsdim,  cov->xdimown, ShapeType,
		    KERNEL, cov->isoown,
		   1, cov->role)) != NOERROR) return err;
  if (cov->vdim[0] != vdim ||  cov->vdim[1] != 1)
    SERR("mismatch on the number of multivariability");

  if (vdim == 1) {
    int i;
    double *c = P(COVARIATE_C);
    for (i=0; i<ntot; i++) if (c[i] != 0.0) break;
    cov->ptwise_definite = pt_zero;
    if (i<ntot) {
      if (c[i] > 0.0) {
	cov->ptwise_definite = pt_posdef;
        for (i++; i<ntot; i++)
	  if (c[i] < 0.0) {
	    cov->ptwise_definite = pt_indef;
	    break;
	  }
      } else { // < 0.0
	cov->ptwise_definite = pt_negdef;
        for (i++; i<ntot; i++)
	  if (c[i] > 0.0) {
	    cov->ptwise_definite = pt_indef;
	    break;
	  }
      }
    }
  } else cov->ptwise_definite = pt_unknown; // to do ?! alle Matrizen ueberpruefen...

  if (Loc(cov)->totalpoints <= 0) SERR("no locations given");
  
  EXTRA_STORAGE;
  return NOERROR;
}

void rangecovariate(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[COVARIATE_C] = RF_NEGINF;
  range->max[COVARIATE_C] = RF_INF;
  range->pmin[COVARIATE_C] = 1e10;
  range->pmax[COVARIATE_C] = -1e10;
  range->openmin[COVARIATE_C] = true;
  range->openmax[COVARIATE_C] = true;

  range->min[COVARIATE_X] = RF_NEGINF;
  range->max[COVARIATE_X] = RF_INF;
  range->pmin[COVARIATE_X] = 1e10;
  range->pmax[COVARIATE_X] = -1e10;
  range->openmin[COVARIATE_X] = true;
  range->openmax[COVARIATE_X] = true;

  range->min[COVARIATE_T] = RF_NEGINF;
  range->max[COVARIATE_T] = RF_INF;
  range->pmin[COVARIATE_T] = 1e10;
  range->pmax[COVARIATE_T] = -1e10;
  range->openmin[COVARIATE_T] = true;
  range->openmax[COVARIATE_T] = true;

  range->min[COVARIATE_GRID] = 0;
  range->max[COVARIATE_GRID] = 1;
  range->pmin[COVARIATE_GRID] = 0;
  range->pmax[COVARIATE_GRID] = 1;
  range->openmin[COVARIATE_GRID] = false;
  range->openmax[COVARIATE_GRID] = false;

  range->min[COVARIATE_FACTOR] = RF_NEGINF;
  range->max[COVARIATE_FACTOR] = RF_INF;
  range->pmin[COVARIATE_FACTOR] = 1e10;
  range->pmax[COVARIATE_FACTOR] = -1e10;
  range->openmin[COVARIATE_FACTOR] = true;
  range->openmax[COVARIATE_FACTOR] = true;
}





/* coxgauss, cmp with nsst1 !! */
// see Gneiting.cc

/* cubic */
void cubic(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y=*x, y2=y * y;
  *v = (y >= 1.0) ? 0.0
    : (1.0 + (((0.75 * y2 - 3.5) * y2 + 8.75) * y - 7) * y2);
}
void Dcubic(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) { 
  double y=*x, y2=y * y;
  *v = (y >= 1.0) ? 0.0 : y * (-14.0 + y * (26.25 + y2 * (-17.5 + 5.25 * y2)));
}

/* cutoff */
// see operator.cc

/* dagum */
#define DAGUM_BETA 0
#define DAGUM_GAMMA 1
#define DAGUM_BETAGAMMA 2
void dagum(double *x, cov_model *cov, double *v){
  double gamma = P0(DAGUM_GAMMA), 
    beta=P0(DAGUM_BETA);
  *v = 1.0 - pow((1 + pow(*x, -beta)), -gamma/beta);
}
void Inversedagum(double *x, cov_model *cov, double *v){ 
  double gamma = P0(DAGUM_GAMMA), 
    beta=P0(DAGUM_BETA);
    *v = *x == 0.0 ? RF_INF 
      : pow(pow(1.0 - *x, - beta / gamma ) - 1.0, 1.0 / beta);
} 
void Ddagum(double *x, cov_model *cov, double *v){
  double y=*x, xd, 
    gamma = P0(DAGUM_GAMMA), 
    beta=P0(DAGUM_BETA);
  xd = pow(y, -beta);
  *v = -gamma * xd / y * pow(1 + xd, -gamma/ beta -1);
}

int checkdagum(cov_model *cov){
  if (PisNULL(DAGUM_GAMMA) || PisNULL(DAGUM_BETA))
    SERR("parameters are not given all");
  double
    gamma = P0(DAGUM_GAMMA), 
    beta = P0(DAGUM_BETA);
  kdefault(cov, DAGUM_BETAGAMMA, beta/gamma);
  FIRST_INIT(initdagum);
 

  cov->monotone =  beta >= gamma ? MONOTONE 
    : gamma <= 1.0 ? COMPLETELY_MON
    : gamma <= 2.0 ? NORMAL_MIXTURE : MISMATCH;
    
  return NOERROR;
}

int initdagum(cov_model *cov, gen_storage *stor) {  
  double
    gamma = P0(DAGUM_GAMMA), 
    beta = P0(DAGUM_BETA),
    betagamma = P0(DAGUM_BETAGAMMA);
  
  if (stor->check) {
    bool tcf = isTcf(cov->typus) || cov->isoown == SPHERICAL_ISOTROPIC;
    bool isna_gamma = tcf && ISNA(gamma);
    if (isna_gamma) {
      if (cov->q == NULL) QALLOC(1); // just as a flag
    } else {
      P(DAGUM_BETAGAMMA)[0] = 1.0;
      // dummy value just to be in the range !
    }
  } else {     
    bool isna_gamma = cov->q != NULL;
    if (isna_gamma) P(DAGUM_GAMMA)[0] = beta / betagamma;
  }
  return NOERROR;
}

sortsofparam paramtype_dagum(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  return k==DAGUM_GAMMA ? IGNOREPARAM : ANYPARAM;
}
 

void rangedagum(cov_model *cov, range_type *range){
  bool tcf = isTcf(cov->typus) || cov->isoown == SPHERICAL_ISOTROPIC;
  range->min[DAGUM_BETA] = 0.0;
  range->max[DAGUM_BETA] = 1.0;
  range->pmin[DAGUM_BETA] = 0.01;
  range->pmax[DAGUM_BETA] = 1.0;
  range->openmin[DAGUM_BETA] = true;
  range->openmax[DAGUM_BETA] = false;

  range->min[DAGUM_GAMMA] = 0.0;
  range->max[DAGUM_GAMMA] = tcf ? 1.0 : 2.0;
  range->pmin[DAGUM_GAMMA] = 0.01;
  range->pmax[DAGUM_GAMMA] = range->max[DAGUM_GAMMA];
  range->openmin[DAGUM_GAMMA] = true;
  range->openmax[DAGUM_GAMMA] = tcf;

  range->min[DAGUM_BETAGAMMA] = 0.0;
  range->max[DAGUM_BETAGAMMA] = tcf ? 1.0 : RF_INF;
  range->pmin[DAGUM_BETAGAMMA] = 0.01;
  range->pmax[DAGUM_BETAGAMMA] = tcf ? 1.0 : 20.0;
  range->openmin[DAGUM_BETAGAMMA] = true;
  range->openmax[DAGUM_BETAGAMMA] = true;
}


/*  damped cosine -- derivative of e xponential:*/
#define DC_LAMBDA 0
void dampedcosine(double *x, cov_model *cov, double *v){
  double y = *x, lambda = P0(DC_LAMBDA);
  *v = (y == RF_INF) ? 0.0 : exp(-y * lambda) * cos(y);
}
void logdampedcosine(double *x, cov_model *cov, double *v, double *Sign){
  double 
    y = *x, 
    lambda = P0(DC_LAMBDA);
  if (y==RF_INF) {
    *v = RF_NEGINF;
    *Sign = 0.0;
  } else {
    double cosy=cos(y);
    *v = -y * lambda + log(fabs(cosy));
    *Sign = cosy > 0.0 ? 1.0 : cosy < 0.0 ? -1.0 : 0.0;
  }
}
void Inversedampedcosine(double *x, cov_model *cov, double *v){ 
  Inverseexponential(x, cov, v);
} 
void Ddampedcosine(double *x, cov_model *cov, double *v){
  double y = *x, lambda = P0(DC_LAMBDA);
  *v = - exp(-lambda*y) * (lambda * cos(y) + sin(y));
}
int checkdampedcosine(cov_model *cov){
   cov->maxdim = ISNAN(P0(DC_LAMBDA)) 
     ? INFDIM : (int) (PIHALF / atan(1.0 / P0(DC_LAMBDA)));
  return NOERROR;
}
void rangedampedcosine(cov_model *cov, range_type *range){
  range->min[DC_LAMBDA]  = 1.0 / tan(PIHALF / cov->tsdim); // ??  1.0 / tan(0.5 * PIHALF / cov->tsdim) ?? to do
  range->max[DC_LAMBDA] = RF_INF;
  range->pmin[DC_LAMBDA] = range->min[DC_LAMBDA] + 1e-10;
  range->pmax[DC_LAMBDA] = 10;
  range->openmin[DC_LAMBDA] = false;
  range->openmax[DC_LAMBDA] = true;
}


/* De Wijsian */
#define DEW_ALPHA 0 // for both dewijsian models
void dewijsian(double *x, cov_model *cov, double *v){
  double alpha = P0(DEW_ALPHA);
  *v = -log(1.0 + pow(*x, alpha));
}
void Ddewijsian(double *x, cov_model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
    p =pow(*x, alpha - 1.0) ;
  *v = - alpha * p / (1.0 + *x * p);
}
void DDdewijsian(double *x, cov_model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
    p = pow(*x, alpha - 2.0),
    ha = p * *x * *x,
    haP1 = ha + 1.0;
  *v = alpha * p * (1.0 - alpha + ha) / (haP1 * haP1);
}

void Inversedewijsian(double *x, cov_model *cov, double *v){ 
  double alpha = P0(DEW_ALPHA);
  *v = pow(exp(*x) - 1.0, 1.0 / alpha);    
} 
int checkdewijsian(cov_model *cov){
  double alpha = P0(DEW_ALPHA);
  cov->logspeed = alpha;
  return NOERROR;
}

void rangedewijsian(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[DEW_ALPHA] = 0.0;
  range->max[DEW_ALPHA] = 2.0;
  range->pmin[DEW_ALPHA] = UNIT_EPSILON;
  range->pmax[DEW_ALPHA] = 2.0;
  range->openmin[DEW_ALPHA] = true;
  range->openmax[DEW_ALPHA] = false; 
}


/* De Wijsian B */
#define DEW_RANGE 1
void DeWijsian(double *x, cov_model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
    range = P0(DEW_RANGE);
  *v = (*x >= range) ? 0.0 : 1.0 -
    log(1.0 + pow(*x, alpha)) / log(1.0 + pow(range, alpha));
}

void InverseDeWijsian(double *x, cov_model *cov, double *v){ 
  double alpha = P0(DEW_ALPHA),
    range = P0(DEW_RANGE);
  *v = *x >= 1.0 ? 0.0 
    : pow(pow(1.0 + pow(range, alpha),  1.0 - *x) - 1.0, 1.0 / alpha);
} 

void rangeDeWijsian(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[DEW_ALPHA] = 0.0;
  range->max[DEW_ALPHA] = 1.0;
  range->pmin[DEW_ALPHA] = UNIT_EPSILON;
  range->pmax[DEW_ALPHA] = 1.0;
  range->openmin[DEW_ALPHA] = true;
  range->openmax[DEW_ALPHA] = false; 

  range->min[DEW_RANGE] = 0.0;
  range->max[DEW_RANGE] = RF_INF;
  range->pmin[DEW_RANGE] = UNIT_EPSILON;
  range->pmax[DEW_RANGE] = 1000;
  range->openmin[DEW_RANGE] = true;
  range->openmax[DEW_RANGE] = true; 
}


/* exponential model */
void exponential(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
   *v = exp(- *x);
 }
void logexponential(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v, double *Sign){
  *v = - *x;
  *Sign = 1.0;
 }
void TBM2exponential(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) 
{
  double y = *x;
  *v = (y==0.0) ?  1.0 : 1.0 - PIHALF * y * I0mL0(y);
}
void Dexponential(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  *v = - exp(- *x);
}
void DDexponential(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
 *v = exp(-*x);
}
void Inverseexponential(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  *v = (*x == 0.0) ? RF_INF : -log(*x);
}

void nonstatLogInvExp(double *x, cov_model *cov,
		      double *left, double *right){
  double 
    z = *x <= 0.0 ? - *x : 0.0;  
  int d,
    dim = cov->tsdim;
 
  for (d=0; d<dim; d++) {
    left[d] = -z;
    right[d] = z;
  }
}


double densityexponential(double *x, cov_model *cov) {
  // to do: ersetzen durch die Familien 

  // spectral density
  int d,
    dim=cov->tsdim;
  double x2 = 0.0,
    dim12 = 0.5 * ((double) (dim + 1));
  for (d=0; d<dim; d++) x2 += x[d] * x[d];
  
  return gammafn(dim12) * pow(M_PI * (1.0 + x2), -dim12);
}

#define HAS_SPECTRAL_ROLE(cov)\
  cov->role == ROLE_GAUSS && cov->method==SpectralTBM

int initexponential(cov_model *cov, gen_storage *s) {
  int dim = cov->tsdim;
  double D = (double) dim;

  if (HAS_SPECTRAL_ROLE(cov)) {
    spec_properties *cs = &(s->spec);
    if (cov->tsdim <= 2) return NOERROR;
    cs->density = densityexponential;
    return search_metropolis(cov, s);
  }
  
  else if (hasAnyShapeRole(cov)) {
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
	i[xi] = dimfak * 2 * pow(M_PI / (double) (xi*xi), dimHalf) / 
	  gammafn(dimHalf);
	
	if (R < RF_INF) {
	  // Gradstein 3.351 fuer n=d-1
	  
	  //printf("here\n");
	  for (sum = 1.0, factor=1.0, d=1; d<dim; d++) {
	    factor *= R / D;
	    sum += factor;
	    //	printf("%d %f %f %f\n", d, sum, factor, R);
	}
	  sum *= dimfak * exp(-R);
	  // printf("%d %f %f %f\n", d, sum, factor, R);
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
  
  else ILLEGAL_ROLE;

  return NOERROR;
}
void do_exp(cov_model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 
}

void spectralexponential(cov_model *cov, gen_storage *S, double *e) {
  spectral_storage *s = &(S->Sspectral);
  if (cov->tsdim <= 2) {
    double A = 1.0 - UNIFORM_RANDOM;
    E12(s, cov->tsdim, sqrt(1.0 / (A * A) - 1.0), e);
  } else {
    metropolis(cov, S, e);
  }
}

int checkexponential(cov_model *cov) {
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (cov->tsdim != 2) cov->pref[Hyperplane] = PREF_NONE;
  return NOERROR;
}

int hyperexponential(double radius, double *center, double *rx,
		     cov_model VARIABLE_IS_NOT_USED *cov, bool simulate, 
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
    err;
  
  assert(cov->tsdim==2);

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
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }  
  
  /* creating the lines; some of the lines are not relevant as they
     do not intersect the rectangle investigated --> k!=4
     (it is checked if all the corners of the rectangle are on one 
     side (?) )
  */
   for(i=0; i<p; i++) {
    phi = UNIFORM_RANDOM * TWOPI;
    hx[q] = cos(phi);     hy[q] = sin(phi);    
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
  return -err;
}

void coinitExp(cov_model VARIABLE_IS_NOT_USED *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0; //  q[CUTOFF_A] 
  li->msg[0] = MSGLOCAL_OK;
}
void ieinitExp(cov_model VARIABLE_IS_NOT_USED *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = 1.0;
  li->msg[0] = MSGLOCAL_OK;
}
void DrawMixExp(cov_model VARIABLE_IS_NOT_USED *cov, double *random) {
  // GR 3.325: int_-infty^infty exp(x^2/2 - b/x^2) = exp(-\sqrt b)
  double x = GAUSS_RANDOM(1.0);
  *random = 1.0 / (x * x);
}
double LogMixDensExp(double VARIABLE_IS_NOT_USED *x, double VARIABLE_IS_NOT_USED logV, cov_model VARIABLE_IS_NOT_USED *cov) {
  // todo: check use of mixdens --- likely to programme it now completely differently 
  return 0.0;
}


// Brownian motion 
void fractalBrownian(double *x, cov_model *cov, double *v) {
  double alpha = P0(BROWN_ALPHA);
  *v = - pow(*x, alpha);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}

void logfractalBrownian(double *x, cov_model *cov, double *v, double *Sign) {
  double alpha = P0(BROWN_ALPHA);
  *v = log(*x) * alpha;//this is an invalid covariance function!
  *Sign = -1.0;
  // keep definition such that the value at the origin is 0
}
/* fractalBrownian: first derivative at t=1 */
void DfractalBrownian(double *x, cov_model *cov, double *v) 
{// FALSE VALUE FOR *x==0 and  alpha < 1
  double alpha = P0(BROWN_ALPHA);
  *v = (*x != 0.0) ? -alpha * pow(*x, alpha - 1.0)
    : alpha > 1.0 ? 0.0 
    : alpha < 1.0 ? RF_NEGINF
    : -1.0;
}
/* fractalBrownian: second derivative at t=1 */
void DDfractalBrownian(double *x, cov_model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  double alpha = P0(BROWN_ALPHA);
  *v = (alpha == 1.0) ? 0.0
    : (*x != 0.0) ? -alpha * (alpha - 1.0) * pow(*x, alpha - 2.0)
    : alpha < 1.0 ? RF_INF 
    : alpha < 2.0 ? RF_NEGINF 
    : -2.0;
}
void D3fractalBrownian(double *x, cov_model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  double alpha = P0(BROWN_ALPHA);
  *v = alpha == 1.0 || alpha == 2.0 ? 0.0 
    : (*x != 0.0) ? -alpha * (alpha - 1.0) * (alpha - 2.0) * pow(*x, alpha-3.0)
    : alpha < 1.0 ? RF_NEGINF 
    : RF_INF;
}

void D4fractalBrownian(double *x, cov_model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  double alpha = P0(BROWN_ALPHA);
  *v = alpha == 1.0 || alpha == 2.0 ? 0.0 
    : (*x != 0.0) ? -alpha * (alpha - 1.0) * (alpha - 2.0) * (alpha - 3.0) *
                     pow(*x, alpha-4.0)
    : alpha < 1.0 ? RF_INF 
    : RF_NEGINF;
}
int checkfractalBrownian(cov_model *cov){
  double alpha = P0(BROWN_ALPHA);
  cov->logspeed = RF_INF;
  cov->full_derivs = alpha <= 1.0 ? 0 : alpha < 2.0 ? 1 : cov->rese_derivs;
  cov->taylor[0][TaylorPow] = cov->tail[0][TaylorPow] = alpha;
  return NOERROR;
}

int initfractalBrownian(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double alpha = P0(BROWN_ALPHA);
  cov->taylor[0][TaylorPow] = cov->tail[0][TaylorPow] = alpha;  
  return NOERROR;
}

void InversefractalBrownian(double *x, cov_model *cov, double *v) {
  double alpha = P0(BROWN_ALPHA);
  *v = - pow(*x, 1.0 / alpha);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
void rangefractalBrownian(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[BROWN_ALPHA] = 0.0;
  range->max[BROWN_ALPHA] = 2.0;
  range->pmin[BROWN_ALPHA] = UNIT_EPSILON;
  range->pmax[BROWN_ALPHA] = 2.0;
  range->openmin[BROWN_ALPHA] = true;
  range->openmax[BROWN_ALPHA] = false;
}
void ieinitBrownian(cov_model *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = (cov->tsdim <= 2)
    ? ((P0(BROWN_ALPHA) <= 1.5) ? 1.0 : 2.0)
    : ((P0(BROWN_ALPHA) <= 1.0) ? 1.0 : 2.0);
  li->msg[0] = cov->tsdim <= 3 ? MSGLOCAL_OK : MSGLOCAL_JUSTTRY;
}



/* FD model */
#define FD_ALPHA 0
void FD(double *x, cov_model *cov, double *v){
  double alpha = P0(FD_ALPHA), y, d, k, skP1;
  static double dold=RF_INF;
  static double kold, sk;
  d = alpha * 0.5;
  y = *x;
  if (y == RF_INF) {*v = 0.0; return;} 
  k = trunc(y);
  if (dold!=d || kold > k) {
    sk = 1;
    kold = 0.0;
  }
  // Sign (-1)^k is (kold+d), 16.11.03, checked. 
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
void rangeFD(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[FD_ALPHA] = -1.0;
  range->max[FD_ALPHA] = 1.0;
  range->pmin[FD_ALPHA] = range->min[FD_ALPHA] + UNIT_EPSILON;
  range->pmax[FD_ALPHA] = range->max[FD_ALPHA] - UNIT_EPSILON;
  range->openmin[FD_ALPHA] = false;
  range->openmax[FD_ALPHA] = true;
}



/* fractgauss */
#define FG_ALPHA 0
void fractGauss(double *x, cov_model *cov, double *v){
  double y = *x, alpha = P0(FG_ALPHA);
  *v = (y == 0.0) ? 1.0 :  (y==RF_INF) ? 0.0 : 
    0.5 *(pow(y + 1.0, alpha) - 2.0 * pow(y, alpha) + pow(fabs(y - 1.0),alpha));
}
void rangefractGauss(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[FG_ALPHA] = 0.0;
  range->max[FG_ALPHA] = 2.0;
  range->pmin[FG_ALPHA] = UNIT_EPSILON;
  range->pmax[FG_ALPHA] = 2.0;
  range->openmin[FG_ALPHA] = true;
  range->openmax[FG_ALPHA] = false;
}


/* Gausian model */
void Gauss(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = exp(- *x * *x);
}
void logGauss(double *x, cov_model VARIABLE_IS_NOT_USED  *cov, double *v, double *Sign) {
  *v = - *x * *x;
  *Sign = 1.0;
}
void DGauss(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x; 
  *v = -2.0 * y * exp(- y * y);
}
void DDGauss(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x * *x; 
  *v = (4.0 * y - 2.0)* exp(- y);
}
void D3Gauss(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x * *x; 
  *v = *x * (12 - 8 * y) * exp(- y);
}
void D4Gauss(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x * *x; 
  *v = ((16 * y - 48) * y + 12) * exp(- y);
}
void InverseGauss(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = sqrt(-log(*x));
}
void nonstatLogInvGauss(double *x, cov_model VARIABLE_IS_NOT_USED *cov, 
		     double *left, double *right) {
  double 
    z = *x <= 0.0 ? sqrt(- *x) : 0.0;  
  int d,
    dim = cov->tsdim;
  for (d=0; d<dim; d++) {
    left[d] = -z;
    right[d] = z;
  }
}
double densityGauss(double *x, cov_model *cov) {  
    int d, dim=cov->tsdim;
    double x2=0.0;
    for (d=0; d<dim; d++) x2 += x[d] * x[d];
    return exp(- 0.25 * x2 - (double) dim * (M_LN2 + M_LN_SQRT_PI));
}
int struct_Gauss(cov_model *cov, cov_model **newmodel) {  
  ASSERT_NEWMODEL_NOT_NULL;

  switch (cov->role) {
  case ROLE_POISSON_GAUSS :// optimierte density fuer den Gauss-Fall
    double invscale;
    addModel(newmodel, GAUSS, cov);       
    addModel(newmodel, DOLLAR);
    kdefault(*newmodel, DSCALE, INVSQRTTWO);
    addModel(newmodel, TRUNCSUPPORT);
    InverseGauss(&GLOBAL.mpp.about_zero, cov, &invscale);
    kdefault(*newmodel, TRUNC_RADIUS, invscale);
    break;  
  case ROLE_MAXSTABLE :   // optimierte density fuer den Gauss-Fall
    // crash();
    addModel(newmodel, GAUSS_DISTR, cov); // to
    kdefault(*newmodel, GAUSS_DISTR_MEAN, 0.0);
    kdefault(*newmodel, GAUSS_DISTR_SD, INVSQRTTWO);
    return NOERROR;
  default : ILLEGAL_ROLE_STRUCT;
  }

  return NOERROR;
}

double IntUdeU2_intern(int d, double R, double expMR2) {
  // int_0^R u^{d-1} exp(-u^2) \D u
  return d == 0 ? (pnorm(R, 0.0, INVSQRTTWO, true, false) - 0.5)  * SQRTPI
    : d == 1 ? 0.5  * (1.0 - expMR2)
    : 0.5 * (expMR2 + (d - 1.0) * IntUdeU2_intern(d - 2, R, expMR2));    
}

double IntUdeU2(int d, double R) {
  // int_0^R u^{d-1} exp(-u^2) \D u
  return IntUdeU2_intern(d, R, exp(-R*R));
}

int initGauss(cov_model *cov, gen_storage *s) {

  if (hasNoRole(cov)) return NOERROR;

  if (HAS_SPECTRAL_ROLE(cov)) {

   spec_properties *cs = &(s->spec);
  
    if (cov->tsdim <= 2) return NOERROR;
    cs->density = densityGauss;
    return search_metropolis(cov, s);
  }

  else if (hasAnyShapeRole(cov)) {
    // int_{b(0,R) e^{-a r^2} dr = d b_d int_0^R r^{d-1} e^{-a r^2} dr
    // where a = 2.0 * xi / sigma^2
    // here : 2.0 is a factor so that the mpp function leads to the
    //            gaussian covariance model exp(-x^2)
    //        xi : 1 : integral ()
    //             2 : integral ()^2

   
   double R = RF_INF;
    int 
      dim = cov->tsdim;
    
    assert(cov->nr == GAUSS);
         
  
    assert(cov->mpp.maxheights[0] = 1.0);
    if (cov->mpp.moments >= 1) {
      cov->mpp.mM[1] = cov->mpp.mMplus[1] = 
	SurfaceSphere(dim - 1, 1.0) * IntUdeU2(dim - 1, R);
      int i;
      for (i=2; i<=cov->mpp.moments; i++) {
	cov->mpp.mM[i] = cov->mpp.mM[1] * pow((double) i, -0.5 * dim);
      }
    }    
    cov->mpp.maxheights[0] = 1.0;
  }

  else ILLEGAL_ROLE;
 
  return NOERROR;

}

void do_Gauss(cov_model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 

}

void spectralGauss(cov_model *cov, gen_storage *S, double *e) {   
  spectral_storage *s = &(S->Sspectral);
  if (cov->tsdim <= 2) {
    E12(s, cov->tsdim, 2.0 * sqrt(-log(1.0 - UNIFORM_RANDOM)), e);

  } else {
    metropolis(cov, S, e);
  }
}
void DrawMixGauss(cov_model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *random) {
  *random = 1.0;
}
double LogMixDensGauss(double VARIABLE_IS_NOT_USED *x, double VARIABLE_IS_NOT_USED logV, cov_model VARIABLE_IS_NOT_USED *cov) {
  return 0.0;
}

/*
void densGauss(double *x, cov_model *cov, double *v) {
  int 
    factor[MAXMPPDIM+1] = {0, 1 / M_SQRT_PI, INVPI, INVPI / M_SQRT_PI, 
			   INVPI * INVPI},
    dim = cov->tsdim;
    *v = factor[dim] * exp(- *x * *x);
}
*/

/*
void getMassGauss(double *a, cov_model *cov, double *kappas, double *m) {
  int i, j, k, kold,
    dim = cov->tsdim;
  double val[MAXMPPDIM + 1],
    sqrt2pi = SQRT2 * SQRTPI,
    factor[MAXMPPDIM+1] = {1, 
			   sqrt(2) / M_SQRT_PI, 
			   2 * INVPI,
			   2 * sqrt(2) * INVPI / M_SQRT_PI, 
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
void simuGauss(cov_model *cov, int dim, double *v) {
  int i;
  double 
    dummy;
  *v = 0.0;
  if (dim <= 2) {
    *v = dim == 1 ? fabs(GAUSS_RANDOM(1.0)) : rexp(1.0); 
  } else {
    for (i=0; i<dim; i++) {
      dummy = GAUSS_RANDOM(1.0);
      *v += dummy * dummy;
    }
    *v = sqrt(*v);
  }
}
*/



/* generalised fractal Brownian motion */
void genBrownian(double *x, cov_model *cov, double *v) {
  double 
    alpha = P0(BROWN_ALPHA),
    beta =  P0(BROWN_GEN_BETA),
    delta = beta / alpha;
  *v = 1.0 - pow(pow(*x, alpha) + 1.0, delta);
  //this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}

void loggenBrownian(double *x, cov_model *cov, double *v, double *Sign) {
  genBrownian(x, cov, v);
  *v = log(-*v);
  *Sign = -1.0;
}
void InversegenBrownian(double *x, cov_model *cov, double *v) {
  double 
    alpha = P0(BROWN_ALPHA),
    beta =  P0(BROWN_GEN_BETA),
    delta = beta / alpha;
  *v = pow(pow(*x + 1.0, 1.0/delta) - 1.0, 1.0/alpha); 
}
static bool genfbmmsg=true;
int checkgenBrownian(cov_model *cov){
  if (genfbmmsg) warning("Note that the definition of 'genfbm' has changed. This warning appears only once per session.");
  genfbmmsg = false;

  cov->logspeed = RF_INF;
  return NOERROR;
}
void rangegenBrownian(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[BROWN_ALPHA] = 0.0;
  range->max[BROWN_ALPHA] = 2.0;
  range->pmin[BROWN_ALPHA] = UNIT_EPSILON;
  range->pmax[BROWN_ALPHA] = 2.0 - UNIT_EPSILON;
  range->openmin[BROWN_ALPHA] = true;
  range->openmax[BROWN_ALPHA] = false;

  range->min[BROWN_GEN_BETA] = 0.0;
  range->max[BROWN_GEN_BETA] = 2.0;
  range->pmin[BROWN_GEN_BETA] = UNIT_EPSILON;
  range->pmax[BROWN_GEN_BETA] = 2.0 - UNIT_EPSILON;
  range->openmin[BROWN_GEN_BETA] = true;
  range->openmax[BROWN_GEN_BETA] = false;
}


/* gencauchy */
#define GENC_ALPHA 0
#define GENC_BETA 1
void generalisedCauchy(double *x, cov_model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA);
  *v = pow(1.0 + pow(*x, alpha), -beta/alpha);
}
void DgeneralisedCauchy(double *x, cov_model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA), ha, y=*x;
  if (y ==0.0) {
    *v = ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? -INFTY : -beta); 
  } else {
    ha=pow(y, alpha - 1.0);
    *v = -beta * ha * pow(1.0 + ha * y, -beta / alpha - 1.0);
  }
}
void DDgeneralisedCauchy(double *x, cov_model *cov, double *v){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA), ha, y=*x;
  if (y ==0.0) {
    *v = ((alpha==2.0) ? beta * (beta + 1.0) : INFTY); 
  } else {
    ha=pow(y, alpha);
    *v = beta * ha / (y * y) * (1.0 - alpha + (1.0 + beta) * ha)
      * pow(1.0 + ha, -beta / alpha - 2.0);
  }
}
void loggeneralisedCauchy(double *x, cov_model *cov, double *v, double *Sign){
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA);
  *v = log(1.0 + pow(*x, alpha)) *  -beta/alpha;
  *Sign = 1.0;
}
void InversegeneralisedCauchy(double *x, cov_model *cov, double *v) {
  double alpha = P0(GENC_ALPHA), beta=P0(GENC_BETA);
  *v =  (*x == 0.0) ? RF_INF : pow(pow(*x, -alpha / beta) - 1.0, 1.0 / alpha);
  // MLE works much better with 0.01 then with 0.05
}
int checkgeneralisedCauchy(cov_model *cov){
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  cov->monotone = P0(GENC_ALPHA) <= 1.0 ? COMPLETELY_MON : NORMAL_MIXTURE;
  return NOERROR;
}

void rangegeneralisedCauchy(cov_model *cov, range_type *range) {
  bool tcf = isTcf(cov->typus) || cov->isoown == SPHERICAL_ISOTROPIC;
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

void coinitgenCauchy(cov_model *cov, localinfotype *li) {
  double thres[2] = {0.5, 1.0}, alpha=P0(GENC_ALPHA); 
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
  if (P0(GENC_ALPHA) <= 1.0) {
    li->value[0] = 1.0;
    li->msg[0] = MSGLOCAL_OK;
  } else {
    li->value[0] = 1.5;
    li->msg[0] = MSGLOCAL_JUSTTRY;
  }
}



/* epsC -> generalised Cauchy : leading 1 is now an eps */
#define EPS_ALPHA 0
#define EPS_BETA 1
#define EPS_EPS 2
void epsC(double *x, cov_model *cov, double *v){
  double 
    alpha = P0(EPS_ALPHA), 
    beta=P0(EPS_BETA),
    eps=P0(EPS_EPS);
  *v = pow(eps + pow(*x, alpha), -beta/alpha);
 }
void logepsC(double *x, cov_model *cov, double *v, double *Sign){
  double 
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  *v = log(eps + pow(*x, alpha)) * -beta/alpha;
  *Sign = 1.0;
 }
void DepsC(double *x, cov_model *cov, double *v){
  double ha, 
    y=*x,
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  if (y ==0.0) {
    *v = (eps == 0.0) ? -INFTY : 
      ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? -INFTY : -beta); 
  } else {
    ha=pow(y, alpha - 1.0);
    *v = -beta * ha * pow(eps + ha * y, -beta / alpha - 1.0);
  }
}
void DDepsC(double *x, cov_model *cov, double *v){
  double ha, 
    y=*x,
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  if (y ==0.0) {
    *v = (eps == 0.0) ? INFTY : ((alpha==2.0) ? beta * (beta + 1.0) : INFTY); 
  } else {
    ha=pow(y, alpha);
    *v = beta * ha / (y * y) * ( (1.0 - alpha) * eps + (1.0 + beta) * ha)
      * pow(eps + ha, -beta / alpha - 2.0);
  }
}
void InverseepsC(double *x, cov_model *cov, double *v){
  double 
    alpha = P0(EPS_ALPHA),
    beta=P0(EPS_BETA), 
    eps=P0(EPS_EPS);
  *v = (*x == 0.0) ? RF_INF : pow(pow(*x, -alpha / beta) - eps, 1.0 / alpha);
}
int checkepsC(cov_model *cov){
  double eps=P0(EPS_ALPHA);
  int i, err;
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if ((err = checkkappas(cov, false)) != NOERROR) return err;
  kdefault(cov, EPS_ALPHA, 1.0); 
  kdefault(cov, EPS_BETA, 1.0); 
  kdefault(cov, EPS_EPS, 0.0); 
  if (ISNAN(eps) || eps == 0.0) {
    //  cov->domown=GENERALISEDCOVARIANCE; // later
    for (i=CircEmbed; i<Nothing; i++) cov->pref[i] = PREF_NONE;
  }
  
  return NOERROR;
}

void rangeepsC(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
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


/* gengneiting */
#define GENGNEITING_K 0
#define GENGNEITING_MU 1
void genGneiting(double *x, cov_model *cov, double *v)
{
  int kk = P0INT(GENGNEITING_K);
  double s,
    mu=P0(GENGNEITING_MU),
    y=*x;
  if (y >= 1.0) {
    *v = 0.0; 
    return;
  }
  s = mu + 2.0 * (double) kk + 0.5;


  switch (kk) {
  case 0:
    *v = 1.0;
    break;
  case 1:
    *v =  1.0 + s*y ;
    break;
  case 2:
    *v = 1.0 + y * (s + y * (s*s-1.0) * ONETHIRD);
    break;
  case 3:
    *v = 1 + y * (s + y * (0.2 * (2.0*s*s - 3 + y * (s*s-4) * s * ONETHIRD)));
    break;
  default : BUG;
  }
  *v *=  pow(1.0 - y, s);
}


// control thanks to http://calc101.com/webMathematica/derivatives.jsp#topdoit

void DgenGneiting(double *x, cov_model *cov, double *v)
{
  int kk = P0INT(GENGNEITING_K);
  double s,
    mu=P0(GENGNEITING_MU), 
    y=*x;
  if (y >= 1.0) {
    *v = 0.0; 
    return;
  }
  s = mu + 2.0 * (double) kk + 0.5;
  
  switch (kk) {
  case 0: 
    *v = s;
    break;
  case 1:
    *v =  y * s  * (s + 1.0);
    break;
  case 2: 
    *v = ONETHIRD * (s * s + 3.0 * s + 2.0) * y * ((s - 1.0) * y + 1.0);
    break;
  case 3: 
    *v = y * (s * (s + 5) + 6) * (y * (s * (s-2) * y + 3 * s - 3) + 3) / 15.0;
    break;
  default : BUG;
  }
  *v *=  -pow(1.0 - y, s - 1.0);
  
}
void DDgenGneiting(double *x, cov_model *cov, double *v){
  int kk = P0INT(GENGNEITING_K);
  double s,
    mu=P0(GENGNEITING_MU), 
    y=*x;
  if (y >= 1.0) {
    *v = 0.0; 
    return;
  }
  s = mu + 2.0 * (double) kk + 0.5;
  switch (kk) {
  case 0: 
    *v = s * (s-1);
    break;
  case 1:
    *v = s  * (s + 1.0) * (s * y - 1) ; 
    break;
  case 2: {
    double s2;
    s2 = s * s;
    *v = ONETHIRD * (s2 + 3.0 * s + 2.0) * ( y * ((s2 - 1) * y - s + 2) - 1);
  }
    break;
  case 3: 
    double s2;
    s2 = s * s;
    *v = (s2  + 5 * s + 6) / 15.0 * 
      (y * (y * ((s2 - 4) * s * y + 6 * s - 3) -3 * s + 6) - 3);
    break;
  default : BUG;
  }
  *v *=  pow(1.0 - y, s - 2.0);
}


int checkgenGneiting(cov_model *cov){
  double mu=P0(GENGNEITING_MU), 
    dim = 2.0 * mu;
  cov->maxdim = (ISNAN(dim) || dim >= INFDIM) ? INFDIM : (int) dim;
  return NOERROR;
}

void rangegenGneiting(cov_model *cov, range_type *range){
  range->min[GENGNEITING_K] = range->pmin[GENGNEITING_K] = 0;
  range->max[GENGNEITING_K] = range->pmax[GENGNEITING_K] = 3; 
  range->openmin[GENGNEITING_K] = false;
  range->openmax[GENGNEITING_K] = false;
  
  range->min[GENGNEITING_MU] = 0.5 * (double) cov->tsdim; 
  range->max[GENGNEITING_MU] =  RF_INF;
  range->pmin[GENGNEITING_MU] = range->min[GENGNEITING_MU];
  range->pmax[GENGNEITING_MU] = range->pmin[GENGNEITING_MU] + 10.0;
  range->openmin[GENGNEITING_MU] = false;
  range->openmax[GENGNEITING_MU] = true;
}

/*
double porcusG(double t, double nu, double mu, double gamma) {
  double t0 = fabs(t);
  if (t0 > mu) return 0.0;
  return pow(t0, nu) * pow(1.0 - t0 / mu, gamma);
}
*/



#define GNEITING_K GENGNEITING_K    // important to keep !
#define GNEITING_MU 1
#define GNEITING_S 2
#define GNEITING_SRED 3
#define GNEITING_GAMMA 4
#define GNEITING_CDIAG 5
#define GNEITING_RHORED 6
#define GNEITING_C 7
double biGneitQuot(double t, double* scale, double *gamma) {
  double t0 = fabs(t);
  return pow(1.0 - t0 / scale[0], gamma[0]) *
    pow(1.0 - t0 / scale[3], gamma[3]) * pow(1.0 - t0 / scale[1], -2 * gamma[1]);
}

void biGneitingbasic(cov_model *cov, 
		     double *scale, 
		     double *gamma, 
		     double *cc
		     ){ 
  double 
    Sign, x12, min, det, a, b, c, sum,
    k = (double) (P0INT(GNEITING_K)),
    kP1 = k + (k >= 1),
    Mu = P0(GNEITING_MU),
    nu = Mu + 0.5,
    *sdiag = P(GNEITING_S), // >= 0 
    s12red = P0(GNEITING_SRED), // in [0,1]
    s12 = s12red * (sdiag[0] <= sdiag[1] ? sdiag[0] : sdiag[1]),

    *tildegamma = P(GNEITING_GAMMA), // 11,22,12

    *Cdiag = P(GNEITING_CDIAG),
    *C = P(GNEITING_C),
    rho = P0(GNEITING_RHORED);

  scale[0] = sdiag[0];
  scale[1] = scale[2] = s12; // = scale[2]
  scale[3] = sdiag[1];

  gamma[0] = tildegamma[0];
  gamma[1] =  gamma[2] = tildegamma[1]; //gamma[2]
  gamma[3] = tildegamma[2];

  sum = 0.0;
  if (scale[0] == scale[1]) sum += gamma[0];
  if (scale[0] == scale[3]) sum += gamma[3];
  if (sum > 2.0 * gamma[1]) error("values of gamma not valid.");

  min = 1.0;
  a = 2 * gamma[1] - gamma[0] - gamma[3];
  b = - 2 * gamma[1] * (scale[0] + scale[3]) + gamma[0] * (scale[1] + scale[3])
    + gamma[3] * (scale[1] + scale[0]);
  c = 2 * gamma[1] * scale[0] * scale[3] - gamma[0] * scale[1] * scale[3]
    - gamma[3] * scale[1] * scale[0]; 
  det = b * b - 4 * a * c;
  
  if (det >= 0) {
    det = sqrt(det);
    for (Sign=-1.0; Sign<=1.0; Sign+=2.0) {
      x12 = 0.5 / a * (-b + Sign * det);
      if (x12>0 && x12<scale[1]) {
	double dummy = biGneitQuot(x12, scale, gamma);
	if (dummy < min) min = dummy;
      }
    }    
  }

  cc[0] = C[0] = Cdiag[0];
  cc[3] = C[2] = Cdiag[1];
  cc[1] = cc[2] = C[1] = rho * sqrt(C[0] * C[2] * min) *
    pow(scale[1] * scale[1] / (scale[0] * scale[3]), 0.5 * (nu + 1 + 2.0 * k)) *
    exp(lgammafn(1.0 + gamma[1]) - lgammafn(2.0 + nu + gamma[1] + kP1)
	+ 0.5 * (  lgammafn(2 + nu + gamma[0] + kP1) - lgammafn(1 + gamma[0])
		 + lgammafn(2 + nu + gamma[3] + kP1) - lgammafn(1 + gamma[3]))
	);
}


int initbiGneiting(cov_model *cov, gen_storage *stor) {
  double  
    *rho = P(GNEITING_RHORED),
    *s = P(GNEITING_S),
    *sred = P(GNEITING_SRED),
    *p_gamma = P(GNEITING_GAMMA),
    *c = P(GNEITING_C),
    *cdiag = P(GNEITING_CDIAG);
  bool check = stor->check;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  
  if (s == NULL) {
    PALLOC(GNEITING_S, 2, 1);
    s = P(GNEITING_S);
    s[0] = s[1] = 1.0;
  }
  
  if (sred == NULL) {
    PALLOC(GNEITING_SRED, 1, 1);
    sred = P(GNEITING_SRED);
    sred[0] = 1.0;
  }
  
  if (sred[0] == 1.0) {
    double sum =0.0;
    if (s[0] <= s[1]) sum += p_gamma[0];
    if (s[1] <= s[0]) sum += p_gamma[2];
    if (sum > 2.0 * p_gamma[1]) {
      SERR7("if '%s'=1, then %s[2] must be greater than min(%s[1], %s[3]) / 2 if%s[1] and %s[3] differ and %s[1] otherwise.", KNAME(GNEITING_SRED),
	    KNAME(GNEITING_GAMMA), KNAME(GNEITING_GAMMA),
	    KNAME(GNEITING_GAMMA), KNAME(GNEITING_GAMMA), 
	    KNAME(GNEITING_GAMMA), KNAME(GNEITING_GAMMA));
    }
  }

  if  (S->cdiag_given) {
    if (cdiag == NULL) {
      PALLOC(GNEITING_CDIAG, 2, 1);
      cdiag = P(GNEITING_CDIAG);
      cdiag[0] = cdiag[1] = 1.0;
    }
    if (rho == NULL) 
      QERRC2(GNEITING_RHORED, 
	     "'%s' and '%s' must be set at the same time ", GNEITING_CDIAG,
	     GNEITING_RHORED);
    if (check && c != NULL) {
      if (cov->nrow[GNEITING_C] != 3 || cov->ncol[GNEITING_C] != 1)
	QERRC(GNEITING_C, "must be a 3 x 1 vector");
      if (fabs(c[i11] - cdiag[0]) > c[i11] * epsilon || 
	  fabs(c[i22] - cdiag[1]) > c[i22] * epsilon ) {
	QERRC1(GNEITING_C, "does not match '%s'.", GNEITING_CDIAG);
      }
      double savec12 = c[i21];
      biGneitingbasic(cov, S->scale, S->gamma, S->c);
      if (fabs(c[i21] - savec12) > fabs(c[i21]) * epsilon)
 	QERRC1(GNEITING_C, "does not match '%s'.", GNEITING_RHORED);
    } else {
      if (PisNULL(GNEITING_C)) PALLOC(GNEITING_C, 3, 1);
      c = P(GNEITING_C);
      biGneitingbasic(cov, S->scale, S->gamma, S->c);
    }
  } else {
    if (c == NULL) 
      QERRC2(GNEITING_C, "either '%s' or '%s' must be set", 
	    GNEITING_C, GNEITING_RHORED);
    if (!ISNAN(c[i11]) && !ISNAN(c[i22]) && (c[i11] < 0.0 || c[i22] < 0.0))
      QERRC2(GNEITING_C, "variance parameter %s[0], %s[2] must be non-negative",
	     GNEITING_C, GNEITING_C);
    if (PisNULL(GNEITING_CDIAG)) PALLOC(GNEITING_CDIAG, 2, 1);
    if (PisNULL(GNEITING_RHORED)) PALLOC(GNEITING_RHORED, 1, 1);
    cdiag = P(GNEITING_CDIAG);
    rho = P(GNEITING_RHORED);
    cdiag[0] = c[i11];
    cdiag[1] = c[i22];
    double savec1 = c[i21];
    if (savec1 == 0.0)  rho[0] = 0.0; // wichtig falls c[0] oder c[2]=NA
    else {
      rho[0] = 1.0;
      biGneitingbasic(cov, S->scale, S->gamma, S->c);
      rho[0] = savec1 / c[i21];
    }
  }

  biGneitingbasic(cov, S->scale, S->gamma, S->c);
  cov->initialised = true;
  return NOERROR;
}


void kappa_biGneiting(int i, cov_model *cov, int *nr, int *nc){
  *nc = *nr = i < CovList[cov->nr].kappas? 1 : -1;
  if (i==GNEITING_S || i==GNEITING_CDIAG) *nr=2; else
    if (i==GNEITING_GAMMA || i==GNEITING_C) *nr=3 ;
}

void biGneiting(double *x, cov_model *cov, double *v) { 
  double z, 
    mu = P0(GNEITING_MU);
  int i;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  // wegen ML aufruf immer neu berechnet
 
  assert(cov->initialised);
   
  for (i=0; i<4; i++) {
    if (i==2) {
      v[2] = v[1];
      continue;
    }
    z = fabs(*x / S->scale[i]);
    P(GENGNEITING_MU)[0] = mu + S->gamma[i] + 1.0;
    genGneiting(&z, cov, v + i);
    v[i] *= S->c[i]; 
  }
  P(GENGNEITING_MU)[0] = mu;
}


void DbiGneiting(double *x, cov_model *cov, double *v){ 
  double z, 
    mu = P0(GENGNEITING_MU);
  int i;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  assert(cov->initialised);
 

  for (i=0; i<4; i++) {
    if (i==2) {
      v[2] = v[1];
      continue;
    }
    z = fabs(*x / S->scale[i]);
    P(GENGNEITING_MU)[0] = mu + S->gamma[i] + 1.0;
    DgenGneiting(&z, cov, v + i);
    v[i] *= S->c[i] / S->scale[i];
  }
  P(GENGNEITING_MU)[0] = mu;
}


void DDbiGneiting(double *x, cov_model *cov, double *v){ 
  double z,
    mu = P0(GENGNEITING_MU);
  int i;
 biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  assert(cov->initialised);

 
   for (i=0; i<4; i++) {
    if (i==2) {
      v[2] = v[1];
      continue;
    }
    z = fabs(*x / S->scale[i]);
    P(GENGNEITING_MU)[0] = mu + S->gamma[i] + 1.0;
    DDgenGneiting(&z, cov, v + i);
    v[i] *= S->c[i] / (S->scale[i] * S->scale[i]);
  }
  P(GENGNEITING_MU)[0] = mu;
}

int checkbiGneiting(cov_model *cov) { 
  int err;
  gen_storage s;
  gen_NULL(&s);
  s.check = true;
  
  if ((err = checkkappas(cov, false)) != NOERROR) return err;
  if (PisNULL(GNEITING_K)) QERRC(GNEITING_K, "must be given.");
  if (PisNULL(GNEITING_MU)) QERRC(GNEITING_MU, "must be given.");
  if (PisNULL(GNEITING_GAMMA)) QERRC(GNEITING_GAMMA,"must be given.");

  NEW_STORAGE(biwm);
  biwm_storage *S = cov->Sbiwm;
  S->cdiag_given = !PisNULL(GNEITING_CDIAG) || !PisNULL(GNEITING_RHORED);
 
  if ((err=initbiGneiting(cov, &s)) != NOERROR) return err;

  int dim = 2.0 * P0(GENGNEITING_MU);
  cov->maxdim = (ISNAN(dim) || dim >= INFDIM) ? INFDIM : (int) dim;
  return NOERROR;
}
  
sortsofparam paramtype_biGneiting(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  return k == GNEITING_S ? SCALEPARAM : 
    k == GNEITING_CDIAG ? VARPARAM : 
    k == GNEITING_C ? DONOTRETURNPARAM :
    (k==GNEITING_MU || k==GNEITING_GAMMA) ? CRITICALPARAM : ANYPARAM;
}


void rangebiGneiting(cov_model *cov, range_type *range){
 // *n = P0(GNEITING_K], 
  range->min[GNEITING_K] = range->pmin[GNEITING_K] = 0;
  range->max[GNEITING_K] = range->pmax[GNEITING_K] = 3;
  range->openmin[GNEITING_K] = false;
  range->openmax[GNEITING_K] = false;
  
 // *mu = P0(GNEITING_MU], 
  range->min[GNEITING_MU] = 0.5 * (double) cov->tsdim; 
  range->max[GNEITING_MU] =  RF_INF;
  range->pmin[GNEITING_MU] = range->min[GNEITING_MU];
  range->pmax[GNEITING_MU] = range->pmin[GNEITING_MU] + 10.0;
  range->openmin[GNEITING_MU] = false;
  range->openmax[GNEITING_MU] = true;
  
 // *scalediag = P0(GNEITING_S], 
  range->min[GNEITING_S] = 0.0;
  range->max[GNEITING_S] = RF_INF;
  range->pmin[GNEITING_S] = 1e-2;
  range->pmax[GNEITING_S] = 6.0;
  range->openmin[GNEITING_S] = true;
  range->openmax[GNEITING_S] = true;
  
  // *scalered12 = P0(GNEITING_SRED], 
  range->min[GNEITING_SRED] = 0;
  range->max[GNEITING_SRED] = 1;
  range->pmin[GNEITING_SRED] = 0.01;
  range->pmax[GNEITING_SRED] = 1;
  range->openmin[GNEITING_SRED] = true;
  range->openmax[GNEITING_SRED] = false;
    
  //   *gamma = P0(GNEITING_GAMMA], 
  range->min[GNEITING_GAMMA] = 0.0;
  range->max[GNEITING_GAMMA] = RF_INF;
  range->pmin[GNEITING_GAMMA] = 1e-5;
  range->pmax[GNEITING_GAMMA] = 100.0;
  range->openmin[GNEITING_GAMMA] = false;
  range->openmax[GNEITING_GAMMA] = true;
   
  //    *c diag = P0(GNEITING_CDIAG]; 
  range->min[GNEITING_CDIAG] = 0;
  range->max[GNEITING_CDIAG] = RF_INF;
  range->pmin[GNEITING_CDIAG] = 1e-05;
  range->pmax[GNEITING_CDIAG] = 1000;
  range->openmin[GNEITING_CDIAG] = true;
  range->openmax[GNEITING_CDIAG] = true; 
  
 //    *rho = P0(GNEITING_RHORED]; 
  range->min[GNEITING_RHORED] = -1;
  range->max[GNEITING_RHORED] = 1;
  range->pmin[GNEITING_RHORED] = -0.95;
  range->pmax[GNEITING_RHORED] = 0.95;
  range->openmin[GNEITING_RHORED] = false;
  range->openmax[GNEITING_RHORED] = false;    

 //    *rho = P0(GNEITING_C]; 
   range->min[GNEITING_C] = RF_NEGINF;
  range->max[GNEITING_C] = RF_INF;
  range->pmin[GNEITING_C] = -1000;
  range->pmax[GNEITING_C] = 1000;
  range->openmin[GNEITING_C] = true;
  range->openmax[GNEITING_C] = true; 
  
 }


/* Gneiting's functions -- alternative to Gaussian */
// #define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
#define GNEITING_ORIG 0
#define gneiting_origmu 1.5
#define NumericalScale 0.301187465825
#define kk_gneiting 3
#define mu_gneiting 2.683509
#define s_gneiting 0.2745640815
void Gneiting(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){ 
  double y = *x * cov->q[0];
  genGneiting(&y, cov, v);
}

 
//void InverseGneiting(cov_model *cov, int scaling) {return 0.5854160193;}
void DGneiting(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){ 
  double y = *x * cov->q[0];
  DgenGneiting(&y, cov, v);
  *v  *= cov->q[0];
}


//void InverseGneiting(cov_model *cov, int scaling) {return 0.5854160193;}
void DDGneiting(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){ 
  double y = *x * cov->q[0];
  DgenGneiting(&y, cov, v);
  *v  *= cov->q[0] * cov->q[0];
}

int checkGneiting(cov_model *cov) { 
  int err;
  kdefault(cov, GNEITING_ORIG, 1);
  if ((err = checkkappas(cov)) != NOERROR) return err;
  bool orig = P0INT(GNEITING_ORIG);
  PFREE(GNEITING_ORIG);

  assert(cov->q == NULL);
  cov->nr = GNEITING_INTERN;
  QALLOC(1);
  
  if (orig) {
    cov->q[0] = NumericalScale;
    kdefault(cov, GENGNEITING_MU, gneiting_origmu);
    kdefault(cov, GENGNEITING_K, kk_gneiting);
  } else {
    cov->q[0] = s_gneiting;
    kdefault(cov, GENGNEITING_MU, mu_gneiting);
    kdefault(cov, GENGNEITING_K, kk_gneiting);
  }
  return checkgenGneiting(cov);
}

void rangeGneiting(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[GNEITING_ORIG] = range->pmin[GNEITING_ORIG] = 0;
  range->max[GNEITING_ORIG] = range->pmax[GNEITING_ORIG] = 1;
  range->openmin[GNEITING_ORIG] = false;
  range->openmax[GNEITING_ORIG] = false;
}



/* hyperbolic */
#define BOLIC_NU 0
#define BOLIC_XI 1
#define BOLIC_DELTA 2
void hyperbolic(double *x, cov_model *cov, double *v){ 
  double Sign;
  loghyperbolic(x, cov, v, &Sign);
  *v = exp(*v);
}
void loghyperbolic(double *x, cov_model *cov, double *v, double *Sign){ 
  double 
    nu = P0(BOLIC_NU),
    xi=P0(BOLIC_XI), 
    delta=P0(BOLIC_DELTA);
  static double nuOld = RF_INF;
  static double xiOld= RF_INF;
  static double deltaOld = RF_INF;
  static double deltasq;
  static double xidelta;
  static double logconst;
  double 
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
    *v = logWM(y * xi, nu, 0.0);
  } else if (xi==0) { //cauchy   => NU2 < 0 !
    y /= delta;
    /* note change in Sign as NU2<0 */
    *v = log(1 + y * y) * 0.5 * nu; 
  } else {
    if ((nu!=nuOld) || (xi!=xiOld) || (delta!=deltaOld)) {
    nuOld = nu; 
    xiOld = xi;
    deltaOld = delta;
    deltasq = deltaOld * deltaOld;
    xidelta = xiOld * deltaOld;
    logconst = xidelta - log(bessel_k(xidelta, nuOld, 2.0)) 
      - nu * log(deltaOld);
    }
    y=sqrt(deltasq + y * y);  
    double xiy = xi * y;
    *v = logconst + nu * log(y) + log(bessel_k(xiy, nu, 2.0)) - xiy;
  }
}
void Dhyperbolic(double *x, cov_model *cov, double *v)
{ 
  double 
    nu = P0(BOLIC_NU), 
    xi=P0(BOLIC_XI), 
    delta=P0(BOLIC_DELTA);
  static double nuOld = RF_INF;
  static double xiOld= RF_INF;
  static double deltaOld = RF_INF;
  static double deltasq;
  static double xidelta;
  static double logconst;
  double s,xi_s, logs, 
    y = *x;
  if (y==0.0) { 
    *v = 1.0;
    return;
  }
  if (delta==0) { // whittle matern
    *v = xi * DWM(y * xi, nu, 0.0);
    *v *= xi;
  } else if (xi==0) { //cauchy
    double ha;
    y /= delta;
    ha = y * y;
    /* note change in Sign as NU2<0 */
    *v = nu * fabs(y) * pow(1.0 + ha, 0.5 * nu - 1.0) / delta;
  } else {
    if ((nu!=nu) || (xi!=xi) || (delta!=delta)) {
      nuOld = nu; 
      xiOld= xi;
      deltaOld = delta;
      deltasq = deltaOld * deltaOld;
      xidelta = xiOld * deltaOld;
      logconst = xidelta - log(bessel_k(xidelta, nuOld, 2.0)) 
	- nu * log(deltaOld);
    }
    s=sqrt(deltasq + y * y);
    xi_s = xi * s;
    logs = log(s);  
    *v = - y * xi * exp(logconst + (nu-1.0) * logs 
		   +log(bessel_k(xi_s, nu-1.0, 2.0)) - xi_s);
  }
}
int checkhyperbolic(cov_model *cov){
  double 
    nu = P0(BOLIC_NU),
    xi=P0(BOLIC_XI),
    delta=P0(BOLIC_DELTA);
  char msg[255];
  int i;
  for (i=0; i<= Nothing; i++)
    cov->pref[i] *= (ISNAN(nu) || nu < BesselUpperB[i]);
  if (nu>0) {
    if ((delta<0) || (xi<=0)) {
      sprintf(msg, "xi>0 and delta>=0 if nu>0. Got nu=%f and delta=%f for xi=%f", nu, delta, xi);
      SERR(msg);
    }
  } else if (nu<0) {
    if ((delta<=0) || (xi<0)) {
      sprintf(msg, "xi>=0 and delta>0 if nu<0. Got nu=%f and delta=%f for xi=%f", nu, delta, xi);
      SERR(msg);
    }
  } else { // nu==0.0
    if ((delta<=0) || (xi<=0)) {
      sprintf(msg, "xi>0 and delta>0 if nu=0. Got nu=%f and delta=%f for xi=%f", nu, delta, xi);
      SERR(msg);
    }
  }
  return NOERROR;
}
void rangehyperbolic(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
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


/* iaco cesare model */
#define CES_NU 0
#define CES_LAMBDA 1
#define CES_DELTA 2
void IacoCesare(double *x, cov_model *cov, double *v){
    double
      nu = P0(CES_NU), 
      lambda=P0(CES_LAMBDA), 
      delta=P0(CES_DELTA);
    *v = pow(1.0 + pow(x[0], nu) + pow(x[1], lambda), - delta); 
}
void rangeIacoCesare(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[CES_NU] = 0.0;
  range->max[CES_NU] = 2.0;
  range->pmin[CES_NU] = 0.0;
  range->pmax[CES_NU] = 2.0;
  range->openmin[CES_NU] = true;
  range->openmax[CES_NU] = false;

  range->min[CES_LAMBDA] = 0.0;
  range->max[CES_LAMBDA] = 2.0;
  range->pmin[CES_LAMBDA] = 0.0;
  range->pmax[CES_LAMBDA] = 2.0;
  range->openmin[CES_LAMBDA] = true;
  range->openmax[CES_LAMBDA] = false;
 
  range->min[CES_DELTA] = 0.0;
  range->max[CES_DELTA] = RF_INF;
  range->pmin[CES_DELTA] = 0.0;
  range->pmax[CES_DELTA] = 10.0;
  range->openmin[CES_DELTA] = true;
  range->openmax[CES_DELTA] = true;
}



/* Kolmogorov model */
void Kolmogorov(double *x, cov_model *cov, double *v){
#define fourthird 1.33333333333333333333333333333333333
#define onethird 0.333333333333333333333333333333333
  int d, i, j,
    dim = cov->tsdim,
    dimP1 = dim + 1,
    dimsq = dim * dim;
  double
    rM23, r23,
    r2 =0.0;

  assert(dim == cov->vdim[0] && dim == cov->vdim[1]);

  for (d=0; d<dimsq; v[d++] = 0.0);
  for (d=0; d<dim; d++) r2 += x[d] * x[d];
  if (r2 == 0.0) return;

  rM23 = onethird / r2; // r^-2 /3

  for (d=0; d<dimsq; d+=dimP1) v[d] = fourthird;
  for (d=i= 0; i<dim ; i++) {
    for(j=0; j<dim; j++) {
      v[d++] -= rM23 * x[i] * x[j];
    }
  }

  r23 = -pow(r2, onethird);  // ! -
  for (d=0; d<dimsq; v[d++] *= r23);

}


int checkKolmogorov(cov_model *cov) { 
   if (cov->xdimown != 3) SERR1("dim (%d) != 3", cov->xdimown);
 if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;
  return NOERROR;
}
 

/* local-global distinguisher */
#define LGD_ALPHA 0
#define LGD_BETA 1
void lgd1(double *x, cov_model *cov, double *v) {
  double y = *x, alpha = P0(LGD_ALPHA), beta=P0(LGD_BETA);
  *v = (y == 0.0) ? 1.0 : (y < 1.0) 
    ? 1.0 - beta / (alpha + beta) * pow(y, alpha)
    : alpha / (alpha + beta) * pow(y,  -beta);
}
void Inverselgd1(double *x, cov_model *cov, double *v) {
  double alpha = P0(LGD_ALPHA), beta=P0(LGD_BETA);
  ERR("scle of lgd1 not programmed yet"); 
   // 19 next line?!
  // 1.0 / .... fehlt auch
  *v = (19 * alpha < beta)
     ? exp(log(1 - *x * (alpha + beta) / beta) / alpha)
     : exp(log(*x * (alpha + beta) / alpha) / beta);
}
void Dlgd1(double *x, cov_model *cov, double *v){
  double y=*x, pp, alpha = P0(LGD_ALPHA), beta=P0(LGD_BETA);
  if (y == 0.0) {
    *v = 0.0;// falscher Wert, aber sonst NAN-Fehler#
    return;
  }
  pp = ( (y < 1.0) ? alpha : -beta ) - 1.0;
  *v = - alpha * beta / (alpha + beta) * exp(pp * y);
}
int checklgd1(cov_model *cov){
  double dim = 2.0 * (1.5 - P0(LGD_ALPHA));
  cov->maxdim = (ISNAN(dim) || dim >= 2.0) ? 2 : (int) dim;
  return NOERROR;
}
void rangelgd1(cov_model *cov, range_type *range) {
  range->min[LGD_ALPHA] = 0.0;
  range->max[LGD_ALPHA] = (cov->tsdim==2) ? 0.5 : 1.0;
  range->pmin[LGD_ALPHA] = 0.01;
  range->pmax[LGD_ALPHA] = range->max[LGD_ALPHA];
  range->openmin[LGD_ALPHA] = true;
  range->openmax[LGD_ALPHA] = false;

  range->min[LGD_BETA] = 0.0;
  range->max[LGD_BETA] = RF_INF;
  range->pmin[LGD_BETA] = 0.01;
  range->pmax[LGD_BETA] = 20.0;
  range->openmin[LGD_BETA] = true;
  range->openmax[LGD_BETA] = true;
 
}


/* mastein */
// see Hypermodel.cc



/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */ 

#define MATERN_NU_THRES 100
double WM(double x, double nu, double factor) {
  // check calling functions, like hyperbolic and gneiting if any changings !!
  return exp(logWM(x, nu, factor));
}

double logWM(double x, double nu, double factor) {
  // check calling functions, like hyperbolic and gneiting if any changings !!

  static double nuOld=RF_INF;
  static double loggamma;
  double v, y, Sign,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    scale = (factor != 0.0) ? factor * sqrt(nuThres) : 1.0;

  if (x > LOW_BESSELK) {
    if (nuThres != nuOld) {
      nuOld = nuThres;
      loggamma = lgammafn(nuThres);
    }
    y = x  * scale;
    v = LOG2 + nuThres * log(0.5 * y) - loggamma + 
		  log(bessel_k(y, nuThres, 2.0)) - y;
  } else v = 0.0;
    
  if (nu > MATERN_NU_THRES) { // factor!=0.0 && 
    double w, 
      g = MATERN_NU_THRES / nu;
    y = x * factor / 2;
    logGauss(&y, NULL, &w, &Sign);
    v = v * g + (1.0 - g) * w;
  }

  return v;
}


double DWM(double x, double nu, double factor) { 
  static double nuOld=RF_INF;
  static double loggamma;
  double   y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    scale = (factor != 0.0) ? factor * sqrt(nuThres) : 1.0;
  
  if (x > LOW_BESSELK) {
    if (nuThres!=nuOld) {
      nuOld = nuThres;
      loggamma = lgammafn(nuThres);
    }
    y = x * scale;  
    v = - 2.0 * exp(nuThres * log(0.5 * y) - loggamma + 
			     log(bessel_k(y, nuThres - 1.0, 2.0)) - y);
  } else {
    v = (nuThres > 0.5) ? 0.0 : (nuThres < 0.5) ? INFTY : 1.253314137;
  }
  v *= scale;

  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    y = x * scale;
    DGauss(&y, NULL, &w);
    w *= scale;
    v = v * g + (1.0 - g) * w;
  }
  return v;
}

double DDWM(double x, double nu, double factor) { 
  static double nuOld=RF_INF;
  static double gamma;
  double  y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    scale = (factor != 0.0) ? factor * sqrt(nuThres) : 1.0,
    scaleSq  = scale * scale;
   
  if (x > LOW_BESSELK) {
    if (nuThres!=nuOld) {
      nuOld = nuThres;
      gamma = gammafn(nuThres);
    }
    y = x * scale;
    v = pow(0.5 * y , nuThres - 1.0) / gamma *
      (- bessel_k(y, nuThres - 1.0, 1.0) + y * bessel_k(y, nuThres - 2.0, 1.0));
  } else {
    v = (nu > 1.0) ? -0.5 / (nu - 1.0) : INFTY;
  }
  v *= scaleSq;

  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    scaleSq = scale * scale;
    y = x * scale;
    DDGauss(&y, NULL, &w);
    w *= scaleSq;
    v = v * g + (1.0 - g) * w;
  }
  return v;
}

double D3WM(double x, double nu, double factor) { 
  static double nuOld=RF_INF;
  static double gamma;
  double y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    scale = (factor != 0.0) ? factor * sqrt(nuThres) : 1.0,
    scaleSq  = scale * scale;
  
  if (x > LOW_BESSELK) {
    if (nuThres!=nuOld) {
      nuOld = nuThres;
      gamma = gammafn(nuThres);
    }
    y = x * scale;
    v = pow(0.5 * y , nuThres - 1.0) / gamma *
      ( 3.0 * bessel_k(y, nuThres - 2.0, 1.0) 
	-y * bessel_k(y, nuThres - 3.0, 1.0)); 
  } else {
    v = 0.0;
  }
  v *= scaleSq * scale;
 
  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    scaleSq = scale * scale;
    y = x * scale;
    D3Gauss(&y, NULL, &w);
    w *= scaleSq * scale;
    v = v * g + (1.0 - g) * w;
  }
  return v;
}

double D4WM(double x,  double nu, double factor) { 
  static double nuOld=RF_INF;
  static double gamma;
  double y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    scale = (factor != 0.0) ? factor * sqrt(nuThres) : 1.0,
    scaleSq  = scale * scale;
  
  if (x > LOW_BESSELK) {
    if (nuThres!=nuOld) {
      nuOld = nuThres;
      gamma = gammafn(nuThres);
    }
    y = x * scale;
    v = 0.25 * pow(0.5 * y , nuThres - 3.0) / gamma *
      (+ 6.0 * (nuThres - 3.0 - y * y) * bessel_k(y, nuThres - 3.0, 1.0)
       + y * (3.0  + y * y) * bessel_k(y, nuThres - 4.0, 1.0)); 
  } else {
    v = (nuThres > 2.0) ? 0.75 / ((nuThres - 1.0) * (nuThres - 2.0)) : INFTY;
  }
  v *= scaleSq * scaleSq;

  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    scaleSq = scale * scale;
    y = x * scale;
    D4Gauss(&y, NULL, &w);
    w *= scaleSq * scaleSq;
     v = v * g + (1.0 - g) * w;
  }
  return v;
}


double ScaleWM(double nu){
  // it is working reasonably well if nu is in [0.001,100]
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
  
  return interpolate(log(nu) * INVLOG2, stuetz, nstuetz, stuetzorigin, 
		     1.5, 5);
}


int checkWM(cov_model *cov) { 
  static double
    spectrallimit=0.17,
    spectralbest=0.4;
  double notinvnu, nu;
  int i, err;
  bool isna_nu;

  if ((err = checkkappas(cov, false)) != NOERROR) return err;
  if (PisNULL(WM_NU)) QERRC(0, "parameter unset"); 

  nu = (PINT(WM_NOTINV) == NULL 
	|| ISNAN(notinvnu = (double) (P0INT(WM_NOTINV))) 
	|| notinvnu != 0.0) ? P0(WM_NU) : 1.0 / P0(WM_NU);
  isna_nu = ISNAN(nu);
  
  for (i=0; i<= Nothing; i++) cov->pref[i] *= isna_nu || nu < BesselUpperB[i];
  if (nu<spectralbest) {
    cov->pref[SpectralTBM] = (nu < spectrallimit) ? PREF_NONE : 3;
  } 
  if (cov->tsdim > 2) 
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (nu > 2.5)  cov->pref[CircEmbed] = 2;

  cov->full_derivs = isna_nu ? -1 : (int) (nu - 1.0);
  cov->monotone = nu <= 0.5 ? COMPLETELY_MON : NORMAL_MIXTURE;

  return NOERROR;
}

void rangeWM(cov_model *cov, range_type *range){
  bool tcf = isTcf(cov->typus) || cov->isoown == SPHERICAL_ISOTROPIC;
  if (tcf) {    
    if (PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)) {
      range->min[WM_NU] = 0.0;
      range->max[WM_NU] = 0.5;
      range->pmin[WM_NU] = 1e-1;
      range->pmax[WM_NU] = 0.5;
      range->openmin[WM_NU] = true;
      range->openmax[WM_NU] = false;
    } else {
      range->min[WM_NU] = 2.0;
      range->max[WM_NU] = RF_INF;
      range->pmin[WM_NU] = 2.0;
      range->pmax[WM_NU] = 10.0;
      range->openmin[WM_NU] = false;
      range->openmax[WM_NU] = true;
    }
  } else { // not tcf
    range->min[WM_NU] = 0.0;
    range->max[WM_NU] = RF_INF;
    range->pmin[WM_NU] = 1e-1;
    range->pmax[WM_NU] = 10.0;
    range->openmin[WM_NU] = true;
    range->openmax[WM_NU] = false;
  }
 
  range->min[WM_NOTINV] = 0.0;
  range->max[WM_NOTINV] = 1.0;
  range->pmin[WM_NOTINV] = 0.0;
  range->pmax[WM_NOTINV] = 1.0;
  range->openmin[WM_NOTINV] = false;
  range->openmax[WM_NOTINV] = false;
}


void coinitWM(cov_model *cov, localinfotype *li) {
  // cutoff_A
  double thres[2] = {0.25, 0.5},
    nu = (PINT(WM_NOTINV) == NULL || P0INT(WM_NOTINV)) 
             ? P0(WM_NU) : 1.0 / P0(WM_NU);
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

void ieinitWM(cov_model *cov, localinfotype *li) {
  double nu = (PINT(WM_NOTINV) == NULL || P0INT(WM_NOTINV))
    ? P0(WM_NU) : 1.0 / P0(WM_NU);
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

double densityWM(double *x, cov_model *cov, double factor) {
  double x2,
    nu = P0(WM_NU),
    powfactor = 1.0;
  int d,
    dim =  cov->tsdim;
  if (nu > 50)
    warning("nu>50 in density of matern class numerically instable. The results cannot be trusted.");
  if (factor == 0.0) factor = 1.0; else factor *= sqrt(nu);
  x2 = x[0] * x[0];
  for (d=1; d<dim; d++) {
    x2 += x[d] * x[d];
    powfactor *= factor;
  }
  x2 /= factor * factor;
  x2 += 1.0;
  
  return powfactor * exp(lgammafn(nu + 0.5 * (double) dim)
			 - lgammafn(nu)
			 - (double) dim * M_LN_SQRT_PI
			 - (nu + 0.5 * (double) dim) * log(x2));
}

  
/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */ 


void Matern(double *x, cov_model *cov, double *v) {
  *v = WM(*x, P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU), SQRT2);
}

void logMatern(double *x, cov_model *cov, double *v, double *Sign) {
  *v = logWM(*x, P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU), SQRT2);
  *Sign = 1.0;
}

void DMatern(double *x, cov_model *cov, double *v) {
  *v =DWM(*x, P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU), SQRT2);
} 

void DDMatern(double *x, cov_model *cov, double *v) {
  *v=DDWM(*x, P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU), SQRT2);
} 

void D3Matern(double *x, cov_model *cov, double *v) {
  *v=D3WM(*x, P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU), SQRT2);
} 

void D4Matern(double *x, cov_model *cov, double *v) {
  *v=D4WM(*x, P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU), SQRT2);
} 

void InverseMatern(double *x, cov_model *cov, double *v) {
  double
    nu = P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU);
  *v =  *x == 0.05 ?  SQRT2 * sqrt(nu) /  ScaleWM(nu) : RF_NA;
}

int checkMatern(cov_model *cov) { 
  if (PINT(WM_NOTINV) == NULL) kdefault(cov, WM_NOTINV, true);
  return checkWM(cov);
}

double densityMatern(double *x, cov_model *cov) {
  return densityWM(x, cov, SQRT2);
}

int initMatern(cov_model *cov, gen_storage *s) {
  if (HAS_SPECTRAL_ROLE(cov)) {
    spec_properties *cs = &(s->spec);
    if (cov->tsdim <= 2) return NOERROR;
    cs->density = densityMatern;
    return search_metropolis(cov, s);
  }

  else ILLEGAL_ROLE;

}

void spectralMatern(cov_model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  /* see Yaglom ! */
  if (cov->tsdim <= 2) {
    double nu = P0INT(WM_NOTINV) ? P0(WM_NU) : 1.0 / P0(WM_NU);
    E12(s, cov->tsdim, 
	sqrt( 2.0 * nu * (pow(1.0 - UNIFORM_RANDOM, -1.0 / nu) - 1.0) ), e);
  } else {
    metropolis(cov, S, e);
  }
}



// Brownian motion 
#define OESTING_BETA 0
void oesting(double *x, cov_model *cov, double *v) {
  // klein beta interagiert mit 'define beta Rf_beta' in Rmath
  double Beta = P0(OESTING_BETA),
    x2 = *x * *x;
  *v = - x2 * pow(1 + x2, -Beta);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
/* oesting: first derivative at t=1 */
void Doesting(double *x, cov_model *cov, double *v) 
{// FALSE VALUE FOR *x==0 and  Beta < 1
  double Beta = P0(OESTING_BETA),
    x2 = *x * *x;
  *v = 2 * *x  * (1 + (1-Beta) * x2) * pow(1 + x2, -Beta-1);
}
/* oesting: second derivative at t=1 */
void DDoesting(double *x, cov_model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  beta < 2
  double Beta = P0(OESTING_BETA),
    x2 = *x * *x;
  *v = 2 * (1 + (2 - 5 * Beta) * x2 + (1 - 3 * Beta + 2 * Beta * Beta) * 
	    x2  *x2) * pow(1 + x2, -Beta-2);
}
int checkoesting(cov_model *cov){
  cov->logspeed = RF_INF;
  cov->full_derivs = cov->rese_derivs;
  return NOERROR;
}
int initoesting(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double Beta = P0(OESTING_BETA);
  cov->taylor[1][TaylorConst] = + Beta;
  cov->taylor[2][TaylorConst] = - 0.5 * Beta * (Beta + 1);  
  cov->tail[0][TaylorPow] = 2 - 2 * Beta;
  return NOERROR;
}
void rangeoesting(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[OESTING_BETA] = 0.0;
  range->max[OESTING_BETA] = 1.0;
  range->pmin[OESTING_BETA] = UNIT_EPSILON;
  range->pmax[OESTING_BETA] = 1.0;
  range->openmin[OESTING_BETA] = true;
  range->openmax[OESTING_BETA] = false;
}


/* penta */
void penta(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v)
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
void Dpenta(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v)
{ 
  double y=*x, y2 = y * y;
  *v = (y >= 1.0) ? 0.0 :
    y * (-14.66666666666666667 + 
	 y2 * (132.0 + 
	       y * (-192.5 + 
		    y2 * (115.5 + 
			  y2 * (-49.5 + 
				y2 * 9.16666666666666667)))));
  
}
void Inversepenta(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = (*x==0.05) ? 1.0 / 1.6552838957365 : RF_NA;
}


/* power model */ 
#define POW_ALPHA 0
void power(double *x, cov_model *cov, double *v){
  double alpha = P0(POW_ALPHA), y = *x;
  *v = (y >= 1.0) ? 0.0 : pow(1.0 - y, alpha);
}
void TBM2power(double *x, cov_model *cov, double *v){
  // only alpha=2 up to now !
  double y = *x;
  if (P0(POW_ALPHA) != 2.0) 
    ERR("TBM2 of power only allowed for alpha=2");
  *v = (y > 1.0)  
    ? (1.0 - 2.0 * y *(asin(1.0 / y) - y + sqrt(y * y - 1.0) ))
    : 1.0 - y * (PI - 2.0 * y);
}
void Dpower(double *x, cov_model *cov, double *v){
  double alpha = P0(POW_ALPHA), y = *x;
  *v = (y >= 1.0) ? 0.0 : - alpha * pow(1.0 - y, alpha - 1.0);
}
int checkpower(cov_model *cov) {
  double alpha = P0(POW_ALPHA);
  double dim = 2.0 * alpha - 1.0;
  cov->maxdim = (ISNAN(dim) || dim >= INFDIM) 
    ? INFDIM-1 : (int) dim;
  cov->monotone = alpha >= (double) ((cov->tsdim / 2) + 1)
    ? COMPLETELY_MON : NORMAL_MIXTURE;
  return NOERROR;
}


// range definition:
// 0: min, theory, 1:max, theory
// 2: min, practically 3:max, practically
void rangepower(cov_model *cov, range_type *range){
  bool tcf = isTcf(cov->typus) || cov->isoown == SPHERICAL_ISOTROPIC;
  range->min[POW_ALPHA] = tcf
    ? (double) ((cov->tsdim / 2) + 1) // Formel stimmt so! Nicht aendern!
    : 0.5 * (double) (cov->tsdim + 1);
  range->max[POW_ALPHA] = RF_INF;
  range->pmin[POW_ALPHA] = range->min[POW_ALPHA];
  range->pmax[POW_ALPHA] = 20.0;
  range->openmin[POW_ALPHA] = false;
  range->openmax[POW_ALPHA] = true;
}


/* qexponential -- derivative of exponential */
#define QEXP_ALPHA 0
void qexponential(double *x, cov_model *cov, double *v){
  double 
    alpha = P0(QEXP_ALPHA),
    y = exp(-*x);
  *v = y * (2.0  - alpha * y) / (2.0 - alpha);
}
void Inverseqexponential(double *x, cov_model *cov, double *v){
  double alpha = P0(QEXP_ALPHA);
  *v = -log( (1.0 - sqrt(1.0 - alpha * (2.0 - alpha) * *x)) / alpha);
} 
void Dqexponential(double *x, cov_model *cov, double *v) {
  double 
    alpha = P0(QEXP_ALPHA), 
    y = exp(-*x);
  *v = y * (alpha * y - 1.0) * 2.0 / (2.0 - alpha);
}
void rangeqexponential(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[QEXP_ALPHA] = 0.0;
  range->max[QEXP_ALPHA] = 1.0;
  range->pmin[QEXP_ALPHA] = 0.0;
  range->pmax[QEXP_ALPHA] = 1.0;
  range->openmin[QEXP_ALPHA] = false;
  range->openmax[QEXP_ALPHA] = false;
}


/* spherical model */ 
void spherical(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x;
  *v = (y >= 1.0) ? 0.0 : 1.0 + y * 0.5 * (y * y - 3.0);
}
// void Inversespherical(cov_model *cov){ return 1.23243208931941;}
void TBM2spherical(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x , y2 = y * y;
  *v = (y>1.0) 
      ? (1.0- 0.75 * y * ((2.0 - y2) * asin(1.0/y) + sqrt(y2 -1.0)))
      : (1.0 - 0.375 * PI * y * (2.0 - y2));
}
void Dspherical(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x;
  *v = (y >= 1.0) ? 0.0 : 1.5 * (y * y - 1.0);
}

void DDspherical(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  *v = (*x >= 1.0) ? 0.0 : 3 * *x;
}

int structspherical(cov_model *cov, cov_model **newmodel) {
  return structCircSph( cov, newmodel, 3);
}
void dospherical(cov_model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 
  // todo diese void do... in Primitive necessary??
  //  mppinfotype *info = &(s->mppinfo);
  //info->radius = cov->mpp.refradius;
}


double alphaIntSpherical(int a) {
  // int r^a C(r) \D r
  double A = (double) a;
  return 1.0 / (A + 1.0) - 1.5 / (A + 2.0) + 0.5 / (A + 4.0);
}

int initspherical(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int  dim = cov->tsdim;

  if (hasNoRole(cov)) {
    if (cov->mpp.moments >= 1) SERR("too high moments required");
    return NOERROR;
  }

  if (hasAnyShapeRole(cov)) {

    assert(cov->mpp.maxheights[0] == 1.0);
    if (cov->mpp.moments >= 1) {
      cov->mpp.mM[1] = cov->mpp.mMplus[1] =  
	SurfaceSphere(dim - 1, 1) * alphaIntSpherical(dim - 1);
    }

  }

  else ILLEGAL_ROLE;

 return NOERROR;
}

/* stable model */
#define STABLE_ALPHA 0
void stable(double *x, cov_model *cov, double *v){
  double y = *x, alpha = P0(STABLE_ALPHA);  
  *v = (y==0.0) ? 1.0 : exp(-pow(y, alpha));
}
void logstable(double *x, cov_model *cov, double *v, double *Sign){
  double y = *x, alpha = P0(STABLE_ALPHA);  
  *v = (y==0.0) ? 0.0 : -pow(y, alpha);
  *Sign = 1.0;
}
void Dstable(double *x, cov_model *cov, double *v){
  double z, y = *x, alpha = P0(STABLE_ALPHA);
  if (y == 0.0) {
    *v = (alpha > 1.0) ? 0.0 : (alpha < 1.0) ? INFTY : 1.0;
  } else {
    z = pow(y, alpha - 1.0);
    *v = -alpha * z * exp(-z * y);
  }
}
/* stable: second derivative at t=1 */
void DDstable(double *x, cov_model *cov, double *v) 
{
  double z, y = *x, alpha = P0(STABLE_ALPHA), xalpha;
  if (y == 0.0) {
      *v = (alpha == 1.0) ? 1.0 : (alpha == 2.0) ? alpha * (1 - alpha) : INFTY;
  } else {
    z = pow(y, alpha - 2.0);
    xalpha = z * y * y;
    *v = alpha * (1.0 - alpha + alpha * xalpha) * z * exp(-xalpha);
  }
}
void Inversestable(double *x, cov_model *cov, double *v){
  double y = *x, alpha = P0(STABLE_ALPHA);  
  *v = y>1.0 ? 0.0 : y == 0.0 ? RF_INF : pow( - log(y), 1.0 / alpha);
}
void nonstatLogInversestable(double *x, cov_model *cov,
			     double *left, double *right){
  double 
    alpha = P0(STABLE_ALPHA),
    z = *x <= 0.0 ? pow( - *x, 1.0 / alpha) : 0.0;
  int d,
    dim = cov->tsdim;
  for (d=0; d<dim; d++) {
    left[d] = -z;
    right[d] = z;
  }
}


int checkstable(cov_model *cov) {
  double alpha = P0(STABLE_ALPHA);
  if (cov->tsdim > 2)
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (alpha == 2.0)
    cov->pref[CircEmbed] = 2;
  
  cov->monotone = alpha <= 1.0 ? COMPLETELY_MON : NORMAL_MIXTURE;

  return NOERROR;
}
void rangestable(cov_model *cov, range_type *range){
  bool tcf = isTcf(cov->typus) || cov->isoown == SPHERICAL_ISOTROPIC;
  range->min[STABLE_ALPHA] = 0.0;
  range->max[STABLE_ALPHA] = tcf ? 1.0 : 2.0;
  range->pmin[STABLE_ALPHA] = 0.06;
  range->pmax[STABLE_ALPHA] =  range->max[STABLE_ALPHA];
  range->openmin[STABLE_ALPHA] = true;
  range->openmax[STABLE_ALPHA] = false;
}
void coinitstable(cov_model *cov, localinfotype *li) {
  coinitgenCauchy(cov, li);
}
void ieinitstable(cov_model *cov, localinfotype *li) {  
  ieinitgenCauchy(cov, li);
}


/* SPACEISOTROPIC stable model for testing purposes only */
void stableX(double *x, cov_model *cov, double *v){
  double y, alpha = P0(STABLE_ALPHA);
  y = x[0] * x[0] + x[1] * x[1];
  *v = (y==0.0) ? 1.0 : exp(-pow(y, 0.5 * alpha));
}
void DstableX(double *x, cov_model *cov, double *v){
  double z, y, alpha = P0(STABLE_ALPHA);
  y = x[0] * x[0] + x[1] * x[1];
  if (y == 0.0) {
    *v = ((alpha > 1.0) ? 0.0 : (alpha < 1.0) ? INFTY : 1.0);
  } else {
    z = pow(y, 0.5 * alpha - 1.0);
    *v = -alpha * x[0] * z * exp(- z * y);
  }
}
/* END SPACEISOTROPIC stable model for testing purposes only */




/* stein space-time model */
#define STEIN_NU 0
#define STEIN_Z 1
void kappaSteinST1(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = (i == STEIN_NU) ? 1 : (i==STEIN_Z) ?  cov->tsdim - 1 : -1;
}
void SteinST1(double *x, cov_model *cov, double *v){
/* 2^(1-nu)/Gamma(nu) [ h^nu K_nu(h) - 2 * tau (x T z) t h^{nu-1} K_{nu-1}(h) /
   (2 nu + d + 1) ]

   \|tau z \|<=1 hence \tau z replaced by z !!
*/
  int d,
    dim = cov->tsdim,
    time = dim - 1;
  double logconst, hz, y,
    nu = P0(STEIN_NU),
    *z=P(STEIN_Z);
  
  static double nuold=RF_INF;
  static double loggamma;
  static int dimold;

  if (nu != nuold || dimold != dim) {
    nuold = nu;
    dimold = dim;
    loggamma = lgammafn(nu);
  }

  hz = 0.0;
  y = x[time] * x[time];
  for (d=0; d<time; d++) {
    y += x[d] * x[d];
    hz += x[d] * z[d];
  }
  
  if ( y==0.0 ) *v = 1.0;
  else {
    y = sqrt(y);
    logconst = (nu - 1.0) * log(0.5 * y)  - loggamma;
    *v =  y * exp(logconst + log(bessel_k(y, nu, 2.0)) - y)
      - 2.0 * hz * x[time] * exp(logconst + log(bessel_k(y, nu - 1.0, 2.0)) -y) 
      / (2.0 * nu + dim);
  }

}

int checkSteinST1(cov_model *cov) {  
  double nu = P0(STEIN_NU), *z= P(STEIN_Z), absz;
  int d, spatialdim=cov->tsdim-1;

  for (d=0; d<= Nothing; d++) cov->pref[d] *= (nu < BesselUpperB[d]);
  if (nu >= 2.5) cov->pref[CircEmbed] = 2;
  if (spatialdim < 1) 
    SERR("dimension of coordinates, including time, must be at least 2");
  	
  for (absz=0.0, d=0;  d<spatialdim; d++)  absz += z[d] * z[d];
  if (ISNAN(absz))
    SERR("currently, components of z cannot be estimated by MLE, so NA's are not allowed");
  if (absz > 1.0 + UNIT_EPSILON && !GLOBAL.general.skipchecks) {
    SERR("||z|| must be less than or equal to 1");
  }
  return NOERROR;
}
double densitySteinST1(double *x, cov_model *cov) {
  double x2, wz, dWM, nu = P0(STEIN_NU), 
    *z=P(STEIN_Z);
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


int initSteinST1(cov_model *cov, gen_storage *s) {
  if (HAS_SPECTRAL_ROLE(cov)) {
    spec_properties *cs = &(s->spec);
    cs->density = densitySteinST1;
    return search_metropolis(cov, s);
  }

  else ILLEGAL_ROLE;

}
void spectralSteinST1(cov_model *cov, gen_storage *S, double *e){
  metropolis(cov, S, e);
}

void rangeSteinST1(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){  
  range->min[STEIN_NU] = 1.0;
  range->max[STEIN_NU] = RF_INF;
  range->pmin[STEIN_NU] = 1e-10;
  range->pmax[STEIN_NU] = 10.0;
  range->openmin[STEIN_NU] = true;
  range->openmax[STEIN_NU] = true;
 
  range->min[STEIN_Z] = RF_NEGINF;
  range->max[STEIN_Z] = RF_INF;
  range->pmin[STEIN_Z] = -10.0;
  range->pmax[STEIN_Z] = 10.0;
  range->openmin[STEIN_Z] = true;
  range->openmax[STEIN_Z] = true;
}


/* wave */
void wave(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x;
  *v = (y==0.0) ? 1.0 : (y==RF_INF) ? 0 : sin(y) / y;
}
void Inversewave(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = (*x==0.05) ? 1.0/0.302320850755833 : RF_NA;
}
int initwave(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  if (HAS_SPECTRAL_ROLE(cov)) {
    return (cov->tsdim <= 2) ? NOERROR : ERRORFAILED;
  } 

  else ILLEGAL_ROLE;

}
void spectralwave(cov_model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  /* see Yaglom ! */
  double x;  
  x = UNIFORM_RANDOM; 
  E12(s, cov->tsdim, sqrt(1.0 - x * x), e);
}



void Whittle(double *x, cov_model *cov, double *v) {
  *v = WM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV) 
	  ? P0(WM_NU) : 1.0 / P0(WM_NU), 0.0);
  assert(!ISNAN(*v));
}


void logWhittle(double *x, cov_model *cov, double *v, double *Sign) {
  *v = logWM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)
	     ? P0(WM_NU)
	     : 1.0 / P0(WM_NU), 0.0);
  assert(!ISNA(*v));
  *Sign = 1.0;
}

void DWhittle(double *x, cov_model *cov, double *v) {
  *v =DWM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)
	  ? P0(WM_NU) : 1.0 / P0(WM_NU), 0.0);
}

void DDWhittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=DDWM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV) 
	  ? P0(WM_NU) : 1.0 / P0(WM_NU), 0.0);
}


void D3Whittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=D3WM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)
	  ? P0(WM_NU) : 1.0 / P0(WM_NU),
	  PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}

void D4Whittle(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=D4WM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV) 
	  ? P0(WM_NU) : 1.0 / P0(WM_NU),
	  PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}

void InverseWhittle(double *x, cov_model *cov, double *v){
  double 
    nu = (PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)) ? P0(WM_NU) : 1.0 / P0(WM_NU);
  *v = (*x == 0.05) ? 1.0 / ScaleWM(nu) : RF_NA;
}

void TBM2Whittle(double *x, cov_model *cov, double *v) 
{
  double nu = P0(WM_NU);
  assert(PisNULL(WM_NOTINV));
  if (nu == 0.5) TBM2exponential(x, cov, v);
  else BUG;
}


double densityWhittle(double *x, cov_model *cov) {
  return densityWM(x, cov, PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}
int initWhittle(cov_model *cov, gen_storage *s) {
  if (HAS_SPECTRAL_ROLE(cov)) {
    if (PisNULL(WM_NU)) {
      spec_properties *cs = &(s->spec);
      if (cov->tsdim <= 2) return NOERROR;
      cs->density = densityWhittle; 
      int err=search_metropolis(cov, s);
      return err;
    } else return initMatern(cov, s);
  }

  else ILLEGAL_ROLE;

}

void spectralWhittle(cov_model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  /* see Yaglom ! */
  if (PisNULL(WM_NOTINV)) {
    if (cov->tsdim <= 2) {
      double nu = P0(WM_NU);
      E12(s, cov->tsdim, sqrt(pow(1.0 - UNIFORM_RANDOM, -1.0 / nu) - 1.0), e);
    } else {
      metropolis(cov, S, e);
    }
  } else spectralMatern(cov, S, e);
}


void DrawMixWM(cov_model VARIABLE_IS_NOT_USED *cov, double *random) { // inv scale
  // V ~ F in PSgen
  *random = -0.25 / log(UNIFORM_RANDOM);
}

double LogMixDensW(double VARIABLE_IS_NOT_USED *x, double logV, cov_model *cov) {
  double
    nu=P0(WM_NU);
  return PisNULL(WM_NOTINV)
    ? (M_LN2  + 0.5 * logV) * (1.0 - nu) - 0.5 *lgammafn(nu)
    // - 0.25 /  cov->mpp[DRAWMIX_V]  - 2.0 * (LOG2 + logV) )
    : RF_NA;  /* !! */
}


void Whittle2(double *x, cov_model *cov, double *v) {
  *v = WM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV) 
	  ? P0(WM_NU) : 1.0 / P0(WM_NU), 0.0);
  assert(!ISNAN(*v));
}


void logWhittle2(double *x, cov_model *cov, double *v, double *Sign) {
  *v = logWM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)
	     ? P0(WM_NU)
	     : 1.0 / P0(WM_NU), 0.0);
  assert(!ISNA(*v));
  *Sign = 1.0;
}

void DWhittle2(double *x, cov_model *cov, double *v) {
  *v =DWM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)
	  ? P0(WM_NU) : 1.0 / P0(WM_NU), 0.0);
}

void DDWhittle2(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=DDWM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV) 
	  ? P0(WM_NU) : 1.0 / P0(WM_NU), 0.0);
}


void D3Whittle2(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=D3WM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)
	  ? P0(WM_NU) : 1.0 / P0(WM_NU),
	  PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}

void D4Whittle2(double *x, cov_model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=D4WM(*x, PisNULL(WM_NOTINV) || P0INT(WM_NOTINV) 
	  ? P0(WM_NU) : 1.0 / P0(WM_NU),
	  PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}

void InverseWhittle2(double *x, cov_model *cov, double *v){
  double 
    nu = (PisNULL(WM_NOTINV) || P0INT(WM_NOTINV)) ? P0(WM_NU) : 1.0 / P0(WM_NU);
  *v = (*x == 0.05) ? 1.0 / ScaleWM(nu) : RF_NA;
}


/* Whittle-Matern or Whittle or Besset */ 
/// only lower parts of the matrices, including the diagonals are used when estimating !!

static bool Bi = !true;

#define BInudiag 0
#define BInured 1
#define BInu 2
#define BIs 3
#define BIcdiag 4
#define BIrhored 5
#define BIc 6
#define BInotinvnu 7


/* Whittle-Matern or Whittle or Besset */ 

void biWM2basic(cov_model *cov, 
		double *a, double *lg,
		double *aorig, double *nunew) {
  // !! nudiag, nured has priority over nu
  // !! cdiag, rhored has priority over c

  double factor, beta, gamma, tsq, t1sq, t2sq, inf, infQ, discr,
    alpha , a2[3],
    dim = (double) cov->tsdim, 
    d2 = dim * 0.5, 
    *nudiag = P(BInudiag),
    nured = P0(BInured),
    *nu = P(BInu),
    *c = P(BIc),
    *s = P(BIs),
    *cdiag = P(BIcdiag),
    rho = P0(BIrhored);
  int i, 
    *notinvnu = PINT(BInotinvnu);

  nu[i11] = nudiag[0];
  nu[i22] = nudiag[1];
  nu[i21] = 0.5 * (nu[i11] + nu[i22]) * nured;
  
  for (i=0; i<3; i++) {
    aorig[i] = 1.0 / s[i];
    if (Bi) print("%d %f %f \n", i, s[i], aorig[i]);
  } 

  if (PisNULL(BInotinvnu)) {
    for (i=0; i<3; i++) {
      a[i] = aorig[i];
      nunew[i] = nu[i];
    }
  } else {
    if (!notinvnu[0]) for (i=0; i<3; i++) nu[i] = 1.0 / nu[i];
    for (i=0; i<3; i++) {
      nunew[i] = nu[i] < MATERN_NU_THRES ? nu[i] : MATERN_NU_THRES;
      a[i] = aorig[i] * sqrt(2.0 * nunew[i]);
    }
  }


  for (i=0; i<3; i++) {
    a2[i] = a[i] * a[i];
    lg[i] = lgammafn(nunew[i]);
  }
  
  alpha = 2 * nunew[i21] - nunew[i11] - nunew[i22];


  // **************** ACHTUNG !!!!! *********
  // nicht gut, sollte in check einmal berechnet werden 
  // dies wiederspricht aber der MLE Maximierung, da dann
  // neu berechnet werden muss; verlg. natsc und cutoff (wo es nicht
  // programmiert ist !!)
  
  factor = exp(lgammafn(nunew[i11] + d2) - lg[i11]
	       + lgammafn(nunew[i22] + d2) - lg[i22]
		   + 2.0 * (lg[i21]  - lgammafn(nunew[i21] + d2)
			    + nunew[i11] * log(a[i11]) + nunew[i22] *log(a[i22])
			    - 2.0 * nunew[i21] * log(a[i21]))
	);

  
  // alpha u^2 + beta u + gamma = 0
  gamma = (2.0 * nunew[i21] + dim) * a2[i11] * a2[i22] 
    - (nunew[i22] + d2) * a2[i11] * a2[i21]
    - (nunew[i11] + d2) * a2[i22] * a2[i21];
      
  beta = (2.0 * nunew[i21] - nunew[i22] + d2) * a2[i11]
      + (2.0 * nunew[i21] - nunew[i11] + d2) * a2[i22]
      - (nunew[i11] + nunew[i22] + dim) * a2[i21];
  
  if (Bi) print("%f %f %f %f %f\n"
		, 2.0 * nunew[i21], - nunew[i11], + d2 , a2[i22]
	    , (nunew[i11] + nunew[i22] + dim) * a2[i22]);

//  
  if (Bi) print("\nalpha=%f beta=%f gamma=%f\n", alpha, beta, gamma);
  if (Bi)  print("\nnu=%f %f %f, a2=%f %f %f\n", 
		  nunew[0], nunew[1], nunew[2], a2[0], a2[1], a2[2]);
     
  if (Bi) print("%d %f %f %f NU22 %f\n", i22, nu[0], nu[1], nu[2], nu[i22]);
     
  if (nured == 1.0) { // lin. eqn only
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


  inf = nured == 1.0 ? 1.0 : RF_INF; // t^2 = infty  ; nudiag[1]>1.0 not
  //                                   allowed by definition
  for (i=0; i<3; i++) {
    tsq = (i == 0) ? 0.0 
      : (i == 1) ? t1sq
      : t2sq;
    
    infQ = pow(a2[i21] + tsq, 2.0 * nunew[i21] + dim) /
      (pow(a2[i11] + tsq, nunew[i11] + d2) 
       * pow(a2[i22] + tsq, nunew[i22] + d2));

    if (infQ < inf) inf = infQ;
  }

  c[i11] = cdiag[0];
  c[i22] = cdiag[1];
  c[i21] = rho * sqrt(factor * inf * c[i11] *  c[i22]);
 
  if (Bi) print("c=%f %f %f rho=%f %f %f\n", c[0], c[1],c[2], rho, factor, inf);
  Bi = false;
 
}

int initbiWM2(cov_model *cov, gen_storage *stor) {
  double a[3], lg[3], aorig[3], nunew[3], 
    *c = P(BIc),
    *cdiag = P(BIcdiag),
    *nu = P(BInu),
    *nudiag = P(BInudiag);
  bool check = stor->check;
  int i,
    *notinvnu = PINT(BInotinvnu);
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  if (S->nudiag_given) {
    kdefault(cov, BInured, 1.0);
    if (check && !PisNULL(BInu)) {
      if (cov->nrow[BInu] != 3 || cov->ncol[BInu] != 1)
	QERRC(BInu, "must be a 3 x 1 vector");
      if (fabs(nu[i11] - nudiag[0]) > nu[i11] * epsilon || 
	  fabs(nu[i22] - nudiag[1]) > nu[i22] * epsilon ||
	  fabs(nu[i21] - 0.5 * (nudiag[i11] + nudiag[1]) * P0(BInured))
	  > nu[i21] * epsilon) {
	QERRC2(BInu, "does not match '%s' and '%s'.", BInudiag, BInured);
      }
    } else {
      if (PisNULL(BInu)) PALLOC(BInu, 3, 1);
      nu = P(BInu); 
      nu[i11] = nudiag[0];
      nu[i22] = nudiag[1];
      nu[i21] = 0.5 * (nu[i11] + nu[i22]) * P0(BInured);
    }
  } else {
    if (check && !PisNULL(BInured)) 
      QERRC1(BInured, "may not be set if '%s' is given", BInu);
    if (PisNULL(BInu)) 
      QERRC2(BInu, "either '%s' or '%s' must be set", BInu, BInured);
    if (PisNULL(BInudiag)) PALLOC(BInudiag, 2, 1);
    nudiag = P(BInudiag);
    nudiag[0] = nu[i11]; 
    nudiag[1] = nu[i22];
    if (PisNULL(BInured)) PALLOC(BInured, 1, 1);
    P(BInured)[0] =  nu[i21] / (0.5 * (nudiag[i11] + nudiag[1]));
  }

  if (!PisNULL(BInotinvnu) && !notinvnu[0]) 
    for (i=0; i<3; i++) nu[i] = 1.0 / nu[i];
 
  cov->full_derivs = cov->rese_derivs = 1; // kann auf 2 erhoeht werden, falls programmiert 
  for (i=0; i<3; i++) {
    int derivs = (int) (nu[i] - 1.0);
    if (cov->full_derivs < derivs) cov->full_derivs = derivs;
  }
  
  if (PisNULL(BIs)) {
    PALLOC(BIs, 3, 1);
    double *s = P(BIs);
    for (i=0; i<3; s[i++] = 1.0);
  }
  
  if  (S->cdiag_given) {
    if (PisNULL(BIrhored))
      QERRC2(BIrhored, "'%s' and '%s' must be set", BIcdiag, BIrhored);
    if (check && !PisNULL(BIc)) {
      if (cov->nrow[BIc] != 3 || cov->ncol[BIc] != 1)
	QERRC(BIc, "must be a 3 x 1 vector");
      if (fabs(c[i11] - cdiag[0]) > c[i11] * epsilon || 
	  fabs(c[i22] - cdiag[1]) > c[i22] * epsilon ) {
	QERRC1(BIc, "does not match '%s'.", BIcdiag);
      }
      double savec12 = c[i21];
      biWM2basic(cov, a, lg, aorig, nunew);
      if (fabs(c[i21] - savec12) > fabs(c[i21]) * epsilon)
 	QERRC1(BIc, "does not match '%s'.", BIrhored);
    } else {
      if (PisNULL(BIc)) PALLOC(BIc, 3, 1);
      c = P(BIc);
      biWM2basic(cov, a, lg, aorig, nunew);
   }
  } else {
    if (check && !PisNULL(BIrhored)) 
      QERRC1(BIrhored, "may not be set if '%s' is not given", BIcdiag);
    if (PisNULL(BIc)) QERRC2(BIc, "either '%s' or '%s' must be set",
			    BIc, BIcdiag);
    if (!ISNAN(c[i11]) && !ISNAN(c[i22]) && (c[i11] < 0.0 || c[i22] < 0.0))
      QERRC2(BIc, "variance parameter %s[0], %s[2] must be non-negative.",
	     BIc, BIc);
    if (PisNULL(BIcdiag)) PALLOC(BIcdiag, 2, 1);
    if (PisNULL(BIrhored)) PALLOC(BIrhored, 1, 1);
    cdiag = P(BIcdiag);
    cdiag[0] = c[i11];
    cdiag[1] = c[i22];
    double savec1 = c[i21];
    if (savec1 == 0.0)  P(BIrhored)[0] = 0.0; // wichtig falls c[0] oder c[2]=NA
    else {
      P(BIrhored)[0] = 1.0;
      biWM2basic(cov, a, lg, aorig, nunew);
      P(BIrhored)[0] = savec1 / c[i21];
    }
  }
  biWM2basic(cov, S->a, S->lg, S->aorig, S->nunew);

  cov->initialised = true;
  return NOERROR;
}


void kappa_biWM(int i, cov_model *cov, int *nr, int *nc){
  *nc = *nr = i < CovList[cov->nr].kappas ? 1 : -1;
  if (i==BInudiag || i==BIcdiag) *nr = 2; else
    if (i==BInu || i==BIs || i==BIc) *nr = 3;
}

void biWM2(double *x, cov_model *cov, double *v) {
  int i;
  double 
    *c = P(BIc),
    *nu = P(BInu),
    xx = *x;
  biwm_storage *S = cov->Sbiwm;

  assert(S != NULL);

  assert(cov->initialised);

  for (i=0; i<3; i++) {
    v[i] = c[i] * WM(fabs(S->a[i] * xx), S->nunew[i], 0.0);
    if (!PisNULL(BInotinvnu) && nu[i] > MATERN_NU_THRES) {
      double w, y;
      y = fabs(S->aorig[i] * xx * INVSQRTTWO);
      Gauss(&y, cov, &w);
      *v = *v * MATERN_NU_THRES / nu[i] + 
	(1 - MATERN_NU_THRES / nu[i]) * w;
    }
  }
  v[3] = v[i22];
  v[2] = v[i21];

}

void biWM2D(double *x, cov_model *cov, double *v) {
  int i;
  double
    *c = P(BIc),
    *nu = P(BInu),
    xx = *x;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  assert(cov->initialised);

  for (i=0; i<3; i++) {
    v[i] = c[i] * S->a[i] * DWM(fabs(S->a[i] * xx), S->nunew[i], 0.0);
    if (!PisNULL(BInotinvnu) && nu[i] > MATERN_NU_THRES) {
      double w, y,
	scale = S->aorig[i] * INVSQRTTWO;
      y = fabs(scale * xx);
      DGauss(&y, cov, &w);
      w *= scale;
      *v = *v * MATERN_NU_THRES / nu[i] + 
	(1 - MATERN_NU_THRES / nu[i]) * w;
    }
  }
  v[3] = v[i22];
  v[2] = v[i21];
}


int checkbiWM2(cov_model *cov) { 

  // !! nudiag, nured has priority over nu
  // !! cdiag, rhored has priority over c
  int err;
  gen_storage s;
  gen_NULL(&s);
  s.check = true;
 
  assert(PisNULL(BIrhored) || ISNAN(P0(BIrhored)) || P0(BIrhored) <= 1);
  if ((err = checkkappas(cov, false)) != NOERROR) return err;
  NEW_STORAGE(biwm);
  biwm_storage *S = cov->Sbiwm;
  S->nudiag_given = !PisNULL(BInudiag);
  S->cdiag_given = !PisNULL(BIcdiag);
  if ((err=initbiWM2(cov, &s)) != NOERROR) return err;

  cov->vdim[0] = cov->vdim[1] = 2;
 
  return NOERROR;
}

sortsofparam paramtype_biWM(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  return  ( k== BInudiag || k==BInured) ? CRITICALPARAM 
    : (k == BInu || k == BIc) ? DONOTRETURNPARAM
    : (k == BIs) ? SCALEPARAM 
    : ( k == BIcdiag) ? VARPARAM : ANYPARAM; // c ignored
}
  
void rangebiWM2(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){

  // *nudiag = P0(BInudiag], 
  range->min[BInudiag] = 0.0;
  range->max[BInudiag] = RF_INF;
  range->pmin[BInudiag] = 1e-1;
  range->pmax[BInudiag] = 4.0; 
  range->openmin[BInudiag] = true;
  range->openmax[BInudiag] = true;
  
  // *nured12 = P0(BInured], 
  range->min[BInured] = 1;
  range->max[BInured] = RF_INF;
  range->pmin[BInured] = 1;
  range->pmax[BInured] = 5;
  range->openmin[BInured] = false;
  range->openmax[BInured] = true;
    
 // *nu = P0(BInu], 
  range->min[BInu] = 0.0;
  range->max[BInu] = RF_INF;
  range->pmin[BInu] = 1e-1;
  range->pmax[BInu] = 4.0;
  range->openmin[BInu] = true;
  range->openmax[BInu] = true;
 
 //   s = P0(BIs], 
  range->min[BIs] = 0.0;
  range->max[BIs] = RF_INF;
  range->pmin[BIs] = 1e-2;
  range->pmax[BIs] = 100.0;
  range->openmin[BIs] = true;
  range->openmax[BIs] = true;
  
  //    *cdiag = P0(BIcdiag]; 
  range->min[BIcdiag] = 0;
  range->max[BIcdiag] = RF_INF;
  range->pmin[BIcdiag] = 0.001;
  range->pmax[BIcdiag] = 1000;
  range->openmin[BIcdiag] = false;
  range->openmax[BIcdiag] = true;  
  
 //    *rho = P0(BIrhored]; 
  range->min[BIrhored] = -1;
  range->max[BIrhored] = 1;
  range->pmin[BIrhored] = -0.99;
  range->pmax[BIrhored] = 0.99;
  range->openmin[BIrhored] = false;
  range->openmax[BIrhored] = false;    

  //    *c = P0(BIc]; 
  range->min[BIc] = RF_NEGINF;
  range->max[BIc] = RF_INF;
  range->pmin[BIc] = 0.001;
  range->pmax[BIc] = 1000;
  range->openmin[BIc] = false;
  range->openmax[BIc] = true;  
  

//    *notinvnu = P0(BInotinvnu]; 
  range->min[BInotinvnu] = 0;
  range->max[BInotinvnu] = 1;
  range->pmin[BInotinvnu] = 0;
  range->pmax[BInotinvnu] = 1;
  range->openmin[BInotinvnu] = false;
  range->openmax[BInotinvnu] = false;    

}





// PARS WM


/* Whittle-Matern or Whittle or Besset */ 
// only lower parts of the matrices, including the diagonals are used when estimating !!

#define PARSnudiag 0
void kappa_parsWM(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i==PARSnudiag) {
    *nr = 0; 
    *nc = 1;
  } else  *nc = *nr = -1;
}

void parsWMbasic(cov_model *cov) {
  // !! nudiag, nured has priority over nu
  // !! cdiag, rhored has priority over c

  double 
    dim = (double) cov->tsdim, 
    d2 = dim * 0.5, 
    *nudiag = P(PARSnudiag);
  int i, j, vdiag,
    vdim = cov->nrow[PARSnudiag], // hier noch nicht cov->vdim, falls parsWMbasic von check aufgerufen wird, und cov->vdim noch nicht gesetzt ist
    vdimP1 = vdim + 1;
 
  // **************** ACHTUNG !!!!! *********
  // nicht gut, sollte in check einmal berechnet werden 
  // dies wiederspricht aber der MLE Maximierung, da dann
  // neu berechnet werden muss; verlg. natsc und cutoff (wo es nicht
  // programmiert ist !!)

  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    cov->q[vdiag] = 1.0;
    for (j=i+1; j<vdim; j++) {
      double half = 0.5 * (nudiag[i] + nudiag[j]);
      int idx = vdiag + j - i;
      cov->q[idx] = cov->q[vdiag + vdim * (j-i)] =
	exp(0.5 * (lgammafn(nudiag[i] + d2) + lgammafn(nudiag[j] + d2)
		   - lgammafn(nudiag[i]) - lgammafn(nudiag[j])
		   + 2.0 * (lgammafn(half) - lgammafn(half + d2))
		   ));
    }
  }
}

void parsWM(double *x, cov_model *cov, double *v) {
  int i, j, vdiag,
    vdim = cov->vdim[0],
    vdimP1 = vdim + 1;
  double 
    *nudiag = P(PARSnudiag);

  parsWMbasic(cov);
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    for (j=i; j<vdim; j++) {
      double half = 0.5 * (nudiag[i] + nudiag[j]);      
      int idx = vdiag + j - i;
      v[idx] = v[vdiag + vdim * (j-i)] = WM(*x, half, 0.0) * cov->q[idx];
    }
  }
}


void parsWMD(double *x, cov_model *cov, double *v) {
  int i, j, vdiag,
    vdim = cov->vdim[0],
    vdimP1 = vdim + 1;
  double 
    *nudiag = P(PARSnudiag);
  parsWMbasic(cov);
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    for (j=i; j<vdim; j++) {
      double half = 0.5 * (nudiag[i] + nudiag[j]);      
      int idx = vdiag + j - i;
      v[idx] = v[vdiag + vdim * (j-i)] = DWM(*x, half, 0.0) * cov->q[idx];
    }
  }
}


int checkparsWM(cov_model *cov) { 
 
  double
    *nudiag = P(PARSnudiag);
  
  int i, err, 
    vdim = cov->nrow[PARSnudiag],
    vdimSq = vdim * vdim;
 
  cov->vdim[0] = cov->vdim[1] = vdim;
  if (vdim == 0) SERR1("'%s' not given", KNAME(PARSnudiag));
  if ((err = checkkappas(cov, true)) != NOERROR) return err;
  if (cov->q == NULL) QALLOC(vdimSq);
  
  cov->full_derivs = cov->rese_derivs = 1; 
  for (i=0; i<vdim; i++) {
    int derivs = (int) (nudiag[i] - 1.0);
    if (cov->full_derivs < derivs) cov->full_derivs = derivs;
  }

  /*
    #define dummyN (5 * ParsWMMaxVDim)
    double value[ParsWMMaxVDim], ivalue[ParsWMMaxVDim], 
    dummy[dummyN], 
    min = RF_INF;
    int j,
    ndummy = dummyN;
  */

  
  return NOERROR;
}

sortsofparam paramtype_parsWM(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  return  ( k== PARSnudiag) ? CRITICALPARAM : ANYPARAM; // c ignored
}
  
void rangeparsWM(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){

  // *nudiag = P0(PARSnudiag], 
  range->min[PARSnudiag] = 0.0;
  range->max[PARSnudiag] = RF_INF;
  range->pmin[PARSnudiag] = 1e-1;
  range->pmax[PARSnudiag] = 4.0; 
  range->openmin[PARSnudiag] = true;
  range->openmax[PARSnudiag] = true;
}





// ************************* NOT A COVARIANCE FCT

double SurfaceSphere(int d, double r) { 
    // d = Hausdorff-Dimension der Oberflaeche der Sphaere
   //  NOT   2 \frac{\pi^{d/2}}{\Gamma(d/2)} r^{d-1}, 
   //  BUT  2 \frac{\pi^{(d+1)/2}}{\Gamma((d+1)/2)} r^{d}, 
   double D = (double) d;
  // printf("r=%f, %f %f %f\n", r, D, pow(SQRTPI * r, D - 1.0), gammafn(0.5 * D));

   return 2.0 * SQRTPI * pow(SQRTPI * r, D) / gammafn(0.5 * (D + 1.0));

}

double VolumeBall(int d, double r) {
  //  V_n(R) = \frac{\pi^{d/2}}{\Gamma(\frac{d}{2} + 1)}R^n, 
 double D = (double) d;
 return pow(SQRTPI * r, D) / gammafn(0.5 * D + 1.0);  
}


void ball(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) { 
  // isotropic function, expecting only distance; BALL_RADIUS=1.0
  assert(*x >= 0);
  *v = (double) (*x <= BALL_RADIUS);
}

void Inverseball(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){
  *v = *x > 1.0 ? 0.0 : BALL_RADIUS;
}


int struct_ball(cov_model *cov, cov_model **newmodel){
  ASSERT_NEWMODEL_NOT_NULL;
  if (hasMaxStableRole(cov)) {
    return addUnifModel(cov, BALL_RADIUS, newmodel);
  } else {
    ILLEGAL_ROLE;
  }

  return NOERROR;
}

int init_ball(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  assert(s != NULL);
  if (hasNoRole(cov)) return NOERROR;
  
  if (hasAnyShapeRole(cov)) {
    cov->mpp.maxheights[0] = 1.0;
    
    if (cov->mpp.moments >= 1) {
      cov->mpp.mM[1] = cov->mpp.mMplus[1] = VolumeBall(cov->tsdim, BALL_RADIUS);
      int i;
      for (i=2; i<=cov->mpp.moments; i++)  
	cov->mpp.mM[i] = cov->mpp.mMplus[i] = cov->mpp.mM[1];
    }
  }
  
  else ILLEGAL_ROLE;


  return NOERROR;
}


void do_ball(cov_model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 
  assert(s != NULL);
 
}



/// Poisson polygons


double meanVolPolygon(int dim, double beta) {
  double kd = VolumeBall(dim, 1.0),
    kdM1 = VolumeBall(dim-1, 1.0);
  return intpow(dim * kd / (kdM1 * beta), dim) / kd;
}

void Polygon(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) { 
  polygon_storage *ps = cov->Spolygon;
  assert(ps != NULL);
  *v = (double) isInside(ps->P, x);
}
 
void Inversepolygon(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v){
  polygon_storage *ps = cov->Spolygon;
  int d,
    dim = cov->tsdim;

  if (ps == NULL) {
    *v = RF_NA;
    return;
  }
  polygon *P = ps->P;
  if (P != NULL) {
    double max = 0.0;
    for (d=0; d<dim; d++) {
      double y = fabs(P->box0[d]);
      if (y > max) max = y;
      y = fabs(P->box1[d]);
      if (y > max) max = y;
    }
  } else {
    BUG;

    *v = pow(meanVolPolygon(dim, P0(POLYGON_BETA)) / VolumeBall(dim, 1.0), 
	     1.0/dim);
    // to do kann man noch mit factoren multiplizieren, siehe
    // strokorb/schlather, max
    // um unabhaengigkeit von der Dimension zu erlangen
  }
}

void InversepolygonNonstat(double VARIABLE_IS_NOT_USED *v, cov_model *cov,
			   double *min, double *max){
  polygon_storage *ps = cov->Spolygon;
  int d,
    dim = cov->tsdim;
  assert(ps != NULL);
  if (ps == NULL) {
    for (d=0; d<dim; d++) min[d] = max[d] = RF_NA;
    return;
  }
  polygon *P = ps->P;
  if (P != NULL) {
     for (d=0; d<dim; d++) {
      min[d] = P->box0[d];
      max[d] = P->box1[d];   
    }
  } else { // gibt aquivalenzradius eines "mittleres" Polygon zurueck

    BUG;

    double size = pow(meanVolPolygon(dim, P0(POLYGON_BETA)) / 
		      VolumeBall(dim, 1.0), 1.0/dim);    
    // to do kann man noch mit factoren multiplizieren, siehe
    // strokorb/schlather, max-stabile Modelle mit gleichen tcf
    for (d=0; d<dim; d++) {
      min[d] = -size;
      max[d] = size;
    }
  }
}

int check_polygon(cov_model *cov) {
  int err;
  assert(cov->tsdim <= MAXDIM_POLY);
  if (cov->tsdim != 2)
    SERR("random polygons only defined for 2 dimensions"); // to do
  assert(cov->tsdim == cov->xdimown);
  kdefault(cov, POLYGON_BETA, 1.0);
  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->deterministic = FALSE;
  return NOERROR;
}

void range_polygon(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[POLYGON_BETA] = 0;
  range->max[POLYGON_BETA] = RF_INF;
  range->pmin[POLYGON_BETA] = 1e-3;
  range->pmax[POLYGON_BETA] = 1e3;
  range->openmin[POLYGON_BETA] = true;
  range->openmax[POLYGON_BETA] = true;
}


int struct_polygon(cov_model VARIABLE_IS_NOT_USED *cov,
		   cov_model VARIABLE_IS_NOT_USED **newmodel){
  /*
  ASSERT_NEWMODEL_NOT_NULL;
  double beta = P0(POLYGON_BETA);
  if (hasAnyShapeRole(cov)) {
    double 
      dim = cov->tsdim,
      safety = P0(POLYGON_SAFETY), // to do : zhou approach !
      equiv_cube_length = pow(beta, -1.0/dim);
    return addUnifModel(cov,  // to do : zhou approach !
			0.5 * equiv_cube_length * safety,
			newmodel);
  } else {
    ILLEGAL_ROLE;
  }
  */
  BUG;

  return NOERROR;
}

polygon_storage *create_polygon() {
  polygon_storage *ps;
  if ((ps = (polygon_storage*) MALLOC(sizeof(polygon_storage))) == NULL)
    return NULL;
  if ((ps->P = (polygon*)  MALLOC(sizeof(polygon))) == NULL) {
    FREE(ps);
    return NULL;
  }
  polygon_NULL(ps);
  return ps;
}

int init_polygon(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int i, err,
    dim = cov->tsdim;
  double beta = P0(POLYGON_BETA),
    lambda = beta; // / VolumeBall(dim - 1, 1.0) * (dim * VolumeBall(dim, 1.0));
  assert(s != NULL);
  if (cov->Spolygon == NULL) {
    if ((cov->Spolygon = create_polygon()) ==NULL) return ERRORMEMORYALLOCATION;
  }
   
  // nicht nur aber auch zur Probe
  struct polygon_storage *ps = cov->Spolygon;
  assert(ps != NULL && ps->P != NULL);

  if (!false) {
    freePolygon(ps->P); 
    if ((err=rPoissonPolygon(ps->P, lambda, true)) != NOERROR)
      SERR1("poisson polygon cannot be simulated (error=%d)", err);
  } else {
    BUG; // valgrind memalloc loss ! + lange Laufzeiten ?!
    if ((err=rPoissonPolygon2(ps, lambda, true)) != NOERROR)
      SERR1("Poisson polygon cannot be simulated (error=%d)", err);
  }
 
  if (hasAnyShapeRole(cov)) {
    double c = meanVolPolygon(dim, beta);
    cov->mpp.maxheights[0] = 1.0; 
   for (i=1; i<=cov->mpp.moments; i++) cov->mpp.mM[i] = cov->mpp.mMplus[i] = c;	
  }  else ILLEGAL_ROLE;

  return NOERROR;
}


#define USER_TYPE 0
#define USER_DOM 1
#define USER_ISO 2
#define USER_VDIM 3
#define USER_BETA 4
#define USER_VARIAB 5
#define USER_FCTN 6
#define USER_FST 7
#define USER_SND 8
#define USER_ENV 9
#define USER_LAST USER_ENV
void evaluateUser(double *x, double *y, bool Time, cov_model *cov,
		  sexp_type *which, double *Res) {
  SEXP  res,
    env = PENV(USER_ENV)->sexp; 
  int i,
    vdim = cov->vdim[0] * cov->vdim[1],
    ncol = cov->ncol[USER_BETA],
    n = cov->xdimown;
  double *beta = P(USER_BETA);  
 
  assert(which != NULL);
  if (cov->nrow[USER_VARIAB] >= 2 && PINT(USER_VARIAB)[1] != -2) {
    if (Time) {
      addVariable( (char *) "T", x + (--n), 1, 1, env); 
    }
    switch (n) {
    case 3 : addVariable( (char *) "z", x+2, 1, 1, env); 		
    case 2 : addVariable( (char *) "y", x+1, 1, 1, env); 		
    case 1 : addVariable( (char *) "x", x+0, 1, 1, env); 		
      break;
    default:
      BUG;
    }
  } else {
    addVariable( (char *) "x", x, n, 1, env);	
    if (y != NULL) addVariable( (char *) "y", y, n, 1, env);
  }

  res = eval(which->sexp, env);
  if (beta == NULL) {
    for (i=0; i<vdim; i++) {
      Res[i] = REAL(res)[i]; 
    }
  } else {
    Ax(beta, REAL(res), vdim, ncol, Res);
  }
}


void kappaUser(int i, cov_model *cov, int *nr, int *nc){
  *nc = *nr = i < CovList[cov->nr].kappas ? 1 : -1;
  if (i == USER_VDIM) *nr = SIZE_NOT_DETERMINED;
  if (i == USER_VARIAB) *nr=SIZE_NOT_DETERMINED;
  if (i == USER_BETA) *nr=*nc=SIZE_NOT_DETERMINED;
}

void User(double *x, cov_model *cov, double *v){
  evaluateUser(x, NULL, Loc(cov)->Time, cov, PSEXP(USER_FCTN), v);
}
void UserNonStat(double *x, double *y, cov_model *cov, double *v){
  evaluateUser(x, y, false, cov, PSEXP(USER_FCTN), v);
}
void DUser(double *x, cov_model *cov, double *v){
  evaluateUser(x, NULL, Loc(cov)->Time, cov, PSEXP(USER_FST), v);
}
void DDUser(double *x, cov_model *cov, double *v){
  evaluateUser(x, NULL, Loc(cov)->Time, cov, PSEXP(USER_SND), v);
}


int checkUser(cov_model *cov){
  cov_fct *C = CovList + cov->nr;

  kdefault(cov, USER_DOM, XONLY);
  if (PisNULL(USER_ISO)) {
    Types type = (Types) P0INT(USER_TYPE);
    if (isVariogram(type)) kdefault(cov, USER_ISO, ISOTROPIC);
    else if (isProcess(type) || isShape(type)) 
      kdefault(cov, USER_ISO, CARTESIAN_COORD);
    else SERR1("type='%s' not allowed", TYPENAMES[type]);
  }
  int
    *dom = PINT(USER_DOM), 
    *iso = PINT(USER_ISO),
    *vdim =PINT(USER_VDIM);
  bool 
    Time,
    fctn = !PisNULL(USER_FCTN),
    fst = !PisNULL(USER_FST),
    snd = !PisNULL(USER_SND);
  int err, 
    *nrow = cov->nrow,
    *variab = PINT(USER_VARIAB), //codierung variab=1:x, 2:y 3:z, 4:T
    nvar = cov->nrow[USER_VARIAB],
    *pref = cov->pref;


  if (nvar < 1) SERR("variables not of required form ('x', 'y', 'z', 'T')");
  
  if (cov->domown != *dom) 
    SERR2("wrong domain (requ='%s'; provided='%s')", 
	  DOMAIN_NAMES[cov->domown], DOMAIN_NAMES[*dom]);
  if (cov->isoown != *iso)
    SERR2("wrong isotropy assumption (requ='%s'; provided='%s')",
	  ISONAMES[cov->isoown], ISONAMES[*iso]);

  if (PENV(USER_ENV) == NULL) BUG;

  if (!PisNULL(USER_BETA)) {
    if (vdim == NULL) kdefault(cov, USER_VDIM, nrow[USER_BETA]);
    else if (*vdim != nrow[USER_BETA]) 
      SERR4("number of columns of '%s' (=%d) does not equal '%s' (=%d)",
	    KNAME(USER_BETA), nrow[USER_BETA], KNAME(USER_VDIM), *vdim);
    cov->vdim[0] = nrow[USER_BETA];
    cov->vdim[1] = 1;
  } else {
    if (PisNULL(USER_VDIM)) kdefault(cov, USER_VDIM, 1);  
    cov->vdim[0] = P0INT(USER_VDIM);
    cov->vdim[1] = cov->nrow[USER_VDIM] == 1 ? 1 : PINT(USER_VDIM)[1];
  }
  if (cov->nrow[USER_VDIM] > 2) SERR1("'%s' must be a scalar or a vector of length 2", KNAME(USER_VDIM));

  if ((err = checkkappas(cov, false)) != NOERROR) return err;

  if (variab[0] != 1 && (variab[0] != 4 || nvar > 1)) SERR("'x' not given");
  if (nvar > 1) {
    variab[1] = abs(variab[1]); // integer!
    //                              it is set to its negativ value
    //                              below, when a kernel is defined
    if (variab[1] == 3) SERR("'z' given but not 'y'"); 
  } 
  Time = variab[nvar-1] == 4;

  if (((nvar >= 3 || variab[nvar-1] == 4)
       && (!ISNA(GLOBAL.coords.xyz_notation) && 
	   !GLOBAL.coords.xyz_notation))
      //  ||
      //  (nrow[USER_VARIAB] == 1 && !ISNA_INT(GLOBAL.coords.xyz_notation) 
      //  && GLOBAL.coords.xyz_notation)
      )
    SERR("mismatch of indicated xyz-notation");

  if (Time xor Loc(cov)->Time)
    SERR("If 'T' is given, it must be given both as coordinates and as variable 'T' in the function 'fctn'");

  if ((nvar > 2 || (nvar == 2 && variab[1] != 2)) && cov->domown == KERNEL)
    SERR1("'%s' not valid for anisotropic models", coords[COORDS_XYZNOTATION]);
  
  if (nvar == 2 && variab[1] == 2) {
    // sowohl nonstat also auch x, y Schreibweise moeglich
    if (ISNA(GLOBAL.coords.xyz_notation))
      SERR1("'%s' equals 'NA', but should be a logical value.", 
	   coords[COORDS_XYZNOTATION]);     
    if (cov->domown == KERNEL && GLOBAL.coords.xyz_notation==2) // RFcov
      SERR1("'%s' is not valid for anisotropic models", 
	   coords[COORDS_XYZNOTATION]);
  }
 
  if (nvar > 1) {
    if (cov->domown == XONLY) {
      if (cov->isoown == ISOTROPIC) {
	SERR("two many variables given for motion invariant function");
      } else if (cov->isoown == SPACEISOTROPIC && nvar != 2)
	SERR("number of variables does not match a space-isotropic model");
    } else {
      if (nvar == 2 && variab[1] == 2) variab[1] = -2;//i.e. non domain (kernel)
    }
  }

  if (!ISNA(GLOBAL.coords.xyz_notation)) {    
    if ((GLOBAL.coords.xyz_notation == 2 && cov->domown == KERNEL)
	||
	((nvar > 2 || (nvar == 2 && cov->domown==XONLY)) && variab[1] == -2)) {
      SERR("domain assumption, model and coordinates do not match.");
    }
  }

  if ((cov->xdimown == 4 && !Loc(cov)->Time) || cov->xdimown > 4)
    SERR("Notation using 'x', 'y', 'z' and 'T' allows only for 3 spatial dimensions and an additional time component.");


  if (fctn) {
    C->F_derivs = C->RS_derivs = 0;
    pref[Direct] = pref[Sequential] = pref[CircEmbed] = pref[Nothing] = 5;
    if (fst) {
      C->F_derivs = C->RS_derivs = 1;
      pref[TBM] = 5;
      if (snd) {
	C->F_derivs = C->RS_derivs = 2;
      }
    } else { // !fst
      if (snd) SERR2("'%s' given but not '%s'",
		     KNAME(USER_SND), KNAME(USER_FST));
    }
  } else { // !fctn
    if (fst || snd) 
      SERR3("'%s' or '%s' given but not '%s'", 
	    KNAME(USER_SND), KNAME(USER_FST), KNAME(USER_FCTN));
  }
   return NOERROR;
}


bool TypeUser(Types required, cov_model *cov, int VARIABLE_IS_NOT_USED depth) {
  int *type = PINT(USER_TYPE);
  if (type == NULL) return false;
  if (!isShape((Types) type[0])) return false;
  return TypeConsistency(required, (Types) type[0]);
}

void rangeUser(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[USER_TYPE] = TcfType;
  range->max[USER_TYPE] = TrendType;
  range->pmin[USER_TYPE] = TcfType;
  range->pmax[USER_TYPE] = TrendType;
  range->openmin[USER_TYPE] = false;
  range->openmax[USER_TYPE] = false;

  range->min[USER_DOM] = XONLY;
  range->max[USER_DOM] = KERNEL;
  range->pmin[USER_DOM] = XONLY;
  range->pmax[USER_DOM] = KERNEL;
  range->openmin[USER_DOM] = false;
  range->openmax[USER_DOM] = false;

  range->min[USER_ISO] = ISOTROPIC;
  range->max[USER_ISO] = CARTESIAN_COORD;
  range->pmin[USER_ISO] = ISOTROPIC;
  range->pmax[USER_ISO] = CARTESIAN_COORD;
  range->openmin[USER_ISO] = false;
  range->openmax[USER_ISO] = false;

  range->min[USER_VDIM] = 1;
  range->max[USER_VDIM] = INFDIM;
  range->pmin[USER_VDIM] = 1;
  range->pmax[USER_VDIM] = 10;
  range->openmin[USER_VDIM] = false;
  range->openmax[USER_VDIM] = true;

  range->min[USER_BETA] = RF_NEGINF;
  range->max[USER_BETA] = RF_INF;
  range->pmin[USER_BETA] = -1e5;
  range->pmax[USER_BETA] = 1e5;
  range->openmin[USER_BETA] = true;
  range->openmax[USER_BETA] = true;

  range->min[USER_VARIAB] = -2;
  range->max[USER_VARIAB] = 4;
  range->pmin[USER_VARIAB] = 1;
  range->pmax[USER_VARIAB] = 4;
  range->openmin[USER_VARIAB] = false;
  range->openmax[USER_VARIAB] = false;
}



// Sigma(x) = diag>0 + A'xx'A
void kappa_Angle(int i, cov_model *cov, int *nr, int *nc){
  int dim = cov->xdimown;
  *nr = i == ANGLE_DIAG ? dim : 1;
  *nc = i <= ANGLE_DIAG && (dim==3 || i != ANGLE_LATANGLE) &&  
    (dim==2 || i != ANGLE_RATIO)  ? 1 : -1;
}

void AngleMatrix(cov_model *cov, double *A) {
  double
    c = cos(P0(ANGLE_ANGLE)),
    s = sin(P0(ANGLE_ANGLE)),
    *diag = P(ANGLE_DIAG);
  int 
    dim = cov->xdimown ;
  assert(dim ==2 || dim == 3);
  
  if (dim == 2) {
    A[0] = c;
    A[1] = s;
    A[2] = -s;
    A[3] = c;
  } else {
    double 
      pc = cos(P0(ANGLE_LATANGLE)),
      ps = sin(P0(ANGLE_LATANGLE));
    /*
    c -s 0    pc 0 -ps   c*pc  -s  -c*ps
    s  c 0  *  0 1  0  = s*pc   c  -s*ps
    0  0 1    ps 0  pc     ps   0   pc 
    */
    A[0] = c * pc;
    A[1] = s * pc;
    A[2] = ps;
    A[3] = -s;
    A[4] = c;
    A[5] = 0.0;
    A[6] = -c * ps;
    A[7] = -s * ps;
    A[8] = pc;
  }

  if (diag!= NULL) {
    int k,d,j;
    for (k=d=0; d<dim; d++) {
      for (j=0; j<dim; j++) {
	A[k++] *= diag[j];
      }
    }
  } else {
    double ratio = 1.0 / P0(ANGLE_RATIO);
    A[1] *= ratio;
    A[3] *= ratio;
  }  
}



void Angle(double *x, cov_model *cov, double *v) { /// to do: extend to 3D!
  double A[9];
  int 
    dim = cov->xdimown ;
  AngleMatrix(cov, A);
  Ax(A, x, dim, dim, v); 
}


void invAngle(double *x, cov_model *cov, double *v) { /// to do: extend to 3D!
  double A[9],
    c = cos(P0(ANGLE_ANGLE)),
    s = sin(P0(ANGLE_ANGLE)),
    *diag = P(ANGLE_DIAG);
  int d, 
      dim = cov->xdimown ;
  
  bool inf = x[0] == RF_INF;
  for (d=1; d<dim; d++) inf &= x[d] == RF_INF;
  if (inf) {
    for (d=0; d<dim; v[d++] = RF_INF);
    return;
  }

  bool neginf = x[0] == RF_NEGINF;
  for (d=1; d<dim; d++) neginf &= x[d] == RF_NEGINF;
  if (neginf) {
    for (d=0; d<dim; v[d++] = RF_NEGINF);
    return;
  }
  
  if (dim == 2) {
    A[0] = c;
    A[1] = -s;
    A[2] = s;
    A[3] = c;
  } else {
    double pc = cos(P0(ANGLE_LATANGLE)),
      ps = sin(P0(ANGLE_LATANGLE));

    A[0] = c * pc;
    A[1] = -s;
    A[2] = -c * ps;
    A[3] = s * pc;
    A[4] = c;
    A[5] = -s * ps;
    A[6] = ps;
    A[7] = 0.0;
    A[8] = pc;
  }

  if (diag!= NULL) {
    int  j, k;
    for (k=d=0; d<dim; d++) {
      for (j=0; j<dim; j++) {
	A[k++] /= diag[d];
      }
    }
  } else {
    double ratio = P0(ANGLE_RATIO);
    A[2] *= ratio;
    A[3] *= ratio;
  }

  Ax(A, x, dim, dim, v);
}

int checkAngle(cov_model *cov){
  int dim = cov->xdimown;

  if (dim != 2 && dim != 3)
    SERR1("'%s' only works for 2 and 3 dimensions", NICK(cov));

  if (PisNULL(ANGLE_DIAG)) {
    if (PisNULL(ANGLE_RATIO)) {
      SERR2("either '%s' or '%s' must be given",
	    KNAME(ANGLE_RATIO), KNAME(ANGLE_DIAG))
    } else if (dim != 2) { 
      SERR1("'%s' may be given only if dim=2",  KNAME(ANGLE_RATIO))
    }
  } else {
    if (!PisNULL(ANGLE_RATIO)) 
      SERR2("'%s' and '%s' may not given at the same time",
	    KNAME(ANGLE_RATIO), KNAME(ANGLE_DIAG));
  }
  cov->vdim[0] = dim;
  cov->vdim[1] = 1;
  cov->mpp.maxheights[0] = RF_NA;
  cov->matrix_indep_of_x = true;
  return NOERROR;
}
 
void rangeAngle(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[ANGLE_ANGLE] = 0.0;
  range->max[ANGLE_ANGLE] = TWOPI;
  range->pmin[ANGLE_ANGLE] = 0.0;
  range->pmax[ANGLE_ANGLE] = TWOPI;
  range->openmin[ANGLE_ANGLE] = false;
  range->openmax[ANGLE_ANGLE] = true;

  range->min[ANGLE_LATANGLE] = 0.0;
  range->max[ANGLE_LATANGLE] = PI;
  range->pmin[ANGLE_LATANGLE] = 0.0;
  range->pmax[ANGLE_LATANGLE] = PI;
  range->openmin[ANGLE_LATANGLE] = false;
  range->openmax[ANGLE_LATANGLE] = true;

  range->min[ANGLE_RATIO] = 0;
  range->max[ANGLE_RATIO] = RF_INF;
  range->pmin[ANGLE_RATIO] = 1e-5;
  range->pmax[ANGLE_RATIO] = 1e5;
  range->openmin[ANGLE_RATIO] = false;
  range->openmax[ANGLE_RATIO] = true;

  range->min[ANGLE_DIAG] = 0;
  range->max[ANGLE_DIAG] = RF_INF;
  range->pmin[ANGLE_DIAG] = 1e-5;
  range->pmax[ANGLE_DIAG] = 1e5;
  range->openmin[ANGLE_DIAG] = false;
  range->openmax[ANGLE_DIAG] = true;
}
 


void idcoord(double *x, cov_model *cov, double *v){
  int d,
    vdim = cov->vdim[0];
  for (d=0; d<vdim; d++) v[d]=x[d];
}
int checkidcoord(cov_model *cov){
  if (cov->isoprev != cov->isoown) return ERRORWRONGISO;
  cov->vdim[0] = cov->xdimown;
  cov->vdim[1] = 1;
  return NOERROR;
}


// obsolete??!!
#define NULL_TYPE 0
void NullModel(double VARIABLE_IS_NOT_USED *x, 
	       cov_model VARIABLE_IS_NOT_USED *cov, 
	       double VARIABLE_IS_NOT_USED *v){ return; }
void logNullModel(double VARIABLE_IS_NOT_USED *x, 
		  cov_model VARIABLE_IS_NOT_USED *cov, 
		  double VARIABLE_IS_NOT_USED *v, 
		  int VARIABLE_IS_NOT_USED *Sign){ return; }
bool TypeNullModel(Types required, cov_model *cov, int VARIABLE_IS_NOT_USED depth) {
  Types type = (Types) P0INT(NULL_TYPE);
  return TypeConsistency(required, type); 
}
void rangeNullModel(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[NULL_TYPE] = TcfType;
  range->max[NULL_TYPE] = OtherType;
  range->pmin[NULL_TYPE] = TcfType;
  range->pmax[NULL_TYPE] = OtherType;
  range->openmin[NULL_TYPE] = false;
  range->openmax[NULL_TYPE] = false;
}


int checkMath(cov_model *cov){
  int i,
    variables = CovList[cov->nr].kappas,
    err = NOERROR;
  if (CovList[cov->nr].kappas > 2) kdefault(cov, 2, 1.0);
  if (cov->typus == TrendType && !isCoordinateSystem(cov->isoown))
    SERR("full coordinates need");

  //printf("checmath: %s %d\n", NAME(cov), variables);

 
  for (i=0; i<variables; i++) {
    cov_model *sub = cov->kappasub[i];
    //printf("   %d %d\n", i, sub == NULL);
    if (sub != NULL) {
      if (i >= 2) SERR("only numerics allowed");
      bool plus = CovList[sub->nr].cov == Mathplus ||
	CovList[sub->nr].check == checkplus ||
	CovList[sub->nr].cov == Mathminus
	;
       if ((err = CHECK(sub, cov->tsdim, cov->xdimown, 
			plus  ? cov->typus : ShapeType, 
			// auch falls cov = TrendType ist
			cov->domown, cov->isoown,
			1, ROLE_COV)) != NOERROR){
	return err;  
      }
      if  (sub->vdim[0] != 1 || sub->vdim[1] != 1)
	SERR("only scalar functions are allowed");
      setbackward(cov, sub); 
    } else if (PisNULL(i)) {
      if (i==0 || (CovList[cov->nr].cov!=Mathplus && 
		   CovList[cov->nr].cov!=Mathminus && 
		   CovList[cov->nr].cov!=Mathbind
		   )) SERR("not enough arguments given")
      else break;
    }
  }

  cov->pref[TrendEval] = 5;
  cov->pref[Direct] = 1;
  return NOERROR;
}


void rangeMath(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  int i,
    variables = CovList[cov->nr].kappas;
  
  cov->maxdim = cov->xdimown;
  assert(cov->maxdim > 0);
  for (i=0; i<variables; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = 1e5;
    range->pmax[i] = -1e5;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}

#define MATH_FACTOR 2
void Mathminus(double *x, cov_model *cov, double *v){
  MATH_DEFAULT;
  double f = P0(MATH_FACTOR); 
  *v = ( (PisNULL(1) && cov->kappasub[1]==NULL) ? -w[0] : w[0 ]- w[1]) * f;
  //printf("minus (%f - %f) * %f  = %f %d\n", w[0], w[1], f, *v, (PisNULL(1) && cov->kappasub[1]==NULL));
  // APMI(cov);
}

void Mathplus(double *x, cov_model *cov, double *v){
  MATH_DEFAULT;
   double f = P0(MATH_FACTOR); 
  *v = ( (PisNULL(1) && cov->kappasub[1]==NULL)? w[0] : w[0] + w[1]) * f;
}

void Mathdiv(double *x, cov_model *cov, double *v){
  MATH_DEFAULT;
   double f = P0(MATH_FACTOR); 
 *v =  (w[0] / w[1]) * f;

}

void Mathmult(double *x, cov_model *cov, double *v){
  MATH_DEFAULT;
   double f = P0(MATH_FACTOR); 
  *v = (w[0] * w[1]) * f;
}
 


void Mathbind(double *x, cov_model *cov, double *v){
  int vdim = cov->vdim[0];
  MATH_DEFAULT_0(vdim); 

  double f = P0(MATH_FACTOR); 
  for (i=0; i<vdim; i++) v[i] = w[i] * f; 
}

int check_bind(cov_model *cov) {  
  int err;
  int 
    variables = CovList[cov->nr].kappas - 1; // last is factor !!
  if ((err = checkMath(cov)) != NOERROR) return err;

  int i;
  for (i =0; i<variables; i++)
    // printf("%d %d %d null: %d %d\n", i, cov->nrow[i], cov->ncol[i], PisNULL(i), 	  cov->kappasub[i] == NULL );


  while (variables > 0 && cov->nrow[variables - 1] == 0 &&
	 cov->kappasub[variables - 1] == NULL) variables--;
  cov->vdim[0] = variables;
  cov->vdim[1] = 1;
  cov->ptwise_definite = pt_mismatch;
  return NOERROR;
}




#define C_C 0
void Mathc(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v) {
  *v = P0(C_C);
}

void rangec(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[C_C] = RF_NEGINF;
  range->max[C_C] = RF_INF;
  range->pmin[C_C] = 1e5;
  range->pmax[C_C] = -1e5;
  range->openmin[C_C] = true;
  range->openmax[C_C] = true;
}

int check_c(cov_model *cov){
  if (cov->kappasub[0] != NULL) SERR("only numerics allowed");
  cov->ptwise_definite =  
    P0(C_C) > 0 ? pt_posdef : P0(C_C) < 0 ? pt_negdef : pt_zero;
  return NOERROR;
}


 

#define PROJ_PROJ 0
#define PROJ_ISO 1
#define PROJ_FACTOR 2

void proj(double *x, cov_model *cov, double *v){
  *v = x[P0INT(PROJ_PROJ) - 1] * P0(PROJ_FACTOR);
}

int checkproj(cov_model *cov){
  int
    isoown = cov->isoown;
  kdefault(cov, PROJ_PROJ, 1);
  kdefault(cov, PROJ_FACTOR, 1.0);
  kdefault(cov, PROJ_ISO, CARTESIAN_COORD);
  int iso = P0INT(PROJ_ISO);
  
  if (isoown != iso && (iso != UNREDUCED || !isCoordinateSystem(isoown))) {
    SERR2("Offered system ('%s') does not match the required one ('%s')",
	  ISONAMES[isoown], ISONAMES[iso]);
  }
  return NOERROR;
}

void rangeproj(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[PROJ_PROJ] = 0;
  range->max[PROJ_PROJ] = cov->xdimown;
  range->pmin[PROJ_PROJ] = 1;
  range->pmax[PROJ_PROJ] = cov->xdimown;
  range->openmin[PROJ_PROJ] = false;
  range->openmax[PROJ_PROJ] = false;

  range->min[PROJ_FACTOR] = RF_NEGINF;
  range->max[PROJ_FACTOR] = RF_INF;
  range->pmin[PROJ_FACTOR] = 1e5;
  range->pmax[PROJ_FACTOR] = -1e5;
  range->openmin[PROJ_FACTOR] = true;
  range->openmax[PROJ_FACTOR] = true;

  range->min[PROJ_ISO] = ISOTROPIC;
  range->max[PROJ_ISO] = UNREDUCED;
  range->pmin[PROJ_ISO] = ISOTROPIC;
  range->pmax[PROJ_ISO] = UNREDUCED;
  range->openmin[PROJ_ISO] = false;
  range->openmax[PROJ_ISO] = false;
}



// fix (matrix); 
void fix(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v){
  location_type *loc = Loc(cov);
  listoftype *list= PLIST(FIX_M);
  int i,j,
    element = P0INT(FIX_ELMNT),
    *nrow = list->nrow; // anzahl listen elemente
  int vdim = cov->vdim[0],
    lx = nrow[element] / vdim, 
    ly = nrow[element] * lx;
  assert(cov->vdim[0] == cov->vdim[1]);
  double **p = list->p,
    *X = (p[element] + loc->i_row + loc->i_col * nrow[element]),
    *Y=X,
    *w = v;
  if (loc->i_row >= lx || loc->i_col >= lx) {
    PRINTF("size=%d current indices=(%d, %d)\n", lx, loc->i_row, loc->i_col);
    ERR("fix model: indices out of range");
  }
  if (element >= cov->nrow[FIX_M])
    ERR("fix model: list index out of range");
  for (i=0; i<vdim; i++, Y+=ly) {
     for (X=Y, j=0; j<vdim; j++, w++, X+=lx) {
       *w = *X;
     }
  }
}


void fix_nonstat(double *x, double VARIABLE_IS_NOT_USED *y, cov_model *cov, double *v){
    fix(x, cov, v);
}

void covmatrix_fix(cov_model *cov, double *v, int *nonzeros) {
  int j, nz=0,
    element = P0INT(FIX_ELMNT);
  listoftype *M= PLIST(FIX_M);
  double *p = M->p[element];
  int n = M->nrow[element] * M->ncol[element];   
  for (j=0; j<n; j++) nz = (v[j] = p[j]) != 0.0;
  *nonzeros = nz;
  loc_set(cov, M->nrow[element]);
}
char iscovmatrix_fix(cov_model VARIABLE_IS_NOT_USED *cov) {  return 2; }

int checkfix(cov_model *cov) {
  listoftype *list = PLIST(FIX_M); 
  long total;
  int info, err, i, vdim, totpts,
    m = cov->nrow[FIX_M],
    *q,
    *ncol = list->ncol,
    *nrow = list->nrow; // anzahl listen elemente
  double *dummy,
    **p = list->p;


  if (cov->q != NULL) {
    cov->vdim[0] = cov->vdim[1] = P0INT(FIX_VDIM);
    return ((int*) cov->q)[0]; 
  } else QALLOC(1);

  q = ((int*) cov->q);
  q[0] = NOERROR;

  if ((err = checkkappas(cov, false)) != NOERROR) return err;
  kdefault(cov, FIX_ELMNT, 0);
  kdefault(cov, FIX_VDIM, 1); // vdim !
  cov->vdim[0] = cov->vdim[1] = vdim = P0INT(FIX_VDIM);
  if (vdim > 1) return q[0] = ERRORVDIMNOTPROGRAMMEDYET;      
  // frag ist hier in welcher Reihenfolde die multivariaten und raeumlichen
  // korrelationen abfolgen, siehe vario** fuer die 2 Moeglichkeiten
   
  cov->ptwise_definite = pt_posdef;
  for (i=0; i<m; i++) {
    if (nrow[i] != ncol[i] || cov->nrow[i] == 0) {
      return q[0] = ERROR_MATRIX_SQUARE;      
    }
    totpts = nrow[i] / vdim;
    if (totpts * vdim != nrow[i]) {
	return q[0] = ERROR_MATRIX_VDIM;      
    }
    
    // check whether positive eigenvalue  
    total = nrow[i] * ncol[i] * sizeof(double); // NICHT weiter vorne
    dummy = (double*) MALLOC(total);
    MEMCOPY(dummy, p[i], total);
    
    if (false) {
      printf("in Primitive.cc; fix\n"); //
      int k, l, mm;
      for (mm=k=0; k<nrow[i]; k++) {
	for (l=0; l<ncol[i]; l++) {
	  printf("%4.0f ", dummy[mm++]); // /* check of false */
	}
	printf("\n");//
      }
      printf("i=%d %d %d\n", i, nrow[i], ncol[i]);//
    }
    
    F77_CALL(dpofa)(dummy, nrow + i, ncol + i, &info); // cholesky
    
    FREE(dummy);
    if (info != 0)  return q[0] = ERROR_MATRIX_POSDEF;
  }
    

  cov->matrix_indep_of_x = true;
  cov->mpp.maxheights[0] = RF_NA;
  err = checkkappas(cov);

  return err;
}


void rangefix(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[FIX_ELMNT] = 0;
  range->max[FIX_ELMNT] = MAXELEMENTS;
  range->pmin[FIX_ELMNT] = 0;
  range->pmax[FIX_ELMNT] = MAXELEMENTS;
  range->openmin[FIX_ELMNT] = false;
  range->openmax[FIX_ELMNT] = false;

  range->min[FIX_M] = RF_NEGINF;
  range->max[FIX_M] = RF_INF;
  range->pmin[FIX_M] = -1e10;
  range->pmax[FIX_M] = 1e10;
  range->openmin[FIX_M] = true;
  range->openmax[FIX_M] = true;
  
  range->min[FIX_VDIM] = 1;
  range->max[FIX_VDIM] = 9999;
  range->pmin[FIX_VDIM] = 1;
  range->pmax[FIX_VDIM] = 9999;
  range->openmin[FIX_VDIM] = false;
  range->openmax[FIX_VDIM] = false;
}
