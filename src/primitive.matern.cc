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



#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

#include "RF.h"
#include "operator.h"
#include "primitive.h"
#include "shape.h"
#include "AutoRandomFields.h"
#include "Processes.h"
#include "questions.h"


//#define LOG2 M_LN2

#define i11 0
#define i21 1
#define i22 2

#define epsilon 1e-15


double WhittleUpperNu[Nothing + 1] =
{80, 80, 80, // circulant
 80, RF_INF, // TBM
 80, 80,    // direct & sequ
 RF_NA, RF_NA, RF_NA,  // GMRF, ave, nugget
 RF_NA, RF_NA, RF_NA,  // exotic
 RF_INF   // Nothing
};


#define nstuetz 19L
#define stuetzorigin 11L
 double stuetz[]=
  {1.41062516176753e-14, 4.41556861847671e-12, 2.22633601732610e-06, 
   1.58113643548649e-03, 4.22181082102606e-02, 2.25024764696152e-01,
   5.70478218148777e-01, 1.03102017706644e+00, 1.57836638352906e+00,
   2.21866372852304e+00, 2.99573229151620e+00, 3.99852231863082e+00,
   5.36837527567695e+00, 7.30561120838150e+00, 1.00809957038601e+01,
   1.40580075785156e+01, 1.97332533513488e+01, 2.78005149402352e+01,
   3.92400265713477e+01};

/* mastein */
// see Hypermodel.cc


#define GET_NU_GEN(NU) (PisNULL(WM_NOTINV) || P0INT(WM_NOTINV) ? NU : 1.0 / NU)
#define GET_NU GET_NU_GEN(P0(WM_NU))


/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */

#define LOGGAMMA(nu) lgammafn(nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES)


#define LOW_MATERN 1e-20
double logWM(double x, double nu, double loggamma, double factor) {
  // check calling functions, like hyperbolic and gneiting if any changings !!
  //  printf("&& %10g\n", x);
  double v, y, 
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
		   scale = 1.0;
  if (factor!=0.0) scale = factor * SQRT(nuThres);
  double bk[MATERN_NU_THRES + 1L];
 
  if (x > LOW_MATERN) {
    //
    y = x  * scale;
    //    printf("$$\n");
    
      v = LOG2 + nuThres * LOG(0.5 * y) - loggamma + 
	LOG(bessel_k_ex(y, nuThres, 2.0, bk)) - y;
  } else v = 0.0;
    
  if (nu > MATERN_NU_THRES) { // factor!=0.0 && 
    //  printf("UU \n");
    double w, sign,
      g = MATERN_NU_THRES / nu;
    y = 0.5 * x * factor;
    logGauss(&y, NULL, &w, &sign);
    v = v * g + (1.0 - g) * w;
  }
  //  printf("^^ %10g\n", x);

  return v;
}


double Ext_WM(double x, double nu, double loggamma, double factor) {
  // check calling functions, like hyperbolic and gneiting if any changings !!
  return EXP(logWM(x, nu, loggamma, factor));
}

double Ext_DWM(double x, double nu, double loggamma, double factor) { 
  double   y, v,
    nuThres = nu <= MATERN_NU_THRES ? nu : MATERN_NU_THRES,
  		   scale = 1.0;
  if (factor!=0.0) scale = factor * SQRT(nuThres);
  double bk[MATERN_NU_THRES + 1L];
  
  if (x > LOW_MATERN) {  
    y = x * scale;  
    v = - 2.0 * EXP(nuThres * LOG(0.5 * y) - loggamma + 
		    LOG(bessel_k_ex(y, nuThres - 1.0, 2.0, bk)) - y);
  } else {
    v = (nuThres > 0.5) ? 0.0 : (nuThres < 0.5) ? INFTY : 1.253314137;
  }
  v *= scale;

  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = 0.5 * factor;
    y = x * scale;
    DGauss(&y, NULL, &w);
    w *= scale;
    v = v * g + (1.0 - g) * w;
  }
  return v;
}

double Ext_DDWM(double x, double nu, double gamma, double factor) { 
  double  y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES, // >
 		   scale = 1.0;
  if (factor!=0.0) scale = factor * SQRT(nuThres);
  double
		   scaleSq  = scale * scale,
		   bk[MATERN_NU_THRES + 1L];
  
  if (x > LOW_MATERN) {  
    y = x * scale;
    v = POW(0.5 * y , nuThres - 1.0) / gamma *
      (- bessel_k_ex(y, nuThres - 1.0, 1.0, bk) +
       y * bessel_k_ex(y, nuThres - 2.0, 1.0, bk));
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

double Ext_D3WM(double x, double nu, double gamma, double factor) { 
  double y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
 		   scale = 1.0;
  if (factor!=0.0) scale = factor * SQRT(nuThres);
  double scaleSq  = scale * scale;
  double bk[MATERN_NU_THRES + 1L];
  
  if (x > LOW_MATERN) {
    y = x * scale;
    v = POW(0.5 * y , nuThres - 1.0) / gamma *
      ( 3.0 * bessel_k_ex(y, nuThres - 2.0, 1.0, bk) 
	-y * bessel_k_ex(y, nuThres - 3.0, 1.0, bk)); 
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

double Ext_D4WM(double x,  double nu, double gamma, double factor) {
  
  double y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
		   scale = 1.0;
  if (factor!=0.0) scale = factor * SQRT(nuThres);
  double scaleSq  = scale * scale;
  double bk[MATERN_NU_THRES + 1L];

 //  printf("x=%10g nu=%10g\n", x, nuThres);
  
  if (x > LOW_MATERN) {
    y = x * scale;
    v = 0.25 * POW(0.5 * y , nuThres - 3.0) / gamma *
      (+ 6.0 * (nuThres - 3.0 - y * y) *
       bessel_k_ex(y, nuThres - 3.0, 1.0, bk)
       + y * (3.0  + y * y) * bessel_k_ex(y, nuThres-4.0, 1.0, bk));
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



  
double logNonStWM(double *x, double *y, model *cov, double factor){
  model *nu = cov->kappasub[WM_NU];
  int 
     dim = OWNLOGDIM(0);
  double nux, nuy, v,
    norm=0.0;
  for (int d=0; d<dim; d++) {
    double delta = x[d] - y[d];
    norm += delta * delta;
  }
  norm = SQRT(norm);
  
  if (nu == NULL) {
    nux = nuy = GET_NU;
  } else {
    FCTN(x, nu, &nux);
    FCTN(y, nu, &nuy);
    if (nux <= 0.0 || nuy <= 0.0)
      ERR1("'%.50s' is not a positive function", KNAME(WM_NU));
    nux = GET_NU_GEN(nux);
    nuy = GET_NU_GEN(nuy);
  }
    v = Ext_logWM(norm, nux, nuy, factor);
  assert(!ISNAN(v));
  return v;
}


double ScaleWM(double nu){
  // it is working reasonably well if nu is in [0.001,100]
  // happy to get any better suggestion!!
  return interpolate(LOG(nu) * INVLOG2, stuetz, nstuetz, stuetzorigin, 
		     1.5, 5);
}


 
bool setWM(model *cov) {
  model *nusub = cov->kappasub[WM_NU];
  //  PMI(cov);  printf("X %d\n", LASTi((PREV)[0]));  printf("%d\n", PREV_INITIALISED);   printf("Y\n");
  isotropy_type iso = CONDPREVISO(0);
  if (!isFixed(iso)) return false;
  if (nusub != NULL && !isRandom(nusub)) {    
    set_dom(OWN, 0, KERNEL);
    set_iso(OWN, 0, isAnySpherical(iso) ? SPHERICAL_SYMMETRIC : SYMMETRIC);
  } else {
    set_dom(OWN, 0, XONLY);
    set_iso(OWN, 0, isAnySpherical(iso) ? SPHERICAL_ISOTROPIC : ISOTROPIC);
 }
  return true;
}

bool allowedDWM(model *cov) {
  model *nusub = cov->kappasub[WM_NU];
  bool *D = cov->allowedD;
  for (int i=FIRST_DOMAIN; i<LAST_DOMAINUSER; D[i++]=false);
  D[nusub != NULL && !isRandom(nusub) ? KERNEL : XONLY] = true;
  return false;
}

bool allowedIWM(model *cov) {
  model *nusub = cov->kappasub[WM_NU];
  bool *I = cov->allowedI;
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  if (nusub != NULL && !isRandom(nusub)) {
    I[SYMMETRIC] = I[SPHERICAL_SYMMETRIC] = true;
  } else {
    I[ISOTROPIC] = I[SPHERICAL_ISOTROPIC] = true;
  }
  return false;
}

int initWM(model *cov, gen_storage VARIABLE_IS_NOT_USED  *s) {
  if (!PisNULL(WM_NU)) {
    double nu = GET_NU;
    if (!ISNA(nu)) {
      double nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES;
      cov->q[WM_LOGGAMMA] = lgammafn(nuThres);
      cov->q[WM_GAMMA] = gammafn(nuThres);      
    }
  } 
  RETURN_NOERROR;
}


int checkWM(model *cov) { 
  model *nusub = cov->kappasub[WM_NU];
#define spectrallimit 0.17
#define spectralbest 0.4
  int i, err,
    dim = OWNLOGDIM(0);
  bool isna_nu;
  
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);  

  set_logdim(OWN, 0, GATTERLOGDIM(0));
  if (nusub != NULL && !isRandom(nusub)) {
    //PMI(cov)
    if (!isKernel(OWN) || !equalsSymmetric(OWNISO(0)))
      SERR2("kernel needed. Got %.50s, %.50s.",
	    DOMAIN_NAMES[OWNDOM(0)], ISO_NAMES[OWNISO(0)]);
    ASSERT_CARTESIAN;
    set_xdim(OWN, 0, GATTERXDIM(0));
    if ((err = CHECK(nusub, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		     SCALAR, cov->frame // VariogramType changed 20.7.14 wg spectral
		     )) != NOERROR) 
      RETURN_ERR(err);
     
    if (LOGDIM(SYSOF(nusub), 0) != dim) RETURN_ERR(ERRORWRONGDIM);
    cov->monotone = NORMAL_MIXTURE;
    RETURN_NOERROR;    
  }


  if (OWNDOM(0) != XONLY || !isAnyIsotropic(OWNISO(0))) 
    // RETURN_ERR(ERRORFAILED);
    SERR2("isotropic function needed. Got %.50s, %.50s.",
	  DOMAIN_NAMES[OWNDOM(0)], ISO_NAMES[OWNISO(0)]);

  if (PisNULL(WM_NU)) QERRC(0, "parameter unset"); 
  double nu = GET_NU;
  isna_nu = ISNAN(nu);

  if (cov->q == NULL) {
    QALLOC(2);
    initWM(cov, NULL);
  }
  
  for (i=0; i<= Nothing; i++) cov->pref[i] *= isna_nu || nu < WhittleUpperNu[i];
  if (nu<spectralbest) {
    cov->pref[SpectralTBM] = (nu < spectrallimit) ? PREF_NONE : 3;
  } 
  if (dim > 2) 
   cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  if (nu > 2.5)  cov->pref[CircEmbed] = 2;

  cov->full_derivs = isna_nu ? 0
    : nu == (int) nu ? 2L * (int) nu - 2L : 2L * (int) nu;
  cov->monotone = nu <= 0.5 ? COMPLETELY_MON : NORMAL_MIXTURE;

  set_xdim(OWN, 0, 1);
  RETURN_NOERROR;
}

void rangeWM(model *cov, range_type *range){
  bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));
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
 
  booleanRange(WM_NOTINV);

}


void coinitWM(model *cov, localinfotype *li) {
  // cutoff_A
  double thres[2] = {0.25, 0.5},
    nu = GET_NU;
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

void ieinitWM(model *cov, localinfotype *li) {
  double nu = GET_NU;
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

double densityWM(double *x, model *cov, double factor) {
  double x2,
    nu = GET_NU,
    powfactor = 1.0;
  int d,
    dim =  PREVLOGDIM(0);
  if (nu > 50)
    warning("nu>50 in density of matern class numerically instable. The results cannot be trusted.");
  if (factor == 0.0) factor = 1.0; else factor *= SQRT(nu);
  x2 = x[0] * x[0];
  for (d=1; d<dim; d++) {
    x2 += x[d] * x[d];
    powfactor *= factor;
  }
  x2 /= factor * factor;
  x2 += 1.0;
  
  return powfactor * EXP(lgammafn(nu + 0.5 * (double) dim)
			 - lgammafn(nu)
			 - (double) dim * M_LN_SQRT_PI
			 - (nu + 0.5 * (double) dim) * LOG(x2));
}


void Matern(double *x, model *cov, double *v) {
  *v = Ext_WM(*x, GET_NU, cov->q[WM_LOGGAMMA], SQRT2);
}

void logMatern(double *x, model *cov, double *v, double *Sign) { 
  double nu = GET_NU;
  *v = logWM(*x, nu, cov->q[WM_LOGGAMMA], SQRT2);
  *Sign = 1.0;
}

void NonStMatern(double *x, double *y, model *cov, double *v){ 
  *v = EXP(logNonStWM(x, y, cov, SQRT2));
}

void logNonStMatern(double *x, double *y, model *cov, double *v, 
		     double *Sign){ 
  *Sign = 1.0;
  *v = logNonStWM(x, y, cov, SQRT2);
}

void DMatern(double *x, model *cov, double *v) {
  *v =Ext_DWM(*x, GET_NU, cov->q[WM_LOGGAMMA], SQRT2);
} 

void DDMatern(double *x, model *cov, double *v) {
  *v=Ext_DDWM(*x, GET_NU, cov->q[WM_GAMMA], SQRT2);
} 

void D3Matern(double *x, model *cov, double *v) {
  *v=Ext_D3WM(*x, GET_NU,  cov->q[WM_GAMMA], SQRT2);
} 

void D4Matern(double *x, model *cov, double *v) {
  *v=Ext_D4WM(*x, GET_NU, cov->q[WM_GAMMA], SQRT2);
} 

void InverseMatern(double *x, model *cov, double *v) {
  double
    nu = GET_NU;
  *v = RF_NA;
  if (*x == 0.05) *v = SQRT2 * SQRT(nu) /  ScaleWM(nu);
}

int checkMatern(model *cov) { 
  kdefault(cov, WM_NOTINV, true);
  return checkWM(cov);
}

Types TypeWM(Types required, model *cov, isotropy_type requ_iso){
  //  printf("requ=%.50s, %.50s\n", TYPE_NAMES[required], ISO_NAMES[requ_iso]);
  //  PMI(cov->calling);
  model *nusub = cov->kappasub[WM_NU];
  if (isCartesian(requ_iso)) {
    if (nusub != NULL) {
      if ((equalsXonly(OWNDOM(0)) && !isRandom(nusub)) ||
	  !isSymmetric(requ_iso))
	return BadType;
      return TypeConsistency(required, PosDefType );   
    }
    double nu = GET_NU;
    return TypeConsistency(required,
			   ISNAN(nu) || nu <= 0.5 ? TcfType : PosDefType );
  } else if (isSpherical(requ_iso)) {
    if (!isSphericalSymmetric(requ_iso) || nusub != NULL) return BadType;
    return TypeConsistency(required, PosDefType );
  } else if (isEarth(requ_iso)) {
    // assert(requ_iso == PREVISO(0));
    //    printf("TC0 %.50s %d %d %.50s\n",  TYPE_NAMES[required], !isEarthSymmetric(requ_iso), nusub != NULL, ISO_NAMES[requ_iso]);
    if (!isEarthSymmetric(requ_iso) || nusub != NULL) return BadType;
    return TypeConsistency(required, PosDefType );
  }
  return BadType;
}

double densityMatern(double *x, model *cov) {
  return densityWM(x, cov, SQRT2);
}

int initMatern(model *cov, gen_storage *s) {
  if (HAS_SPECTRAL_FRAME(cov)) {
    spec_properties *cs = &(s->spec);
    if (OWNLOGDIM(0) <= 2) RETURN_NOERROR;
    cs->density = densityMatern;
    return search_metropolis(cov, s);
  }

  else if (hasRandomFrame(cov)) { RETURN_NOERROR; }
  
  else ILLEGAL_FRAME;

}

void spectralMatern(model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  /* see Yaglom ! */
  int dim = PREVLOGDIM(0);
  if (dim <= 2) {
    double nu = GET_NU;
    E12(s, dim, 
	SQRT( 2.0 * nu * (POW(1.0 - UNIFORM_RANDOM, -1.0 / nu) - 1.0) ), e);
  } else {
    metropolis(cov, S, e);
  }
}


void Whittle(double *x, model *cov, double *v) {
  //  printf("!! %10g\n", *x);
  double nu = GET_NU;
  assert(cov->q != NULL);
  *v = Ext_WM(*x, nu, cov->q[WM_LOGGAMMA], 0.0);
  assert(!ISNAN(*v));
  //  printf("<< %10g\n", *x);
}


void logWhittle(double *x, model *cov, double *v, double *Sign) {
  double nu = GET_NU;
  *v = logWM(*x, nu, cov->q[WM_LOGGAMMA], 0.0);
  assert(!ISNAN(*v));
  *Sign = 1.0;
}


void NonStWhittle(double *x, double *y, model *cov, double *v){ 
  *v = EXP(logNonStWM(x, y, cov, 0.0));
}

void logNonStWhittle(double *x, double *y, model *cov, double *v, 
		     double *Sign){ 
  *Sign = 1.0;
  *v = logNonStWM(x, y, cov, 0.0);
}

void DWhittle(double *x, model *cov, double *v) {
  *v =Ext_DWM(*x, GET_NU, cov->q[WM_LOGGAMMA], 0.0);
}

void DDWhittle(double *x, model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=Ext_DDWM(*x, GET_NU, cov->q[WM_GAMMA], 0.0);
}


void D3Whittle(double *x, model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=Ext_D3WM(*x, GET_NU,  cov->q[WM_GAMMA], PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}

void D4Whittle(double *x, model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{
  *v=Ext_D4WM(*x, GET_NU, cov->q[WM_GAMMA], PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}

void InverseWhittle(double *x, model *cov, double *v){
  double nu = GET_NU;
  *v = (*x == 0.05) ? 1.0 / ScaleWM(nu) : RF_NA;
}

void TBM2Whittle(double *x, model *cov, double *v) 
{
  double nu = GET_NU;
  assert(PisNULL(WM_NOTINV));
  if (nu == 0.5) TBM2exponential(x, cov, v);
  else BUG;
}


double densityWhittle(double *x, model *cov) {
  return densityWM(x, cov, PisNULL(WM_NOTINV) ? 0.0 : SQRT2);
}
int initWhittle(model *cov, gen_storage *s) {
  if (HAS_SPECTRAL_FRAME(cov)) {
    if (PisNULL(WM_NU)) {
      spec_properties *cs = &(s->spec);
      if (OWNLOGDIM(0) <= 2) RETURN_NOERROR;
      cs->density = densityWhittle; 
      int err=search_metropolis(cov, s);
      RETURN_ERR(err);
    } else return initMatern(cov, s);
  }

  else if (hasRandomFrame(cov)) { RETURN_NOERROR;}

  else ILLEGAL_FRAME;

}

void spectralWhittle(model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  /* see Yaglom ! */
  if (PisNULL(WM_NOTINV)) {
    int dim = PREVLOGDIM(0);
    if (dim <= 2) {
      double nu = P0(WM_NU);
      E12(s, dim, SQRT(POW(1.0 - UNIFORM_RANDOM, -1.0 / nu) - 1.0), e);
    } else {
      metropolis(cov, S, e);
    }
  } else spectralMatern(cov, S, e);
}


void DrawMixWM(model VARIABLE_IS_NOT_USED *cov, double *random) { // inv scale
  // V ~ F in PSgen
  *random = -0.25 / LOG(UNIFORM_RANDOM);
}

double LogMixDensW(double VARIABLE_IS_NOT_USED *x, double logV, model *cov) {
  double
    nu=GET_NU;
  return PisNULL(WM_NOTINV)
    ? (M_LN2  + 0.5 * logV) * (1.0 - nu) - 0.5 *lgammafn(nu)
    // - 0.25 /  cov->mpp[DRAWMIX_V]  - 2.0 * (LOG2 + logV) )
    : RF_NA;  /* !! */
}




// using nu^(-1-nu+a)/2 for g and v^-a e^{-1/4v} as density instead of frechet
// the bound 1/3 can be dropped
// s tatic double eM025 = EXP(-0.25);
//void DrawMixNonStWM(model *cov, double *random) { // inv scale
//  // V ~ F in stp
//  model *nu = cov->sub[WM_NU];  
//  double minnu;
//  double alpha;
//
//  if (nu == NULL) {
//    minnu = P(WM_NU][0];
//  } else {
//    double minmax[2];
//    DefList[nu->nr].minmaxeigenvalue(nu, minmax);
//    minnu = minmax[0];
//  }
//  alpha = 1.0 + 0.5 /* 0< . < 1*/ * (3.0 * minnu - 0.5 * cov->tsdim);
//  if (alpha > 2.0) alpha = 2.0; // original choice
//  if (alpha <= 1.0) ERR("minimal nu too low or dimension too high");
//
// ERR("logmixdensnonstwm not programmed yet");
 /* 
  double beta = GLOBAL.mpp.beta,
    p = GLOBAL.mpp.p,
    logU;
  if (UNIFORM_RANDOM < p){
    cov_a->WMalpha = beta;
    logU =  LOG(UNIFORM_RANDOM * eM025);
    cov_a->WMfactor = -0.5 * LOG(0.25 * p * (beta - 1.0)) + 0.25;
  } else {
    cov_a->WMalpha = alpha;
    logU = LOG(eM025 + UNIFORM_RANDOM * (1.0 - eM025));
    cov_a->WMfactor = -0.5 * LOG(0.25 * (1.0 - p) * (alpha - 1.0));
  } 
  
  logmix!!

  *random = LOG(-0.25 / logU) / (cov_a->WMalpha - 1.0); //=invscale
  */
//}




//double LogMixDensNonStWM(double *x, double logV, model *cov) {
//  // g(v,x) in stp
//  double z = 0.0;
//  ERR("logmixdensnonstwm not programmed yet");
  // wmfactor ist kompletter unsinn; die 2 Teildichten muessen addiert werden
  /*
  model *calling = cov->calling,
    *Nu = cov->sub[0];
  if (calling == NULL) BUG;
   double nu,
    alpha = cov_a->WMalpha,
    logV = cov_a->logV,
    V = cov_a->V;
  
  if (Nu == NULL) 
    nu = P(WM_NU][0];
  else 
     FCTN(x, Nu, &nu);


   z = - nu  * M_LN2 // in g0  // eine 2 kuerzt sich raus
    + 0.5 * ((1.0 - nu) // in g0
    + alpha // lambda
    - 2.0 //fre*
    ) * logV
    - 0.5 * lgammafn(nu)  // in g0
    + cov_a->WMfactor // lambda
    - 0.125 / V   // g: Frechet
    + 0.125 * POW(V, 1.0 - alpha); // lambda: frechet

  if (!(z < 7.0)) {
    s tatic double storage = 0.0; 
    if (gen_storage != logV) {
      if (PL >= PL_DETAILS) 
	PRINTF("alpha=%10g, is=%10g, cnst=%10g logi=%10g lgam=%10g loga=%10g invlogs=%10g pow=%10g z=%10g\n",
	       alpha,V,
	       (1.0 - nu) * M_LN2 
	       , + ((1.0 - nu) * 0.5 + alpha - 2.0) * logV
	       ,- 0.5 * lgammafn(nu) 
	       , -cov_a->WMfactor
	       ,- 0.25 / V 
	       , + 0.25 * POW(V, - alpha)
	       , z);
      storage = logV;
    }
    //assert(z < 10.0);
  }
*/
//  return z;
//				      
//}







void Whittle2(double *x, model *cov, double *v) {
  *v =Ext_WM(*x, GET_NU,  cov->q[WM_LOGGAMMA], 0.0);
  assert(!ISNAN(*v));
}


void logWhittle2(double *x, model *cov, double *v, double *Sign) {
  double nu = GET_NU;
  *v = logWM(*x, nu, cov->q[WM_LOGGAMMA], 0.0);
  assert(!ISNAN(*v));
  *Sign = 1.0;
}

void DWhittle2(double *x, model *cov, double *v) {
  *v =Ext_DWM(*x, GET_NU, cov->q[WM_LOGGAMMA], 0.0);
}

void DDWhittle2(double *x, model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=Ext_DDWM(*x, GET_NU, cov->q[WM_GAMMA], 0.0);
}


void D3Whittle2(double *x, model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=Ext_D3WM(*x, GET_NU, cov->q[WM_GAMMA], 0.0);
}

void D4Whittle2(double *x, model *cov, double *v)
// check calling functions, like hyperbolic and gneiting if any changings !!
{ 
  *v=Ext_D4WM(*x, GET_NU, cov->q[WM_GAMMA], 0.0);
}

void InverseWhittle2(double *x, model *cov, double *v){
  *v = (*x == 0.05) ? 1.0 / ScaleWM(GET_NU) : RF_NA;
}


/* Whittle-Matern or Whittle or Besset */ 
/// only lower parts of the matrices, including the diagonals are used when estimating !!

/* Whittle-Matern or Whittle or Besset */ 
void biWM2basic(model *cov, 
		double *a, double *lg,
		double *aorig, double *nunew) {
  // !! nudiag, nured has priority over nu
  // !! cdiag, rhored has priority over c

  double factor, beta, gamma, tsq, t1sq, t2sq, inf, infQ, discr,
    alpha , a2[3],
    dim = (double) OWNLOGDIM(0), 
    d2 = dim * 0.5, 
    *nudiag = P(BInudiag),
    nured = P0(BInured),
    *nu = P(BInu),
    *c = P(BIc),
    *s = P(BIs),
    *cdiag = P(BIcdiag),
    rho = P0(BIrhored);
  int i;
  bool notinvnu = PisNULL(BInotinvnu) || (bool) P0INT(BInotinvnu);

  nu[i11] = nudiag[0];
  nu[i22] = nudiag[1];
  nu[i21] = 0.5 * (nu[i11] + nu[i22]) * nured;
  
  for (i=0; i<3; i++) {
    aorig[i] = 1.0 / s[i];
  } 

  if (PisNULL(BInotinvnu)) {
    for (i=0; i<3; i++) {
      a[i] = aorig[i];
      nunew[i] = nu[i];
    }
  } else {
    if (!notinvnu) for (i=0; i<3; i++) nu[i] = 1.0 / nu[i];
    for (i=0; i<3; i++) {
      nunew[i] = nu[i] < MATERN_NU_THRES ? nu[i] : MATERN_NU_THRES;
      a[i] = aorig[i] * SQRT(2.0 * nunew[i]);
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
  
  factor = EXP(lgammafn(nunew[i11] + d2) - lg[i11]
	       + lgammafn(nunew[i22] + d2) - lg[i22]
		   + 2.0 * (lg[i21]  - lgammafn(nunew[i21] + d2)
			    +nunew[i11] * LOG(a[i11]) + nunew[i22] * LOG(a[i22])
			    - 2.0 * nunew[i21] * LOG(a[i21]))
	       );
  
  
  // alpha u^2 + beta u + gamma = 0
  gamma = (2.0 * nunew[i21] + dim) * a2[i11] * a2[i22] 
    - (nunew[i22] + d2) * a2[i11] * a2[i21]
    - (nunew[i11] + d2) * a2[i22] * a2[i21];
      
  beta = (2.0 * nunew[i21] - nunew[i22] + d2) * a2[i11]
      + (2.0 * nunew[i21] - nunew[i11] + d2) * a2[i22]
      - (nunew[i11] + nunew[i22] + dim) * a2[i21];
  
  if (nured == 1.0) { // lin. eqn only
    t2sq = (beta == 0.0) ? 0.0 : -gamma / beta;
    if (t2sq < 0.0) t2sq = 0.0;
    t1sq =  t2sq;
  } else { // quadratic eqn.
    discr = beta * beta - 4.0 * alpha * gamma;
    if (discr < 0.0) {
      t1sq = t2sq = 0.0;
    } else {
      discr = SQRT(discr);
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
    
    infQ = POW(a2[i21] + tsq, 2.0 * nunew[i21] + dim) /
      (POW(a2[i11] + tsq, nunew[i11] + d2) 
       * POW(a2[i22] + tsq, nunew[i22] + d2));

    if (infQ < inf) inf = infQ;
  }

  c[i11] = cdiag[0];
  c[i22] = cdiag[1];
  c[i21] = rho * SQRT(factor * inf * c[i11] *  c[i22]);
 
}


int initbiWM2(model *cov, gen_storage *stor) {
  double a[3], lg[3], aorig[3], nunew[3], 
    *c = P(BIc),
    *cdiag = P(BIcdiag),
    *nu = P(BInu),
    *nudiag = P(BInudiag);
  bool check = stor->check;
  bool notinvnu = PisNULL(BInotinvnu) || (bool) P0INT(BInotinvnu);
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  if (S->nudiag_given) {
    kdefault(cov, BInured, 1.0);
    if (check && !PisNULL(BInu)) {
      if (cov->nrow[BInu] != 3 || cov->ncol[BInu] != 1)
	QERRC(BInu, "must be a 3 x 1 vector");
      if (FABS(nu[i11] - nudiag[0]) > nu[i11] * epsilon || 
	  FABS(nu[i22] - nudiag[1]) > nu[i22] * epsilon ||
	  FABS(nu[i21] - 0.5 * (nudiag[i11] + nudiag[1]) * P0(BInured))
	  > nu[i21] * epsilon) {
	QERRC2(BInu, "does not match '%.50s' and '%.50s'.", BInudiag, BInured);
      }
      if ((bool) ISNA(nu[i11]) + (bool) ISNA(nu[i22]) + (bool)ISNA(nu[i21])
	  != (bool) ISNA(nudiag[0]) + (bool) ISNA(nudiag[1]) +
	  (bool) ISNA(P0(BInured)) ){
	QERRC2(BInu, "does not have matching NAs with '%.50s' and '%.50s'.",
	       BInudiag, BInured);
      }
    } else {
      if (PisNULL(BInu)) PALLOC(BInu, 3, 1);
      nu = P(BInu); 
      nu[i11] = nudiag[0];
      nu[i22] = nudiag[1];
      nu[i21] = 0.5 * (nu[i11] + nu[i22]) * P0(BInured);
    }
  } else {
    if (check && !PisNULL(BInured)) {
      QERRC1(BInured, "may not be set if '%.50s' is given", BInu);
    }
    if (PisNULL(BInu)) 
      QERRC2(BInu, "either '%.50s' or '%.50s' must be set", BInu, BInured);
    if (PisNULL(BInudiag)) PALLOC(BInudiag, 2, 1);
    nudiag = P(BInudiag);
    nudiag[0] = nu[i11]; 
    nudiag[1] = nu[i22];
    if (PisNULL(BInured)) PALLOC(BInured, 1, 1);
    P(BInured)[0] =  nu[i21] / (0.5 * (nudiag[i11] + nudiag[1]));
  }

  if (!notinvnu) 
    for (int i=0; i<3; i++) nu[i] = 1.0 / nu[i];
 
  cov->full_derivs = cov->rese_derivs = 1; // kann auf 2 erhoeht werden, falls programmiert 
  for (int i=0; i<3; i++) {
    int derivs = nu[i] - 1.0 > MAXINT ? MAXINT : (int) (nu[i] - 1.0);
    if (cov->full_derivs < derivs) cov->full_derivs = derivs;
  }
  
  if (PisNULL(BIs)) {
    PALLOC(BIs, 3, 1);
    double *s = P(BIs);
    for (int i=0; i<3; s[i++] = 1.0);
  }
  
 
  if  (S->cdiag_given) {
    if (PisNULL(BIrhored))
      QERRC2(BIrhored, "'%.50s' and '%.50s' must be set", BIcdiag, BIrhored);
    if (check && !PisNULL(BIc)) {
      if (cov->nrow[BIc] != 3 || cov->ncol[BIc] != 1)
	QERRC(BIc, "must be a 3 x 1 vector");
      if (FABS(c[i11] - cdiag[0]) > c[i11] * epsilon || 
	  FABS(c[i22] - cdiag[1]) > c[i22] * epsilon ) {
	QERRC1(BIc, "does not match '%.50s'.", BIcdiag);
      }
      double savec12 = c[i21];
      biWM2basic(cov, a, lg, aorig, nunew);
      if (FABS(c[i21] - savec12) > FABS(c[i21]) * epsilon)
 	QERRC1(BIc, "does not match '%.50s'.", BIrhored);
    } else {
      if (PisNULL(BIc)) PALLOC(BIc, 3, 1);
      c = P(BIc);
      biWM2basic(cov, a, lg, aorig, nunew);
   }
  } else {
    if (check && !PisNULL(BIrhored))  {
      QERRC1(BIrhored, "may not be set if '%.50s' is not given", BIcdiag);
    }
    if (PisNULL(BIc))
      QERRC2(BIc, "either '%.50s' or '%.50s' must be set", BIc, BIcdiag);
    if (!ISNAN(c[i11]) && !ISNAN(c[i22]) && (c[i11] < 0.0 || c[i22] < 0.0))
      QERRC2(BIc, "variance parameter %.50s[0], %.50s[2] must be non-negative.",
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
  for (int i=0; i<3; i++) {
    double nuThres = nu[i] < MATERN_NU_THRES ? nu[i] : MATERN_NU_THRES;
     cov->q[WM_LOGGAMMA + 2 * i] =  lgammafn(nuThres);
    cov->q[WM_GAMMA + 2 * i] = gammafn(nuThres);
  }
  

  cov->initialised = true;

  RETURN_NOERROR;
}


void kappa_biWM(int i, model *cov, int *nr, int *nc){
  *nc = *nr = i < DefList[COVNR].kappas ? 1 : -1;
  if (i==BInudiag || i==BIcdiag) *nr = 2; else
    if (i==BInu || i==BIs || i==BIc) *nr = 3;
}

void biWM2(double *x, model *cov, double *v) {
  int i;
  double 
    *c = P(BIc),
    *nu = P(BInu),
    xx = *x;
  biwm_storage *S = cov->Sbiwm;
  
  assert(S != NULL);

  assert(cov->initialised);

  for (i=0; i<3; i++) {
    v[i] = c[i] * Ext_WM(FABS(S->a[i] * xx), S->nunew[i],
			 cov->q[WM_LOGGAMMA + 2 * i], 0.0);
    if (!PisNULL(BInotinvnu) && nu[i] > MATERN_NU_THRES) {
      double w, y;
      y = FABS(S->aorig[i] * xx * INVSQRTTWO);
      Gauss(&y, cov, &w);
      *v = *v * MATERN_NU_THRES / nu[i] + 
	(1 - MATERN_NU_THRES / nu[i]) * w;
    }
  }
  v[3] = v[i22];
  v[2] = v[i21];

  
  assert(R_FINITE(v[1]));

}

void biWM2D(double *x, model *cov, double *v) {
  int i;
  double
    *c = P(BIc),
    *nu = P(BInu),
    xx = *x;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  assert(cov->initialised);

  for (i=0; i<3; i++) {
    v[i] = c[i] * S->a[i] * Ext_DWM(FABS(S->a[i] * xx), S->nunew[i],
				    cov->q[WM_LOGGAMMA + 2 * i], 0.0);
    if (!PisNULL(BInotinvnu) && nu[i] > MATERN_NU_THRES) {
      double w, y,
	scale = S->aorig[i] * INVSQRTTWO;
      y = FABS(scale * xx);
      DGauss(&y, cov, &w);
      w *= scale;
      *v = *v * MATERN_NU_THRES / nu[i] + 
	(1 - MATERN_NU_THRES / nu[i]) * w;
    }
  }
  v[3] = v[i22];
  v[2] = v[i21];
}


void biWM2DD(double *x, model *cov, double *v) {
  int i;
  double
    *c = P(BIc),
    *nu = P(BInu),
    xx = *x;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  assert(cov->initialised);

  for (i=0; i<3; i++) {
    v[i] = c[i] * S->a[i] * S->a[i] *
      Ext_DDWM(FABS(S->a[i] * xx), S->nunew[i], cov->q[WM_GAMMA + 2 * i],0.0);
    if (!PisNULL(BInotinvnu) && nu[i] > MATERN_NU_THRES) {
      double w, y,
    scale = S->aorig[i] * INVSQRTTWO;
      y = FABS(scale * xx);
      DDGauss(&y, cov, &w);
      w *= scale;
      *v = *v * MATERN_NU_THRES / nu[i] +
    (1 - MATERN_NU_THRES / nu[i]) * w;
    }
  }
  v[3] = v[i22];
  v[2] = v[i21];
}

void biWM2D3(double *x, model *cov, double *v) {
  int i;
  double
    *c = P(BIc),
    *nu = P(BInu),
    xx = *x;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  assert(cov->initialised);

  for (i=0; i<3; i++) {
    v[i] = c[i] * S->a[i] * S->a[i] * S->a[i] *
      Ext_D3WM(FABS(S->a[i] * xx), S->nunew[i], cov->q[WM_GAMMA + 2 * i], 0.0);
    if (!PisNULL(BInotinvnu) && nu[i] > MATERN_NU_THRES) {
      double w, y,
    scale = S->aorig[i] * INVSQRTTWO;
      y = FABS(scale * xx);
      D3Gauss(&y, cov, &w);
      w *= scale;
      *v = *v * MATERN_NU_THRES / nu[i] +
    (1 - MATERN_NU_THRES / nu[i]) * w;
    }
  }
  v[3] = v[i22];
  v[2] = v[i21];
}

void biWM2D4(double *x, model *cov, double *v) {
  
  int i;
  double
    *c = P(BIc),
    *nu = P(BInu),
    xx = *x;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  assert(cov->initialised);

  for (i=0; i<3; i++) {
    double Sa2 = S->a[i] * S->a[i];
    v[i] = c[i] * Sa2 * Sa2 * Ext_D4WM(FABS(S->a[i] * xx), S->nunew[i],
				       cov->q[WM_GAMMA + 2 * i],0.0);
    if (!PisNULL(BInotinvnu) && nu[i] > MATERN_NU_THRES) {
      double w, y,
    scale = S->aorig[i] * INVSQRTTWO;
      y = FABS(scale * xx);
      D4Gauss(&y, cov, &w);
      w *= scale;
      *v = *v * MATERN_NU_THRES / nu[i] +
    (1 - MATERN_NU_THRES / nu[i]) * w;
    }
  }
  v[3] = v[i22];
  v[2] = v[i21];
  
}


int checkbiWM2(model *cov) { 

  // !! nudiag, nured has priority over nu
  // !! cdiag, rhored has priority over c
  int err;
  gen_storage s;
  gen_NULL(&s);
  s.check = true;

  assert(PisNULL(BIrhored) || ISNAN(P0(BIrhored)) || P0(BIrhored) <= 1);

  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);

  if (cov->Sbiwm == NULL) {
    ONCE_NEW_STORAGE(biwm);
    biwm_storage *S = cov->Sbiwm;
    S->nudiag_given = !PisNULL(BInudiag);
    S->cdiag_given = !PisNULL(BIcdiag);
  }

  //  PMI(cov);  printf("%d\n", cov->Sbiwm->nudiag_given);
  if (cov->q == NULL) QALLOC(6);
  if ((err=initbiWM2(cov, &s)) != NOERROR) {
    biwm_storage *S = cov->Sbiwm;
    if (S->nudiag_given) {
      PFREE(BInu);
    } else {
      PFREE(BInured);
      PFREE(BInudiag);
     }
    if (S->cdiag_given) {
      PFREE(BIc);
    } else {
      PFREE(BIrhored);
      PFREE(BIcdiag);
    }
  }

  VDIM0 = VDIM1 = 2;  

  RETURN_ERR(err);
}

sortsofparam sortof_biwm2(model *cov, int k, int VARIABLE_IS_NOT_USED row,
			  int VARIABLE_IS_NOT_USED col,
			  sort_origin origin){
  biwm_storage *S = cov->Sbiwm;
  if (S == NULL) return UNKNOWNPARAM;
  switch(k) {
  case BInudiag : case BInured : 
    return S->nudiag_given || origin != original ? CRITICALPARAM : IGNOREPARAM;
  case BInu :
    return S->nudiag_given || origin == mle_conform ? IGNOREPARAM : ONLYRETURN;
  case BIcdiag :  
    return S->cdiag_given || origin != original ? VARPARAM : IGNOREPARAM;
  case BIrhored : 
    return S->cdiag_given || origin != original ? ANYPARAM : IGNOREPARAM;
  case BIc :
    return S->cdiag_given || origin == mle_conform ? IGNOREPARAM : ONLYRETURN;
  case BIs :      return SCALEPARAM;
  case BInotinvnu :return ONLYRETURN;
  default : BUG;
  }
}
  
 
sortsofparam sortof_biwm2_INisOUT(model *cov, int k, int VARIABLE_IS_NOT_USED row,
			  int VARIABLE_IS_NOT_USED col){
  biwm_storage *S = cov->Sbiwm;
  if (S == NULL) return UNKNOWNPARAM;
  switch(k) {
  case BInudiag : return S->nudiag_given ? CRITICALPARAM : CRITONLYMLE;
  case BInured :  return S->nudiag_given ? CRITICALPARAM : CRITONLYMLE;
  case BInu :     return S->nudiag_given ? IGNOREPARAM : ONLYRETURN;
  case BIcdiag :  return S->cdiag_given ? VARPARAM : VARONLYMLE;
  case BIrhored : return S->cdiag_given ? ANYPARAM : ONLYMLE;
  case BIc :      return S->cdiag_given ? IGNOREPARAM : ONLYRETURN;
  case BIs :      return SCALEPARAM;
  case BInotinvnu :return ONLYRETURN;
  default : BUG;
  }
}
  
 
void rangebiWM2(model VARIABLE_IS_NOT_USED *cov, range_type *range){

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
  range->min[BIrhored] = - 1;
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
  
  booleanRange(BInotinvnu);

}



void coinitbiWM2(model *cov, localinfotype *li) {
  double
    thres = 1.5,
    //         *c = P(BIc),
    *nu = P(BInu);
  
  if ( ( nu[0] <= thres ) &&  ( nu[1] <= thres ) && ( nu[2] <= thres ) ) {
    li->instances = 1;
    li->value[0] = 1; //  q[CUTOFF_A]
    li->msg[0] =MSGLOCAL_OK;
  }
  else {
    li->instances = 1;
    li->value[0] = 1 ; //  q[CUTOFF_A]
    li->msg[0] = MSGLOCAL_JUSTTRY;
  }
}



// PARS WM


/* Whittle-Matern or Whittle or Besset */ 
// only lower parts of the matrices, including the diagonals are used when estimating !!

#define PARSWM_V vdim2
void kappa_parsWM(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i==PARSnudiag) {
    *nr = 0; 
    *nc = 1;
  } else  *nc = *nr = OUT_OF_RANGE;
}


int initparsWM(model *cov, gen_storage VARIABLE_IS_NOT_USED  *s) {
  // !! nudiag, nured has priority over nu
  // !! cdiag, rhored has priority over c

  double 
    dim = (double) OWNLOGDIM(0), 
    d2 = dim * 0.5, 
    *nudiag = P(PARSnudiag);
  int i, j, vdiag,
    vdim = cov->nrow[PARSnudiag], // hier noch nicht cov->vdim, falls initparsWM von check aufgerufen wird, und cov->vdim noch nicht gesetzt ist
    vdim2 = vdim * vdim, 
    vdimP1 = vdim + 1;
 
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    for (j=i; j<vdim; j++) {
      cov->q[vdiag + j - i] = cov->q[vdiag + vdim * (j-i)] =
	lgammafn(0.5 * (nudiag[i] + nudiag[j]));
    }
  }
 
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    cov->q[vdiag + PARSWM_V] = 1.0;
    for (j=i+1; j<vdim; j++) {
      double half = 0.5 * (nudiag[i] + nudiag[j]);
      int idx = vdiag + j - i;
      cov->q[idx + PARSWM_V] = cov->q[vdiag + vdim * (j-i) + PARSWM_V] =
	EXP(0.5 * (lgammafn(nudiag[i] + d2) + lgammafn(nudiag[j] + d2)
		   -cov->q[vdiag] /* i */ -cov->q[j * vdimP1] /* j */
		   + 2.0 * (cov->q[idx] /* half */ - lgammafn(half + d2))
		   ));
     }
  }
  RETURN_NOERROR;
}

void parsWM(double *x, model *cov, double *v) {
  int i, j, vdiag,
    vdim = VDIM0,
    vdim2 = vdim * vdim,
    vdimP1 = vdim + 1;
  double 
    *nudiag = P(PARSnudiag);

   
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    for (j=i; j<vdim; j++) {
      double half = 0.5 * (nudiag[i] + nudiag[j]);      
      int idx = vdiag + j - i;
      v[idx] = v[vdiag + vdim * (j-i)] = cov->q[PARSWM_V + idx] *
	Ext_WM(*x, half, cov->q[idx], 0.0);
    }
  }
}


void parsWMD(double *x, model *cov, double *v) {
  int i, j, vdiag,
    vdim = VDIM0,
    vdim2 = vdim * vdim,
    vdimP1 = vdim + 1;
  double 
    *nudiag = P(PARSnudiag);
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    for (j=i; j<vdim; j++) {
      double half = 0.5 * (nudiag[i] + nudiag[j]);      
      int idx = vdiag + j - i;
      v[idx] = v[vdiag + vdim * (j-i)] = cov->q[PARSWM_V + idx] *
	Ext_DWM(*x, half, cov->q[idx], 0.0);
    }
  }
}



int checkparsWM(model *cov) { 
 
  double
    *nudiag = P(PARSnudiag);
  
  int i, err, 
    vdim = cov->nrow[PARSnudiag],
    vdimSq = vdim * vdim;
 
  VDIM0 = VDIM1 = vdim;
  if (vdim == 0) SERR1("'%.50s' not given", KNAME(PARSnudiag));
  if ((err = checkkappas(cov, true)) != NOERROR) RETURN_ERR(err);
  
  cov->full_derivs = cov->rese_derivs = 1; 
  for (i=0; i<vdim; i++) {
    int derivs = nudiag[i] - 1.0 > MAXINT ? MAXINT : (int) (nudiag[i] - 1.0);
    if (cov->full_derivs < derivs) cov->full_derivs = derivs;
  }

  if (cov->q == NULL) {
    QALLOC(vdimSq * 2);
    initparsWM(cov, NULL);
  }

  /*
    #define dummyN (5 * ParsWMMaxVDim)
    double value[ParsWMMaxVDim], ivalue[ParsWMMaxVDim], 
    dummy[dummyN], 
    min = RF_INF;
    int j,
    ndummy = dummyN;
  */

  
  RETURN_NOERROR;
}

void rangeparsWM(model VARIABLE_IS_NOT_USED *cov, range_type *range){

  // *nudiag = P0(PARSnudiag], 
  range->min[PARSnudiag] = 0.0;
  range->max[PARSnudiag] = RF_INF;
  range->pmin[PARSnudiag] = 1e-1;
  range->pmax[PARSnudiag] = 4.0; 
  range->openmin[PARSnudiag] = true;
  range->openmax[PARSnudiag] = true;
}
