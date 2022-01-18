/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

 Copyright (C) 2005 -- 2017 Martin Schlather
 
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
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "extern.h"
#include "questions.h"
#include "operator.h"
#include "Processes.h"
#include "variogramAndCo.h"
#include "rf_interfaces.h"
#include "kleinkram.h"
#include "startGetNset.h"

#define BINARY_P 0
#define BINARY_CORR 1
#define BINARY_CENTRED 2

#define MAXCOMPONENTS 6

/*
  Felix' integral 
  + trafo  v = (1-r) / (1+r)
  + taylor 
  + GR 2.211/2

  V ~ 25 : Grenze zur Asymptoic
  GR 8.357


fctn = function(rho, t) EXP(-t^2 / (1+rho)) / SQRT(1.0-rho^2)
C = function(r, t) (2*pi)^(-1) * integrate(fctn, 0, r, t=t)$value
fctn1 = function(rho, t) EXP(-t^2 * rho / 2) / SQRT(rho) / (1 + rho)
C1 = function(r, t)  EXP(-t^2 / 2) / (2*pi) * integrate(fctn1, (1-r)/(1+r), 1, t=t)$value
deltaC = function(r, t) C(r,t) - C1(r,t)
deltaC(0.5, 0.5)

for (r in seq(0, 1, len=10)) for (t in seq(0.1, 2, len=10)) print(deltaC(r,t))#  //
                      


for (r in seq(0.1, 1, len=10)) for (t in seq(0.1, 2, len=10)) {
  print(deltaC(r,t))#  //
 
  a = t^2 /2

  n = 10
  d = 0:n
  fak = c(1, cumprod(d[-1]))
  s = cumsum(a^d / fak)
  
  v = (1-r) / (1+r)
  (U = 2 * SQRT(v) * sum((-v)^d / (2 * d + 1) * s))
  v =  1
  (O = 2 * SQRT(v) * sum((-v)^d / (2 * d + 1) * s))
  (O - U) / (2 * pi)
  C(r, t)
  
  v = (1-r) / (1+r)
  (U = 2 * ATAN(SQRT(v)) + 2 * SQRT(v) * sum((-v)^d / (2 * d + 1) * (EXP(-a) * s - 1)) )
  v =  1
  (O = 2 * ATAN(SQRT(v)) + 2 * SQRT(v) * sum((-v)^d / (2 * d + 1) * (EXP(-a) * s - 1)) )
  (O - U) / (2 * pi) 
  C(r, t)
  print( (O - U) / (2 * pi)  - C1(r,t))#  //


model = RM e x p()
x = seq(0.1, 1, len=10)
t = 1.67
for (r in x) print(C1(r,t=t))#  //
RFcov(m odel, LOG(x))
RFcov(RMbernoulli(m odel, t=t), LOG(x))
for (r in x) print(C1(r,t=t))  #  //

 */


#define Binary_Eps 1e-13
void binary(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  double a, var, V, r, expMa, Vd,
    s, Oned, sum, sumOne, summand, summandOne, 
    d,
    factor,
    t = P0(BINARY_P),
    p = pnorm(t, 0, 1.0, true, false);
  
  COV(ZERO(next), next, &var);
  COV(x, next, &r);
  
  if (t == 0.0) {
    *v = ((M_1_PI * ASIN(r / var) + 0.5) - p) * p;
  } else {      
    a = 0.5 * t * t / var;
    expMa = EXP(-a);
    r /= var;
    
    // als BA-Arbeit
    if (r < -0.9)
      ERR0("correlation of submodel must be >= -0.9 for numerical reasons");
    
    V = (1 - r) / (1 + r);
    
    double ad = expMa;
    d = 0.0;
    sum = sumOne = 0.0;
    Vd = 1;
    Oned = 1;
    s = ad;
    factor = (s - 1.0) / (2.0 * d + 1.0);
    summand = Vd * factor;
    summandOne = Oned * factor;
    while (FABS(summand) > Binary_Eps || FABS(summandOne) > Binary_Eps) {
      sum += summand;
      sumOne += summandOne;
      d += 1.0;
      ad *= a / d;
      s += ad;
      Vd *= -V;
      Oned *= - 1;
      factor =  (s - 1.0) / (2.0 * d + 1.0);
      summand = Vd * factor;
      summandOne = Oned * factor;
    }
    sum += summand;
    sumOne += summandOne;

    // 0.25 = 2 ATAN(SQRT(1.0)) / 2pi
    *v = 0.25 + INVPI * (sumOne - (ATAN(SQRT(V)) + SQRT(V) * sum));
  }

  if (!P0INT(BINARY_CENTRED)) {
    *v += p * p;
  }

  if (P0INT(BINARY_CORR)) {
    *v /= p;
  }
}


int checkbinary(model *cov) {
  // to do extend to multivariate
  
  model
    *next = cov->sub[0];
  double v;
  int i,
    vdim = VDIM0,
    err = NOERROR;
  if (VDIM0 != VDIM1) BUG;
  kdefault(cov, BINARY_P, 0.0);
  //  if (P0(BINARY_P) != 0.0) SERR("currently only threshold=0 is possible");//todo
  kdefault(cov, BINARY_CORR, 1);
  kdefault(cov, BINARY_CENTRED, 1);
  if ((err = CHECK_PASSTYPE(next, PosDefType )) != NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, cov->tsdim,  cov->xdimprev, PosDefType,
  //		     OWNDOM(0), OWNISO(0),
  //		     SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
  setbackward(cov, next);
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = 1.0;
  COV(ZERO(next), next, &v);
  RETURN_NOERROR;
}


void rangebinary(model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[BINARY_P] = RF_NEGINF; 
  range->max[BINARY_P] = RF_INF;
  range->pmin[BINARY_P] = -4.0;
  range->pmax[BINARY_P] = 4.0;
  range->openmin[BINARY_P] = false;
  range->openmax[BINARY_P] = false;

  booleanRange(BINARY_CORR); 
  booleanRange(BINARY_CENTRED); 
}



#define MASTEIN_NU 0
#define MASTEIN_DELTA 1
#define MASTEIN_MAXNU 80
void MaStein(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  double  nuG, loggammas, v1, v2,
    nu = P0(MASTEIN_NU),
    delta = P0(MASTEIN_DELTA);

  COV(ZERO(next), next, &v1);
  COV(x + 1, next, &v2);
  nuG = nu + v1 - v2;
  if (nuG >= (double) MASTEIN_MAXNU) {
    ERR0("Whittle Matern function cannot be evaluated with parameter value b+g(t) greater than 80.");
  }
  loggammas = lgammafn(nu + delta) - lgammafn(nu) -lgammafn(nuG + delta);
  double s = *x,
    bk[MASTEIN_MAXNU + 1L];
  if (s == 0.0) *v = EXP(lgammafn(nuG) + loggammas);
  else *v = 2.0 * EXP(nuG * LOG(0.5 * s) + loggammas
		      + LOG(bessel_k_ex(s, nuG, 2.0, bk)) - s);
}

int check_MaStein(model *cov) {
  model *next = cov->sub[0];
  int err;
  ASSERT_ONESYSTEM;

  if (OWNTOTALXDIM != 2) SERR("reduced dimension must be 2");

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  if ((err = CHECK(next, 1, 1, VariogramType, XONLY,
		     SYMMETRIC, SCALAR, EvaluationType)) != NOERROR) RETURN_ERR(err);
  if (cov->ncol[MASTEIN_NU] != 1 || cov->nrow[MASTEIN_NU] != 1) 
    SERR("nu not scalar");
  if (cov->ncol[MASTEIN_DELTA] != 1 || cov->nrow[MASTEIN_DELTA] != 1) 
    SERR("d not scalar");
  // no setbackward !
  set_maxdim(OWN, 0, MAXDIM(NEXT, 0));
  RETURN_NOERROR;
  // no setbackward !
}

void range_MaStein(model *cov, range_type *range){
  range->min[MASTEIN_NU] = 0.0; // 
  range->max[MASTEIN_NU] = RF_INF;
  range->pmin[MASTEIN_NU] = 1e-2;
  range->pmax[MASTEIN_NU] = 10.0;
  range->openmin[MASTEIN_NU] = true;
  range->openmax[MASTEIN_NU] = true;

  range->min[MASTEIN_DELTA] = 0.5 * (double) (OWNLOGDIM(0) - 1); // d
  range->max[MASTEIN_DELTA] = RF_INF;
  range->pmin[MASTEIN_DELTA] = range->min[MASTEIN_DELTA];
  range->pmax[MASTEIN_DELTA] = 10;
  range->openmin[MASTEIN_DELTA] = false;
  range->openmax[MASTEIN_DELTA] = true;
}



/* bivariate shift model */
#define SHIFT_DELAY 0
void kappashift(int i, model *cov, int *nr, int *nc){
  *nc = 0;
  *nr = i < DefList[COVNR].kappas  ? OWNLOGDIM(0) : -1;
}
void shift(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  double y[ShiftMaxDim],
    z[ShiftMaxDim] = { RF_NAN },
    *jh, *ih, 
    *pv = v,
    *h = P(SHIFT_DELAY);
  int i, j, d,
    logicaldim = OWNLOGDIM(0),
    vdim = VDIM0,
    vdimM1 = vdim - 1,
    vdimSq = vdim * vdim,
    vdimP1 = vdim + 1;
  
  COV(x, next, v);
  for (i=vdimP1; i<vdimSq; i+=vdimP1) v[i] = v[0];
  
  for (jh=h-logicaldim, j=-1; j<vdimM1; j++, jh+=logicaldim) {
    if (j < 0) for (d=0; d<logicaldim; d++) z[d] = x[d];
    else for (d=0; d<logicaldim; d++) z[d] = x[d] + jh[d];
    for (ih=h-logicaldim, i=-1; i<vdimM1; i++, ih+=logicaldim, pv++) {
      if (i==j) continue;
      if (i < 0) for (d=0; d<logicaldim; d++) y[d] = z[d];
      else for (d=0; d<logicaldim; d++) y[d] = z[d] - ih[d];  
      COV(y, next, pv);
    }
  }
}

int checkshift(model *cov) {
  model *next = cov->sub[0];
  int err;
      
  if (OWNTOTALXDIM > ShiftMaxDim)
   SERR2("For technical reasons max. dimension for ave is %d. Got %d.", 
	  StpMaxDim, OWNXDIM(0));

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  COPYALLSYSTEMS(PREVSYSOF(next), OWN, false);

  if ((err = CHECK_GEN(next, SCALAR, SCALAR, EvaluationType, true))
      != NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, dim, dim, PosDefType, XONLY,
  //		     (dim > 1) ? SYMMETRIC : ISOTROPIC, 
  //		     SCALAR, EvaluationType)) != NOERROR) 
  //    RETURN_ERR(err);
  setbackward(cov, next);
  VDIM0 = VDIM1 = cov->ncol[SHIFT_DELAY] + 1;
  RETURN_NOERROR;
}

void rangeshift(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SHIFT_DELAY] = RF_NEGINF;
  range->max[SHIFT_DELAY] = RF_INF;
  range->pmin[SHIFT_DELAY] = - 1000;
  range->pmax[SHIFT_DELAY] = 1000;
  range->openmin[SHIFT_DELAY] = true;
  range->openmax[SHIFT_DELAY] = true;
}



/* Vector - Delta Delta^T */
#define VECTOR_A 0
#define VECTOR_D 1
void vector(double *x, model *cov, double *v) {
  /*
    r = \| x\|
    d/d x_i f(r) = f'(r) x_i / r
    d^2/d x_i^2 f(r) = f''(r) x_i^2/r^2 + f'(r) (r^2 - x_i^2) / r^3
    d^2/d x_i x_j f(r) = f''(r) x_i x_j/r^2 + f'(r) (- x_i x_j) / r^3

    r = 0:
    d/d x_i f(r) |_0 = 0
    d^2/d x_i^2 f(r) |_0 = f''(0)
    d^2/d x_i x_j f(r) = 0

    Vector (n = dimension of spatial field)
    r=0
    L = n * f''(0)

    operator:  -0.5 * (a + 1) laplace + a * hessian,  a in [-1,1]; a=1
  */
  model *next = cov->sub[0];
  double norm[2], normSq0, normL2, normT2, D, D2, diag,
    a = P0(VECTOR_A),
    b = - 0.5 * (1.0 + a);
  int i,
    logicaldim = OWNLOGDIM(0),
    dim = P0INT(VECTOR_D),
    dimP1 = dim + 1,
    dimsq = dim * dim;
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<logicaldim; i++) normT2 += x[i] * x[i];
  if (isIsotropic(NEXT)) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = SQRT(normT2); 
  }
  norm[0] = SQRT(normSq0);
  Abl1(norm, next, &D);
  Abl2(norm, next, &D2);

  if (normSq0 == 0.0) {
    diag = (b * dim + a) * D2;
    for (i=0; i<dimsq; i++) {
      v[i] = (i % dimP1 == 0) ? diag : 0.0; 
    }
  } else {
    double 
      D2nsq = D2 / normSq0,
      D1n = D / norm[0],
      D1n3 = D / (norm[0] * normSq0),
      delta = a * (D2nsq - D1n3),
      Vector = b * ((D2nsq - D1n3) * normL2 + dim * D1n);
    int k,l;
    
    diag = Vector + a * D1n;    
    for (i=k=0; k<dim; k++) {
      for (l=0; l<dim; l++, i++) {
	v[i] = (i % dimP1 == 0) ? diag : 0.0; 
	v[i] += delta * x[k] * x[l];
      }
    }
  }
}


/* Vector - Delta Delta^T */
void vectorAniso(double *x, model *cov, double *v) {
  /*
    operator:  -0.5 * (a + 1) laplace + a * hessian,  a in [-1,1]; a=1
  */

  model *next = cov->sub[0];
  double laplace,
    a = P0(VECTOR_A),
    b = - 0.5 * (1.0 + a);
  int i, j, k, endfor,
    dim = P0INT(VECTOR_D),    
    dimP1 = dim + 1,
    dimsq = dim * dim,
    xdim = OWNXDIM(0),
    xdimsq = xdim * xdim,
    xdimP1 = xdim + 1,
    dimxdim = dim * xdim;

  TALLOC_XX1(D, xdimsq);
 
  // should pass back the hessian
  HESSE(x, next, D);
  laplace = 0.0;

  for (i=0; i<dimxdim; i+=xdimP1) {
    laplace += D[i];
  }
  laplace *= b;

  k = 0;
  for (i=0; i<dimxdim; i+=xdim) {  
   endfor = i + dim; 
    for (j=i; j<endfor; j++) {
      v[k++] = D[j] * a;
    }
  }
  
  for (i=0; i<dimsq; i+=dimP1) {
    v[i] += laplace;
  }
  END_TALLOC_XX1;
}

int checkvector(model *cov) {
  model *next = cov->sub[0];
  int err, i,
   dim = OWNLOGDIM(0);
   ASSERT_UNREDUCED;
   // isotropy_type isotropy = OWNDOM(0);

  kdefault(cov, VECTOR_A, 0.5); // a
  kdefault(cov, VECTOR_D, equalsSpaceIsotropic(OWN) ? dim - 1 : dim); // Dspace
  //  kdefault(cov, VECTOR_D, 
  //	   (isotropy==DOUBLEISOTROPIC  ) 
  //	   ? dim - 1 : dim); // Dspace
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  //  if ( (isotropy==DOUBLEISOTROPIC) 
  //      && P0INT(VECTOR_D) != dim - 1) {
  if (equalsSpaceIsotropic(OWN) && P0INT(VECTOR_D) != dim - 1) {
    SERR1("for spatiotemporal submodels '%.50s' must be applied to spatial part",
	 NICK(cov));
  }
   
  set_nr(OWN, VECTOR); // wegen nr++ unten; JA NICHT SET_NR !!!
  if ((err = CHECK(next, dim, 1, PosDefType, OWNDOM(0), ISOTROPIC,
		     SCALAR, EvaluationType)) != NOERROR) {
    if ((err = CHECK(next, dim, dim, PosDefType, OWNDOM(0),
		     SYMMETRIC, SCALAR, EvaluationType)) != NOERROR) {
	RETURN_ERR(err);
    }
  }

  setbackward(cov, next);
  int diffpref = MIN(2, PREF_BEST - cov->pref[CircEmbed]);
  if (diffpref > 0) cov->pref[CircEmbed] += diffpref;

  for (i=0; i<dim; i++) cov->mpp.maxheights[i] = RF_NA;

  if (next->full_derivs < 2 && !next->hess) {
       SERR("2nd derivative of submodel not defined (for the given paramters)");
  }
  if (!isSpaceIsotropic(NEXT)) {
    if (!next->hess) SERR("hess matrix not defined");
    set_nr(OWN, COVNR + 1); // JA NICHT SET_NR !!!
  }

  VDIM0 = VDIM1 = P0INT(VECTOR_D);

  EXTRA_STORAGE;
  RETURN_NOERROR;
}

void rangevector(model *cov, range_type *range){
  range->min[VECTOR_A] = - 1.0;
  range->max[VECTOR_A] = 1.0;
  range->pmin[VECTOR_A] = - 1.0;
  range->pmax[VECTOR_A] = 1.0;
  range->openmin[VECTOR_A] = false;
  range->openmax[VECTOR_A] = false;

  range->min[VECTOR_D] = 1.0;
  range->max[VECTOR_D] = OWNLOGDIM(0);
  range->pmin[VECTOR_D] = 1.0;
  range->pmax[VECTOR_D] = range->max[VECTOR_D];
  range->openmin[VECTOR_D] = false;
  range->openmax[VECTOR_D] = false;
}


int zzz=0;
void kappadivcurl(int i, model VARIABLE_IS_NOT_USED *cov,
		  int *nr, int *nc){
  *nc = 1;
  *nr = i == 0  ? SIZE_NOT_DETERMINED : -1;
} 


/* Rot - Delta Delta^T */
void curl(double *x, model *cov, double *v) {
  /*
    r = \| x\|
    d/d x_i f(r) = f'(r) x_i / r
    d^2/d x_i^2 f(r) = f''(r) x_i^2/r^2 + f'(r) (r^2 - x_i^2) / r^3
    d^2/d x_i x_j f(r) = f''(r) x_i x_j/r^2 + f'(r) (- x_i x_j) / r^3

    r = 0:
    d/d x_i f(r) |_0 = 0
    d^2/d x_i^2 f(r) |_0 = f''(0)
    d^2/d x_i x_j f(r) = 0

    Rot (n = dimension of spatial field)
    r=0
    L = n * f''(0)

    operator:  -0.5 * (a + 1) rot + a * hessian,  a in [-1,1]; a=1/2
  */
 
  model *next = cov->sub[0];
  defn *N = DefList + NEXTNR; // !!! nicht gatternr !!!!
  double norm[2], normSq0, normL2, normT2, D, D2, D3, diag,
      a = - 1.0, // curl free
      b = - 0.5 * (1.0 + a);
  int i,
     logicaldim = OWNLOGDIM(0),
      dim = logicaldim,
       // currently no space-time allowed -> BA thesis
      dimP1 = dim + 1,
      dimP2 = dim + 2,
      dimP3 = dim + 3,
      dimP2sqM1 =  dimP2 * dimP2 -1; // for three dimensions much
  // more complicated
 
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<logicaldim; i++) normT2 += x[i] * x[i];
  if (isIsotropic(NEXT)) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = SQRT(normT2); 
  }
  norm[0] = SQRT(normSq0);
  N->D(norm, next, &D);
  N->D2(norm, next, &D2);

  N->D3(norm, next, &D3);

  if (normSq0 == 0.0) {
    for (i=0; i<=dimP2sqM1; i++) v[i] = 0.0;

    N->cov(norm, next, v); // v[0]
    
    diag = (b * dim + a) * D2;
    for (i=dimP3; i<dimP2sqM1; i+=dimP3) v[i] = diag; // diag 2,3

    N->D2(norm, next, v + dimP1);
    v[dimP1] *= 2.0;
    v[dimP2 * dimP1] = v[dimP1];
    N->D4(norm, next, v + dimP2sqM1);
    v[dimP2sqM1] *= (8.0 / 3.0); 
  } else {
    double 
      z[2],
      D2nsq = D2 / normSq0,
      D1n = D / norm[0],
      D3n = D3 / norm[0],
      D1n3 = D / (norm[0] * normSq0),
      delta = a * (D2nsq - D1n3),
      div = b * ((D2nsq - D1n3) * normL2 + dim * D1n);
    int k,l;
    
    N->cov(norm, next, v); // v[0]
    
    z[0] = x[0];
    z[1] = x[1];
    for (i=1; i<=dim; i++) { // kante links und oben
	v[i] = -(v[i * dimP2] = z[i-1] * D1n);
    }
 
    diag = div + a * D1n;    
    for (i= dimP3, k=0; k<dim; k++, i+=2) {
      for (l=0; l<dim; l++, i++) {
	v[i] = (i % dimP3 == 0) ? diag : 0.0; 
	v[i] += delta * x[k] * x[l];
      }
    }

    v[dimP1 * dimP2] = v[dimP1] = - v[dimP3] - v[2 * dimP3]; // ecken

    for (i=1; i<=dim; i++) { // kanten unten und hinten
	v[dimP2 * dimP1 + i] =
	    - (v[dimP2 * (i+1) - 1] = z[i-1] * (D3n + D2nsq - D1n3));
    }

    N->D4(norm, next, v + dimP2sqM1); // rechts unten 
    v[dimP2sqM1] += 2.0 * D3n - D2nsq + D1n3;
 
  } 

// for (i=0; i<=15; i++) {
//      if (i!=15) v[i] = 0.0;
//  }

}

#define DIVCURL_WHICH 0
/* Div - Delta Delta^T */
void diverge(double *x, model *cov, double *w) { // div -free !!
  /*
    r = \| x\|
    d/d x_i f(r) = f'(r) x_i / r
    d^2/d x_i^2 f(r) = f''(r) x_i^2/r^2 + f'(r) (r^2 - x_i^2) / r^3
    d^2/d x_i x_j f(r) = f''(r) x_i x_j/r^2 + f'(r) (- x_i x_j) / r^3

    r = 0:
    d/d x_i f(r) |_0 = 0
    d^2/d x_i^2 f(r) |_0 = f''(0)
    d^2/d x_i x_j f(r) = 0

    Div (n = dimension of spatial field)
    r=0
    L = n * f''(0)

    operator:  -0.5 * (a + 1) div + a * hessian,  a in [-1,1]; a=1
  */
  model *next = cov->sub[0];
  
  defn *N = DefList + NEXTNR; // gatter oder keine gatter??
  double
    norm[2], normSq0, normL2, normT2, D, D2, D3, diag,
      a = 1.0, // divergenz free
      b = - 0.5 * (1.0 + a);
  int i,
      logicaldim = OWNLOGDIM(0),
      dim = logicaldim, 
      // currently no space-time allowed -> BA thesis
      dimP1 = dim + 1,
      dimP2 = dim + 2,
      dimP3 = dim + 3,
      dimP2sqM1 =  dimP2 * dimP2 -1; // for three dimensions much
  // more complicated
  double extra_a[MAXCOMPONENTS * MAXCOMPONENTS];
  double *v = PisNULL(DIVCURL_WHICH) ? w : extra_a;
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<logicaldim; i++) normT2 += x[i] * x[i];
  if (isIsotropic(NEXT)) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = SQRT(normT2); 
  }
  norm[0] = SQRT(normSq0);
  N->D(norm, next, &D);
  N->D2(norm, next, &D2);
  N->D3(norm, next, &D3);


  if (normSq0 == 0.0) {
    for (i=0; i<=dimP2sqM1; i++) v[i] = 0.0;

    N->cov(norm, next, v); // v[0]
    
    diag = (b * dim + a) * D2;
    for (i=dimP3; i<dimP2sqM1; i+=dimP3) v[i] = diag;  // diag 2,3

    N->D2(norm, next, v + dimP1);
    v[dimP1] *= 2.0;
    v[dimP2 * dimP1] = v[dimP1];
    N->D4(norm, next, v + dimP2sqM1);
    v[dimP2sqM1] *= (8.0 / 3.0);    
  } else {
    double 
      z[2],
      D2nsq = D2 / normSq0,
      D1n = D / norm[0],
      D3n = D3 / norm[0],
      D1n3 = D / (norm[0] * normSq0),
      delta = a * (D2nsq - D1n3),
      div = b * ((D2nsq - D1n3) * normL2 + dim * D1n);
    int k,l;
    
    N->cov(norm, next, v); // v[0]
    
    z[0] = -x[1];
    z[1] = x[0];
    for (i=1; i<=dim; i++) { // kante links und oben
      v[i] = -(v[i * dimP2] = z[i-1] * D1n);
    }
     
    diag = div + a * D1n;    
    for (i= dimP3, k=0; k<dim; k++, i+=2) {
      for (l=0; l<dim; l++, i++) {
	v[i] = delta * x[k] * x[l] + (k == l ? diag : 0.0);
      }
    }

    v[dimP1 * dimP2] = v[dimP1] = - v[dimP3] - v[2 * dimP3]; // ecken

    for (i=1; i<=dim; i++) { // kanten unten und hinten
	v[dimP2 * dimP1 + i] =
	    -(v[dimP2 * (i+1) - 1] =  z[i-1] * (D3n + D2nsq - D1n3));
    }

    N->D4(norm, next, v + dimP2sqM1); // rechts unten 
    v[dimP2sqM1] += 2.0 * D3n - D2nsq + D1n3;
  }

  if (!PisNULL(DIVCURL_WHICH)) {
    int len = NROW(DIVCURL_WHICH),
      size = (int) cov->q[0];
    for (i=0; i<len; i++) {
      int ii = PINT(DIVCURL_WHICH)[i] - 1;
      for (int j=0; j<len; j++) {
	w[i + j * len] = v[ii + (PINT(DIVCURL_WHICH)[j] - 1) * size];
      }
    }
  }
 
  /*
  printf("\n\n%10g %10g \n", x[0], x[1]);//
  for (i=0; i<=15; i++) {
    printf("%10g ", v[i]);//
  }
  */


//  for (i=0; i<=15; i++) {
//      if (i!=15) v[i] = 0.0;
//  }

}

int checkdivcurl(model *cov) {
  model  *next = cov->sub[0];
  // been programmed yet 
  int err, i,
    totaldim = OWNLOGDIM(0),
    spacedim = GetLocspatialdim(cov);

  // statt STATIONARY : VARIOGRAM ?? s. Paper mit Scheuerer

  if ((err = CHECK(next, totaldim, 1, PosDefType, OWNDOM(0), ISOTROPIC,
		   SCALAR, EvaluationType)) != NOERROR) {
    if ((err = CHECK(next, totaldim, 1, PosDefType, OWNDOM(0), DOUBLEISOTROPIC,
		       SCALAR, EvaluationType)) != NOERROR)
    RETURN_ERR(err);
  }

  if (next->full_derivs < 4) SERR("4th derivative of submodel not defined");
  if (totaldim != 2) SERR("currently coded only for dim=2");
  if (!isSpaceIsotropic(NEXT)) SERR("submodel must be spaceisotropic");
  if (spacedim != 2)
    SERR1("model '%.50s' currently coded only for dim=2", NAME(cov));

  
  setbackward(cov, next);
  int diffpref = MIN(2, PREF_BEST - cov->pref[CircEmbed]);
  if (diffpref > 0) cov->pref[CircEmbed] += diffpref;

  int components =  spacedim + 2; // only correct for spacedim == 2!!!
  if (components > MAXCOMPONENTS) SERR("dimension of space too large");
  int nwhich = NROW(DIVCURL_WHICH);
  if (nwhich > 0) {
    for (i=0; i<nwhich; i++) {
      if (PINT(DIVCURL_WHICH)[i] <= 0 || PINT(DIVCURL_WHICH)[i] > components) 
	SERR4("value %.50s[%d]=%d outside range 1,...,%d.",
	      KNAME(i), i+1, PINT(DIVCURL_WHICH)[i], components);
    }
  } else {
    nwhich = components;
  }
  for (i=0; i<totaldim; i++) cov->mpp.maxheights[i] = RF_NA;
  VDIM0 = VDIM1 = nwhich;
  assert(spacedim == 2);
  if (cov->q == NULL) {
    QALLOC(1);
    cov->q[0] = components;
  }
  
  RETURN_NOERROR;
}


void rangedivcurl(model *cov, range_type *range){
  model  *next = cov->sub[0];
  int  dim = OWNLOGDIM(0),
    spacedim = dim - equalsSpaceIsotropic(NEXT),
    max = spacedim + 2;
  if (spacedim != 2)
    ERR0("div and curl currently programmed only for spatial dimension 2.");
  range->min[DIVCURL_WHICH] = 1;
  range->max[DIVCURL_WHICH] = max;
  range->pmin[DIVCURL_WHICH] = 1;
  range->pmax[DIVCURL_WHICH] = max;
  range->openmin[DIVCURL_WHICH] = false;
  range->openmax[DIVCURL_WHICH] = false;
}


/*

#define LP_P 0
void lp(double *x, model *cov, double *v){
  model *next = cov->sub[0];  
  double z,
    p = P0(LP_P);
  int d,
    dim = cov->tsdim;
  for (z=0.0, d=0; d<dim; d++) {
    z += POW(FABS(x[d]), p);
  }
  z = POW(z, 1 / p);
  
  COV(&z, next, v);
}
int checklp(model *cov) {
  bool skipchecks = GLOBAL.general.skipchecks;
  model *next = cov->sub[0];  
  int err;

  kdefault(cov, 0, 1.0);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
 
  if ((err = CHECK(next, cov->tsdim + 1, 1, PosDefType, XONLY,ISOTROPIC,
		     SCALAR, EvaluationType)) 
	!= NOERROR) RETURN_ERR(err);

  // no setbackward
  double p =P0(LP_P);
  if (p==1) {
    cov->maxdim = next->maxdim - 1;  // true for p==1
    if (cov->maxdim == 0) cov->maxdim = 1;
  } else if (p==2) {
    cov->maxdim = next->maxdim;
  } else if (!skipchecks) SERR("p must be 1 or 2.");
  cov->monotone = next->monotone;
  cov->semiseparatelast = false;
  updatepref(cov, next);

  assert(false); // nicht verwendete Fkt

  RETURN_NOERROR;
}
void rangelp(model *cov, range_type *range){
  range->min[LP_P] = 0.0;
  range->max[LP_P] = 2.0;
  range->pmin[LP_P] = 0.01;
  range->pmax[LP_P] = 2.0;
  range->openmin[LP_P] = true;
  range->openmax[LP_P] = false;
}

*/

#define MA1_ALPHA 0
#define MA1_BETA 1
void ma1(double *x, model *cov, double *v){
  model *next = cov->sub[0];  
  double z,
    alpha = P0(MA1_ALPHA),
    theta = P0(MA1_BETA);
  COV(x, next, &z);
  *v = POW(theta / (1 - (1-theta) * z), alpha);  
}
int checkma1(model *cov) {
//  bool skipchecks = GLOBAL.general.skipchecks;
  model *next = cov->sub[0];  
//  double p;
  int err;

  kdefault(cov, 0, 1.0);
  kdefault(cov, 1, 0.5);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
       
  if ((err = CHECK_PASSFRAME(next, EvaluationType)) != NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, OWNDOM(0),
  //		     OWNISO(0), SCALAR, EvaluationType)) 
  //	!= NOERROR) RETURN_ERR(err);

  //cov->semiseparatelast = cov->separatelast =false;
  cov->logspeed = 0.0;
  //  updatepref(cov, next);
  setbackward(cov, next);
  cov->mpp.maxheights[0] = 1.0;
  RETURN_NOERROR;
}
void rangema1(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[MA1_ALPHA] = 0.0;
  range->max[MA1_ALPHA] = RF_INF;
  range->pmin[MA1_ALPHA] = 0.01;
  range->pmax[MA1_ALPHA] = 10;
  range->openmin[MA1_ALPHA] = true;
  range->openmax[MA1_ALPHA] = true;

  range->min[MA1_BETA] = 0.0;
  range->max[MA1_BETA] = 1.0;
  range->pmin[MA1_BETA] = 0.0;
  range->pmax[MA1_BETA] = 1.0;
  range->openmin[MA1_BETA] = true;
  range->openmax[MA1_BETA] = true;
}


void ma2(double *x, model *cov, double *v){
  model *next = cov->sub[0];  
  double z, z0;
  COV(ZERO(next), next, &z0);
  COV(x, next, &z);
  z = z0 - z;
  if (z == 0) *v = 1.0; else *v = (1.0 - EXP(-z)) / z;  
}
int checkma2(model *cov) {
//  bool skipchecks = GLOBAL.general.skipchecks;
  model *next = cov->sub[0];  
//  double p;
  int err;

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
    
  if ((err = CHECK_PASSTYPE(next, VariogramType )) != NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, cov->tsdim, cov->xdimown, VariogramType, OWNDOM(0),
  //		     OWNISO(0), SCALAR, EvaluationType)) 
  //	!= NOERROR) RETURN_ERR(err);

  // cov->semiseparatelast = cov->separatelast =false;
  cov->logspeed = 0.0;
 //  updatepref(cov, next);
  setbackward(cov, next);
  cov->mpp.maxheights[0] = 1.0;

  RETURN_NOERROR;
}




void kappaM(int i, model *cov, int *nr, int *nc){
  *nc = SIZE_NOT_DETERMINED;
  *nr = i < DefList[COVNR].kappas  ? SIZE_NOT_DETERMINED : -1;
}

void M(model *cov, double *M1, double *Z,  double *M2, double *V) {
  // M C M^t
  assert(M2 != NULL);
  model *next = cov->sub[0];
  int
    n = cov->nsub,
    ncol = cov->ncol[M_M],
    nrow = cov->nrow[M_M];

//  printf("nrow=%d ncol=%d\n", nrow, ncol);
//      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K, ALPHA,A,LDA,B,LDB, BETA,C,LDC)
// C := alpha*op( A )*op( B ) + beta*C,
// M : rows of op(A)
// N : cols of op(B)
// K : cols of op(A)= rows of op(B)
// LD-X : first dimension of X   \command{\link{RMflat}} \tab constant in space \cr
      
  if (next->vdim[0] == 1) {
    if (cov->kappasub[M_M] == NULL && n == 1) {
      // V = Z * M1 %*% t(M2)
      //      printf("M:%10g %10g   %10g %10g (%d)\n", M1[0], M1[1], M2[0], M2[1], cov->zaehler);
      double z = *Z;
      assert(nrow <= MAXVDIM);      
      int nrow2 = nrow * nrow;   
      for (int i=0; i<nrow2; i++) V[i] = cov->q[i] * z;
      //printf("\n"); for (int i=0; i<nrow2; i++) printf("V(%d)=%10g ", i, V[i]);
    } else {
      double Mz[MAXVDIM * MAXVDIM];
      double *p = Mz;
      for (int i=0; i<ncol; i++) {
	double zi = Z[i % n];
	for (int j=0; j<nrow; j++) *(p++) = *(M1++) * zi;
      }
      matmult_2ndtransp(Mz, M2, V, nrow, ncol, nrow);
    }
  } else {
    assert(nrow <= MAXVDIM);      
    assert(ncol <= MAXVDIM);
    double Mz[MAXVDIM * MAXVDIM];
    matmult(M1, Z, Mz, nrow, ncol, ncol);
    matmult_2ndtransp(Mz, M2, V, nrow, ncol, nrow);
  }
}

void NoM(double *z, int ncol, int n, double *v) {
  double *p = v;
  *p = z[0];
  for (int i=1; i<ncol; i++) {
    for (int j=0; j<ncol; j++) *(++p) = 0;
    *(++p) = z[i % n];
  }
}

// plus + M + whittle

void Mstat(double *x, model *cov, double *v) { // RMmatrix
  assert(cov->kappasub[M_M] == NULL);
  int
    n = cov->nsub,
    ncol = cov->ncol[M_M];
  assert(ncol <= MAXVDIM);
  double z[MAXVDIM * MAXVDIM];
  for (int i=0; i<n; i++) COV(x, cov->sub[i], z + i);

  if (PisNULL(M_M)) NoM(z, ncol, n, v);
  else M(cov, P(M_M), z, P(M_M), v);
}

void Mnonstat(double *x, double *y, model *cov, double *v){
  int
    n = cov->nsub,
    ncol = cov->ncol[M_M];
  assert(ncol <= MAXVDIM);
  double z[MAXVDIM * MAXVDIM];
   for (int i=0; i<n; i++) NONSTATCOV(x, y, cov->sub[i], z + i);
  if (cov->kappasub[M_M] == NULL) {
    if (PisNULL(M_M)) NoM(z, ncol, n, v);
    else M(cov, P(M_M), z,  P(M_M), v);
  } else {
    assert(cov->nrow[M_M] <= MAXVDIM);
    double M1[MAXVDIM * MAXVDIM], M2[MAXVDIM * MAXVDIM];
    COV(x, cov->kappasub[M_M], M1);
    COV(y, cov->kappasub[M_M], M2);
    M(cov, M1, z, M2, v);
  }
}

int initM(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  if (!PisNULL(M_M) && cov->nsub==1) {
    int vdim = cov->nrow[M_M],
      ncol = cov->ncol[M_M];
    assert(cov->q != NULL);
    matmult_2ndtransp(P(M_M), P(M_M), cov->q, vdim, ncol, vdim);
  }
  RETURN_NOERROR;
}

int checkM(model *cov) {
  assert(COVNR == M_PROC || cov->Splus == NULL || !cov->Splus->keys_given);
  int err, 
    dim = OWNTOTALXDIM;
  //  printf("check\n");

  ASSERT_ONESYSTEM;
 
  model *sub = cov->sub[0];
  if (SUBNR == BIND) {
    if (cov->nsub > 1) SERR("only 1 vector of submodels might be given");
    int nsub = DefList[SUBNR].kappas,
      j = 0;
    for (int i=0; i<nsub; i++) j += sub->kappasub[i] != NULL;
    if (j >= MAXSUB) SERR("too many submodels");
    j = 0;
    for (int i=0; i<nsub; i++) {
       if (sub->kappasub[i] != NULL) {
	cov->sub[j] = sub->kappasub[i];
	SET_CALLING(cov->sub[j], cov);
	sub->kappasub[i] = NULL;
	j++;
       }
    }
    COV_DELETE_WITHOUT_LOC(&sub, cov);
    cov->nsub = j;
  }

  int vdim,
    ncol = 0,
    subvdim = !PisNULL(M_VDIM) ? P0INT(M_VDIM) : cov->nsub > 1 ? cov->nsub : 0;
  model *k = cov->kappasub[M_M];

  if (k != NULL) {
    // for a rfere
    if ((err = CHECK(k, dim, dim, ShapeType, XONLY,
		     CoordinateSystemOf(OWNISO(0)),
		     SUBMODEL_DEP, EvaluationType)) != NOERROR) RETURN_ERR(err);
    vdim = k->vdim[0];
    ncol = k->vdim[1];

    //  printf("vdim=%d %d\n", vdim, ncol);
    
  } else if (!PisNULL(M_M)) {
    vdim = cov->nrow[M_M];
    ncol = cov->ncol[M_M];
    if (cov->q == NULL) {
      int vdim2 = vdim * vdim;
      QALLOC(vdim2);
      initM(cov, NULL);
    }
  } else {
    vdim = subvdim;
  }
  if (subvdim == 0) subvdim = ncol;
  else if (ncol != 0 && subvdim != ncol)
    SERR("number of submodel is different from the number of columns of 'M'")

      // Achtung: nicht am Anfang!
      //  PMI(cov);
      int required_vdim = cov->nsub > 1 ? 1 : subvdim != 0 ? subvdim : VDIM0;
  if ((err = CHECKPOS2NEG(cov->sub[0], required_vdim, EvaluationType,
			  OWNDOM(0)))!=NOERROR)
    RETURN_ERR(err);
  int vdimback = cov->sub[0]->vdim[0];
  if (subvdim == 0) {
    subvdim = vdimback;
    if (vdim == 0) vdim = subvdim;
  } else if (cov->nsub == 1 && subvdim != vdimback) {
    SERR("submodel does not have correct subdimensionality")
      assert(vdim > 0 && subvdim > 0);
  }
  if (vdim > MAXVDIM) {
    SERR2("the maximum multivariate dimension is %d, but %d is given by the user", MAXVDIM, vdim);
  }
  if (subvdim > MAXVDIM) RETURN_ERR(ERRORMAXVDIM);

  bool variogram =  equalsnowVariogram(cov->sub[0]);
  if (isnowPosDef(cov) && variogram) SERR("not a positive definite function");
  setbackward(cov, cov->sub[0]);
  for (int i=1; i<cov->nsub; i++) {
    err = CHECK_PASSTF(cov->sub[i], variogram ? VariogramType : PosDefType,
		       1, EvaluationType);
    if (err != NOERROR) SERR("submodels are not of the same type");
    setbackward(cov, cov->sub[i]);
  }

  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);

  VDIM0 = VDIM1 = vdim; // after setbackward
  cov->ncol[M_M] = subvdim;
  for (int i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_NA;
  cov->pref[Specific] =  vdim < subvdim ? 2 :
    (vdim == subvdim && cov->kappasub[M_M] == NULL) ? 4 : 5;

  if (cov->kappasub[M_M] != NULL)
    cov->pref[CircEmbed] = cov->pref[Sequential] = PREF_NONE;

  //PMI(cov->calling == NULL ? cov : cov->calling->calling == NULL ? cov->calling : cov->calling->calling );
   
  RETURN_NOERROR;
}



bool allowedDM(model *cov) {
  if (cov->kappasub[M_M] != NULL) {
    bool *D = cov->allowedD;
    for (int i=(int) FIRST_DOMAIN; i<=(int) LAST_DOMAINUSER; D[i++] = false);
    D[KERNEL] = true;
    return false;
  }
  return allowedDstandard(cov);
}


bool allowedIM(model *cov) {
  if (cov->kappasub[M_M] != NULL) {
    bool *I = cov->allowedI;
    for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
    I[SYMMETRIC] = I[SPHERICAL_SYMMETRIC] = I[EARTH_SYMMETRIC] = true;
    return false;
  }
  return allowedIstandard(cov);
}


sortsofparam sortof_M(model VARIABLE_IS_NOT_USED *cov,
		      int VARIABLE_IS_NOT_USED  k, int row, int col,
		      sort_origin VARIABLE_IS_NOT_USED original){
  return k == M_M ? (row==col ? SDPARAM : SIGNEDSDPARAM) : INTEGERPARAM;
}
  
void rangeM(model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[M_M] = RF_NEGINF;
  range->max[M_M] = RF_INF;
  range->pmin[M_M] = - 1e5;
  range->pmax[M_M] = 1e5;
  range->openmin[M_M] = true;
  range->openmax[M_M] = true; 
  
  range->min[M_VDIM] = 1;
  range->max[M_VDIM] = MAXVDIM;
  range->pmin[M_VDIM] = 1;
  range->pmax[M_VDIM] = MAXVDIM;
  range->openmin[M_VDIM] = false;
  range->openmax[M_VDIM] = false; 
}



Types TypeM(Types required, model *cov, isotropy_type required_iso){
  if (!(isShape(required) || isTrend(required) || isProcess(required)))
    return BadType;

  model **Sub = cov->sub,
    *sub = Sub[0];
  int nsub = cov->nsub;
  
  if (SUBNR == BIND) {
    Sub = sub->kappasub;
    nsub = DefList[SUBNR].kappas;
  }

  for (int i=0; i<nsub; i++) {
    if (Sub[i] != NULL &&
	TypeConsistency(required, Sub[i], required_iso) == BadType)
      return BadType;
  }
  return required;
}



int structMproc(model *cov, model **newmodel) {
  int err;
  ASSERT_NEWMODEL_NULL;
  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
  if (PrevLoc(cov)->distances) 
    SERR("distances do not allow for more sophisticated simulation methods");
     
  NEW_STORAGE_WITH_SAVE(plus);
  plus_storage *s =cov->Splus;

  int newdim = PREVLOGDIM(0);
  for (int i=0; i<cov->nsub; i++) {
    if ((err = covcpy(s->keys + i, cov->sub[i])) != NOERROR) RETURN_ERR(err);
    assert(!isGaussMethod(s->keys[i]));
    addModel(s->keys + i, GAUSSPROC);
    model *key = s->keys[i]; // hier und nicht frueher!
    assert(key->calling == cov);    
    
    if ((err = CHECK_NO_TRAFO(key, newdim, newdim, ProcessType, XONLY,
			      CoordinateSystemOf(OWNISO(0)),
			      cov->sub[i]->vdim[1], cov->frame)) != NOERROR) {
      RETURN_ERR(err);
    }

    //    printf("structM A %d\n", i);
   
    if ((err = STRUCT(key, NULL)) != NOERROR) RETURN_ERR(err);
    //printf("structM B\n");
  }
  cov->Splus->keys_given = true;
  if ((err = ReturnOwnField(cov)) != NOERROR) RETURN_ERR(err);
  cov->simu.active = true;
  //printf("structM C\n");

  RETURN_NOERROR;
}



int initMproc(model *cov, gen_storage *S){
  int 
    //    prevdim = prevloc->timespacedim,
    //  dim = GetLoctsdim(cov),
    err = NOERROR;

  //  PMI(cov); //crash();
  
  // assert(cov->key != NULL);
  plus_storage *s =cov->Splus;
  assert(s != NULL);

  if ((err = alloc_pgs(cov)) != NOERROR) RETURN_ERR(err);

  for (int i=0; i<cov->nsub; i++) {
    if ((err = INIT(s->keys[i], 0, S)) != NOERROR) RETURN_ERR(err);
    s->keys[i]->simu.active = true;
  }

  EXTRA_STORAGE;

  RETURN_NOERROR;
}


void doMproc(model *cov, gen_storage *S){
  // assert(cov->key != NULL);
  //  printf("mvdim=%d\n", M_VDIM);
  //  PMI0(cov)
  int 
    vdim = VDIM0,
    keyvdim = cov->sub[0]->vdim[0];
  double z[MAXVDIM], zz[MAXVDIM],
    *rf = cov->rf,
    *subrf[MAXSUB] = {NULL};
  assert(keyvdim < MAXVDIM);
  assert(vdim < MAXVDIM);

 
  bool
    subM = cov->kappasub[M_M] != NULL,
    anyM = subM || !PisNULL(M_M),
    close = GLOBAL.general.vdim_close_together;

  plus_storage *s =cov->Splus;
  for (int i=0; i<cov->nsub; i++) {
    assert(s->keys[i] != NULL);
    DO(s->keys[i], S);
    subrf[i] = s->keys[i]->rf;
    assert(subrf[i] != NULL);
  }

  if (!anyM) {
    int tot = Gettotalpoints(cov);
    if (!close) {
      long bytes = tot * sizeof(double);
      for (int j=0; j<vdim; j++, rf+=tot) MEMCOPY(rf, subrf[j], bytes);
    } else {
      for (int i=0; i<tot; i++) for (int j=0; j<vdim; *(rf++) = subrf[j++][i]);
    }
    return;
  }
  
  bool kernel = false;  
  double M0[MAXVDIM * MAXVDIM],
    *M = subM ? M0 : P(M_M),
    *y = cov->Spgs->supportmin,
    *sub = subrf[0];
  assert(sub != NULL);
  double
    *w = close && cov->nsub == 1 ? sub - keyvdim : z; 
  FINISH_START(cov, cov, true, 0);
  
#define todo								\
  if (cov->nsub > 1) for (int j=0; j<keyvdim; j++) w[j] = *(++subrf[j]); \
  else if (!close) { for (int j=0; j<keyvdim; j++) w[j] = sub[i_row +j*tot]; } \
  else w += keyvdim;							\
  if (subM) COV(x, cov->kappasub[M_M], M);				\
  if (close) {								\
    Ax(M, w, vdim, keyvdim, rf);					\
    rf += vdim;								\
  } else {								\
    Ax(M, w, vdim, keyvdim, zz);					\
    for (int j=0; j<vdim; j++) rf[i_row + j * tot] = zz[j];		\
  }

  trafo &= subM;
  vdimSq = 0; // do as if always we have multivariate case
  double *value = NULL, *zero = NULL;

  //PMI(cov);
  
  PERFORM({}, todo, {}, {});

  

  /*
  if (grid) {								
    if (ygiven || kernel) {						
      STANDARDSTART_Y_SUPPL;						
      if (vdimSq == 1) { GRIDCYCLES(; value+=vdimSq);	
      } else { GRIDCYCLES(); }				
    } else { /// grid, y not given
      if (vdimSq == 1) {						
	GRIDCYCLE_X(; value+=vdimSq);			
      } else {
	//	GRIDCYCLE_X(todo);
     while (true) {			       
       if (cov->nsub > 1) for (int j=0; j<keyvdim; j++) w[j] = *(++subrf[j]); 
       else if (!close) for (int j=0; j<keyvdim; j++)   {pr intf("j=%d %d\n", j, i_row); w[j] = sub[i_row +j*tot];}
       else w += keyvdim;							
       if (subM) COV(x, cov->kappasub[M_M], M);				
       if (close) {								
	 Ax(M, w, vdim, keyvdim, rf);					
	 rf += vdim;								
       } else {								
	 Ax(M, w, vdim, keyvdim, zz);					
	 for (int j=0; j<vdim; j++) rf[i_row + j * tot] = zz[j];		
       }
       
    STANDARDINKREMENT_X;		       
  }
   }				
    }									
  } else { // not a grid 					
    int localdim = tsxdim;						
    double *xxx, *yyy = zero;						
    if (trafo) {							
      localdim = TransformLoc(cov, &xx, &yy, false);			
      xxx = xx;								
      if (ygiven) yyy = yy;						
    } else {								
      xxx=loc->x;							
      if (ygiven) yyy=loc->y;						
    }									
    int //  xMem = sizeof(double) * (tsxdim + 1), 
      xCpy = sizeof(double) * localdim;					
    assert(ygiven xor (yyy==zero));					
    if (ygiven || kernel) {						
      double *y0 = yyy,							
	*yend = ygiven ? yyy + tsxdim * loc->ly : yyy;			
      NONGRIDCYCLE(DO_INCREMENTY, PREPAREY, ; value+=vdimSq,	
		   );					
    } else {								
      NONGRIDCYCLE(EMPTY, PREPAREX, ; value+=vdimSq,	
		   todo);					
    }									
    if (err != NOERROR) XERR(err);					
  }

  */

  
  STANDARD_ENDE_X;
}




void kappaScatter(int i, model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i < DefList[COVNR].kappas ? 0 : -1;
}

void Scatter(double *xx, model *cov, double *v){
  model *next = cov->sub[0];
  int
    vdim = VDIM0 * VDIM1,
    tsxdim = OWNTOTALXDIM;
  scatter_storage *s = cov->Sscatter;
   int
     i_row = 0,
    *start = s->min,   				
    *end = s->max;						
   double						
    *inc = s->step;						
  TALLOC_X1(xstart, tsxdim);
  TALLOC_X2(x, tsxdim + 1);
  TALLOC_X3(value, vdim);
  TALLOC_L1(nx, tsxdim);
     

  for (int i=0; i<vdim; i++) v[i]=0.0;
  for (int d=0; d<tsxdim; d++) {
    if (P(SCATTER_STEP)[d] <= 0) {
      BUG;
    } else {
      xstart[d] = (double) s->min[d] * s->step[d] + xx[d]; // S
      // printf("xstart =%10g %10g %10g %10g\n", xstart[d], (double) s->min[d], s->step[d], xx[d]);
    }
  }
   
							
  for (int d=0; d<tsxdim; d++) {				
    nx[d] = start[d];					
    x[d] = xstart[d];					
  }

  while (true) {
    int d;
    // printf("Scatter[vdim=%d] =%10g\n", vdim, *x);    
    COV(x, next, value);
    for (int i=0; i<vdim; i++) v[i] += value[i];
    STANDARDINKREMENT_X;
  }
  END_TALLOC_X1;
  END_TALLOC_X2;
  END_TALLOC_X3;
  END_TALLOC_L1;
}


int TaylorScatter(model *cov) {
  model *next = cov->sub[0];

  if (hasRandomFrame(cov)) {
    for (int i=0; i<=cov->mpp.moments; i++) 
      cov->mpp.mM[i] = cov->mpp.mMplus[i] = RF_NA;

    //    PMI0(cov);    printf("prevxtot=%d \n", PREVTOTALXDIM);
    double *zero = ZERO(cov);
    //    printf("zero = %10g\n", zero[0]);
    Scatter(zero, cov, cov->mpp.maxheights); 
    
    if (next->taylor[0][TaylorPow] < 0.0) {
      cov->taylorN = next->taylorN;
      for (int i=0; i < next->taylorN; i++) 
	for (int t=TaylorConst; t <= TaylorPow; t++)
	  cov->taylor[i][t] = next->taylor[i][t];
    } else {
      cov->taylorN = 1;
      cov->taylor[0][TaylorConst] = cov->mpp.maxheights[0];
      cov->taylor[0][TaylorPow] = 0.0;
    } 
    
    cov->tailN = next->tailN;
    for (int i=0; i < next->tailN; i++) 
      for (int t=TaylorConst; t <= TaylorExpPow; t++)
	cov->tail[i][t] = next->tail[i][t];    
  }
  
  else ILLEGAL_FRAME;

  RETURN_NOERROR;
}

int checkScatter(model *cov) {
  model *next = cov->sub[0];
  mpp_param *gp = &(GLOBAL.mpp);
  int d, nr, err,
    dim = OWNTOTALXDIM,
    vdim = VDIM0 * VDIM1;
  
  if (PisNULL(SCATTER_MAX)) {
    PALLOC(SCATTER_MAX, dim, 1);
    for (d=0; d<dim; d++) PINT(SCATTER_MAX)[d] = gp->scatter_max[d];
  } else if ((nr = cov->nrow[SCATTER_MAX]) < dim) {
    int j, *m = PINT(SCATTER_MAX);
    assert(nr > 0);
    PtoNULL(SCATTER_MAX);
    PALLOC(SCATTER_MAX, dim, 1);
    for (j=d=0; d<dim; d++)  {
      PINT(SCATTER_MAX)[d] = m[j++];
      j = j % nr;
    }
  }
  if (PisNULL(SCATTER_STEP)) {
    PALLOC(SCATTER_STEP, dim, 1);
    for (d=0; d<dim; d++) P(SCATTER_STEP)[d] = gp->scatter_step[d];
  } else if ((nr = cov->nrow[SCATTER_STEP]) < dim) {
    int j;
    double *m = P(SCATTER_STEP);
    assert(nr > 0);
    PtoNULL(SCATTER_MAX);
    PALLOC(SCATTER_MAX, dim, 1);
    for (j=d=0; d<dim; d++)  {
      PINT(SCATTER_MAX)[d] = m[j++];
      j = j % nr;
    }
  }
  COPYALLSYSTEMS(PREVSYSOF(next), OWN, false);	
  set_system_type(PREVSYSOF(next), ShapeType);
  if ((err = CHECK_GEN(next, VDIM0, VDIM1, cov->frame, true)) != NOERROR)
    RETURN_ERR(err);
  // if ((err = CHECK_VDIM(next, cov->tsdim, cov->xdimown, ShapeType, OWNDOM(0),  
  //			OWNISO(0), VDIM0, VDIM1,
  //			cov->frame)) != NOERROR)
  //RETURN_ERR(err);
  
  setbackward(cov, next);
  EXTRA_STORAGE;
  
  if (cov->Sscatter == NULL || cov->Sscatter->dim != dim ||
      cov->Sscatter->vdim != vdim) {
    NEW_STORAGE(scatter);
    cov->Sscatter->vdim = vdim;
    cov->Sscatter->dim = dim;
    ALLC_NEWINT(Sscatter, min, dim, min);
    ALLC_NEWINT(Sscatter, max, dim, max);
    ALLC_NEW(Sscatter, step, dim, step);
    TALLOC_X1(x, dim);
    TALLOC_X2(xmin, dim);
    location_type *loc = Loc(cov);
    for (d=0; d<dim; d++) {
      if (R_finite(P(SCATTER_STEP)[d])) {
 	step[d] = P(SCATTER_STEP)[d];
      } else {
	if (loc->grid) {
	  step[d] = loc->xgr[d][XSTEP];
	} else SERR1("non-positive '%.50s' only allowed for grids", 
		     KNAME(SCATTER_STEP));
      } 
      if (PINT(SCATTER_MAX)[d] < 0) {
	if (P(SCATTER_STEP)[d] > 0) {
	  NONSTATINVERSE(&(GLOBAL.mpp.about_zero), next, xmin, x);
	  BUG; // ehem xmin vorige Zeile kann nicht stimmen?! x doch auch nicht
	} else {
	  min[d] = 0;
	  max[d] = loc->xgr[d][XLENGTH];
	}
      } else {
	if (P(SCATTER_STEP)[d] > 0) {
	  min[d] = -(max[d] = PINT(SCATTER_MAX)[d]);
	} else {
	  min[d] = max[d] = NA_INTEGER;
	}
      }
      max[d]++; // da auf '<' getestet wird
    }
    END_TALLOC_X1;
    END_TALLOC_X2;
  }


  if ((err = TaylorScatter(cov)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}
  
void rangeScatter(model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[SCATTER_STEP] = RF_NEGINF;
  range->max[SCATTER_STEP] = RF_INF;
  range->pmin[SCATTER_STEP] = 1e-5;
  range->pmax[SCATTER_STEP] = 1e5;
  range->openmin[SCATTER_STEP] = true;
  range->openmax[SCATTER_STEP] = true; 

  range->min[SCATTER_MAX] = - 1;
  range->max[SCATTER_MAX] = MAXINT;
  range->pmin[SCATTER_MAX] = 0;
  range->pmax[SCATTER_MAX] = 999;
  range->openmin[SCATTER_MAX] = false;
  range->openmax[SCATTER_MAX] = true;
}



int struct_scatter(model VARIABLE_IS_NOT_USED *cov, 
		   model VARIABLE_IS_NOT_USED  **newmodel){
  RETURN_NOERROR;
}


int init_scatter(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int err;
  //vdim = VDIM0;

  if (VDIM1 != 1) SERR("matrix-valued shape functions cannot be initialised");
  if ((err = TaylorScatter(cov)) != NOERROR) RETURN_ERR(err);
 
  RETURN_NOERROR;
}


void do_scatter(model *cov, gen_storage *s){
  BUG;
  model *key = cov->key;
  DO(key, s);   
}






#define SCHUR_M 0
#define SCHUR_DIAG 1
#define SCHUR_RED 2
void kappaSchur(int i, model *cov, int *nr, int *nc){
  double *M = P(SCHUR_M);
  int vdim = cov->nrow[M != NULL ? SCHUR_M : SCHUR_DIAG]; 
  *nc = i == SCHUR_M ? vdim : 1;
  *nr = i == SCHUR_RED ? vdim * (vdim-1) / 2 
    : i < DefList[COVNR].kappas ? vdim : -1;
}

void SchurMult(double VARIABLE_IS_NOT_USED *x, model *cov, double *v){
  double *M = P(SCHUR_M);
  int i,
    vdim = VDIM0;
  if (!PisNULL(SCHUR_M)) {
    int nrow2 = cov->nrow[SCHUR_M] * cov->nrow[SCHUR_M];
    for (i=0; i<nrow2; i++) v[i] *= M[i]; 
  } else {
    int j, k;
    TALLOC_X1(q, vdim);
    double 
      *diag=P(SCHUR_DIAG),
      *red=P(SCHUR_RED);
    for (i=0; i<vdim; i++) q[i] = SQRT(diag[i]);
    for (k=i=0; i<vdim; i++) {
      for (j=0; j<vdim; j++, k++) {
	v[k] *= q[i] * q[j];
      }
    }
    for (k=i=0; i<vdim; i++) {
      for (j=0; j<vdim; j++, k++) {
	v[i + vdim * j] *= red[k];	
  	v[j + vdim * i] *= red[k];
      }
    }
    END_TALLOC_X1;
  }
}

void Schurstat(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  COV(x, next, v); 
  SchurMult(x, cov, v);
}
void DSchur(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  Abl1(x, next, v); 
  SchurMult(x, cov, v);
}
void D2Schur(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  Abl2(x, next, v); 
  SchurMult(x, cov, v);
}
void D3Schur(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  Abl3(x, next, v); 
  SchurMult(x, cov, v);
}
void D4Schur(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  Abl4(x, next, v); 
  SchurMult(x, cov, v);
}

void Schurnonstat(double *x, double *y, model *cov, double *v){
  model *next = cov->sub[0];
  NONSTATCOV(x, y, next, v);  
  SchurMult(x, cov, v);
}
int checkSchur(model *cov) {
  model *next = cov->sub[0];
  double *C = NULL,
    *M=P(SCHUR_M),
    *diag=P(SCHUR_DIAG),
    *red=P(SCHUR_RED);
  int i,j, k, l,
    err = NOERROR, 
    vdim = cov->nrow[M != NULL ? SCHUR_M : SCHUR_DIAG],
    vdimP1 = vdim + 1,
    bytes = vdim * vdim * sizeof(double),
    *nrow = cov->nrow,
    *ncol = cov->ncol;

  VDIM0 = VDIM1 = vdim;
  if ((err = CHECK_PASSTF(next, PosDefType, nrow[SCHUR_M], EvaluationType))
      != NOERROR) goto ErrorHandling;
  //  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, OWNDOM(0), 
  //		     OWNISO(0), nrow[SCHUR_M], EvaluationType)) != NOERROR) 
  //RETURN_ERR(err);
  setbackward(cov, next);

  if ((M != NULL) xor (diag == NULL || red == NULL))
    GERR3("either '%.50s' and '%.50s' or '%.50s' must be given", 
	  KNAME(SCHUR_DIAG), KNAME(SCHUR_RED), KNAME(SCHUR_M));
  C = (double*) MALLOC(bytes);
  if (M == NULL) {
    for (i=0; i<vdim; i++) 
      if (diag[i]<0) GERR1("elements of '%.50s' negative.", KNAME(SCHUR_DIAG));
    for (k=l=i=0; i<vdim; i++, l+=vdimP1) {
      for (j=0; j<vdim; j++, k++) {
	C[i + vdim * j] = C[i * vdim + j] = red[k];
      }
      C[l] = 1.0;
    }
  } else MEMCOPY(C, M, bytes);
  
  if (!Ext_is_positive_definite(C, ncol[SCHUR_M])) 
    GERR3("%d x %d matrix '%.50s' is not (strictly) positive definite",
	  nrow[SCHUR_M], ncol[SCHUR_M], KNAME(SCHUR_M));
    
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = 1.0;

  /*

    Check symmetry??

 for (vdiag = i=0; i<vdim; i++, vdiag+=vdimP1) {
    if (rhored[vdiag] != 1.0) GERR("diagonal elements do not equal 1.0");
    if (value[i]  <  -GLOBAL.direct.svdtolerance) {
    GERR("' r  hored' has negative eigenvalues");
    }
  }
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    for (j=i+1; j<vdim; j++) {
      if (rhored[ vdiag + j - i] != rhored[vdiag + vdim * (j-i)]) 
	GERR("' rh ored' not s ymmetric");
    }
  }
  */
 ErrorHandling:
   FREE(C);

   EXTRA_STORAGE;
  RETURN_ERR(err);
 }
  
void rangeSchur(model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[SCHUR_M] = RF_NEGINF;
  range->max[SCHUR_M] = RF_INF;
  range->pmin[SCHUR_M] = - 1e5;
  range->pmax[SCHUR_M] = 1e5;
  range->openmin[SCHUR_M] = true;
  range->openmax[SCHUR_M] = true; 

  range->min[SCHUR_DIAG] = 0;
  range->max[SCHUR_DIAG] = RF_INF;
  range->pmin[SCHUR_DIAG] = 0;
  range->pmax[SCHUR_DIAG] = 1e5;
  range->openmin[SCHUR_DIAG] = false;
  range->openmax[SCHUR_DIAG] = true; 

  range->min[SCHUR_RED] = 0;
  range->max[SCHUR_RED] = 1;
  range->pmin[SCHUR_RED] =  0;
  range->pmax[SCHUR_RED] = 1;
  range->openmin[SCHUR_RED] = false;
  range->openmax[SCHUR_RED] = false; 
}


#define ID_VDIM 0
void IdStat(double *x, model *cov, double *v) {
    model *next = cov->sub[0];
    COV(x, next, v);
}
void IdNonStat(double *x, double *y, model *cov, double *v){
  model *next = cov->sub[0];
  NONSTATCOV(x, y, next, v);
}

int checkId(model *cov) {
  model *next = cov->sub[0];
  int err;

  VDIM0 = VDIM1 = !PisNULL(ID_VDIM) ? P0INT(ID_VDIM) : SUBMODEL_DEP;
  if ((err = CHECK_NOPASS(next)) !=NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, OWNDOM(0),
  //		     OWNISO(0), cov->vdim, cov->frame)) !=NOERROR) RETURN_ERR(err);
  if (VDIM0 == SUBMODEL_DEP) {
    VDIM0 = next->vdim[0];
    VDIM1 = next->vdim[1];
  }
  cov->logspeed = next->logspeed;
  setbackward(cov, next);
   
  RETURN_NOERROR;
}

void DId(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  Abl1(x, next, v);
}

void DDId(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  Abl2(x, next, v);
}
void TBM2Id(double *x, model *cov, double *v){  
  model *next = cov->sub[0];
  TBM2CALL(x, next, v)
}
void IdInverse(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  INVERSE(x, next, v); 
}

int initId(model *cov, gen_storage *S) {
  model *next = cov->sub[0];
  return INIT(next, cov->mpp.moments, S);
}
void spectralId(model *cov, gen_storage *S, double *e) { 
  //  spectral_storage *s = &(S->Sspectral);
  model *next = cov->sub[0];
  SPECTRAL(next, S, e); // nicht nr
}
void coinitId(model *cov, localinfotype *li) {
  model *next = cov->sub[0];
  DefList[NEXTNR].coinit(next, li); // nicht nr
}
void ieinitId(model *cov, localinfotype *li) {
  model *next = cov->sub[0];
  DefList[NEXTNR].ieinit(next, li); // nicht nr
}
void rangeId(model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[ID_VDIM] = 1.0;
  range->max[ID_VDIM] = RF_INF;
  range->pmin[ID_VDIM] = 1.0;
  range->pmax[ID_VDIM] = 10;
  range->openmin[ID_VDIM] = false;
  range->openmax[ID_VDIM] = true;
}
Types TypeId(Types required, model *cov, isotropy_type required_iso) {
  model *next = cov->sub[0];
  return TypeConsistency(required, next, required_iso);
}


#define EXP_N 0
#define EXP_STANDARDISED 1
void Exp(double *x, model *cov, double *v, int n, bool standardize){
  double v0, 
    s = 0.0, 
    w = 1.0;
  model *next = cov->sub[0];
  int k,
    vdim = VDIM0,
    vdimq = vdim * vdim;

  COV(x, next, v);
  if (vdim == 1) {
    for (k=0; k<=n; k++) {
      s += w;
      w *= *v / (double) (k+1);
    }    
    *v = EXP(*v) - s;
    if (standardize) {
      Exp(ZERO(cov), cov, &v0, n, false);
      *v /= v0;
    }
  } else {
    int i;
    BUG;
    // fehlt die multiplication von links und rechts mit C(0)^{-1/2}
    for (i=0; i<vdimq; i++) v[i] = EXP(v[i]); 
  }
}

void Exp(double *x, model *cov, double *v){	
  Exp(x, cov, v, P0INT(EXP_N), P0INT(EXP_STANDARDISED));
}


void nonstatExp(double *x, double *y, model *cov, double *v, int n, 
		bool standardised){
  model *next = cov->sub[0];
  double v0, 
    s = 0.0, 
    w = 1.0;
  int k,
    vdim = VDIM0,
    vdimq = vdim * vdim;

  NONSTATCOV(x, y, next, v);
  if (vdim == 1) {
    for (k=0; k<=n; k++) {
      s += w;
      w *= *v / (double) (k+1);
    }    
    *v = EXP(*v) - s;
    if (standardised) {
      double *zero = ZERO(cov);
      nonstatExp(zero, zero, cov, &v0, n, false);
       *v /= v0;
    }
  } else {
    int i;
    BUG;
    // fehlt die multiplication von links und rechts mit C(0)^{-1/2}
    for (i=0; i<vdimq; i++) v[i] = EXP(v[i]); 
  }
}

void nonstatExp(double *x, double *y, model *cov, double *v){	
  nonstatExp(x, y, cov, v, P0INT(EXP_N), P0INT(EXP_STANDARDISED));
}

void DExp(double *x, model *cov, double *v){
  double D;
  model *next = cov->sub[0];
  int n = P0(EXP_N);
  assert(VDIM0 == 1);
  
  Abl1(x, next, &D);
  Exp(x, cov, v, n - 1, false);
  *v *= -D;

  if (P0INT(EXP_STANDARDISED)) {
    double v0;
    Exp(ZERO(cov), cov, &v0, n, false);
    *v /= v0;
  }
}

  // Hier fehlt die multidim variante (Hesse matrix)
  

void DDExp(double *x, model *cov, double *v){
  double D, Abl2, w;
  model *next = cov->sub[0];
  int n = P0INT(EXP_N);
  assert(VDIM0 == 1);

  Abl1(x, next, &D);
  Abl2(x, next, &Abl2);
  Exp(x, cov, v, n - 2, false);
  Exp(x, cov, &w, n - 1, false);

  *v =  D * D * *v + Abl2 * w; // Achtung + D2 nicht - D2, da D2 Ableitung einer Cov.

  if (P0INT(EXP_STANDARDISED)) {
    double v0;
    Exp(ZERO(cov), cov, &v0, n, false);
    *v /= v0;
  }
}

int checkExp(model *cov) {
  model *next = cov->sub[0];
  int err, i,
      vdim = VDIM0;

  kdefault(cov, EXP_N, -1);
  kdefault(cov, EXP_STANDARDISED, 1);
 
  if ((err = CHECKPOS2VAR(next, SCALAR, cov->frame, OWNDOM(0))) != NOERROR)
    RETURN_ERR(err);
  
  //  if ((err = CHECKPD2ND(next, cov->tsdim,  cov->xdimown, OWNISO(0),
  //			SCALAR, EvaluationType)) != NOERROR) RETURN_ERR(err);
  if (!isnowPosDef(next) && P0INT(EXP_N) != - 1)
    SERR("for variograms only n=-1 allowed");
  
  setbackward(cov, next);
  if (VDIM0 > 1 && P0INT(EXP_N) != - 1)
    SERR1("'%.50s' must be '-1' in the multivariate case", KNAME(EXP_N));
  if (VDIM0 > 1) SERR("multivariate case not programmed yet");
 
  if (isXonly(NEXT)) {
   defn *C = DefList + COVNR;
    cov->pref[CircEmbed] = C->pref[CircEmbed];
    cov->pref[Direct] = C->pref[Direct];
    cov->pref[Sequential] = C->pref[Sequential];    
    if (!isnowVariogram(cov))
      SERR1("negative definite function expected -- got '%.50s'",
	    TYPE_NAMES[SYSTYPE(OWN, 0)]);
  } else {
    if (!isnowPosDef(cov)) 
      SERR1("positive definite function expected -- got '%.50s'",
	    TYPE_NAMES[SYSTYPE(OWN, 0)]); 
  }

  double height= isnowVariogram(next) && !isnowPosDef(next) ? 1.0 : RF_NA;
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = height;

  cov->monotone = 
    (isBernstein(next)) ? NORMAL_MIXTURE : 
    isMonotone(next->monotone) ? MONOTONE : NOT_MONOTONE;
  cov->logspeed = 0.0;
  RETURN_NOERROR;
}


void rangeExp(model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[EXP_N] = - 1;
  range->max[EXP_N] = RF_INF;
  range->pmin[EXP_N] = - 1;
  range->pmax[EXP_N] = 5;
  range->openmin[EXP_N] = false;
  range->openmax[EXP_N] = true;

  booleanRange(EXP_STANDARDISED);
}



void Pow(double *x, model *cov, double *v){
  double v0, v1,
    alpha = P0(POW_ALPHA);
  model *next = cov->sub[0];
  if (isnowShape(cov)) {
    //printf("PowShape x=%10g\n", *x);
    COV(x, next, v);
    *v = POW(*v, alpha);
  } else {
   COV(ZERO(next), next, &v0);
    COV(x, next, &v1);
    *v = POW(v0, alpha) - POW(v0 - v1, alpha);
  }
}

void DPow(double *x, model *cov, double *v){
  double v0, v1, gamma,
    alpha = P0(POW_ALPHA);
  model *next = cov->sub[0];

  Abl1(x, next, v);
  if (alpha == 1.0) return;
  COV(ZERO(next), next, &v0);
  if (isnowShape(cov)) {
    *v *= - alpha * POW(v0, alpha -1.0);
  } else {  
    COV(x, next, &v1);
    gamma =  v0 - v1;
    *v *= - alpha * POW(gamma, alpha -1.0);
  }
  // Achtung  "-" , da covarianzfunktion angenommen
}

void DDPow(double *x, model *cov, double *v){
  double D, v0, v1, gamma,
    alpha = P0(POW_ALPHA);
  model *next = cov->sub[0];

  Abl2(x, next, v);
  if (alpha == 1.0) return;
  Abl1(x, next, &D);
  if (isnowShape(cov)) {
    COV(x, next, &v0);
    *v *= alpha *  POW(v0, alpha - 2.0) * ((alpha - 1.0) * D + v0 * (*v));
  } else {
    COV(ZERO(next), next, &v0);
    COV(x, next, &v1);
    gamma =  v0 - v1;
    *v *= - alpha *  POW(gamma, alpha-2.0) * ((alpha - 1.0) * D + gamma * (*v));
 } // Achtung  "-" , da covarianzfunktion angenommen
}

void InversePow(double *x, model *cov, double *v) {
  // invPow only defined for submodel having variance 1 !!!!
  model *next = cov->sub[0];
  double v0,
    alpha = P0(POW_ALPHA);
    
  if (isnowShape(cov)) {
    INVERSE(x, next, v);
    *v = POW(*v, 1.0 / alpha);
  } else {
    COV(ZERO(next), next, &v0);
    v0 = v0 - POW(POW(v0, alpha) - *x, 1 / alpha);
    INVERSE(&v0, next, v);

    /*
      COV(x, next, v);
      double y = v0 - *v;
      if (y < 0.0 || y > 1.0) {
      if (y > -1e-14 && y < 0.0)  y=0.0;
      else if (y < 1.0 + 1e-14)  y=1.0;
      else {
      //PRINTF("covariance value %10e (1 - %10e = %10e) at %10e\n", *v, *v, 1.0-*v, *x);
      ERR0("invPow valid only for non-negative covariance models with variance 1");
      }
      }
      
      *v = 1.0 - POW(y, 1.0 / alpha);
      //  print("%10g %10g %10g\n", y, *v, P0(POW_ALPHA));
      //  print("Inv %10g %10g %10g v=%10g\n", *x,  P0(POW_ALPHA),y, *v);
      */
  }
}

int initPow(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  model *next = cov->sub[0];
  int null = 0;
  double 
    alpha = P0(POW_ALPHA),
    saveconst = RF_NAN,
    eps = 1e-14;
  cov->taylorN = next->taylorN;
  cov->tailN = next->tailN;
  if (!ISNA(alpha)) {
    if (cov->taylorN > 0) {
      if ((null = (int) !isnowShape(cov))) {
	assert(next->taylor[0][TaylorPow] == 0.0);
	assert(next->taylor[null][TaylorConst] > 0);
	cov->taylor[0][TaylorConst] = POW(next->taylor[0][TaylorConst], alpha);
	saveconst = next->taylor[null][TaylorConst];
	next->taylor[null][TaylorConst] *= -1.0;
      }
    }

    if (cov->taylorN > null) {
     cov->taylor[null][TaylorConst]=POW(next->taylor[null][TaylorConst],alpha);
      cov->taylor[null][TaylorPow] = next->taylor[null][TaylorPow] * alpha;
      double factor = alpha * cov->taylor[null][TaylorConst] /
	next->taylor[null][TaylorConst],
	summand = cov->taylor[null][TaylorPow] - next->taylor[null][TaylorPow];
      for (int i=null + 1; i<next->taylorN; i++) {      
	cov->taylor[i][TaylorConst] = factor * next->taylor[i][TaylorConst];	
	cov->taylor[i][TaylorPow] = summand + next->taylor[i][TaylorPow]; 
      }
      double factor2 = 0.5 * alpha * (alpha - 1.0) *
	cov->taylor[null][TaylorConst] /
	(next->taylor[null][TaylorConst] * next->taylor[null][TaylorConst]),
	summand2 = cov->taylor[null][TaylorPow] -
	           2.0 * next->taylor[null][TaylorPow];
      for (int i = null + 1; i < next->taylorN; i++) {
	for (int j=i; i<next->taylorN; i++) {
	  double p = summand2 +
	    next->taylor[i][TaylorPow] + next->taylor[j][TaylorPow];
	  int k = null + 1;
	  while (k < cov->taylorN && p > cov->taylor[k][TaylorPow] - eps) k++;
	  if (k >= cov->taylorN) break;
	  if (p > cov->taylor[k][TaylorPow] + eps) {
	    for (int m=next->taylorN - 2; m >=k; m--) {
	      cov->taylor[m + 1][TaylorConst] = cov->taylor[m][TaylorConst];
	      cov->taylor[m + 1][TaylorPow] = cov->taylor[m][TaylorPow];
	    }
	    cov->taylor[k][TaylorPow] = p;
	    cov->taylor[k][TaylorConst] = 0.0;
	  }
	  cov->taylor[k][TaylorConst] += (i==j ? 1.0 : 2.0) * factor2 *
	    next->taylor[i][TaylorConst] * next->taylor[j][TaylorConst];
	}
      }
      int i=null + 1;
      double summand3 = cov->taylor[null][TaylorPow] -
	3.0 * next->taylor[null][TaylorPow],
	pmin = summand3 + 3.0 * next->taylor[i][TaylorPow];
      
      while (pmin < cov->taylor[cov->taylorN - 1][TaylorPow] + eps) {
	assert(false); // obige quadratische Approx nicht gut genug
	cov->taylorN--;
      }

      if (!isnowShape(cov)) {
	for (i=null; i < cov->taylorN; i++) cov->taylor[i][TaylorConst] *= -1.0;
      }
    }
    if (null) next->taylor[null][TaylorConst] = saveconst;
          
    if (cov->tailN > 0) {
      assert(cov->tailN == 1);
      cov->tailN = 1; // falls mal tailN groessere Werte annehmen darf
      if (isnowShape(cov)) {
	cov->tail[0][TaylorConst] = POW(next->tail[0][TaylorConst], alpha);
	cov->tail[0][TaylorPow] = next->tail[0][TaylorPow] * alpha;
	cov->tail[0][TaylorExpConst] = next->tail[0][TaylorExpConst] * alpha;
	cov->tail[0][TaylorExpPow] = next->tail[0][TaylorExpPow];
      } else {
	assert(next->taylor[0][TaylorPow] == 0.0); // !! nicht tail
	cov->tail[0][TaylorConst] = next->tail[0][TaylorConst]
	  * alpha * POW(next->taylor[0][TaylorConst], alpha - 1);//!! nicht tail
	cov->tail[0][TaylorPow] = next->tail[0][TaylorPow];
	cov->tail[0][TaylorExpConst] = next->tail[0][TaylorExpConst];
	cov->tail[0][TaylorExpPow] = next->tail[0][TaylorExpPow];	
      }
    }
  }

  RETURN_NOERROR;
}


int checkPow(model *cov) {
  model *next = cov->sub[0];
  int err;
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  if (!isXonly(OWN)) RETURN_ERR(ERRORSTATVARIO);
  //   APMI0(cov);
  if ((err = CHECK_PASSFRAME(next, //EvaluationType
			      cov->frame
			     )) != NOERROR){
    RETURN_ERR(err);  
  }
  //  if ((err = CHECK(next, cov->tsdim, cov->xdimown, cov->typus, //PosDefType,
  //		   OWNDOM(0), OWNISO(0), SCALAR, EvaluationType)) != NOERROR){
  //    // print("OK %d\n",err);
  //    RETURN_ERR(err);  
  //  }
  setbackward(cov, next);      
  assert(VDIM0 == 1);
  cov->mpp.maxheights[0] = RF_NA;
  cov->monotone =
    isMonotone(next->monotone) && P0(POW_ALPHA) > 0 ? MONOTONE : NOT_MONOTONE;
  
  if ((err = initPow(cov, NULL)) != NOERROR) RETURN_ERR(err);
  
  RETURN_NOERROR;
}

void rangePow(model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  if (isnowVariogram(cov)) {
    range->min[POW_ALPHA] = 0.0;
    range->max[POW_ALPHA] = 1.0;
    range->pmin[POW_ALPHA] = 0.01;
    range->pmax[POW_ALPHA] = 1.0;
    range->openmin[POW_ALPHA] = true;
    range->openmax[POW_ALPHA] = false;
  } else {
    range->min[POW_ALPHA] = RF_NEGINF;
    range->max[POW_ALPHA] = RF_INF;
    range->pmin[POW_ALPHA] = - 10;
    range->pmax[POW_ALPHA] = 10;
    range->openmin[POW_ALPHA] = true;
    range->openmax[POW_ALPHA] = true;
  }
}





////////////////////////////////////////////////////////////////////


/* qam */
void kappaqam(int i, model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i < DefList[COVNR].kappas  ? cov->nsub-1 : -1;
}


void qam(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  int i, 
    nsub = cov->nsub;
  double sum, s, w,
    *theta = P(QAM_THETA);
  
  sum = 0.0;
  for (i=1; i<nsub; i++) {
    model * sub = cov->sub[i];
    COV(x, sub, &s);
    INVERSE(&s, next, &w);
    sum += theta[i - 1] * w * w;
  }

  sum = SQRT(sum);
  COV(&sum, next, v); 
}

int checkqam(model *cov) {
  model *next = cov->sub[0],
    *sub;
  int i, err, 
    nsub = cov->nsub,
    nsubM1 = nsub - 1;
  double v, sum;
      
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  
  //  cov->monotone = NORMAL_MIXTURE;
  sum = 0.0;
  for (i=0; i<nsubM1; i++) {
    sum += P(QAM_THETA)[i];
  }
  if (FABS(sum - 1.0) > 1e-14) SERR("theta must sum up to 1");    

  if ((err = CHECK(next, 1,  1, PosDefType, OWNDOM(0), OWNISO(0), 
		     SCALAR, EvaluationType)) != NOERROR)
    RETURN_ERR(err);
  if (!isNormalMixture(next->monotone))
    SERR("phi is not a normal mixture");   
 
  for (i=1; i<nsub; i++) {
    sub = cov->sub[i];
    if ((err = CHECK_PASSFRAME(sub, EvaluationType)) != NOERROR) RETURN_ERR(err);
    //    if ((err = CHECK(sub, cov->tsdim, cov->tsdim, PosDefType,
    //		       OWNDOM(0), OWNISO(0),  SCALAR, EvaluationType)) != NOERROR)
    //RETURN_ERR(err);
    COV(ZERO(sub), sub, &v);
    if (v != 1.0) SERR("unit variance required");
    setbackward(cov, sub);
  } 

  INVERSE(ZERO(next), next, &v);

  if (ISNAN(v)) SERR1("inverse function of '%.50s' unknown", NICK(next));

  cov->logspeed = 0.0;   
  RETURN_NOERROR;
}


void rangeqam(model *cov, range_type *range){
  double pmax = 2.0 / (double)  (cov->nsub-1);
  range->min[QAM_THETA] = 0.0;
  range->max[QAM_THETA] = 1.0;
  range->pmin[QAM_THETA] = 0.0;
  range->pmax[QAM_THETA] = pmax;
  range->openmin[QAM_THETA] = false;
  range->openmax[QAM_THETA] = false;
}



/* mqam */
void kappamqam(int i, model *cov, int *nr, int *nc) {
  int nsub = cov->nsub - 1;
  *nc = (i==QAM_THETA) ? 1 : -1;
  *nr = i == QAM_THETA ? nsub : -1;
}
void mqam(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  int i, j, k, l,
    vdim = VDIM0,
    vdimP1 = vdim + 1;
  double s0,
    *theta = P(QAM_THETA),
    s[MAXSUB];
  
  for (i=0; i<vdim; i++) {
    model * sub = cov->sub[i+1];
    COV(x, sub, &s0);
    INVERSE(&s0, next, s + i);
    s[i] *= theta[i] * s[i];
  }
    
  for (j=0; j<vdim; j++) {
    l = k = j * vdimP1; // nicht vdim
    for (i=j; i<vdim; i++, k++, l+=vdim) {
      s0 = SQRT(s[i] + s[j]);
      COV(&s0, next, v + k);
      v[l] = v[k]; 
      // if (k != l) v[k] *= rho[k];
      // assert (rho[l] == rho[k]) ;
            
    }
  }

}

int checkmqam(model *cov) {
  int err,
    nsub = cov->nsub,
    vdim = nsub - 1;     

  //  NotProgrammedYet("");

  if ((err = checkqam(cov)) != NOERROR) RETURN_ERR(err);
  VDIM0 = VDIM1 = vdim;
  RETURN_NOERROR;
}

void rangemqam(model *cov, range_type *range){
  //double pmax = 2.0 / (double)  (cov->nsub-1);

  rangeqam(cov, range);
}



/////////////// NATSC

void natsc(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  double invscale, y;

  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = *x * invscale;
  // letzteres darf nur passieren wenn dim = 1!!
  COV(&y, next, v);
}

void Dnatsc(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  double invscale, y;
 
  assert(DefList[NEXTNR].inverse != NULL);
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = *x * invscale;

  Abl1(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= invscale; 
}

void DDnatsc(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  double invscale, y, invScSq;

  assert(DefList[NEXTNR].inverse != NULL);
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = *x * invscale;
  invScSq = invscale * invscale;

  Abl2(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= invScSq; 
}

void Inversenatsc(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  double invscale, modelinv;
  assert(DefList[NEXTNR].inverse != NULL);
  INVERSE(x, next, &modelinv);    
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  *v = modelinv / invscale;
}

int checknatsc(model *cov) {
  model *next = cov->sub[0];
  int err;


  assert(isNatsc(cov));

   if ((err = CHECK_PASSFRAME(next, EvaluationType)) != NOERROR) RETURN_ERR(err);
   // if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, OWNDOM(0),
   //		   OWNISO(0), SUBMODEL_DEP, EvaluationType))
   // != NOERROR) {
   //    RETURN_ERR(err);
   //  }

  if (DefList[NEXTNR].inverse == NULL)
    SERR1("natural scaling is not defined for %.50s", NICK(next));
 
  double invscale;
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);

  if (ISNAN(invscale))
    SERR1("inverse function of '%.50s' unknown", NICK(next));


  cov->logspeed = 0.0;
  setbackward(cov, next);

  VDIM0 = next->vdim[0];
  VDIM1 = next->vdim[1];

  RETURN_NOERROR;
}


int initnatsc(model *cov, gen_storage *s){
  //  location_type *loc = Loc(cov);

  if (hasGaussMethodFrame(cov)) {
    model *next = cov->sub[0];
    return INIT(next, cov->mpp.moments, s);
  }

  else if (hasMaxStableFrame(cov) || hasAnyPoissonFrame(cov)) {
   
    SERR("natsc for max-stable processes and poisson process not programmed yet");

  }

  else ILLEGAL_FRAME;

  RETURN_NOERROR;
}

void donatsc(model *cov, gen_storage *s){
  model *next = cov->sub[0];
  DO(next, s);
}

void spectralnatsc(model *cov, gen_storage *S, double *e) {
  //  spectral_storage *s = &(S->Sspectral);
  model *next = cov->sub[0];
  int d,
    dim = OWNLOGDIM(0);
  double invscale;

  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);  
  SPECTRAL(next, S, e);
  for (d=0; d<dim; d++) e[d] *=  invscale;
}

void coinitnatsc(model *cov, localinfotype *li) {
  model *next = cov->sub[0];
  defn *C = DefList + NEXTNR; // not gatternr
  if ( C->coinit == NULL) 
    ERR0("# cannot find coinit -- please inform author");
  C->coinit(next, li); // not gatternr
}

void ieinitnatsc(model *cov, localinfotype *li) {
  model *next = cov->sub[0];
  defn *C = DefList + NEXTNR; // not gatternr
  if ( C->ieinit == NULL) // not gatternr
    ERR0("# cannot find ieinit -- please inform author");
  C->ieinit(next, li); // not gatternr
}

void tbm2natsc(double *x, model *cov, double *v){
  model *next = cov->sub[0];
  double invscale, y;

  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = x[0] * invscale;
     
  TBM2CALL(&y, next, v)
}

//////////////////////////////////////





void mult_inverse(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  FCTN(x, next, v);
  *v = 1.0 / *v;
}


void mult_inverseNonstat(double *x, double *y, model *cov, double *v) {
  model *next = cov->sub[0];
  NONSTATCOV(x, y, next, v);
  *v = 1.0 / *v;
}


int checkmult_inverse(model *cov) {
  model
    *next = cov->sub[0];
  assert(VDIM0 == 1 && VDIM1 == 1 );
   int err = NOERROR;
  if ((err = CHECK_PASSTF(next,  ShapeType, SUBMODEL_DEP, cov->frame)) != NOERROR)
    RETURN_ERR(err);
  //  if ((err = CHECK(next, cov->tsdim,  cov->xdimprev, ShapeType,
  //		   OWNDOM(0), OWNISO(0),
  //		   SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
  setbackward(cov, next);
  cov->mpp.maxheights[0] = RF_NA;
  RETURN_NOERROR;
}




////////////////////

void addSetParam(model **newmodel, model* remote, 
		 param_set_fct set, bool performdo, int variant, int nr) {
  set_storage *S;
  assert(SETPARAM_LOCAL == 0);
  addModel(newmodel, nr, remote);
  kdefault(*newmodel, SET_PERFORMDO, performdo);
  NEW_COV_STORAGE(*newmodel, set);
  S = (*newmodel)->Sset;
  // S->from wird unten gesetzt !
  S->remote = remote;
  S->set = set;  
  S->variant = variant;   
}   
void addSetParam(model **newmodel, model * remote, 
		 param_set_fct set,  bool performdo, 
		 int variant) {
  addSetParam(newmodel, remote, set, performdo, variant, SETPARAM);
}
void addSetDistr(model **newmodel, model * remote, 
		 param_set_fct set,  bool performdo, int variant) {
  addSetParam(newmodel, remote, set, performdo, variant, SET_DISTR);
}


void setparamStat(double *x, model *cov, double *v){
  COV(x, cov->sub[SETPARAM_LOCAL], v);
}

void setparamNonStat(double *x, double *y, model *cov, double *v){
  NONSTATCOV(x, y, cov->sub[SETPARAM_LOCAL], v);
}

void Dsetparam(double *x, model *cov, double *v){
  Abl1(x, cov->sub[SETPARAM_LOCAL], v);
}

void DDsetparam(double *x, model *cov, double *v){
  Abl2(x, cov->sub[SETPARAM_LOCAL], v);
}

void D3setparam(double *x, model *cov, double *v){
  Abl3(x, cov->sub[SETPARAM_LOCAL], v);
}
void D4setparam(double *x, model *cov, double *v){
  Abl4(x, cov->sub[SETPARAM_LOCAL], v);
}

void Inverse_setparam(double *v, model *cov, double *x){
  model *next = cov->sub[SETPARAM_LOCAL];
  INVERSE(v, next, x);
}

void NonstatInverse_setparam(double *v, model *cov, double *x, double *y){
  model *next = cov->sub[SETPARAM_LOCAL];
  NONSTATINVERSE(v, next, x, y);
}
void LogNonstatInverse_setparam(double *v, model *cov, double *x, double *y){
  model *next = cov->sub[SETPARAM_LOCAL];
  NONSTATLOGINVERSE(v, next, x, y);
}

int checksetparam(model *cov) {
  model *next = cov->sub[SETPARAM_LOCAL];
  int err;
  Types type = OWNTYPE(0);
 
  kdefault(cov, SET_PERFORMDO, true);
  ASSERT_ONESYSTEM;
  if (equalsRandom(type) || equalsRandom(SYSTYPE(NEXT, 0))) {
    BUG;
  }
   if ((err = CHECK_PASSTYPE(next, OWNTYPE(0))) != NOERROR) RETURN_ERR(err);
   // if ((err = CHECK(next, dim, xdim, type, dom, iso, 
   //		   SUBMODEL_DEP, frame)) != NOERROR) RETURN_ERR(err);
  setbackward(cov, next);
  VDIM0 = next->vdim[0];
  VDIM1 = next->vdim[1];
  cov->randomkappa = true;

  TaylorCopy(cov, next);
  
  // ACHTUNG ! weder SETPARAM_FROM (da nur link) noch SETPARAM_SET (da 
  //           i.a. keine echten Modelle) werden ueberprueft!
 
  RETURN_NOERROR;
}


Types Typesetparam(Types required, model *cov, isotropy_type i) {
  return TypeConsistency(required, cov->sub[SETPARAM_LOCAL], i);
}

void spectralsetparam(model *cov, gen_storage *s, double *e){
   SPECTRAL(cov->sub[SETPARAM_LOCAL], s, e);  // nicht gatternr
}

int initsetparam(model *cov, gen_storage *s){
  model *next= cov->sub[SETPARAM_LOCAL];
  set_storage *X = cov->Sset;

  assert(X != NULL);
  int err, i,
    vdim = VDIM0;
  if (VDIM0 != VDIM1) BUG;

  if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR)
    RETURN_ERR(err);
  if (X->remote != NULL) {
    assert(X->set != NULL);
    X->set(cov->sub[0], X->remote, X->variant);
  }

  TaylorCopy(cov, next);
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = next->mpp.maxheights[i];
  RETURN_NOERROR;
}


void dosetparam(model *cov, gen_storage *s) {
  bool performDo = P0INT(SET_PERFORMDO);
  if (performDo) DO(cov->sub[SETPARAM_LOCAL], s); 
}

void covmatrix_setparam(model *cov, double *v) {
  model *next = cov->sub[SETPARAM_LOCAL];
  DefList[NEXTNR].covmatrix(next, v);
}

char iscovmatrix_setparam(model *cov) {
  model *next = cov->sub[SETPARAM_LOCAL];
  return DefList[NEXTNR].is_covmatrix(next);
}

void range_setparam(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  booleanRange(SET_PERFORMDO);
}



void kappatrafo(int i, model  VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  *nc = 1;
  *nr =  (i == TRAFO_ISO) ? SIZE_NOT_DETERMINED : -1;
}

void trafo(double *x, model *cov, double *v){  
  model *next = cov->sub[0];
  assert(next != NULL);
  COV(x, next, v);
}
void logtrafo(double *x, model *cov, double *v, 
		     double *Sign){  
  model *next = cov->sub[0];
  assert(next != NULL);
  LOGCOV(x, next, v, Sign);
  assert(P0INT(TRAFO_ISO) == NEXTISO(0)); 
}
void nonstattrafo(double *x, double *y, model *cov, double *v){  
  model *next = cov->sub[0];
  assert(next != NULL);
  NONSTATCOV(x, y, next, v);
}
void lognonstattrafo(double *x, double *y, model *cov, double *v, 
		     double *Sign){  
  model *next = cov->sub[0];
  assert(P0INT(TRAFO_ISO) == ISO(PREVSYSOF(next), 0)); 
  assert(next != NULL);
  LOGNONSTATCOV(x, y, next, v, Sign);
}


int checktrafo(model *cov){
  ASSERT_ONESYSTEM;
  if (PisNULL(TRAFO_ISO)) SERR("parameter not given");
  if (cov->nsub == 0) addModel(cov, 0, IDCOORD);
  assert(cov->nsub == 1);
  model *next = cov->sub[0];
  int err = NOERROR;
  isotropy_type iso = (isotropy_type) P0INT(TRAFO_ISO);
  //PMI(cov);
  if (isAnyIsotropic(iso)) set_xdim(OWN, 0, 1);
  else set_xdim(OWN, 0, isSpaceIsotropic(iso) ? 2 : GATTERXDIM(0));
  set_logdim(OWN, 0, GATTERLOGDIM(0));

  isotropy_type newiso = OWNISO(0),
    previso = PREVISO(0);
  
  if ((equalsCoordinateSystem(newiso) || equalsAnySymmetric(newiso) ||
       isEarthProjection(newiso)) &&
      newiso != CoordinateSystemOf(previso)) {
    if (!isCartesian(newiso)) SERR("Only transformations from earth systems to cartesian systems are currently programmed.");
    if (isAnyIsotropic(previso)) newiso=ISOTROPIC;
    else if (equalsEarthSymmetric(previso) || equalsSphericalSymmetric(previso))
      newiso = SYMMETRIC;
    set_iso(OWN, 0, newiso);
  }

  if (next == NULL) {
    addModel(cov, 0, IDCOORD);
    next = cov->sub[0];
  }
  if ((err = CHECK_NOPASS(next)) != NOERROR) RETURN_ERR(err);  

  setbackward(cov, next);
  // no setbackward ?!!
  
  if (VDIM0 == SUBMODEL_DEP || VDIM0 == PARAM_DEP) {
    VDIM0 = next->vdim[0];
    VDIM1 = next->vdim[1];
  } else {
    if (VDIM0 != next->vdim[0] || VDIM1 != next->vdim[1]) {
      PMI(cov); //
      BUG;
    }
  }
  
  RETURN_NOERROR;
}

bool settrafo(model *cov) {
  isotropy_type iso = CONDPREVISO(0); 
  if (!isFixed(iso)) return false;
  
  if (PisNULL(TRAFO_ISO))
    ERR2("argument '%.50s' in '%.50s' not given.", KNAME(TRAFO_ISO), NICK(cov));
  isotropy_type newiso = (isotropy_type) P0INT(TRAFO_ISO);
  if (PREV_INITIALISED) {
    if (equalsCoordinateSystem(newiso) && equalsAnySymmetric(PREVISO(0)))
      newiso = SymmetricOf(newiso);
  }
  if (cov->nsub == 0 || MODELNR(cov->sub[0]) == IDCOORD) {
    VDIM0 = SUBMODEL_DEP; // OWNTOTALXDIM ist noch nicht definiert bei Aufruf
    VDIM1 = 1;
  } else {
    VDIM1 = VDIM0 = SUBMODEL_DEP;
  }

  set_iso(OWN, 0, newiso);
  set_dom(OWN, 0, isAnyIsotropic(newiso) ? XONLY : PREVDOM(0));
  
  return true;
}

bool allowedDtrafo(model *cov) {
  model *next = cov->sub[0];
  bool *D = cov->allowedD;
  isotropy_type previso = CONDPREVISO(0);      
  if (equalsIsoMismatch(previso)) {
    cov->IallowedDone = false;
    return allowedDfalse(cov);
  }
  if (isNegDef(PREVTYPE(0)) &&
      (EssentialCoordinateSystemOf((isotropy_type) P0INT(TRAFO_ISO)) !=
       EssentialCoordinateSystemOf(previso))){
    for (int i=FIRST_DOMAIN; i<=LAST_DOMAINUSER; D[i++]=false);
    D[KERNEL] = true;
    return false;
  } // else domain is given by subsequent model, see lines below

  if (next == NULL) {
    for (int i=FIRST_DOMAIN; i<=LAST_DOMAINUSER; D[i++]=false);
    D[XONLY] = true;
  } else return allowedDstandard(cov);
  return false;
}

bool allowedItrafo(model *cov) {
  bool *I = cov->allowedI;
  isotropy_type iso = (isotropy_type) P0INT(TRAFO_ISO);
  MEMSET(I, 0, sizeof(allowedI_type)); // default is false
  if (isCartesian(iso)) {
    for (int i=(int) iso; i<=LAST_CARTESIAN; I[i++] = true);    
    if (equalsVectorIsotropic(iso)) {
      I[SYMMETRIC] = false;
    } else {
      if (isEarthProjection(iso)) 
	for (int i=FIRST_PROJECTION; i<=LAST_PROJECTION; i++) I[i] = i == iso;
      I[EARTH_SYMMETRIC] = isSymmetric(iso);
    }
    I[EARTH_COORD] = true;
  } else if (isSpherical(iso)) {
    for (int i=iso; i<=LAST_TRUESPHERICAL; i++)
      I[i] = I[i + FIRST_EARTH - FIRST_SPHERICAL] = true;
  } else if (isEarth(iso)) {
    for (int i=(int) iso; i<=LAST_EARTH; I[i++] = true);    
  } else BUG;

  I[CYLINDER_COORD] = false; // not programmed yet
  I[UNREDUCED] = true;

  return false;
}

Types Typetrafo(Types required, model *cov, isotropy_type i) {
  if (cov->sub[0] == NULL) return required==ShapeType ? required : BadType;
  return TypeConsistency(required, cov->sub[0], i);
}

void rangetrafo(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[TRAFO_ISO] = FIRST_ISOUSER;
  range->max[TRAFO_ISO] = LAST_ISOUSER;
  range->pmin[TRAFO_ISO] = FIRST_ISOUSER;
  range->pmax[TRAFO_ISO] = LAST_ISOUSER;
  range->openmin[TRAFO_ISO] = false;
  range->openmax[TRAFO_ISO] = false;
}


int checktrafoproc(model *cov) { // auch fuer Trendproc
  model *next = cov->sub[0],
    *key = cov->key;  
  location_type *loc = Loc(cov);
  int err;  
  FRAME_ASSERT_GAUSS;
  ASSERT_ONESYSTEM;
  if (PisNULL(TRAFO_ISO)) SERR("parameter not given");
  assert(cov->nsub == 1);

  
  if (cov->key != NULL) {
    if ((err = CHECK(key, 3 + loc->Time, 3+ loc->Time, ProcessType,
		     XONLY, CARTESIAN_COORD,
		     SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
  } else {
    if ((err = CHECK(next, OWNLOGDIM(0), OWNXDIM(0), OWNTYPE(0), OWNDOM(0), 
		     (isotropy_type) P0INT(TRAFO_ISO), -1, EvaluationType))
	!= NOERROR){
      RETURN_ERR(err);  
    }
    if (!isnowVariogram(cov)) SERR("definite function needed");
  }

  // no setbackward??!!
  
  VDIM0 = next->vdim[0];
  VDIM1 = next->vdim[1];
  
  RETURN_NOERROR;
}


int structtrafoproc(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  model *sub = cov->sub[0];
  int err,
    n = OWNSYSTEMS;
  location_type *loc = Loc(cov);
  bool 
    Time = GetTime(cov); // muss ganz oben stehen
  FRAME_ASSERT_GAUSS;

   if (n > 2 || (n == 2 && !Time) ||
      !equalsCartCoord((isotropy_type) P0INT(TRAFO_ISO)) ||
      cov->calling == NULL || !equalsEarthCoord(ISO(CALLING, 0)))
     SERR("currently only earth-to-cartesian allowed");

  if (cov->key != NULL) BUG;
 
  // ACHTUNG In der Reihenfolge, da (i) Loc einzuverwendet ist
  // (ii) umgehaengt werden muss da TransformLoc ownloc verwendet, ownloc
  // aber durch die Erdtrafo ueberschreiben wird (siehe weiter unten).
  // Somit muss ownloc gesichert werden. Wird eingehaengt, so dass bei Fehlern
  // ownloc dennoch zum Schluss geloescht wird.

  TransformLoc(cov, true, True, false);
  // SetLoc2NewLoc(sub, PLoc(cov));  ist hier zu frueh, da Bindung
  // durch loc_set unten zerstoert wird
  
  loc = Loc(cov);
  assert(!loc->grid);
  if (loc->len != 1)
    SERR("trafo currently only possible for a single data set");

  int 
    spatialdim = loc->spatialdim,
    newdim =  spatialdim < 2 ? NA_INTEGER : spatialdim > 3 ? spatialdim -1 : 3,
    spatialpts = loc->spatialtotalpoints;
  double aequ, pol, T[3],
    *y = (double *) MALLOC(sizeof(double) * newdim * spatialpts);
  if (y == NULL) RETURN_ERR(ERRORMEMORYALLOCATION);

  if (Time) MEMCOPY(T, loc->T, sizeof(double) * 3);
  if (STRCMP(GLOBAL.coords.newunits[0], UNITS_NAMES[units_km]) == 0) {
    aequ = 6378.1;
    pol =  6356.8;
  } else {
    aequ = 3963.17;
    pol = 3949.93;
  }

  Earth2Cart(cov, aequ, pol, y);       
  loc_set(y, NULL, T, newdim, newdim, spatialpts, 0, Time, false, false, cov);
  SetLoc2NewLoc(sub, PLoc(cov));  

  if ((err = covcpy(&(cov->key), sub)) != NOERROR) {
    if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
    goto ErrorHandling;
  }

  addModel(&(cov->key), GAUSSPROC);

#if MAXSYSTEMS > 1
  BUG; logdim und xdim passen nicht
#endif
 
  if ((err = CHECK(cov->key, 3 + Time, 3 + Time, ProcessType, 
		   XONLY, CARTESIAN_COORD,
		   SUBMODEL_DEP, cov->frame)) != NOERROR) goto ErrorHandling;

  if ((err = STRUCT(cov->key, NULL)) != NOERROR) goto ErrorHandling;
 
 ErrorHandling :
  FREE(y);

  RETURN_ERR(err);
}


int inittrafoproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){// auch fuer Trendproc
  model *key = cov->key;
  int err;

  if (VDIM0 != 1) NotProgrammedYet("");
  if ((err = INIT(key, 0, s)) != NOERROR) RETURN_ERR(err);

  ReturnOtherField(cov, key);
  cov->simu.active = true;
  RETURN_NOERROR;
}


void dotrafoproc(model *cov, gen_storage *s){
  model *key = cov->key;
  DO(key, s);   
}




void nonstatprod(double *x, double *y, model *cov, double *v){  
  model *next = cov->sub[0];
  int
    rows = next->vdim[0],
    cols = next->vdim[1],
    vdimsq = rows * cols;
  
  TALLOC_XX1(w, vdimsq);
  FCTN(y, next, w);

  if (vdimsq == 1) {
    FCTN(x, next, v);
    *v *= *w;
  } else {
    TALLOC_XX2(u, vdimsq);
    FCTN(x, next, u); 
    matmulttransposed(u, w, v, cols, rows, cols);
    END_TALLOC_XX2;
  }
  END_TALLOC_XX1;
}

int checkprod(model *cov) {
  //  location_type *loc = Loc(cov);
  if (cov->sub[0]  == NULL) {
    addModel(cov, 0, IDCOORD);
  }
  model *next = cov->sub[0];
  int err;  
  if ((err = CHECK(next, OWNLOGDIM(0), OWNXDIM(0), ShapeType, XONLY, OWNISO(0),
		   SUBMODEL_DEP, EvaluationType)) !=NOERROR) RETURN_ERR(err);

  assert((OWNDOM(0) == KERNEL) xor (COVNR == PROD_PROC));
  
  setbackward(cov, next);
  VDIM0 = next->vdim[0]; 
  VDIM1 = next->vdim[1];
  cov->pref[Direct] = 1;
  cov->pref[Specific] = cov->pref[Nothing]= 5;
  
  EXTRA_STORAGE;
  RETURN_NOERROR;  
}


int checkprodproc(model *cov) { // auch fuer Trendproc
  int err;  
  if ((err = checkprod(cov))) RETURN_ERR(err);
  if (VDIM0 != 1) NotProgrammedYet("");
  if (cov->q == NULL) {QALLOC(1); cov->q[PRODPROC_RANDOM] = (double) true;}
  RETURN_NOERROR;
}


int structprodproc(model  VARIABLE_IS_NOT_USED *cov,
		   model VARIABLE_IS_NOT_USED **newmodel){// auch fuer Trendproc
   RETURN_NOERROR;
}


int initprodproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){// auch fuer Trendproc
  int err;
  if (VDIM0 != 1) NotProgrammedYet("");

  if ((err = check_fctn(cov)) != NOERROR) RETURN_ERR(err);
  
  FRAME_ASSERT_GAUSS;
  err = ReturnOwnField(cov);
  cov->simu.active = err == NOERROR;
  if (PL>= PL_STRUCTURE) {PRINTF("\n'%.50s' is now initialized.\n", NAME(cov));}
  RETURN_ERR(err);
}


void doprodproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  location_type *loc = Loc(cov);					
  long i,
    vdim = VDIM0,
    totvdim = loc->totalpoints * vdim;
  double *res = cov->rf;
  assert(res != NULL);

  Fctn(NULL, cov, res);

  if (cov->q[PRODPROC_RANDOM]) {
    double g = GAUSS_RANDOM(1.0);  
    for(i=0; i<totvdim; i++) res[i] *= g;
  }
}


void nonstatsum(double *x, double *y, model *cov, double *v){  
  model *next = cov->sub[0];
  int 
    rows = next->vdim[0],
    cols = next->vdim[1],
    vdimsq = rows * cols;
  
  TALLOC_XX1(w, vdimsq);
  FCTN(y, next, w);
  FCTN(x, next, v);
  for (int i=0; i<vdimsq; i++) v[i] += w[i];
  END_TALLOC_XX1;
}

int checksum(model *cov) {
  //  location_type *loc = Loc(cov);
  if (cov->sub[0] == NULL) {
    addModel(cov, 0, IDCOORD);
  }
  model *next = cov->sub[0];
  int err;

  if ((err = CHECK(next, OWNLOGDIM(0), OWNXDIM(0),  ShapeType,
		     XONLY, OWNISO(0),
		     SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
  
  setbackward(cov, next);
  if (VDIM0 != VDIM1) 
    SERR("sub model must return symmetric a square matrix");
  EXTRA_STORAGE;
  RETURN_NOERROR;  
}

//  LocalWords:  NOERROR



/* g(x)+g(y)-g(x,y) - g(0,0) */
void kappavariogram2cov(int i, model VARIABLE_IS_NOT_USED *cov,
			int *nr, int *nc){
  *nc = *nr = i <= VAR2COV_C ? SIZE_NOT_DETERMINED : -1;
}


void variogram2cov(double *x, double *y, model *cov, double *v){
  model *next = cov->sub[0];
  int vdim = VDIM0,
    vdimSq = vdim * vdim,
    dim = OWNLOGDIM(0),
    n = NROW(VAR2COV_C);
  assert(cov->Scovariate != NULL && cov->Scovariate->x);
  double 
    *c = P(VAR2COV_C),
    *z = cov->Scovariate->x,
    *zbase = z;

  TALLOC_X1(resx, vdimSq);
  TALLOC_X2(resy, vdimSq);
  NONSTATCOV(x, y, next, v);

  // printf("n=%d %d %d %10g %10g\n", n, dim, cov->Scovariate->pts, z[0], z[1]); assert(n == 2);

  for (int k=0; k<n; k++, z += dim) {
    double ck = c[k];
    NONSTATCOV(x, z, next, resx);
    NONSTATCOV(z, y, next, resy);
    //  printf("k=%d x=%2.1f %2.1f y=%2.1f %2.1f ", k, *x, *resx, *y, *resy);
    for (int d=0; d<vdimSq; d++) v[d] -= ck * (resx[d] + resy[d]);

    for (int m=0; m<n; m++) {
      double ckm = ck * c[m];      
      NONSTATCOV(z, zbase + m * dim, next, resx);
      for (int d=0; d<vdimSq; d++) {
	v[d] += ckm * resx[d];
	//	printf("m=%d z=%2.1f %2.1f %2.1f", m, *z, zbase[m*dim], *resx);
      }
    }
    // printf("v=%2.3f\n", *v);
  }
  END_TALLOC_X1;
  END_TALLOC_X2;
}


int checkvariogram2cov(model *cov) {  
  model *next = cov->sub[0];
  int err,
    dim = ANYOWNDIM; // ANOWNTOTALDIM
  defn *C = DefList + COVNR;
  SEXP xlist = PVEC(VAR2COV_X);

  if (cov->Scovariate == NULL) {
    NEW_STORAGE(covariate);
    covariate_storage *s = cov->Scovariate;
    if (length(xlist) == 1) {      
      int method = (int) REAL(VECTOR_ELT(xlist, 0))[0] - 1;
      location_type *loc = Loc(cov);
      switch(method) {
      case 0: case 1: { // origin or center 
	if ((s->x = (double*) CALLOC(dim, sizeof(double))) == NULL)
	  RETURN_ERR(ERRORMEMORYALLOCATION); 
	double *xcond = s->x;
	if (method == 1) {
	  if (Getgrid(cov)) {
	    assert(loc != NULL);
	    assert(loc->xgr != NULL);
	    for (int d=0; d<dim; d++) 
	      xcond[d] = loc->xgr[d][XSTART] +
		loc->xgr[d][XSTEP] * 0.5 * (loc->xgr[d][XLENGTH] - 1.0);
	  } else {
	    int spatial = Getspatialpoints(cov),
	      spatial_dim =  GetLocspatialdim(cov);
	    double *x = loc->x,
	      dummy = (double) spatial;
	    for (int i=0; i<spatial; i++, x+=spatial_dim) {
	      for (int d=0; d<spatial_dim; d++) xcond[d] += x[d];
	    }
	    for (int d=0; d<spatial_dim; d++) xcond[d] /= dummy;
	    if (GetTime(cov)) {
	      double *T = loc->T;
	      xcond[spatial_dim] = T[XSTART] + T[XSTEP] * 0.5*(T[XLENGTH]-1.0);
	    }
	  }
	}
	s->pts = 1;
      }
	break;
      case 2 : // extremals
	if (Getgrid(cov)) {
	  int pts = s->pts = 1 << dim;
	  if ((s->x = (double*) CALLOC(pts * dim, sizeof(double))) ==NULL)
	    RETURN_ERR(ERRORMEMORYALLOCATION); 
	  double *xcond = s->x,
	    *factor = (double*) CALLOC(dim, sizeof(double));
	  for (int i=0; i<pts; i++, xcond+=dim) {
	    for (int d=0; d<dim; d++) {
	      xcond[d] = loc->xgr[d][XSTART] + factor[d] * loc->xgr[d][XSTEP];
	      //printf("i=%d %10g %10g\n", i, *xcond, *factor);
	    }
	    for (int d=0; d<dim; d++) {
	      if (factor[d] == 0.0) {
		factor[d] = loc->xgr[d][XLENGTH] - 1.0;
		break;
	      } else {
		factor[d] = 0.0;
	      }
	    }
	  }
	  //	  BUG;
	  FREE(factor);
	} else {
	  covariate_DELETE(&(cov->Scovariate));
	  SERR2("currently only grids are allowed for option '%.50s' of '%.50s'.\n",
		RMCOV_X[method], C->kappanames[VAR2COV_X]);
	}
   	break;
      case 3: // all
	if (Getgrid(cov)) {	  
	  int pts = s->pts = Gettotalpoints(cov);
	  if ((s->x = (double*) CALLOC(pts * dim, sizeof(double))) ==NULL)
	    RETURN_ERR(ERRORMEMORYALLOCATION);
	  double *xcond = s->x,
	    *factor = (double*) CALLOC(dim, sizeof(double));

	  for (int i=0; i<pts; i++, xcond+=dim) {
	    for (int d=0; d<dim; d++) {
	      xcond[d] = loc->xgr[d][XSTART] + factor[d] * loc->xgr[d][XSTEP];
	    }
	    for (int d=0; d<dim; d++) {
	      factor[d]++;
	      if (factor[d] < loc->xgr[d][XLENGTH]) break;
	      else factor[d] = 0.0;
	    }
	  }
	  //
	  FREE(factor);
	} else {
	  covariate_DELETE(&(cov->Scovariate));
	  SERR2("currently only grids are allowed for option '%.50s' of '%.50s'.\n",
		RMCOV_X[method], C->kappanames[VAR2COV_X]);
	}
	break;
      default: BUG;
      }
    } else { // x, y, z, T given explicitely
      int store = GLOBAL.general.set;
      location_type *loc = LocLoc(s->loc);
      GLOBAL.general.set = 0;
      s->loc = loc_set(PVEC(VAR2COV_X), true);
      assert(s->loc != NULL);
      GLOBAL.general.set = store;
      if (loc->timespacedim != GetLoctsdim(cov)) {
	covariate_DELETE(&(cov->Scovariate));
	ERR0("dimension of the conditioning coordinates does not match dimension of the field");
      }
      GLOBAL.general.set = 0;
      TransformLoc(cov, loc, &(s->x));
      GLOBAL.general.set = store;            
      if (loc->timespacedim != OWNLOGDIM(0))  {
	covariate_DELETE(&(cov->Scovariate));
	SERR1("number of columns of '%.50s' do not equal the dimension",
	      KNAME(VAR2COV_X));
      }
      GLOBAL.general.set = 0;
      s->pts = LocLoc(s->loc)->totalpoints;
    }
  } // Scovariate == NULL


  int pts = cov->Scovariate->pts;
  if (pts == 0) BUG;
  if (PisNULL(VAR2COV_C)) {
    PALLOC(VAR2COV_C, pts, 1); // automatisch 0
    double *p = P(VAR2COV_C),
      ninv = 1.0 / (double) pts;    
    for (int i=0; i<pts; p[i++] = ninv);
  } else {
    if (pts != cov->nrow[VAR2COV_C]) SERR1("length of '%.50s' does not match the number of given locations", C->kappanames[VAR2COV_C]);
    double
      sum = 0.0,
      *c = P(VAR2COV_C);
    for (int i=0; i<pts; i++) {
      if (ISNAN(c[i]))
	SERR1("currently, components of '%.50s' cannot be estimated by MLE, so NA's are not allowed.", C->kappanames[VAR2COV_C]);
      sum += c[i];
    }
    if (FABS(sum - 1.0) > 1e-15)
      SERR1("Values of '%.50s' do not sum up to 1.", C->kappanames[VAR2COV_C]);
  }

  if ((err = CHECK_PASSTF(next, VariogramType, VDIM0, EvaluationType))
      != NOERROR) RETURN_ERR(err);
  setbackward(cov, next);

  EXTRA_STORAGE;
  RETURN_NOERROR;
}

void rangevariogram2cov(model VARIABLE_IS_NOT_USED *cov, range_type *range){  
  range->min[VAR2COV_X] = RF_NEGINF;
  range->max[VAR2COV_X] = RF_INF;
  range->pmin[VAR2COV_X] = - 10.0;
  range->pmax[VAR2COV_X] = 10.0;
  range->openmin[VAR2COV_X] = true;
  range->openmax[VAR2COV_X] = true;

  range->min[VAR2COV_C] = RF_NEGINF;
  range->max[VAR2COV_C] = RF_INF;
  range->pmin[VAR2COV_C] = - 10.0;
  range->pmax[VAR2COV_C] = 10.0;
  range->openmin[VAR2COV_C] = true;
  range->openmax[VAR2COV_C] = true;

  range->min[VAR2COV_SPECIAL] = 0;
  range->max[VAR2COV_SPECIAL] = nVAR2COV_METHODS + 1;
  range->pmin[VAR2COV_SPECIAL] =  0;
  range->pmax[VAR2COV_SPECIAL] = range->max[VAR2COV_SPECIAL] ;
  range->openmin[VAR2COV_SPECIAL] = false;
  range->openmax[VAR2COV_SPECIAL] = false;
}

//  LocalWords:  NOERROR




int checkvar2covproc(model *cov) { // auch fuer Trendproc
  int err;
  if (isProcess(cov)) {
    RETURN_ERR(ERRORNOTPROGRAMMEDYET);
  } else {
    if ((err = checkvariogram2cov(cov))) RETURN_ERR(err);
  }
  RETURN_NOERROR;
}


int structvar2covproc(model *cov,
		   model VARIABLE_IS_NOT_USED **newmodel){// auch fuer Trendproc
  model
    *next = cov->sub[0];
  int err;
  if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
  if (!isGaussMethod(cov->key)) addModel(&(cov->key), GAUSSPROC);
  ASSERT_ONESYSTEM;
  if ((err = CHECK_NO_TRAFO(cov->key, OWNLOGDIM(0), OWNXDIM(0),
			    ProcessType, XONLY,
			    OWNISO(0), VDIM0, GaussMethodType)) != NOERROR) {
    RETURN_ERR(err);
  }
  err = STRUCT(cov->key, NULL);
  RETURN_ERR(err);
}


int initvar2covproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){// auch fuer Trendproc
  int err;
  model 
    *key = cov->key;
  assert(key != NULL);
  
  FRAME_ASSERT_GAUSS;
  if ((err = INIT(key, 0, s)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  ReturnOtherField(cov, key);
  cov->simu.active = true;
  if (PL>= PL_STRUCTURE) {PRINTF("\n'%.50s' is now initialized.\n", NAME(cov));}
  RETURN_NOERROR;
}


void dovar2covproc(model *cov, gen_storage *s){
  //location_type *loc = Loc(cov);					
  //ong
  // vdim = VDIM0,
  //totvdim = loc->totalpoints * vdim;
  assert(cov->rf != NULL);
  //  double *res = cov->rf;


  DO(cov->key, s); 

  BUG; // to do

}



#define SCALE_COV 0
#define SCALE_SCALE 1
#define SCALE_PENALTY 2
void nonstatscale(double *x, double *y, model *cov, double *v){  
  model
    *next = cov->sub[SCALE_COV],
    *scale = cov->sub[SCALE_SCALE],
    *penalty = cov->sub[SCALE_PENALTY];
  double dummy, sx, sy,
    normSq = 0.0;
  int dim = OWNXDIM(0);
  assert(scale != NULL);
  FCTN(x, scale, &sx);
  FCTN(y, scale, &sy);
  if (sx <= 0.0 || sy  <= 0.0)
     ERR1("'%.50s' must be a positive function", KNAME(DSCALE));
  for (int i=0; i<dim; i++) {
    dummy = x[i] / sx - y[i] / sy;
    normSq += dummy * dummy;
  }
  if (penalty == NULL) {
     dummy = sx - sy;   
  } else {
    double fx, fy;
    FCTN(&sx, penalty, &fx);
    FCTN(&sy, penalty, &fy);
    dummy = fx - fy;   
  }
  double z = SQRT(normSq + dummy * dummy);  
  COV(&z, next, v);       
}



int checkscale(model *cov){
  ASSERT_ONESYSTEM;
  model   *next = cov->sub[SCALE_COV],
    *scale = cov->sub[SCALE_SCALE],
    *penalty = cov->sub[SCALE_PENALTY];
  int err,
    dim = OWNXDIM(0);

  if (next==NULL || scale == NULL) SERR("submodel(s) missing");

  ASSERT_QUASIONESYSTEM;
  if ((err = CHECK(next, OWNLOGDIM(0) + 1, 1, PosDefType,
		   XONLY, ISOTROPIC, SCALAR, cov->frame)) != NOERROR) {
    RETURN_ERR(err);
  }	
  if ((err = CHECK(scale, dim, dim, TrendType, XONLY, SYMMETRIC, SCALAR,
		   TrendType)) != NOERROR) RETURN_ERR(err);
  if (penalty!= NULL &&
      (err = CHECK(penalty, 1, 1, TrendType, XONLY, SYMMETRIC, SCALAR,
		   TrendType)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}




void kappablend(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc) {
   *nc = i == BLEND_THRES ? 1 : -1;
   *nr = SIZE_NOT_DETERMINED;
}

void nonstatblend(double *x, double *y, model *cov, double *v){
  model *blend = cov->sub[BLEND_BLEND],
    *multi = cov->sub[BLEND_MULTI];
  double bx, by,
    *thres =  P(BLEND_THRES);
  int i, j,
    n_thresholds = cov->nrow[BLEND_THRES],
    vdim = multi->vdim[0];
  
  FCTN(x, blend, &bx);
  for (i=0; i<n_thresholds; i++) if (bx <= thres[i]) break;
  i = i % vdim;
  
  FCTN(y, blend, &by);
  for (j=0; j<n_thresholds; j++) if (by <= thres[j]) break;
  j = j % vdim;

  // printf("vdim=%d %d\n", vdim, VDIM0);
  TALLOC_X2(w, vdim * vdim);
  NONSTATCOV(x, y, multi, w);
  *v = w[i + vdim * j];
  END_TALLOC_X2;
}

int checkblend(model *cov){
  ASSERT_ONESYSTEM;
  model *blend = cov->sub[BLEND_BLEND],
    *multi = cov->sub[BLEND_MULTI];
  int err,
    dim = OWNXDIM(0);

  kdefault(cov, BLEND_THRES, 0.5);
  int n=cov->vdim[BLEND_THRES];
  double *thres = P(BLEND_THRES);
  for (int i=1; i<n; i++)
    if (thres[i] <= thres[0])
      RFERROR("Threshold numbers must be given in strictly isotone ordering.");
  if ((err = CHECK(blend, dim, dim, TrendType, XONLY, CARTESIAN_COORD, SCALAR,
		   TrendType)) != NOERROR) RETURN_ERR(err);
  if ((err = CHECK(multi, dim, dim, PosDefType, KERNEL, SYMMETRIC,
		   SUBMODEL_DEP, GaussMethodType)) != NOERROR) RETURN_ERR(err);
  EXTRA_STORAGE;
  RETURN_NOERROR;
}

void rangeblend(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[BLEND_THRES] = RF_NEGINF;
  range->max[BLEND_THRES] = RF_INF;
  range->pmin[BLEND_THRES] = -10;
  range->pmax[BLEND_THRES] = 10;
  range->openmin[BLEND_THRES] = false;
  range->openmax[BLEND_THRES] = false;
}



void kappabubble(int i, model *cov, int *nr, int *nc) {
  if (i==BUBBLE_Z) {
    *nc = SIZE_NOT_DETERMINED;
    *nr = OWNXDIM(0);
    //    printf("nr=%d %d\n", *nr, *nc);
  } else  if (i==BUBBLE_BARYCENTRE) {
    *nc = *nr = 1;
  } else if (i == BUBBLE_MINSCALE || i == BUBBLE_WEIGHT) {
    if (PisNULL(BUBBLE_Z)) *nc = *nr = 0;
    else {
      *nc = 1;
      *nr = SIZE_NOT_DETERMINED;
    }
  } else *nc = *nr = -1;
}

#define tol_bubble 0.0
void nonstatbubble(double *X, double *Y, model *cov, double *v){
  model
    *next = cov->sub[BUBBLE_COV];
  bubble_storage *s = cov->Sbubble;
  assert(s != NULL);
  int
    dim = OWNXDIM(0),
    ix = (int) X[dim],
    iy = (int) Y[dim];
  
  double *x, *y, 
    *w = P(BUBBLE_WEIGHT),
    *z = P(BUBBLE_Z);
  int rankx = s->rank[ix],
    ranky = s->rank[iy];

  // printf("ix=%d %d %d %d %10g %10g \n", ix, iy, rankx, ranky, s->tau[0], s->tau[1]);
  
  if (rankx < ranky) {
    int k = rankx; rankx = ranky, ranky = k;
    k = ix; ix = iy; iy = k;
    x = Y;
    y = X;
  } else {
    x = X;
    y = Y;
  }

  //  printf("x=%3.1f y=%3.1f ", *x, *y);

  double
    // tau_x = s->tau[rankx], 
    tau_y = s->tau[ranky];
  // printf("tau_x %3.1f %3.1f ", tau_x, tau_y);
  assert(s->tau[rankx] >= tau_y);
  assert(rankx >= ranky);
  //#define max_tau tau_x
#define min_tau tau_y
  double dist = 0.0;
  for (int i=0; i<dim; i++) {
    double dummy = x[i] - y[i];
    dist += dummy * dummy;
  }

  //  printf("dist=%3.1f mintau=%3.1f ", dist, min_tau);
  dist *= min_tau; // first part of the distance
  
  
  for (int i=ranky + 1; i <= rankx; i++) {
    int endfor = s->end[i];
    double D = 0.0;
    for (int within=s->start[i]; within < endfor; within++) {
      double 
	*xx = x,
	*u = z + dim * within;
      for (int j=0; j<dim; j++) {
	double dummy = xx[j] - u[j];
	D += dummy * dummy * w[within];
      }      
    }
    double d_xi = s->tau[i] - s->tau[i-1];
    dist += D * d_xi; 
  }

  dist = SQRT(dist);
  // printf("dist=%3.1f\n",  dist);  
  COV(&dist, next, v);
}



int checkbubble(model *cov){
  //  PMI0(cov);
  
  ASSERT_ONESYSTEM;
  model   *next = cov->sub[BUBBLE_COV],
    *scale = cov->sub[BUBBLE_SCALE];
  int err,
    total = Gettotalpoints(cov),
    dim = OWNXDIM(0);
  kdefault(cov, BUBBLE_BARYCENTRE, false);
  
  if (next==NULL || scale == NULL) SERR("submodel(s) missing");

  ASSERT_QUASIONESYSTEM;
  if ((err = CHECK(next, INFDIM, 1, PosDefType,
		   XONLY, ISOTROPIC, SCALAR, cov->frame)) != NOERROR) {
    RETURN_ERR(err);
  }	
  if ((err = CHECK(scale, dim, dim, ShapeType,
		   XONLY, CARTESIAN_COORD, SCALAR, TrendType)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  if (cov->Sbubble == NULL) {
    NEW_STORAGE(bubble);
    bubble_storage *s = cov->Sbubble;
   
    double *tau = (double*) MALLOC(sizeof(double) * total);
    s->tau = (double*) MALLOC(sizeof(double) * total); // kann viel zu viel sein
    //  int *rk = (int*) MALLOC(sizeof(int) * total);
    s->rank = (int*) MALLOC(sizeof(int) * total);
    int *ordr = (int*) MALLOC(sizeof(int) * (total + 1));
 
    FctnExtern(cov, cov, scale, tau, true);
    
    for (int i=0; i<total; i++) tau[i] = 1 / (tau[i] * tau[i]);
    Ext_ordering(tau, total, 1, ordr);
    //   Ext_orderingInt(ordr, total, 1, rk); //rank

    //    printf("tau:");for (int i=0;i<total;i++) printf("%10g ",tau[i]);printf("\n");
    //   printf("ord:");for (int i=0;i<total;i++) printf("%d ",ordr[i]);printf("\n");
   //printf("rk:");for (int i=0;i<total;i++) printf("%d ",rk[i]);printf("\n");

   int cur_rk = 0,
     idx = ordr[0];
    s->tau[0] = tau[idx];
    s->rank[idx] = cur_rk;
    for (int i=1; i<total; i++) {
      idx = ordr[i];
      double value = tau[idx];
      if (FABS(value - s->tau[cur_rk]) > tol_bubble) s->tau[++cur_rk] = value;
      s->rank[idx] = cur_rk;
      //   printf("i=%d ord=%d val=%10g; %d %10g\n", i, idx, value, cur_rk, s->tau[cur_rk]);
    }
    //  printf("rank:");for (int i=0;i<total;i++) printf("%d ",s->rank[i]);printf("\n");


    s->start = (int*) CALLOC(cur_rk + 1, sizeof(int));
    s->end = (int*) CALLOC(cur_rk + 1, sizeof(int));

    if (PisNULL(BUBBLE_Z)) {
      bool barycentre = P0INT(BUBBLE_BARYCENTRE);
      if (!PisNULL(BUBBLE_WEIGHT) || !PisNULL(BUBBLE_MINSCALE)) {
	ERR3("'%.50s' and '%.50s' may only be given if '%.50s' is given",
	     KNAME(BUBBLE_WEIGHT), KNAME(BUBBLE_MINSCALE), KNAME(BUBBLE_Z) );
      }
      if (!Getgrid(cov))
	ERR1("'%.50s' must be always given for arbitrary locations.",
	    KNAME(BUBBLE_Z))
      PALLOC(BUBBLE_Z, dim, total);
      PALLOC(BUBBLE_WEIGHT, total, 1);
      PALLOC(BUBBLE_MINSCALE, total, 1);
      int 
	twodim = 2 * dim,
	i = 0,	// (current) number of groups of same scale
	end = 0, // eins mehr als unten
	points = 0, // number of z's and weights
	*len = (int*) MALLOC(sizeof(int) * dim), // grid len
	*ix = (int*) MALLOC(sizeof(int) * dim), // grid coordinates
	*cumlen = (int*) MALLOC(sizeof(int) * (dim + 1));
      double
	*z = P(BUBBLE_Z),
	*w = P(BUBBLE_WEIGHT),
	*ms = P(BUBBLE_MINSCALE),
	*left = (double*) MALLOC(sizeof(double) * dim), // grid start
	*step = (double*) MALLOC(sizeof(double) * dim); // grid step

      location_type *loc = Loc(cov);
      cumlen[0] = 1;
      for (int d=0; d<dim; d++) {
	len[d] = (int) loc->xgr[d][XLENGTH];
	step[d] = loc->xgr[d][XSTEP];
	left[d] = loc->xgr[d][XSTART];
	cumlen[d+1] = cumlen[d] * len[d];
      }
      
      while (end < total) {
	int
	  o = ordr[end],
	  start = end,
	  cmp_rk = s->rank[o];
	double sum = 0.0,
	  *z_start = z,
	  minscale = 1.0 / SQRT(tau[o] + tol_bubble);
	s->start[i] = points;
	// printf("%d s->rank[start]=%d\n", start, s->rank[o]);
	while (end < total && s->rank[o] == cmp_rk) {
	  // printf("\n end=%d %d %3.1f ", end, o, tau[o]);
	  int oo = o;
	  for (int d=0; d<dim; d++) {
	    ix[d] = oo % len[d];
	    oo /= len[d];
	    // printf("%d ", ix[d]);
	  }
	  double tau0 = tau[o],
	    gradient = 0.0;
	  for (int k=0; k<twodim; k++) {
	    idx = k / 2;
	    int
	      sign = 2 * (k % 2) - 1,
	      i_neigh = ix[idx] + sign;
	    if (i_neigh < 0 || i_neigh >= len[idx]) continue;
	    double
	      tau_n = tau[o + sign * cumlen[idx]],
	      delta_tau = FABS((tau_n - tau0) / step[idx]);	    
	    gradient = MAX(gradient, delta_tau);
	  }
	  //	  printf("grd=%3.1f\n", gradient);
	  if (gradient != 0.0) {
	    for (int d=0; d<dim; d++) {
	      z[d] = left[d] + (double) ix[d] * step[d];
	      // printf("%d %10g %10g %d %10g\n", d, z[d], left[d], ix[d], step[d]);
	    }
	    z += dim;
	    w[points] = gradient;
	    sum += gradient;
	    points++;
	  }
	  end++;
	  o = ordr[end]; // not ok
	}

	// printf("\n sum=%10g\n", sum);
	
	if (sum == 0.0 && (start > 0 || end < total))
	  ERR1("strange field -- contact author and give '%.50s' explicitely",
	       KNAME(BUBBLE_Z));

	for (int k=s->start[i]; k<points; k++) {
	  w[k] /= sum;
	  ms[k] = minscale; // if (barycentre) only ms[s->start[i]] is used!
	}

	if (barycentre) {
	  int s_start = s->start[i];
	  z = z_start;
	  double w_cur = w[s_start];
	  for (int d=0; d<dim; d++) z_start[d] *= w_cur;
	  for (int k=s_start+1; k<points; k++, z+=dim) {
	    w_cur = w[k];
	    for (int d=0; d<dim; d++) z_start[d] += w_cur * z[d];
	  }
	  z = z_start + dim;
	  w[s_start] = 1.0;
	  points = s_start + 1;

	  //	  for(int d=0; d<dim; d++) printf("d=%d %10g\n", d, z_start[d]);
	} 

	s->end[i++] = points; // ganz zum schluss;
	// nur bei z=NULL gilt s->end[i] = s->start[i+1] !!
	

	// printf("i=%d start/end = %d %d\n", i-1, s->start[i-1], s->end[i-1]);
      } // while end < total

      
      //      z = P(BUBBLE_Z);  for(i=0; i<points; i++) { printf("%d (%3.1f %3.1f %5.5f %5.5f)\n", points, z[(i)*2], z[(i)*2+1], w[i], ms[i]); }
      
      cov->nrow[BUBBLE_WEIGHT] = cov->ncol[BUBBLE_Z] =
	cov->nrow[BUBBLE_MINSCALE] = points;
      FREE(len);
      FREE(ix);
      FREE(cumlen);
      FREE(step);
      FREE(left);

      //APMI(cov);      
    } else {
      if (PisNULL(BUBBLE_MINSCALE)) {
	for (int j=0; j<=cur_rk; j++) s->end[j] = cov->ncol[BUBBLE_Z];
      } else {
	double *minscale = P(BUBBLE_MINSCALE);
	if (cov->nrow[BUBBLE_MINSCALE] != cov->ncol[BUBBLE_Z])
	  ERR2("lengths of '%.50s' and '%.50s' do not match.", KNAME(BUBBLE_MINSCALE),
	       KNAME(BUBBLE_Z));
	int start,
	  points = 0, // minscales for z
	  i = 0, // classes of scales of the field
	  len = cov->nrow[BUBBLE_MINSCALE];
#define bubble_eps 1e-15	

	while (i <= cur_rk) {
	  //	  printf("pts=%d, minscale=%10g %d st=%d len=%d\n", points, minscale[points], cur_rk, start, len);
	  if (points >= len) {
	    ERR1("lowest actual scale value passes lowest value of '%.50s'",
		 KNAME(BUBBLE_MINSCALE));
	  }
	  double
	    ms = minscale[points],
	    maxtau = 1.0 / (ms * ms) + bubble_eps;
	  start = points++;
	  for ( ; points < len; points++) {
	    if (minscale[points] > ms)
	      ERR1("'%.50s' not in descending order.", KNAME(BUBBLE_MINSCALE));
	    if (minscale[points] < ms) break;
	    // "==" : same class of scales
	  }
	  while (i <= cur_rk && s->tau[i] <= maxtau) {
	    s->start[i] = start;
	    s->end[i] = points;
	    i++;
	  }
	  //	  printf("i=%d s->tau=%10g %10g %d pts=%d\n", i, s->tau[i],  maxtau, s->tau[i] <= maxtau, points);
	}
      }

      
      int i = 0;// classes of scales of the field
      //      PMI0(cov);
      if (PisNULL(BUBBLE_WEIGHT)) {
	PALLOC(BUBBLE_WEIGHT, cov->ncol[BUBBLE_Z], 1);
	double *w = P(BUBBLE_WEIGHT);
	while(i <= cur_rk) {
	  int s_end = s->end[i];

	  //  printf("start=%d %d %d\n", s->start[i], s->end[i], s_end);
	  
	  double weight = 1.0 / (s_end - s->start[i]);
	  for (int points=s->start[i]; points < s_end; w[points++] = weight);
	  i++;
	  while(i <= cur_rk && s->start[i] == s->start[i-1]) {
	    if (s->end[i] != s->end[i-1])
	      ERR0("Identical starting position must have identical end position");
	    i++;
	  } // while i <= cur_rk
	} // weight not NULL
      } else { // bubble not NULL
	if (cov->nrow[BUBBLE_WEIGHT] != cov->ncol[BUBBLE_Z])
	  ERR2("lengths of '%.50s' and '%.50s' do not match.", KNAME(BUBBLE_WEIGHT),
	       KNAME(BUBBLE_Z));
	double *w = P(BUBBLE_WEIGHT);
	while(i <= cur_rk) {
	  int s_end = s->end[i];
	  double sum=0.0;

	  //printf("i=%d start/end = %d %d\n", i, s->start[i], s_end);
	  
	  for (int points=s->start[i]; points<s_end; points++) sum+=w[points];
	  
	  // printf("delta %10e\n", sum - 1.0	  );
	  if (FABS(sum - 1.0) > 1e-14)
	    ERR0("weights do not partially sum up to 1.");
	  i++;
	  while(i <= cur_rk && s->start[i] == s->start[i-1]) {
	    if (s->end[i] != s->end[i-1])
	      ERR0("identical starting position must have identical end position");
	    i++;
	  }
	} // while i <= cur_rk
      } // weight not NULL
    } // z not NULL    

    FREE(tau);
    FREE(ordr);
  } // end create Sbubble
  //  PMI0(cov);
  

  RETURN_NOERROR;
}



void rangebubble(model VARIABLE_IS_NOT_USED *cov, range_type *range){

  range->min[BUBBLE_Z] = RF_NEGINF;
  range->max[BUBBLE_Z] = RF_INF;
  range->pmin[BUBBLE_Z] = -1000;
  range->pmax[BUBBLE_Z] = 1000;
  range->openmin[BUBBLE_Z] = true;
  range->openmax[BUBBLE_Z] = true;

  range->min[BUBBLE_MINSCALE] = 0;
  range->max[BUBBLE_MINSCALE] = RF_INF;
  range->pmin[BUBBLE_MINSCALE] = 1e-5;
  range->pmax[BUBBLE_MINSCALE] = 10;
  range->openmin[BUBBLE_MINSCALE] = true;
  range->openmax[BUBBLE_MINSCALE] = true;
  
  range->min[BUBBLE_WEIGHT] = 0;
  range->max[BUBBLE_WEIGHT] = RF_INF;
  range->pmin[BUBBLE_WEIGHT] = 1e-5;
  range->pmax[BUBBLE_WEIGHT] = 10;
  range->openmin[BUBBLE_WEIGHT] = true;
  range->openmax[BUBBLE_WEIGHT] = true;

  booleanRange(BUBBLE_BARYCENTRE); 

}




#define DERIV_WHICH 0
void kappaderivative(int i, model VARIABLE_IS_NOT_USED *cov,
		  int *nr, int *nc){
  *nc = 1;
  *nr = i == DERIV_WHICH ? SIZE_NOT_DETERMINED : -1;
} 

/* Div - Delta Delta^T */
void derivative(double *x, model *cov, double *v) { // div -free !!

  
  /*
    r = \| h\| 
    mit h = x-y 
    d/d h_i f(r) = f'(r) h_i / r
    d^2/d h_i^2 f(r) =   f''(r) h_i^2/r^2 + f'(r) (r^2 - h_i^2) / r^3
    d^2/d h_i h_j f(r) = f''(r) h_i h_j/r^2 + f'(r) (- h_i h_j) / r^3

    r = 0:
    d/d h_i f(r) |_0 = 0
    d^2/d h_i^2 f(r) |_0 = f''(0)
    d^2/d h_i h_j f(r) = 0
  */
  model *next = cov->sub[0];  
  defn *N = DefList + NEXTNR; // gatter oder keine gatter??
  double
    r, rSq, D, D2;
  int
    dim = OWNLOGDIM(0), 
      // currently no space-time allowed -> BA thesis
    dimP1 = dim + 1,
    dimP2 = dim + 2,
    dimPsq = dimP1 * dimP1; 
  double extra_a[MAXVDIM * MAXVDIM];
  double *zz = PisNULL(DERIV_WHICH) ? v : extra_a;

  
  rSq = 0.0;
  for (int i=0; i<dim; i++) rSq += x[i] * x[i];
  r = SQRT(rSq);
  N->D(&r, next, &D);
  N->D2(&r, next, &D2);
  N->cov(&r, next, zz); // v[0]
  
   if (rSq == 0.0) {
    for (int i=1; i< dimPsq; i++) zz[i] = 0.0;
    double diag = -D2;
    for (int i=dimP2; i<dimPsq; i+=dimP2) zz[i] = diag;  // diag 2,3
  } else {
    double 
      D1n =  D / r,
      D1n3 = D / (r * rSq),      
      D2nsq = - D2 / rSq,
      DX = D1n3 + D2nsq;
    
    for (int i=1; i<=dim; i++)  // kante links und oben
      zz[i * dimP1] = -(zz[i] = x[i-1] * D1n);
    
    for (int i=dimP2, k=0; k<dim; k++, i++) {
      for (int l=0; l<dim; l++) {
	zz[i++] = x[k] * x[l] * DX - (k == l ? D1n : 0.0);
      }
    }
  }

  if (!PisNULL(DERIV_WHICH)) {
    int len = NROW(DERIV_WHICH);
    for (int i=0; i<len; i++) {
      int ii = PINT(DERIV_WHICH)[i] - dimP2; // immer PINT() - 1 da R -> C
      for (int j=0; j<len; j++) {
	v[i + j * len] = zz[ii + PINT(DERIV_WHICH)[j] * dimP1]; // -1 s.o.
      }
    }
  }

}


int checkderivative(model *cov) {
  model  *next = cov->sub[0];
  // been programmed yet 
  int err, i,
    dim = OWNLOGDIM(0);

  // statt STATIONARY : VARIOGRAM ?? s. Paper mit Scheuerer

  if ((err = CHECK(next, dim, 1, PosDefType, OWNDOM(0), ISOTROPIC,
		   SCALAR, EvaluationType)) != NOERROR) {
    //if ((err =CHECK(next, dim, 1, PosDefType, OWNDOM(0), DOUBLEISOTROPIC,
    // SCALAR, EvaluationType)) != NOERROR)
    RETURN_ERR(err);
  }

  if (next->full_derivs < 2) SERR("2nd derivative of submodel not defined");
  if (dim >= MAXVDIM) SERR("too high dimensions");
  
  setbackward(cov, next);
  int diffpref = MIN(2, PREF_BEST - cov->pref[CircEmbed]);
  if (diffpref > 0) cov->pref[CircEmbed] += diffpref;

  int nwhich = NROW(DERIV_WHICH),
    components = dim + 1;
  if (nwhich > 0) {
    for (i=0; i<nwhich; i++) {
      if (PINT(DERIV_WHICH)[i] <= 0 || PINT(DERIV_WHICH)[i] > components) {
	SERR4("value %.50s[%d]=%d outside range 1:%d.",
	      KNAME(i), i+1, PINT(DERIV_WHICH)[i], components);
      }
    }
  } else {
    nwhich = components;
  }
  for (i=0; i<dim; i++) cov->mpp.maxheights[i] = RF_NA;
  VDIM0 = VDIM1 = nwhich;

  RETURN_NOERROR;
}

void rangederivative(model *cov, range_type *range) {
  int dim = OWNLOGDIM(0),
    max = dim + 1;
  range->min[DERIV_WHICH] = 1;
  range->max[DERIV_WHICH] = max;
  range->pmin[DERIV_WHICH] = 1;
  range->pmax[DERIV_WHICH] = max;
  range->openmin[DERIV_WHICH] = false;
  range->openmax[DERIV_WHICH] = false;
}
