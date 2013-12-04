/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

 Copyright (C) 2005 -- 2013 Martin Schlather

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
#include "Covariance.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

#define BINARY_P 0


/*
  Felix' integral 
  + trafo  v = (1-r) / (1+r)
  + taylor 
  + GR 2.211/2

  V ~ 25 : Grenze zur Asymptoic
  GR 8.357


fctn = function(rho, t) exp(-t^2 / (1+rho)) / sqrt(1-rho^2)
C = function(r, t) (2*pi)^(-1) * integrate(fctn, 0, r, t=t)$value
fctn1 = function(rho, t) exp(-t^2 * rho / 2) / sqrt(rho) / (1 + rho)
C1 = function(r, t)  exp(-t^2 / 2) / (2*pi) * integrate(fctn1, (1-r)/(1+r), 1, t=t)$value
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
  (U = 2 * sqrt(v) * sum((-v)^d / (2 * d + 1) * s))
  v =  1
  (O = 2 * sqrt(v) * sum((-v)^d / (2 * d + 1) * s))
  (O - U) / (2 * pi)
  C(r, t)
  
  v = (1-r) / (1+r)
  (U = 2 * atan(sqrt(v)) + 2 * sqrt(v) * sum((-v)^d / (2 * d + 1) * (exp(-a) * s - 1)) )
  v =  1
  (O = 2 * atan(sqrt(v)) + 2 * sqrt(v) * sum((-v)^d / (2 * d + 1) * (exp(-a) * s - 1)) )
  (O - U) / (2 * pi) 
  C(r, t)
  print( (O - U) / (2 * pi)  - C1(r,t))#  //


model = RMexp()
x = seq(0.1, 1, len=10)
t = 1.67
for (r in x) print(C1(r,t=t))#  //
RFcov(model, log(x))
RFcov(RMbernoulli(model, t=t), log(x))
for (r in x) print(C1(r,t=t))  #  //

 */


#define Binary_Eps 1e-13
void binary(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double a, var, V, r, expMa, Vd,
    s, Oned, sum, sumOne, summand, summandOne, 
    d,
    factor,
    t = cov->p[BINARY_P][0];
  
  COV(ZERO, next, &var);
  a = 0.5 * t * t / var;
  expMa = exp(-a);
 
  COV(x, next, &r);
  r /= var;

  // als BA-Arbeit
  if (r < -0.9)
    ERR("correlation of submodel must be >= -0.9 for numerical reasons");

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
  while (fabs(summand) > Binary_Eps || fabs(summandOne) > Binary_Eps) {
    sum += summand;
    sumOne += summandOne;
    d += 1.0;
    ad *= a / d;
    s += ad;
    Vd *= -V;
    Oned *= -1;
    factor =  (s - 1.0) / (2.0 * d + 1.0);
    summand = Vd * factor;
    summandOne = Oned * factor;
  }
  sum += summand;
  sumOne += summandOne;

  // 0.25 = 2 atan(sqrt(1)) / 2pi
  *v = 0.25 + INVPI * (sumOne - (atan(sqrt(V)) + sqrt(V) * sum));
}


int checkbinary(cov_model *cov) {
  cov_model
    *next = cov->sub[0];
  double v;
  int err = NOERROR;
  kdefault(cov, BINARY_P, 0.0);
  if ((err = CHECK(next, cov->tsdim,  cov->xdimprev, PosDefType,
		     cov->domown, cov->isoown,
		     SUBMODEL_DEP, cov->role)) != NOERROR) return err;
  setbackward(cov, next);
  cov->mpp.maxheight = 1.0;
  COV(ZERO, next, &v);
  return NOERROR;
}


void rangebinary(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[BINARY_P] = RF_NEGINF; 
  range->max[BINARY_P] = RF_INF;
  range->pmin[BINARY_P] = -4.0;
  range->pmax[BINARY_P] = 4.0;
  range->openmin[BINARY_P] = false;
  range->openmax[BINARY_P] = false;
}



//////////////////////////////////////////////////////////////////////
// extremal gaussian 

void extrgauss(double *x, cov_model *cov, double *v) {
  // BrownResnick to binary Gaussian
  cov_model *next = cov->sub[0];
  double var, z;
  
  COV(ZERO, next, &var);
  COV(x, next, &z);
  *v = 1.0 - sqrt(0.5 * (1.0 - z / var));
}


int check_extrgauss(cov_model *cov) {
  cov_model
    *next = cov->sub[0];
  double v;
  int err = NOERROR;
  if ((err = CHECK(next, cov->tsdim,  cov->xdimprev, PosDefType,
		     cov->domown, cov->isoown,
		     SUBMODEL_DEP, cov->role)) != NOERROR) return err;
  setbackward(cov, next);
  cov->mpp.maxheight = 1.0;
  COV(ZERO, next, &v);
  if (v != 1.0) SERR("only correlation functions allowed");
  return NOERROR;
}


//////////////////////////////////////////////////////////////////////
// Brown Resnick

#define BR_FACTOR 0.5
void brownresnick(double *x, cov_model *cov, double *v) {
  // BrownResnick to binary Gaussian
  cov_model *next = cov->sub[0];
  double z;  

  assert(cov->role == ROLE_COV);

  COV(ZERO, next, &z);
  COV(x, next, v);
  *v = 2.0 * pnorm(sqrt((z - *v) * BR_FACTOR), 0.0, 1.0, 0, 0);
}

void Dbrownresnick(double *x, cov_model *cov, double *v) {
  // b = BR_FACTOR
  // gamma(h) = C(0) - C(h) 
  // varphi(sqrt(gamma * b)) * (b C') / sqrt(gamma * b)
  // varphi density standard normal

  // BrownResnick to binary Gaussian
  cov_model *next = cov->sub[0];
  double z, abl, s;
  
  //  

  if (cov->role != ROLE_COV && cov->role != ROLE_MAXSTABLE) BUG;
  if (cov->taylorN < 1) BUG;

  if (cov->taylor[1][TaylorPow] == 0.0) {
    *v = 0.0;
    return;
  }
  
  if (*x != 0.0) {
    COV(ZERO, next, &z);
    COV(x, next, v);
    
    //     A  PMI(cov);
    //  printf("%ld %ld %s\n", CovList[next->gatternr].D, D_2, 
    //	   CovList[next->gatternr].name);
    assert(CovList[next->gatternr].D != NULL);
    assert(CovList[next->gatternr].D != ErrCov);

    Abl1(x, next, &abl);   // Vorsicht: abl = -\gamma' 


    abl *= BR_FACTOR;
    s = sqrt((z - *v) * BR_FACTOR); // sqrt(c * gamma)
    *v = dnorm(s, 0.0, 1.0, 0) * abl / s; // -\varphi * \gamma' / sqrt\gamma
    //                                      =-2 \varphi * (-C') / (2 sqrt\gamma)
    //                              mit C= c_0 -\gamma
 
    //     printf("dbrown %f %f %f %e\n", *x, abl, s, v); //assert(false);
    // 
    assert(*v <= 0);

  } else {
    if (cov->taylor[1][TaylorPow] < 1.0) {     
      *v = -RF_INF;
    } else if (cov->taylor[1][TaylorPow] == 1.0) {
      *v = fabs(cov->taylor[1][TaylorConst]);
      assert(*v > 0.0);
    } else BUG;
  }
  // 2 * c * varphi(s) * gamma' / s
}


void DDbrownresnick(double *x, cov_model *cov, double *v) {
  // D = \varphi * b C' / sqrt(\gamma b)

  // b = BR_FACTOR
  // gamma(h) = C(0) - C(h) 
  // varphi(sqrt(gamma * b)) [(b C')^2 / [2 sqrt (gamma * b)]
  //                          +(b C'') / sqrt (gamma * b)
  //                          +1/2 * (b C')^2 / sqrt (gamma * b)^3 ]
 // varphi density standard normal

  cov_model *next = cov->sub[0];
  double z, abl, abl2, s0, s;
  
  if (cov->role != ROLE_COV && cov->role != ROLE_MAXSTABLE) BUG;
  if (cov->taylor[1][TaylorPow] == 0.0) {
    *v = 0.0;
    return;
  }
   
  if (*x != 0.0) {
    COV(ZERO, next, &z);
    COV(x, next, v);
    Abl1(x, next, &abl);
    Abl2(x, next, &abl2);
    s0 = (z - *v) * BR_FACTOR;
    s = sqrt(s0); // sqrt(c * gamma)
    abl  *= BR_FACTOR;
    abl2 *= BR_FACTOR;
    *v = dnorm(s, 0.0, 1.0, 0) / s * (abl2 + abl * abl * 0.5 * (1/s0 + 1));
    
    //printf("br x=%f v=%e s=%f abl=%f s0=%f abl2=%f\n", *x, *v, s, abl, s0, abl2);
    assert(*v >= 0);

  } else {
    *v = cov->taylor[1][TaylorPow]==1 ? 0.0 : RF_INF;
  }
}


void D3brownresnick(double *x, cov_model *cov, double *v) {
  // D = \varphi * b C' / sqrt(\gamma b)

  // b = BR_FACTOR
  // gamma(h) = C(0) - C(h) 
  // varphi(sqrt(gamma * b)) [(b C')^3 / [4 sqrt (gamma * b)]
  //                         +3 (b C') (b C'')/ [2 sqrt (gamma * b)]
  //                         + (b C''') / [sqrt (gamma * b)]
  //                         + (b C')^2 / [2 sqrt (gamma * b)]
  //                         +3(b C') (b C'') / [2 sqrt (gamma * b)^3]
  //                         +3(b C')^3 / [4 sqrt (gamma * b)^5]

  cov_model *next = cov->sub[0];
  double z, abl, abl2, abl3, s0, s;
  
  if (cov->role != ROLE_COV && cov->role != ROLE_MAXSTABLE) BUG;
  if (cov->taylor[1][TaylorPow] == 0.0) {
    *v = 0.0;
    return;
  }
   
  if (*x != 0.0) {
    COV(ZERO, next, &z);
    COV(x, next, v);
    Abl1(x, next, &abl);
    Abl2(x, next, &abl2);
    Abl3(x, next, &abl3);
    s0 = (z - *v) * BR_FACTOR;
    s = sqrt(s0); // sqrt(c * gamma)
    abl  *= BR_FACTOR;
    abl2 *= BR_FACTOR;
    abl3 *= BR_FACTOR;
    *v = dnorm(s, 0.0, 1.0, 0) / s * 
      (abl3 + 
       1.5 * abl2 * abl * (1/s0 + 1) +
       abl * abl * abl * (0.25 + 0.5 / s0 + 0.75 / (s0 * s0)));

    //printf("br x=%f v=%e s=%f abl=%f s0=%f abl2=%f\n", *x, *v, s, abl, s0, abl2);
    assert(*v >= 0);

  } else {
    *v = cov->taylor[1][TaylorPow]==1 ? 0.0 : RF_NEGINF; // da zweite Ableitung im Ursprung Pol +\infty hat.
  }
}


int checkbrownresnick(cov_model *cov) {
  cov_model  
   *next = cov->sub[0];
  int err, idx, dim = cov->tsdim;

  if ((err = CHECK(next, dim,  dim, NegDefType, cov->domown, 
		     cov->isoown, SUBMODEL_DEP, ROLE_COV)) != NOERROR) {
    return err;
  } 
  
  idx = isPosDef(next->typus); // Taylorentw. in 2 Glieder falls pos def.
  assert(idx == 0);

  //  assert(!idx);

  if (next->taylor[idx][TaylorPow] >= 2.0) {//for pos/neg def only == 2 possible
    cov->full_derivs = 1;
  } else cov->full_derivs = 0;
  cov->rese_derivs = next->rese_derivs;
  if (cov->rese_derivs > 3) cov->rese_derivs = 3;
  if (next->taylorN < 1 + idx) cov->rese_derivs = 0;
  // else if (next->taylor[idx][TaylorPow] == 2) {
  //  if (next->taylorN < 2 + idx) cov->rese_derivs = 0;
    // else if (cov->rese_derivs > 2) cov->rese_derivs = 2; 
  //  }
    

  setbackward(cov, next);
  cov->monotone = isBernstein(next) ? GNEITING_MON : 
    isMonotone(next) ? MONOTONE : NOT_MONOTONE;

  if (next->taylorN >= 1 + idx &&  next->taylor[idx][TaylorConst] < 0.0) {
    // 2 \Phi(sqrt(gamma * br_f)) =
    //1+ 2 * phi(0) * sqrt(gamma * br_f) - 2 phi(0) / 6 *sqrt(gamma*br_f)^3
    //sei gamma(x) = c x^alpha + d x^beta, beta > alpha. Dann gilt
    // also  2 \Phi(sqrt(gamma * br_f)) ~
    // 1 + 2 * phi(0) * sqrt((c x^alpha + d x^beta)* br_f)
    //   - 2 * phi(0) / 6 * sqrt((c x^alpha + d x^beta) *br_f)^3
    // = 1 + 2 * phi(0)[ (c f x^a)^{1/2} (1 + 0.5 d / c * x^{b-a})
    //                   (c f x^a)^{3/2} (1 + 1.5 d / c * x^{b-a}) ]
    // ~ 1 + 2 * phi(0) (c f x^a)^{1/2}
    //     + phi(0)(c f)^{1/2} d / c * x^{b - a/2})            (*)
    //     + 2 * phi(0) (c f x^a)^{3/2}                        (**)
    // 
    //
    // da a <= 2 und b > a gilt:
    // 1. Abltg ~ phi(0) (c f)^{1/2} x^{a/2 - 1}
    //           + phi(0) (c f)^{1/2} d / c * (b - a/2) x^{b - a/2 - 1} 
    //           + 3 phi(0 (c f)^{3/2} a x^{3/2 * a - 1} 
    // 2. Abltg ~ phi(0) (c f)^{1/2} (a/2 - 1) x^{a/2 - 2} ~ infty,   a != 2
    // und fuer a=2 gilt folgendes: da gamma neg def ist e^{-\gamma} pos def
    // und mit sasvari folgt aus der Taylor entwicklung dass (a=2) < b <= 4).
    // Andererseits dominiert x{3/2 a} den Term x^{b - a/2} falls 
    //  3/2 a < b - a/2  <=> b > 2a = 4. 
    // Also ist b=4 ein Grenzfall, wo beide (*) und (**) eine Rolle spielen.
    // Ansonsten nur (*).
    //
 
    assert(next->taylor[idx][TaylorConst] <= 0.0);
     
    cov->taylorN = 2;
    cov->taylor[0][TaylorConst] = 1.0;
    cov->taylor[0][TaylorPow] = 0.0;    
    double
      next_taylor_const =  - next->taylor[idx][TaylorConst],
      g = sqrt(next_taylor_const  * BR_FACTOR * 0.5 / M_PI);
    cov->taylor[1][TaylorConst] = - 2 * g;
    cov->taylor[1][TaylorPow] = 0.5 * next->taylor[idx][TaylorPow];
    if (next->taylor[idx][TaylorPow] == 2) {
      if (next->taylorN >= 2 + idx) {
	cov->taylorN = 3;
	if (next->taylor[idx + 1][TaylorConst] != 0) {
	  cov->taylor[2][TaylorConst] = 
	    g * next->taylor[idx + 1][TaylorConst] / next_taylor_const;
	  cov->taylor[2][TaylorPow] = 
	    next->taylor[idx + 1][TaylorPow] - 0.5*next->taylor[idx][TaylorPow];
	} else {	// Spezialfall:L fractionelle BR
	  cov->taylor[2][TaylorConst] = 0.0;
	  cov->taylor[2][TaylorPow] = 4.0;
	}
	if (next->taylor[idx + 1][TaylorPow] == 4) {
	  cov->taylor[1][TaylorConst] += 2 * g * next_taylor_const * BR_FACTOR;
	}
      } else cov->taylorN = 0;
    }
  } else cov->taylorN = 0;

  
  if (next->tailN >= 1 && next->taylorN >= 1) {
    cov->tailN = 1;    
    cov->tail[0][TaylorPow] = -0.5 * next->tail[0][TaylorPow];
    if (next->tail[0][TaylorPow] > 0) {  
      assert( next->tail[0][TaylorConst] < 0.0);
      double next_tail_const = - next->tail[0][TaylorConst];
      cov->tail[0][TaylorConst] = 
	2.0 / sqrt(2.0 * M_PI * BR_FACTOR * next_tail_const);
      cov->tail[0][TaylorExpConst] = 0.5 * BR_FACTOR * next_tail_const;
      cov->tail[0][TaylorExpPow] = next->tail[0][TaylorPow];
    } else {
      cov->tail[0][TaylorConst] = 
	2.0 / sqrt(2.0 * M_PI * BR_FACTOR * next->tail[0][TaylorConst])
	* exp(-0.5 * BR_FACTOR * next->tail[0][TaylorConst]);
      cov->tail[0][TaylorPow] = cov->tail[0][TaylorExpConst] = 
	cov->tail[0][TaylorExpPow] = 0.0;
    }     
  } else cov->tailN = 0;

  if (cov->taylorN < 1) cov->rese_derivs = 0;

  cov->mpp.maxheight = 1.0;
 
  return NOERROR;  
}



int struct_brownresnick(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {  
  cov_model *next = cov->sub[0];

  if (cov->role == ROLE_SMITH) {
    if (1 > next->taylorN  || 1 > next->tailN) 
      SERR2("role '%s' not possible for submodel '%s'", 
	    ROLENAMES[cov->role], NICK(next));

    BUG;

    assert(false);

    // shape ist somit die Ableitung, falls d=1 und i.w. die 
    // zweifache Ableitung, falls d=3

    // hier ist auch Taylor zu setztn fuer den neuen Shape

    // STRUCT(next);
  } else ILLEGAL_ROLE;


  return NOERROR;
}

int init_brownresnick(cov_model VARIABLE_IS_NOT_USED *cov, storage VARIABLE_IS_NOT_USED *s) {
  // cov_model *next = cov->sub[0];
  return NOERROR;
}

void do_brownresnick(cov_model *cov, storage *s) {
  cov_model *next = cov->sub[0];
  DO(next, s); // nicht gatternr
}


// #define UNIF_LOGREFAREA dim
#define SIGN_P 0
void randomsign(double *x, cov_model *cov, double *v) { 
  cov_model *next = cov->sub[0];
  assert(next != NULL);
  COV(x, next, v);
  (*v) *= cov->q[0];
  // printf("%f ", cov->q[0]);
}

void lograndomsign(double *x, cov_model *cov, double *v, double *sign) { 
  cov_model *next = cov->sub[0];
  assert(next != NULL);
  LOGCOV(x, next, v, sign);
  (*sign) *= cov->q[0];
}

void randomsignInverse(double *v, cov_model *cov, double *x){
  cov_model *next = cov->sub[0];
  INVERSE(v, next, x);
}
void randomsignNonstatInverse(double *v, cov_model *cov, double *x, double *y){
  cov_model *next = cov->sub[0];
  NONSTATINVERSE(v, next, x, y);
}

int check_randomsign(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err, 
    size = 1;
  if (cov->q == NULL) {
    if ((cov->q  = (double*) CALLOC(sizeof(double),  size)) == NULL)
      return ERRORMEMORYALLOCATION;
    cov->qlen = size;
  }
  
  kdefault(cov, 0, 0.5);
  if ((err = checkkappas(cov)) != NOERROR) return err;
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, ShapeType,
		     cov->domown, cov->isoown, 
		     SCALAR, cov->role)) != NOERROR)  return err;
  setbackward(cov, next);
  return NOERROR;
}

void range_randomsign(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SIGN_P] = 0; 
  range->max[SIGN_P] = 1;
  range->pmin[SIGN_P] = 0.0001;
  range->pmax[SIGN_P] = 0.9999;
  range->openmin[SIGN_P] = false;
  range->openmax[SIGN_P] = false;
}

int init_randomsign(cov_model *cov, storage *s) {
  cov_model *next = cov->sub[0];
  int err;
  double Eminus;
  assert(next != NULL);
  if ((err=INIT(next, cov->mpp.moments, s))!= NOERROR) return err;
  if (next->fieldreturn && next->loggiven) 
    SERR("log return is incompatible with random sign");
  
  if (cov->mpp.moments >= 1) {
    cov->mpp.M[0] = next->mpp.M[0];
    cov->mpp.Mplus[0] = next->mpp.Mplus[0];
    Eminus = cov->mpp.Mplus[1] - cov->mpp.M[1];
    cov->mpp.Mplus[1] = 
      cov->p[SIGN_P][0] * (cov->mpp.Mplus[1] - Eminus) + Eminus; 
    cov->mpp.M[1] = 0.0;
  }
  cov->mpp.maxheight = next->mpp.maxheight;
  cov->fieldreturn = next->fieldreturn;
  cov->origrf = false;
  cov->rf = next->rf;  

  // assert(cov->mpp.maxheight == 1.00);

  return err;
}

void do_randomsign(cov_model *cov, storage *s) {
  cov_model *next = cov->sub[0];
  assert(next != NULL);
  DO(next, s); // nicht gatternr
  cov->q[0] = 2.0 * (UNIFORM_RANDOM <= cov->p[SIGN_P][0]) - 1.0;
  if (cov->q[0] != 1.0 && next->fieldreturn) { 
    assert(cov->q[0] == -1.0);
    if (next->loggiven) ERR("log return is incompatible with random sign");
    int i,
      endfor = Loc(next)->totalpoints;
    double *rf = cov->rf;
    for (i=0; i<endfor; i++) rf[i] = -rf[i];
  }
}


int struct_randomsign(cov_model *cov, cov_model **newmodel) {  
  int err = STRUCT(cov->sub[0], newmodel);
  //  assert(cov->sub[0]->mpp.maxheight == 1.0);
  return err;
}



void BR2BG(double *x, cov_model *cov, double *v) {
  // BrownResnick to binary Gaussian
  cov_model *next = cov->sub[0];
  double z;
  COV(ZERO, next, &z);
  COV(x, next, v);
  *v = cos(M_PI * (1.0 - 2.0 * pnorm(sqrt( (z - *v) / 8.0), 0.0, 1.0, 0, 0 )));
}

int check_BR2BG(cov_model *cov) {
  cov_model *next = cov->sub[0];
  double v, t;
  int err;
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType,
		     cov->domown, cov->isoown, 
		     SCALAR, cov->role)) != NOERROR)  return err;
  setbackward(cov, next);
  cov->mpp.maxheight = 1.0;
  if (next->pref[Nothing] == PREF_NONE) return ERRORPREFNONE;

  COV(ZERO, next, &v);
  t = qnorm(0.5 * (1.0 - 0.5), 0.0, 1.0, 0, 0);

  if (v > 8.0 * t * t)
    SERR("variance must be less than 8(erf^{-1}(1/sqrt 2))^2 = 1.8197");
  return NOERROR;  
}


void BR2EG(double *x, cov_model *cov, double *v) {
  // BrownResnick to binary Gaussian
  cov_model *next = cov->sub[0];
  double z;
  COV(ZERO, next, &z);
  COV(x, next, v);
  z = 1.0 - 2.0 * pnorm(sqrt( (z - *v) / 8.0), 0.0, 1.0, 0, 0 );
  *v = 1.0 - 2.0 * z * z;
}

int check_BR2EG(cov_model *cov) {
  cov_model *next = cov->sub[0];
  double v, t;
  int err;
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, 
		     cov->domown, cov->isoown, 
		     SCALAR, cov->role)) != NOERROR)  return err;
  setbackward(cov, next);
  cov->mpp.maxheight = 1.0;
  if (next->pref[Nothing] == PREF_NONE) return ERRORPREFNONE;
  COV(ZERO, next, &v);
  t = qnorm(0.5 * (1.0 - INVSQRTTWO), 0.0, 1.0, 0, 0);
  if (v > 8.0 * t * t)
    SERR("variance must be less than 8(erf^{-1}(1/2))^2 = 4.425098");
  return NOERROR;  
}



#define MASTEIN_NU 0
#define MASTEIN_DELTA 1
void MaStein(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double  nuG, gammas, v1, v2,
    **p = cov->p,
    nu = p[MASTEIN_NU][0],
    delta = p[MASTEIN_DELTA][0];

  COV(ZERO, next, &v1);
  COV(x + 1, next, &v2);
  nuG = nu + v1 - v2;
  if (nuG >= 80.0) {
    ERR("Whittle Matern function cannot be evaluated with parameter value b+g(t) greater than 80.");
  }
  gammas = lgammafn(nu + delta) - lgammafn(nu) -lgammafn(nuG + delta);
  double s = *x;
  *v = (s==0.0) ? exp(lgammafn(nuG) + gammas)
    : 2.0 * exp(nuG * log(0.5 * s) + gammas + log(bessel_k(s, nuG, 2.0)) - s);
}

int check_MaStein(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;

  if (cov->xdimown != 2) SERR("reduced dimension must be 2");

  if ((err = checkkappas(cov)) != NOERROR) return err;
  if ((err = CHECK(next, 1, 1, NegDefType, XONLY,
		     SYMMETRIC, SCALAR, ROLE_COV)) != NOERROR) return err;
  if (cov->ncol[MASTEIN_NU] != 1 || cov->nrow[MASTEIN_NU] != 1) 
    SERR("nu not scalar");
  if (cov->ncol[MASTEIN_DELTA] != 1 || cov->nrow[MASTEIN_DELTA] != 1) 
    SERR("d not scalar");
  // no setbackward
  cov->maxdim = next->maxdim;
  return NOERROR;
  // no setbackward !
}

void range_MaStein(cov_model *cov, range_type *range){
  range->min[MASTEIN_NU] = 0.0; // 
  range->max[MASTEIN_NU] = RF_INF;
  range->pmin[MASTEIN_NU] = 1e-2;
  range->pmax[MASTEIN_NU] = 10.0;
  range->openmin[MASTEIN_NU] = true;
  range->openmax[MASTEIN_NU] = true;

  range->min[MASTEIN_DELTA] = 0.5 * (double) (cov->tsdim - 1); // d
  range->max[MASTEIN_DELTA] = RF_INF;
  range->pmin[MASTEIN_DELTA] = range->min[MASTEIN_DELTA];
  range->pmax[MASTEIN_DELTA] = 10;
  range->openmin[MASTEIN_DELTA] = false;
  range->openmax[MASTEIN_DELTA] = true;
}



/* bivariate shift model */
#define SHIFT_DELAY 0
void kappashift(int i, cov_model *cov, int *nr, int *nc){
  *nc = 0;
  *nr = i < CovList[cov->nr].kappas  ? cov->tsdim : -1;
}
void shift(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double y[ShiftMaxDim], z[ShiftMaxDim],
    *jh, *ih, 
    *pv = v,
    *h = cov->p[SHIFT_DELAY];
  int i, j, d,
    tsdim = cov->tsdim,
    vdim = cov->vdim,
    vdimM1 = vdim - 1,
    vdimSq = vdim * vdim,
    vdimP1 = vdim + 1;
  
  COV(x, next, v);
  for (i=vdimP1; i<vdimSq; i+=vdimP1) v[i] = v[0];

  //printf("%d %d\n", vdim, vdimM1);
  // APMI(cov->calling);
  
  for (jh=h-tsdim, j=-1; j<vdimM1; j++, jh+=tsdim) {
    if (j < 0) for (d=0; d<tsdim; d++) z[d] = x[d];
    else for (d=0; d<tsdim; d++) z[d] = x[d] + jh[d];
    for (ih=h-tsdim, i=-1; i<vdimM1; i++, ih+=tsdim, pv++) {
      if (i==j) continue;
      if (i < 0) for (d=0; d<tsdim; d++) y[d] = z[d];
      else for (d=0; d<tsdim; d++) y[d] = z[d] - ih[d];  
      COV(y, next, pv);
    }
  }
}

int checkshift(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err,
    dim = cov->tsdim;
      
  if (cov->xdimown > ShiftMaxDim)
    SERR2("For technical reasons max. dimension for ave is %d. Got %d.", 
	  StpMaxDim, cov->xdimown);

  if ((err = checkkappas(cov)) != NOERROR) return err;
  if ((err = CHECK(next, dim, dim, PosDefType, XONLY,
		     (dim > 1) ? SYMMETRIC : ISOTROPIC, 
		     SCALAR, ROLE_COV)) != NOERROR) 
    return err;
  setbackward(cov, next);
  cov->vdim = cov->ncol[SHIFT_DELAY] + 1;
  return NOERROR;
}

void rangeshift(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SHIFT_DELAY] = -RF_INF;
  range->max[SHIFT_DELAY] = RF_INF;
  range->pmin[SHIFT_DELAY] = -1000;
  range->pmax[SHIFT_DELAY] = 1000;
  range->openmin[SHIFT_DELAY] = true;
  range->openmax[SHIFT_DELAY] = true;
}



/* Vector - Delta Delta^T */
#define VECTOR_A 0
#define VECTOR_D 1
void vector(double *x, cov_model *cov, double *v) {
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
  cov_model *next = cov->sub[0];
  double norm[2], normSq0, normL2, normT2, D, D2, diag,
    a = cov->p[VECTOR_A][0],
    b = - 0.5 * (1.0 + a);
  int i,
    td = cov->tsdim,
    dim = ((int *)cov->p[VECTOR_D])[0],
    dimP1 = dim + 1,
    dimsq = dim * dim;
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<td; i++) normT2 += x[i] * x[i];
  if (next->isoown == ISOTROPIC) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = sqrt(normT2); 
  }
  norm[0] = sqrt(normSq0);
  Abl1(norm, next, &D);
  Abl2(norm, next, &D2);

  if (normSq0 == 0.0) {
    diag = (b * dim + a) * D2;
    for (i=0; i<dimsq; i++) {
      //    print("%f D2 %f\n", diag, D2);
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
void vectorAniso(double *x, cov_model *cov, double *v) {
  /*
    operator:  -0.5 * (a + 1) laplace + a * hessian,  a in [-1,1]; a=1
  */

  cov_model *next = cov->sub[0];
  double *D,  laplace,
    a = cov->p[VECTOR_A][0],
    b = - 0.5 * (1.0 + a);
  int i, j, k, endfor,
    dim = ((int *)cov->p[VECTOR_D])[0],    
    dimP1 = dim + 1,
    dimsq = dim * dim,
    tsxdim = cov->tsdim,
    xdimsq = tsxdim * tsxdim,
    xdimP1 = tsxdim + 1,
    dimxdim = dim * tsxdim;
  
  D = (double*) MALLOC(sizeof(double) * xdimsq); // alloc in check schieben, wie in $; auch bei allen anderen Funktionen

  // should pass back the hessian
  HESSE(x, next, D);
  laplace = 0.0;


  //  print("return matrix %d\n", xdimsq);
  //  for (i=0; i<xdimsq; i+=tsxdim) {
  //      for (j=i; j <i+tsxdim; j++) {
  //	  print("%f;%d  ", D[j], j);
  //      }
  //      print("\n");
  //  }


  for (i=0; i<dimxdim; i+=xdimP1) {
    laplace += D[i];
    //      print("lplace %d %f\n", i, D[i]);
  }
  laplace *= b;

  k = 0;
  for (i=0; i<dimxdim; i+=tsxdim) {  
    endfor = i + dim; 
    for (j=i; j<endfor; j++) {
      //	  print(">> %d %d %f %f %f\n", k, j, D[j], a, laplace );
      v[k++] = D[j] * a;
    }
  }
  
  for (i=0; i<dimsq; i+=dimP1) {
    v[i] += laplace;
  }

  free(D);
}

int checkvector(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err,
    dim = cov->tsdim;
  isotropy_type isotropy = cov->domown;

  kdefault(cov, 0, 0.5); // a
  kdefault(cov, 1, 
	   (isotropy==SPACEISOTROPIC  || isotropy==ZEROSPACEISO) 
	   ? dim - 1 : dim); // Dspace
  if ((err = checkkappas(cov)) != NOERROR) return err;

  if ( (isotropy==SPACEISOTROPIC || isotropy==ZEROSPACEISO) 
       && ((int*) (cov->p[VECTOR_D]))[0] != dim - 1) {
    //      print(".... %d %d\n", ((int*) (cov->p[VECTOR_D]))[0], dim);
    SERR("for spatiotemporal models 'vector' must be applied to spatial part");
  }
  
  if (cov->tsdim != cov->xdimown || cov->tsdim != cov->xdimprev)
    return ERRORDIM;

  cov->nr = VECTOR; // wegen nr++ unten;
  if ((err = CHECK(next, dim, 1, PosDefType, cov->domown, ISOTROPIC,
		     SCALAR, ROLE_COV)) != NOERROR) {
    //     print("err=%d\n", err);
    if ((err = CHECK(next, dim, dim, PosDefType, cov->domown, ZEROSPACEISO,
		       SCALAR, ROLE_COV)) != NOERROR) {
      //	print("err2=%d\n", err);
      if ((err = CHECK(next, dim, dim, PosDefType, cov->domown, SYMMETRIC,
			 SCALAR, ROLE_COV)) != NOERROR) {
	//	 print("err3=%d\n", err);
	return err;
      }
    }
  }

  //  APMI(cov);

  setbackward(cov, next);
  cov->mpp.maxheight = RF_NAN;

  if (next->full_derivs < 2 && !next->hess) {
    //    print("%s (%d): full_derivs %d  %d\n", NICK(to), 
    //	  next->nr, next->full_derivs, CovList[next->nr].F_derivs);
    //  PMI(next);
    
       SERR("2nd derivative of submodel not defined (for the given paramters)");
  }
  isotropy = next->isoown;

  if (isotropy != ISOTROPIC && isotropy != SPACEISOTROPIC) {
      if (!next->hess) SERR("hess matrix not defined");
      cov->nr++;
  }

  cov->vdim = ((int*) (cov->p[VECTOR_D]))[0];
  next->delflag = DEL_COV;
  return NOERROR;
}

void rangevector(cov_model *cov, range_type *range){
  range->min[VECTOR_A] = -1.0;
  range->max[VECTOR_A] = 1.0;
  range->pmin[VECTOR_A] = -1.0;
  range->pmax[VECTOR_A] = 1.0;
  range->openmin[VECTOR_A] = false;
  range->openmax[VECTOR_A] = false;

  range->min[VECTOR_D] = 1.0;
  range->max[VECTOR_D] = cov->tsdim;
  range->pmin[VECTOR_D] = 1.0;
  range->pmax[VECTOR_D] = cov->tsdim;
  range->openmin[VECTOR_D] = false;
  range->openmax[VECTOR_D] = false;
}



int zzz=0;

/* Rot - Delta Delta^T */
void curl(double *x, cov_model *cov, double *v) {
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
 
  cov_model *next = cov->sub[0];
  cov_fct *N = CovList + next->nr; // !!! nicht gatternr !!!!
  double norm[2], normSq0, normL2, normT2, D, D2, D3, diag,
      a = -1.0, // curl free
      b = - 0.5 * (1.0 + a);
  int i,
      td = cov->tsdim,
      dim = td,
      // currently no space-time allowed -> BA thesis
      dimP1 = dim + 1,
      dimP2 = dim + 2,
      dimP3 = dim + 3,
      dimP2sqM1 =  dimP2 * dimP2 -1; // for three dimensions much
  // more complicated
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<td; i++) normT2 += x[i] * x[i];
  if (next->isoown == ISOTROPIC) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = sqrt(normT2); 
  }
  norm[0] = sqrt(normSq0);
  N->D(norm, next, &D);
  N->D2(norm, next, &D2);

  //  printf("next %s\n", N->name);

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


/* Div - Delta Delta^T */
void div(double *x, cov_model *cov, double *v) { // div -free !!
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
  cov_model *next = cov->sub[0];
  cov_fct *N = CovList + next->nr; // gatter oder keine gatter??
  double norm[2], normSq0, normL2, normT2, D, D2, D3, diag,
      a = 1.0, // divergenz free
      b = - 0.5 * (1.0 + a);
  int i,
      td = cov->tsdim,
      dim = td, 
      // currently no space-time allowed -> BA thesis
      dimP1 = dim + 1,
      dimP2 = dim + 2,
      dimP3 = dim + 3,
      dimP2sqM1 =  dimP2 * dimP2 -1; // for three dimensions much
  // more complicated
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<td; i++) normT2 += x[i] * x[i];
  if (next->isoown == ISOTROPIC) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = sqrt(normT2); 
  }
  norm[0] = sqrt(normSq0);
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
	v[i] = (i % dimP3 == 0) ? diag : 0.0; 
	v[i] += delta * x[k] * x[l];
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

//  for (i=0; i<=15; i++) {
//      if (i!=15) v[i] = 0.0;
//  }

}

int checkdivcurl(cov_model *cov) {
  cov_model  *next = cov->sub[0];
  // been programmed yet 
  int err,
      dim = cov->tsdim,
      spacedim;
  isotropy_type isotropy;

  // statt STATIONARY : VARIOGRAM ?? s. Paper mit Scheuerer

  if ((err = CHECK(next, dim, 1, PosDefType, cov->domown, ISOTROPIC,
		     SCALAR, ROLE_COV)) != NOERROR) {
    if ((err = CHECK(next, dim, 1, PosDefType, cov->domown, SPACEISOTROPIC,
		       SCALAR, ROLE_COV)) != NOERROR)
    return err;
  }

  //  PMI(next);
  
  if (next->full_derivs < 4) SERR("4th derivative of submodel not defined");
  if (cov->tsdim != 2) SERR("currently coded only for dim=2");
  isotropy = next->isoown;
  spacedim = dim - (isotropy == SPACEISOTROPIC);

  if (isotropy != ISOTROPIC && isotropy != SPACEISOTROPIC)
      SERR("submodel must be spaceisotropic");
  if (spacedim != 2)
       SERR("model currently only defined for the plane");
  
  setbackward(cov, next);
  cov->mpp.maxheight = RF_NAN;
  cov->vdim = spacedim + 2; // only correct for spacedim == 2!!!
  assert(spacedim == 2);
  next->delflag = DEL_COV;
  return NOERROR;
}

// void rangedivcurl(cov_model *cov, range_type *range){ }

/*

#define LP_P 0
void lp(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];  
  double z,
    p = cov->p[LP_P][0];
  int d,
    dim = cov->tsdim;
  for (z=0.0, d=0; d<dim; d++) {
    z += pow(fabs(x[d]), p);
  }
  z = pow(z, 1 / p);
  
  COV(&z, next, v);
}
int checklp(cov_model *cov) {
  bool skipchecks = GLOBAL.general.skipchecks;
  cov_model *next = cov->sub[0];  
  int err;

  kdefault(cov, 0, 1.0);
  if ((err = checkkappas(cov)) != NOERROR) return err;
 
  if ((err = CHECK(next, cov->tsdim + 1, 1, PosDefType, XONLY,ISOTROPIC,
		     SCALAR, ROLE_COV)) 
	!= NOERROR) return err;
  next->delflag = DEL_COV;

  // no setbackward
  double p =cov->p[LP_P][0];
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

  return NOERROR;
}
void rangelp(cov_model *cov, range_type *range){
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
void ma1(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];  
  double z,
    alpha = cov->p[MA1_ALPHA][0],
    theta = cov->p[MA1_BETA][0];
  COV(x, next, &z);
  *v = pow(theta / (1 - (1-theta) * z), alpha);  
}
int checkma1(cov_model *cov) {
//  bool skipchecks = GLOBAL.general.skipchecks;
  cov_model *next = cov->sub[0];  
//  double p;
  int err;

  kdefault(cov, 0, 1.0);
  kdefault(cov, 1, 0.5);
  if ((err = checkkappas(cov)) != NOERROR) return err;
       
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, cov->domown,
		     cov->isoown, SCALAR, ROLE_COV)) 
	!= NOERROR) return err;

  cov->semiseparatelast = false;
  cov->logspeed = 0.0;
  //  updatepref(cov, next);
  setbackward(cov, next);
  cov->mpp.maxheight = 1.0;
  return NOERROR;
}
void rangema1(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
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


void ma2(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];  
  double z, z0;
  COV(ZERO, next, &z0);
  COV(x, next, &z);
  z = z0 - z;
  *v = (z == 0) ? 1.0 : (1.0 - exp(-z)) / z;  
}
int checkma2(cov_model *cov) {
//  bool skipchecks = GLOBAL.general.skipchecks;
  cov_model *next = cov->sub[0];  
//  double p;
  int err;

  if ((err = checkkappas(cov)) != NOERROR) return err;
    
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, NegDefType, cov->domown,
		     cov->isoown, SCALAR, ROLE_COV)) 
	!= NOERROR) return err;

  cov->semiseparatelast = false;
  cov->logspeed = 0.0;
 //  updatepref(cov, next);
  setbackward(cov, next);
  cov->mpp.maxheight = 1.0;

  return NOERROR;
}


#define M_M 0
void kappaM(int i, cov_model *cov, int *nr, int *nc){
  *nc = SIZE_NOT_DETERMINED;
  *nr = i < CovList[cov->nr].kappas  ? SIZE_NOT_DETERMINED : -1;
}

void M(cov_model *cov, double *Z, double *V) {
  // M C M^t
  cov_model *next = cov->sub[0];
  int 
    *ncol = cov->ncol,
    *nrow = cov->nrow;
  double *MZ, 
    alpha = 1.0,
    beta = 0.0,
    *M = cov->p[M_M];
//      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K, ALPHA,A,LDA,B,LDB, BETA,C,LDC)
// C := alpha*op( A )*op( B ) + beta*C,
// M : rows of op(A)
// N : cols of op(B)
// K : cols of op(A)= rows of op(B)
// LD-X : first dimension of X

  // printf("dim=%d %d %d", next->vdim, nrow[0], ncol[0]);


  if (next->vdim == 1) {
    F77_CALL(dgemm)("N", "T", nrow, nrow, ncol, Z, M, nrow, M, nrow, 
		    &beta, V, nrow); 
  } else {
    MZ = (double*) MALLOC(sizeof(double) * *nrow * *ncol);
    F77_CALL(dgemm)("N", "N", nrow, ncol, ncol, &alpha, M, nrow, Z, ncol, 
		    &beta, MZ, nrow); 
    F77_CALL(dgemm)("N", "T", nrow, nrow, ncol, &alpha, MZ, nrow, M, nrow, 
		    &beta, V, nrow);
    free(MZ);
  }
}
void Mstat(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int    ncol = cov->ncol[M_M];
  double *z = (double*) MALLOC(sizeof(double) * ncol * ncol);
  COV(x, next, z);  
  M(cov, z, v);
  free(z);
}
void Mnonstat(double *x, double *y, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int    ncol = cov->ncol[M_M];
  double *z = (double*) MALLOC(sizeof(double) * ncol * ncol);
  NONSTATCOV(x, y, next, z);  
  M(cov, z, v);
  free(z);
}
int checkM(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;

  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->vdim = cov->nrow[M_M]; // zwingend vor CHECK
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, cov->domown, 
		     cov->isoown, cov->ncol[M_M], ROLE_COV)) != NOERROR) 
    return err;
  setbackward(cov, next);
  //  if (cov->vdim != cov->nrow[M_M]) APMI(cov);
  cov->mpp.maxheight = RF_NAN;
  return NOERROR;
}

sortsofparam paramtype_M(int k, int row, int col) {
  return k<0 ? VARPARAM : (row==col) ? SDPARAM : SIGNEDSDPARAM;
}
  
void rangeM(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[M_M] = RF_NEGINF;
  range->max[M_M] = RF_INF;
  range->pmin[M_M] = -1e5;
  range->pmax[M_M] = 1e5;
  range->openmin[M_M] = true;
  range->openmax[M_M] = true; 
}



#define SCHUR_M 0
#define SCHUR_DIAG 1
#define SCHUR_RED 2
void kappaSchur(int i, cov_model *cov, int *nr, int *nc){
  double *M = cov->p[SCHUR_M];
  int vdim = cov->nrow[M != NULL ? SCHUR_M : SCHUR_DIAG]; 
  *nc = i == SCHUR_M ? vdim : 1;
  *nr = i == SCHUR_RED ? vdim * (vdim-1) / 2 
    : i < CovList[cov->nr].kappas ? vdim : -1;
}

void SchurMult(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v){
  double *M = cov->p[SCHUR_M];
  int i,
    vdim = cov->vdim;
  if (cov->p[SCHUR_M] != NULL) {
    int nrow2 = cov->nrow[SCHUR_M] * cov->nrow[SCHUR_M];
    for (i=0; i<nrow2; i++) v[i] *= M[i]; 
  } else {
    int j, k;
    double 
      *q = cov->q,
      *diag=cov->p[SCHUR_DIAG],
      *red=cov->p[SCHUR_RED];
    for (i=0; i<vdim; i++) q[i] = sqrt(diag[i]);
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
  }
}

void Schurstat(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  COV(x, next, v); 
  SchurMult(x, cov, v);
}
void Schurnonstat(double *x, double *y, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  NONSTATCOV(x, y, next, v);  
  SchurMult(x, cov, v);
}
int checkSchur(cov_model *cov) {
  cov_model *next = cov->sub[0];
  double *C,
    *M=cov->p[SCHUR_M],
    *diag=cov->p[SCHUR_DIAG],
    *red=cov->p[SCHUR_RED];
  int i,j, k, l,
    err = NOERROR, 
    vdim = cov->nrow[M != NULL ? SCHUR_M : SCHUR_DIAG],
    vdimP1 = vdim + 1,
    bytes = vdim * vdim * sizeof(double),
    *nrow = cov->nrow,
    *ncol = cov->ncol;

  cov->vdim =  vdim;
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, cov->domown, 
		     cov->isoown, nrow[SCHUR_M], ROLE_COV)) != NOERROR) 
    return err;
  setbackward(cov, next);

  if (M == NULL) {
    if (diag == NULL || red == NULL)
      SERR("either 'diag' and 'red' or 'M' must be given");
    for (i=0; i<vdim; i++) if (diag[i]<0) SERR("elements of 'diag' negative.");
    C = (double*) MALLOC(bytes);
    for (k=l=i=0; i<vdim; i++, l+=vdimP1) {
      for (j=0; j<vdim; j++, k++) {
	C[i + vdim * j] = C[i * vdim + j] = red[k];
      }
      C[l] = 1.0;
    }
    F77_CALL(dpofa)(C, ncol, ncol, &err); // C i s now cholesky
    if (err != 0)
      SERR2("%d x %d matrix M is not (strictly) positive definite",
	    nrow[SCHUR_M] * vdim, ncol[SCHUR_M] * vdim);   
    cov->q = (double*)  MALLOC(vdim * sizeof(double));
    cov->qlen = vdim;
  } else {
    if (diag != NULL || red != NULL)
      SERR("if 'M' is given, neither 'diag' nor 'red' might be given.")
    C = (double*) MALLOC(bytes);
    MEMCOPY(C, M, bytes);
    F77_CALL(dpofa)(C, ncol, ncol, &err); // C i s now cholesky
    if (err != 0)
      SERR2("%d x %d matrix M is not (strictly) positive definite",
	    nrow[SCHUR_M] * vdim, ncol[SCHUR_M] * vdim);
  }
    
  free(C);
  cov->mpp.maxheight = RF_NAN;

  /*

    Check symmetry??

 for (vdiag = i=0; i<vdim; i++, vdiag+=vdimP1) {
    if (rhored[vdiag] != 1.0) SERR("diagonal elements do not equal 1.0");
    if (value[i]  <  -GLOBAL.direct.svdtolerance) {
      SERR("'rhored' has negative eigenvalues");
    }
  }
  for (vdiag=i=0; i<vdim; i++, vdiag+=vdimP1) {
    for (j=i+1; j<vdim; j++) {
      if (rhored[ vdiag + j - i] != rhored[vdiag + vdim * (j-i)]) 
	SERR("'rhored' not symmetric");
    }
  }
  */


  return err;
}
  
void rangeSchur(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[SCHUR_M] = RF_NEGINF;
  range->max[SCHUR_M] = RF_INF;
  range->pmin[SCHUR_M] = -1e5;
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
void IdStat(double *x, cov_model *cov, double *v) {
  //  cov_model 
  //    *key = cov->key;
  // if (key != NULL) CovList[key->nr].cov(x, key, v);
  // else {
    cov_model *next = cov->sub[0];
    COV(x, next, v);
    // }
}
void IdNonStat(double *x, double *y, cov_model *cov, double *v){
  // cov_model 
  //    *key = cov->key;
  //  if (key != NULL) CovList[key->nr].nonstat_cov(x, y, key, v);
  // else {
    cov_model *next = cov->sub[0];
    NONSTATCOV(x, y, next, v);
    // }
}

int checkId(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;

  cov->vdim = cov->p[ID_VDIM] != NULL  ?
      ((int*) cov->p[ID_VDIM])[0] : SUBMODEL_DEP;
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, cov->domown,
		     cov->isoown, cov->vdim, cov->role)) != NOERROR) return err;
  if (cov->vdim == SUBMODEL_DEP) cov->vdim = next->vdim;
  cov->logspeed = next->logspeed;
  setbackward(cov, next);
   
  return NOERROR;
}

void DId(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  Abl1(x, next, v);
}

void DDId(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  Abl2(x, next, v);
}
void TBM2Id(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  TBM2CALL(x, next, v);
}
void IdInverse(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  INVERSE(x, next, v);
}

int initId(cov_model *cov, storage *S) {
  cov_model *next = cov->sub[0];
  return INIT(next, cov->mpp.moments, S);
}
void spectralId(cov_model *cov, storage *S, double *e) { 
  //  spectral_storage *s = &(S->Sspectral);
  cov_model *next = cov->sub[0];
  return SPECTRAL(next, S, e); // nicht nr
}
void coinitId(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0];
  CovList[next->nr].coinit(next, li); // nicht nr
}
void ieinitId(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0];
  CovList[next->nr].ieinit(next, li); // nicht nr
}
void rangeId(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[ID_VDIM] = 1.0;
  range->max[ID_VDIM] = RF_INF;
  range->pmin[ID_VDIM] = 1.0;
  range->pmax[ID_VDIM] = 10;
  range->openmin[ID_VDIM] = false;
  range->openmax[ID_VDIM] = true;
}


#define EXP_N 0
#define EXP_STANDARDISED 1
void Exp(double *x, cov_model *cov, double *v, int n, bool standardize){
  double v0, 
    s = 0.0, 
    w = 1.0;
  cov_model *next = cov->sub[0];
  int k,
    vdim = cov->vdim,
    vdim2q = vdim * vdim;

  COV(x, next, v);
  if (vdim == 1) {
    for (k=0; k<=n; k++) {
      s += w;
      w *= *v / (double) (k+1);
    }    
    *v = exp(*v) - s;
    if (standardize) {
      Exp(ZERO, cov, &v0, n, false);
      *v /= v0;
    }
  } else {
    int i;
    BUG;
    // fehlt die multiplication von links und rechts mit C(0)^{-1/2}
    for (i=0; i<vdim2q; i++) v[i] = exp(v[i]); 
  }
}

void Exp(double *x, cov_model *cov, double *v){	
  Exp(x, cov, v, ((int*)cov->p[EXP_N])[0], ((int*)cov->p[EXP_STANDARDISED])[0]);
}


void nonstatExp(double *x, double *y, cov_model *cov, double *v, int n, 
		bool standardised){
  cov_model *next = cov->sub[0];
  double v0, 
    s = 0.0, 
    w = 1.0;
  int k,
    vdim = cov->vdim,
    vdim2q = vdim * vdim;

  NONSTATCOV(x, y, next, v);
  if (vdim == 1) {
    for (k=0; k<=n; k++) {
      s += w;
      w *= *v / (double) (k+1);
    }    
    *v = exp(*v) - s;
    if (standardised) {
      nonstatExp(ZERO, ZERO, cov, &v0, n, false);
      *v /= v0;
    }
  } else {
    int i;
    BUG;
    // fehlt die multiplication von links und rechts mit C(0)^{-1/2}
    for (i=0; i<vdim2q; i++) v[i] = exp(v[i]); 
  }
}

void nonstatExp(double *x, double *y, cov_model *cov, double *v){	
  nonstatExp(x, y, cov, v, ((int*)cov->p[EXP_N])[0], 
	     ((int*)cov->p[EXP_STANDARDISED])[0]);
}

void DExp(double *x, cov_model *cov, double *v){
  double D;
  cov_model *next = cov->sub[0];
  int n = ((int*) cov->p[EXP_N])[0];
  assert(cov->vdim == 1);
  
  Abl1(x, next, &D);
  Exp(x, cov, v, n - 1, false);
  *v *= -D;

  if (((int*) cov->p[EXP_STANDARDISED])[0]) {
    double v0;
    Exp(ZERO, cov, &v0, n, false);
    *v /= v0;
  }
}

  // Hier fehlt die multidim variante (Hesse matrix)
  

void DDExp(double *x, cov_model *cov, double *v){
  double D, Abl2, w;
  cov_model *next = cov->sub[0];
  int n = ((int*) cov->p[EXP_N])[0];
  assert(cov->vdim == 1);

  Abl1(x, next, &D);
  Abl2(x, next, &Abl2);
  Exp(x, cov, v, n - 2, false);
  Exp(x, cov, &w, n - 1, false);

  *v =  D * D * *v + Abl2 * w; // Achtung + D2 nicht - D2, da D2 Ableitung einer Cov.

  if (((int*) cov->p[EXP_STANDARDISED])[0]) {
    double v0;
    Exp(ZERO, cov, &v0, n, false);
    *v /= v0;
  }
}

int checkExp(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;

  kdefault(cov, EXP_N, -1);
  if (!isPosDef(next->typus) && ((int*) cov->p[EXP_N])[0] != -1)
    SERR("for variograms only n=-1 allowed");
  kdefault(cov, EXP_STANDARDISED, 1);
 
  if ((err = CHECKPD2ND(next, cov->tsdim,  cov->xdimown, cov->isoown,
			SCALAR, ROLE_COV)) != NOERROR) return err;

  next->delflag = DEL_COV - 10;

  setbackward(cov, next);
  if (cov->vdim > 1 && ((int*) cov->p[EXP_N])[0] != -1)
    SERR("'n' must be '-1' in the multivariate case");
  if (cov->vdim > 1) SERR("multivariate case not programmed yet");
 
  if (next->domown == XONLY) {
    cov_fct *C = CovList + cov->nr;
    cov->pref[CircEmbed] = C->pref[CircEmbed];
    cov->pref[Direct] = C->pref[Direct];
    cov->pref[Sequential] = C->pref[Sequential];    
    if (!isNegDef(cov->typus))
      SERR1("negative definite function expected -- got '%s'",
	    TYPENAMES[cov->typus]);
  } else {
    if (!isPosDef(cov)) SERR1("positive definite function expected -- got '%s'",
			      TYPENAMES[cov->typus]);
 
  }
  cov->mpp.maxheight = isNegDef(next->typus) && !isPosDef(next->typus)
    ? 1.0 : RF_NAN;

  cov->monotone = 
    (isBernstein(next)) ? NORMAL_MIXTURE : 
    isMonotone(next->monotone) ? MONOTONE : NOT_MONOTONE;
  cov->logspeed = 0.0;
  return NOERROR;
}


void rangeExp(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[EXP_N] = -1;
  range->max[EXP_N] = RF_INF;
  range->pmin[EXP_N] = -1;
  range->pmax[EXP_N] = 5;
  range->openmin[EXP_N] = false;
  range->openmax[EXP_N] = true;

  range->min[EXP_STANDARDISED] = 0;
  range->max[EXP_STANDARDISED] = 1;
  range->pmin[EXP_STANDARDISED] = 0;
  range->pmax[EXP_STANDARDISED] = 1;
  range->openmin[EXP_STANDARDISED] = false;
  range->openmax[EXP_STANDARDISED] = false;
}



#define POW_ALPHA 0
void Pow(double *x, cov_model *cov, double *v){
  double v0, v1,
    alpha = cov->p[POW_ALPHA][0];
  cov_model *next = cov->sub[0];
 
  COV(ZERO, next, &v0);
  COV(x, next, &v1);
  *v = pow(v0, alpha) - pow(v0 - v1, alpha); 
  //  print("Pow %f %f %f v=%f\n", v0, v1, alpha, *v);
}

void DPow(double *x, cov_model *cov, double *v){
  double v0, v1, gamma,
    alpha = cov->p[POW_ALPHA][0];
  cov_model *next = cov->sub[0];

  Abl1(x, next, v);
  if (alpha == 1.0) return;
  COV(ZERO, next, &v0);
  COV(x, next, &v1);
  gamma =  v0 - v1;
  *v *= - alpha * pow(gamma, alpha -1.0);
  // Achtung  "-" , da covarianzfunktion angenommen
}

void DDPow(double *x, cov_model *cov, double *v){
  double D, v0, v1, gamma,
    alpha = cov->p[POW_ALPHA][0];
  cov_model *next = cov->sub[0];

  Abl2(x, next, v);
  if (alpha == 1.0) return;
  Abl1(x, next, &D);
  COV(ZERO, next, &v0);
  COV(x, next, &v1);
  gamma =  v0 - v1;
  *v *= - alpha *  pow(gamma, alpha - 2.0) * ((alpha - 1.0) * D + gamma * (*v));
  // Achtung  "-" , da covarianzfunktion angenommen
}

void InversePow(double *x, cov_model *cov, double *v) {
  // invPow only defined for submodel having variance 1 !!!!
  cov_model *next = cov->sub[0];
  double alpha = cov->p[POW_ALPHA][0];
  
  COV(x, next, v);
  double y = 1.0 - *v;
  if (y < 0.0 || y > 1.0) {
    if (y > -1e-14 && y < 0.0)  y=0.0;
    else if (y < 1.0 + 1e-14)  y=1.0;
    else {
      //PRINTF("covariance value %e (1 - %e = %e) at %e\n", *v, *v, 1.0-*v, *x);
      ERR("invPow valid only for non-negative covariance models with variance 1");
    }
  }
  
  *v = 1.0 - pow(y, 1.0 / alpha);
  //  print("%f %f %f\n", y, *v, cov->p[POW_ALPHA][0]);
  //  print("Inv %f %f %f v=%f\n", *x,  cov->p[POW_ALPHA][0],y, *v);
}

int checkPow(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;
  //print("OK %d %d\n", cov->domown, cov->isoown);
  if ((err = checkkappas(cov)) != NOERROR) return err;
  if (!isNegDef(cov) || cov->domown != XONLY) return ERRORSTATVARIO;
  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, cov->domown, 
		     cov->isoown, SCALAR, ROLE_COV)) != NOERROR){
    // print("OK %d\n",err);
    return err;  
  }
  // print("xOK %d\n",err);
  setbackward(cov, next);      
  cov->mpp.maxheight = RF_NAN;
   cov->monotone = isMonotone(next->monotone) ? MONOTONE : NOT_MONOTONE;
  
  return NOERROR;
}

void rangePow(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[POW_ALPHA] = 0.0;
  range->max[POW_ALPHA] = 1.0;
  range->pmin[POW_ALPHA] = 0.01;
  range->pmax[POW_ALPHA] = 1.0;
  range->openmin[POW_ALPHA] = true;
  range->openmax[POW_ALPHA] = false;
}





/* qam */
#define QAM_THETA 0
void kappaqam(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i < CovList[cov->nr].kappas  ? cov->nsub-1 : -1;
}


void qam(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int i, 
    nsub = cov->nsub;
  double sum, s, w,
    *theta = cov->p[QAM_THETA];
  
  sum = 0.0;
  for (i=1; i<nsub; i++) {
    cov_model * sub = cov->sub[i];
    COV(x, sub, &s);
    INVERSE(&s, next, &w);
    sum += theta[i - 1] * w * w;
  }

  sum = sqrt(sum);
  COV(&sum, next, v); 
}

int checkqam(cov_model *cov) {
  cov_model *next = cov->sub[0],
    *sub;
  int i, err, 
    nsub = cov->nsub,
    nsubM1 = nsub - 1;
  double v, sum;
      
  if ((err = checkkappas(cov)) != NOERROR) return err;
  
  //  cov->monotone = NORMAL_MIXTURE;
  sum = 0.0;
  for (i=0; i<nsubM1; i++) {
    sum += cov->p[QAM_THETA][i];
  }
  if (fabs(sum - 1.0) > 1e-14) SERR("theta must sum up to 1");    

  if ((err = CHECK(next, 1,  1, PosDefType, cov->domown, cov->isoown, 
		     SCALAR, ROLE_COV)) != NOERROR)
    return err;
  if (!isNormalMixture(next->monotone))
    SERR("phi is not a normal mixture");   
 
  for (i=1; i<nsub; i++) {
    sub = cov->sub[i];
    if ((err = CHECK(sub, cov->tsdim, cov->tsdim, PosDefType,
		       cov->domown, cov->isoown,  SCALAR, ROLE_COV)) != NOERROR)
      return err;
    COV(ZERO, sub, &v);
    if (v != 1.0) SERR("unit variance required");
    setbackward(cov, sub);
  } 

  INVERSE(ZERO, next, &v);
  if (!ISNA(v)) SERR1("inverse function of '%s' unknown", NICK(next));

  cov->logspeed = 0.0;   
  return NOERROR;
}

sortsofparam paramtype_qam(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  return k==0 ? CRITICALPARAM : ANYPARAM;
}

void rangeqam(cov_model *cov, range_type *range){
  double pmax = 2.0 / (double)  (cov->nsub-1);
  range->min[QAM_THETA] = 0.0;
  range->max[QAM_THETA] = 1.0;
  range->pmin[QAM_THETA] = 0.0;
  range->pmax[QAM_THETA] = pmax;
  range->openmin[QAM_THETA] = false;
  range->openmax[QAM_THETA] = false;
}



/* mqam */
void kappamqam(int i, cov_model *cov, int *nr, int *nc) {
  int nsub = cov->nsub - 1;
  *nc = (i==QAM_THETA) ? 1 : -1;
  *nr = i == QAM_THETA ? nsub : -1;
  
  // print("kappa %d nsub=%d %d %d\n", i, nsub, *nc, *nr);

}
void mqam(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int i, j, k, l,
    vdim = cov->vdim,
    vdimP1 = vdim + 1;
  double s0,
    *theta = cov->p[QAM_THETA],
    s[MAXSUB];
  
  for (i=0; i<vdim; i++) {
    cov_model * sub = cov->sub[i+1];
    COV(x, sub, &s0);
    INVERSE(&s0, next, s + i);
    s[i] *= theta[i] * s[i];
  }
    
  for (j=0; j<vdim; j++) {
    l = k = j * vdimP1; // nicht vdim
    for (i=j; i<vdim; i++, k++, l+=vdim) {
      s0 = sqrt(s[i] + s[j]);
      COV(&s0, next, v + k);
      v[l] = v[k]; 
      // if (k != l) v[k] *= rho[k];
      // assert (rho[l] == rho[k]) ;
            
    }
  }

}

int checkmqam(cov_model *cov) {
  int err,
    nsub = cov->nsub,
    vdim = nsub - 1;     

  return ERRORNOTPROGRAMMED;

  if ((err = checkqam(cov)) != NOERROR) return err;
  cov->vdim = vdim;
  return NOERROR;
}

void rangemqam(cov_model *cov, range_type *range){
  //double pmax = 2.0 / (double)  (cov->nsub-1);

  rangeqam(cov, range);
}



/////////////// NATSC

void natsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double invscale, y;

  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = *x * invscale;
  // letzteres darf nur passieren wenn dim = 1!!
  COV(&y, next, v);
}

void Dnatsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim;
  double invscale, y;
 
  assert(CovList[next->nr].inverse != NULL);
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = *x * invscale;

  Abl1(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= invscale; 
}

void DDnatsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim;
  double invscale, y, invScSq;

  assert(CovList[next->nr].inverse != NULL);
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = *x * invscale;
  invScSq = invscale * invscale;

  Abl2(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= invScSq; 
}

void Inversenatsc(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim;
  double invscale;

  assert(CovList[next->nr].inverse != NULL);
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);    
  Abl1(x, next, v);
  for (i=0; i<vdimSq; i++) v[i] /= invscale; 

}

int checknatsc(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;


  assert(cov->nr == NATSC);

  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, cov->domown,
		   cov->isoown, SUBMODEL_DEP, ROLE_COV))
      != NOERROR) {
    return err;
  }

  if (next->domown == cov->domown && next->isoown == cov->isoown) {
    next->delflag = DEL_COV - 12;
  }

  if (CovList[next->nr].inverse == NULL) {
    sprintf(ERRORSTRING, "natural scaling is not defined for %s", 
	    NICK(next));
    return ERRORFAILED;
  }
 
  double invscale;
  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  if (R_FINITE(invscale) && invscale == INVERSE_UNKNOWN)
    SERR1("inverse function of '%s' unknown", NICK(next));


  cov->logspeed = 0.0;
  setbackward(cov, next);
  return NOERROR;
}


int initnatsc(cov_model *cov, storage *s){
  //  location_type *loc = Loc(cov);

  if (cov->role == ROLE_GAUSS) {
    cov_model *next = cov->sub[0];
    return INIT(next, cov->mpp.moments, s);
  }

  else if (cov->role == ROLE_BROWNRESNICK || cov->role == ROLE_SMITH ||
	   cov->role == ROLE_SCHLATHER || cov->role == ROLE_POISSON
	   || cov->role == ROLE_POISSON_GAUSS) {
   
    SERR("natsc for max-stable processes and poisson process not programmed yet");

  }

  else ILLEGAL_ROLE;

  return NOERROR;
}

void donatsc(cov_model *cov, storage *s){
  cov_model *next = cov->sub[0];
  DO(next, s);
}

void spectralnatsc(cov_model *cov, storage *S, double *e) {
  //  spectral_storage *s = &(S->Sspectral);
  cov_model *next = cov->sub[0];
  int d,
    dim = cov->tsdim;
  double invscale;

  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);  
  SPECTRAL(next, S, e);
  for (d=0; d<dim; d++) e[d] *=  invscale;
}

void coinitnatsc(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr; // not gatternr
  if ( C->coinit == NULL) 
    ERR("# cannot find coinit -- please inform author");
  C->coinit(next, li); // not gatternr
}

void ieinitnatsc(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr; // not gatternr
  if ( C->ieinit == NULL) // not gatternr
    ERR("# cannot find ieinit -- please inform author");
  C->ieinit(next, li); // not gatternr
}

void tbm2natsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double invscale, y;

  INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
  y = x[0] * invscale;
     
  TBM2CALL(&y, next, v);
}

////////////////////////////////////////////////////////////////////



void truncsupport(double *x, cov_model *cov, double *v){
  // truncsupport erwartet dass vorher $ kommt, sofern skalenmischung
  cov_model *next = cov->sub[0];
  int xdimown = cov->xdimown;
  double dist,
    radius = cov->p[TRUNC_RADIUS][0]; // default -1
  if (xdimown > 1) {
    int i;
    dist = 0.0;
    for (i=0; i<xdimown; i++) dist += x[i] * x[i];
    dist = sqrt(dist);
  } else dist = fabs(*x);

  if (radius>=0 && dist > radius) { *v=0.0; return; }
  FCTN(x, next, v);

  // printf("ts %f %f\n", x, v);
}

int checktruncsupport(cov_model *cov) {
  cov_model *next=cov->sub[0];
  int err,
    dim = cov->tsdim; // taken[MAX DIM],
 
  cov->maxdim=INFDIM;
  cov->monotone = isMonotone(next->monotone) ? MONOTONE : NOT_MONOTONE;

  if (cov->tsdim != cov->xdimown || cov->tsdim != cov->xdimprev)
    return ERRORDIM;

  if ((err = CHECK(next, dim, dim, ShapeType, cov->domown,
		     cov->isoown, SUBMODEL_DEP, cov->role)) != NOERROR) {
    //       print("error !!\n");
    return err;
  }

  next->delflag = DEL_COV - 20;
  setbackward(cov, next);

  return NOERROR;
}


void truncsupportInverse(double VARIABLE_IS_NOT_USED *x,
			 cov_model *cov, double *v){
  *v = cov->p[TRUNC_RADIUS][0];
}

void rangetruncsupport(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[TRUNC_RADIUS] = 0;      // < 0 : unset
  range->max[TRUNC_RADIUS] = RF_INF;
  range->pmin[TRUNC_RADIUS] = 0.00001;
  range->pmax[TRUNC_RADIUS] = 100.0;
  range->openmin[TRUNC_RADIUS] = true;
  range->openmax[TRUNC_RADIUS] = true;
}

int struct_truncsupport(cov_model *cov, cov_model **newmodel) {  
  int err;

  ASSERT_NEWMODEL_NOT_NULL;

  if (cov->role == ROLE_POISSON_GAUSS || cov->role == ROLE_MAXSTABLE || 
      cov->role == ROLE_POISSON) {
    if ((err = addUnifModel(cov, cov->p[TRUNC_RADIUS][0], newmodel)) != NOERROR)
      return err;
  } else ILLEGAL_ROLE_STRUCT;
 

  switch (cov->role) {
  case ROLE_POISSON_GAUSS :
    double invscale;
    addModel(newmodel, GAUSS);       
    addModel(newmodel, DOLLAR);
    kdefault(*newmodel, DSCALE, INVSQRTTWO);
    addModel(newmodel, TRUNCSUPPORT);
    InverseGauss(&GLOBAL.mpp.about_zero, cov, &invscale);
    kdefault(*newmodel, TRUNC_RADIUS, invscale);
    break;      
  case ROLE_POISSON : // optimierte density
    return addUnifModel(cov, 1.0, newmodel);
  case ROLE_MAXSTABLE :    
    return addUnifModel(cov, 1.0, newmodel);
  default : ILLEGAL_ROLE_STRUCT;
  }

  return NOERROR;
}

int init_truncsupport(cov_model *cov,  storage *s) {
  int err;
  if (cov->role == ROLE_BROWNRESNICK || cov->role == ROLE_SMITH ||
      cov->role == ROLE_SCHLATHER || cov->role == ROLE_POISSON
      || cov->role == ROLE_POISSON_GAUSS) {
    
    cov_model *next = cov->sub[0];
    //    double
    //      radius = cov->p[TRUNC_RADIUS][0]; // default -1
    
    if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) return err;
    //   if (radius>=0 && radius < cov->mpp.refradius) cov->mpp.refradius = radius;
    
    // Eplus, M2 are assumed to be still precise !!
    
    cov->mpp.maxheight = next->mpp.maxheight;

    return NOERROR;
  }

  else ILLEGAL_ROLE;
}


void do_truncsupport(cov_model *cov, storage *s) {
  //  mppinfotype *info = &(s->mppinfo);
  cov_model *next = cov->sub[0];
  //  double
  //    radius = cov->p[TRUNC_RADIUS][0]; // default -1
  
  DO(next, s);  
  cov->mpp.maxheight = next->mpp.maxheight;

  // if (radius>=0 && radius < info->radius) info->radius = radius;
}





void tbm3(double *x, cov_model *cov, double *v, double tbmdim){
  cov_model *next = cov->sub[TBMOP_COV];
  int i,
    vdim2 = cov->vdim * cov->vdim;
  double v1[MAXTBMVDIM * MAXTBMVDIM];
  COV(x, next, v); // x has dim 2, if turning planes
  //  print(" cov=%4.4f ", *v);
  if (x[0] != 0.0) {
    Abl1(x, next, v1);
    //  print(" D=%4.4f ", v1);
    for (i=0; i<vdim2; i++) v[i] += *x * v1[i] / tbmdim;
  }
}

typedef struct TBM2_integr{
  cov_model *cov;
  double *x;
  // bool layers;
} TBM2_integr;

void TBM2NumIntegrFct(double *u,  int n, void *ex) {
  int i;
  double z[2];
  TBM2_integr *info;
  info = (TBM2_integr*) ex;
  cov_model *cov=info->cov;
  double *x = info->x;
    
  for (i=0; i<n; i++) {
    z[0] = x[0] * sqrt(1 - u[i] * u[i]);
    tbm3(z, cov, u + i, 1.0);
  }
}

void tbm2num(double *x, cov_model *cov, double *v){
  TBM2_integr info;
#define MAXSUBDIVISIONS 100
  static double a = 0, 
      b = 1, 
      eps = 1e-10; // 1e-15
  static int maxsubdivisions = MAXSUBDIVISIONS, 
      lenw = 4 * MAXSUBDIVISIONS;
  double abserr, work[4 * MAXSUBDIVISIONS];
  int subdivisions, integralevaluations, err, iwork[MAXSUBDIVISIONS];

  info.cov = cov; // a bit strange: this is tbm2, but will call tbm3...
  //                 in TBM2NumIntegrFct ...
  info.x = x;
  Rdqags(TBM2NumIntegrFct, (void *) &info, &a, &b, &eps, &eps,
	 v, &abserr, &integralevaluations, &err,
	 &maxsubdivisions, &lenw, &subdivisions, iwork, work);
}


void tbm(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[TBMOP_COV];
  int 
    fulldim = ((int *) cov->p[TBMOP_FULLDIM])[0],
    tbmdim = ((int *) cov->p[TBMOP_TBMDIM])[0];
    //vdim2 = cov->vdim * cov->vdim;
  //  print("x=%4.4f t=%4.4f ", x[0], x[1]); 

  if (cov->role != ROLE_COV) COV(x, next, v)

  else if (fulldim == tbmdim + 2) tbm3(x, cov, v, (double) tbmdim);

  else if (fulldim == 2 && tbmdim == 1) {
    if (CovList[next->nr].tbm2 != NULL) TBM2CALL(x, next, v);
    else tbm2num(x, cov, v);
  }

  else XERR(ERRORTBMCOMBI);

}

void Dtbm(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[TBMOP_COV];
  int i,
    fulldim = ((int *) cov->p[TBMOP_FULLDIM])[0],
    vdim2 = cov->vdim * cov->vdim;
  double v1[MAXTBMVDIM * MAXTBMVDIM],
    tbmdim = (double) ((int *) cov->p[TBMOP_TBMDIM])[0],
    f = 1.0 + 1.0 / tbmdim;

  BUG;// to do : Taylor-Entwicklung im Ursprung

  if (fulldim == tbmdim + 2) {
    Abl1(x, next, v); // x has dim 2, if turning planes
    Abl2(x, next, v1);
    for (i=0; i<vdim2; i++) {
      v[i] = v[i] * f + *x * v1[i] / tbmdim;
    }
  }

  else XERR(ERRORTBMCOMBI);
}


int checktbmop(cov_model *cov) {
  cov_model *next=cov->sub[TBMOP_COV];
  tbm_param *gp  = &(GLOBAL.tbm);
  int err; 

  kdefault(cov, TBMOP_FULLDIM, cov->p[TBMOP_TBMDIM] == NULL || gp->tbmdim >= 0
	   ? gp->fulldim : ((int*)cov->p[TBMOP_TBMDIM])[0] - gp->tbmdim);
  kdefault(cov, TBMOP_TBMDIM, gp->tbmdim > 0 
	   ? gp->tbmdim : ((int *) cov->p[TBMOP_FULLDIM])[0] + gp->tbmdim);
  kdefault(cov, TBMOP_LAYERS, gp->layers);
  if ((err = checkkappas(cov)) != NOERROR)  return err;

  int 
    tbmdim = ((int *) cov->p[TBMOP_TBMDIM])[0],
    fulldim = ((int *) cov->p[TBMOP_FULLDIM])[0];
  double
    storedlayer = cov->p[TBMOP_LAYERS][0];
  bool layers = !ISNA(storedlayer) ? storedlayer :
    cov->xdimown == tbmdim + 1 && cov->isoown == SPACEISOTROPIC;
  if (tbmdim >= fulldim)
     SERR2("'reduceddim (=%d)' must be less than 'fulldim' (=%d)", 
	   tbmdim, fulldim);
  if (cov->tsdim > fulldim + layers) return ERRORWRONGDIM;
  //  printf("%d %d %d %d\n",  !ISNA_INT(storedlayer), storedlayer, cov->xdimown == tbmdim + 1 , cov->isoown == SPACEISOTROPIC );
  //printf("%d %d %d %d\n",cov->xdimown,tbmdim,layers, cov->isoown);
  if (cov->xdimown > tbmdim + layers) {
    //    APMI(cov);
    SERR("dimension of coordinates does not match reduced dimension of tbm");
  }

  if ((err = CHECK(next,  cov->tsdim, cov->xdimown, PosDefType, cov->domown,
		     cov->isoown, SUBMODEL_DEP, ROLE_COV)) != NOERROR) {
    // PMI(cov); printf("errr=%d %d %d\n", err, cov->tsdim, cov->xdimown);
    //XERR(err);
    return err;
  }
  if (next->pref[TBM] == PREF_NONE) return ERRORPREFNONE;

  if (cov->isoown != ISOTROPIC && cov->isoown != SPACEISOTROPIC) {
     return ERRORANISO;
  }
  if (!isNegDef(cov->typus) || cov->domown != XONLY) {
    return ERRORSTATVARIO;
  }


  cov->maxdim = 0;  
  setbackward(cov, next);
  cov->monotone=NOT_MONOTONE;
  cov->maxdim = fulldim + layers;  
  cov->rese_derivs = next->rese_derivs - 1;
  cov->finiterange = ((fulldim - tbmdim) % 2 == 0) && next->finiterange == true;

  if (cov->vdim > MAXTBMVDIM) 
    SERR2("vdim (%d) exceeds max. value of vdim in tbm3 (%d)", 
	  cov->vdim, MAXTBMVDIM);

  //  printf("gp=%f sto=%d %d %d\n", gp->layers, storedlayer, 
  //	 layers,   ((int*) cov->p[TBMOP_LAYERS])[0]);
  // APMI(cov);
  
  // only after being sure that the subsequent model does not cause
  // problems. So the time dimension should be fixed.
  // This is not absolutely safe programming, but should be OK.
  // But difficult to get around for MLE calls that do not allow 
  // for NAs values in integer variables.
  cov->p[TBMOP_LAYERS][0] = layers;
  return NOERROR;
}


void rangetbm_common(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range, bool tbmop ){ 
  int 
    TBMDIM = tbmop ? TBMOP_TBMDIM : TBM_TBMDIM,
    FULLDIM = tbmop ? TBMOP_FULLDIM : TBM_FULLDIM,
    LAYERS = tbmop ? TBMOP_LAYERS : TBM_LAYERS;
 
  range->min[FULLDIM] = 1.0;
  range->max[FULLDIM] = RF_INF;
  range->pmin[FULLDIM] = 1.0;
  range->pmax[FULLDIM] = 100;
  range->openmin[FULLDIM] = false;
  range->openmax[FULLDIM] = true;

  range->min[TBMDIM] = RF_NEGINF;
  range->max[TBMDIM] = RF_INF;
  range->pmin[TBMDIM] = RF_NEGINF;
  range->pmax[TBMDIM] = 100;
  range->openmin[TBMDIM] = false;
  range->openmax[TBMDIM] = true;

  range->min[LAYERS] = 0.0;
  range->max[LAYERS] = 1.0;
  range->pmin[LAYERS] = 0.0;
  range->pmax[LAYERS] = 1.0;
  range->openmin[LAYERS] = false;
  range->openmax[LAYERS] = false;
}

void rangetbmop(cov_model *cov, range_type *range){ 
  rangetbm_common(cov, range, true);
}





int set_stein_q(cov_model *next, double r, double d, double *q) {
  double C0, phi0, phi1, phi2, 
    zero = 0.0, 
    rP1 = r + 1.0,
    rM1 = r - 1.0, 
    dsq = d * d;

  COV(&zero, next, &C0);
  COV(&d, next, &phi0);
  Abl1(&d, next, &phi1); // derivative of
  phi1 *= d;             // scaled fctn at 1
  Abl2(&d, next, &phi2); // 2nd derivative of
  phi2 *= dsq;            // scaled fctn at 1
  q[LOCAL_R] =  r * d;   
  q[INTRINSIC_A2] = (phi2 - phi1) / (3.0 * r * rP1) ;
  q[INTRINSIC_B]  = (r == 1.0) ? 0.0 : q[INTRINSIC_A2] / (rM1 * dsq);
  q[INTRINSIC_A2] = (q[INTRINSIC_A2] - phi1 / 3.0 - phi2 / 6.0) / dsq;
  q[INTRINSIC_A0] = 0.5 * rM1 / rP1 * phi2 + phi1 / rP1 - phi0;


  //print("check: %f %e r=%f intr:%f %f  %f\n", 
  //	  phi2, phi1, r, q[INTRINSIC_B],  q[INTRINSIC_A0],  q[INTRINSIC_A2]);


//   assert(false);

  if ((q[INTRINSIC_B]  < 0.0) || (q[INTRINSIC_A2] < 0.0) ||
      (q[INTRINSIC_A0] + C0 < 0.0)) 
    return MSGLOCAL_INITINTRINSIC;

  return NOERROR;
}




int set_cutoff_q(cov_model *cov, double a, double d, double *q) {
  // auf modell ebene, d.h. co->sub[0] oder stein->sub[0]
  double phi0, phi1, a2; //, one=1.0;
  a2 = a * a;
  //  cov_model *neu = NULL;
  //int err;

  //  if ((err = covcpy(&neu, cov)) != NOERROR) return err;
  // CovList[cov->gatternr].cov(&d, neu, &phi0);
  // CovList[cov->gatternr].D(&d, neu, &phi1);
  //free(neu);
  COV(&d, cov, &phi0);
  Abl1(&d, cov, &phi1);
 
  phi1 *= d;

  //  print("dphi0ph1 %e %e %e\n", d, phi0, phi1);

  if (phi0 <= 0.0) return MSGLOCAL_SIGNPHI;
  if (phi1 >= 0.0) return MSGLOCAL_SIGNPHIFST;

//  print("check_co phi0=%f phi1=%f a=%f d=%f a2=%f\n", phi0, phi1, a, d, a2);
//  assert(false);
//  print("p %f\n", phi1);
//  print("a2 %f\n", a2);
//  print("p0 %f\n", phi0);
//  print("a %f\n", a);
//  print("d %f\n", d);
  q[CUTOFF_B] =//sound qconstant even if variance of submodel is not 1
    pow(- phi1 / (2.0 * a2 * phi0), 2.0 * a) * phi0 / pow(d, 2.0 * a2);
//  assert(false);
  q[CUTOFF_THEOR] = pow(1.0 - 2.0 * a2 * phi0 / phi1, 1.0 / a);
  q[LOCAL_R] = d * q[CUTOFF_THEOR];
  q[CUTOFF_ASQRTR] = pow(q[LOCAL_R], a);

//  print("phi0=%f %f %f %f\n", phi0, phi1, a, d);
//   print("%f\n", d); assert(false);
//  print("%f %f %f %f\n", q[CUTOFF_B], q[CUTOFF_THEOR], q[LOCAL_R], q[CUTOFF_ASQRTR]); //assert(false);
  return NOERROR;
}


int check_local(cov_model *cov,  
		 Methods method, int maxq, 
		 getlocalparam init, // next->coinit
		 set_local_q_type set_local) {
  location_type *loc = Loc(cov);
  int i, msg, 
    dim = cov->tsdim,
    err=NOERROR;
  //dim = cov->tsdim; // timespacedim -- needed ?;
  double  *q, q2[LOCAL_MAX],
    d=RF_NAN, **p = cov->p;
  cov_model *next = cov->sub[0];
  localinfotype li;
  //ce_param *gp  = &(GLOBAL.localce); // ok

  //  print("entering check local from %d:%s\n", cov->calling->nr,
  //	 CovList[cov->calling->nr].name);

  if ((err = CHECK(next, dim,  1,
		     method == CircEmbedCutoff ? PosDefType : NegDefType, 
		     cov->domown, cov->isoown, SCALAR, ROLE_COV)) != NOERROR)
    return err;
 


  // no setbackward ?!
  setbackward(cov, next); 
  if (next->pref[method] == PREF_NONE) return ERRORPREFNONE;

  //  PMI(cov);

  if (init == NULL){
    return ERRORUNKNOWNMETHOD;	  
  }

  //  PMI(cov);
  // assert(p[pLOC_DIAM] != NULL) ;

  
  // cov->variogram = false; 
  // no changing in pref by submodel !!
  if (cov->q != NULL) {
      free(cov->q);
      // ERR("q not NULL in check_local -- ask author");
  }
  cov->qlen = maxq;
  q = cov->q = (double*) CALLOC(maxq, sizeof(double));
  for (i = 0; i < maxq; i++) q2[i] = RF_NAN; // q2 will be copied to cov->q;
  //  {int i; print("%d %f\n", cov->qlen); for(i=0;i<9;i++)print("%f ", cov->q[i]); print("\n");  assert(false);}
  
  
  if (p[pLOC_DIAM] == NULL) {  
    double 
      diameter = GetDiameter(loc);
    if (PL>=PL_DETAILS) { LPRINT("diameter %f\n", diameter); }
    kdefault(cov, pLOC_DIAM, diameter);
  } else {
    d = p[pLOC_DIAM][0];
  }

  if (p[pLOC_A] == NULL) {
    if (CovList[next->nr].implemented[method]) {
      assert(init != NULL);

      init(next, &li);
      if (li.instances == 0) {
	SERR("parameter values do not allow for finding second parameter");
      }
      q[LOCAL_R] = R_PosInf;

      // print("%f\n", d); assert(false);

      msg = MSGLOCAL_FAILED;
      for (i=0; i<li.instances; i++) {
	// 	print("li %d %d %d %f %d\n", i, li.instances, li.msg[i], li.value[i],  msg);
	if (li.msg[i] <= msg) {
	  //	
	  err = set_local(next, li.value[i], d, q2);
	  //
	  if (err == NOERROR && (li.msg[i] < msg || q2[LOCAL_R] < q[LOCAL_R])){
	    MEMCOPY(q, q2, sizeof(double) * maxq);
	    msg = li.msg[i];
	    if (p[pLOC_A] == NULL) {
	      p[pLOC_A] = (double*) MALLOC(sizeof(double));
	      cov->ncol[pLOC_A] = cov->nrow[pLOC_A] = 1;
 	    }
	    cov->p[pLOC_A][0] = li.value[i]; // co: a; ie: rawr
	  }
	}
      }
      cov->q[LOCAL_MSG] = msg;

      //  print("msg %d %d\n", msg, MSGLOCAL_FAILED);

      if (msg == MSGLOCAL_FAILED) err = ERRORFAILED;
    } else {
	// PRINTF("%s\n", C->name);
      SERR("2nd parameter is neither given nor can be found automatically");
    }

    // APMI(cov);

  } else {
    if (cov->ncol[pLOC_A] != 1 || cov->nrow[pLOC_A] != 1) 
      SERR("`a' must be a scale");
    //   print("here 2\n");
    err = set_local(next, cov->p[pLOC_A][0], d, q2);
    MEMCOPY(q, q2, sizeof(double) * maxq);
  }

  cov->pref[CircEmbed] = 5;
  // APMI(cov);

// print("pdq %f %f %f %d\n", p[pLOC_A]==NULL ?RF_NAN :p[pLOC_A][0],d,q2, rr);
//
  //if (err != NOERROR) { X ERR(err)
    //   char Msg[255];
    //   err orMessage(err);
    //  err or(Msg);
  //}
  return err;
}


void co(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double  y=*x, *q=cov->q,
    diameter = cov->p[pLOC_DIAM][0],  
    a = cov->p[pLOC_A][0];
 
  assert(cov->role == ROLE_COV);

  if (y <= diameter) COV(x, next, v)
  else {
    *v = (y >= q[LOCAL_R]) ? 0.0
      : q[CUTOFF_B] * pow(q[CUTOFF_ASQRTR] - pow(y, a), 2.0 * a);
  }
}


int check_co(cov_model *cov) {
  cov_model *next = cov->sub[0];
  return check_local(cov, CircEmbedCutoff, CUTOFF_MAX, 
		     CovList[next->nr].coinit,
		     set_cutoff_q);
}

bool alternativeparam_co(cov_model VARIABLE_IS_NOT_USED *cov){
  return false;
}

void range_co(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){

  range->min[pLOC_DIAM] = 0.0; //  CUTOFF_DIAM
  range->max[pLOC_DIAM] = RF_INF;
  range->pmin[pLOC_DIAM] = 1e-10;
  range->pmax[pLOC_DIAM] = 1e10;
  range->openmin[pLOC_DIAM] = true;
  range->openmax[pLOC_DIAM] = true;

  range->min[pLOC_A] = 0.0; // cutoff_a
  range->max[pLOC_A] = RF_INF;
  range->pmin[pLOC_A] = 0.5;
  range->pmax[pLOC_A] = 2.0;
  range->openmin[pLOC_A] = true;
  range->openmax[pLOC_A] = true;
}


void Stein(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double y=*x, z, *q=cov->q,
    diameter = cov->p[pLOC_DIAM][0];

  // printf("diameter %f\n", diameter);

  assert(cov->role == ROLE_COV);

  if (y <= diameter) {
    COV(x, next, v);
    *v += q[INTRINSIC_A0] + q[INTRINSIC_A2] * y * y;

    //       print("%f %f %f %f %d %d \n ", 
    //	y, v[0], q[INTRINSIC_A0], q[INTRINSIC_A2], INTRINSIC_A0, INTRINSIC_A2);
  } else {
     z = q[LOCAL_R] - y;
     *v = (z <= 0.0) ? 0.0 : q[INTRINSIC_B] * z * z * z / y;
  }
}

int check_Stein(cov_model *cov)
{  
  cov_model *next = cov->sub[0]; // dito
  return check_local(cov, CircEmbedIntrinsic, INTRINSIC_MAX, 
		     CovList[next->nr].ieinit,
		     set_stein_q);
}  

bool alternativeparam_Stein(cov_model *cov) {
  cov->p[pLOC_A][0] *= 2.0;
  return true;
}

void range_Stein(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[pLOC_DIAM] = 0.0; 
  range->max[pLOC_DIAM] = RF_INF;
  range->pmin[pLOC_DIAM] = 0.01;
  range->pmax[pLOC_DIAM] = 100;
  range->openmin[pLOC_DIAM] = true;
  range->openmax[pLOC_DIAM] = true;

  range->min[pLOC_R] = 1.0; // stein_r
  range->max[pLOC_R] = RF_INF;
  range->pmin[pLOC_R] = 1.0;
  range->pmax[pLOC_R] = 20.0;
  range->openmin[pLOC_R] = false;
  range->openmax[pLOC_R] = true;
}



void strokorb(double *x, cov_model *cov, double *v) {
  // BrownResnick to strokorb Gaussian
  cov_model *next = cov->sub[0];
  int dim = cov->tsdim;
  double u;
      int idx = 0;
  
  u = 2 * *x;
  switch (dim) {
  case 1 :
    Abl1(&u, next, v);
    *v = - *v;
    //printf("u=%f %e\n", u, *v);
    break;
  case 3 :
     if (*x == 0.0) {
      while (idx < next->taylorN && (next->taylor[idx][TaylorPow] ==0.0 || 
				     next->taylor[idx][TaylorPow] ==1.0)) idx++;
      double p = next->taylor[idx][TaylorPow];
      if (idx >= next->taylorN || p > 3.0) BUG;
      // > 3 ist mathematisch vermutlich nicht moeglich, da bei strokorb
      // das submodel die entwicklung 1 - c x (oder rauher) hat
      *v = p < 3.0 ? RF_INF 
	: next->taylor[idx][TaylorConst] * p * (p - 1) / M_PI; // 3.0
    } else {
      Abl2(&u, next, v);
      *v /= (M_PI * *x);
    }
    break;
  default: BUG;
  }

  if (*v<0) {
     BUG;
  }
  
  assert(*v >= 0);
  assert(*x >= 0);
  //  assert(false);
}


int checkstrokorb(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int 
    dim = cov->tsdim,
    err = NOERROR;

  if ((err = CHECK(next, cov->tsdim, cov->xdimprev, TcfType,
		     cov->domown, cov->isoown,
		     SCALAR, ROLE_COV)) != NOERROR) return err;
  if (!isGneiting(next)) 
    SERR("member of the Gneiting-Schaback class as submodel needed");
  switch(dim) {
  case 1: 
    if (next->rese_derivs < 1) SERR("submodel must be once differentiable");
    break;
  case 3: 
    if (next->rese_derivs < 2) SERR("submodel must be twice differentiable");
    break;
  default:
    SERR("only dimensions 1 and 3 are allowed");
  }

  //  PMI(cov);

  if (!(hasMaxStableRole(cov) || hasNoRole(cov)))
    SERR1("'%s' may be used only as a shape function with max-stable field simulation", NICK(cov));

  if (next->tailN < 1)
    SERR2("%d members of the Taylor expansion at infinity of '%s', but at least order 1 required.", next->tailN, NICK(next));


  cov->taylorN = cov->tailN = 1;
  cov->tail[0][TaylorExpConst] = next->tail[0][TaylorExpConst];
  cov->tail[0][TaylorExpPow] = next->tail[0][TaylorExpPow];
  
  switch(dim) {
  case 1 : 
    cov->taylor[0][TaylorConst] =
      -next->taylor[1][TaylorConst] * next->taylor[1][TaylorPow];
    cov->taylor[0][TaylorPow] = next->taylor[1][TaylorPow] - 1.0;
    
    if (next->tail[0][TaylorExpPow] != 0.0) {
      if (next->tail[0][TaylorExpConst] == 0.0) BUG;
      cov->tail[0][TaylorConst] = next->tail[0][TaylorConst]
	* next->tail[0][TaylorExpConst] * next->tail[0][TaylorExpPow];
      cov->tail[0][TaylorPow] = 
	next->tail[0][TaylorPow] + next->tail[0][TaylorExpPow] - 1;
    } else {
      if (next->tail[0][TaylorExpConst] != 0.0) BUG;
      if (next->tail[0][TaylorPow] >= 0) {
	if (next->tail[0][TaylorConst] == 0) SERR("trivial submodel");
	//APMI(next);
	SERR1("'%s' is not integrable", NICK(next));
      }
      cov->tail[0][TaylorConst] = 
	- next->tail[0][TaylorConst] * next->tail[0][TaylorPow];
      cov->tail[0][TaylorPow] = next->tail[0][TaylorPow] - 1.0;
    }
    break;
  case 3 : 
    int idx;
    idx = 1;
    
    //PMI(next, -1);
    
    if (next->taylorN <= 1)
      SERR2("%d members of the Taylor expansion of '%s' known, but at least 2 members required.", next->taylorN, NICK(next));
    
    if (next->taylor[idx][TaylorPow] == 1.0) {	
      if (next->taylorN <= idx + 1)
	SERR3("%d members of the Taylor expansion of '%s' known, but at least %d members required.", next->taylorN, NICK(next), idx + 1);
	idx++;
    } else assert(next->taylor[idx][TaylorPow] < 1.0);
    
    cov->taylor[0][TaylorPow] = (next->taylor[idx][TaylorPow] - 2.0) - 1.0;
    cov->taylor[0][TaylorConst] = next->taylor[idx][TaylorConst] / M_PI * 
      next->taylor[idx][TaylorPow] * (next->taylor[idx][TaylorPow] - 1) *
      pow(2.0,  cov->taylor[0][TaylorPow]); // 
    

    if (next->tail[0][TaylorExpPow] != 0.0) {
      if (next->tail[0][TaylorExpConst] == 0.0) BUG;
      // Achtung ! Reihenfolge der Berechnung!
      double f = next->tail[0][TaylorExpConst] * next->tail[0][TaylorExpPow];
      cov->tail[0][TaylorPow] = next->tail[0][TaylorPow] + 
	 2.0 * (next->tail[0][TaylorExpPow] - 1.0);
      cov->tail[0][TaylorConst]= next->tail[0][TaylorConst] * f * f *
	pow(2.0, cov->tail[0][TaylorPow]);
      cov->tail[0][TaylorPow] -= 1.0;
      cov->tail[0][TaylorExpConst] *= pow(2.0, cov->tail[0][TaylorExpPow]);
    } else {
      // Achtung ! Reihenfolge der Berechnung!
      if (next->tail[0][TaylorExpConst] != 0.0) BUG;
      cov->tail[0][TaylorPow] = next->tail[0][TaylorPow] - 2.0;      
      cov->tail[0][TaylorConst] = next->tail[0][TaylorConst] / M_PI
	* next->tail[0][TaylorPow] * (next->tail[0][TaylorPow] - 1.0) 
	* pow(2.0, next->tail[0][TaylorPow] - 2.0);
      cov->tail[0][TaylorPow] -= 1.0;
    }
    break;
  default: BUG;
  }

  setbackward(cov, next);

  return NOERROR;
}



int init_strokorb(cov_model *cov,  storage VARIABLE_IS_NOT_USED *s) {
  if (cov->role == ROLE_MAXSTABLE) {
    
    //cov_model *next = cov->sub[0];
    //    double
    //      radius = cov->p[TRUNC_RADIUS][0]; // default -1
    
    // if ((err = INIT(next, 0, s)) != NOERROR) return err;

    //if (radius>=0 && radius < cov->mpp.refradius) cov->mpp.refradius = radius;
    // Eplus, M2 are assumed to be still precise !!
    
    cov->mpp.maxheight = RF_NAN;
    
  }
  
  else ILLEGAL_ROLE;

  cov->mpp.maxheight = 1.0; // all with be covered by pgs
  if (cov->mpp.moments >= 1) {
    cov->mpp.M[1] = cov->mpp.Mplus[1] = 1;     
  }

  return NOERROR;
}


void do_strokorb(cov_model VARIABLE_IS_NOT_USED *cov, storage VARIABLE_IS_NOT_USED *s) {
  assert(false);
}


////////////////// strokorbBall
/*
void strokorbBall(double *x, cov_model *cov, double *v) {
  // BrownResnick to strokorbBall Gaussian
  cov_model *key = cov->key;
  if (key == NULL) { *v = RF_NAN; return; }
  assert(key->nr >= DOLLAR && key->nr < LASTDOLLAR);
  int dim = cov->tsdim;  
  FCTN(x, key, v);
}
*/


#define STROKORBBALL_DIM 0 // Inner
int checkstrokorbBall(cov_model *cov) {
  cov_model
    //*key = cov->key,
    *next = cov->sub[0];
  int 
    dim = cov->tsdim,
    err = NOERROR;

  assert(cov->key == NULL);
  if ((err = CHECK(next, dim, cov->xdimprev, TcfType,
		   cov->domown, cov->isoown,
		   SCALAR, ROLE_COV)) != NOERROR) return err;
  if (!isGneiting(next)) 
    SERR("member of the Gneiting-Schaback class as submodel needed");

  switch(dim) {
  case 1: 
    if (next->rese_derivs < 2) SERR("submodel must be twice differentiable");
    break;
  case 3: 
    if (next->rese_derivs < 3) 
      SERR("submodel must be three times differentiable");
    break;
  default:
    SERR("only dimensions 1 and 3 are allowed");
  }
  
  
  //  PMI(cov);
  
  if (!(hasMaxStableRole(cov) || hasNoRole(cov)))
    SERR1("'%s' may be used only as a shape function with max-stable field simulation", NICK(cov));
  
  if (next->tailN < 1)
    SERR2("%d members of the Taylor expansion at infinity of '%s', but at least order 1 required.", next->tailN, NICK(next));
  
  if (next->taylorN < 2)
    SERR2("%d members of the Taylor expansion at infinity of '%s', but at least order 1 required.", next->tailN, NICK(next));

  setbackward(cov, next);

  return NOERROR;
}



int struct_strokorbBall(cov_model *cov, cov_model **newmodel) {  
  int 
    dim = cov->tsdim;
  
  ASSERT_NEWMODEL_NOT_NULL;
  assert(*newmodel == NULL);
  //PMI(cov, "strokorbBall");
  
  if (cov->role == ROLE_MAXSTABLE) {
    addModel(newmodel, BALL);
    addModel(newmodel, DOLLAR);
    covcpy((*newmodel)->kappasub + DSCALE, cov);
    // !! inverse scale gegenueber paper
    (*newmodel)->kappasub[DSCALE]->nr = STROKORB_BALL_INNER;
    kdefault((*newmodel)->kappasub[DSCALE], STROKORBBALL_DIM, dim);
    addModel( (*newmodel)->kappasub + DSCALE, RECTANGULAR); 
    kdefault((*newmodel)->kappasub[DSCALE], RECT_APPROX, false);
    kdefault((*newmodel)->kappasub[DSCALE], RECT_ONESIDED, true);
    addModel( (*newmodel)->kappasub + DSCALE, DOLLAR); 
    kdefault((*newmodel)->kappasub[DSCALE], DVAR, 0.5);
  } else ILLEGAL_ROLE_STRUCT;
  
  return NOERROR;
}

/*
int init_strokorbBall(cov_model *cov,  storage VARIABLE_IS_NOT_USED *s) {
  if (cov->role == ROLE_MAXSTABLE) {
    
    //cov_model *next = cov->sub[0];
    //    double
    //      radius = cov->p[TRUNC_RADIUS][0]; // default -1
    
    // if ((err = INIT(next, 0, s)) != NOERROR) return err;

    //if (radius>=0 && radius < cov->mpp.refradius) cov->mpp.refradius = radius;
    // Eplus, M2 are assumed to be still precise !!
    
    cov->mpp.maxheight = RF_NAN;
    
  }
  
  else ILLEGAL_ROLE;

  cov->mpp.maxheight = 1.0; // all with be covered by pgs
  if (cov->mpp.moments >= 1) {
    cov->mpp.M[1] = cov->mpp.Mplus[1] = 1;     
  }

  return NOERROR;
}


void do_strokorbBall(cov_model VARIABLE_IS_NOT_USED *cov, storage VARIABLE_IS_NOT_USED *s) {
  assert(false);
}

*/



void strokorbBallInner(double *x, cov_model *cov, double *v) {
  // only proportional to a density !!!
  cov_model *next = cov->sub[0];
  int dim = ((int *) cov->p[STROKORBBALL_DIM])[0];
  if (*x <= 0) {
    *v = 0.0;
    return;
  }
  switch (dim) {
  case 1 :
    // chi''(1/s)/s^3  --> \chi''(s)*s
    Abl2(x, next, v);
    *v *= intpow(*x, dim + 1);
    //printf("u=%f %e\n", u, *v);
    break;
  case 3 :
    double w;
    //[chi''(1/s) - \chi'''(1/s)/s] / (3s^3) --> [chi''(s)-chi'''(s) * s] * s/3
    Abl2(x, next, v);
    Abl3(x, next, &w);
    *v = (*v - w * *x) * intpow(*x, dim + 1) / 3.0;
    break;
  default: BUG;
  }
  if (*v<0) {
     BUG;
  }  
  //  assert(*v >= 0);
}


int check_strokorbBallInner(cov_model *cov) {
  cov_model *next = cov->key;
  int err, idx;
  double TP;
  
  ROLE_ASSERT(ROLE_MAXSTABLE);
  
  if ((err = checkkappas(cov)) != NOERROR) return err;
  if (cov->tsdim != 1) SERR("only dimension 1 allowed");
  if ((err = checkstrokorbBall(cov)) != NOERROR) return err;

  cov->taylorN = cov->tailN = 1;
  cov->tail[0][TaylorExpConst] = next->tail[0][TaylorExpConst];
  cov->tail[0][TaylorExpPow] = next->tail[0][TaylorExpPow];
  
  double tep = next->tail[0][TaylorExpPow],
    tp = next->tail[0][TaylorPow];
  
  switch(((int *) cov->p[STROKORBBALL_DIM])[0]) {
  case 1 :
    if (tep != 0.0) {
      cov->tail[0][TaylorConst] = next->tail[0][TaylorExpConst] * tep;
      cov->tail[0][TaylorConst] *= cov->tail[0][TaylorConst];
      cov->tail[0][TaylorPow] = tp + 2 * tep - 1.0;
    } else {
      cov->tail[0][TaylorConst] = tp * (tp -1.0);
      cov->tail[0][TaylorPow] = tp - 1.0; // !! nicht -2.0, da noch ein x drauf-multipliziert wird
    }
    
    idx = 1;
    if (next->taylor[1][TaylorPow] == (int) next->taylor[1][TaylorPow]) {
	// 2 sollte nie auftreten, laeuft aber gleich
      if (cov->taylorN >= 3) idx++;
    } 
    
    TP = next->taylor[idx][TaylorPow];
    cov->taylor[0][TaylorConst] =
      next->taylor[idx][TaylorConst] * TP * (TP - 1.0);
    cov->taylor[0][TaylorPow] = next->taylor[idx][TaylorPow] - 1.0;
    break;
  case 3 : 
    idx = 1;
    if (tep != 0.0) {
      double dummy = next->tail[0][TaylorExpConst] * tep;
      cov->tail[0][TaylorConst] = dummy * dummy * dummy;
      cov->tail[0][TaylorPow] = tp + 3 * tep - 2.0;
    } else {
      cov->tail[0][TaylorConst] = tp * (tp - 1.0) * (3.0 - tp);
      cov->tail[0][TaylorPow] = tp - 2.0; // !! nicht -2.0, da noch ein x drauf-multipliziert wird
    }
    cov->tail[0][TaylorConst] *= next->tail[0][TaylorConst] / 3.0;
    
    if (next->taylor[1][TaylorPow] == (int) next->taylor[1][TaylorPow]) {
      // 2 sollte nie auftreten, laeuft aber gleich
      if (cov->taylorN >= 3) idx++;
    } 
    
    TP = next->taylor[idx][TaylorPow];
    cov->taylor[0][TaylorConst] =
      next->taylor[idx][TaylorConst] * TP * (TP - 1.0) * (3.0 - TP) / 3.0;
    cov->taylor[0][TaylorPow] = next->taylor[idx][TaylorPow] - 2.0;
    break;
  default: BUG;
  }
  return NOERROR;
}

void range_strokorbBallInner(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[STROKORBBALL_DIM] = 1; 
  range->max[STROKORBBALL_DIM] = 3;
  range->pmin[STROKORBBALL_DIM] = 1;
  range->pmax[STROKORBBALL_DIM] = 3;
  range->openmin[STROKORBBALL_DIM] = false;
  range->openmax[STROKORBBALL_DIM] = false;
}




int checkstrokorbPoly(cov_model *cov) {
  cov_model
    //    *key = cov->key,
    *next = cov->sub[0];
  int 
    dim = cov->tsdim,
    err = NOERROR;

  assert(cov->key == NULL);
  if ((err = CHECK(next, dim, cov->xdimprev, TcfType,
		   cov->domown, cov->isoown,
		   SCALAR, ROLE_COV)) != NOERROR) return err;
  if (!isGneiting(next)) 
    SERR("member of the Gneiting-Schaback class as submodel needed");
  switch(dim) {
  case 1: 
    if (next->rese_derivs < 2) SERR("submodel must be twice differentiable");
    break;
  case 3: 
    if (next->rese_derivs < 3) 
      SERR("submodel must be three times differentiable");
    break;
  default:
    SERR("only dimensions 1 and 3 are allowed");
  }
  
  
  //  PMI(cov);
  
  if (!(hasMaxStableRole(cov) || hasNoRole(cov)))
    SERR1("'%s' may be used only as a shape function with max-stable field simulation", NICK(cov));
  
  if (next->tailN < 1)
    SERR2("%d members of the Taylor expansion at infinity of '%s', but at least order 1 required.", next->tailN, NICK(next));
  
  if (next->taylorN < 2)
    SERR2("%d members of the Taylor expansion at infinity of '%s', but at least order 1 required.", next->tailN, NICK(next));

  setbackward(cov, next);

  return NOERROR;
}



int struct_strokorbPoly(cov_model *cov, cov_model **newmodel) {  
  //  int err,
  //    dim = cov->tsdim;
  BUG;

  ASSERT_NEWMODEL_NOT_NULL;
  assert(*newmodel == NULL);
  //PMI(cov, "strokorbPoly");
  
  if (cov->role == ROLE_MAXSTABLE) {
    addModel(newmodel, POLYGON);
    addModel(newmodel, DOLLAR);
    covcpy((*newmodel)->kappasub + DSCALE, cov);
    (*newmodel)->kappasub[DSCALE]->nr = ARCSQRT_DISTR;
    //  addModel( (*newmodel)->kappasub + DSCALE, MULT_INVERSE); 
  } else ILLEGAL_ROLE_STRUCT;
  return NOERROR;
}




void mult_inverse(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  FCTN(x, next, v);
  *v = 1.0 / *v;
}


void mult_inverseNonstat(double *x, double *y, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  NONSTATCOV(x, y, next, v);
  *v = 1.0 / *v;
}


int checkmult_inverse(cov_model *cov) {
  cov_model
    *next = cov->sub[0];
   int err = NOERROR;
  if ((err = CHECK(next, cov->tsdim,  cov->xdimprev, ShapeType,
		   cov->domown, cov->isoown,
		   SUBMODEL_DEP, cov->role)) != NOERROR) return err;
  setbackward(cov, next);
  cov->mpp.maxheight = RF_NAN;
  return NOERROR;
}




////////////////////


void setparamStat(double *x, cov_model *cov, double *v){
  COV(x, cov->sub[SETPARAM_TO], v);
}

void setparamNonStat(double *x, double *y, cov_model *cov, double *v){
  NONSTATCOV(x, y, cov->sub[SETPARAM_TO], v);
}

void Dsetparam(double *x, cov_model *cov, double *v){
  Abl1(x, cov->sub[SETPARAM_TO], v);
}

void DDsetparam(double *x, cov_model *cov, double *v){
  Abl2(x, cov->sub[SETPARAM_TO], v);
}

void D3setparam(double *x, cov_model *cov, double *v){
  Abl3(x, cov->sub[SETPARAM_TO], v);
}
void D4setparam(double *x, cov_model *cov, double *v){
  Abl4(x, cov->sub[SETPARAM_TO], v);
}

void Inverse_setparam(double *v, cov_model *cov, double *x){
  cov_model *next = cov->sub[0];
  INVERSE(v, next, x);
}


int checksetparam(cov_model *cov) {
  cov_model *next = cov->sub[SETPARAM_TO];
  int err,
    dim = cov->tsdim,
    xdim = cov->xdimown,
    role = cov->role;
  Types type = cov->typus;
  domain_type dom = cov->domown;
  isotropy_type iso = cov->isoown;

  kdefault(cov, SET_PERFORMDO, 0);
  if ((err = CHECK(next, dim, xdim, type, dom, iso, 
		   SUBMODEL_DEP, role)) == NOERROR) return err;
  setbackward(cov, next);
   
  // ACHTUNG ! weder SETPARAM_FROM (da nur link) noch SETPARAM_SET (da 
  //           i.a. keine echten Modelle) werden ueberprueft!
 
  return NOERROR;
}


bool Typesetparam(Types required, cov_model *cov) {
  return TypeConsistency(required, cov->sub[SETPARAM_TO]);
}

void spectralsetparam(cov_model *cov, storage *s, double *e){
   SPECTRAL(cov->sub[SETPARAM_TO], s, e);  // nicht gatternr
}

//int initsetparam(cov_model *cov, storage *s){
//  err = INIT(cov->sub[SETPARAM_TO], cov->mpp.moments, s);
//}


void dosetparam(cov_model *cov, storage *s) {
  set_storage *X = cov->Sset;
  int performDo = ((int *) cov->p[SET_PERFORMDO])[0];
  assert(X != NULL);
  if (performDo < 0) DO(cov->sub[SETPARAM_TO], s); 
  if (X->from != NULL) {
    assert(X->set != NULL);
    X->set(cov->sub[0], X->from, X->variant);
  }
  if (performDo < 0) DO(cov->sub[SETPARAM_TO], s); 
}

void covmatrix_setparam(cov_model *cov, double *v) {
  cov_model *next = cov->sub[SETPARAM_TO];
  CovList[next->nr].covmatrix(next, v);
}

char iscovmatrix_setparam(cov_model *cov) {
  cov_model *next = cov->sub[SETPARAM_TO];
  return CovList[next->nr].is_covmatrix(next);
}

void range_setparam(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SET_PERFORMDO] = -1; 
  range->max[SET_PERFORMDO] = 1;
  range->pmin[SET_PERFORMDO] = -1;
  range->pmax[SET_PERFORMDO] = 1;
  range->openmin[SET_PERFORMDO] = false;
  range->openmax[SET_PERFORMDO] = false;
}

