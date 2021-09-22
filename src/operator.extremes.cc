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
#include <R_ext/Applic.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "operator.h"
#include "shape.h"


//////////////////////////////////////////////////////////////////////
// extremal gaussian 

void extrgauss(double *x, model *cov, double *v) {
  // BrownResnick to binary Gaussian
  model *next = cov->sub[0];
  double var, z;
  
  COV(ZERO(next), next, &var);
  COV(x, next, &z);
  *v = 1.0 - SQRT(0.5 * (1.0 - z / var));
}


int check_extrgauss(model *cov) {
  // to do extend to multivariate
  model
    *next = cov->sub[0];
  double v;
  int i,
    vdim = VDIM0,
    err = NOERROR;
  if (VDIM0 != VDIM1) BUG;
  if ((err = CHECK_PASSTYPE(next, PosDefType)) !=NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, cov->tsdim,  cov->xdimprev, PosDefType,
  //		     OWNDOM(0), OWNISO(0),
  //		     SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
  setbackward(cov, next);
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = 1.0;
  COV(ZERO(next), next, &v);
  if (v != 1.0) SERR("only correlation functions allowed");
  RETURN_NOERROR;
}


//////////////////////////////////////////////////////////////////////
// Brown Resnick

#define BR_FACTOR 0.25  //  KSH: '/ 4' 
#define BR_SEMI_FACTOR (2 * BR_FACTOR) // hier Semi-Variogram
void brownresnick(double *x, model *cov, double *v) {
  // BrownResnick to binary Gaussian
  model *next = cov->sub[0];
  double z;  

  COV(ZERO(next), next, &z);
  COV(x, next, v);
  *v = 2.0 * pnorm(SQRT((z - *v) * BR_SEMI_FACTOR), 0.0, 1.0, false, false);
}

void Dbrownresnick(double *x, model *cov, double *v) {
  // b = BR_SEMI_FACTOR
  // gamma(h) = C(0) - C(h) 
  // varphi(SQRT(gamma * b)) * (b C') / SQRT(gamma * b)
  // varphi density standard normal

  // BrownResnick to binary Gaussian
  model *next = cov->sub[0];
  double z, abl, s;
  
  //  

  assert((hasMaxStableFrame(cov) || hasGaussMethodFrame(cov) ||
	  hasAnyEvaluationFrame(cov)) && cov->taylorN > 1);

  if (cov->taylor[1][TaylorPow] == 0.0) { // const??
    *v = 0.0;
    return;
  }
  
  if (*x != 0.0) {
    COV(ZERO(next), next, &z);
    COV(x, next, v);
    
    assert(DefList[NEXTNR].D != NULL);
    assert(DefList[NEXTNR].D != ErrD);

    Abl1(x, next, &abl);   // Vorsicht: abl = -\gamma' 

    //    printf("x=%10g abl=%10g %10g z=%10g %10g\n", *x, abl, BR_SEMI_FACTOR, z, *v);
    //    PMI0(next);
    //    crash();

    abl *= BR_SEMI_FACTOR;
    s = SQRT((z - *v) * BR_SEMI_FACTOR); // SQRT(c * gamma)
    *v = dnorm(s, 0.0, 1.0, false) * abl / s; // -\varphi * \gamma' / sqrt\gamma
    //                                      =-2 \varphi * (-C') / (2 sqrt\gamma)
    //                              mit C= c_0 -\gamma
    //printf("v=%10g %10g %10g\n", *v, abl, s);
 
    assert(*v <= 0);

  } else {
    if (cov->taylor[1][TaylorPow] < 1.0) {     
      *v = RF_NEGINF;
    } else if (cov->taylor[1][TaylorPow] == 1.0) {
      *v = FABS(cov->taylor[1][TaylorConst]);
      assert(*v > 0.0);
    } else BUG;
  }
  // 2 * c * varphi(s) * gamma' / s
}


void DDbrownresnick(double *x, model *cov, double *v) {
  // D = \varphi * b C' / SQRT(\gamma b)

  // b = BR_SEMI_FACTOR
  // gamma(h) = C(0) - C(h) 
  // varphi(SQRT(gamma * b)) [(b C')^2 / [2 sqrt (gamma * b)]
  //                          +(b C'') / sqrt (gamma * b)
  //                          +1/2 * (b C')^2 / sqrt (gamma * b)^3 ]
 // varphi density standard normal

  model *next = cov->sub[0];
  double z, abl, abl2, s0, s;
  
  if (cov->taylor[1][TaylorPow] == 0.0) {
    *v = 0.0;
    return;
  }
   
  if (*x != 0.0) {
    COV(ZERO(next), next, &z);
    COV(x, next, v);
    Abl1(x, next, &abl);
    Abl2(x, next, &abl2);
    s0 = (z - *v) * BR_SEMI_FACTOR;
    s = SQRT(s0); // SQRT(c * gamma)
    abl  *= BR_SEMI_FACTOR;
    abl2 *= BR_SEMI_FACTOR;
    *v = dnorm(s, 0.0, 1.0, false) / s * (abl2 + abl * abl * 0.5 * (1/s0 + 1));

    assert(*v >= 0);

  } else {
    *v = cov->taylor[1][TaylorPow]==1 ? 0.0 : RF_INF;
  }
}


void D3brownresnick(double *x, model *cov, double *v) {
  // D = \varphi * b C' / SQRT(\gamma b)

  // b = BR_SEMI_FACTOR
  // gamma(h) = C(0) - C(h) 
  // varphi(SQRT(gamma * b)) [(b C')^3 / [4 sqrt (gamma * b)]
  //                         +3 (b C') (b C'')/ [2 sqrt (gamma * b)]
  //                         + (b C''') / [sqrt (gamma * b)]
  //                         + (b C')^2 / [2 sqrt (gamma * b)]
  //                         +3(b C') (b C'') / [2 sqrt (gamma * b)^3]
  //                         +3(b C')^3 / [4 sqrt (gamma * b)^5]

  model *next = cov->sub[0];
  double z, abl, abl2, abl3, s0, s;
  
   if (cov->taylor[1][TaylorPow] == 0.0) {
    *v = 0.0;
    return;
  }
   
  if (*x != 0.0) {
    COV(ZERO(next), next, &z);
    COV(x, next, v);
    Abl1(x, next, &abl);
    Abl2(x, next, &abl2);
    Abl3(x, next, &abl3);
    s0 = (z - *v) * BR_SEMI_FACTOR;
    s = SQRT(s0); // SQRT(c * gamma)
    abl  *= BR_SEMI_FACTOR;
    abl2 *= BR_SEMI_FACTOR;
    abl3 *= BR_SEMI_FACTOR;
    *v = dnorm(s, 0.0, 1.0, false) / s * 
      (abl3 + 
       1.5 * abl2 * abl * (1/s0 + 1) +
       abl * abl * abl * (0.25 + 0.5 / s0 + 0.75 / (s0 * s0)));

    // printf("br x=%10g v=%10e s=%10g abl=%10g s0=%10g abl2=%10g\n", *x, *v, s, abl, s0, abl2);
    //  assert(*v >= 0);

  } else {
    *v = cov->taylor[1][TaylorPow]==1 ? 0.0 : RF_NEGINF; // da zweite Ableitung im Ursprung Pol +\infty hat.
  }
}

int TaylorBrownresnick(model *cov) {
  model  
   *next = cov->sub[0];
  int idx = isnowPosDef(next); // Taylorentw. in 2 Glieder falls pos def.
 assert(idx == 0);

  //  assert(!idx);

  if (next->taylor[idx][TaylorPow] >= 2.0) {//for pos/neg def only == 2 possible
    cov->full_derivs = 1;
  } else cov->full_derivs = 0;
  cov->rese_derivs = MIN(3, next->rese_derivs); 
   
   // else if (next->taylor[idx][TaylorPow] == 2) {
  //  if (next->taylorN < 2 + idx) cov->rese_derivs = 0;
    // else if (cov->rese_derivs > 2) cov->rese_derivs = 2; 
  //  }
   
  if (next->taylorN >= 1 + idx &&  next->taylor[idx][TaylorConst] < 0.0) {
    // 2 \Phi(SQRT(gamma * br_f)) =
    //1+ 2 * phi(0) * SQRT(gamma * br_f) - 2 phi(0) / 6 *SQRT(gamma*br_f)^3
    //sei gamma(x) = c x^alpha + d x^beta, beta > alpha. Dann gilt
    // also  2 \Phi(SQRT(gamma * br_f)) ~
    // 1 + 2 * phi(0) * SQRT((c x^alpha + d x^beta)* br_f)
    //   - 2 * phi(0) / 6 * SQRT((c x^alpha + d x^beta) *br_f)^3
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
      g = SQRT(next_taylor_const  * BR_SEMI_FACTOR * 0.5 / M_PI);
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
	  cov->taylor[1][TaylorConst] += 2 * g * next_taylor_const * BR_SEMI_FACTOR;
	}
      } else cov->taylorN = 0;
    }
  } else cov->taylorN = 0;

  
  if (next->tailN >= 1) {
    cov->tailN = 1;    
    cov->tail[0][TaylorPow] = -0.5 * next->tail[0][TaylorPow];
    if (next->tail[0][TaylorPow] > 0) {  
      assert( next->tail[0][TaylorConst] < 0.0);
      double next_tail_const = - next->tail[0][TaylorConst];
      cov->tail[0][TaylorConst] = 
	2.0 / SQRT(2.0 * M_PI * BR_SEMI_FACTOR * next_tail_const);
      cov->tail[0][TaylorExpConst] = 0.5 * BR_SEMI_FACTOR * next_tail_const;
      cov->tail[0][TaylorExpPow] = next->tail[0][TaylorPow];
    } else {
      cov->tail[0][TaylorConst] = 
	2.0 / SQRT(2.0 * M_PI * BR_SEMI_FACTOR * next->tail[0][TaylorConst])
	* EXP(-0.5 * BR_SEMI_FACTOR * next->tail[0][TaylorConst]);
      cov->tail[0][TaylorPow] = cov->tail[0][TaylorExpConst] = 
	cov->tail[0][TaylorExpPow] = 0.0;
    }     
  } else cov->tailN = 0;

  if (cov->taylorN < 1) cov->rese_derivs = 0;

  RETURN_NOERROR;
}

int checkbrownresnick(model *cov) {
  // to do extend to multivariate
  model  
   *next = cov->sub[0];
  int i, err, 
    vdim = VDIM0;
  if (vdim != VDIM1) BUG;

   if ((err = CHECK_PASSTF(next, VariogramType, SUBMODEL_DEP, 
			   // hasMaxStableFrame(cov) ? BrMethodType  :
			   EvaluationType))
       // if ((err = CHECK(next, dim,  dim, VariogramType, OWNDOM(0), 
       //	   OWNISO(0), SUBMODEL_DEP, 
       //	   hasMaxStableFrame(cov) ? Ma xStableType : EvaluationType))
      != NOERROR) {
    RETURN_ERR(err);
  } 
  setbackward(cov, next);
  cov->monotone = isBernstein(next) ? GNEITING_MON : 
    isMonotone(next) ? MONOTONE : NOT_MONOTONE;

  if ((err = TaylorBrownresnick(cov)) != NOERROR) RETURN_ERR(err);

  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = 1.0;
  MEMCOPY(cov->pref, DefList[COVNR].pref, sizeof(pref_shorttype)); 
 
  RETURN_NOERROR;  
}



int struct_brownresnick(model *cov, model VARIABLE_IS_NOT_USED **newmodel) {  
  model *next = cov->sub[0];

  if (hasSmithFrame(cov)) {
    if (1 > next->taylorN  || 1 > next->tailN) 
      SERR2("frame '%.50s' not possible for submodel '%.50s'", 
	    TYPE_NAMES[cov->frame], NICK(next));

    BUG;

    // shape ist somit die Ableitung, falls d=1 und i.w. die 
    // zweifache Ableitung, falls d=3

    // hier ist auch Taylor zu setztn fuer den neuen Shape

  } else ILLEGAL_FRAME;


  RETURN_NOERROR;
}

int init_brownresnick(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int err;
  // model *next = cov->sub[0];
  if ((err = TaylorBrownresnick(cov)) != NOERROR) RETURN_ERR(err);
  RETURN_NOERROR;
}

void do_brownresnick(model *cov, gen_storage *s) {
  model *next = cov->sub[0];
  DO(next, s); // nicht gatternr
}


void BR2EG(double *x, model *cov, double *v) {
  // BrownResnick to binary Gaussian
  model *next = cov->sub[0];
  double z;
  COV(ZERO(next), next, &z);
  COV(x, next, v);
  z = 2.0 * pnorm(SQRT( (z - *v) * BR_SEMI_FACTOR), 0.0, 1.0, true, false) -1.0;
  *v = 1.0 - 2.0 * z * z;
}

int check_BR2EG(model *cov) {
  // to do extend to multivariate
  model  
   *next = cov->sub[0];
  double v, t, alpha;
  int err, i,
    vdim = VDIM0;
   if (VDIM0 != VDIM1) BUG;
  
  if ((err = CHECK_PASSTYPE(next, PosDefType)) !=NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType, 
  //		     OWNDOM(0), OWNISO(0), 
  //		     SCALAR, cov->frame)) != NOERROR)  RETURN_ERR(err);
  setbackward(cov, next);
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = 1.0;
  if (next->pref[Nothing] == PREF_NONE) RETURN_ERR(ERRORPREFNONECOV);

  // Erfc(x) = 2(1 - Phi(x * SQRT(2)))
  // Erf(x) = 2 Phi(x * SQRT(2)) - 1
  // Erf^{-1}(x) = Phi^{-1}( (y + 1) / 2 ) / SQRT(2) 

  // Sei c = 1-2 * Erf(SQRT(semivario / 4))^2 .
  // Also c = 1 - 2 [ 2 Phi(SQRT(semivario / 2)) - 1]^2
  // Umgekehrt semivario = 4 * { Erf^{-1} (sqrt[0.5 (1 - c)]) } ^2
  // mit c = 0 folgt SQRT(0.5 (1-c)) = 1 / SQRT(2)
  // semivario = 2 * {Phi^{-1}( (1 / SQRT(2) + 1) / 2 ) } ^2

  alpha = 0.0;
  COV(ZERO(next), next, &v);
  t = qnorm(0.5 * (1.0 + INVSQRTTWO), 0.0, 1.0, true, false);
  t *=  t / (BR_SEMI_FACTOR * (1.0 - alpha)); // 1/2 wegen Erf->qnorm
  if (v > t)
    SERR2("variance equals %10g, but must be at most 4(Erf^{-1}(1/2))^2 = %10g",
	  v, t);
  RETURN_NOERROR;  
}



void BR2BG(double *x, model *cov, double *v) {
  // BrownResnick to binary Gaussian
  model *next = cov->sub[0];
  double z;
  COV(ZERO(next), next, &z);
  COV(x, next, v); 
  z = 2.0 * pnorm(SQRT( (z - *v) * BR_SEMI_FACTOR), 0.0, 1.0, true, false) -1;
  *v = COS(M_PI * z);
}

int check_BR2BG(model *cov) {
  // to do extend to multivariate
  model *next = cov->sub[0];
  double v, t, alpha;
  int err, i,
    vdim = VDIM0;
  if (VDIM0 != VDIM1) BUG;
  //  if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType,
  //		     OWNDOM(0), OWNISO(0), 
  //		     SCALAR, cov->frame)) != NOERROR)  RETURN_ERR(err);
  if ((err = CHECK_PASSTF(next, PosDefType, vdim, //cov->frame
			  //hasMaxStableFrame(cov) ? BrMethodType :
			  EvaluationType
			  )) != NOERROR)
    RETURN_ERR(err);
   setbackward(cov, next);
   for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = 1.0;
   if (next->pref[Nothing] == PREF_NONE) RETURN_ERR(ERRORPREFNONECOV);

  // Erfc(x) = 2(1 - Phi(x * SQRT(2)))
  // Erf(x) = 2 Phi(x * SQRT(2)) - 1
  // Erf^{-1}(x) = Phi^{-1}( (y + 1) / 2 ) / SQRT(2) 


  // Sei c = COS(pi * Erf(SQRT(semivario / 4))) .
  // Also c = COS(pi * [2 * Phi(SQRT(semivario / 2)) - 1] )
  // Umgekehrt semivario = 2 * { Phi^{-1}(0.5 * [arcCOS( c ) / pi + 1]) }^2
  // mit c = 0 folgt arcCOS(c)/ pi + 1 = 3/2
  // semivario = 2 * { Phi^{-1}( 3 / 4) }^2
  
  COV(ZERO(next), next, &v);
  alpha = 0.0; // to do
  t = qnorm(0.75, 0.0, 1.0, false, false);
  t *= t / (BR_SEMI_FACTOR * (1.0 - alpha)); // 1/2 wegen Erf->qnorm
 
  if (v > t) { 
     SERR2("variance equals %10g, but must be at most 4(Erf^{-1}(1 / 2))^2 = %10g", 
	 v, t);
    }
  RETURN_NOERROR;   
}



void strokorb(double *x, model *cov, double *v) {
  // BrownResnick to strokorb Gaussian
  model *next = cov->sub[0];
  int dim = OWNLOGDIM(0);
  double u;
  int idx = 0;
  
  u = 2 * *x;
  switch (dim) {
  case 1 :
    Abl1(&u, next, v);
    *v = - *v;
    break;
  case 3 :
     if (*x == 0.0) {
      while (idx < next->taylorN && (next->taylor[idx][TaylorPow] ==0.0 || 
				     next->taylor[idx][TaylorPow] ==1.0)) idx++;
      if (idx >= next->taylorN) BUG;
      double p = next->taylor[idx][TaylorPow];
      if (p > 3.0) BUG;
      // > 3 ist mathematisch vermutlich nicht moeglich, da bei strokorb
      // das submodel die entwicklung 1 - c x (oder rauher) hat
      if (p == 3.0)
	*v = next->taylor[idx][TaylorConst] * p * (p - 1) * POW(2, p-2)/ M_PI;
      else *v = RF_INF;
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
  //   printf("x=%10g %10g ", *x, *v);    
  
  assert(*v >= 0);
  assert(*x >= 0);
  //  assert(false);
}


int TaylorStrokorb(model *cov) {
  model *next = cov->sub[0];
  int dim = OWNLOGDIM(0);
  
  if (next->tailN < 1)
    SERR2("%d members of the Taylor expansion at infinity of '%.50s', but at least order 1 required.", next->tailN, NICK(next));


  cov->taylorN = cov->tailN = 1;
  cov->tail[0][TaylorExpConst] = next->tail[0][TaylorExpConst];
  cov->tail[0][TaylorExpPow] = next->tail[0][TaylorExpPow];
  
  switch(dim) {
  case 1 : 
    cov->taylor[0][TaylorConst] =
      -next->taylor[1][TaylorConst] * next->taylor[1][TaylorPow];
    cov->taylor[0][TaylorPow] = next->taylor[1][TaylorPow] - 1.0;
    cov->taylorN = 1;

    cov->tailN = 1;
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
	SERR1("'%.50s' is not integrable", NICK(next));
      }
      cov->tail[0][TaylorConst] = 
	- next->tail[0][TaylorConst] * next->tail[0][TaylorPow];
      cov->tail[0][TaylorPow] = next->tail[0][TaylorPow] - 1.0;
    }
    break;
  case 3 : {
    int idx = 1;
    
    if (next->taylorN <= 1)
      SERR2("%d members of the Taylor expansion of '%.50s' known, but at least 2 members required.", next->taylorN, NICK(next));
    
    if (next->taylor[idx][TaylorPow] == 1.0) {	
      if (next->taylorN <= idx + 1)
	SERR3("%d members of the Taylor expansion of '%.50s' known, but at least %d members required.", next->taylorN, NICK(next), idx + 1);
	idx++;
    } else assert(next->taylor[idx][TaylorPow] < 1.0);
    
    cov->taylor[0][TaylorPow] = (next->taylor[idx][TaylorPow] - 2.0) - 1.0;
    cov->taylor[0][TaylorConst] = next->taylor[idx][TaylorConst] / M_PI * 
      next->taylor[idx][TaylorPow] * (next->taylor[idx][TaylorPow] - 1) *
      POW(2.0,  cov->taylor[0][TaylorPow]); // 
    cov->taylorN = 1;

    cov->tailN = 1;
    if (next->tail[0][TaylorExpPow] != 0.0) {
      if (next->tail[0][TaylorExpConst] == 0.0) BUG;
      // Achtung ! Reihenfolge der Berechnung!
      double f = next->tail[0][TaylorExpConst] * next->tail[0][TaylorExpPow];
      cov->tail[0][TaylorPow] = next->tail[0][TaylorPow] + 
	 2.0 * (next->tail[0][TaylorExpPow] - 1.0);
      cov->tail[0][TaylorConst]= next->tail[0][TaylorConst] * f * f *
	POW(2.0, cov->tail[0][TaylorPow]);
      cov->tail[0][TaylorPow] -= 1.0;
      cov->tail[0][TaylorExpConst] *= POW(2.0, cov->tail[0][TaylorExpPow]);
    } else {
      // Achtung ! Reihenfolge der Berechnung!
      if (next->tail[0][TaylorExpConst] != 0.0) BUG;
      cov->tail[0][TaylorPow] = next->tail[0][TaylorPow] - 2.0;      
      cov->tail[0][TaylorConst] = next->tail[0][TaylorConst] / M_PI
	* next->tail[0][TaylorPow] * (next->tail[0][TaylorPow] - 1.0) 
	* POW(2.0, next->tail[0][TaylorPow] - 2.0);
      cov->tail[0][TaylorPow] -= 1.0;
    }
  }
    break;
  default: BUG;
  }
  RETURN_NOERROR;
}
  
int checkstrokorb(model *cov) {
  model *next = cov->sub[0];
  int 
    dim = OWNLOGDIM(0),
    err = NOERROR;

 if ((err = CHECK_PASSTF(next, TcfType, SCALAR, EvaluationType)) != NOERROR)
   RETURN_ERR(err);
 //  if ((err = CHECK(next, cov->tsdim, cov->xdimprev, TcfType,
 //		     OWNDOM(0), OWNISO(0),
 //		     SCALAR, EvaluationType)) != NOERROR) RETURN_ERR(err);
 if (next->randomkappa) RETURN_ERR(ERRORRANDOMKAPPA);

  if (!isGneiting(next)) 
    SERR("member of the Gneiting-Schaback class as submodel needed");
  switch(dim) {
  case 1:
    if (next->rese_derivs < 1) {
      //APMI(cov);
      SERR("submodel must be once differentiable");
    }
    break;
  case 3: 
    if (next->rese_derivs < 2) SERR("submodel must be twice differentiable");
    break;
  default:
    SERR("only dimensions 1 and 3 are allowed");
  }

  if (false && !(hasMaxStableFrame(cov) || hasRandomFrame(cov))) {
    //APMI(cov);    crash();    crash();
     SERR1("'%.50s' may be used only as a shape function with max-stable field simulation", NICK(cov));
  }

  if ((err = TaylorStrokorb(cov)) != NOERROR) RETURN_ERR(err);
  setbackward(cov, next);

  RETURN_NOERROR;
}



int init_strokorb(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int err;
  if (hasSmithFrame(cov) || hasRandomFrame(cov)) {
    
    //model *next = cov->sub[0];
    //    double
    //      radius = P(TRUNC_RADIUS)[0]; // default -1
    
    // if ((err = INIT(next, 0, s)) != NOERROR) RETURN_ERR(err);

    //if (radius>=0 && radius < cov->mpp.refradius) cov->mpp.refradius = radius;
    // Eplus, M2 are assumed to be still precise !!
    
    assert(VDIM0 == 1 && VDIM1 == 1);
    cov->mpp.maxheights[0] = RF_NA;
    
  }
  
  else ILLEGAL_FRAME;

  /*
    if ((err = CHECK(cov->sub[0], cov->tsdim, cov->xdimprev, 
		   ShapeType, // 
		   OWNDOM(0), OWNISO(0),
		   SCALAR, 
		   M axStableType//  otherwise do will fail
		   )) != NOERROR) RETURN_ERR(err);
  */

  cov->mpp.maxheights[0] = 1.0; // all with be covered by pgs
  if (cov->mpp.moments >= 1) {
    cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1;     
  }

  if ((err = TaylorStrokorb(cov)) != NOERROR) RETURN_ERR(err);
  RETURN_NOERROR;
}


void do_strokorb(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  BUG;
}


////////////////// strokorbBall

#define STROKORBBALL_DIM 0 // Inner
int checkstrokorbBall(model *cov) {
  model
    //*key = cov->key,
    *next = cov->sub[0];
  int 
    dim = OWNLOGDIM(0),
    err = NOERROR;

  assert(cov->key == NULL);
  if ((err = CHECK_PASSTF(next, TcfType, SCALAR, EvaluationType)) !=NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, dim, cov->xdimprev, TcfType,
  //		   OWNDOM(0), OWNISO(0),
  //		   SCALAR, EvaluationType)) != NOERROR) RETURN_ERR(err);
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
  
  if (false && !(hasMaxStableFrame(cov) || hasRandomFrame(cov)))
    SERR1("'%.50s' may be used only as a shape function with max-stable field simulation", NICK(cov));
  
  if (next->tailN < 1)
    SERR2("%d members of the Taylor expansion at infinity of '%.50s' found, but at least 1 is required.", next->tailN, NICK(next));
  
  if (next->taylorN < 2)
    SERR2("%d members of the Taylor expansion of '%.50s' found, but at least 2 is required.", next->taylorN, NICK(next));

  setbackward(cov, next);

  RETURN_NOERROR;
}


void ScaleToVar(model *local, model *remote, 
		int VARIABLE_IS_NOT_USED variant) {
  // ACHTUNG!! from and to sind hier vertauscht!! 
  assert(MODELNR(local)==POWER_DOLLAR && MODELNR(remote)==LOC);
  //  int dim = local->tsdim;
  double scale = PARAM0(local, POWSCALE);
  // scale im entfernten model setzen
  PARAM(remote, LOC_SCALE)[0] = scale;
  remote->randomkappa = true;
}
  

    
int struct_strokorbBall(model *cov, model **newmodel) {  
  int  err,
    dim = OWNLOGDIM(0);
  
  ASSERT_NEWMODEL_NOT_NULL;
  
  if (hasSmithFrame(cov)) {
    addModel(newmodel, BALL, cov);    
    addModel(newmodel, POWER_DOLLAR);
    kdefault(*newmodel, POWSCALE, 1.0);    
    kdefault(*newmodel, POWPOWER, -dim);    
    kdefault(*newmodel, POWVAR,  1.0 / VolumeBall(dim, BALL_RADIUS));    
    model *pts=NULL, *scale=NULL;
    if ((err = covcpy(&pts, *newmodel)) != NOERROR) RETURN_ERR(err);
    

    if (DefList[COVNR].kappas < 2) {
      if ((err = covcpy(&scale, cov)) != NOERROR) RETURN_ERR(err);
      // !! inverse scale gegenueber paper
      SET_NR(scale, STROKORB_BALL_INNER);
      kdefault(scale, STROKORBBALL_DIM, dim);
      addModel(&scale, RECTANGULAR, *newmodel); 
      kdefault(scale, RECT_APPROX, false);
      kdefault(scale, RECT_ONESIDED, true);
      (*newmodel)->kappasub[POWSCALE] = scale;
    } else { // for testing only
      addModelKappa(*newmodel, POWSCALE, UNIF); 
      kdefault((*newmodel)->kappasub[POWSCALE], UNIF_MIN, P0(0));
      kdefault((*newmodel)->kappasub[POWSCALE], UNIF_MAX, P0(1));
    }
      
    addModel(&pts, RECTANGULAR, *newmodel);
    addModel(&pts, LOC, *newmodel);
    kdefault(pts, LOC_SCALE, 1.0);
    kdefault(pts, LOC_POWER, -dim);
    addModelKappa(pts, LOC_SCALE, NULL_MODEL); 
    kdefault(pts->kappasub[LOC_SCALE], NULL_TYPE, RandomType);

    addSetParam(newmodel, pts, ScaleToVar, true, 0);
    addModel(newmodel, ZHOU); // to do : unif better ?!
    (*newmodel)->sub[PGS_LOC] = pts;
    SET_CALLING(pts, *newmodel);

  } else ILLEGAL_FRAME_STRUCT;
  
  RETURN_NOERROR;
}


void rangestrokorbball(model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[0] = RF_NEGINF; 
  range->max[0] = RF_INF;
  range->pmin[0] = -4.0;
  range->pmax[0] = 4.0;
  range->openmin[0] = false;
  range->openmax[0] = false;

  range->min[1] = RF_NEGINF; 
  range->max[1] = RF_INF;
  range->pmin[1] = -4.0;
  range->pmax[1] = 4.0;
  range->openmin[1] = false;
  range->openmax[1] = false;
}
/*
int init_strokorbBall(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  if (cov->frame = = M axStableType) {
    
    //model *next = cov->sub[0];
    //    double
    //      radius = P(TRUNC_RADIUS)[0]; // default -1
    
    // if ((err = INIT(next, 0, s)) != NOERROR) RETURN_ERR(err);

    //if (radius>=0 && radius < cov->mpp.refradius) cov->mpp.refradius = radius;
    // Eplus, M2 are assumed to be still precise !!
    
    cov->mpp.maxheight = RF_NA;
    
  }
  
  else ILLEGAL_FRAME;

  cov->mpp.maxheight = 1.0; // all with be covered by pgs
  if (cov->mpp.moments >= 1) {
    cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1;     
  }

  RETURN_NOERROR;
}


void do_strokorbBall(model VARIABLE_IS_NOT_USED *cov, gen_storage VA// RIABLE_IS_NOT_USED *s) {
//   assert(false);
// }

// */



void strokorbBallInner(double *x, model *cov, double *v) {
  // only proportional to a density !!!
  model *next = cov->sub[0];

  int dim = COVNR != STROKORB_BALL_INNER || PisNULL(STROKORBBALL_DIM) 
    ? OWNLOGDIM(0) : P0INT(STROKORBBALL_DIM);
  if (*x <= 0) {
    *v = 0.0;
    return;
  }

  double y = 2 * *x;
  switch (dim) {
  case 1 :
    // \chi''(s)*s [ diameter ] -> 4\chi''(2 s)*s
    Abl2(&y, next, v);
    *v *= 2.0 * y;
    break;
  case 3 :
    double w;
    Abl2(&y, next, v);
    Abl3(&y, next, &w);
    *v = 2.0 * (*v - w * y) * y / 3.0;
    break;
  default: BUG;
  }
  if (*v<0) {
     BUG;
  }  
  //  assert(*v >= 0);  
}


int TaylorBall(model *cov) {
  model *next = cov->sub[0];
 
  if (next->tailN < 1 || next->taylorN < 2) 
    SERR1("taylor expansions of '%.50s' not programmed yet", NICK(next));
  double tep = next->tail[0][TaylorExpPow],
    tp = next->tail[0][TaylorPow];
 
  cov->taylorN = cov->tailN = 1;
  cov->tail[0][TaylorExpConst] = POW(2.0, tep) * next->tail[0][TaylorExpConst];
  cov->tail[0][TaylorExpPow] = tep;
  
  int idx = 1;
  if (next->taylor[1][TaylorPow] == (int) next->taylor[1][TaylorPow]) {
    assert(next->taylor[1][TaylorPow] == 1.0);
    // 2 sollte nie auftreten, laeuft aber gleich
    if (next->taylorN >= 3) idx++;
    else SERR1("%.50s does not have a long enough taylor development programmed", 
	      NICK(next));
  }
  double TP = next->taylor[idx][TaylorPow];

  switch(P0INT(STROKORBBALL_DIM)) {
  case 1 :
    if (tep != 0.0) {
      cov->tail[0][TaylorPow] = tp + 2.0 * (tep - 1.0) + 1.0;// !! +1, da noch ein x drauf-multipliziert wird
      cov->tail[0][TaylorConst] = next->tail[0][TaylorExpConst] * tep;
      cov->tail[0][TaylorConst] *= cov->tail[0][TaylorConst];
     } else {
      cov->tail[0][TaylorConst] = tp * (tp -1.0);
      cov->tail[0][TaylorPow] = tp - 1.0; // !! nicht -2.0, da noch ein x drauf-multipliziert wird
    }
         
    cov->taylor[0][TaylorConst] = TP * (TP - 1.0);
    cov->taylor[0][TaylorPow] = TP - 1.0;
    break;
  case 3 : 
    if (tep != 0.0) {
      double dummy = next->tail[0][TaylorExpConst] * tep;
      cov->tail[0][TaylorConst] = dummy * dummy * dummy / 3.0;
      cov->tail[0][TaylorPow] = tp + 3 * tep - 1.0;
    } else {
      cov->tail[0][TaylorConst] = tp * (tp - 1.0) * (3.0 - tp) / 3.0;
      cov->tail[0][TaylorPow] = tp - 1.0; // !! nicht -2.0, da noch ein x drauf-multipliziert wird
    }
        
    cov->taylor[0][TaylorConst] = TP * (TP - 1.0) * (3.0 - TP) / 3.0;
    cov->taylor[0][TaylorPow] = TP - 2.0;
    break;
  default: BUG;
  }

  cov->tail[0][TaylorConst] *= 2.0 * next->tail[0][TaylorConst] * 
    POW(2.0, cov->tail[0][TaylorPow]);
  cov->taylor[0][TaylorConst] *= 2.0 * next->taylor[idx][TaylorConst] * 
    POW(2.0, cov->taylor[0][TaylorPow]);

  RETURN_NOERROR;
}
  
int check_strokorbBallInner(model *cov) {
  model *next = cov->sub[0];
  int err;
  
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  if (OWNLOGDIM(0) != 1) SERR("only dimension 1 allowed");
  if ((err = checkstrokorbBall(cov)) != NOERROR) RETURN_ERR(err);
  switch(P0INT(STROKORBBALL_DIM)) {
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

  if ((err = TaylorBall(cov)) != NOERROR) RETURN_ERR(err);
  
  RETURN_NOERROR;
}

int init_strokorbBallInner(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  // 
  model *next = cov->sub[0];
  int err;

  if (next->randomkappa) RETURN_ERR(ERRORRANDOMKAPPA);// u.a. Taylor fehlt

  assert(VDIM0 == 1 && VDIM1 == 1);
  cov->mpp.maxheights[0] = 1.0; 
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1;     
  if (cov->mpp.moments >= 1) {
    cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1;     
  }

  if ((err = TaylorBall(cov)) != NOERROR) RETURN_ERR(err);
  RETURN_NOERROR;
}

void do_strokorbBallInner(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) {  
  //model *next = cov->sub[0];  
  //DO(next, s); // nicht gatternr
}


void range_strokorbBallInner(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[STROKORBBALL_DIM] = 1; 
  range->max[STROKORBBALL_DIM] = 3;
  range->pmin[STROKORBBALL_DIM] = 1;
  range->pmax[STROKORBBALL_DIM] = 3;
  range->openmin[STROKORBBALL_DIM] = false;
  range->openmax[STROKORBBALL_DIM] = false;
}





void poly2unif(model *local, model *remote, 
		int VARIABLE_IS_NOT_USED variant) {
  assert(MODELNR(local)==POLYGON && MODELNR(remote)==UNIF &&
	 local->Spolygon != NULL && local->Spolygon->P != NULL);
  polygon *P = local->Spolygon->P;
  assert(P->e != NULL);
  int d,
    dim = LOGDIM(SYSOF(local), 0);
  assert(dim == 2);

  for (d=0; d<dim; d++) {
    PARAM(remote, UNIF_MIN)[d] = P->box0[d];
    PARAM(remote, UNIF_MAX)[d] = P->box1[d];
  }
  remote->randomkappa = true;
}
  


void strokorbPoly(double *x, model *cov, double *v) {
  // only proportional to a density !!!
  model *next = cov->sub[0];
  COV(x, next, v);
}


int checkstrokorbPoly(model *cov) {
  model
    //    *key = cov->key,
    *next = cov->sub[0];
  int 
    dim = OWNLOGDIM(0),
    err = NOERROR;

  assert(cov->key == NULL);
  if ((err= CHECK_PASSTF(next, TcfType, SCALAR, EvaluationType)) != NOERROR) RETURN_ERR(err);
  //  if ((err = CHECK(next, dim, cov->xdimprev, TcfType,
  //		   OWNDOM(0), OWNISO(0),
  //		   SCALAR, EvaluationType)) != NOERROR) RETURN_ERR(err);
  if (!isGneiting(next)) 
    SERR("member of the Gneiting-Schaback class as submodel needed");
  
  if (dim != 2) SERR("only dimension 2 currently programmed");
  if (!(hasSmithFrame(cov))) {
    SERR1("'%.50s' may be used only as a shape function of a Smith field", NICK(cov));
  }

  setbackward(cov, next);

  RETURN_NOERROR;
}


int struct_strokorbPoly(model *cov, model **newmodel) {  
  model *sub = cov->sub[0];
  int 
    dim = OWNLOGDIM(0);
  double var = 1.0;
  model *pts=NULL, *shape=NULL;

  ASSERT_NEWMODEL_NOT_NULL;
    
  if (hasSmithFrame(cov)) {
    if (SUBNR != BROWNRESNICK) 
      SERR1("only tcf '%.50s' allowed", DefList[BROWNRESNICK].nick);
    sub = sub->sub[0];
    if (isDollar(sub)) {
      var = PARAM0(sub, DVAR);    
      sub = sub->sub[0];
    }
    if (SUBNR != BROWNIAN || PARAM0(sub, BROWN_ALPHA) != 1.0) {
      SERR2("Numerical inverse Laplace transform has not been implemented yet. Currently, only '%.50s' with parameter %.50s=1 is a valid submodel", DefList[BROWNIAN].nick,
	    DefList[BROWNIAN].kappanames[BROWN_ALPHA]);
    }
    
    addModel(&pts, UNIF, NULL, true);
    kdefault(pts, UNIF_NORMED, (int) false);
    PARAMALLOC(pts, UNIF_MIN, dim, 1);
    PARAMALLOC(pts, UNIF_MAX, dim, 1);
   
   
    addModel(&shape, POLYGON, NULL, true);
    addModelKappa(shape, POLYGON_BETA, ARCSQRT_DISTR);
    kdefault(shape->kappasub[POLYGON_BETA], ARCSQRT_SCALE, 1.0 / var);
    addSetParam(&shape, pts, poly2unif, true, 0);
   
    addModel(newmodel, ZHOU);
    kdefault(*newmodel, ZHOU_NORMED, false);
    kdefault(*newmodel, ZHOU_ISOTROPIC, false);
    SET_CALLING(shape, *newmodel);
    SET_CALLING(pts, *newmodel);
    (*newmodel)->sub[PGS_LOC] = pts;
    (*newmodel)->sub[PGS_FCT] = shape;

  } else ILLEGAL_FRAME_STRUCT;

  RETURN_NOERROR;
}

