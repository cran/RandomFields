/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
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
extern const char *coords[coordsN];

//#define LOG2 M_LN2

#define i11 0
#define i21 1
#define i22 2

#define epsilon 1e-15


#define GOLDENR 0.61803399
#define GOLDENC (1.0 -GOLDENR)
#define GOLDENTOL 1e-8

#define GOLDENSHIFT2(a, b, c) (a) = (b); (b) = (c);
#define GOLDENSHIFT3(a, b, c, d) (a) = (b); (b) = (c); (c) = (d);

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
    a    = EXP( -lambda * diff * diff);
    weights += a;
    sum += a * stuetz[i];  // or  a/stuetz[i]
  }
  return (double) (weights/sum); // and then   sum/weights       
}




/* BCW */
#define BCW_ALPHA 0
#define BCW_BETA 1
#define BCW_C 2
#define BCW_EPS 1e-7
#define BCW_TAYLOR_ZETA \
  (- LOG2 * (1.0 + 0.5 * zetalog2 * (1.0 + ONETHIRD * zetalog2)))
#define BCW_CAUCHY (POW(1.0 + POW(*x, alpha), zeta) - 1.0)
void bcw(double *x, model *cov, double *v){
  double alpha = P0(BCW_ALPHA), beta=P0(BCW_BETA),
    zeta = beta / alpha,
    absZeta = FABS(zeta);
  if (absZeta > BCW_EPS) 
    *v = BCW_CAUCHY / (1.0 - POW(2.0, zeta));
  else {
    double dewijs = LOG(1.0 + POW(*x, alpha)),
      zetadewijs = zeta * dewijs,
      zetalog2 = zeta * LOG2
      ;
    if (FABS(zetadewijs) > BCW_EPS)
      *v = BCW_CAUCHY / (zeta * BCW_TAYLOR_ZETA);
    else {
      *v = dewijs * (1.0 + 0.5 * zetadewijs * (1.0 + ONETHIRD *zetadewijs))
	/ BCW_TAYLOR_ZETA;
    }
  }
  if (!PisNULL(BCW_C)) *v += P0(BCW_C);
}


void Dbcw(double *x, model *cov, double *v){
  double ha,
    alpha = P0(BCW_ALPHA), 
    beta=P0(BCW_BETA),
    zeta=beta / alpha,
    y=*x;
  
  if (y ==0.0) {
    *v = (alpha > 1.0) ? 0.0 : (alpha < 1.0) ? RF_INF : alpha; 
  } else {
    ha = POW(y, alpha - 1.0);
    *v = alpha * ha * POW(1.0 + ha * y, zeta - 1.0);
  }
  
  if (FABS(zeta) > BCW_EPS) *v *= zeta / (1.0 - POW(2.0, zeta));
  else {
    double zetalog2 = zeta * LOG2;
    *v /= BCW_TAYLOR_ZETA;
  }
}


void DDbcw(double *x, model *cov, double *v){
  double 
  alpha = P0(BCW_ALPHA), 
  beta = P0(BCW_BETA),
  zeta = beta / alpha,
  y=*x;

  if (y == 0.0) {
    *v = alpha==2.0 ? alpha * (alpha - 1.0) :
      alpha==1.0 ? beta - 1.0 :
      alpha > 1.0 ? RF_INF : RF_NEGINF; 
  } else {
    double ha2 = POW(y, alpha - 2),
      ha = ha2 * y * y;
    *v = alpha * ha2 * (alpha - 1.0 + (beta - 1.0) * ha) *
      POW(1.0 + ha, zeta- 2.0);
  }
  if (FABS(zeta) > BCW_EPS) *v *= zeta / (1.0 - POW(2.0, zeta));
  else {
    double zetalog2 = zeta * LOG2;
    *v /= BCW_TAYLOR_ZETA;
  }
}

void D3bcw(double *x, model *cov, double *v){
  double
    alpha = P0(BCW_ALPHA),
    beta=P0(BCW_BETA),
    zeta = beta / alpha,
    absZeta = FABS(zeta),
     y=*x;

  if (y == 0.0) {
    *v =  RF_INF;
  } else {
    double
    ha3 = POW(y, alpha - 3),
    ha = ha3 * y * y * y;
    *v = alpha * POW(1.0 + ha, zeta - 3.0) * ha3 *
      (  (beta - 1.0) * (beta - 2.0) * ha * ha + 
          (alpha -1.0) * (3.0 * beta - alpha - 4.0) * ha +
	 (alpha - 2.0) * (alpha - 1.0) );
  }
  if (absZeta > BCW_EPS) *v *= zeta / (1.0 - POW(2.0, zeta));
  else {
    double zetalog2 = zeta * LOG2;
    *v /= BCW_TAYLOR_ZETA;
  }
}

 void D4bcw(double *x, model *cov, double *v){
  double
    alpha = P0(BCW_ALPHA),
    beta=P0(BCW_BETA),
    zeta = beta / alpha,
    absZeta = FABS(zeta),
  y = *x;

  if (y == 0.0) {
    *v =  RF_INF;
  } else {
  // olga please counter check, also for cauchy
   double
     ha4 = POW(y, alpha - 4),
     y2q = y * y,
     //alphaSq = alpha * alpha,
     a123 =  (alpha - 1.0) * (alpha - 2.0) * (alpha - 3.0),
     a1 = alpha - 1,
     coefha = -a1* (alpha * (4*alpha - 7*beta + 4) + 11 * beta - 18  ),
     coefha2  = a1 * (-4*alpha*beta + alpha*(alpha + 7) + 6*beta*beta
                      -22*beta +18),
     b1b2b3  = (beta - 1.0) * (beta - 2.0) * (beta - 3.0),
     ha = ha4 * y2q * y2q;

     *v = alpha * ha4 *
     ( a123 +
        coefha * ha +
       coefha2 * ha * ha +
       b1b2b3 * ha * ha * ha
       );
//   *v = alpha * ha / (y * y * y * y) * (alpha*alpha*alpha*(ha*ha - 4*ha +1)
//         +alpha*alpha*(-4*beta*ha*ha + 7*beta*ha + 6*ha*ha -6 )
//          -2*(3*beta*beta - 11*beta +9)*ha*ha
//           +alpha*((6*beta*beta -18*beta + 11)*ha*ha + (22-18*beta)*ha +11 )
//            +(beta*beta*beta - 6*beta*beta +11*beta - 6)*ha*ha*ha
//            + (11*beta - 18)*ha - 6)
//      * POW(1.0 + ha, beta / alpha - 4.0);
     }
  if (absZeta > BCW_EPS) *v *= zeta / (1.0 - POW(2.0, zeta));
  else {
    double zetalog2 = zeta * LOG2;
    *v /= BCW_TAYLOR_ZETA;
  }
}


void Inversebcw(double *x, model *cov, double *v) {  
  double 
    alpha = P0(BCW_ALPHA), 
    beta=P0(BCW_BETA),
    zeta = beta / alpha,
    y = *x;
  if (y == 0.0) {
    *v = beta < 0.0 ? RF_INF : 0.0; 
    return;
  }
  if (!PisNULL(BCW_C)) y = P0(BCW_C) - y;
  if (zeta != 0.0)
    *v = POW(POW(y * (POW(2.0, zeta) - 1.0) + 1.0, 1.0/zeta) - 1.0,
	       1.0/alpha); 
  else 
    *v =  POW(EXP(y * LOG2) - 1.0, 1.0 / alpha);   
}

int checkbcw(model *cov) {
  double
    alpha = P0(BCW_ALPHA), 
    beta=P0(BCW_BETA);
  if (OWNLOGDIM(0) > 2)
     cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = PREF_NONE;
  cov->logspeed = beta > 0 ? RF_INF : beta < 0.0 ? 0 : alpha * INVLOG2;
  RETURN_NOERROR;
}

void rangebcw(model *cov, range_type *range) {
  bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));
  bool posdef = isnowPosDef(cov) && PisNULL(BCW_C); // tricky programming
  range->min[BCW_ALPHA] = 0.0;
  range->max[BCW_ALPHA] = tcf ? 1.0 : 2.0;
  range->pmin[BCW_ALPHA] = 0.05;
  range->pmax[BCW_ALPHA] = range->max[BCW_ALPHA];
  range->openmin[BCW_ALPHA] = true;
  range->openmax[BCW_ALPHA] = false;

  range->min[BCW_BETA] = RF_NEGINF;
  range->max[BCW_BETA] = posdef ? 0.0 : 2.0;
  range->pmin[BCW_BETA] = - 10.0;
  range->pmax[BCW_BETA] = 2.0;
  range->openmin[BCW_BETA] = true;
  range->openmax[BCW_BETA] = posdef;

  range->min[BCW_C] = 0;
  range->max[BCW_C] = RF_INF;
  range->pmin[BCW_C] = 0;
  range->pmax[BCW_C] = 1000;
  range->openmin[BCW_C] = false;
  range->openmax[BCW_C] = true;
}



void coinitbcw(model *cov, localinfotype *li) {
  double
    beta=P0(BCW_BETA); 
  if (beta < 0) coinitgenCauchy(cov, li);
  else {
    if (beta == 0) coinitdewijsian(cov, li);
    else {
      double thres[2] = {0.5, 1.0},
	thresBeta[2] = {0.5, 1.0},
        alpha=P0(BCW_ALPHA);
      if (beta <= thresBeta[0] && alpha <= thres[0]) {
	li->instances = 2;
	li->value[0] = 0.5;
	li->value[1] = 1.0;
	li->msg[0] = li->msg[1] = MSGLOCAL_OK;
      } else {
	if (beta <= thresBeta[0] && alpha > thres[0] && alpha <= thres[1])  {
	  li->instances = 1;
	  li->value[0] = 1.0;
	  li->msg[0] =  MSGLOCAL_OK;
	} else {
	  if (beta>thresBeta[0] && beta <= thresBeta[1] && alpha <= thres[1]) {
	    li->instances = 1;
	    li->value[0] = 1.0;
	    li->msg[0] =  MSGLOCAL_OK;
	  } else {
	    if (beta > thresBeta[1] && alpha <= thres[0]  ) {
	      li->instances = 1;
	      li->instances = 1;
	      li->value[0] = CUTOFF_THIRD_CONDITION ;
	      li->msg[0] = MSGLOCAL_JUSTTRY;
	    } else {
	      li->instances = 0;
	    }
	  }
	}
      }
    }
  }
}

//void coinitbcw(model *cov, localinfotype *li) {
//  double
//    beta=P0(BCW_BETA);
//  if (beta < 0) coinitgenCauchy(cov, li);
//  else {
//    li->instances = 0;
//  }
//}
void ieinitbcw(model *cov, localinfotype *li) {
  if (P0(BCW_BETA) < 0) ieinitgenCauchy(cov, li);
  else {
    ieinitBrownian(cov, li); // to do: can nicht passen!!
    // todo: nachrechnen, dass OK!
  }
}



#define LOW_BESSELJ 1e-20
#define LOW_BESSELK 1e-20
#define BESSEL_NU 0
void Bessel(double *x, model *cov, double *v){
  // printf("bessel chec \n");
   /*
 GR 8.402:  Bessel(2 x SQRT(nu), nu) -> e^(-x^2)
 Differenz des k-ten koeffizienten in Taylor zu Gauss(ohne gemeinsame Faktoren):
 1 - 1 / [ (1+k/nu)(1+(k-1)/nu) ... (1 + 1/nu) ] \approx 0.5 k(k+1) / nu.
 Also ist der Abstand zu ungefaehr
 e^(-x^2) - Bessel(2 x SQRT(nu), nu) \approx F / nu
 fuer irgendein F > 0.
 Also
 e^(-x^2) - Bessel(2 x SQRT(nu), nu)  \approx
             nu0 nu^{-1} [ e^(-x^2) -  Bessel(2 x SQRT(nu0), nu0) ]
  */
  
  double nu = P0(BESSEL_NU),
    nuThres = nu <= BESSEL_NU_THRES ? nu : BESSEL_NU_THRES, 
    y=*x,
    bk[BESSEL_NU_THRES + 1L]; 
  if  (y <= LOW_BESSELJ) {*v = 1.0; return;}
  if (y == RF_INF)  {*v = 0.0; return;} // also for cosine i.e. nu=-1/2
  assert(cov->qlen != 0 && cov->q != NULL);

  //PMI(cov);
  //  printf("%10g %10g\n", nu, QVALUE);
  
  *v = QVALUE * POW(2.0 / y, nuThres) * bessel_j_ex(y, nuThres, bk);

  if (nu > BESSEL_NU_THRES) { // factor!=0.0 && 
    //  printf("UU \n");
    double w,
      g = BESSEL_NU_THRES / nu;
    y = 0.5 * *x / SQRT(nuThres);
    Gauss(&y, NULL, &w);
    *v = *v * g + (1.0 - g) * w;
  }
}
int initBessel(model *cov, gen_storage 
	       VARIABLE_IS_NOT_USED *s) {
  //  printf("init chec \n");
  assert(cov->q != NULL);
  double nu = P0(BESSEL_NU),
    nuThres = nu <= BESSEL_NU_THRES ? nu : BESSEL_NU_THRES;
   QVALUE = gammafn(nuThres+1.0);
   //  printf("intit end chec \n");
  ASSERT_GAUSS_METHOD(SpectralTBM);
   RETURN_NOERROR; 
}

int checkBessel(model *cov) {
  //  printf("chec \n");
  // Whenever TBM3Bessel exists, add further check against too small nu! 
  double nu = P0(BESSEL_NU);
  int i;
  double dim = (2.0 * P0(BESSEL_NU) + 2.0);

  for (i=0; i<= Nothing; i++)
    cov->pref[i] *= ((bool) ISNAN(nu)) || nu < BesselUpperB[i];
  if (OWNLOGDIM(0)>2) cov->pref[SpectralTBM] = PREF_NONE; // do to d > 2 !
  set_maxdim(OWN, 0, ISNAN(dim) || dim >= INFDIM ? INFDIM : (int) dim);
  if (cov->q == NULL) {
    EXTRA_Q;
    initBessel(cov, NULL);
  }
  //  printf("chec e\n");
  RETURN_NOERROR;
}
void spectralBessel(model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  double 
    nu =  P0(BESSEL_NU);
/* see Yaglom ! */
  // nu==0.0 ? 1.0 : // not allowed anymore;
	// other wise densityBessel (to define) will not work
  if (nu >= 0.0) {
    E12(s, OWNLOGDIM(0), nu > 0 ? SQRT(1.0 - POW(UNIFORM_RANDOM, 1/nu)) : 1, e);
  } else {
    double A;
    assert(OWNLOGDIM(0) == 1);
    if (nu == -0.5) A = 1.0;
    else { // siehe private/bessel.pdf, p. 6, Remarks
      // http://homepage.tudelft.nl/11r49/documents/wi4006/bessel.pdf
      while (true) {
	A = 1.0 - POW(UNIFORM_RANDOM, 1.0 / ( P0(BESSEL_NU) + 0.5));
	if (UNIFORM_RANDOM <= POW(1 + A, nu - 0.5)) break;
      }
    }
    E1(s, A, e);
  }
}
void rangeBessel(model *cov, range_type *range){
  range->min[BESSEL_NU] = 0.5 * ((double) OWNLOGDIM(0) - 2.0);
  range->max[BESSEL_NU] = RF_INF;
  range->pmin[BESSEL_NU] = 0.0001 + range->min[BESSEL_NU];
  range->pmax[BESSEL_NU] = range->pmin[BESSEL_NU] + 10.0;
  range->openmin[BESSEL_NU] = false;
  range->openmax[BESSEL_NU] = true;
}



/* circular model */
void circular(double *x, model VARIABLE_IS_NOT_USED  *cov, double *v) {
  double y = *x;
  *v = 0.0;
  if (y < 1.0) *v = 1.0 - (2.0 * (y * SQRT(1.0 - y * y) + ASIN(y))) * INVPI;
}
void Dcircular(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x * *x;
  *v = 0.0;
  if (y < 1.0) *v =-4 * INVPI * SQRT(1.0 - y);
}
int structCircSph(model *cov, model **newmodel, int dim) { 
  ASSERT_NEWMODEL_NOT_NULL;

  switch (cov->frame) {
  case PoissonGaussType :
    {
      addModel(newmodel, BALL, cov);      
      addModel(newmodel, DOLLAR);
      addModelKappa(*newmodel, DSCALE, SCALESPHERICAL);
      kdefault((*newmodel)->kappasub[DSCALE], SPHERIC_SPACEDIM,
	       (double) OWNLOGDIM(0));
      kdefault((*newmodel)->kappasub[DSCALE], SPHERIC_BALLDIM, (double) dim);
     } 
    break;
  case PoissonType : return addUnifModel(cov, 1.0, newmodel); // to do: felix
  case SmithType : return addUnifModel(cov, 1.0, newmodel);
  default:
    BUG;
  }
  RETURN_NOERROR;
}

int structcircular(model *cov, model **newmodel) {
  return structCircSph(cov, newmodel, 2);
}





/* coxgauss, cmp with nsst1 !! */
// see Gneiting.cc

/* cubic */
void cubic(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y=*x, y2=y * y;
  *v = 0.0;
  if (y < 1.0) *v = (1.0 + (((0.75 * y2 - 3.5) * y2 + 8.75) * y - 7) * y2);
}
void Dcubic(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) { 
  double y=*x, y2=y * y;
  *v =  0.0;
  if (y < 1.0) *v = y * (-14.0 + y * (26.25 + y2 * (-17.5 + 5.25 * y2)));
}

/* cutoff */
// see operator.cc

/* dagum */
#define DAGUM_BETA 0
#define DAGUM_GAMMA 1
#define DAGUM_BETAGAMMA 2



#define FIRST_INIT(init)				\
  int xerr;						\
  gen_storage s;					\
  gen_NULL(&s);						\
  s.check = true;					\
  if ((xerr=init(cov, &s)) != NOERROR) RETURN_ERR(xerr);


void dagum(double *x, model *cov, double *v){
  double gamma = P0(DAGUM_GAMMA), 
    beta=P0(DAGUM_BETA);
  *v = 1.0 - POW((1 + POW(*x, -beta)), -gamma/beta);
}
void Inversedagum(double *x, model *cov, double *v){ 
  double gamma = P0(DAGUM_GAMMA), 
    beta=P0(DAGUM_BETA);
  *v = RF_INF;
  if (*x != 0.0) *v = POW(POW(1.0 - *x, - beta / gamma ) - 1.0, 1.0 / beta);
} 
void Ddagum(double *x, model *cov, double *v){
  double y=*x, xd, 
    gamma = P0(DAGUM_GAMMA), 
    beta=P0(DAGUM_BETA);
  xd = POW(y, -beta);
  *v = -gamma * xd / y * POW(1 + xd, -gamma/ beta -1);
}

int checkdagum(model *cov){
  if (PisNULL(DAGUM_GAMMA) || PisNULL(DAGUM_BETA))
    SERR("parameters are not given all");
  double
    gamma = P0(DAGUM_GAMMA), 
    beta = P0(DAGUM_BETA);
  kdefault(cov, DAGUM_BETAGAMMA, beta/gamma);
  FIRST_INIT(initdagum);
 

  cov->monotone =  beta >= gamma ? MONOTONE 
    : gamma <= 1.0 ? COMPLETELY_MON
    : gamma <= 2.0 ? NORMAL_MIXTURE : MON_MISMATCH;
    
  RETURN_NOERROR;
}

int initdagum(model *cov, gen_storage *s) {  
  double
    gamma = P0(DAGUM_GAMMA), 
    beta = P0(DAGUM_BETA),
    betagamma = P0(DAGUM_BETAGAMMA);
  
  if (s->check) {
    bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));
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
  RETURN_NOERROR;
}

void rangedagum(model *cov, range_type *range){
  bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));
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
void dampedcosine(double *x, model *cov, double *v){
  double y = *x, lambda = P0(DC_LAMBDA);
  *v = (y == RF_INF) ? 0.0 : EXP(-y * lambda) * COS(y);
}
void logdampedcosine(double *x, model *cov, double *v, double *Sign){
  double 
    y = *x, 
    lambda = P0(DC_LAMBDA);
  if (y==RF_INF) {
    *v = RF_NEGINF;
    *Sign = 0.0;
  } else {
    double cosy=COS(y);
    *v = -y * lambda + LOG(FABS(cosy));
    *Sign = cosy > 0.0 ? 1.0 : cosy < 0.0 ? -1.0 : 0.0;
  }
}
void Inversedampedcosine(double *x, model *cov, double *v){ 
  Inverseexponential(x, cov, v);
} 
void Ddampedcosine(double *x, model *cov, double *v){
  double y = *x, lambda = P0(DC_LAMBDA);
  *v = - EXP(-lambda*y) * (lambda * COS(y) + SIN(y));
}
int checkdampedcosine(model *cov){
  int dim = INFDIM;
  if (!ISNAN(P0(DC_LAMBDA))) dim = (int) (PIHALF / ATAN(1.0 / P0(DC_LAMBDA)));
  set_maxdim(OWN, 0, dim);
  RETURN_NOERROR;
}
void rangedampedcosine(model *cov, range_type *range){
  int dim = OWNLOGDIM(0);
  if (dim <= 3)
    range->min[DC_LAMBDA] = dim==1 ? 0 : dim==2 ? 1 : 1.7320508075688771932;
  else range->min[DC_LAMBDA] = 1.0 / TAN(PIHALF / dim); 
  range->max[DC_LAMBDA] = RF_INF;
  range->pmin[DC_LAMBDA] = range->min[DC_LAMBDA] +  1e-10;
  range->pmax[DC_LAMBDA] = 10;
  range->openmin[DC_LAMBDA] = false;
  range->openmax[DC_LAMBDA] = true;
}


/* De Wijsian */
#define DEW_ALPHA 0 // for both dewijsian models
#define DEW_D 1
void dewijsian(double *x, model *cov, double *v){
  double alpha = P0(DEW_ALPHA);
  *v = -LOG(1.0 + POW(*x, alpha));
}
void Ddewijsian(double *x, model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
      p  = POW(*x, alpha - 1.0) ;
  *v = - alpha * p / (1.0 + *x * p);
}
void DDdewijsian(double *x, model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
      p = POW(*x, alpha - 2.0),
      ha = p * *x * *x,
      haP1 = ha + 1.0;
  *v = alpha * p * (1.0 - alpha + ha) / (haP1 * haP1);
}

void D3dewijsian(double *x, model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
    p = POW(*x, alpha - 3.0),
    ha = p * *x * *x * *x,
    haP1 = ha + 1.0,
    haP12 = haP1*haP1;
  *v = alpha * p * (alpha*alpha*(ha - 1) +3*alpha*haP1 -2*haP12  )
      / (haP1 * haP12);
}

void D4dewijsian(double *x, model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
    alpha2 = alpha*alpha,
    p = POW(*x, alpha - 4.0),
    ha = p * *x * *x * *x* *x,
    haSq = ha * ha,
    haP1 = ha + 1.0,
    haP12 = haP1*haP1;
  *v = -alpha * p * (alpha2*alpha*(haSq - 4*ha + 1) +
                     6*alpha2*(haSq - 1) + 11*alpha*haP12 -
                     6*haP1*haP12  ) / (haP12*haP12);
}

void Inversedewijsian(double *x, model *cov, double *v){ 
  double alpha = P0(DEW_ALPHA);
  *v = POW(EXP(*x) - 1.0, 1.0 / alpha);    
} 
int checkdewijsian(model *cov){
  double alpha = P0(DEW_ALPHA);
  cov->logspeed = alpha;
  RETURN_NOERROR;
}

void rangedewijsian(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[DEW_ALPHA] = 0.0;
  range->max[DEW_ALPHA] = 2.0;
  range->pmin[DEW_ALPHA] = UNIT_EPSILON;
  range->pmax[DEW_ALPHA] = 2.0;
  range->openmin[DEW_ALPHA] = true;
  range->openmax[DEW_ALPHA] = false; 
}

void coinitdewijsian(model *cov, localinfotype *li) {
  double
      thres[3] = {0.5, 1.0, 1.8},
      alpha=P0(DEW_ALPHA);

  if (alpha <= thres[0]) {
      li->instances = 2;
      li->value[0] = 0.5;
      li->value[1] = 1.0;
      li->msg[0] = li->msg[1] = MSGLOCAL_OK;
    } else {
      if (alpha <= thres[1]) {
          li->instances = 1;
          li->value[0] = 1.0; //  q[CUTOFF_A]
          li->msg[0] =MSGLOCAL_OK;
        } else {
          if (alpha <= thres[2]) {
              li->instances = 1;
              li->value[0] = CUTOFF_THIRD_CONDITION ; //  q[CUTOFF_A]
              li->msg[0] = MSGLOCAL_OK;
            }
          else {
	    //TO DO: this is copied from Cauchy model and must be understood and changed
              li->instances = 1;
              li->value[0] = CUTOFF_THIRD_CONDITION ; //  q[CUTOFF_A]
              li->msg[0] = MSGLOCAL_JUSTTRY;
            }
        }
    }
}

/* De Wijsian B */
#define DEW_RANGE 1
void DeWijsian(double *x, model *cov, double *v){
  double alpha = P0(DEW_ALPHA),
    range = P0(DEW_RANGE);
  *v = 0.0;
  if (*x < range) *v = 1.0-LOG(1.0 + POW(*x, alpha)) / LOG(1.0 + POW(range, alpha));
}

void InverseDeWijsian(double *x, model *cov, double *v){ 
  double alpha = P0(DEW_ALPHA),
    range = P0(DEW_RANGE);
  *v = 0.0;
  if (*x < 1.0)
    *v = POW(POW(1.0 + POW(range, alpha),  1.0 - *x) - 1.0, 1.0 / alpha);
} 

void rangeDeWijsian(model VARIABLE_IS_NOT_USED *cov, range_type *range){
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



// Brownian motion 

int initlsfbm(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {   
  int d = OWNLOGDIM(0);
  double alpha = P0(LOCALLY_BROWN_ALPHA),
    a2 = 0.5 * alpha,
    d2 = 0.5 * d;
  if (PisNULL(LOCALLY_BROWN_C)) {
    QVALUE = EXP(-alpha * LOG2 + lgammafn(a2 + d2) + lgammafn(1.0 - a2) 
		 - lgammafn(d2));
    if (PL >= PL_DETAILSUSER) {
      PRINTF("'%.50s' in '%.50s' equals %10g for '%.50s'=%10g\n", KNAME(LOCALLY_BROWN_C),
	     NICK(cov), QVALUE, KNAME(LOCALLY_BROWN_ALPHA), alpha);
    }
  } else {
    QVALUE =  P0(LOCALLY_BROWN_C);
  }
  cov->taylor[0][TaylorPow] = cov->tail[0][TaylorPow] = alpha;
  RETURN_NOERROR;
}

int checklsfbm(model *cov){
  int err;
  if (P(BROWN_ALPHA) == NULL) ERR("alpha must be given");
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  double alpha = P0(BROWN_ALPHA);
  cov->logspeed = RF_INF;
  cov->full_derivs = alpha <= 1.0 ? 0 : alpha < 2.0 ? 1 : cov->rese_derivs;
  if (cov->q == NULL) {
    EXTRA_Q;
    if ((err = initlsfbm(cov, NULL)) != NOERROR) RETURN_ERR(err);
  }
  RETURN_NOERROR;
}
void rangelsfbm(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[LOCALLY_BROWN_ALPHA] = 0.0;
  range->max[LOCALLY_BROWN_ALPHA] = 2.0;
  range->pmin[LOCALLY_BROWN_ALPHA] = UNIT_EPSILON;
  range->pmax[LOCALLY_BROWN_ALPHA] = 2.0;
  range->openmin[LOCALLY_BROWN_ALPHA] = true;
  range->openmax[LOCALLY_BROWN_ALPHA] = false;

  range->min[LOCALLY_BROWN_C] = 0.5;
  range->max[LOCALLY_BROWN_C] = RF_INF;
  range->pmin[LOCALLY_BROWN_C] = 1.0;
  range->pmax[LOCALLY_BROWN_C] = 1000;
  range->openmin[LOCALLY_BROWN_C] = true;
  range->openmax[LOCALLY_BROWN_C] = P(LOCALLY_BROWN_ALPHA)==NULL;
}

void lsfbm(double *x, model *cov, double *v) {
  if (*x > 1.0)
    ERR1("the coordinate distance in '%.50s' must be at most 1.", NICK(cov));
  double alpha = P0(LOCALLY_BROWN_ALPHA);
  *v = QVALUE - POW(*x, alpha);
}
/* lsfbm: first derivative at t=1 */
void Dlsfbm(double *x, model *cov, double *v) 
{// FALSE VALUE FOR *x==0 and  alpha < 1
  if (*x > 1.0)
    ERR1("the coordinate distance in '%.50s' must be at most 1.", NICK(cov));
  double alpha = P0(LOCALLY_BROWN_ALPHA);
  *v = (*x != 0.0) ? -alpha * POW(*x, alpha - 1.0)
    : alpha > 1.0 ? 0.0 
    : alpha < 1.0 ? RF_NEGINF
	      : -1.0;
}
/* lsfbm: second derivative at t=1 */
void DDlsfbm(double *x, model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  if (*x > 1.0)
    ERR1("the coordinate distance in '%.50s' must be at most 1.", NICK(cov));
  double alpha = P0(LOCALLY_BROWN_ALPHA);
  *v = (alpha == 1.0) ? 0.0
    : (*x != 0.0) ? -alpha * (alpha - 1.0) * POW(*x, alpha - 2.0)
    : alpha < 1.0 ? RF_INF 
	      : alpha < 2.0 ? RF_NEGINF 
			: -2.0;
}
void D3lsfbm(double *x, model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  if (*x > 1.0)
    ERR1("the coordinate distance in '%.50s' must be at most 1.", NICK(cov));
  double alpha = P0(LOCALLY_BROWN_ALPHA);
  *v = alpha == 1.0 || alpha == 2.0 ? 0.0 
   : (*x != 0.0) ? -alpha * (alpha - 1.0) * (alpha - 2.0) * POW(*x, alpha-3.0)
	      : alpha < 1.0 ? RF_NEGINF 
			: RF_INF;
}

void D4lsfbm(double *x, model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  if (*x > 1.0)
    ERR1("the coordinate distance in '%.50s' must be at most 1.", NICK(cov));
  double alpha = P0(LOCALLY_BROWN_ALPHA);
  *v = alpha == 1.0 || alpha == 2.0 ? 0.0 
    : (*x != 0.0) 
    ? -alpha * (alpha - 1.0) * (alpha - 2.0) * (alpha - 3.0) * POW(*x,alpha-4.0)
    : alpha < 1.0 ? RF_INF 
    : RF_NEGINF;
}

void Inverselsfbm(double *x, model *cov, double *v) {
  if (*x > 1.0)
    ERR1("the coordinate distance in '%.50s' must be at most 1.", NICK(cov));
  double alpha = P0(LOCALLY_BROWN_ALPHA);
  *v = POW(QVALUE - *x, 1.0 / alpha);
}

// Brownian motion 
void fractalBrownian(double *x, model *cov, double *v) {
  double alpha = P0(BROWN_ALPHA);
  *v = - POW(*x, alpha);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}


void logfractalBrownian(double *x, model *cov, double *v, double *Sign) {
  double alpha = P0(BROWN_ALPHA);
  *v = LOG(*x) * alpha;//this is an invalid covariance function!
  *Sign = - 1.0;
  // keep definition such that the value at the origin is 0
}

/* fractalBrownian: first derivative at t=1 */
void DfractalBrownian(double *x, model *cov, double *v) 
{// FALSE VALUE FOR *x==0 and  alpha < 1
  double alpha = P0(BROWN_ALPHA);
  *v = (*x != 0.0) ? -alpha * POW(*x, alpha - 1.0)
    : alpha > 1.0 ? 0.0 
    : alpha < 1.0 ? RF_NEGINF
    : -1.0;
}
/* fractalBrownian: second derivative at t=1 */
void DDfractalBrownian(double *x, model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  double alpha = P0(BROWN_ALPHA);
  *v = (alpha == 1.0) ? 0.0
    : (*x != 0.0) ? -alpha * (alpha - 1.0) * POW(*x, alpha - 2.0)
    : alpha < 1.0 ? RF_INF 
    : alpha < 2.0 ? RF_NEGINF 
    : -2.0;
}
void D3fractalBrownian(double *x, model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  double alpha = P0(BROWN_ALPHA);
  *v = alpha == 1.0 || alpha == 2.0 ? 0.0 
    : (*x != 0.0) ? -alpha * (alpha - 1.0) * (alpha - 2.0) * POW(*x, alpha-3.0)
    : alpha < 1.0 ? RF_NEGINF 
    : RF_INF;
}

void D4fractalBrownian(double *x, model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  alpha < 2
  double alpha = P0(BROWN_ALPHA);
  *v = alpha == 1.0 || alpha == 2.0 ? 0.0 
    : (*x != 0.0) ? -alpha * (alpha - 1.0) * (alpha - 2.0) * (alpha - 3.0) *
                     POW(*x, alpha-4.0)
    : alpha < 1.0 ? RF_INF 
    : RF_NEGINF;
}

int initfractalBrownian(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double alpha = P0(BROWN_ALPHA);
  cov->taylor[0][TaylorPow] = cov->tail[0][TaylorPow]=alpha;
  //PMI(cov);
  RETURN_NOERROR;
}

int checkfractalBrownian(model *cov){
  int err;
  double alpha = P0(BROWN_ALPHA);
  cov->logspeed = RF_INF;
  cov->full_derivs = alpha <= 1.0 ? 0 : alpha < 2.0 ? 1 : cov->rese_derivs;
  if ((err = initfractalBrownian(cov, NULL)) != NOERROR) RETURN_ERR(err);
  RETURN_NOERROR;
}


void InversefractalBrownian(double *x, model *cov, double *v) {
  double alpha = P0(BROWN_ALPHA);
  *v = POW(*x, 1.0 / alpha);
}

void rangefractalBrownian(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[BROWN_ALPHA] = 0.0;
  range->max[BROWN_ALPHA] = 2.0;
  range->pmin[BROWN_ALPHA] = UNIT_EPSILON;
  range->pmax[BROWN_ALPHA] = 2.0;
  range->openmin[BROWN_ALPHA] = true;
  range->openmax[BROWN_ALPHA] = false;
}
void ieinitBrownian(model *cov, localinfotype *li) {
  li->instances = 1;
  li->value[0] = (OWNLOGDIM(0) <= 2)
    ? ((P0(BROWN_ALPHA) <= 1.5) ? 1.0 : 2.0)
    : ((P0(BROWN_ALPHA) <= 1.0) ? 1.0 : 2.0);
  li->msg[0] = OWNLOGDIM(0) <= 3 ? MSGLOCAL_OK : MSGLOCAL_JUSTTRY;
}



/* FD model */
#define FD_ALPHA 0
void FD(double *x, model *cov, double *v){
  double
    alpha = P0(FD_ALPHA),
    d = alpha * 0.5,
    y = *x;
  if (y == RF_INF) {*v = 0.0; return;}
  double
    k = TRUNC(y);
#ifdef DO_PARALLEL
  double sk = 1.0,
    kold = 0.0;
#else
  static double kold = RF_INF, // only if not parallel
    sk = RF_INF,
    dold=RF_INF;
  if (dold!=d || kold > k) {
    sk = 1;
    kold = 0.0;
  }
#endif	    
  // Sign (-1)^k is (kold+d), 16.11.03, checked. 
  for (; kold<k; kold += 1.0) sk =  sk * (kold + d) / (kold + 1.0 - d);
#ifndef DO_PARALLEL
  dold = d;
  kold = k;
#endif	    
  if (k == y) {
    *v = sk;
  } else {
    double skP1 = sk * (kold + d) / (kold + 1.0 - d);
    *v = sk + (y - k) * (skP1 - sk);
  }
}

void rangeFD(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[FD_ALPHA] = - 1.0;
  range->max[FD_ALPHA] = 1.0;
  range->pmin[FD_ALPHA] = range->min[FD_ALPHA] + UNIT_EPSILON;
  range->pmax[FD_ALPHA] = range->max[FD_ALPHA] - UNIT_EPSILON;
  range->openmin[FD_ALPHA] = false;
  range->openmax[FD_ALPHA] = true;
}



/* fractgauss */
#define FG_ALPHA 0
void fractGauss(double *x, model *cov, double *v){
  double y = *x, alpha = P0(FG_ALPHA);
  *v = (y == 0.0) ? 1.0 :  (y==RF_INF) ? 0.0 : 
    0.5 *(POW(y + 1.0, alpha) - 2.0 * POW(y, alpha) + POW(FABS(y - 1.0),alpha));
}
void rangefractGauss(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[FG_ALPHA] = 0.0;
  range->max[FG_ALPHA] = 2.0;
  range->pmin[FG_ALPHA] = UNIT_EPSILON;
  range->pmax[FG_ALPHA] = 2.0;
  range->openmin[FG_ALPHA] = true;
  range->openmax[FG_ALPHA] = false;
}



/* generalised fractal Brownian motion */
void genBrownian(double *x, model *cov, double *v) {
  double 
    alpha = P0(BROWN_ALPHA),
    beta =  P0(BROWN_GEN_BETA),
    delta = beta / alpha;
  *v = 1.0 - POW(POW(*x, alpha) + 1.0, delta);
  //this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}

void loggenBrownian(double *x, model *cov, double *v, double *Sign) {
  genBrownian(x, cov, v);
  *v = LOG(-*v);
  *Sign = - 1.0;
}
void InversegenBrownian(double *x, model *cov, double *v) {
  double 
    alpha = P0(BROWN_ALPHA),
    beta =  P0(BROWN_GEN_BETA),
    delta = beta / alpha;
  *v = POW(POW(*x + 1.0, 1.0/delta) - 1.0, 1.0/alpha); 
}

int checkgenBrownian(model *cov){
  cov->logspeed = RF_INF;
  RETURN_NOERROR;
}
void rangegenBrownian(model VARIABLE_IS_NOT_USED *cov, range_type *range){
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



/* gengneiting */
void genGneiting(double *x, model *cov, double *v)
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
  *v *=  POW(1.0 - y, s);
}


// control thanks to http://calc101.com/webMathematica/derivatives.jsp#topdoit

void DgenGneiting(double *x, model *cov, double *v)
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
  *v *=  -POW(1.0 - y, s - 1.0);
  
}
void DDgenGneiting(double *x, model *cov, double *v){
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
  *v *=  POW(1.0 - y, s - 2.0);
}


int checkgenGneiting(model *cov){
  double mu=P0(GENGNEITING_MU), 
    dim = 2.0 * mu;
  set_maxdim(OWN, 0, ISNAN(dim) || dim >= INFDIM ? INFDIM : (int) dim);
  RETURN_NOERROR;
}

void rangegenGneiting(model *cov, range_type *range){
  range->min[GENGNEITING_K] = range->pmin[GENGNEITING_K] = 0;
  range->max[GENGNEITING_K] = range->pmax[GENGNEITING_K] = 3; 
  range->openmin[GENGNEITING_K] = false;
  range->openmax[GENGNEITING_K] = false;
  
  range->min[GENGNEITING_MU] = 0.5 * (double) OWNLOGDIM(0); 
  range->max[GENGNEITING_MU] =  RF_INF;
  range->pmin[GENGNEITING_MU] = range->min[GENGNEITING_MU];
  range->pmax[GENGNEITING_MU] = range->pmin[GENGNEITING_MU] + 10.0;
  range->openmin[GENGNEITING_MU] = false;
  range->openmax[GENGNEITING_MU] = true;
}

/*
double porcusG(double t, double nu, double mu, double gamma) {
  double t0 = FABS(t);
  if (t0 > mu) return 0.0;
  return POW(t0, nu) * POW(1.0 - t0 / mu, gamma);
}
*/


double biGneitQuot(double t, double* scale, double *gamma) {
  double t0 = FABS(t);
  return POW(1.0 - t0 / scale[0], gamma[0]) *
    POW(1.0 - t0 / scale[3], gamma[3]) * POW(1.0 - t0 / scale[1], -2 * gamma[1]);
}

void biGneitingbasic(model *cov, 
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
  if (sum > 2.0 * gamma[1]) ERR("values of gamma not valid.");

  min = 1.0;
  a = 2 * gamma[1] - gamma[0] - gamma[3];
  b = - 2 * gamma[1] * (scale[0] + scale[3]) + gamma[0] * (scale[1] + scale[3])
    + gamma[3] * (scale[1] + scale[0]);
  c = 2 * gamma[1] * scale[0] * scale[3] - gamma[0] * scale[1] * scale[3]
    - gamma[3] * scale[1] * scale[0]; 
  det = b * b - 4 * a * c;
  
  if (det >= 0) {
    det = SQRT(det);
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
  cc[1] = cc[2] = C[1] = rho * SQRT(C[0] * C[2] * min) *
    POW(scale[1] * scale[1] / (scale[0] * scale[3]), 0.5 * (nu + 1 + 2.0 * k)) *
    EXP(lgammafn(1.0 + gamma[1]) - lgammafn(2.0 + nu + gamma[1] + kP1)
	+ 0.5 * (  lgammafn(2 + nu + gamma[0] + kP1) - lgammafn(1 + gamma[0])
		 + lgammafn(2 + nu + gamma[3] + kP1) - lgammafn(1 + gamma[3]))
	);
}


int initbiGneiting(model *cov, gen_storage *stor) {
  double  
    *rho = P(GNEITING_RHORED),
    *scale = P(GNEITING_S),
    *sred = P(GNEITING_SRED),
    *p_gamma = P(GNEITING_GAMMA),
    *c = P(GNEITING_C),
    *cdiag = P(GNEITING_CDIAG);
  bool check = stor->check;
  biwm_storage *S = cov->Sbiwm;
  assert(S != NULL);
  
  if (scale == NULL) {
    PALLOC(GNEITING_S, 2, 1);
    scale = P(GNEITING_S);
    scale[0] = scale[1] = 1.0;
  }
  
  if (sred == NULL) {
    PALLOC(GNEITING_SRED, 1, 1);
    sred = P(GNEITING_SRED);
    sred[0] = 1.0;
  }
  
  if (sred[0] == 1.0) {
    double sum =0.0;
    if (scale[0] <= scale[1]) sum += p_gamma[0];
    if (scale[1] <= scale[0]) sum += p_gamma[2];
    if (sum > 2.0 * p_gamma[1]) {
      SERR7("if '%.50s'=1, then %.50s[2] must be greater than min(%.50s[1], %.50s[3]) / 2 if%.50s[1] and %.50s[3] differ and %.50s[1] otherwise.", KNAME(GNEITING_SRED),
	    KNAME(GNEITING_GAMMA), KNAME(GNEITING_GAMMA),
	    KNAME(GNEITING_GAMMA), KNAME(GNEITING_GAMMA), 
	    KNAME(GNEITING_GAMMA), KNAME(GNEITING_GAMMA));
    }
  }

  if  (S->cdiag_given) { // cdiag or rho given
    if (cdiag == NULL) {
      PALLOC(GNEITING_CDIAG, 2, 1);
      cdiag = P(GNEITING_CDIAG);
      cdiag[0] = cdiag[1] = 1.0;
    }
    if (rho == NULL) 
      QERRC2(GNEITING_RHORED, 
	     "'%.50s' and '%.50s' must be set at the same time ", GNEITING_CDIAG,
	     GNEITING_RHORED);
    if (check && c != NULL) {
      if (cov->nrow[GNEITING_C] != 3 || cov->ncol[GNEITING_C] != 1)
	QERRC(GNEITING_C, "must be a 3 x 1 vector");
      if (FABS(c[i11] - cdiag[0]) > c[i11] * epsilon || 
	  FABS(c[i22] - cdiag[1]) > c[i22] * epsilon ) {
	QERRC1(GNEITING_C, "does not match '%.50s'.", GNEITING_CDIAG);
      }
      double savec12 = c[i21];
      biGneitingbasic(cov, S->scale, S->gamma, S->c);
      if (FABS(c[i21] - savec12) > FABS(c[i21]) * epsilon)
 	QERRC1(GNEITING_C, "does not match '%.50s'.", GNEITING_RHORED);
    } else {
      if (PisNULL(GNEITING_C)) PALLOC(GNEITING_C, 3, 1);
      c = P(GNEITING_C);
      biGneitingbasic(cov, S->scale, S->gamma, S->c);
    }
  } else {
    if (c == NULL) 
      QERRC2(GNEITING_C, "either '%.50s' or '%.50s' must be set", 
	    GNEITING_C, GNEITING_RHORED);
    if (!ISNAN(c[i11]) && !ISNAN(c[i22]) && (c[i11] < 0.0 || c[i22] < 0.0))
      QERRC2(GNEITING_C, "variance parameter %.50s[0], %.50s[2] must be non-negative",
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
  RETURN_NOERROR;
}


void kappa_biGneiting(int i, model *cov, int *nr, int *nc){
  *nc = *nr = i < DefList[COVNR].kappas? 1 : -1;
  if (i==GNEITING_S || i==GNEITING_CDIAG) *nr=2; else
    if (i==GNEITING_GAMMA || i==GNEITING_C) *nr=3 ;
}

void biGneiting(double *x, model *cov, double *v) { 
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
    z = FABS(*x / S->scale[i]);
    P(GENGNEITING_MU)[0] = mu + S->gamma[i] + 1.0;
    genGneiting(&z, cov, v + i);
    v[i] *= S->c[i]; 
  }
  P(GENGNEITING_MU)[0] = mu;
}


void DbiGneiting(double *x, model *cov, double *v){ 
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
    z = FABS(*x / S->scale[i]);
    P(GENGNEITING_MU)[0] = mu + S->gamma[i] + 1.0;
    DgenGneiting(&z, cov, v + i);
    v[i] *= S->c[i] / S->scale[i];
  }
  P(GENGNEITING_MU)[0] = mu;
}


void DDbiGneiting(double *x, model *cov, double *v){ 
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
    z = FABS(*x / S->scale[i]);
    P(GENGNEITING_MU)[0] = mu + S->gamma[i] + 1.0;
    DDgenGneiting(&z, cov, v + i);
    v[i] *= S->c[i] / (S->scale[i] * S->scale[i]);
  }
  P(GENGNEITING_MU)[0] = mu;
}

int checkbiGneiting(model *cov) { 
  int err;
  gen_storage s;
  gen_NULL(&s);
  s.check = true;
  
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  if (PisNULL(GNEITING_K)) QERRC(GNEITING_K, "must be given.");
  if (PisNULL(GNEITING_MU)) QERRC(GNEITING_MU, "must be given.");
  if (PisNULL(GNEITING_GAMMA)) QERRC(GNEITING_GAMMA,"must be given.");

  if (cov->Sbiwm == NULL) {
    NEW_STORAGE(biwm);
    biwm_storage *S = cov->Sbiwm;
    S->cdiag_given = !PisNULL(GNEITING_CDIAG) || !PisNULL(GNEITING_RHORED);
  }
 
  if ((err=initbiGneiting(cov, &s)) != NOERROR) RETURN_ERR(err);

  int dim = 2.0 * P0(GENGNEITING_MU);
  set_maxdim(OWN, 0, ISNAN(dim) || dim >= INFDIM ? INFDIM : (int) dim);
  RETURN_NOERROR;
}
  
sortsofparam sortof_biGneiting(model *cov, int k, int VARIABLE_IS_NOT_USED row,
			       int VARIABLE_IS_NOT_USED col,
			       sort_origin origin){
  biwm_storage *S = cov->Sbiwm;
  if (S == NULL) return UNKNOWNPARAM;
  switch(k) {
  case GNEITING_K :     return ONLYRETURN;
  case GNEITING_MU :    return CRITICALPARAM;
  case GNEITING_S :     return SCALEPARAM;
  case GNEITING_SRED :  return ANYPARAM;
  case GNEITING_GAMMA : return ANYPARAM;
  case GNEITING_CDIAG :
    //    printf("xxxxxx %d\n", S ==NULL);
    //    printf("%d\n",  S->cdiag_given);
    //printf("%d\n", origin != original);
    return S->cdiag_given || origin != original ? VARPARAM : IGNOREPARAM;
  case GNEITING_RHORED :
    return S->cdiag_given || origin != original ? ANYPARAM : IGNOREPARAM;
  case GNEITING_C:
    return S->cdiag_given || origin == mle_conform ? IGNOREPARAM : ONLYRETURN;
  default : BUG;
  }
}
 
sortsofparam sortof_biGneiting_INisOUT(model *cov, int k, int VARIABLE_IS_NOT_USED row,
			  int VARIABLE_IS_NOT_USED col){
  biwm_storage *S = cov->Sbiwm;
  if (S == NULL) return UNKNOWNPARAM;
  switch(k) {
  case GNEITING_K :     return ONLYRETURN;
  case GNEITING_MU :    return CRITICALPARAM;
  case GNEITING_S :     return SCALEPARAM;
  case GNEITING_SRED :  return ANYPARAM;
  case GNEITING_GAMMA : return ANYPARAM;
  case GNEITING_CDIAG : return S->cdiag_given ? VARPARAM : VARONLYMLE;
  case GNEITING_RHORED :return S->cdiag_given ? ANYPARAM : ONLYMLE;
  case GNEITING_C:      return S->cdiag_given ? IGNOREPARAM : ONLYRETURN;
  default : BUG;
  }
}

void rangebiGneiting(model *cov, range_type *range){
 // *n = P0(GNEITING_K], 
  range->min[GNEITING_K] = range->pmin[GNEITING_K] = 0;
  range->max[GNEITING_K] = range->pmax[GNEITING_K] = 3;
  range->openmin[GNEITING_K] = false;
  range->openmax[GNEITING_K] = false;
  
 // *mu = P0(GNEITING_MU], 
  range->min[GNEITING_MU] = 0.5 * (double) OWNLOGDIM(0); 
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
  range->min[GNEITING_RHORED] = - 1.0;
  range->max[GNEITING_RHORED] = 1;
  range->pmin[GNEITING_RHORED] = -0.95;
  range->pmax[GNEITING_RHORED] = 0.95;
  range->openmin[GNEITING_RHORED] = false;
  range->openmax[GNEITING_RHORED] = false;    

 //    *rho = P0(GNEITING_C]; 
   range->min[GNEITING_C] = RF_NEGINF;
  range->max[GNEITING_C] = RF_INF;
  range->pmin[GNEITING_C] = - 1000;
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
void Gneiting(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){ 
  double y = *x * cov->q[0];
  genGneiting(&y, cov, v);
}

 
void DGneiting(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){ 
  double y = *x * cov->q[0];
  DgenGneiting(&y, cov, v);
  *v  *= cov->q[0];
}


void DDGneiting(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){ 
  double y = *x * cov->q[0];
  DgenGneiting(&y, cov, v);
  *v  *= cov->q[0] * cov->q[0];
}

int checkGneiting(model *cov) { 
  int err;
  kdefault(cov, GNEITING_ORIG, 1);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  bool orig = P0INT(GNEITING_ORIG);
  PFREE(GNEITING_ORIG);

  assert(cov->q == NULL);
  set_nr(OWN, GNEITING_INTERN);
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


void rangeGneiting(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  booleanRange(GNEITING_ORIG);
}



/* iaco cesare model */
#define CES_NU 0
#define CES_LAMBDA 1
#define CES_DELTA 2
void IacoCesare(double *x, model *cov, double *v){
    double
      nu = P0(CES_NU), 
      lambda=P0(CES_LAMBDA), 
      delta=P0(CES_DELTA);
     *v = POW(1.0 + POW(x[0], nu) + POW(x[1], lambda), - delta); 
     //   printf("x=%10g %10g  nu=%10g %10g %10g 1-v=%10g\n", x[0], x[1], nu, lambda, delta, 1-*v);
   //   APMI(cov);
}
void rangeIacoCesare(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[CES_NU] = 0.0;
  range->max[CES_NU] = 2.0;
  range->pmin[CES_NU] = 1e-3;
  range->pmax[CES_NU] = 2.0;
  range->openmin[CES_NU] = true;
  range->openmax[CES_NU] = false;

  range->min[CES_LAMBDA] = 0.0;
  range->max[CES_LAMBDA] = 2.0;
  range->pmin[CES_LAMBDA] = 1e-3;
  range->pmax[CES_LAMBDA] = 2.0;
  range->openmin[CES_LAMBDA] = true;
  range->openmax[CES_LAMBDA] = false;
 
  range->min[CES_DELTA] = 0.0;
  range->max[CES_DELTA] = RF_INF;
  range->pmin[CES_DELTA] = 1e-3;
  range->pmax[CES_DELTA] = 10.0;
  range->openmin[CES_DELTA] = true;
  range->openmax[CES_DELTA] = true;
}



/* Kolmogorov model */
void Kolmogorov(double *x, model *cov, double *v){
#define fourthird 1.33333333333333333333333333333333333
#define onethird 0.333333333333333333333333333333333
  int d, i, j,
    dim = OWNLOGDIM(0),
    dimP1 = dim + 1,
    dimsq = dim * dim;
  double
    rM23, r23,
    r2 =0.0;

  assert(dim == VDIM0 && dim == VDIM1);

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

  r23 = -POW(r2, onethird);  // ! -
  for (d=0; d<dimsq; v[d++] *= r23);

}


int checkKolmogorov(model *cov) { 
    if (ANYDIM != 3) SERR1("dim (%d) != 3", ANYDIM);
   RETURN_NOERROR;
}
 

/* local-global distinguisher */
#define LGD_ALPHA 0
#define LGD_BETA 1
void lgd1(double *x, model *cov, double *v) {
  double y = *x, alpha = P0(LGD_ALPHA), beta=P0(LGD_BETA);
  *v = 1.0;
  if (y != 0.0) {
    if (y < 1.0) *v = 1.0 - beta / (alpha + beta) * POW(y, alpha);
    else *v = alpha / (alpha + beta) * POW(y,  -beta);
  }
}
void Inverselgd1(double *x, model *cov, double *v) {
  double alpha = P0(LGD_ALPHA), beta=P0(LGD_BETA);
  ERR("scle of lgd1 not programmed yet"); 
   // 19 next line?!
  // 1.0 / .... fehlt auch
  if (19*alpha < beta) *v = EXP(LOG(1.0 - *x * (alpha + beta) / beta) / alpha);
  else *v = EXP(LOG(*x * (alpha + beta) / alpha) / beta);
}
void Dlgd1(double *x, model *cov, double *v){
  double y=*x, pp, alpha = P0(LGD_ALPHA), beta=P0(LGD_BETA);
  if (y == 0.0) {
    *v = 0.0;// falscher Wert, aber sonst NAN-Fehler#
    return;
  }
  pp = ( (y < 1.0) ? alpha : -beta ) - 1.0;
  *v = - alpha * beta / (alpha + beta) * EXP(pp * y);
}
int checklgd1(model *cov){
  double dim = 2.0 * (1.5 - P0(LGD_ALPHA));
  set_maxdim(OWN, 0, ISNAN(dim) || dim >= 2.0 ? 2 : (int) dim);
  RETURN_NOERROR;
}
void rangelgd1(model *cov, range_type *range) {
  range->min[LGD_ALPHA] = 0.0;
  range->max[LGD_ALPHA] = (OWNLOGDIM(0)==2) ? 0.5 : 1.0;
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



#define MULTIQUAD_DELTA 0
#define MULTIQUAD_TAU 1
#define MULTIQUAD_EPS 1e-7
void multiquad(double *x, model *cov, double *v){
  
  double delta = P0(MULTIQUAD_DELTA), // Auslesen der Parameter aus cov
    deltaM1 = delta - 1.0,
    tau = P0(MULTIQUAD_TAU); 
  double y = *x;
  if (y < 0.0 || y >= PI) BUG;
  
  *v = POW( deltaM1 * deltaM1 / (1.0 + delta * (delta - 2.0 * COS(y))), tau);
  return;
}


/* range der Parameter delta & tau */
void rangemultiquad(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[MULTIQUAD_DELTA] = 0.0;  // true range
  range->max[MULTIQUAD_DELTA] = 1.0;   
  range->pmin[MULTIQUAD_DELTA] = MULTIQUAD_EPS; // practical range (assumed)
  range->pmax[MULTIQUAD_DELTA] = 1.0 - MULTIQUAD_EPS; 
  range->openmin[MULTIQUAD_DELTA] = true;  // lower endpoint included 
  range->openmax[MULTIQUAD_DELTA] = true; // upper endpoint excluded 

  range->min[MULTIQUAD_TAU] = 0;  // true range
  range->max[MULTIQUAD_TAU] = RF_INF;   
  range->pmin[MULTIQUAD_TAU] = 0.0001; // practical range (assumed)
  range->pmax[MULTIQUAD_TAU] = 500; 
  range->openmin[MULTIQUAD_TAU] = true;  // lower endpoint included 
  range->openmax[MULTIQUAD_TAU] = true; // upper endpoint excluded  
}



// Brownian motion 
#define OESTING_BETA 0
void oesting(double *x, model *cov, double *v) {
  // klein beta interagiert mit 'define beta Rf_beta' in Rmath
  double Beta = P0(OESTING_BETA),
    x2 = *x * *x;
  *v = - x2 * POW(1 + x2, -Beta);//this is an invalid covariance function!
  // keep definition such that the value at the origin is 0
}
/* oesting: first derivative at t=1 */
void Doesting(double *x, model *cov, double *v) 
{// FALSE VALUE FOR *x==0 and  Beta < 1
  double Beta = P0(OESTING_BETA),
    x2 = *x * *x;
  *v = 2 * *x  * (1 + (1-Beta) * x2) * POW(1 + x2, -Beta-1);
}
/* oesting: second derivative at t=1 */
void DDoesting(double *x, model *cov, double *v)  
{// FALSE VALUE FOR *x==0 and  beta < 2
  double Beta = P0(OESTING_BETA),
    x2 = *x * *x;
  *v = 2 * (1 + (2 - 5 * Beta) * x2 + (1 - 3 * Beta + 2 * Beta * Beta) * 
	    x2  *x2) * POW(1 + x2, -Beta-2);
}
  
int initoesting(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double Beta = P0(OESTING_BETA);
  cov->taylor[1][TaylorConst] = + Beta;
  cov->taylor[2][TaylorConst] = - 0.5 * Beta * (Beta + 1);  
  cov->tail[0][TaylorPow] = 2 - 2 * Beta;
  RETURN_NOERROR;
}

int checkoesting(model *cov){
  int err;
  cov->logspeed = RF_INF;
  cov->full_derivs = cov->rese_derivs;
  if ((err = initoesting(cov, NULL)) != NOERROR) RETURN_ERR(err);
 RETURN_NOERROR;
}
void rangeoesting(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[OESTING_BETA] = 0.0;
  range->max[OESTING_BETA] = 1.0;
  range->pmin[OESTING_BETA] = UNIT_EPSILON;
  range->pmax[OESTING_BETA] = 1.0;
  range->openmin[OESTING_BETA] = true;
  range->openmax[OESTING_BETA] = false;
}


/* penta */
void penta(double *x, model VARIABLE_IS_NOT_USED *cov, double *v)
{ ///
  double y=*x, y2 = y * y;
  *v = 0.0;
  if (y < 1.0) 
    *v = (1.0 + y2 * (-7.333333333333333 
		      + y2 * (33.0 +
			      y * (-38.5 + 
				   y2 * (16.5 + 
					 y2 * (-5.5 + 
					       y2 * 0.833333333333333))))));
}
void Dpenta(double *x, model VARIABLE_IS_NOT_USED *cov, double *v)
{ 
  double y=*x, y2 = y * y;
  *v = 0.0;
  if (y < 1.0) 
    *v = y * (-14.66666666666666667 + 
	      y2 * (132.0 + 
		    y * (-192.5 + 
			 y2 * (115.5 + 
			       y2 * (-49.5 + 
				     y2 * 9.16666666666666667)))));
  
}
void Inversepenta(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = (*x==0.05) ? 1.0 / 1.6552838957365 : RF_NA;
}


/* power model */ 
#define POW_ALPHA 0
void power(double *x, model *cov, double *v){
  double alpha = P0(POW_ALPHA), y = *x;
  *v = 0.0;
  if (y < 1.0) *v = POW(1.0 - y, alpha);
}
void TBM2power(double *x, model *cov, double *v){
  // only alpha=2 up to now !
  double y = *x;
  if (P0(POW_ALPHA) != 2.0) 
    ERR("TBM2 of power only allowed for alpha=2");
  if (y > 1.0) *v = (1.0 - 2.0 * y *(ASIN(1.0 / y) - y + SQRT(y * y - 1.0) ));
  else *v = 1.0 - y * (PI - 2.0 * y);
}
void Dpower(double *x, model *cov, double *v){
  double alpha = P0(POW_ALPHA), y = *x;
  *v = 0.0;
  if (y < 1.0) *v = - alpha * POW(1.0 - y, alpha - 1.0);
}
int checkpower(model *cov) {
  double alpha = P0(POW_ALPHA);
  double dim = 2.0 * alpha - 1.0;
  set_maxdim(OWN, 0, ISNAN(dim) || dim >= INFDIM ? INFDIM-1 : (int) dim);
  cov->monotone = alpha >= (double) ((OWNLOGDIM(0) / 2) + 1)
    ? COMPLETELY_MON : NORMAL_MIXTURE;
  RETURN_NOERROR;
}


// range definition:
// 0: min, theory, 1:max, theory
// 2: min, practically 3:max, practically
void rangepower(model *cov, range_type *range){
 bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));
  int dim = OWNLOGDIM(0);
   range->min[POW_ALPHA] = tcf
    ? (double) ((dim / 2) + 1) // Formel stimmt so! Nicht aendern!
    : 0.5 * (double) (dim + 1);
  range->max[POW_ALPHA] = RF_INF;
  range->pmin[POW_ALPHA] = range->min[POW_ALPHA];
  range->pmax[POW_ALPHA] = 20.0;
  range->openmin[POW_ALPHA] = false;
  range->openmax[POW_ALPHA] = true;
}


/* qexponential -- derivative of exponential */
#define QEXP_ALPHA 0
void qexponential(double *x, model *cov, double *v){
  double 
    alpha = P0(QEXP_ALPHA),
    y = EXP(-*x);
  *v = y * (2.0  - alpha * y) / (2.0 - alpha);
}
void Inverseqexponential(double *x, model *cov, double *v){
  double alpha = P0(QEXP_ALPHA);
  *v = -LOG( (1.0 - SQRT(1.0 - alpha * (2.0 - alpha) * *x)) / alpha);
} 
void Dqexponential(double *x, model *cov, double *v) {
  double 
    alpha = P0(QEXP_ALPHA), 
    y = EXP(-*x);
  *v = y * (alpha * y - 1.0) * 2.0 / (2.0 - alpha);
}
void rangeqexponential(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[QEXP_ALPHA] = 0.0;
  range->max[QEXP_ALPHA] = 1.0;
  range->pmin[QEXP_ALPHA] = 0.0;
  range->pmax[QEXP_ALPHA] = 1.0;
  range->openmin[QEXP_ALPHA] = false;
  range->openmax[QEXP_ALPHA] = false;
}



#define SINEPOWER_ALPHA 0
void sinepower(double *x, model *cov, double *v){
  double alpha = P0(SINEPOWER_ALPHA);  // Auslesen des Parameters aus cov  
  double y = *x;
  if (y < 0.0 || y >= PI) BUG;
  *v = 1.0 - POW(SIN(y * 0.5), alpha);
  return;
}


/* range des Parameters alpha */
void rangesinepower(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SINEPOWER_ALPHA] = 0.0;  // true range
  range->max[SINEPOWER_ALPHA] = 2.0;   
  range->pmin[SINEPOWER_ALPHA] = 0.001; // practical range (assumed)
  range->pmax[SINEPOWER_ALPHA] = 2.0; 
  range->openmin[SINEPOWER_ALPHA] = true;  // lower endpoint included 
  range->openmax[SINEPOWER_ALPHA] = false; // upper endpoint excluded  
}


/* spherical model */ 
void spherical(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x;
  *v = 0.0;
  if (y < 1.0) *v = 1.0 + y * 0.5 * (y * y - 3.0);
}
// void Inversespherical(model *cov){ return 1.23243208931941;}
void TBM2spherical(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x , y2 = y * y;
  if (y>1.0) *v = (1.0- 0.75 * y * ((2.0 - y2) * ASIN(1.0/y) + SQRT(y2 -1.0)));
  else *v = (1.0 - 0.375 * PI * y * (2.0 - y2));
}
void Dspherical(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  double y = *x;
  *v = (y >= 1.0) ? 0.0 : 1.5 * (y * y - 1.0);
}

void DDspherical(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  *v = (*x >= 1.0) ? 0.0 : 3 * *x;
}

int structspherical(model *cov, model **newmodel) {
  return structCircSph( cov, newmodel, 3);
}
void dospherical(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 
  // todo diese void do... in Primitive necessary??
  //  mppinfotype *info = &(s->mppinfo);
  //info->radius = cov->mpp.refradius;
}


double alphaIntSpherical(int a) {
  // int r^a C(r) \D r
  double A = (double) a;
  return 1.0 / (A + 1.0) - 1.5 / (A + 2.0) + 0.5 / (A + 4.0);
}

int initspherical(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int  dim = OWNLOGDIM(0);

  if (hasAnyEvaluationFrame(cov)) {
    if (cov->mpp.moments >= 1) SERR("too high moments required");
  }

  else if (hasSmithFrame(cov)) {    
    assert(cov->mpp.maxheights[0] == 1.0);
    if (cov->mpp.moments >= 1) {
      cov->mpp.mM[1] = cov->mpp.mMplus[1] =  
	SurfaceSphere(dim - 1, 1) * alphaIntSpherical(dim - 1);
    }
  }

  else if (hasRandomFrame(cov)) { RETURN_NOERROR; }

  else ILLEGAL_FRAME;

 RETURN_NOERROR;
}


/* stein space-time model */
#define STEIN_NU 0
#define STEIN_Z 1
void kappaSteinST1(int i, model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = (i == STEIN_NU) ? 1 : (i==STEIN_Z) ? OWNLOGDIM(0) - 1 : -1;
}
void SteinST1(double *x, model *cov, double *v){
/* 2^(1-nu)/Gamma(nu) [ h^nu K_nu(h) - 2 * tau (x T z) t h^{nu-1} K_{nu-1}(h) /
   (2 nu + d + 1) ]

   \|tau z \|<=1 hence \tau z replaced by z !!
*/
  int d,
   dim = OWNLOGDIM(0),
    time = dim - 1;
  double logconst, hz, y,
    nu = P0(STEIN_NU),
    *z=P(STEIN_Z);

  hz = 0.0;
  y = x[time] * x[time];
  for (d=0; d<time; d++) {
    y += x[d] * x[d];
    hz += x[d] * z[d];
  }
  
  if ( y==0.0 ) *v = 1.0;
  else {
    double bk[MATERN_NU_THRES];
    y = SQRT(y);
    logconst = (nu - 1.0) * LOG(0.5 * y)  - QVALUE;
    *v =  y * EXP(logconst + LOG(bessel_k_ex(y, nu, 2.0, bk)) - y)
      - 2.0 * hz * x[time] * EXP(logconst +
				 LOG(bessel_k_ex(y, nu - 1.0, 2.0, bk)) - y)
      / (2.0 * nu + dim);
  }

}


int checkSteinST1(model *cov) {  
  double nu = P0(STEIN_NU), *z= P(STEIN_Z), absz;
   int d, 
    spatialdim=OWNLOGDIM(0) - 1;

  for (d=0; d<= Nothing; d++) cov->pref[d] *= (nu < BesselUpperB[d]);
  if (nu >= 2.5) cov->pref[CircEmbed] = 2;
  if (spatialdim < 1) 
    SERR("dimension of coordinates, including time, must be at least 2");
  if (nu > MATERN_NU_THRES) SERR1("'nu'>%d is too large for precise returns",
				  MATERN_NU_THRES)
  	
  for (absz=0.0, d=0;  d<spatialdim; d++)  absz += z[d] * z[d];
  if (ISNAN(absz))
    SERR("currently, components of z cannot be estimated by MLE, so NA's are not allowed");
  if (absz > 1.0 + UNIT_EPSILON && !GLOBAL_UTILS->basic.skipchecks) {
    SERR("||z|| must be less than or equal to 1");
  }
  if (cov->q == NULL) {
    EXTRA_Q;
    initSteinST1(cov, NULL);
  }
  RETURN_NOERROR;
}
double densitySteinST1(double *x, model *cov) {
  double x2, wz, dWM, nu = P0(STEIN_NU), 
    *z=P(STEIN_Z);
  int d,
    dim = PREVLOGDIM(0),
    spatialdim = dim - 1;

  x2 = x[spatialdim] * x[spatialdim]; // v^2
  wz = 0.0;
  for (d=0; d<spatialdim; d++) {
    x2 += x[d] * x[d];  // w^2 + v^2
    wz += x[d] * z[d];
  }
  dWM = EXP(QVALUE2 - QVALUE3 * LOG(x2 + 1.0));
  return (1.0 + 2.0 * wz * x[spatialdim] + x2) * dWM
    / (2.0 * nu + (double) dim + 1.0);
}


int initSteinST1(model *cov, gen_storage *s) {
  int dim = PREVLOGDIM(0);
  double nu = P0(STEIN_NU);
  QVALUE = lgammafn(nu);
  QVALUE2 = QVALUE - lgammafn(nu +  0.5 * dim) - dim * M_LN_SQRT_PI;
  QVALUE3 = nu + dim; // factor

  if (HAS_SPECTRAL_FRAME(cov)) {
    spec_properties *cs = &(s->spec);
    cs->density = densitySteinST1;
    return search_metropolis(cov, s);
  }

  RETURN_NOERROR;
}

void spectralSteinST1(model *cov, gen_storage *S, double *e){
  metropolis(cov, S, e);
}

void rangeSteinST1(model VARIABLE_IS_NOT_USED *cov, range_type *range){  
  range->min[STEIN_NU] = 1.0;
  range->max[STEIN_NU] = RF_INF;
  range->pmin[STEIN_NU] = 1e-10;
  range->pmax[STEIN_NU] = 10.0;
  range->openmin[STEIN_NU] = true;
  range->openmax[STEIN_NU] = true;
 
  range->min[STEIN_Z] = RF_NEGINF;
  range->max[STEIN_Z] = RF_INF;
  range->pmin[STEIN_Z] = - 10.0;
  range->pmax[STEIN_Z] = 10.0;
  range->openmin[STEIN_Z] = true;
  range->openmax[STEIN_Z] = true;
}


/* wave */
void wave(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  double y = *x;
  *v = (y==0.0) ? 1.0 : (y==RF_INF) ? 0 : SIN(y) / y;
}
void Inversewave(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  *v = (*x==0.05) ? 1.0/0.302320850755833 : RF_NA;
}
int initwave(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  if (HAS_SPECTRAL_FRAME(cov)) {
    return (OWNLOGDIM(0) <= 2) ? NOERROR : ERRORFAILED;
  } 

  else if (hasRandomFrame(cov)) { RETURN_NOERROR; }

  else ILLEGAL_FRAME;

}
void spectralwave(model *cov, gen_storage *S, double *e) { 
  spectral_storage *s = &(S->Sspectral);
  /* see Yaglom ! */
  double x;  
  x = UNIFORM_RANDOM; 
  E12(s, PREVLOGDIM(0), SQRT(1.0 - x * x), e);
}



#define USER_TYPE 0
#define USER_DOM 1
#define USER_ISO 2
#define USER_VDIM 3
#define USER_BETA 4
#define USER_COORD 5
#define USER_FCTN 6
#define USER_FST 7
#define USER_SND 8
#define USER_ENV 9
#define USER_LAST USER_ENV
void evaluateUser(double *x, double *y, bool Time, model *cov,
		  sexp_type *which, double *Res) {
  SEXP  res,
    env = PENV(USER_ENV)->sexp; 
  int i,
    vdim = VDIM0 * VDIM1,
    ncol = cov->ncol[USER_BETA],
    n = OWNXDIM(0);
  double *beta = P(USER_BETA);  
 
  assert(which != NULL);
  if (cov->nrow[USER_COORD] >= 2 && PINT(USER_COORD)[1] != -2) {
    if (Time) {
      addVariable( (char *) "T", x + (--n), 1, 1, env); 
    }
    switch (n) {
    case 3 : addVariable( (char *) "z", x+2, 1, 1, env);
      FALLTHROUGH_OK;
    case 2 : addVariable( (char *) "y", x+1, 1, 1, env);
      FALLTHROUGH_OK;
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


void kappaUser(int i, model *cov, int *nr, int *nc){
  *nc = *nr = i < DefList[COVNR].kappas ? 1 : -1;
  if (i == USER_VDIM) *nr = SIZE_NOT_DETERMINED;
  if (i == USER_COORD) *nr=SIZE_NOT_DETERMINED;
  if (i == USER_BETA) *nr=*nc=SIZE_NOT_DETERMINED;
}

void User(double *x, model *cov, double *v){
  evaluateUser(x, NULL, Loc(cov)->Time, cov, PLANG(USER_FCTN), v);
}
void UserNonStat(double *x, double *y, model *cov, double *v){
  evaluateUser(x, y, false, cov, PLANG(USER_FCTN), v);
}
void DUser(double *x, model *cov, double *v){
  evaluateUser(x, NULL, Loc(cov)->Time, cov, PLANG(USER_FST), v);
}
void DDUser(double *x, model *cov, double *v){
  evaluateUser(x, NULL, Loc(cov)->Time, cov, PLANG(USER_SND), v);
}


int checkUser(model *cov){
  defn *C = DefList + COVNR;

  kdefault(cov, USER_DOM, XONLY);
  if (PisNULL(USER_ISO)) {
    Types type = (Types) P0INT(USER_TYPE);
    if (isVariogram(type)) kdefault(cov, USER_ISO, ISOTROPIC);
    else if (isProcess(type) || isShape(type)) 
      kdefault(cov, USER_ISO, CARTESIAN_COORD);
    else SERR1("type='%.50s' not allowed", TYPE_NAMES[type]);
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
    *variab = PINT(USER_COORD), //codierung variab=1:x, 2:y 3:z, 4:T
    nvar = cov->nrow[USER_COORD],
    *pref = cov->pref;


  if (nvar < 1) SERR("variables not of the required form ('x', 'y', 'z', 'T')");
  
  if (OWNDOM(0) != (domain_type) *dom) 
    SERR2("wrong domain (requ='%.50s'; provided='%.50s')", 
	  DOMAIN_NAMES[OWNDOM(0)], DOMAIN_NAMES[*dom]);
  if (OWNISO(0) != (isotropy_type) *iso)
    SERR2("wrong isotropy assumption (requ='%.50s'; provided='%.50s')",
	  ISO_NAMES[OWNISO(0)], ISO_NAMES[*iso]);

  if (PENV(USER_ENV) == NULL) BUG;

  if (!PisNULL(USER_BETA)) {
    if (vdim == NULL) kdefault(cov, USER_VDIM, nrow[USER_BETA]);
    else if (*vdim != nrow[USER_BETA]) 
      SERR4("number of columns of '%.50s' (=%d) does not equal '%.50s' (=%d)",
	    KNAME(USER_BETA), nrow[USER_BETA], KNAME(USER_VDIM), *vdim);
    VDIM0 = nrow[USER_BETA];
    VDIM1 = 1;
  } else {
    if (PisNULL(USER_VDIM)) kdefault(cov, USER_VDIM, 1);  
    VDIM0 = P0INT(USER_VDIM);
    VDIM1 = cov->nrow[USER_VDIM] == 1 ? 1 : PINT(USER_VDIM)[1];
  }
  if (cov->nrow[USER_VDIM] > 2) SERR1("'%.50s' must be a scalar or a vector of length 2", KNAME(USER_VDIM));

  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);

  if (variab[0] != 1 && (variab[0] != 4 || nvar > 1)) SERR("'x' not given");
  if (nvar > 1) {
    variab[1] = std::abs(variab[1]); // integer ABS!
    //                              it is set to its negativ value
    //                              below, when a kernel is defined
    if (variab[1] == 3) SERR("'z' given but not 'y'"); 
  } 
  Time = variab[nvar-1] == 4;

  if (((nvar >= 3 || variab[nvar-1] == 4) && GLOBAL.coords.xyz_notation==False)
      //  ||
      //  (nrow[USER_COORD] == 1 && !ISNA_INT(GLOBAL.coords.xyz_notation) 
      //  && GLOBAL.coords.xyz_notation)
      )
    SERR("mismatch of indicated xyz-notation");

  if (Time xor Loc(cov)->Time)
    SERR("If 'T' is given, it must be given both as coordinates and as variable 'T' in the function 'fctn'");

  if ((nvar > 2 || (nvar == 2 && variab[1] != 2)) && isKernel(OWN))
    SERR1("'%.50s' not valid for anisotropic models", coords[COORDS_XYZNOTATION]);
  
  if (nvar == 2 && variab[1] == 2) {
    // sowohl nonstat also auch x, y Schreibweise moeglich
    if (GLOBAL.coords.xyz_notation == Nan)
      SERR1("'%.50s' equals 'NA', but should be a logical value.", 
	   coords[COORDS_XYZNOTATION]);     
    //if (isKernel(OWN) && GLOBAL.coords.xyz_notation==2) // RFcov
	// SERR1("'%.50s' is not valid for anisotropic models", 
    //  coords[COORDS_XYZNOTATION]);
  }
 
  if (nvar > 1) {
    if (isXonly(OWN)) {
      if (equalsIsotropic(OWNISO(0))) {
	SERR("two many variables given for motion invariant function");
      } else if (isSpaceIsotropic(OWN) && nvar != 2)
	SERR("number of variables does not match a space-isotropic model");
    } else {
      if (nvar == 2 && variab[1] == 2) variab[1] = -2;//i.e. non domain (kernel)
    }
  }

  if (GLOBAL.coords.xyz_notation != Nan) {    
    if (//(GLOBAL.coords.xyz_notation == 2 && isKernel(OWN)) ||
	((nvar > 2 || (nvar == 2 && isXonly(OWN))) && variab[1] == -2)) {
      SERR("domain assumption, model and coordinates do not match.");
    }
  }

  if ((OWNXDIM(0) == 4 && !Loc(cov)->Time) || OWNXDIM(0) > 4)
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
      if (snd) SERR2("'%.50s' given but not '%.50s'",
		     KNAME(USER_SND), KNAME(USER_FST));
    }
  } else { // !fctn
    if (fst || snd) 
      SERR3("'%.50s' or '%.50s' given but not '%.50s'", 
	    KNAME(USER_SND), KNAME(USER_FST), KNAME(USER_FCTN));
  }
   RETURN_NOERROR;
}


Types TypeUser(Types required, model *cov,
	       isotropy_type VARIABLE_IS_NOT_USED i) {
  if (PisNULL(USER_TYPE)) return BadType;
  Types type = (Types) P0INT(USER_TYPE);
  if (! (isShape(type) || equalsRandom(type)) ) return BadType;
  return TypeConsistency(required, type);
}


bool setUser(model *cov){
  // A  PMI0(cov);
  Types
    type = PisNULL(USER_TYPE) ? BadType : (Types) P0INT(USER_TYPE);
  domain_type
    dom = PisNULL(USER_DOM) ? DOMAIN_MISMATCH : (domain_type) P0INT(USER_DOM);
  isotropy_type
    iso = PisNULL(USER_ISO) ? ISO_MISMATCH : (isotropy_type) P0INT(USER_ISO);
  int
    xdim = cov->nrow[USER_COORD];

  //    vdim = PisNULL(USER_VDIM) ? 0 : P0INT(USER_VDIM);
  //printf("isset = %d %d\n", PREV_INITIALISEDXX, PREV_INITIALISED);
 
  isotropy_type previso = CONDPREVISO(0); 
  set_system(OWN, 0, isFixed(previso) ? PREVLOGDIM(0) : xdim /* xdim only dummy! */,
	     xdim, xdim, type, dom, iso);
  return true;
}


bool allowedDuser(model *cov) {
  if (PisNULL(USER_DOM)) return allowedDtrue(cov);
  
  bool *D = cov->allowedD;
  for (int i=FIRST_DOMAIN; i<LAST_DOMAINUSER; D[i++]=false);
  D[(domain_type) P0INT(USER_DOM)] = true;
  return false;
}

bool allowedIuser(model *cov) {
  if (PisNULL(USER_ISO)) return allowedDtrue(cov);
 
  bool *I = cov->allowedI;
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  I[(domain_type) P0INT(USER_ISO)] = true;
  return false;
}



void rangeUser(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[USER_TYPE] = TcfType;
  range->max[USER_TYPE] = TrendType;
  range->pmin[USER_TYPE] = TcfType;
  range->pmax[USER_TYPE] = TrendType;
  range->openmin[USER_TYPE] = false;
  range->openmax[USER_TYPE] = false;

  range->min[USER_DOM] = FIRST_DOMAIN;
  range->max[USER_DOM] = LAST_DOMAINUSER;
  range->pmin[USER_DOM] = FIRST_DOMAIN;
  range->pmax[USER_DOM] = LAST_DOMAINUSER;
  range->openmin[USER_DOM] = false;
  range->openmax[USER_DOM] = false;

  range->min[USER_ISO] = FIRST_CARTESIAN;
  range->max[USER_ISO] = LAST_CARTESIAN;
  range->pmin[USER_ISO] = FIRST_CARTESIAN;
  range->pmax[USER_ISO] = LAST_CARTESIAN;
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
  range->pmin[USER_BETA] = - 1e5;
  range->pmax[USER_BETA] = 1e5;
  range->openmin[USER_BETA] = true;
  range->openmax[USER_BETA] = true;

  range->min[USER_COORD] = -2;
  range->max[USER_COORD] = 4;
  range->pmin[USER_COORD] = 1;
  range->pmax[USER_COORD] = 4;
  range->openmin[USER_COORD] = false;
  range->openmax[USER_COORD] = false;
  
#define user_n 4
  int idx[user_n] = {USER_FCTN, USER_FST, USER_SND, USER_ENV};
  for (int i=0; i<user_n; i++) {
    int j = idx[i];
    range->min[j] = RF_NAN;
    range->max[j] = RF_NAN;
    range->pmin[j] = RF_NAN;
    range->pmax[j] = RF_NAN;
    range->openmin[j] = false;
    range->openmax[j] = false; 
  }

  int kappas = DefList[COVNR].kappas;
  for (int i=USER_LAST + 1; i<kappas; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = 1e10;
    range->pmax[i] = - 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true; 
  }
}

